#' Create Seurat ATAC Objects
#'
#' This function creates multiple Seurat objects. It takes a list of directories
#' as input. In each directory, there should be at least the following files:
#' peaks.bed, singlecell.csv and fragments.tsv.gz. You can directly give the
#' output folder from cellranger into this function. While reading in the data,
#' the function creates a common peak list from the samples. Following this,
#' peaks on scaffolds and not the major chromosomes are removed.
#'
#'
#' @param data_dirs Path to the directories with ATAC data
#' @param add_treatment Whether to add a treatment column to metadata
#' @param treatment Treatment metadata value (vector)
#' @param peakwidths_max Max peak width for finding combined peaks
#' @param peakwidths_min Min peak width for finding combined peaks
#' @param passed_filters_value Min value for filtering cells based on
#' passed_filters column
#' @return A list of Seurat objects
#' @export

CreateATACObjects <-
  function(data_dirs, add_treatment = FALSE, treatment = NULL,
           peakwidths_max = 10000, peakwidths_min = 20,
           passed_filters_value = 500) {

    if(add_treatment == FALSE & is.null(treatment) == FALSE) {
      stop('\n\n  Error: Treatment vector was added, but add_treatment set to FALSE.\nSet add_treatment to TRUE before proceeding.')
    }

    # Use lapply to read the peak sets for each sample
    peak_data_list <- lapply(data_dirs, function(dir) {
      # Read peaks
      peak_data <- read.table(file = paste(dir, '/outs/peaks.bed', sep = ''),
                              col.names = c("chr", "start", "end"))
      # Make GRanges objects
      gr <- makeGRangesFromDataFrame(peak_data)
    })

    # Create combined peak set
    suppressWarnings(for (i in 2:length(peak_data_list)) {
      combined.peaks <- GenomicRanges::reduce(c(peak_data_list[[1]], peak_data_list[[i]]))
    })
    peakwidths <- width(combined.peaks)
    combined.peaks <- combined.peaks[peakwidths < peakwidths_max &
                                       peakwidths > peakwidths_min]
    #remove scaffolds not in genome
    main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
    keep.peaks <- as.logical(seqnames(granges(combined.peaks)) %in% main.chroms)
    combined.peaks <- combined.peaks[keep.peaks, ]

    # extract gene annotations from EnsDb
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

    # change to UCSC style since the data was mapped to hg19
    seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
    genome(annotations) <- "mm10"

    # Create Seurat objects
    seurat_objects <- lapply(data_dirs, function(dir) {

      # Load metadata for each sample
      md <- read.table(file = paste(dir, "/outs/singlecell.csv", sep = ''), sep = ",", header = TRUE, row.names = 1)[-1, ] # remove the first row
      md <- md[md$passed_filters > passed_filters_value, ]

      # Create fragment objects
      frag.obj <- CreateFragmentObject(path = paste(dir, '/outs/fragments.tsv.gz', sep = ''),
                                       cells = rownames(md))

      # Create Feature matrix objects
      counts <- FeatureMatrix(
        fragments = frag.obj,
        features = combined.peaks,
        cells = rownames(md)
      )

      #Create chromatin assay and final object with QC metics
      assay <- CreateChromatinAssay(counts, fragments = frag.obj)
      seurat.obj <- CreateSeuratObject(assay, assay = "ATAC", meta.data=md,
                                       project = basename(dir))

      # add the gene information to the object
      Annotation(seurat.obj) <- annotations

      seurat.obj <- NucleosomeSignal(seurat.obj)
      seurat.obj$nucleosome_group <- ifelse(seurat.obj$nucleosome_signal > 4,
                                            'NS > 4', 'NS < 4')
      seurat.obj <- TSSEnrichment(seurat.obj)
      seurat.obj$pct_reads_in_peaks <- seurat.obj$peak_region_fragments /
        seurat.obj$passed_filters * 100
      seurat.obj$blacklist_ratio <- seurat.obj$blacklist_region_fragments /
        seurat.obj$peak_region_fragments

      return(seurat.obj)
    })

    names(seurat_objects) <- basename(data_dirs)

    obj <- merge(seurat_objects[[1]], seurat_objects[-1])

    pct_reads_in_peaks.plot <- ggplot(obj@meta.data,
                                      aes(orig.ident, pct_reads_in_peaks)) + geom_boxplot()
    peak_region_fragments.plot <- ggplot(obj@meta.data,
                                         aes(orig.ident, peak_region_fragments)) + geom_boxplot()
    TSS.enrichment.plot <- ggplot(obj@meta.data,
                                  aes(orig.ident, TSS.enrichment)) + geom_boxplot()
    blacklist_ratio.plot <- ggplot(obj@meta.data,
                                   aes(orig.ident, blacklist_ratio)) + geom_boxplot()
    nucleosome_signal.plot <- ggplot(obj@meta.data,
                                     aes(orig.ident, nucleosome_signal)) + geom_boxplot()

    print(pct_reads_in_peaks.plot + peak_region_fragments.plot +
            TSS.enrichment.plot +
            blacklist_ratio.plot + nucleosome_signal.plot +
            patchwork::plot_layout(ncol = 3))

    # Add a column to metadata to specify treatment
    if (is.null(treatment) == FALSE ){
      seurat_objects <- setNames(lapply(seq_along(seurat_objects), function(i) {
        seurat_obj <- seurat_objects[[i]]
        seurat_obj[["Treatment"]] <- treatment[i]
        return(seurat_obj)
      }), names(seurat_objects))
    }

    return(seurat_objects)

  }
