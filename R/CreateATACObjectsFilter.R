#' Create and Filter Seurat ATAC Objects
#'
#' This function creates multiple Seurat objects and filters them. It takes a
#' list of directories as input. In each directory, there should be at least the
#' following files: peaks.bed, singlecell.csv and fragments.tsv.gz. While
#' reading in the data, the function creates a common peak list from the
#' samples. Following this, peaks on scaffolds and not the major chromosomes are
#' removed. You can directly give the output folder from cellranger into this
#' function. The objects are then automatically filtered based on hard cutoffs.
#' This can function can be run interactively to choose cutoffs without having
#' to re-run it.
#'
#' @param data_dirs Path to the directories
#' @return A list of filtered Seurat objects
#' @export
CreateATACObjectsFilter <-
  function(data_dirs,  treatment = NULL, filter = TRUE, interactive = FALSE,
           peak_region_fragments_min = 3000, peak_region_fragments_max = 100000,
           pct_reads_in_peaks_min = 40, blacklist_ratio_max = 0.025,
           nucleosome_signal_max = 4, TSS.enrichment_min = 2,
           peak_region_max = 3000, peakwidths_max = 10000,
           peakwidths_min = 20, passed_filters_value = 500) {

    if (filter == FALSE & interactive == TRUE) {
      stop("Error: Set filter=TRUE to use interactiver, otherwise set
  interactive=FALSE to skip filtering")
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
      combined.peaks <- reduce(c(peak_data_list[[1]], peak_data_list[[i]]))
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
    if (is.null(treatment) == FALSE) {
      seurat_objects <- setNames(lapply(seq_along(seurat_objects), function(i) {
        seurat_obj <- seurat_objects[[i]]
        seurat_obj[["Treatment"]] <- treatment[i]
        return(seurat_obj)
      }), names(seurat_objects))
    }

    if (filter == TRUE & interactive == FALSE) {
      subsetted_objs <- lapply(seurat_objects, function(obj) {
        # Subset based on the threshold
        subset(obj, subset =
                 peak_region_fragments > peak_region_fragments_min &
                 peak_region_fragments < peak_region_fragments_max &
                 pct_reads_in_peaks > pct_reads_in_peaks_min &
                 blacklist_ratio < blacklist_ratio_max &
                 nucleosome_signal < nucleosome_signal_max &
                 TSS.enrichment > TSS.enrichment_min)
      })

      obj <- merge(subsetted_objs[[1]], subsetted_objs[-1])

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

      return(subsetted_objs)

    } else if (filter == TRUE & interactive == TRUE) {
      # Ask user for thresholds interactively
      for (param in c("min pct_reads_in_peaks", "min peak_region_fragments", "max peak_region_fragments", "min TSS.enrichment", "max blacklist_ratio", "max nucleosome_signal")) {
        use_quantile <- readline(prompt = paste("Do you want to use quantile for subsetting", param, "? (yes/no): "))
        if (tolower(use_quantile) == "yes") {
          threshold <- as.numeric(readline(prompt = paste("Enter quantile threshold for", param, " (0 to 1): ")))
          seurat_objects <- lapply(seurat_objects, function(obj) {
            if (param == "min pct_reads_in_peaks") {
              subset(obj, subset = pct_reads_in_peaks > quantile(obj$pct_reads_in_peaks, threshold))
            } else if (param == "min peak_region_fragments") {
              subset(obj, subset = peak_region_fragments > quantile(obj$peak_region_fragments, threshold))
            } else if (param == "max peak_region_fragments") {
              subset(obj, subset = peak_region_fragments < quantile(obj$peak_region_fragments, threshold))
            } else if (param == "min TSS.enrichment") {
              subset(obj, subset = TSS.enrichment > quantile(obj$TSS.enrichment, threshold))
            } else if (param == "max blacklist_ratio") {
              subset(obj, subset = blacklist_ratio < quantile(obj$blacklist_ratio, threshold))
            } else if (param == "max nucleosome_signal") {
              subset(obj, subset = nucleosome_signal < quantile(obj$nucleosome_signal, threshold))
            }
          })
        } else if (tolower(use_quantile) == "no") {
          threshold <- as.numeric(readline(prompt = paste("Enter threshold for", param, ": ")))
          seurat_objects <- lapply(seurat_objects, function(obj) {
            if (param == "min pct_reads_in_peaks") {
              subset(obj, subset = pct_reads_in_peaks > threshold)
            } else if (param == "min peak_region_fragments") {
              subset(obj, subset = peak_region_fragments > threshold)
            } else if (param == "max peak_region_fragments") {
              subset(obj, subset = peak_region_fragments < threshold)
            } else if (param == "min TSS.enrichment") {
              subset(obj, subset = TSS.enrichment > threshold)
            } else if (param == "max blacklist_ratio") {
              subset(obj, subset = blacklist_ratio < threshold)
            } else if (param == "max nucleosome_signal") {
              subset(obj, subset = nucleosome_signal < threshold)
            }
          })
        }
      }


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

      print('Saving Filtered Objects')

      saveRDS(seurat_objects, 'seurat_objects_filtered.rds')

      return(seurat_objects)
    } else {
      return(seurat_objects)
    }
  }
