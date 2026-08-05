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
#' @param genome One of \code{"mm10"} (default) or \code{"hg38"}. Selects the gene
#'   annotation (\code{EnsDb.Mmusculus.v79} / \code{EnsDb.Hsapiens.v86}) and genome
#'   sequence (\code{BSgenome.Mmusculus.UCSC.mm10} / \code{BSgenome.Hsapiens.UCSC.hg38})
#'   used for scaffold-chromosome filtering and the \code{ChromatinAssay} annotation.
#'   Both are optional (Suggests) packages -- only the one matching your chosen
#'   \code{genome} needs to be installed, and you'll get an actionable error naming
#'   it if it's missing.
#' @param object_names Optional character vector of names for the returned list,
#'   the same length as \code{data_dirs}. \code{NULL} (default) uses
#'   \code{basename(data_dirs)}.
#' @param peakwidths_max Max peak width for finding combined peaks
#' @param peakwidths_min Min peak width for finding combined peaks
#' @param passed_filters_value Min value for filtering cells based on
#' passed_filters column
#' @return A list of Seurat objects
#' @export

CreateATACObjects <-
  function(data_dirs, add_treatment = FALSE, treatment = NULL,
           genome = c("mm10", "hg38"), object_names = NULL,
           peakwidths_max = 10000, peakwidths_min = 20,
           passed_filters_value = 500) {

    genome <- match.arg(genome)

    if(add_treatment == FALSE & is.null(treatment) == FALSE) {
      stop('\n\n  Error: Treatment vector was added, but add_treatment set to FALSE.\nSet add_treatment to TRUE before proceeding.')
    }
    if (!is.null(object_names) && length(object_names) != length(data_dirs)) {
      stop('`object_names` must be the same length as `data_dirs` (', length(data_dirs),
          '), got ', length(object_names), '.')
    }

    # ---- Resolve genome-specific annotation/sequence packages -----------------
    # Both mouse and human packages are optional (Suggests): only the one
    # matching the requested `genome` needs to actually be installed.
    genome_pkgs <- switch(
      genome,
      mm10 = list(ensdb = 'EnsDb.Mmusculus.v79', bsgenome = 'BSgenome.Mmusculus.UCSC.mm10'),
      hg38 = list(ensdb = 'EnsDb.Hsapiens.v86',   bsgenome = 'BSgenome.Hsapiens.UCSC.hg38')
    )
    if (!requireNamespace(genome_pkgs$ensdb, quietly = TRUE)) {
      stop("'", genome_pkgs$ensdb, "' is required for genome = '", genome, "'. ",
          "Install with: BiocManager::install('", genome_pkgs$ensdb, "')")
    }
    if (!requireNamespace(genome_pkgs$bsgenome, quietly = TRUE)) {
      stop("'", genome_pkgs$bsgenome, "' is required for genome = '", genome, "'. ",
          "Install with: BiocManager::install('", genome_pkgs$bsgenome, "')")
    }
    # Both packages export an object of the same name as the package itself
    # (the standard Bioconductor EnsDb/BSgenome annotation-package convention).
    ensdb_obj    <- getExportedValue(genome_pkgs$ensdb, genome_pkgs$ensdb)
    bsgenome_obj <- getExportedValue(genome_pkgs$bsgenome, genome_pkgs$bsgenome)

    message(sprintf('--- Reading peak sets (%d directories) ---', length(data_dirs)))
    # Use lapply to read the peak sets for each sample
    peak_data_list <- lapply(data_dirs, function(dir) {
      # Read peaks
      peak_data <- read.table(file = paste(dir, '/outs/peaks.bed', sep = ''),
                              col.names = c("chr", "start", "end"))
      # Make GRanges objects
      gr <- makeGRangesFromDataFrame(peak_data)
    })

    message('--- Building combined peak set ---')
    # Create combined peak set
    suppressWarnings(for (i in 2:length(peak_data_list)) {
      combined.peaks <- GenomicRanges::reduce(c(peak_data_list[[1]], peak_data_list[[i]]))
    })
    peakwidths <- width(combined.peaks)
    combined.peaks <- combined.peaks[peakwidths < peakwidths_max &
                                       peakwidths > peakwidths_min]
    message(sprintf('  Peaks within width range [%d, %d]: %d',
                    peakwidths_min, peakwidths_max, length(combined.peaks)))

    message('--- Removing peaks on scaffolds (keeping main chromosomes) ---')
    #remove scaffolds not in genome
    main.chroms <- standardChromosomes(bsgenome_obj)
    keep.peaks <- as.logical(seqnames(granges(combined.peaks)) %in% main.chroms)
    combined.peaks <- combined.peaks[keep.peaks, ]
    message(sprintf('  Peaks after scaffold removal: %d', length(combined.peaks)))

    message(sprintf('--- Loading gene annotations from %s ---', genome_pkgs$ensdb))
    # extract gene annotations from EnsDb
    annotations <- Signac::GetGRangesFromEnsDb(ensdb = ensdb_obj)

    # change to UCSC style since the data was mapped to hg19
    seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
    genome(annotations) <- genome

    message('--- Building Seurat ATAC objects per sample ---')
    # Create Seurat objects
    seurat_objects <- lapply(seq_along(data_dirs), function(idx) {
      dir <- data_dirs[[idx]]
      message(sprintf('  Building object %d of %d: %s',
                      idx, length(data_dirs), basename(dir)))

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

    names(seurat_objects) <- if (!is.null(object_names)) object_names else basename(data_dirs)

    message('--- Generating ATAC QC plots ---')
    obj <- merge(seurat_objects[[1]], seurat_objects[-1])

    pct_reads_in_peaks.plot <- ggplot(obj@meta.data,
                                      aes(orig.ident, pct_reads_in_peaks)) + geom_boxplot() + Ol_Reliable()
    peak_region_fragments.plot <- ggplot(obj@meta.data,
                                         aes(orig.ident, peak_region_fragments)) + geom_boxplot() + Ol_Reliable()
    TSS.enrichment.plot <- ggplot(obj@meta.data,
                                  aes(orig.ident, TSS.enrichment)) + geom_boxplot() + Ol_Reliable()
    blacklist_ratio.plot <- ggplot(obj@meta.data,
                                   aes(orig.ident, blacklist_ratio)) + geom_boxplot() + Ol_Reliable()
    nucleosome_signal.plot <- ggplot(obj@meta.data,
                                     aes(orig.ident, nucleosome_signal)) + geom_boxplot() + Ol_Reliable()

    print(pct_reads_in_peaks.plot + peak_region_fragments.plot +
            TSS.enrichment.plot +
            blacklist_ratio.plot + nucleosome_signal.plot +
            patchwork::plot_layout(ncol = 3))

    # Add a column to metadata to specify treatment
    if (is.null(treatment) == FALSE ){
      message('--- Adding Treatment metadata column ---')
      seurat_objects <- setNames(lapply(seq_along(seurat_objects), function(i) {
        seurat_obj <- seurat_objects[[i]]
        seurat_obj[["Treatment"]] <- treatment[i]
        return(seurat_obj)
      }), names(seurat_objects))
    }

    return(seurat_objects)

  }
