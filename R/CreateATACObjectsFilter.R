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
#' @param data_dirs Path to the directories with ATAC data
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
#' @param interactive Whether to run the filtering step in interactive mode
#' @param peak_region_fragments_min Min peak region fragments value for
#' filtering
#' @param peak_region_fragments_max Max peak region fragments value for
#' filtering
#' @param pct_reads_in_peaks_min Min percent reads in peaks value for filtering
#' @param blacklist_ratio_max Max blacklist ratio for filtering
#' @param nucleosome_signal_max Max nucleosome signal for filtering
#' @param TSS.enrichment_min Min TSS enrichment for filtering
#' @param peakwidths_max Max peak width for finding combined peaks
#' @param peakwidths_min Min peak width for finding combined peaks
#' @param passed_filters_value Min value for filtering cells based on
#' passed_filters column
#' @return A list of filtered Seurat objects
#' @export
CreateATACObjectsFilter <-
  function(data_dirs,  treatment = NULL, filter = TRUE, interactive = FALSE,
           genome = c("mm10", "hg38"), object_names = NULL,
           peak_region_fragments_min = 3000, peak_region_fragments_max = 100000,
           pct_reads_in_peaks_min = 40, blacklist_ratio_max = 0.025,
           nucleosome_signal_max = 4, TSS.enrichment_min = 2,
           peak_region_max = 3000, peakwidths_max = 10000,
           peakwidths_min = 20, passed_filters_value = 500) {

    genome <- match.arg(genome)

    if (filter == FALSE & interactive == TRUE) {
      stop("Error: Set filter=TRUE to use interactiver, otherwise set
  interactive=FALSE to skip filtering")
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
      gr <- GenomicRanges::makeGRangesFromDataFrame(peak_data)
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
    main.chroms <- GenomeInfoDb::standardChromosomes(bsgenome_obj)
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

    message('--- Generating unfiltered ATAC QC plots ---')
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
    if (is.null(treatment) == FALSE) {
      message('--- Adding Treatment metadata column ---')
      seurat_objects <- setNames(lapply(seq_along(seurat_objects), function(i) {
        seurat_obj <- seurat_objects[[i]]
        seurat_obj[["Treatment"]] <- treatment[i]
        return(seurat_obj)
      }), names(seurat_objects))
    }

    if (filter == TRUE & interactive == FALSE) {
      message('--- Filtering ATAC objects (non-interactive) ---')
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

      message('--- Generating filtered ATAC QC plots ---')
      obj <- merge(subsetted_objs[[1]], subsetted_objs[-1])

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

      return(subsetted_objs)

    } else if (filter == TRUE & interactive == TRUE) {
      message('--- Filtering ATAC objects (interactive) ---')
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


      message('--- Generating filtered ATAC QC plots ---')
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

      message('--- Saving filtered Seurat objects ---')
      saveRDS(seurat_objects, 'seurat_objects_filtered.rds')

      return(seurat_objects)
    } else {
      return(seurat_objects)
    }
  }
