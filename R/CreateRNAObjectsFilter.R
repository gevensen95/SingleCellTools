#' Create and Filter Seurat RNA Objects
#'
#' This function creates multiple Seurat objects and filters them. It takes a
#' list of directories as input. In each directory, there should be the contents
#' of the filtered feature matrix folder, a filtered feature matrix folder, or
#' a .h5 file. You can directly give the output folder from cellranger into this
#' function. It will preferentially choose the filtered_feature_matrix folder.
#' The objects are then automatically filtered based on the 15th
#' and 85th percentile for nFeature_RNA and percent.mt. This can function can
#' be run interactively to choose cutoffs without having to re-run it.
#'
#' @param data_dirs Path to directories containing matrix.mtx, features.tsv, and
#'  barcodes.tsv or .h5 files.
#' @param cells Features must be expressed in at least this many cells
#' @param features Cells must have at least this many features
#' @param mt_pattern Pattern for calculating percent mitochondrial reads
#'   (\code{percent.mt}), which this function's own quantile/threshold
#'   filtering thresholds on directly. Default \code{"^mt-"} (mouse gene
#'   symbol convention, e.g. \code{"mt-Nd1"}); pass \code{"^MT-"} for human
#'   data.
#' @param rb_pattern Pattern for calculating percent ribosomal-protein reads
#'   (\code{percent.rb}), computed alongside \code{percent.mt} for reference
#'   (informational only -- not used by this function's own filtering
#'   logic, which only ever thresholds on \code{percent.mt}). Default
#'   \code{"^(Rp[sl]|RP[SL])"} matches both mouse and human gene symbol
#'   conventions.
#' @param hb_pattern Pattern for calculating percent hemoglobin reads
#'   (\code{percent.hb}), computed alongside \code{percent.mt} for reference
#'   (informational only, like \code{rb_pattern} above). Default
#'   \code{"^(Hb[^p]|HB[^P])"} excludes the unrelated \code{"Hbp1"}/
#'   \code{"HBP1"} gene, a well-known false positive for naive
#'   \code{"^Hb"} patterns.
#' @param treatment Treatment metadata column (e.g., Age, chemical, etc.)
#' @param use_quantile Use quantile filtering method for nFeature_RNA and
#' percent.mt
#' @param quantile_value_min Minimum percentile for filtering (max will be 1-min)
#' @param feature_min Minimum cutoff for features
#' @param feature_max Maximum cutoff for features
#' @param percent_mt_max Maximum cutoff for percent.mt
#' @param interactive Run filtering interactively to pick cutoffs. A figure will
#' be shown of the the features and the percent.mt, or you can open the
#' automatically saved file
#' @param object_names Names for the Seurat objects
#' @return A list of filtered Seurat objects
#' @export

CreateRNAObjectsFilter <-
  function(data_dirs, cells = 3, features = 200,
           mt_pattern = '^mt-',
           rb_pattern = '^(Rp[sl]|RP[SL])',
           hb_pattern = '^(Hb[^p]|HB[^P])',
           treatment = NULL,
           use_quantile = TRUE, quantile_value_min = 0.15,
           feature_min = NA, feature_max = NA,
           percent_mt_max = NA, interactive = FALSE,
           object_names = NULL) {
    # Ensure thresholds are specified if not using quantiles
    if (!use_quantile) {
      if (!is.numeric(feature_min)) stop("Error: Did not specify threshold for feature_min")
      if (!is.numeric(feature_max)) stop("Error: Did not specify threshold for feature_max")
      if (!is.numeric(percent_mt_max)) stop("Error: Did not specify threshold for percent_mt_max")
    }

    message(sprintf('--- Reading data and creating Seurat objects (%d directories) ---',
                    length(data_dirs)))
    # Use lapply to read the data and create Seurat objects

    seurat_objects <- lapply(data_dirs, function(dir) {
      if (rlang::is_empty(
        list.files(dir, 'barcodes.tsv.gz|features.tsv.gz|matrix.mtx.gz')) == FALSE) {
        # Read 10X data
        seurat_data <- Seurat::Read10X(data.dir = dir)
        Seurat::CreateSeuratObject(counts = seurat_data,
                                   min.cells = cells,
                                   min.features = features,
                                   project = basename(dir))
      } else if (rlang::is_empty(
        list.files(paste(dir,'filtered_feature_bc_matrix', sep = '/'),
                   'barcodes.tsv.gz|features.tsv.gz|matrix.mtx.gz')) == FALSE) {

        seurat_data <- Seurat::Read10X(data.dir =
                                         paste(dir,'filtered_feature_bc_matrix',
                                               sep = '/'))

        # Create Seurat object
        Seurat::CreateSeuratObject(counts = seurat_data,
                                   min.cells = cells,
                                   min.features = features,
                                   project = basename(dir))
      } else if (rlang::is_empty(
        list.files(paste(paste(dir, 'outs', sep = '/'),
                         'filtered_feature_bc_matrix', sep = '/'),
                   'barcodes.tsv.gz|features.tsv.gz|matrix.mtx.gz')) == FALSE) {

        seurat_data <- Seurat::Read10X(data.dir =
                                         paste(paste(dir, 'outs', sep = '/'),
                                               'filtered_feature_bc_matrix', sep = '/'))

        # Create Seurat object
        Seurat::CreateSeuratObject(counts = seurat_data,
                                   min.cells = cells,
                                   min.features = features,
                                   project = basename(dir))
      } else if (length(list.files(dir, pattern = '\\.h5$')) > 0) {
        # Anchored on ".h5" (excludes AnnData ".h5ad" files, which end in
        # "ad" not "h5" -- the old unanchored grepl(".h5", x) treated
        # "sample.h5ad" as a match). Prefer a file with "filtered" in its
        # name if one exists (CellRanger's filtered_feature_bc_matrix.h5
        # convention); otherwise fall back to whatever single .h5 file is
        # present. Same "exactly one candidate, else error" discipline as
        # CreateRNAObjects.R's .read_10x_triplet() -- silently picking one
        # of several risks reading the wrong sample's data.
        h5_files      <- list.files(dir, pattern = '\\.h5$')
        filtered_h5   <- grep('filtered', h5_files, value = TRUE, ignore.case = TRUE)
        h5_candidates <- if (length(filtered_h5) > 0) filtered_h5 else h5_files
        if (length(h5_candidates) > 1) {
          stop("Found more than one candidate .h5 file in '", dir, "': ",
              paste(h5_candidates, collapse = ", "),
              ". Expected exactly one -- please split these into separate ",
              "sample directories.")
        }
        seurat_data <- Seurat::Read10X_h5(file.path(dir, h5_candidates))

        # Create Seurat object
        Seurat::CreateSeuratObject(counts = seurat_data,
                                   min.cells = cells,
                                   min.features = features,
                                   project = dirname(dir))
      } else {
        stop("Could not find a barcodes/features/matrix triplet or a .h5 ",
            "file in '", dir, "' (or its filtered_feature_bc_matrix ",
            "subdirectories).")
      }

    })
    # Name the list elements with the base names of the directories
    if (is.null(object_names) == TRUE) {
      names(seurat_objects) <- basename(data_dirs)
    } else {names(seurat_objects) <- object_names}

    # Add a Treatment metadata column (e.g., Age, chemical, etc.), matching
    # each data_dirs[i] to treatment[i] positionally.
    if (!is.null(treatment)) {
      if (length(treatment) != length(data_dirs)) {
        stop(sprintf(
          "`treatment` has length %d but there are %d `data_dirs` -- these must match one-to-one, or samples would silently get wrong/NA treatment labels via recycling.",
          length(treatment), length(data_dirs)))
      }
      message('--- Adding Treatment metadata column ---')
      seurat_objects <- setNames(lapply(seq_along(seurat_objects), function(i) {
        seurat_obj <- seurat_objects[[i]]
        seurat_obj[["Treatment"]] <- treatment[i]
        return(seurat_obj)
      }), names(seurat_objects))
    }

    message('--- Calculating percent mitochondrial / ribosomal / hemoglobin reads ---')
    # Add percent mitochondrial DNA to each Seurat object. percent.rb/percent.hb
    # are informational only here -- this function's own filtering logic below
    # only ever thresholds on percent.mt, matching its existing behavior;
    # they're computed so GenerateQCReport()/QCComparePlots()/CellSuiteSummary()
    # (which already treat percent.rb/percent.hb as standard metrics) can
    # report on them.
    seurat_objects <- lapply(seurat_objects, function(obj) {
      obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = mt_pattern)
      obj[["percent.rb"]] <- Seurat::PercentageFeatureSet(obj, pattern = rb_pattern)
      obj[["percent.hb"]] <- Seurat::PercentageFeatureSet(obj, pattern = hb_pattern)
      return(obj)
    })

    message('--- Saving unfiltered Seurat objects ---')
    saveRDS(seurat_objects, 'seurat_objects_unfiltered.rds')

    message('--- Generating unfiltered QC plots ---')
    # Row-bind just the metadata rather than merge()-ing the objects --
    # merge() would also combine the count matrices, which these QC plots
    # never touch (matching the fix already applied to CreateRNAObjects.R).
    meta <- dplyr::bind_rows(lapply(seurat_objects, function(x) x@meta.data))
    orig.ident <- nFeature_RNA <- percent.mt <- NULL  # silence R CMD check NSE notes

    # Create plots
    gene.plot <- ggplot2::ggplot(meta, ggplot2::aes(orig.ident, nFeature_RNA)) +
      ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered') + Ol_Reliable()
    mt.plot <- ggplot2::ggplot(meta, ggplot2::aes(orig.ident, percent.mt)) +
      ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered') + Ol_Reliable()
    print(gene.plot + mt.plot)

    message('--- Filtering cells ---')
    # Interactive mode
    if (interactive) {
      message('  Mode: interactive')
      # Ask user for thresholds interactively
      for (param in c("min nFeature_RNA", "max nFeature_RNA", "max percent.mt")) {
        use_quantile <- readline(prompt = paste("Do you want to use quantile for subsetting", param, "? (yes/no): "))
        if (tolower(use_quantile) == "yes") {
          threshold <- as.numeric(readline(prompt = paste("Enter quantile threshold for", param, " (0 to 1): ")))
          seurat_objects <- lapply(seurat_objects, function(obj) {
            if (param == "min nFeature_RNA" ) {
              subset(obj, subset = nFeature_RNA > quantile(obj$nFeature_RNA, threshold))
            } else if (param == "max nFeature_RNA") {
              subset(obj, subset = nFeature_RNA < quantile(obj$nFeature_RNA, threshold))
            } else if (param == "max percent.mt") {
              subset(obj, subset = percent.mt < quantile(obj$percent.mt, threshold))
            }
          })
        } else if (tolower(use_quantile) == "no") {
          threshold <- as.numeric(readline(prompt = paste("Enter threshold for", param, ": ")))
          seurat_objects <- lapply(seurat_objects, function(obj) {
            if (param == "min nFeature_RNA") {
              subset(obj, subset = nFeature_RNA > threshold)
            } else if (param == "max nFeature_RNA") {
              subset(obj, subset = nFeature_RNA < threshold)
            } else if (param == "max percent.mt") {
              subset(obj, subset = percent.mt < threshold)
            }
          })
        }
      }

      message('--- Generating filtered QC plots ---')
      meta <- dplyr::bind_rows(lapply(seurat_objects, function(x) x@meta.data))
      gene.plot <- ggplot2::ggplot(meta, ggplot2::aes(orig.ident, nFeature_RNA)) +
        ggplot2::geom_boxplot() + ggplot2::labs(title = 'Filtered') + Ol_Reliable()
      mt.plot <- ggplot2::ggplot(meta, ggplot2::aes(orig.ident, percent.mt)) +
        ggplot2::geom_boxplot() + ggplot2::labs(title = 'Filtered') + Ol_Reliable()
      print(gene.plot + mt.plot)

      message('--- Saving filtered Seurat objects ---')
      saveRDS(seurat_objects, 'seurat_objects_filtered.rds')

      return(seurat_objects)

    } else {
      message(sprintf('  Mode: non-interactive (use_quantile = %s)', use_quantile))
      # Non-interactive mode
      subsetted_objs <- lapply(seurat_objects, function(obj) {
        if (use_quantile) {
          feature_threshold_min <- quantile(obj$nFeature_RNA, quantile_value_min)
          feature_threshold_max <- quantile(obj$nFeature_RNA, 1 - quantile_value_min)
          mt_threshold_max <- quantile(obj$percent.mt, 1 - quantile_value_min)

          subset(obj, subset = nFeature_RNA > feature_threshold_min &
                   nFeature_RNA < feature_threshold_max &
                   percent.mt < mt_threshold_max)
        } else {
          subset(obj, subset = nFeature_RNA > feature_min &
                   nFeature_RNA < feature_max &
                   percent.mt < percent_mt_max)
        }
      })

      # Name the subsetted list with original names for clarity
      names(subsetted_objs) <- names(seurat_objects)

      message('--- Generating filtered QC plots ---')
      meta <- dplyr::bind_rows(lapply(subsetted_objs, function(x) x@meta.data))
      gene.plot <- ggplot2::ggplot(meta, ggplot2::aes(orig.ident, nFeature_RNA)) +
        ggplot2::geom_boxplot() + ggplot2::labs(title = 'Filtered') + Ol_Reliable()
      mt.plot <- ggplot2::ggplot(meta, ggplot2::aes(orig.ident, percent.mt)) +
        ggplot2::geom_boxplot() + ggplot2::labs(title = 'Filtered') + Ol_Reliable()
      print(gene.plot + mt.plot)

      message('--- Saving filtered Seurat objects ---')
      saveRDS(subsetted_objs, 'seurat_objects_filtered.rds')

      return(subsetted_objs)
    }
  }
