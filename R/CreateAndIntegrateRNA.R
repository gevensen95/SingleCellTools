#' Create and Integrate Seurat RNA Objects
#'
#' This function creates multiple Seurat objects and then integrates them. It
#' takes a list of directories as input. In each directory, there should be
#' the contents of the filtered feature matrix folder, a filtered feature
#' matrix folder, or a .h5 file. You can directly give the output folder from
#' cellranger into this function. It will preferentially choose the
#' filtered_feature_matrix folder. The function then automatically filters each
#' objects and integrates using Harmony. The user can change these parameters.
#' This function also save multiple Seurat objects along the way and
#' runs FindAllMarkers on the clusters.
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
#' @param cell_IDs Adds cell identities to the beginning of each barcode
#' @param to_regress Variable to regress out in normalization
#' @param cluster_resolution Resolution for cell clustering
#' @param max_dims Maximum dimensions for FindNeighbors and UMAP
#' @param use_SCT Use SCTransform method or LogNormalize (FALSE)
#' @param sct_assay Assay for SCTransform normalization
#' @param save_rds_file Save final Seurat object
#' @param file_name Name of files
#' @param use_elbow_plot Use the ElbowPlot to determine number of componenets
#' for FindNeighbors and UMAP
#' @param spatial If your data is a spatial (i.e., Visium)
#' @param integration Method for integrating data, see IntegrateLayers
#' @param integration_normalization Normalization method used
#' @param integration_assay Assay that has normalized data
#' @param integration_reduction Reduction to use for integration
#' @param new_reduction Name of new integrated reduction
#' @param k_anchor How many neighbors (k) to use when picking anchors
#' @param k_weight Number of neighbors to consider when weighting anchors
#' @return An integrated Seurat object
#' @param markers Find all markers
#' @param workers Number of parallel workers to use (via \code{future.apply})
#'   for reading/creating each sample's Seurat object -- fully independent
#'   across samples. Defaults to \code{length(data_dirs)} (one worker per
#'   sample); errors up front if that (or an explicit value) exceeds
#'   \code{parallel::detectCores()}, naming the number of cores actually
#'   available. Pass \code{workers = 1} to run sequentially instead.
#'   \code{workers > 1} spins up that many background R sessions via
#'   \code{future::plan(multisession)}, restored on exit. Note each worker
#'   holds its own copy of that sample's data, so peak memory scales with
#'   \code{workers}.
#' @export
CreateAndIntegrateRNA <-
  function(data_dirs, cells = 3, features = 200,
           mt_pattern = '^mt-',
           rb_pattern = '^(Rp[sl]|RP[SL])',
           hb_pattern = '^(Hb[^p]|HB[^P])',
           treatment = NULL,
           use_quantile = TRUE, quantile_value_min = 0.15,
           feature_min = NA, feature_max = NA,
           percent_mt_max = NA, interactive = FALSE,
           object_names = NULL,
           # Relies on R's lazy default-argument evaluation: this default
           # isn't resolved until cell_IDs is first used (at the real merge
           # below, well after `seurat_objects` has been built/filtered in
           # this function's own frame), so it correctly picks up the final
           # sample names -- but only because nothing forces `cell_IDs`
           # earlier in the body. Pass an explicit vector if that ordering
           # ever changes.
           cell_IDs = names(seurat_objects), to_regress = 'percent.mt',
           cluster_resolution = 0.3, max_dims = 15, use_SCT = TRUE,
           sct_assay = 'RNA',
           save_rds_file = TRUE, file_name = NULL,
           use_elbow_plot = FALSE, spatial = FALSE,
           integration = 'HarmonyIntegration',
           integration_normalization = 'SCT', integration_assay = 'SCT',
           integration_reduction = 'pca', new_reduction = 'harmony',
           k_anchor = NULL, k_weight = NULL,
           markers = TRUE, workers = length(data_dirs)) {
    workers <- .resolve_workers(workers, n_samples = length(data_dirs),
                                was_default = missing(workers))
    # Ensure thresholds are specified if not using quantiles
    if (!use_quantile) {
      if (!is.numeric(feature_min)) stop("Error: Did not specify threshold for feature_min")
      if (!is.numeric(feature_max)) stop("Error: Did not specify threshold for feature_max")
      if (!is.numeric(percent_mt_max)) stop("Error: Did not specify threshold for percent_mt_max")
    }

    if (use_quantile) {
      if (is.numeric(feature_min)) stop("Error: Set quantile=TRUE and specificed hard cut off for min nFeature_RNA Pick only one.")
      if (is.numeric(feature_max)) stop("Error: Set quantile=TRUE and specificed hard cut off for max nFeature_RNA Pick only one.")
      if (is.numeric(percent_mt_max)) stop("Error: Set quantile=TRUE and specificed hard cut off for percent.mt. Pick only one.")
    }

    if (workers > 1) {
      if (!requireNamespace("future.apply", quietly = TRUE)) {
        stop("Package 'future.apply' is required for workers > 1. ",
            "install.packages('future.apply')")
      }
      old_plan <- future::plan(future::multisession, workers = workers)
      on.exit(future::plan(old_plan), add = TRUE)
    }

    message(sprintf('--- Reading data and creating Seurat objects (%d directories)%s ---',
                    length(data_dirs),
                    if (workers > 1) sprintf(', %d parallel workers', workers) else ''))
    # Reading + CreateSeuratObject is fully independent per directory, so
    # this runs in parallel when workers > 1; otherwise a plain sequential
    # lapply, unchanged from before.

    .read_one <- function(dir) {
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

    }

    seurat_objects <- if (workers > 1) {
      future.apply::future_lapply(data_dirs, .read_one, future.seed = TRUE)
    } else {
      lapply(data_dirs, .read_one)
    }
    # Name the list elements with the base names of the directories
    if (is.null(object_names) == TRUE) {
      names(seurat_objects) <- basename(data_dirs)
    } else {names(seurat_objects) <- object_names}

    message('--- Calculating percent mitochondrial / ribosomal / hemoglobin reads ---')
    # Add percent mitochondrial DNA to each Seurat object. percent.rb/percent.hb
    # are informational only here -- this function's own filtering logic below
    # (interactive and non-interactive blocks) only ever thresholds on
    # percent.mt, matching its existing behavior; they're computed so
    # GenerateQCReport()/QCComparePlots()/CellSuiteSummary() (which already
    # treat percent.rb/percent.hb as standard metrics) can report on them.
    seurat_objects <- lapply(seurat_objects, function(obj) {
      obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = mt_pattern)
      obj[["percent.rb"]] <- Seurat::PercentageFeatureSet(obj, pattern = rb_pattern)
      obj[["percent.hb"]] <- Seurat::PercentageFeatureSet(obj, pattern = hb_pattern)
      return(obj)
    })

    # Add a to_regress to metadata to specify treatment
    if (is.null(treatment) == FALSE){
      message('--- Adding Treatment metadata column ---')
      seurat_objects <- setNames(lapply(seq_along(seurat_objects), function(i) {
        seurat_obj <- seurat_objects[[i]]
        seurat_obj[["Treatment"]] <- treatment[i]
        return(seurat_obj)
      }), names(seurat_objects))
    }

    message('--- Saving unfiltered Seurat objects ---')
    saveRDS(seurat_objects, file = 'seurat_objects_unfiltered.rds')

    message('--- Generating unfiltered QC plots ---')
    # Row-bind just the metadata rather than merge()-ing the objects here --
    # this merge only ever fed two boxplots and was thrown away immediately
    # after (the real merge happens further down, right before
    # normalization, where the full counts are actually needed).
    qc_meta <- dplyr::bind_rows(lapply(seurat_objects, function(x) x@meta.data))
    orig.ident <- nFeature_RNA <- percent.mt <- NULL  # silence R CMD check NSE notes

    # Create plots
    gene.plot <- ggplot2::ggplot(qc_meta, ggplot2::aes(orig.ident, nFeature_RNA)) +
      ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered') + Ol_Reliable()
    mt.plot <- ggplot2::ggplot(qc_meta, ggplot2::aes(orig.ident, percent.mt)) +
      ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered') + Ol_Reliable()
    print(gene.plot + mt.plot)
    ggplot2::ggsave('unfiltered_features_percentMT.pdf', height = 5, width = 7)

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

      meta <- dplyr::bind_rows(lapply(seurat_objects, function(x) x@meta.data))
      gene.plot <- ggplot2::ggplot(meta, ggplot2::aes(orig.ident, nFeature_RNA)) +
        ggplot2::geom_boxplot() + ggplot2::labs(title = 'Filtered') + Ol_Reliable()
      mt.plot <- ggplot2::ggplot(meta, ggplot2::aes(orig.ident, percent.mt)) +
        ggplot2::geom_boxplot() + ggplot2::labs(title = 'Filtered') + Ol_Reliable()
      print(gene.plot + mt.plot)

    } else {
      message(sprintf('  Mode: non-interactive (use_quantile = %s)', use_quantile))
      # Non-interactive mode
      seurat_objects <- lapply(seurat_objects, function(obj) {
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

      meta <- dplyr::bind_rows(lapply(seurat_objects, function(x) x@meta.data))
      gene.plot <- ggplot2::ggplot(meta, ggplot2::aes(orig.ident, nFeature_RNA)) +
        ggplot2::geom_boxplot() + ggplot2::labs(title = 'Filtered') + Ol_Reliable()
      mt.plot <- ggplot2::ggplot(meta, ggplot2::aes(orig.ident, percent.mt)) +
        ggplot2::geom_boxplot() + ggplot2::labs(title = 'Filtered') + Ol_Reliable()
      print(gene.plot + mt.plot)
    }

    message('--- Saving filtered Seurat objects ---')
    saveRDS(seurat_objects, file = 'seurat_objects_filtered.rds')

    if(integration != 'HarmonyIntegration' & new_reduction == 'harmony') {
      stop('\n\n  Error: Integration method is not the default (HarmonyIntegration).\nChange new_reduction to match integration method')
    }

    # `&` binds tighter than `|` in R, so the old un-parenthesized version of
    # these checks (`A | B | C & D`) parsed as `A | B | (C & D)` -- the
    # is.null() check only ever attached to the JointPCAIntegration clause,
    # so RPCAIntegration/CCAIntegration unconditionally hit stop() below
    # regardless of whether k_anchor/k_weight were actually supplied.
    # Explicit %in% + a single & fixes both.
    if (integration %in% c('RPCAIntegration', 'CCAIntegration', 'JointPCAIntegration') &
        is.null(k_anchor)) {
      stop('\n\n  Error: Integration method is RPCAIntegration/CCAIntegration/JointPCAIntegration.\nSpecifiy k_anchor to an integar value (recommend 20)')
    }

    if (integration %in% c('RPCAIntegration', 'CCAIntegration', 'JointPCAIntegration') &
        is.null(k_weight)) {
      stop('\n\n  Error: Integration method is RPCAIntegration/CCAIntegration/JointPCAIntegration.\nSpecifiy k_weight to an integar value (recommend 100)')
    }

    message('--- Merging Seurat objects ---')
    if (spatial){
      obj <- suppressWarnings(merge(seurat_objects[[1]], seurat_objects[-1],
                                    add.cell.ids = cell_IDs))
    } else {
      obj <- merge(seurat_objects[[1]], seurat_objects[-1],
                   add.cell.ids = cell_IDs)
    }

    message('--- Normalizing data ---')
    if (use_SCT){
      message('  Running SCTransform')
      calculate_median <- function(data, column_name) {
        data |>
          dplyr::group_by(orig.ident) |>
          dplyr::summarise(Median = median(.data[[column_name]], na.rm = TRUE)) |>
          dplyr::arrange(Median)
      }
      med_counts <- calculate_median(
        obj@meta.data,
        colnames(obj@meta.data)[stringr::str_detect(colnames(obj@meta.data), 'nCount')][1]
      )

      obj <- Seurat::SCTransform(obj, vars.to.regress = to_regress,
                                 assay = sct_assay, scale_factor = med_counts$Median[1])
    }
    else{
      message('  Running NormalizeData / FindVariableFeatures / ScaleData')
      obj <- Seurat::NormalizeData(obj)
      obj <- Seurat::FindVariableFeatures(obj)
      obj <- Seurat::ScaleData(obj, vars.to.regress = to_regress)
    }

    message('--- Running PCA ---')
    obj <- Seurat::RunPCA(obj)
    if (use_elbow_plot) {

      elbow_plot <- Seurat::ElbowPlot(obj)
      print(elbow_plot)

      pcs <- as.numeric(readline(prompt = 'Enter # of PCs: '))

      message(sprintf('--- Integrating layers (method: %s) ---', integration))
      obj <- Seurat::IntegrateLayers(object = obj, method = integration,
                             orig.reduction = integration_reduction,
                             assay = integration_assay,
                             normalization.method = integration_normalization,
                             new.reduction = new_reduction,
                             k.anchor = k_anchor, k.weight = k_weight)
      message('--- Finding neighbors and clusters ---')
      obj <- Seurat::FindNeighbors(obj, reduction = new_reduction, dims = 1:pcs)
      obj <- Seurat::FindClusters(obj, resolution = cluster_resolution)
      message('--- Running UMAP ---')
      obj <- Seurat::RunUMAP(obj, reduction = new_reduction, dims = 1:pcs)

      if (save_rds_file == TRUE & is.null(file_name) == TRUE) {
        message('--- Saving integrated Seurat object to RDS ---')
        saveRDS(obj, paste(new_reduction, 'merged_seurat_objects.rds',
                           sep = '_'))
      } else if (save_rds_file == TRUE & is.null(file_name) == FALSE) {
        message('--- Saving integrated Seurat object to RDS ---')
        saveRDS(obj, paste(file_name, 'merged_seurat_objects.rds',
                           sep = '_'))
      }

    } else {
      message(sprintf('--- Integrating layers (method: %s) ---', integration))
      obj <- Seurat::IntegrateLayers(object = obj, method = integration,
                             orig.reduction = integration_reduction,
                             assay = integration_assay,
                             normalization.method = integration_normalization,
                             new.reduction = new_reduction,
                             k.anchor = k_anchor, k.weight = k_weight)
      message('--- Finding neighbors and clusters ---')
      obj <- Seurat::FindNeighbors(obj, reduction = new_reduction, dims = 1:max_dims)
      obj <- Seurat::FindClusters(obj, resolution = cluster_resolution)
      message('--- Running UMAP ---')
      obj <- Seurat::RunUMAP(obj, reduction = new_reduction, dims = 1:max_dims)
      if (use_SCT == TRUE & markers == TRUE){
        if (save_rds_file == TRUE & is.null(file_name) == TRUE) {
          message('--- Saving integrated Seurat object to RDS ---')
          saveRDS(obj, paste(new_reduction, 'merged_seurat_objects.rds',
                             sep = '_'))
        } else if (save_rds_file == TRUE & is.null(file_name) == FALSE) {
          message('--- Saving integrated Seurat object to RDS ---')
          saveRDS(obj, paste(file_name, 'merged_seurat_objects.rds',
                             sep = '_'))
        }
        else {
          return(obj)
        }
      }
    }

    message('--- Saving cluster DimPlot ---')
    Seurat::DimPlot(obj, label = T)
    ggplot2::ggsave('dimplot_seurat_clusters.pdf', height = 5, width = 7)

    if (markers == TRUE){
      message('--- Running FindAllMarkers ---')
      markers <- Seurat::FindAllMarkers(obj, logfc.threshold = 1, only.pos = TRUE,
                                        min.pct = 0.25)
      write.csv(markers, 'markers_all.csv')
      return(obj)
    } else {return(obj)}
  }
