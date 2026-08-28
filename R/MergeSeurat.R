#' Merge and Integrate Seurat Objects
#'
#' Merges a list of Seurat objects, normalizes (SCT or LogNormalize), runs PCA,
#' integrates layers, finds neighbors/clusters, runs UMAP, and (optionally)
#' computes per-cluster markers.
#'
#' \strong{ATAC (Signac) input} is detected automatically -- checked via
#' whether the first object's default assay is a \code{ChromatinAssay}; every
#' object in \code{seurat_objects} must agree (mixing ATAC and non-ATAC
#' objects in one call errors). When detected, \code{use_SCT}/normalization
#' and \code{RunPCA()} are replaced entirely by \code{RunATACWrapper()}
#' (TF-IDF -> top features -> LSI, Signac's standard ATAC recipe -- peak
#' accessibility counts aren't gene-expression counts, so neither
#' log-normalization/SCT nor PCA are the right tool), producing an
#' \code{"lsi"} reduction instead of \code{"pca"}. \code{integration_reduction}
#' defaults to \code{"lsi"} automatically in this case (still overridable).
#' By LSI convention the first component usually correlates with sequencing
#' depth rather than biology, so downstream \code{FindNeighbors()}/
#' \code{RunUMAP()} use \code{dims = atac_lsi_first_dim:max_dims} (default
#' \code{2:max_dims}) rather than starting at 1. Only \code{spatial = 'no'}
#' and \code{integration = 'HarmonyIntegration'} are supported for ATAC input
#' (errors otherwise) -- the merge-time Visium/Xenium assay coercion doesn't
#' apply to a ChromatinAssay, and anchor-based integration (RPCA/CCA/JointPCA)
#' on LSI space isn't a validated combination here, mirroring the same
#' constraint already placed on \code{banksy = TRUE}. \code{markers = TRUE}
#' (the default) is silently skipped for ATAC input rather than run with
#' RNA-tuned thresholds -- call \code{Signac::FindMarkers()}/
#' \code{FindAllMarkers()} directly with ATAC-appropriate settings instead.
#'
#' @param seurat_objects A named list of Seurat objects.
#' @param cell_IDs Character vector of cell ID prefixes (defaults to
#'   `names(seurat_objects)`).
#' @param to_regress Variables to regress out during normalization.
#' @param cluster_resolution Resolution for `FindClusters`.
#' @param max_dims Number of PCs to use when `use_elbow_plot = FALSE`.
#' @param use_SCT Logical; use SCTransform (TRUE) or LogNormalize pipeline (FALSE).
#' @param sct_assay Assay to pass to `SCTransform`.
#' @param save_rds_file Logical; save merged object to disk.
#' @param file_name Optional prefix for the saved RDS (defaults to
#'   `new_reduction`).
#' @param use_elbow_plot Logical; show elbow plot and prompt for number of PCs.
#' @param spatial One of `'no'`, `'Visium'`, `'Xenium'`.
#' @param integration Integration method (e.g. `'HarmonyIntegration'`,
#'   `'RPCAIntegration'`, `'CCAIntegration'`, `'JointPCAIntegration'`).
#'   `'HarmonyIntegration'` is always run via `harmony::RunHarmony()` called
#'   directly rather than `Seurat::IntegrateLayers()`, since Seurat's
#'   built-in `HarmonyIntegration` wrapper depends on
#'   `harmony::HarmonyMatrix()`, which was removed from the harmony package
#'   in its 1.0 rewrite (current CRAN harmony only exports `RunHarmony()`).
#'   Every other method still goes through `IntegrateLayers()` unchanged.
#' @param integration_normalization,integration_assay,integration_reduction
#'   Arguments passed through to `IntegrateLayers`.
#' @param new_reduction Name of the reduction created by integration.
#' @param k_anchor,k_weight Arguments required for RPCA/CCA/JointPCA integration.
#' @param markers Logical; if TRUE (default), run
#'   \code{\link[Seurat]{FindAllMarkers}}, write the full result to
#'   \code{markers_all.csv}, and draw the top \code{marker_n} markers per
#'   cluster to \code{marker_plot.pdf} via \code{\link{TopMarkerPlot}} --
#'   so the panel matches the annotated style and auto-sizing of every other
#'   marker plot in this package. If no gene survives filtering the plot is
#'   skipped with a message; the merged object is still returned either way.
#' @param marker_n Number of top marker genes per cluster to draw when
#'   \code{markers = TRUE}. Default \code{10}. Ignored otherwise.
#' @param group_column Metadata column used to compute median nCount per group
#'   for `SCTransform`'s scale_factor.
#' @param common_genes_only Logical; if TRUE, subset all input objects to genes
#'   present in every object before merging.
#' @param common_genes_assay Assay to inspect for gene names when
#'   `common_genes_only = TRUE`.
#' @param banksy Logical; if TRUE, run BANKSY spatial-aware clustering
#'   (via `RunBanksyWrapper()`) on the merged, normalized object instead of
#'   plain PCA. Requires `spatial` to be `'Visium'` or `'Xenium'` (BANKSY
#'   needs spatial coordinates) and `integration` to be
#'   `'HarmonyIntegration'` -- RunBanksy() produces a single consolidated
#'   (non-split-layer) assay, which `Seurat::IntegrateLayers()`'s
#'   anchor-based methods (RPCA/CCA/JointPCA) can't work with. Since
#'   `'HarmonyIntegration'` already bypasses `IntegrateLayers()` entirely
#'   (see `integration`), this isn't a special case any more -- Harmony
#'   runs directly on the `pca_banksy` embedding the same way it runs on
#'   `pca` otherwise (same as the BANKSY-Seurat vignette's own Harmony
#'   workflow). Default FALSE.
#' @param banksy_lambda Numeric in `[0,1]`; BANKSY's spatial weight parameter.
#'   Low values (~0.2, default) favor cell-typing, high values (~0.8) favor
#'   spatial domain segmentation. Only used when `banksy = TRUE`.
#' @param banksy_k_geom Number of spatial neighbors for BANKSY. Only used
#'   when `banksy = TRUE`.
#' @param banksy_assay Assay to compute BANKSY on. Defaults to `'SCT'` when
#'   `use_SCT = TRUE`, else the technology-native assay actually normalized
#'   above (`'Spatial'` or `'Xenium'`, matching `spatial`). Only used when
#'   `banksy = TRUE`.
#' @param atac_lsi_first_dim First LSI component to include in
#'   `FindNeighbors()`/`RunUMAP()`'s `dims` when ATAC input is detected.
#'   Default `2`, dropping the first component -- by LSI convention it
#'   usually correlates with sequencing depth rather than biology (see
#'   Signac's own vignette). Pass `1` to include it. Ignored for non-ATAC
#'   input.
#' @return A merged, integrated, clustered Seurat object.
#' @export
MergeSeurat <- function(seurat_objects,
                        cell_IDs = names(seurat_objects),
                        to_regress = 'percent.mt',
                        cluster_resolution = 0.3,
                        max_dims = 15,
                        use_SCT = TRUE,
                        sct_assay = 'RNA',
                        save_rds_file = TRUE,
                        file_name = NULL,
                        use_elbow_plot = FALSE,
                        spatial = 'no',
                        integration = 'HarmonyIntegration',
                        integration_normalization = 'SCT',
                        integration_assay = 'SCT',
                        integration_reduction = 'pca',
                        new_reduction = 'harmony',
                        k_anchor = NULL,
                        k_weight = NULL,
                        markers = TRUE,
                        marker_n = 10,
                        group_column = 'orig.ident',
                        common_genes_only = FALSE,
                        common_genes_assay = NULL,
                        banksy = FALSE,
                        banksy_lambda = 0.2,
                        banksy_k_geom = 15,
                        banksy_assay = NULL,
                        atac_lsi_first_dim = 2) {

  # ---- Detect ATAC (Signac ChromatinAssay) input --------------------------
  # Every object in one call must be the same modality -- checked up front so
  # a mixed list fails here with an actionable message rather than partway
  # through NormalizeData()/RunPCA() on a ChromatinAssay object it was never
  # designed for. Tolerant of non-Seurat elements (e.g. the plain numeric
  # placeholders some argument-validation tests pass in to exercise the
  # `stop()`s below without building real Seurat objects) -- those are
  # treated as non-ATAC rather than erroring out here.
  .is_atac_obj <- function(o) {
    if (!inherits(o, 'Seurat')) return(FALSE)
    a <- tryCatch(o[[Seurat::DefaultAssay(o)]], error = function(e) NULL)
    inherits(a, 'ChromatinAssay')
  }
  is_atac <- .is_atac_obj(seurat_objects[[1]])
  if (length(seurat_objects) > 1) {
    other_is_atac <- vapply(seurat_objects[-1], .is_atac_obj, logical(1))
    if (any(other_is_atac != is_atac)) {
      stop('\n\n  Error: `seurat_objects` mixes ATAC (ChromatinAssay) and ',
           'non-ATAC objects in one call -- MergeSeurat() expects one ',
           'modality per call. Merge each modality separately.')
    }
  }
  # `integration_reduction` defaults to 'pca', which doesn't exist for ATAC
  # input (only 'lsi' does) -- switch the default when the caller didn't
  # explicitly set it themselves.
  if (is_atac && missing(integration_reduction)) integration_reduction <- 'lsi'

  # ---- Argument validation ------------------------------------------------
  if (isTRUE(is_atac) && spatial != 'no') {
    stop('\n\n  Error: ATAC (ChromatinAssay) input was detected, but ',
         'spatial = "', spatial, '". Only spatial = "no" is supported for ',
         'ATAC input -- the Visium/Xenium merge-time assay coercion doesn\'t ',
         'apply to a ChromatinAssay.')
  }
  if (isTRUE(is_atac) && integration != 'HarmonyIntegration') {
    stop('\n\n  Error: ATAC (ChromatinAssay) input was detected, but ',
         'integration = "', integration, '". Only integration = ',
         '"HarmonyIntegration" is currently supported for ATAC input -- ',
         'anchor-based integration (RPCA/CCA/JointPCA) on LSI space isn\'t ',
         'a validated combination here (the same constraint already placed ',
         'on banksy = TRUE, for the same underlying reason).')
  }
  if (integration != 'HarmonyIntegration' & new_reduction == 'harmony') {
    stop('\n\n  Error: Integration method is not the default (HarmonyIntegration).\n  Change new_reduction to match integration method')
  }
  if (integration %in% c('RPCAIntegration', 'CCAIntegration', 'JointPCAIntegration') &&
      is.null(k_anchor)) {
    stop('\n\n  Error: Integration method is ', integration,
         '.\n  Specify k_anchor to an integer value (recommend 20)')
  }
  if (integration %in% c('RPCAIntegration', 'CCAIntegration', 'JointPCAIntegration') &&
      is.null(k_weight)) {
    stop('\n\n  Error: Integration method is ', integration,
         '.\n  Specify k_weight to an integer value (recommend 100)')
  }
  if (isTRUE(banksy) && !spatial %in% c('Visium', 'Xenium')) {
    stop('\n\n  Error: banksy = TRUE requires spatial = "Visium" or "Xenium"',
         ' (BANKSY needs spatial coordinates on the merged object).',
         '\n  For other spatial technologies (e.g. FOV-based Xenium/CosMx',
         ' handled outside this function), call RunBanksyWrapper() directly',
         ' on the merged object instead.')
  }
  if (isTRUE(banksy) && integration != 'HarmonyIntegration') {
    stop('\n\n  Error: banksy = TRUE currently only supports',
         ' integration = "HarmonyIntegration".',
         '\n  RunBanksy() builds a single consolidated (non-split-layer)',
         ' assay, which Seurat::IntegrateLayers()\'s anchor-based methods',
         ' (RPCA/CCA/JointPCA) require split per-sample layers to work',
         ' with. This mirrors the BANKSY-Seurat vignette\'s own workflow,',
         ' which integrates via harmony::RunHarmony() directly rather than',
         ' IntegrateLayers().')
  }

  # ---- Optional: restrict to genes common to every object before merging --
  if (isTRUE(common_genes_only)) {
    message('--- Computing common genes across objects ---')
    if (is.null(common_genes_assay)) {
      common_genes_assay <- switch(
        spatial,
        'Visium' = 'Spatial',
        'Xenium' = 'Xenium',
        Seurat::DefaultAssay(seurat_objects[[1]])
      )
    }
    gene_lists <- lapply(seurat_objects, function(o) {
      rownames(SeuratObject::LayerData(o, assay = common_genes_assay,
                                       layer = 'counts'))
    })
    per_object_n <- vapply(gene_lists, length, integer(1))
    obj_labels <- names(seurat_objects)
    if (is.null(obj_labels)) obj_labels <- paste0('obj_', seq_along(seurat_objects))
    union_genes <- Reduce(union, gene_lists)
    common      <- Reduce(intersect, gene_lists)
    if (length(common) == 0) {
      stop('\n\n  Error: No genes are common across all objects in assay "',
           common_genes_assay, '". Check that gene identifiers ',
           '(symbols vs Ensembl IDs) are harmonized across objects.')
    }
    message(sprintf('  Genes per object (assay: %s):', common_genes_assay))
    for (i in seq_along(per_object_n)) {
      message(sprintf('    %s: %d genes', obj_labels[i], per_object_n[i]))
    }
    pct_kept <- round(100 * length(common) / length(union_genes), 1)
    message(sprintf('  Total unique genes across all objects (union): %d',
                    length(union_genes)))
    message(sprintf('  Genes shared across all objects (intersection): %d (%.1f%% of union)',
                    length(common), pct_kept))
    message(sprintf('  common_genes_only = TRUE: subsetting all %d objects to the %d shared genes.',
                    length(seurat_objects), length(common)))
    seurat_objects <- lapply(seurat_objects, function(o) o[common, ])

    # After dropping non-common genes, some cells may now have zero total
    # counts. Recompute per-cell totals on the shared-gene set and drop those
    # empty cells, then clean up the temporary metadata column.
    message('--- Recomputing per-cell counts on shared genes and dropping empty cells ---')

    # Helper: robust per-cell count summation that tolerates the quirks of
    # SeuratObject::LayerData() in v5. Specifically:
    #   * If the assay still has split layers (e.g. 'counts.1', 'counts.2'
    #     left over from a previous merge), `layer = 'counts'` may not match
    #     anything and return NULL, *or* it may return a joined matrix with
    #     unexpected dim. Iterate over every counts-like layer instead.
    #   * If the result is reduced to a single row (one shared gene), the
    #     v5 method drops the matrix dimension and returns a plain numeric
    #     vector — base::colSums then errors with
    #     "'x' must be an array of at least two dimensions".
    .cell_counts_on_assay <- function(o, assay) {
      lyrs <- tryCatch(SeuratObject::Layers(o, assay = assay),
                       error = function(e) NULL)
      if (is.null(lyrs) || !length(lyrs)) lyrs <- 'counts'
      counts_lyrs <- grep('^counts(\\.|$)', lyrs, value = TRUE)
      if (!length(counts_lyrs)) counts_lyrs <- 'counts'

      cells <- colnames(o)
      totals <- setNames(numeric(length(cells)), cells)

      for (lyr in counts_lyrs) {
        m <- tryCatch(
          SeuratObject::LayerData(o, assay = assay, layer = lyr),
          error = function(e) NULL
        )
        if (is.null(m)) next

        # Coerce 1D returns (single-feature subset) back to a 1-row matrix
        if (is.null(dim(m))) {
          if (length(m) != length(cells)) next
          m <- matrix(as.numeric(m), nrow = 1,
                      dimnames = list(NULL, cells))
        }

        # Matrix::colSums dispatches correctly for both dense and sparse
        cs <- Matrix::colSums(m)
        if (is.null(names(cs))) names(cs) <- colnames(m)
        # Layers in v5 can be cell-subsets of the assay; only add overlap
        idx <- intersect(names(cs), names(totals))
        totals[idx] <- totals[idx] + cs[idx]
      }
      totals
    }

    seurat_objects <- lapply(seq_along(seurat_objects), function(i) {
      o   <- seurat_objects[[i]]
      lab <- obj_labels[i]

      cell_counts <- .cell_counts_on_assay(o, common_genes_assay)

      n_before   <- ncol(o)
      keep_cells <- names(cell_counts)[cell_counts > 0]
      n_drop     <- n_before - length(keep_cells)

      if (n_drop > 0) {
        # Use cells = ... rather than subset = Counts > 0 so Seurat reconciles
        # FOVs / images against an explicit cell list. The expression form can
        # leave stale images attached on spatial objects.
        if (spatial %in% c('Visium', 'Xenium') && exists('subset_opt', mode = 'function')) {
          o <- subset_opt(o, cells = keep_cells)
        } else {
          o <- subset(o, cells = keep_cells)
        }
      }

      message(sprintf('    %s: dropped %d / %d cells with 0 counts on shared genes (%d remaining)',
                      lab, n_drop, n_before, ncol(o)))
      o
    })
    names(seurat_objects) <- obj_labels
  }

  # ---- Merge --------------------------------------------------------------
  message('--- Merging Seurat objects ---')
  if (spatial == 'Visium') {
    obj <- suppressWarnings(merge(seurat_objects[[1]], seurat_objects[-1],
                                  add.cell.ids = cell_IDs))
    obj[["RNA"]] <- methods::as(object = obj[["Spatial"]], Class = "Assay5")
  } else if (spatial == 'Xenium') {
    obj <- suppressWarnings(merge(seurat_objects[[1]], seurat_objects[-1],
                                  add.cell.ids = cell_IDs))
    obj[["RNA"]] <- methods::as(object = obj[["Xenium"]], Class = "Assay5")
  } else if (spatial == 'no') {
    obj <- merge(seurat_objects[[1]], seurat_objects[-1],
                 add.cell.ids = cell_IDs)
  }

  if (isTRUE(is_atac)) {
    # ---- ATAC: TF-IDF -> top features -> LSI (replaces Normalize + PCA) ----
    # Peak accessibility counts aren't gene-expression counts, so neither
    # log-normalization/SCT nor PCA is the right tool here -- RunATACWrapper()
    # runs Signac's standard TF-IDF -> FindTopFeatures -> RunSVD recipe and
    # produces an "lsi" reduction in place of "pca".
    message('--- Running RunATACWrapper() (TF-IDF -> top features -> LSI) ---')
    obj <- RunATACWrapper(obj, n_components = max_dims)
  } else {
    # ---- Normalize ----------------------------------------------------------
    message('--- Normalizing data ---')
    if (use_SCT) {
      calculate_median <- function(data, group_column, column_name) {
        data %>%
          dplyr::group_by(.data[[group_column]]) %>%
          dplyr::summarise(Median = stats::median(.data[[column_name]], na.rm = TRUE)) %>%
          dplyr::arrange(Median)
      }
      nCount_col <- colnames(obj@meta.data)[stringr::str_detect(colnames(obj@meta.data),
                                                                'nCount')][1]
      med_counts <- calculate_median(data = obj@meta.data,
                                     group_column = group_column,
                                     column_name = nCount_col)
      message('  Running SCTransform')
      obj <- Seurat::SCTransform(obj, vars.to.regress = to_regress,
                                 assay = sct_assay,
                                 scale_factor = med_counts$Median[1])
    } else {
      message('  Running NormalizeData / FindVariableFeatures / ScaleData')
      obj <- Seurat::NormalizeData(obj)
      obj <- Seurat::FindVariableFeatures(obj)
      obj <- Seurat::ScaleData(obj, vars.to.regress = to_regress)
    }

    # ---- BANKSY (optional spatial-aware clustering) --------------------------
    if (isTRUE(banksy)) {
      if (is.null(banksy_assay)) {
        # SCTransform always names its output assay "SCT" regardless of
        # `sct_assay`. Without SCT, NormalizeData()/ScaleData() above ran on
        # whatever assay was DefaultAssay(obj) after merging -- for spatial
        # objects that's the technology-native assay ("Spatial"/"Xenium"),
        # *not* the "RNA" copy created during the merge step above.
        banksy_assay <- if (use_SCT) {
          'SCT'
        } else if (spatial == 'Visium') {
          'Spatial'
        } else {
          'Xenium'
        }
      }
      # merge() leaves Seurat v5 assay layers split per input object (e.g.
      # "data.1"/"data.2" rather than one unified "data" layer) until
      # JoinLayers() is called -- NormalizeData()/ScaleData() above tolerate
      # that split fine, but RunBanksy()'s GetAssayData(layer = slot) call
      # does not ("GetAssayData doesn't work for multiple layers in v5
      # assay."). Join defensively; a no-op if the assay is already joined
      # (e.g. a single pre-merged object, or spatial = 'no').
      obj <- tryCatch(SeuratObject::JoinLayers(obj, assay = banksy_assay),
                      error = function(e) obj)
      message('--- Running BANKSY (spatial-aware clustering) ---')
      obj <- RunBanksyWrapper(obj, lambda = banksy_lambda, assay = banksy_assay,
                              k_geom = banksy_k_geom, group = group_column,
                              run_pca = TRUE, npcs = max_dims)
    } else {
      # ---- PCA ----------------------------------------------------------------
      message('--- Running PCA ---')
      obj <- Seurat::RunPCA(obj)
    }
  }

  # ---- Choose dims (elbow prompt or fixed max_dims) -----------------------
  if (use_elbow_plot) {
    elbow_reduction <- if (isTRUE(is_atac)) {
      'lsi'
    } else if (isTRUE(banksy)) {
      'pca_banksy'
    } else {
      'pca'
    }
    elbow_plot <- Seurat::ElbowPlot(obj, reduction = elbow_reduction)
    print(elbow_plot)
    dims_to_use <- as.numeric(readline(prompt = 'Enter # of PCs: '))
  } else {
    dims_to_use <- max_dims
  }

  # ATAC's first LSI component usually correlates with sequencing depth
  # rather than biology (Signac convention) -- start FindNeighbors()/
  # RunUMAP() at atac_lsi_first_dim (default 2) instead of 1 in that case.
  dims_start <- if (isTRUE(is_atac)) atac_lsi_first_dim else 1

  # ---- Integrate ------------------------------------------------------------
  # HarmonyIntegration always goes through harmony::RunHarmony() directly
  # rather than Seurat::IntegrateLayers(method = "HarmonyIntegration").
  # Seurat's built-in HarmonyIntegration wrapper calls the legacy
  # harmony::HarmonyMatrix(), which was removed from the harmony package as
  # part of its 1.0 rewrite (current CRAN harmony only exports
  # RunHarmony()) -- so IntegrateLayers() errors with "'HarmonyMatrix' is
  # not an exported object from 'namespace:harmony'" on any current harmony
  # install. Calling RunHarmony() ourselves sidesteps Seurat's broken
  # wrapper entirely and works regardless of harmony version.
  #
  # This also covers banksy = TRUE (validated above that integration must
  # be 'HarmonyIntegration' in that case): RunBanksy() produces a single
  # consolidated (non-split-layer) assay that IntegrateLayers()'s
  # anchor-based methods can't work with anyway, and RunHarmony() operates
  # directly on the PCA cell embeddings, so the split-layer requirement
  # doesn't apply to it (same as the BANKSY-Seurat vignette's own Harmony
  # workflow).
  if (identical(integration, 'HarmonyIntegration')) {
    if (!requireNamespace('harmony', quietly = TRUE)) {
      stop("'harmony' is required for integration = 'HarmonyIntegration'. ",
           "Install with: install.packages('harmony')")
    }
    harmony_reduction <- if (isTRUE(banksy)) 'pca_banksy' else integration_reduction
    message(sprintf('--- Running Harmony directly on the %s embedding ---',
                    harmony_reduction))
    obj <- harmony::RunHarmony(obj, group.by.vars = group_column,
                               reduction.use = harmony_reduction,
                               reduction.save = new_reduction)
  } else {
    message(sprintf('--- Integrating layers (method: %s) ---', integration))
    obj <- Seurat::IntegrateLayers(object = obj,
                                   method = integration,
                                   orig.reduction = integration_reduction,
                                   assay = integration_assay,
                                   normalization.method = integration_normalization,
                                   new.reduction = new_reduction,
                                   k.anchor = k_anchor,
                                   k.weight = k_weight)
  }

  # ---- Cluster + UMAP -----------------------------------------------------
  message('--- Finding neighbors and clusters ---')
  obj <- Seurat::FindNeighbors(obj, reduction = new_reduction, dims = dims_start:dims_to_use)
  obj <- Seurat::FindClusters(obj, resolution = cluster_resolution)

  message('--- Running UMAP ---')
  obj <- Seurat::RunUMAP(obj, reduction = new_reduction, dims = dims_start:dims_to_use)

  # ---- Save RDS (single, deduplicated path) -------------------------------
  if (isTRUE(save_rds_file)) {
    rds_prefix <- if (is.null(file_name)) new_reduction else file_name
    rds_path   <- paste(rds_prefix, 'merged_seurat_objects.rds', sep = '_')
    message(sprintf('--- Saving merged Seurat object to RDS (%s) ---', rds_path))
    saveRDS(obj, rds_path)
  }

  # ---- Cluster DimPlot ----------------------------------------------------
  message('--- Saving cluster DimPlot ---')
  Seurat::DimPlot(obj, label = TRUE, raster = FALSE)
  ggplot2::ggsave('dimplot_seurat_clusters.pdf', height = 5, width = 7)

  # ---- Markers (only place that early-returned before; now always reached)
  if (isTRUE(markers) && isTRUE(is_atac)) {
    message('--- Skipping FindAllMarkers(): ATAC input detected and RNA-tuned ',
            'thresholds (logfc.threshold, min.pct) don\'t apply to peak ',
            'accessibility. Call Signac::FindMarkers()/FindAllMarkers() ',
            'directly with ATAC-appropriate settings (e.g. ',
            'test.use = "LR", latent.vars = "peak_region_fragments") if you ',
            'need differential accessibility. ---')
  } else if (isTRUE(markers)) {
    message('--- Running FindAllMarkers ---')
    marker_results <- Seurat::FindAllMarkers(obj,
                                             logfc.threshold = 1,
                                             only.pos = TRUE,
                                             min.pct = 0.1,
                                             recorrect_umi = FALSE)
    utils::write.csv(marker_results, 'markers_all.csv')

    # Hand the table we already have to TopMarkerPlot() rather than
    # hand-rolling a DotPlot here. Passing `markers = marker_results` means
    # FindAllMarkers() is NOT re-run, and the panel comes out in the same
    # annotated style as every other marker plot in the package: right-edge
    # cluster labels, shared auto-sizing, and the assay-presence /
    # zero-expression / no-variance filtering that keeps blank and all-grey
    # rows out of the figure.
    #
    # Wrapped in tryCatch() because this is the LAST step of a long and
    # expensive function -- merging, integration, clustering and the RDS save
    # have all already happened by now. A marker panel that legitimately
    # cannot be drawn (nothing clears p_val_adj, or every candidate gene is
    # filtered out) must not throw away the object this function exists to
    # return.
    p_markers <- tryCatch(
      TopMarkerPlot(obj,
                    n         = marker_n,
                    markers   = marker_results,
                    save_path = 'marker_plot.pdf'),
      error = function(e) {
        message('--- Skipping marker_plot.pdf: ', conditionMessage(e), ' ---')
        NULL
      }
    )
  }

  return(obj)
}
