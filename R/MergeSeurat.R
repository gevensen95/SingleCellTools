#' Merge and Integrate Seurat Objects
#'
#' Merges a list of Seurat objects, normalizes (SCT or LogNormalize), runs PCA,
#' integrates layers, finds neighbors/clusters, runs UMAP, and (optionally)
#' computes per-cluster markers.
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
#' @param integration_normalization,integration_assay,integration_reduction
#'   Arguments passed through to `IntegrateLayers`.
#' @param new_reduction Name of the reduction created by integration.
#' @param k_anchor,k_weight Arguments required for RPCA/CCA/JointPCA integration.
#' @param markers Logical; if TRUE, run `FindAllMarkers` and save outputs.
#' @param group_column Metadata column used to compute median nCount per group
#'   for `SCTransform`'s scale_factor.
#' @param common_genes_only Logical; if TRUE, subset all input objects to genes
#'   present in every object before merging.
#' @param common_genes_assay Assay to inspect for gene names when
#'   `common_genes_only = TRUE`.
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
                        group_column = 'orig.ident',
                        common_genes_only = FALSE,
                        common_genes_assay = NULL) {

  # ---- Argument validation ------------------------------------------------
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
    seurat_objects <- lapply(seq_along(seurat_objects), function(i) {
      o <- seurat_objects[[i]]
      lab <- obj_labels[i]

      o$Counts <- colSums(SeuratObject::LayerData(o,
                                                  assay = common_genes_assay,
                                                  layer = 'counts'))
      n_before <- ncol(o)
      n_drop   <- sum(o$Counts == 0)
      if (n_drop > 0) {
        o <- subset(o, subset = Counts > 0)
      }
      o$Counts <- NULL

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
    obj[["RNA"]] <- as(object = obj[["Spatial"]], Class = "Assay5")
  } else if (spatial == 'Xenium') {
    obj <- suppressWarnings(merge(seurat_objects[[1]], seurat_objects[-1],
                                  add.cell.ids = cell_IDs))
    obj[["RNA"]] <- as(object = obj[["Xenium"]], Class = "Assay5")
  } else if (spatial == 'no') {
    obj <- merge(seurat_objects[[1]], seurat_objects[-1],
                 add.cell.ids = cell_IDs)
  }

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

  # ---- PCA ----------------------------------------------------------------
  message('--- Running PCA ---')
  obj <- Seurat::RunPCA(obj)

  # ---- Choose dims (elbow prompt or fixed max_dims) -----------------------
  if (use_elbow_plot) {
    elbow_plot <- Seurat::ElbowPlot(obj)
    print(elbow_plot)
    dims_to_use <- as.numeric(readline(prompt = 'Enter # of PCs: '))
  } else {
    dims_to_use <- max_dims
  }

  # ---- Integrate ----------------------------------------------------------
  message(sprintf('--- Integrating layers (method: %s) ---', integration))
  obj <- Seurat::IntegrateLayers(object = obj,
                                 method = integration,
                                 orig.reduction = integration_reduction,
                                 assay = integration_assay,
                                 normalization.method = integration_normalization,
                                 new.reduction = new_reduction,
                                 k.anchor = k_anchor,
                                 k.weight = k_weight)

  # ---- Cluster + UMAP -----------------------------------------------------
  message('--- Finding neighbors and clusters ---')
  obj <- Seurat::FindNeighbors(obj, reduction = new_reduction, dims = 1:dims_to_use)
  obj <- Seurat::FindClusters(obj, resolution = cluster_resolution)

  message('--- Running UMAP ---')
  obj <- Seurat::RunUMAP(obj, reduction = new_reduction, dims = 1:dims_to_use)

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
  if (isTRUE(markers)) {
    message('--- Running FindAllMarkers ---')
    marker_results <- Seurat::FindAllMarkers(obj,
                                             logfc.threshold = 1,
                                             only.pos = TRUE,
                                             min.pct = 0.1,
                                             recorrect_umi = FALSE)
    utils::write.csv(marker_results, 'markers_all.csv')

    markers_filtered <- marker_results %>%
      dplyr::filter(p_val_adj < 0.05) %>%
      dplyr::arrange(-avg_log2FC) %>%
      dplyr::group_by(cluster) %>%
      dplyr::slice_max(avg_log2FC, n = 10)

    Seurat::DotPlot(obj, features = unique(markers_filtered$gene)) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = '', y = '')
    ggplot2::ggsave('marker_plot.pdf', width = 10, height = 10, units = 'in')
  }

  return(obj)
}
