#' Pseudotime inference (slingshot / monocle3 / destiny)
#'
#' Fits a pseudotemporal ordering of cells and writes it into the Seurat
#' object's metadata, using one of three backends. Every method also
#' writes a standardized structure to \code{obj@misc$pseudotime} --
#' \code{list(method, pseudotime, weights)}, both matrices cells x lineages
#' -- so downstream tools (\code{\link{RunTradeSeqDE}}) can consume the
#' result without caring which method produced it.
#'
#' \describe{
#'   \item{\code{"slingshot"} (default)}{Fits lineage curves through a
#'     reduced-dimensional embedding using cluster labels as anchors. Good
#'     for continuous developmental/differentiation processes; handles
#'     branching natively and returns one pseudotime column per detected
#'     lineage, with real per-cell/per-lineage soft weights
#'     (\code{slingshot::slingCurveWeights}).}
#'   \item{\code{"monocle3"}}{Learns a principal graph over the embedding
#'     (\code{monocle3::learn_graph}) and orders cells along it from a root
#'     defined by \code{start_cluster}. Handles branching structurally (the
#'     graph itself branches), but this wrapper reports one combined
#'     pseudotime column (\code{"Lineage1"}) rather than per-branch curves --
#'     inspect \code{obj@misc$monocle3} directly for the full graph.
#'     Requires the \code{monocle3} and \code{SeuratWrappers} packages.}
#'   \item{\code{"destiny"}}{Diffusion pseudotime (DPT) on a diffusion map
#'     of the expression data. Lightweight and fast, best for simple
#'     continuous trajectories without complex branching; reports a single
#'     \code{"Lineage1"} pseudotime column. Requires the \code{destiny}
#'     package.}
#' }
#'
#' @param obj A Seurat object with a reduction and clusters computed.
#' @param method \code{"slingshot"} (default), \code{"monocle3"}, or
#'   \code{"destiny"}.
#' @param reduction Reduction to run pseudotime on (\code{"slingshot"}/
#'   \code{"monocle3"}). UMAP is common for visualization; PCA/harmony can
#'   give more stable topology. Default \code{"umap"}.
#' @param dims Number of dimensions from \code{reduction} to use. Default 2
#'   for UMAP; bump to 10-30 for PCA/harmony.
#' @param cluster_col Metadata column holding cluster labels. Default
#'   \code{"seurat_clusters"}.
#' @param start_cluster Optional cluster id to fix as the trajectory root
#'   (\code{"slingshot"}/\code{"monocle3"}). \code{NULL} lets the method
#'   choose.
#' @param end_cluster \code{"slingshot"} only: optional cluster id(s) to
#'   constrain as terminal endpoint(s).
#' @param assay \code{"destiny"} only: assay to compute the diffusion map
#'   from. Default \code{DefaultAssay(obj)}.
#' @param prefix Column name prefix for the per-lineage pseudotime metadata
#'   columns. \code{NULL} (default) uses \code{method}.
#' @param verbose Message progress. Default \code{TRUE}.
#' @return The Seurat object with new \code{<prefix>_<Lineage>} metadata
#'   columns, \code{obj@misc$pseudotime} (the standardized
#'   \code{list(method, pseudotime, weights)} described above), and the
#'   method's native fitted object under \code{obj@misc$slingshot} /
#'   \code{obj@misc$monocle3} / \code{obj@misc$destiny}.
#' @examples
#' \dontrun{
#' obj <- PseudotimeWrapper(obj, method = "slingshot",
#'                          cluster_col = "seurat_clusters",
#'                          start_cluster = "3")
#' FeaturePlot(obj, features = "slingshot_Lineage1")
#'
#' # Feed straight into RunTradeSeqDE() -- no extra plumbing needed
#' obj <- PseudotimeWrapper(obj, method = "monocle3", start_cluster = "3")
#' de  <- RunTradeSeqDE(obj)
#' }
#' @importFrom Seurat DefaultAssay Embeddings GetAssayData
#' @export
PseudotimeWrapper <- function(obj,
                              method        = c("slingshot", "monocle3", "destiny"),
                              reduction     = "umap",
                              dims          = 2,
                              cluster_col   = "seurat_clusters",
                              start_cluster = NULL,
                              end_cluster   = NULL,
                              assay         = NULL,
                              prefix        = NULL,
                              verbose       = TRUE) {

  method <- match.arg(method)
  .assert_seurat(obj)
  if (!cluster_col %in% colnames(obj@meta.data)) {
    stop("Cluster column '", cluster_col, "' not found in metadata.")
  }
  px <- if (is.null(prefix)) method else prefix

  if (method == "slingshot") {
    if (!requireNamespace("slingshot", quietly = TRUE)) {
      stop("'slingshot' is required. Install with ",
           "BiocManager::install('slingshot').")
    }
    if (!(reduction %in% names(obj@reductions))) {
      stop("Reduction '", reduction, "' not found.")
    }
    emb <- Seurat::Embeddings(obj, reduction = reduction)[, seq_len(dims), drop = FALSE]
    clusters <- as.character(obj@meta.data[[cluster_col]])

    if (isTRUE(verbose)) {
      message(sprintf("--- Running slingshot on '%s' (%d dims, %d clusters) ---",
                      reduction, dims, length(unique(clusters))))
    }
    sds <- slingshot::slingshot(data = emb, clusterLabels = clusters,
                                start.clus = start_cluster, end.clus = end_cluster)

    pt <- slingshot::slingPseudotime(sds)
    wt <- tryCatch(slingshot::slingCurveWeights(sds),
                   error = function(e) !is.na(pt) * 1)
    colnames(pt) <- colnames(wt) <- paste0("Lineage", seq_len(ncol(pt)))
    pt <- pt[colnames(obj), , drop = FALSE]
    wt <- wt[colnames(obj), , drop = FALSE]

    obj@misc$slingshot <- sds
    n_lineages <- ncol(pt)

  } else if (method == "monocle3") {
    if (!requireNamespace("monocle3", quietly = TRUE) ||
        !requireNamespace("SeuratWrappers", quietly = TRUE)) {
      stop("'monocle3' and 'SeuratWrappers' are required for method = ",
           "'monocle3'. Install with remotes::install_github(",
           "c('cole-trapnell-lab/monocle3', 'satijalab/seurat-wrappers')).")
    }
    if (!(reduction %in% names(obj@reductions))) {
      stop("Reduction '", reduction, "' not found.")
    }
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("'SingleCellExperiment' is required for method = 'monocle3'. ",
           "Install with BiocManager::install('SingleCellExperiment').")
    }
    if (isTRUE(verbose)) message("--- Building monocle3 cell_data_set ---")
    cds <- SeuratWrappers::as.cell_data_set(obj)
    emb <- Seurat::Embeddings(obj, reduction = reduction)[, seq_len(dims), drop = FALSE]
    SingleCellExperiment::reducedDim(cds, "UMAP") <- emb

    if (isTRUE(verbose)) message("--- Learning principal graph ---")
    cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP")
    cds <- monocle3::learn_graph(cds, use_partition = FALSE)

    root_cells <- NULL
    if (!is.null(start_cluster)) {
      root_cells <- colnames(obj)[as.character(obj@meta.data[[cluster_col]]) == start_cluster]
      if (length(root_cells) == 0) {
        stop("No cells found in start_cluster = '", start_cluster, "'.")
      }
    }
    cds <- monocle3::order_cells(cds, reduction_method = "UMAP",
                                 root_cells = root_cells)

    pt_vec <- monocle3::pseudotime(cds)
    pt <- matrix(pt_vec, ncol = 1, dimnames = list(names(pt_vec), "Lineage1"))
    pt[!is.finite(pt)] <- NA
    wt <- matrix(as.numeric(!is.na(pt)), ncol = 1, dimnames = dimnames(pt))
    pt <- pt[colnames(obj), , drop = FALSE]
    wt <- wt[colnames(obj), , drop = FALSE]

    obj@misc$monocle3 <- cds
    n_lineages <- 1

  } else {
    if (!requireNamespace("destiny", quietly = TRUE)) {
      stop("'destiny' is required for method = 'destiny'. Install with ",
           "BiocManager::install('destiny').")
    }
    a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay
    expr <- t(as.matrix(Seurat::GetAssayData(obj, assay = a, layer = "data")))

    if (isTRUE(verbose)) message("--- Fitting diffusion map (destiny) ---")
    dm <- destiny::DiffusionMap(expr)

    tip <- NULL
    if (!is.null(start_cluster)) {
      root_cells <- rownames(expr)[as.character(obj@meta.data[[cluster_col]]) == start_cluster]
      if (length(root_cells) == 0) {
        stop("No cells found in start_cluster = '", start_cluster, "'.")
      }
      tip <- which(rownames(expr) == root_cells[1])
    }
    dpt <- if (is.null(tip)) destiny::DPT(dm) else destiny::DPT(dm, tips = tip)
    pt_vec <- dpt$dpt
    names(pt_vec) <- rownames(expr)
    pt <- matrix(pt_vec, ncol = 1, dimnames = list(names(pt_vec), "Lineage1"))
    wt <- matrix(as.numeric(!is.na(pt)), ncol = 1, dimnames = dimnames(pt))
    pt <- pt[colnames(obj), , drop = FALSE]
    wt <- wt[colnames(obj), , drop = FALSE]

    obj@misc$destiny <- list(diffusion_map = dm, dpt = dpt)
    n_lineages <- 1
  }

  for (col in colnames(pt)) {
    obj@meta.data[[paste0(px, "_", col)]] <- pt[colnames(obj), col]
  }
  obj@misc$pseudotime <- list(method = method, pseudotime = pt, weights = wt)

  if (isTRUE(verbose)) {
    message(sprintf("  %d lineage(s) fit (%s): %s", n_lineages, method,
                    paste0(px, "_", colnames(pt), collapse = ", ")))
  }
  obj
}
