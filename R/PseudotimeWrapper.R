#' Pseudotime inference with slingshot
#'
#' Wraps \code{slingshot::slingshot} to fit lineage curves through a
#' reduced-dimensional embedding using cluster labels as anchors, and
#' writes the resulting per-cell pseudotime values into the Seurat
#' object's metadata (one column per lineage).
#'
#' Slingshot works well for continuous developmental / differentiation
#' processes and doesn't require a graph pre-processing step. For
#' branching topologies with multiple end states, it returns one
#' pseudotime column per detected lineage (\code{Lineage1},
#' \code{Lineage2}, ...).
#'
#' @param obj A Seurat object with a reduction and clusters computed.
#' @param reduction Reduction to run pseudotime on. UMAP is common for
#'   visualization but PCA / harmony can give more stable topology.
#'   Default \code{"umap"}.
#' @param dims Number of dimensions from \code{reduction} to use. Default 2
#'   for UMAP; bump to 10-30 for PCA / harmony.
#' @param cluster_col Metadata column holding cluster labels. Default
#'   \code{"seurat_clusters"}.
#' @param start_cluster Optional cluster id to fix as the trajectory root.
#'   \code{NULL} lets slingshot choose (usually the endpoint of a
#'   minimum-spanning tree).
#' @param end_cluster Optional cluster id(s) to constrain as terminal
#'   endpoint(s).
#' @param prefix Column name prefix for the pseudotime metadata. Default
#'   \code{"slingshot"}. Columns \code{<prefix>_Lineage1}, ... are added.
#' @return The Seurat object with new metadata columns and
#'   \code{obj@misc$slingshot} containing the full SlingshotDataSet.
#' @examples
#' \dontrun{
#' obj <- PseudotimeWrapper(obj,
#'                          reduction     = "umap",
#'                          cluster_col   = "seurat_clusters",
#'                          start_cluster = "3")
#' # Plot one lineage on UMAP
#' FeaturePlot(obj, features = "slingshot_Lineage1")
#'
#' # Access the trajectory itself
#' obj@misc$slingshot
#' }
#' @importFrom Seurat Embeddings
#' @export
PseudotimeWrapper <- function(obj,
                              reduction     = "umap",
                              dims          = 2,
                              cluster_col   = "seurat_clusters",
                              start_cluster = NULL,
                              end_cluster   = NULL,
                              prefix        = "slingshot") {

  if (!requireNamespace("slingshot", quietly = TRUE)) {
    stop("'slingshot' is required. Install with ",
         "BiocManager::install('slingshot').")
  }
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (!(reduction %in% names(obj@reductions))) {
    stop("Reduction '", reduction, "' not found.")
  }
  if (!cluster_col %in% colnames(obj@meta.data)) {
    stop("Cluster column '", cluster_col, "' not found in metadata.")
  }

  emb <- Seurat::Embeddings(obj, reduction = reduction)[, seq_len(dims),
                                                        drop = FALSE]
  clusters <- as.character(obj@meta.data[[cluster_col]])

  message(sprintf("--- Running slingshot on '%s' (%d dims, %d clusters) ---",
                  reduction, dims, length(unique(clusters))))
  sds <- slingshot::slingshot(
    data          = emb,
    clusterLabels = clusters,
    start.clus    = start_cluster,
    end.clus      = end_cluster
  )

  # slingshot pseudotime: cells x lineages
  pt <- slingshot::slingPseudotime(sds)
  colnames(pt) <- paste0(prefix, "_", colnames(pt))
  # Some cells may have NA for lineages they aren't part of; that's expected.
  for (col in colnames(pt)) {
    obj@meta.data[[col]] <- pt[colnames(obj), col]
  }

  obj@misc$slingshot <- sds
  message(sprintf("  %d lineage(s) fit: %s", ncol(pt),
                  paste(colnames(pt), collapse = ", ")))
  obj
}
