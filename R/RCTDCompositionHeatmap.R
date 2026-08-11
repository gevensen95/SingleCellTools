#' Heatmap of cluster-level RCTD composition
#'
#' Visualizes the cluster x cell-type proportion matrix from
#' \code{\link{SummarizeRCTDByCluster}} as a heatmap. In liver (and other
#' tissue with one very abundant cell type), the dominant type's proportion
#' is often high and roughly stable across every cluster -- real biology,
#' not an artifact -- which buries the actually differentiating minority
#' cell types under that dominant type's raw magnitude. \code{scale_cols =
#' TRUE} (default) z-scores each cell type across clusters before plotting
#' so a cell type's \emph{relative} enrichment in a given cluster is what
#' shows up, regardless of its absolute proportion.
#'
#' @param x A Seurat object, OR the list returned by
#'   \code{\link{SummarizeRCTDByCluster}} (i.e. has a \code{composition}
#'   element), OR a cluster x cell-type numeric matrix/data frame directly
#'   (e.g. \code{res$composition} itself).
#' @param cluster_col,weight_cols,min_score,min_margin,unassigned_label
#'   Used only when \code{x} is a Seurat object -- passed through to
#'   \code{\link{SummarizeRCTDByCluster}} to compute the composition matrix.
#' @param scale_cols Logical; z-score each cell type (column) across
#'   clusters before plotting. Default \code{TRUE}. Columns with zero
#'   variance (e.g. a cell type absent everywhere) become \code{0} rather
#'   than \code{NaN}.
#' @param colors Diverging color vector for the gradient. Default RdBu.
#' @param cluster_rows,cluster_cols Logical; hierarchical clustering of
#'   rows / columns. Defaults: rows = TRUE, cols = FALSE -- matching
#'   \code{\link{MarkerHeatmap}}'s defaults.
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' visium <- RunRCTD(visium, reference = ref, celltype_col = "cell_type")
#' RCTDCompositionHeatmap(visium)
#'
#' # Or reuse an already-computed result instead of re-aggregating
#' res <- SummarizeRCTDByCluster(visium)
#' RCTDCompositionHeatmap(res)
#' RCTDCompositionHeatmap(res$composition, scale_cols = FALSE)  # raw proportions
#' }
#' @importFrom ggplot2 ggplot aes geom_tile theme element_text scale_fill_gradientn labs
#' @importFrom RColorBrewer brewer.pal
#' @export
RCTDCompositionHeatmap <- function(x,
                                   cluster_col      = "seurat_clusters",
                                   weight_cols      = NULL,
                                   min_score        = NA,
                                   min_margin       = NA,
                                   unassigned_label = "Unknown",
                                   scale_cols       = TRUE,
                                   colors           = NULL,
                                   cluster_rows     = TRUE,
                                   cluster_cols     = FALSE) {

  # ---- Resolve `x` to a plain cluster x cell-type matrix -------------------
  if (inherits(x, "Seurat")) {
    res <- SummarizeRCTDByCluster(x, cluster_col = cluster_col,
                                  weight_cols = weight_cols,
                                  min_score = min_score, min_margin = min_margin,
                                  unassigned_label = unassigned_label)
    comp <- res$composition
  } else if (is.list(x) && !is.data.frame(x) && "composition" %in% names(x)) {
    comp <- x$composition
  } else {
    comp <- x
  }
  if (!is.matrix(comp)) comp <- as.matrix(comp)
  if (!is.numeric(comp)) {
    stop("Could not resolve `x` to a numeric cluster x cell-type matrix. ",
         "Pass a Seurat object, the list from SummarizeRCTDByCluster(), or ",
         "a numeric matrix/data frame directly.")
  }
  if (nrow(comp) < 2) {
    stop("Need at least 2 clusters to plot a heatmap; got ", nrow(comp), ".")
  }

  if (isTRUE(scale_cols)) {
    comp <- scale(comp)
    comp[!is.finite(comp)] <- 0
  }

  # Optional clustering for row / col order -- same pattern as MarkerHeatmap()
  if (isTRUE(cluster_rows) && nrow(comp) > 2) {
    comp <- comp[stats::hclust(stats::dist(comp))$order, , drop = FALSE]
  }
  if (isTRUE(cluster_cols) && ncol(comp) > 2) {
    comp <- comp[, stats::hclust(stats::dist(t(comp)))$order, drop = FALSE]
  }

  # Long form for ggplot
  long <- data.frame(
    cell_type = rep(colnames(comp), each = nrow(comp)),
    cluster   = rep(rownames(comp), times = ncol(comp)),
    value     = as.vector(comp),
    stringsAsFactors = FALSE
  )
  long$cell_type <- factor(long$cell_type, levels = colnames(comp))
  long$cluster   <- factor(long$cluster,   levels = rownames(comp))

  if (is.null(colors)) {
    colors <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
  }

  cluster_ <- cell_type_ <- value <- NULL  # NSE
  ggplot2::ggplot(long, ggplot2::aes(x = cluster, y = cell_type, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colors = colors) +
    Ol_Reliable() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid  = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = NULL,
                  fill = if (isTRUE(scale_cols)) "Z-score" else "Proportion")
}
