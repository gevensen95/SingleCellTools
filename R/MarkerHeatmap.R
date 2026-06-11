#' Heatmap of top markers per cluster
#'
#' Picks the top \code{n} markers per cluster by \code{avg_log2FC}, averages
#' expression per cluster, and draws a heatmap. Works either from an existing
#' \code{FindAllMarkers} table or computes one on the fly.
#'
#' @param obj A Seurat object.
#' @param markers Optional data frame from \code{FindAllMarkers} (must have
#'   columns \code{gene}, \code{cluster}, \code{avg_log2FC}, \code{p_val_adj}).
#'   If NULL, \code{FindAllMarkers} is run.
#' @param n Number of top markers per cluster to display. Default 10.
#' @param assay Assay to read from. Default DefaultAssay.
#' @param scale_rows Logical; z-score each gene across clusters. Default TRUE.
#' @param colors Diverging color vector for the gradient. Default RdBu.
#' @param cluster_rows,cluster_cols Logical; hierarchical clustering of
#'   rows / columns. Defaults: rows = TRUE, cols = FALSE (clusters in factor
#'   order).
#' @return A \code{ggplot} object.
#' @importFrom Seurat DefaultAssay AverageExpression FindAllMarkers
#' @importFrom dplyr group_by slice_max ungroup arrange
#' @importFrom ggplot2 ggplot aes geom_tile theme_minimal theme element_text scale_fill_gradientn labs
#' @importFrom RColorBrewer brewer.pal
#' @export
MarkerHeatmap <- function(obj,
                          markers       = NULL,
                          n             = 10,
                          assay         = NULL,
                          scale_rows    = TRUE,
                          colors        = NULL,
                          cluster_rows  = TRUE,
                          cluster_cols  = FALSE) {

  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay

  if (is.null(markers)) {
    message("--- Running FindAllMarkers ---")
    markers <- Seurat::FindAllMarkers(obj, only.pos = TRUE,
                                      min.pct = 0.25,
                                      logfc.threshold = 0.25,
                                      verbose = FALSE)
  }
  required <- c("gene", "cluster", "avg_log2FC", "p_val_adj")
  if (!all(required %in% colnames(markers))) {
    stop("`markers` is missing required columns: ",
         paste(setdiff(required, colnames(markers)), collapse = ", "))
  }

  cluster <- avg_log2FC <- p_val_adj <- NULL  # NSE
  top <- markers %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(avg_log2FC, n = n, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(cluster, dplyr::desc(avg_log2FC))

  top_genes <- unique(top$gene)
  if (!length(top_genes)) stop("No genes passed the marker filter.")

  # Average expression per cluster
  avg <- Seurat::AverageExpression(obj, features = top_genes,
                                   assays = a, return.seurat = FALSE)[[a]]
  avg <- as.matrix(avg)

  if (isTRUE(scale_rows)) {
    avg <- t(scale(t(avg)))
    avg[!is.finite(avg)] <- 0
  }

  # Optional clustering for row / col order
  if (isTRUE(cluster_rows) && nrow(avg) > 2) {
    avg <- avg[stats::hclust(stats::dist(avg))$order, , drop = FALSE]
  }
  if (isTRUE(cluster_cols) && ncol(avg) > 2) {
    avg <- avg[, stats::hclust(stats::dist(t(avg)))$order, drop = FALSE]
  }

  # Long form for ggplot
  long <- data.frame(
    gene    = rep(rownames(avg), times = ncol(avg)),
    cluster = rep(colnames(avg), each = nrow(avg)),
    value   = as.vector(avg),
    stringsAsFactors = FALSE
  )
  long$gene    <- factor(long$gene,    levels = rownames(avg))
  long$cluster <- factor(long$cluster, levels = colnames(avg))

  if (is.null(colors)) {
    colors <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
  }

  cluster_ <- gene_ <- value <- NULL  # NSE
  ggplot2::ggplot(long, ggplot2::aes(x = cluster, y = gene, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colors = colors) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid  = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = NULL,
                  fill = if (isTRUE(scale_rows)) "Z-score" else "Avg. expr.")
}
