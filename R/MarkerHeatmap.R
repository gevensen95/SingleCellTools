#' Heatmap of top markers per cluster
#'
#' Picks the top \code{n} markers per cluster by \code{avg_log2FC}, averages
#' expression per cluster, and draws a heatmap. Works either from an existing
#' \code{FindAllMarkers} table, computes one on the fly, or -- via
#' \code{genes} -- plots an explicit, caller-supplied gene list instead of
#' any marker selection.
#'
#' @param obj A Seurat object.
#' @param markers Optional data frame from \code{FindAllMarkers} (must have
#'   columns \code{gene}, \code{cluster}, \code{avg_log2FC}, \code{p_val_adj}).
#'   If NULL, \code{FindAllMarkers} is run. Mutually exclusive with
#'   \code{genes}.
#' @param genes Optional character vector of gene names to plot directly,
#'   bypassing marker selection entirely (\code{n}, \code{p_val_adj} /
#'   \code{avg_log2FC} filtering, and \code{FindAllMarkers} are all skipped).
#'   Useful for a curated panel rather than a data-driven top-N marker list.
#'   Mutually exclusive with \code{markers}. Genes not found in \code{assay}
#'   are dropped with a message.
#' @param n Number of top markers per cluster to display. Ignored when
#'   \code{genes} is supplied. Default 10.
#' @param assay Assay to read from. Default DefaultAssay.
#' @param pseudobulk Logical; if TRUE, sum raw counts per cluster first
#'   (\code{Seurat::AggregateExpression}), normalize once, and z-score that
#'   instead of \code{AverageExpression}'s per-cell-averaged values. Default
#'   \code{FALSE}. See Details.
#' @param scale_rows Logical; z-score each gene across clusters. Default TRUE.
#' @param colors Diverging color vector for the gradient. Default RdBu.
#' @param cluster_rows,cluster_cols Logical; hierarchical clustering of
#'   rows / columns. Defaults: rows = TRUE, cols = FALSE (clusters in factor
#'   order).
#' @return A \code{ggplot} object.
#' @details
#' \strong{\code{pseudobulk = TRUE}.} \code{AverageExpression()} (the
#' default) averages already-normalized per-cell values, which is fast but
#' still exposed to per-cell/per-spot sampling noise -- a real concern on
#' sparse data (Visium spots routinely carry only a few hundred to a few
#' thousand UMIs). Pseudobulking sums \emph{raw} counts across every cell
#' in a cluster into one profile first, then normalizes once (CPM-style,
#' scale factor 1e4, \code{log1p}) -- the statistically preferred order of
#' operations (summing before normalizing, not averaging already-normalized
#' values), and the same noise-collapsing idea \code{\link{PseudobulkDE}}
#' uses for differential expression. Z-scoring that pseudobulk profile
#' across clusters is what actually protects against a broadly/ambiently
#' expressed gene dominating every cluster -- pseudobulking alone reduces
#' \emph{stochastic} noise, not a systematic contamination signal, so it
#' doesn't replace the z-scoring step, it makes the values being z-scored
#' more reliable. Slower and more memory-heavy than the default: the whole
#' assay is aggregated (via \code{Seurat::AggregateExpression}, which
#' doesn't reliably support a \code{features} subset across Seurat
#' versions) before being subset down to the marker genes actually needed.
#' @importFrom Seurat DefaultAssay AverageExpression AggregateExpression FindAllMarkers
#' @importFrom SeuratObject LayerData
#' @importFrom dplyr group_by slice_max ungroup arrange
#' @importFrom ggplot2 ggplot aes geom_tile theme_minimal theme element_text scale_fill_gradientn labs
#' @importFrom RColorBrewer brewer.pal
#' @export
MarkerHeatmap <- function(obj,
                          markers       = NULL,
                          genes         = NULL,
                          n             = 10,
                          assay         = NULL,
                          pseudobulk    = FALSE,
                          scale_rows    = TRUE,
                          colors        = NULL,
                          cluster_rows  = TRUE,
                          cluster_cols  = FALSE) {

  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (!is.null(markers) && !is.null(genes)) {
    stop("Supply either `markers` or `genes`, not both.")
  }
  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay

  if (!is.null(genes)) {
    if (!is.character(genes) || !length(genes)) {
      stop("`genes` must be a non-empty character vector of gene names.")
    }
    top_genes <- unique(genes)
  } else {
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
  }

  # Drop any requested gene not present in the target assay -- applies to
  # both marker-derived and caller-supplied `genes` lists (a hand-typed
  # `markers` table or `genes` vector can easily contain a typo/absent gene;
  # FindAllMarkers()-derived tables never do, so this is a no-op there).
  assay_genes <- rownames(obj[[a]])
  keep_genes  <- intersect(top_genes, assay_genes)
  if (!length(keep_genes)) {
    stop("None of the requested genes are present in assay '", a, "'.")
  }
  if (length(keep_genes) < length(top_genes)) {
    message(sprintf(
      "  MarkerHeatmap: %d/%d gene(s) not found in assay '%s': %s",
      length(top_genes) - length(keep_genes), length(top_genes), a,
      paste(setdiff(top_genes, keep_genes), collapse = ", ")))
  }
  top_genes <- keep_genes

  if (isTRUE(pseudobulk)) {
    # True pseudobulk: sum raw counts per cluster first, then normalize once.
    # AggregateExpression() can return a plain matrix, a dgCMatrix, or an
    # Assay/Assay5 depending on Seurat version -- same coercion PseudobulkDE()
    # already relies on for the same function.
    agg_raw <- Seurat::AggregateExpression(obj, assays = a,
                                           return.seurat = FALSE,
                                           slot = "counts")[[a]]
    if (inherits(agg_raw, c("Assay", "Assay5"))) {
      agg_raw <- SeuratObject::LayerData(agg_raw, layer = "counts")
    }
    agg_raw <- if (is.matrix(agg_raw)) agg_raw else as.matrix(agg_raw)

    # top_genes is already filtered to genes present in `a`'s rownames above,
    # but AggregateExpression()'s output can in principle be missing a gene
    # that was present pre-aggregation (e.g. it got dropped for having zero
    # counts everywhere) -- re-intersect defensively rather than assume.
    keep_genes <- intersect(top_genes, rownames(agg_raw))
    if (!length(keep_genes)) {
      stop("None of the requested genes are present in the pseudobulk ",
           "assay '", a, "'.")
    }
    if (length(keep_genes) < length(top_genes)) {
      message(sprintf(
        "  pseudobulk: %d/%d gene(s) dropped during aggregation from assay '%s': %s",
        length(top_genes) - length(keep_genes), length(top_genes), a,
        paste(setdiff(top_genes, keep_genes), collapse = ", ")))
    }
    agg_raw <- agg_raw[keep_genes, , drop = FALSE]

    # CPM-style normalize per pseudobulk cluster, then log1p -- matches this
    # package's other log-normalization conventions (scale.factor = 1e4),
    # just applied once per cluster instead of once per cell.
    lib_size <- Matrix::colSums(agg_raw)
    avg <- sweep(agg_raw, 2, pmax(lib_size, 1), "/") * 1e4
    avg <- log1p(avg)
  } else {
    # Average expression per cluster
    avg <- Seurat::AverageExpression(obj, features = top_genes,
                                     assays = a, return.seurat = FALSE)[[a]]
    # Guard rather than unconditional as.matrix(): AverageExpression() output
    # here is already small (clusters x marker genes), so this isn't a real
    # memory concern, but skip the coercion/copy when it's already a base
    # matrix rather than always paying for one.
    if (!is.matrix(avg)) avg <- as.matrix(avg)
  }

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
    Ol_Reliable() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid  = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = NULL,
                  fill = if (isTRUE(scale_rows)) {
                    "Z-score"
                  } else if (isTRUE(pseudobulk)) {
                    "Pseudobulk log-CPM"
                  } else {
                    "Avg. expr."
                  })
}
