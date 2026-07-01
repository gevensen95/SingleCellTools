#' Per-sample cell-type composition
#'
#' Computes per-sample proportions of each cluster / cell type, returning
#' a tidy data frame and (optionally) a ggplot. Pairs naturally with
#' \code{CompositionalTest} for statistical inference.
#'
#' @param obj A Seurat object with per-cell metadata for sample and cluster.
#' @param cluster_col Metadata column holding cluster / cell-type ids.
#'   Default \code{"seurat_clusters"}.
#' @param sample_col Metadata column identifying the biological replicate
#'   (typically \code{"orig.ident"}).
#' @param group_col Optional metadata column for a grouping / condition
#'   (e.g. \code{"treatment"}). If supplied, plot styles that support it
#'   (\code{"box"}, \code{"line"}) show per-group distributions.
#' @param style Plot style. One of \code{"stack"} (per-sample stacked
#'   proportions, default), \code{"box"} (proportion per group as
#'   boxplots, one panel per cluster), \code{"line"} (mean +/- SE per
#'   group per cluster), \code{"none"} (no plot).
#' @param normalize Denominator for proportions. \code{"sample"} (default)
#'   normalizes each sample to sum to 1. \code{"cluster"} normalizes each
#'   cluster to sum to 1 across samples (rarer, useful when you want to
#'   compare sample contribution to a cluster).
#' @param colors Optional vector of colors for the cluster fill. \code{NULL}
#'   uses ggplot2 default.
#' @return If \code{style = "none"}, a data frame with columns
#'   \code{sample}, \code{cluster}, \code{n_cells},
#'   \code{n_sample_total}, \code{prop}, (and \code{group} if supplied).
#'   Otherwise, a list \code{list(df, plot)}.
#' @examples
#' \dontrun{
#' out <- CellComposition(obj,
#'                        sample_col = "orig.ident",
#'                        group_col  = "treatment")
#' out$plot
#' head(out$df)
#'
#' # Just the numbers
#' comp_df <- CellComposition(obj, sample_col = "orig.ident", style = "none")
#' }
#' @importFrom ggplot2 aes element_text facet_wrap geom_boxplot geom_col geom_line geom_pointrange ggplot labs position_dodge scale_fill_manual stat_summary theme theme_bw
#' @export
CellComposition <- function(obj,
                            cluster_col = "seurat_clusters",
                            sample_col  = "orig.ident",
                            group_col   = NULL,
                            style       = c("stack", "box", "line", "none"),
                            normalize   = c("sample", "cluster"),
                            colors      = NULL) {

  style     <- match.arg(style)
  normalize <- match.arg(normalize)
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  for (col in c(cluster_col, sample_col)) {
    if (!col %in% colnames(obj@meta.data)) {
      stop("Column '", col, "' not found in obj@meta.data.")
    }
  }
  if (!is.null(group_col) &&
      !group_col %in% colnames(obj@meta.data)) {
    stop("`group_col` '", group_col, "' not found in obj@meta.data.")
  }

  md <- obj@meta.data
  df <- as.data.frame(table(
    sample  = as.character(md[[sample_col]]),
    cluster = as.character(md[[cluster_col]])
  ), responseName = "n_cells", stringsAsFactors = FALSE)

  if (normalize == "sample") {
    tot <- stats::aggregate(n_cells ~ sample, data = df, FUN = sum)
    names(tot)[2] <- "n_sample_total"
    df <- merge(df, tot, by = "sample")
    df$prop <- df$n_cells / pmax(1, df$n_sample_total)
  } else {
    tot <- stats::aggregate(n_cells ~ cluster, data = df, FUN = sum)
    names(tot)[2] <- "n_sample_total"   # reused col name for symmetry
    df <- merge(df, tot, by = "cluster")
    df$prop <- df$n_cells / pmax(1, df$n_sample_total)
  }

  if (!is.null(group_col)) {
    grp_map <- unique(md[, c(sample_col, group_col)])
    grp_map <- setNames(as.character(grp_map[[group_col]]),
                        as.character(grp_map[[sample_col]]))
    df$group <- unname(grp_map[df$sample])
  }

  df <- df[order(df$sample, df$cluster), ]
  rownames(df) <- NULL

  if (style == "none") return(df)

  p <- switch(
    style,
    stack = .plot_composition_stack(df, group_col, colors),
    box   = .plot_composition_box(df, group_col, colors),
    line  = .plot_composition_line(df, group_col, colors)
  )
  list(df = df, plot = p)
}


# ============================================================================
# Style: per-sample stacked bar of cluster proportions
# ============================================================================
#' @keywords internal
#' @noRd
.plot_composition_stack <- function(df, group_col, colors) {
  sample <- prop <- cluster <- group <- NULL  # NSE
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = sample, y = prop, fill = cluster)) +
    ggplot2::geom_col(color = "black", linewidth = 0.15) +
    ggplot2::labs(x = NULL, y = "Proportion", fill = "Cluster") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  if (!is.null(group_col)) {
    p <- p + ggplot2::facet_wrap(~ group, scales = "free_x")
  }
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }
  p
}


# ============================================================================
# Style: per-cluster boxplot of per-sample proportion, grouped by condition
# ============================================================================
#' @keywords internal
#' @noRd
.plot_composition_box <- function(df, group_col, colors) {
  if (is.null(group_col)) {
    stop("`group_col` is required for style = 'box'.")
  }
  cluster <- prop <- group <- NULL  # NSE
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = group, y = prop, fill = group)) +
    ggplot2::geom_boxplot(color = "black", outlier.size = 1) +
    ggplot2::facet_wrap(~ cluster, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Proportion", fill = "Group") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }
  p
}


# ============================================================================
# Style: mean +/- SE per group per cluster
# ============================================================================
#' @keywords internal
#' @noRd
.plot_composition_line <- function(df, group_col, colors) {
  if (is.null(group_col)) {
    stop("`group_col` is required for style = 'line'.")
  }
  group <- prop <- cluster <- NULL  # NSE
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = group, y = prop, color = cluster,
                                    group = cluster)) +
    ggplot2::stat_summary(fun.data = "mean_se", geom = "pointrange") +
    ggplot2::stat_summary(fun = "mean", geom = "line") +
    ggplot2::labs(x = NULL, y = "Proportion", color = "Cluster") +
    ggplot2::theme_bw()
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(values = colors)
  }
  p
}
