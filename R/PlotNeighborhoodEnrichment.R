#' Diverging heatmap for NeighborhoodEnrichment() results
#'
#' Cell-type x cell-type diverging heatmap of \code{\link{NeighborhoodEnrichment}}'s
#' z-score matrix, with significance stars drawn from \code{padj} -- the
#' same visual convention as \code{\link{plotNicheCoExpress}}. Red =
#' neighboring more than the permutation null predicts, blue = less.
#'
#' @param enrich The list returned by \code{\link{NeighborhoodEnrichment}}.
#' @param sig_levels Named thresholds for significance stars. Default
#'   \code{c(`***` = 0.001, `**` = 0.01, `*` = 0.05)}.
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' enrich <- NeighborhoodEnrichment(obj, group.by = "cell_type", k = 10)
#' PlotNeighborhoodEnrichment(enrich)
#' }
#' @seealso \code{\link{NeighborhoodEnrichment}}
#' @importFrom ggplot2 aes element_text geom_text geom_tile ggplot labs scale_fill_gradient2 theme
#' @export
PlotNeighborhoodEnrichment <- function(enrich,
                                       sig_levels = c(`***` = 0.001,
                                                      `**`  = 0.01,
                                                      `*`   = 0.05)) {

  if (!is.list(enrich) || is.null(enrich$results)) {
    stop("`enrich` must be the list returned by NeighborhoodEnrichment() ",
         "(missing `results`).")
  }
  df <- enrich$results

  star <- function(p) {
    vapply(p, function(pp) {
      if (is.na(pp)) return("")
      hit <- names(sig_levels)[pp <= sig_levels]
      if (length(hit)) hit[which.min(sig_levels[hit])] else ""
    }, character(1))
  }
  df$label <- star(df$padj)

  lim <- max(abs(df$z), na.rm = TRUE)
  focal <- neighbor <- z <- label <- NULL  # NSE silencing
  ggplot2::ggplot(df, ggplot2::aes(x = neighbor, y = focal, fill = z)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 4,
                       vjust = 0.78, color = "black") +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, limits = c(-lim, lim), name = "z-score") +
    ggplot2::labs(x = "neighbor cell type", y = "focal cell type",
                 title = "Spatial neighborhood enrichment",
                 subtitle = "red = more neighboring than chance  |  * padj<.05  ** <.01  *** <.001") +
    Ol_Reliable() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
