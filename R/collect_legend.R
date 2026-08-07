#' Extract a plot's legend as a standalone grob
#'
#' Pulls the legend/guide-box out of a `ggplot` object so it can be
#' recombined separately from the plot it came from -- for example, to
#' collect the (per-layer) legends that [stack_polygons()] strips out while
#' building a multi-layer overlay, and lay them out alongside the combined
#' plot instead of losing them.
#'
#' Implemented with base `ggplot2`/`grid`/`gtable` machinery (all already
#' guaranteed available anywhere `ggplot2` is), rather than pulling in
#' `cowplot` or `ggpubr` as an extra dependency just for this.
#'
#' @param plot A `ggplot` object with a visible legend (i.e. not built with
#'   `legend.position = "none"`, and with at least one mapped aesthetic that
#'   produces a guide).
#'
#' @return A `patchwork`-composable object (from `patchwork::wrap_elements()`)
#'   wrapping the legend grob -- combine it with other plots/legends using
#'   `patchwork`'s `+` / `plot_layout()`.
#'
#' @seealso [stack_polygons()], [PlotPolygons()]
#'
#' @examples
#' \dontrun{
#' tcell <- PlotPolygons(subset(poly, cell_type == "T cell"), feature = "Cd3e")
#' bcell <- PlotPolygons(subset(poly, cell_type == "B cell"), feature = "Cd19")
#' legends <- patchwork::wrap_plots(collect_legend(tcell), collect_legend(bcell), ncol = 1)
#' }
#'
#' @export
collect_legend <- function(plot) {
  if (!inherits(plot, "ggplot")) {
    stop("`plot` must be a ggplot object.")
  }

  gt         <- ggplot2::ggplotGrob(plot)
  grob_names <- vapply(gt$grobs, function(g) g$name, character(1))
  guide_idx  <- grep("guide-box", grob_names)

  if (length(guide_idx) == 0) {
    stop("`plot` has no legend to collect -- check that legend.position ",
        "isn't 'none' and that at least one aesthetic (e.g. fill/color) is ",
        "mapped.")
  }

  patchwork::wrap_elements(gt$grobs[[guide_idx[1]]])
}
