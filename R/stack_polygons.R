#' Stack multiple PlotPolygons() layers into one overlay
#'
#' Prepares a [PlotPolygons()] plot to be layered into a single composite
#' figure with other polygon plots -- pins it to a shared coordinate range
#' (from `poly`, typically the full, unsubset data frame every layer's cells
#' were drawn from) and strips its legend/titles so stacking doesn't shift
#' the panel. The first layer (`first = TRUE`) anchors the composite as an
#' ordinary `ggplot`; every later layer (`first = FALSE`) is wrapped as a
#' transparent `patchwork` inset that can be added directly on top of it
#' with `+`.
#'
#' Legends are dropped here because each layer typically has its own,
#' independent color scale (e.g. one cell type on a red gradient, another on
#' blue) -- stacking full plots as transparent insets, rather than merging
#' them into one `ggplot` object, is what lets each layer keep its own scale
#' without the layers fighting over a single shared legend/scale. Recover
#' the dropped legends afterward with [collect_legend()] on the *original*
#' (pre-`stack_polygons()`) plots.
#'
#' Unlike a plain `ggplot2::coord_cartesian(xlim =, ylim =)` (which does not
#' preserve aspect ratio), this reapplies `coord_fixed(xlim =, ylim =)` so
#' every stacked layer keeps the same 1:1 spatial aspect ratio
#' [PlotPolygons()] sets by default -- important for segmentation polygons,
#' where a stretched aspect ratio would visibly distort cell shapes.
#'
#' @param plot A `ggplot` object, typically from [PlotPolygons()].
#' @param poly The data.frame the shared coordinate range should come from --
#'   pass the same (unsubset) data.frame for every layer being stacked
#'   together, even though each layer's own plot may have been built from a
#'   subset of it, so all layers end up on identical axis limits.
#' @param first Logical; `TRUE` for exactly one layer -- the one at the
#'   bottom of the stack. Default `FALSE`.
#'
#' @return For `first = TRUE`, a `ggplot` object. For `first = FALSE`, a
#'   `patchwork` inset element -- add it to the `first = TRUE` plot (or a
#'   previous inset) with `+`.
#'
#' @seealso [PlotPolygons()], [collect_legend()]
#'
#' @examples
#' \dontrun{
#' poly <- get_polygon_coords(xenium, "fov1", meta_cols = "cell_type")
#'
#' bg    <- PlotPolygons(poly, background = "grey90")
#' tcell <- PlotPolygons(subset(poly, cell_type == "T cell"),
#'                       feature = "Cd3e", colors = c("white", "red"))
#' bcell <- PlotPolygons(subset(poly, cell_type == "B cell"),
#'                       feature = "Cd19", colors = c("white", "blue"))
#'
#' overlay <- stack_polygons(bg, poly, first = TRUE) +
#'   stack_polygons(tcell, poly) +
#'   stack_polygons(bcell, poly)
#'
#' legends <- patchwork::wrap_plots(collect_legend(tcell), collect_legend(bcell), ncol = 1)
#' overlay + legends + patchwork::plot_layout(widths = c(1, 0.25))
#' }
#'
#' @importFrom ggplot2 coord_fixed theme element_blank
#' @export
stack_polygons <- function(plot, poly, first = FALSE) {
  if (!inherits(plot, "ggplot")) {
    stop("`plot` must be a ggplot object (e.g. from PlotPolygons()).")
  }
  if (!is.data.frame(poly) || !all(c("x", "y") %in% colnames(poly))) {
    stop("`poly` must be a data.frame with 'x' and 'y' columns (e.g. from ",
        "get_polygon_coords()).")
  }

  xlim <- range(poly$x, na.rm = TRUE)
  ylim <- range(poly$y, na.rm = TRUE)

  prepped <- plot +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim, expand = TRUE) +
    ggplot2::theme(legend.position = "none",
                  plot.title    = ggplot2::element_blank(),
                  plot.subtitle = ggplot2::element_blank(),
                  plot.caption  = ggplot2::element_blank())

  if (isTRUE(first)) {
    prepped
  } else {
    patchwork::inset_element(prepped, left = 0, right = 1, top = 1, bottom = 0,
                             align_to = "panel")
  }
}
