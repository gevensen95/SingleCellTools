#' Flexible polygon plot for spatial segmentation boundaries
#'
#' Draws per-cell segmentation polygons (from [get_polygon_coords()]) as
#' filled shapes, colored by any continuous or discrete column -- a more
#' flexible alternative to `Seurat::ImageFeaturePlot()` /
#' `Seurat::ImageDimPlot()` for building custom or layered spatial figures.
#' Returns a plain `ggplot` object, so it's fully chainable with any
#' `ggplot2` scale/theme/title call afterward.
#'
#' @param poly A data.frame from [get_polygon_coords()] (or anything with the
#'   same `cell`/`x`/`y` columns) -- optionally subset to a particular set of
#'   cells first (e.g. `subset(poly, cell_type == "T cell")`) to build one
#'   layer of a multi-layer overlay with [stack_polygons()].
#' @param feature Optional column in `poly` to color polygons by. `NULL`
#'   (default) draws every polygon in a single flat `background` color
#'   instead -- useful as a base/context layer underneath colored layers in
#'   an overlay.
#' @param type One of `"auto"` (default; numeric columns get a continuous
#'   gradient, everything else gets a discrete palette), `"continuous"`, or
#'   `"discrete"`.
#' @param colors For `type = "continuous"`: a vector of colors passed to
#'   `ggplot2::scale_fill_gradientn(colors = ...)`. For `type = "discrete"`:
#'   a named or unnamed vector passed to `ggplot2::scale_fill_manual(values =
#'   ...)`. `NULL` (default) uses each scale's own default palette.
#' @param fill_limits Optional length-2 numeric vector fixing the color scale
#'   range for `type = "continuous"` (passed to `scale_fill_gradientn(limits
#'   = ...)`). Useful when comparing the same feature across several
#'   separately-built plots (e.g. stacked layers, or different samples) that
#'   should all share one color scale rather than each auto-scaling to its
#'   own subset's range.
#' @param background Fill color used when `feature = NULL`. Default
#'   `"grey80"`.
#' @param border_color Polygon outline color, passed to
#'   `ggplot2::geom_polygon(color = ...)`. Default `"black"`; set `NA` for no
#'   outline (recommended at high cell density, where borders can dominate
#'   the plot).
#' @param legend_name Legend title. `NULL` (default) uses `feature`.
#' @param legend_position Passed to `ggplot2::theme(legend.position = ...)`.
#'   Default `"right"`.
#'
#' @return A `ggplot` object.
#'
#' @seealso [get_polygon_coords()], [stack_polygons()], [collect_legend()]
#'
#' @examples
#' \dontrun{
#' poly <- get_polygon_coords(xenium, image = "fov1", meta_cols = "cell_type")
#' PlotPolygons(poly, feature = "cell_type")
#'
#' # Chainable like any ggplot
#' PlotPolygons(poly, feature = "cell_type") +
#'   ggplot2::ggtitle("Cell types in fov1")
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_polygon scale_fill_gradientn scale_fill_manual coord_fixed labs theme
#' @export
PlotPolygons <- function(poly,
                         feature          = NULL,
                         type             = c("auto", "continuous", "discrete"),
                         colors           = NULL,
                         fill_limits      = NULL,
                         background       = "grey80",
                         border_color     = "black",
                         legend_name      = NULL,
                         legend_position  = "right") {
  type <- match.arg(type)
  if (!is.data.frame(poly)) {
    stop("`poly` must be a data.frame (e.g. from get_polygon_coords()).")
  }
  if (!all(c("cell", "x", "y") %in% colnames(poly))) {
    stop("`poly` must have columns 'cell', 'x', 'y'.")
  }

  x <- y <- cell <- NULL  # NSE

  if (is.null(feature)) {
    p <- ggplot2::ggplot(poly, ggplot2::aes(x = x, y = y, group = cell)) +
      ggplot2::geom_polygon(fill = background, color = border_color)
  } else {
    if (!feature %in% colnames(poly)) {
      stop("Column '", feature, "' not found in `poly`.")
    }
    vals <- poly[[feature]]
    if (type == "auto") type <- if (is.numeric(vals)) "continuous" else "discrete"
    used_name <- if (is.null(legend_name)) feature else legend_name

    p <- ggplot2::ggplot(poly, ggplot2::aes(x = x, y = y, group = cell,
                                            fill = .data[[feature]])) +
      ggplot2::geom_polygon(color = border_color)

    if (type == "continuous") {
      grad_colors <- if (is.null(colors)) {
        c("lightgrey", "#FEB24C", "#F03B20", "#BD0026")
      } else colors
      p <- p + ggplot2::scale_fill_gradientn(colors = grad_colors,
                                             limits = fill_limits,
                                             name   = used_name)
    } else {
      p <- p + ggplot2::scale_fill_manual(values = colors, name = used_name)
    }
  }

  p +
    ggplot2::coord_fixed() +
    Ol_Reliable() +
    ggplot2::theme(legend.position = legend_position) +
    ggplot2::labs(x = NULL, y = NULL)
}
