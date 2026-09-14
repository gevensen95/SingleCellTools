#' Build a named color vector: known categories get real colors, the rest a catch-all
#'
#' Common plotting need: highlight a handful of known categories with real
#' colors and gray/black out everything else, without hand-rolling
#' \code{setNames()}/overwrite logic every time (e.g. one cell type of
#' interest against every other cluster, or a small set of named zonation
#' calls against every other cell-type label a classifier left untouched).
#'
#' @param values The categorical values that will actually be plotted --
#'   typically \code{unique(obj$some_column)} or its factor levels. Every
#'   value present here ends up as a name in the result; values not present
#'   here are ignored even if they appear in \code{known}.
#' @param known Named character vector (or named list) mapping category ->
#'   color, e.g. \code{c(pericentral = "#3B82F6", periportal = "#F59E0B")}.
#' @param other Fallback color for every value in \code{values} not named in
#'   \code{known}. Default \code{"black"}.
#' @return A named character vector of colors, one per unique value in
#'   \code{values}, suitable for \code{cols =} in \code{DimPlot()},
#'   \code{SpatialDimPlot()}, \code{\link{SpatialDimPlotFixed}}, etc.
#' @examples
#' \dontrun{
#' zone_colors <- AssignColors(
#'   unique(cosmx$Zone_final),
#'   known = c(pericentral = "#3B82F6", periportal = "#F59E0B",
#'            midlobular = "#10B981", Unclassified = "grey80")
#' )
#' SpatialDimPlotFixed(cosmx, group.by = "Zone_final", cols = zone_colors)
#' }
#' @export
AssignColors <- function(values, known, other = "black") {
  if (is.list(known)) known <- unlist(known)
  if (is.null(names(known)) || any(names(known) == "")) {
    stop("`known` must be a named vector/list, e.g. c(groupA = \"red\").")
  }
  values <- unique(values)
  colors <- setNames(rep(other, length(values)), values)
  hit <- intersect(names(known), names(colors))
  colors[hit] <- known[hit]
  colors
}
