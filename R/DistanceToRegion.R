#' Distance from every cell to the nearest of a set of polygon regions
#'
#' Generalizes the polygon-distance machinery that
#' \code{\link{AnnotateRegions}} (inside/outside a polygon),
#' \code{\link{detect_fov_edges}} (distance to the FOV's own outer
#' boundary), and \code{detect_tissue_holes} (distance to an internal gap)
#' each compute in their own narrower way, into one reusable "distance from
#' each cell to the nearest boundary of a named region" primitive. Where
#' \code{AnnotateRegions()} answers "which region is this cell in" (a
#' categorical call), this answers "how far is this cell from region X" (a
#' continuous gradient) -- useful whenever the biology is defined by
#' physical distance from a structure rather than a hard inside/outside
#' line, e.g. distance from a portal tract or central vein as a check on a
#' marker-score-based zonation call (see \code{\link{ClassifyByReferenceCutoff}}),
#' or distance from a tumor margin, vessel, or any other hand-drawn/
#' segmented boundary.
#'
#' @param obj A Seurat object with spatial data.
#' @param polygons Named list of polygon data frames, same convention as
#'   \code{\link{AnnotateRegions}}: each needs \code{x}/\code{y} columns,
#'   and the list names become the region labels recorded in
#'   \code{nearest_region_col}.
#' @param image_name Image / FOV to use for cell coordinates (same
#'   coordinate space the polygons were drawn in).
#' @param distance_col Name of the new numeric metadata column holding each
#'   cell's distance to its nearest region. Units match your coordinate
#'   space (pixels, microns, ...), not a fixed physical unit. Default
#'   \code{"distance_to_region"}.
#' @param nearest_region_col Optional name of an additional metadata column
#'   recording which region was nearest for each cell. \code{NULL}
#'   (default) skips this.
#' @param to_boundary If \code{FALSE} (default), distance is computed to
#'   each polygon as a filled area, so a cell already inside a region gets
#'   distance \code{0} to it. If \code{TRUE}, distance is computed to the
#'   polygon's boundary line instead, so cells inside a region also get a
#'   positive "how far from the edge" distance (useful for e.g. "distance
#'   from the nearest vessel wall, whether you're inside or outside it").
#' @param cells Optional character vector of cell names to restrict to.
#'   \code{NULL} (default) uses every cell with coordinates in
#'   \code{image_name}; cells not covered get \code{NA}.
#' @return \code{obj} with \code{distance_col} (and \code{nearest_region_col},
#'   if requested) added.
#' @examples
#' \dontrun{
#' portal_polys <- parse_polygons(portal_tract_coords)
#' xenium <- DistanceToRegion(xenium, portal_polys, image_name = "fov1",
#'                            distance_col = "dist_to_portal",
#'                            to_boundary = TRUE)
#' }
#' @importFrom Seurat GetTissueCoordinates
#' @export
DistanceToRegion <- function(obj,
                             polygons,
                             image_name,
                             distance_col       = "distance_to_region",
                             nearest_region_col = NULL,
                             to_boundary        = FALSE,
                             cells              = NULL) {

  .assert_seurat(obj)
  if (!is.list(polygons) || is.null(names(polygons)) || any(names(polygons) == "")) {
    stop("`polygons` must be a NAMED list of polygon data frames (see ",
         "?AnnotateRegions), e.g. list(portal = df1, central = df2).")
  }
  for (nm in names(polygons)) {
    if (!all(c("x", "y") %in% names(polygons[[nm]]))) {
      stop("Polygon '", nm, "' must have columns 'x' and 'y'.")
    }
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Install with: install.packages('sf')")
  }

  message(sprintf("--- Building %d polygon region(s)%s ---", length(polygons),
                  if (isTRUE(to_boundary)) " (boundary distance)" else ""))
  poly_geoms <- lapply(polygons, function(p) {
    pm <- as.matrix(p[, c("x", "y")])
    if (!all(pm[1, ] == pm[nrow(pm), ])) pm <- rbind(pm, pm[1, ])
    g <- sf::st_sfc(sf::st_polygon(list(pm)))
    if (isTRUE(to_boundary)) sf::st_boundary(g) else g
  })

  message(sprintf("--- Pulling tissue coordinates (image: %s) ---", image_name))
  coords <- as.data.frame(Seurat::GetTissueCoordinates(obj, image = image_name))
  if (all(c("imagecol", "imagerow") %in% names(coords))) {
    coords$x <- coords$imagecol
    coords$y <- coords$imagerow
  } else if (!all(c("x", "y") %in% names(coords))) {
    stop("Couldn't find x/y or imagecol/imagerow in tissue coordinates.")
  }
  if (!"cell" %in% colnames(coords)) coords$cell <- rownames(coords)
  if (!is.null(cells)) {
    coords <- coords[coords$cell %in% cells, , drop = FALSE]
  }
  if (nrow(coords) == 0) {
    stop("No cells with coordinates found for image '", image_name, "'",
         if (!is.null(cells)) " matching `cells`" else "", ".")
  }

  cells_sf <- sf::st_as_sf(coords, coords = c("x", "y"), crs = NA)

  # Distance to every region, then take the nearest per cell.
  region_names <- names(poly_geoms)
  dist_mat <- vapply(region_names, function(nm) {
    as.numeric(sf::st_distance(cells_sf, poly_geoms[[nm]]))
  }, numeric(nrow(cells_sf)))
  if (length(region_names) == 1) {
    dist_mat <- matrix(dist_mat, ncol = 1, dimnames = list(NULL, region_names))
  }

  nearest_idx    <- apply(dist_mat, 1, which.min)
  nearest_dist   <- dist_mat[cbind(seq_len(nrow(dist_mat)), nearest_idx)]
  nearest_region <- region_names[nearest_idx]

  message(sprintf(
    "  %d cell(s); distance to nearest region ranges %.3g to %.3g (your coordinate units)",
    nrow(dist_mat), min(nearest_dist), max(nearest_dist)))

  # Same full-length, fully-named-vector approach as AnnotateRegions() --
  # subsetting a named vector by a key it doesn't contain returns an NA
  # *name* as well as an NA value, which breaks [[<-'s name matching.
  dist_vals <- setNames(rep(NA_real_, ncol(obj)), colnames(obj))
  to_add    <- setNames(nearest_dist, coords$cell)
  common    <- intersect(names(to_add), names(dist_vals))
  dist_vals[common] <- to_add[common]
  obj[[distance_col]] <- unname(dist_vals)

  if (!is.null(nearest_region_col)) {
    region_vals <- setNames(rep(NA_character_, ncol(obj)), colnames(obj))
    to_add_region <- setNames(nearest_region, coords$cell)
    region_vals[common] <- to_add_region[common]
    obj[[nearest_region_col]] <- unname(region_vals)
  }

  obj
}
