#' Assign cells to anatomical regions defined by polygons
#'
#' Given a named list of polygon data frames (each with \code{x} and \code{y}
#' columns) and a Seurat object with spatial coordinates, labels every cell
#' with the name of the polygon it falls inside (or \code{"unassigned"} if it
#' falls outside all of them). Pairs naturally with \code{get_cells_in_polygon}
#' and \code{parse_polygons}.
#'
#' If polygons overlap, the cell is assigned to whichever polygon was checked
#' first (the order of \code{polygons}).
#'
#' @param obj A Seurat object with spatial data.
#' @param polygons Named list of polygon data frames. Each must have columns
#'   \code{x} and \code{y}. The list names become the region labels.
#' @param image_name Image / FOV to use for cell coordinates.
#' @param region_col Name of the new metadata column. Default \code{"region"}.
#' @param unassigned_label Label for cells not in any polygon. Default
#'   \code{"unassigned"}.
#' @return The Seurat object with a new metadata column.
#' @importFrom Seurat GetTissueCoordinates
#' @export
AnnotateRegions <- function(obj,
                            polygons,
                            image_name,
                            region_col       = "region",
                            unassigned_label = "unassigned") {

  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (!is.list(polygons) || is.null(names(polygons))) {
    stop("`polygons` must be a NAMED list of polygon data frames.")
  }
  for (nm in names(polygons)) {
    if (!all(c("x", "y") %in% names(polygons[[nm]]))) {
      stop("Polygon '", nm, "' must have columns 'x' and 'y'.")
    }
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Install with: install.packages('sf')")
  }

  message(sprintf("--- Building %d polygon(s) ---", length(polygons)))
  poly_sfs <- lapply(polygons, function(p) {
    pm <- as.matrix(p[, c("x", "y")])
    if (!all(pm[1, ] == pm[nrow(pm), ])) pm <- rbind(pm, pm[1, ])
    sf::st_sfc(sf::st_polygon(list(pm)))
  })

  message(sprintf("--- Pulling tissue coordinates (image: %s) ---", image_name))
  coords <- Seurat::GetTissueCoordinates(obj, image = image_name)
  coords <- as.data.frame(coords)
  if (all(c("imagecol", "imagerow") %in% names(coords))) {
    coords$x <- coords$imagecol
    coords$y <- coords$imagerow
  } else if (!all(c("x", "y") %in% names(coords))) {
    stop("Couldn't find x/y or imagecol/imagerow in tissue coordinates.")
  }
  if (!"cell" %in% colnames(coords)) {
    coords$cell <- rownames(coords)
  }

  cells_sf <- sf::st_as_sf(coords, coords = c("x", "y"), crs = NA)

  # Assign by checking polygons in order; first hit wins.
  assignment <- rep(unassigned_label, nrow(cells_sf))
  for (region_name in names(poly_sfs)) {
    inside <- sf::st_within(cells_sf, poly_sfs[[region_name]], sparse = FALSE)[, 1]
    # Only fill unassigned cells (preserve first-hit behavior)
    keep <- inside & assignment == unassigned_label
    assignment[keep] <- region_name
    message(sprintf("  %-20s %d cells", region_name, sum(keep)))
  }
  message(sprintf("  %-20s %d cells", unassigned_label,
                  sum(assignment == unassigned_label)))

  # Push back to obj metadata, aligning by cell barcode. `coords$cell` only
  # covers the cells in `image_name` -- when the object has other
  # FOVs/images too, those cells aren't in `to_add`'s names at all.
  # Subsetting a named vector by a key it doesn't contain (to_add[colnames(obj)])
  # returns NA for *both* the value and the name at that position (not the
  # queried cell ID), which then breaks `[[<-`'s internal name-matching
  # check (comparing a real colname against an NA name yields NA, not
  # TRUE/FALSE, and `if (NA)` errors with "missing value where TRUE/FALSE
  # needed"). Build a full-length, fully-named vector explicitly instead so
  # every cell -- covered or not -- gets a real, non-NA name.
  to_add <- setNames(assignment, coords$cell)
  region_vals <- setNames(rep(unassigned_label, ncol(obj)), colnames(obj))
  common <- intersect(names(to_add), names(region_vals))
  region_vals[common] <- to_add[common]
  obj[[region_col]] <- unname(region_vals)
  obj
}
