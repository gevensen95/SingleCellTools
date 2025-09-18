#' Get Cells in Polygon
#'
#' This function will take a polygon data frame (x & y coordinates) and parse out which cells are in that polygon.
#'
#' @param seurat.obj Seurat object with spatial data
#' @param poly_df Data frame of polygon coordinates (x & y as columns)
#' @param image_name Name of image to look for cells in
#' @return List of data frames
#' @export
#'
get_cells_in_polygon <- function(seurat.obj, poly_df, image_name) {
  # deps
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required. Please install it.")
  }
  if (!inherits(seurat.obj, "Seurat")) {
    stop("`seurat.obj` must be a Seurat object.")
  }
  if (!all(c("x","y") %in% names(poly_df))) {
    stop("`poly_df` must be a data frame with columns named 'x' and 'y'.")
  }

  # ---- 1) Polygon (ensure closed ring) ----
  polygon_matrix <- as.matrix(poly_df)
  if (!all(polygon_matrix[1, ] == polygon_matrix[nrow(polygon_matrix), ])) {
    polygon_matrix <- rbind(polygon_matrix, polygon_matrix[1, ])
  }
  poly_sf <- sf::st_sfc(sf::st_polygon(list(polygon_matrix)))

  # ---- 2) Cell coordinates from Seurat ----
  coords <- Seurat::GetTissueCoordinates(seurat.obj, image = image_name)
  coords <- as.data.frame(coords)

  # Standardize to x/y
  if (all(c("imagecol", "imagerow") %in% names(coords))) {
    coords <- dplyr::rename(coords, x = imagecol, y = imagerow)
  } else if (!all(c("x", "y") %in% names(coords))) {
    stop("Couldn't find 'x'/'y' or 'imagecol'/'imagerow' in tissue coordinates.")
  }

  # Points → sf
  cells_sf <- sf::st_as_sf(coords, coords = c("x", "y"), crs = NA)

  # ---- 3) Point-in-polygon ----
  inside <- sf::st_within(cells_sf, poly_sf, sparse = FALSE)[, 1]

  # Map back to Seurat metadata (logical vector aligned to Cells(seurat.obj))
  inside_named <- setNames(inside, cells_sf$cell)

  # Helpful return values
  cells_inside <- inside_named[inside_named==TRUE]
  return(cells_inside)
}
