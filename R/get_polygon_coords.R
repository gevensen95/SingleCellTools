#' Extract per-vertex cell segmentation polygon coordinates
#'
#' Pulls the segmentation-boundary polygon vertices for every cell in a
#' spatial FOV/image (Xenium, CosMx, MERFISH -- any Seurat object whose
#' image has a `"segmentation"` boundary set, as opposed to only cell
#' centroid points), optionally joined with metadata columns for downstream
#' coloring. Companion to [PlotPolygons()], which plots the result.
#'
#' `Seurat::GetTissueCoordinates(obj[[image]], which = "segmentation")`
#' returns one row per polygon VERTEX (many rows per cell), unlike the
#' one-row-per-cell centroid coordinates [get_all_coords()] collects -- so
#' this is deliberately a separate function rather than a
#' `which = "segmentation"` option bolted onto that one. Joining metadata by
#' cell also has to replicate each cell's metadata across all of its vertex
#' rows, which `get_all_coords()`'s row-per-cell logic doesn't do.
#'
#' @param obj A Seurat object.
#' @param image Name of the image/FOV in `obj@images` to extract from.
#' @param meta_cols Optional character vector of `obj@meta.data` columns to
#'   join onto the output (e.g. a cell-type column, or a gene's expression
#'   pulled in beforehand with `obj$gene <- Seurat::FetchData(obj,
#'   "gene")[[1]]`). `NULL` (default) returns just the coordinates.
#'
#' @return A data.frame with columns `cell`, `x`, `y` (one row per polygon
#'   vertex, in vertex order), plus any `meta_cols` requested.
#'
#' @seealso [PlotPolygons()], [stack_polygons()]
#'
#' @examples
#' \dontrun{
#' poly <- get_polygon_coords(xenium, image = "fov1", meta_cols = "cell_type")
#' PlotPolygons(poly, feature = "cell_type")
#' }
#'
#' @export
get_polygon_coords <- function(obj, image, meta_cols = NULL) {
  .assert_seurat(obj)
  if (!image %in% Seurat::Images(obj)) {
    stop("Image '", image, "' not found in obj@images. Available: ",
        paste(Seurat::Images(obj), collapse = ", "))
  }

  img       <- obj@images[[image]]
  # Boundaries() is exported by SeuratObject, not Seurat (confirmed:
  # Seurat::Boundaries() errors with "not an exported object").
  available <- SeuratObject::Boundaries(img)
  if (!"segmentation" %in% available) {
    stop("Image '", image, "' has no 'segmentation' boundary set ",
        "(available: ", paste(available, collapse = ", "), "). ",
        "get_polygon_coords() needs actual cell-boundary polygons, not ",
        "centroid points -- use a loader that stores segmentation ",
        "boundaries for this platform (e.g. LoadXenium2(..., type = ",
        "c('centroids', 'segmentations'))).")
  }

  message(sprintf('--- Extracting segmentation polygons from image "%s" ---', image))
  poly <- as.data.frame(Seurat::GetTissueCoordinates(img, which = "segmentation"))

  if (!"cell" %in% colnames(poly)) {
    stop("GetTissueCoordinates(which = 'segmentation') did not return a ",
        "'cell' column -- got: ", paste(colnames(poly), collapse = ", "))
  }
  if (!all(c("x", "y") %in% colnames(poly))) {
    stop("GetTissueCoordinates(which = 'segmentation') did not return ",
        "'x'/'y' columns -- got: ", paste(colnames(poly), collapse = ", "))
  }
  poly <- poly[, c("cell", "x", "y")]

  if (!is.null(meta_cols)) {
    missing_cols <- setdiff(meta_cols, colnames(obj@meta.data))
    if (length(missing_cols) > 0) {
      stop("Column(s) not found in obj@meta.data: ",
          paste(missing_cols, collapse = ", "))
    }
    meta <- obj@meta.data[poly$cell, meta_cols, drop = FALSE]
    poly <- cbind(poly, meta)
  }

  poly
}
