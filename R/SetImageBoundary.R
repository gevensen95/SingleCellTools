#' Set Boundary of Images
#'
#' This function changes the default boundary of a spatial Seurat object with multiple images.
#'
#' @param seurat_obj Seurat object (spatial)
#' @param boundary Boundary name (e.g. "segmentation" or "centroids")
#' @param verbose Message progress. Default \code{TRUE}.
#' @return a Seurat object
#' @examples
#' \dontrun{
#' xenium <- SetImageBoundary(xenium, "segmentation")
#' }
#' @export
SetImageBoundary <- function(seurat_obj, boundary, verbose = TRUE) {
  .assert_seurat(seurat_obj, "seurat_obj")
  if (!is.character(boundary) || length(boundary) != 1 || !nzchar(boundary)) {
    stop("`boundary` must be a single non-empty string (e.g. 'segmentation' or 'centroids').")
  }

  fovs <- Seurat::Images(seurat_obj)
  if (length(fovs) == 0) {
    warning("`seurat_obj` has no images in seurat_obj@images -- nothing to set.")
    return(seurat_obj)
  }
  if (isTRUE(verbose)) {
    message(sprintf('--- Setting default boundary to "%s" across %d FOVs ---',
                    boundary, length(fovs)))
  }

  for (fov in fovs) {
    img <- seurat_obj@images[[fov]]
    # Boundaries() is exported by SeuratObject, not Seurat (confirmed:
    # Seurat::Boundaries() errors with "not an exported object").
    available <- SeuratObject::Boundaries(img)

    if (!boundary %in% available) {
      warning(sprintf(
        "Boundary '%s' not found in FOV '%s'. Available: %s",
        boundary, fov, paste(available, collapse = ", ")
      ))
      next
    }

    Seurat::DefaultBoundary(img) <- boundary
    seurat_obj@images[[fov]] <- img
  }

  return(seurat_obj)
}
