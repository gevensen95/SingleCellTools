#' Set Boudnary of Images
#'
#' This function changes the default boundary of a spatial Seurat object with multiple images.
#'
#' @param obj Seurat object (spatial)
#' @param boundary Boundry (segmentation or centroids)
#' @return a Seurat object
#' @export
SetImageBoundary <- function(seurat_obj, boundary) {
  fovs <- Images(seurat_obj)
  message(sprintf('--- Setting default boundary to "%s" across %d FOVs ---',
                  boundary, length(fovs)))

  for (fov in fovs) {
    img <- seurat_obj@images[[fov]]
    available <- Boundaries(img)

    if (!boundary %in% available) {
      warning(sprintf(
        "Boundary '%s' not found in FOV '%s'. Available: %s",
        boundary, fov, paste(available, collapse = ", ")
      ))
      next
    }

    DefaultBoundary(img) <- boundary
    seurat_obj@images[[fov]] <- img
  }

  return(seurat_obj)
}
