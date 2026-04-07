#' Set Boudnary of Images
#'
#' This function changes the default boundary of a spatial Seurat object.
#'
#' @param obj Seurat object (spatial)
#' @param boundary Boundry (segmentation or centroids)
#' @return a Seurat object
#' @export
SetImageBoundary <- function(seurat_obj, boundary) {
  for (fov in Images(seurat_obj)) {
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
