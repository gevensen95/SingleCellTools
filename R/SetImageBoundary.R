#' Set Boundary of Images
#'
#' This function changes the default boundary of a spatial Seurat object with multiple images.
#'
#' @param seurat_obj Seurat object (spatial)
#' @param boundary Boundary name (e.g. "segmentation" or "centroids")
#' @return a Seurat object
#' @export
SetImageBoundary <- function(seurat_obj, boundary) {
  fovs <- Seurat::Images(seurat_obj)
  message(sprintf('--- Setting default boundary to "%s" across %d FOVs ---',
                  boundary, length(fovs)))

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
