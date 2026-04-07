#' Get coordinates from all FOVs
#'
#' This function collect the coordinates from each FOV and combine them into a single data frame
#'
#' @param obj Seurat object (spatial)
#' @return a data frame
#' @export
get_all_coords <- function(seurat_obj, meta_cols = NULL) {
  fov_names <- Images(seurat_obj)

  coords_list <- lapply(fov_names, function(fov) {
    coords <- GetTissueCoordinates(seurat_obj, image = fov)
    coords$fov <- fov

    # Optionally join selected metadata columns
    if (!is.null(meta_cols)) {
      meta <- seurat_obj@meta.data[rownames(coords), meta_cols, drop = FALSE]
      coords <- cbind(coords, meta)
    }

    coords
  })

  do.call(rbind, coords_list)
}
