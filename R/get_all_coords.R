#' Get coordinates from all FOVs
#'
#' This function collect the coordinates from each FOV and combine them into a single data frame
#'
#' @param obj Seurat object (spatial)
#' @return a data frame
#' @export
get_all_coords <- function(seurat_obj, meta_cols = NULL) {
  fov_names <- Images(seurat_obj)
  message(sprintf('--- Collecting tissue coordinates from %d FOVs ---', length(fov_names)))

  coords_list <- lapply(fov_names, function(fov) {
    coords <- GetTissueCoordinates(seurat_obj, image = fov)
    coords$fov <- fov

    # Optionally join selected metadata columns. GetTissueCoordinates() on
    # an FOV/centroids image returns cell identity in a "cell" column, not
    # as rownames (rownames are typically just 1, 2, 3, ...), so indexing
    # meta.data by rownames(coords) silently looks up the wrong keys and
    # returns all-NA columns. Use the "cell" column when present, matching
    # the same fallback AnnotateRegions() already uses for this.
    if (!is.null(meta_cols)) {
      cell_ids <- if ("cell" %in% colnames(coords)) coords$cell else rownames(coords)
      meta <- seurat_obj@meta.data[cell_ids, meta_cols, drop = FALSE]
      coords <- cbind(coords, meta)
    }

    coords
  })

  message('--- Combining coordinates into single data frame ---')
  do.call(rbind, coords_list)
}
