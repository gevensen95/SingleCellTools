#' Get coordinates from all FOVs
#'
#' This function collect the coordinates from each FOV and combine them into a single data frame
#'
#' @param seurat_obj Seurat object (spatial)
#' @param meta_cols Optional character vector of \code{seurat_obj@meta.data}
#'   columns to join onto each row, matched by cell ID. \code{NULL} (default)
#'   doesn't join any metadata.
#' @return a data frame with the original per-FOV coordinate columns
#'   plus canonical `x`/`y` columns (aliased from `imagecol`/`imagerow`
#'   when that's what the FOV's image class provides)
#' @export
get_all_coords <- function(seurat_obj, meta_cols = NULL) {
  fov_names <- Seurat::Images(seurat_obj)
  message(sprintf('--- Collecting tissue coordinates from %d FOVs ---', length(fov_names)))

  coords_list <- lapply(fov_names, function(fov) {
    coords <- Seurat::GetTissueCoordinates(seurat_obj, image = fov)
    coords$fov <- fov

    # Canonicalize to x/y columns regardless of which image class this FOV
    # uses. GetTissueCoordinates() output varies: `imagecol`/`imagerow`
    # (older Seurat, VisiumV1) vs. `x`/`y` (current Visium/FOV image
    # classes) -- the same split `.get_fov_coords()` (spatial_utils.R)
    # already normalizes for NeighborhoodEnrichment()/detect_fov_edges()/
    # BuildMultipleNicheAssays(). A single object whose FOVs mix image
    # classes (e.g. an older Visium round merged with a newer one) would
    # otherwise rbind() into one frame with x/y populated for some rows
    # and NA for others (or vice versa for imagecol/imagerow). Original
    # columns are kept alongside, not replaced, so nothing downstream that
    # already reads imagecol/imagerow directly is affected.
    if (!all(c("x", "y") %in% colnames(coords)) &&
        all(c("imagecol", "imagerow") %in% colnames(coords))) {
      coords$x <- coords$imagecol
      coords$y <- coords$imagerow
    }

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
