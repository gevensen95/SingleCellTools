# Internal helpers shared across spatial functions. Not exported.
#
# .assert_seurat() / .as_seurat_list() dedupe the "is this actually a
# Seurat object (or list of them)" checks that were previously copy-pasted
# near-identically across SpatialObjectInfo, DropSpatialImage,
# GetHiresVisiumImage, QueryXeniumMolecules, NeighborhoodEnrichment,
# detect_fov_edges, AnnotateRegions, and get_cells_in_polygon.
#
# .get_fov_coords() dedupes the "normalize Seurat::GetTissueCoordinates()
# output into an x/y matrix keyed by cell ID" logic that was previously
# copy-pasted (with slight, silently-drifting variations -- see below) into
# NeighborhoodEnrichment, detect_fov_edges, and BuildMultipleNicheAssays.
# Column names in GetTissueCoordinates() output vary by Seurat/image-class
# version: `imagecol`/`imagerow` (older Seurat, VisiumV1 image class, pixel
# coordinates, rownames = barcodes) vs. `x`/`y` (current Seurat/SeuratData
# Visium image classes and all FOV-based platforms -- Xenium/CosMx/MERFISH
# -- with a `cell` column holding barcodes instead of rownames).
# BuildMultipleNicheAssays already handled both cases; NeighborhoodEnrichment
# and detect_fov_edges only handled the x/y case, so consolidating onto
# BuildMultipleNicheAssays's more complete logic also fixes a latent bug in
# the other two (they'd have hard-errored on an object using the older
# imagecol/imagerow layout).

#' @keywords internal
#' @noRd
.assert_seurat <- function(obj, arg_name = "obj") {
  if (!inherits(obj, "Seurat")) {
    stop("`", arg_name, "` must be a Seurat object.")
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.as_seurat_list <- function(obj, arg_name = "obj") {
  was_single <- inherits(obj, "Seurat")
  objs <- if (was_single) list(obj) else obj
  if (!is.list(objs) || length(objs) == 0L ||
      !all(vapply(objs, inherits, logical(1), "Seurat"))) {
    stop("`", arg_name, "` must be a Seurat object, or a list of Seurat objects.")
  }
  list(objs = objs, was_single = was_single)
}

#' @keywords internal
#' @noRd
.get_fov_coords <- function(obj, fov, cells = NULL, which = "centroids") {
  coords_df <- as.data.frame(
    Seurat::GetTissueCoordinates(obj[[fov]], which = which)
  )
  if (all(c("imagecol", "imagerow") %in% colnames(coords_df))) {
    coord_cols <- c("imagecol", "imagerow")
  } else if (all(c("x", "y") %in% colnames(coords_df))) {
    coord_cols <- c("x", "y")
  } else {
    stop("Could not find coordinate columns in GetTissueCoordinates() ",
         "output for FOV '", fov, "' (got: ",
         paste(colnames(coords_df), collapse = ", "), ").")
  }
  if ("cell" %in% colnames(coords_df)) {
    rownames(coords_df) <- coords_df[["cell"]]
  }
  coords_mat <- as.matrix(coords_df[, coord_cols, drop = FALSE])
  colnames(coords_mat) <- c("x", "y")
  if (!is.null(cells)) {
    keep <- intersect(rownames(coords_mat), cells)
    coords_mat <- coords_mat[keep, , drop = FALSE]
  }
  coords_mat
}
