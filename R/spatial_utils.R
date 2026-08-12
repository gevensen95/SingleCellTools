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


# ============================================================================
# .spatial_panel_base() / .combine_spatial_plots() -- shared by
# SpatialDimPlotFixed() and SpatialFeaturePlotFixed(). Both exist because
# Seurat's own SpatialDimPlot()/SpatialFeaturePlot() draw everything (tissue
# image + spots) inside a single custom GeomSpatial layer whose point-size
# calculation has had real, documented regressions across Seurat 5.x
# releases (pt.size.factor silently doing nothing -- see
# https://github.com/satijalab/seurat/issues/9491,
# https://github.com/satijalab/seurat/issues/6179,
# https://github.com/satijalab/seurat/issues/4272). These build the
# equivalent plot from ordinary ggplot2 layers instead -- a real
# geom_point(size = ...) that's guaranteed to respond to its size argument
# -- so they don't depend on that code path at all.
# ============================================================================

#' @keywords internal
#' @noRd
.spatial_panel_base <- function(obj, image_name, crop, image.alpha, pad,
                                cells = NULL) {
  if (!image_name %in% names(obj@images)) {
    stop("Image '", image_name, "' not found in obj@images. Available: ",
         paste(names(obj@images), collapse = ", "))
  }
  img_arr <- obj@images[[image_name]]@image
  h <- dim(img_arr)[1]
  w <- dim(img_arr)[2]

  coords <- .get_fov_coords(obj, image_name, cells = cells)
  if (nrow(coords) == 0) {
    stop("No cells with coordinates found for image '", image_name, "'.")
  }
  coords <- as.data.frame(coords)
  coords$cell <- rownames(coords)
  # GetTissueCoordinates()'s row coordinate increases downward (standard
  # image pixel convention); negate it so "up" agrees between the points and
  # the raster below, rather than relying on scale_y_reverse() interacting
  # correctly with annotation_raster() (a real ggplot2/grob-level ambiguity
  # this sidesteps entirely by just doing the flip with arithmetic).
  coords$y <- -coords$y

  if (isTRUE(crop)) {
    xr <- range(coords$x)
    yr <- range(coords$y)
    xpad <- max(diff(xr) * pad, 1)
    ypad <- max(diff(yr) * pad, 1)
    xlim <- xr + c(-xpad, xpad)
    ylim <- yr + c(-ypad, ypad)
  } else {
    xlim <- c(0, w)
    ylim <- c(-h, 0)
  }

  p <- ggplot2::ggplot() +
    ggplot2::annotation_raster(
      grDevices::as.raster(img_arr),
      xmin = 0, xmax = w, ymin = -h, ymax = 0,
      interpolate = TRUE
    )

  if (image.alpha < 1) {
    # annotation_raster() has no alpha argument -- fake dimming by overlaying
    # a semi-transparent white rectangle over the whole image extent.
    p <- p + ggplot2::annotate(
      "rect", xmin = 0, xmax = w, ymin = -h, ymax = 0,
      fill = "white", alpha = 1 - image.alpha
    )
  }

  p <- p +
    ggplot2::coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::theme_void() +
    Ol_Reliable() +
    ggplot2::theme(
      panel.border     = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.grid       = ggplot2::element_blank(),
      axis.text        = ggplot2::element_blank(),
      axis.ticks       = ggplot2::element_blank(),
      axis.title       = ggplot2::element_blank()
    )

  list(plot = p, coords = coords)
}

#' @keywords internal
#' @noRd
.combine_spatial_plots <- function(plots, ncol, combine) {
  if (length(plots) == 1) return(plots[[1]])
  if (!isTRUE(combine)) return(plots)
  if (is.null(ncol)) ncol <- min(length(plots), 4)
  patchwork::wrap_plots(plots, ncol = ncol, guides = "collect")
}
