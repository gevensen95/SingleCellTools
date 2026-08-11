#' Spatial feature plot, without Seurat's GeomSpatial point-size bug
#'
#' A from-scratch \code{ggplot2} replacement for
#' \code{Seurat::SpatialFeaturePlot()}: plots the tissue image with spots
#' colored by one or more continuous features (genes, module scores, or any
#' other numeric metadata/assay value \code{Seurat::FetchData()} can pull) on
#' top. Built with plain \code{geom_point()} rather than Seurat's custom
#' \code{GeomSpatial} layer -- see \code{\link{SpatialDimPlotFixed}} for why.
#'
#' @param obj A Seurat object with at least one image in \code{obj@images}.
#' @param features Character vector of features to plot (genes, module
#'   scores, numeric metadata columns -- anything \code{Seurat::FetchData()}
#'   accepts). One row of panels per feature when multiple are given and
#'   there's also more than one image.
#' @param images Character vector of image names in \code{obj@images} to
#'   plot. Defaults to all of them.
#' @param assay Assay to pull \code{features} from. Defaults to
#'   \code{Seurat::DefaultAssay(obj)}.
#' @param layer Layer to pull \code{features} from (e.g. \code{"data"},
#'   \code{"counts"}, \code{"scale.data"}). Default \code{"data"}.
#' @param pt.size,pt.alpha,pt.stroke Passed to \code{geom_point()} as
#'   \code{size}, \code{alpha}, \code{stroke}.
#' @param image.alpha Dims the tissue image toward white; \code{1} = full
#'   opacity, \code{0} = fully white. Default \code{1}.
#' @param crop Logical; if \code{TRUE} (default), zooms to the bounding box
#'   of the plotted spots (plus \code{pad}). If \code{FALSE}, shows the full
#'   image canvas.
#' @param pad Fractional padding added around the spot bounding box when
#'   \code{crop = TRUE}. Default \code{0.05}.
#' @param ncol Number of columns in the combined panel grid. Defaults to
#'   \code{length(images)} when multiple features and multiple images are
#'   both given (so each row is one feature across samples), otherwise
#'   \code{min(n_panels, 4)}.
#' @param combine Logical; if \code{TRUE} (default), returns a single
#'   combined \code{patchwork} object. Same-feature panels across images
#'   share one color scale/legend; different features get their own.
#'   If \code{FALSE}, returns a flat list of per-panel \code{ggplot} objects.
#' @param colors Color gradient to use. Defaults to a reversed
#'   \code{RColorBrewer} "Spectral" palette.
#' @param cells Optional character vector of cell/spot barcodes to restrict
#'   plotting to. Defaults to all cells present in each image.
#' @return A combined \code{patchwork} object (default), or a list of
#'   \code{ggplot} objects if \code{combine = FALSE}.
#' @examples
#' \dontrun{
#' SpatialFeaturePlotFixed(visium, features = "Alb")
#' SpatialFeaturePlotFixed(visium, features = c("Alb", "Cyp2e1"), pt.size = 2)
#' }
#' @importFrom ggplot2 aes geom_point labs scale_color_gradientn
#' @importFrom Seurat FetchData
#' @export
SpatialFeaturePlotFixed <- function(obj,
                                    features,
                                    images      = NULL,
                                    assay       = NULL,
                                    layer       = "data",
                                    pt.size     = 1.5,
                                    pt.alpha    = 1,
                                    pt.stroke   = 0,
                                    image.alpha = 1,
                                    crop        = TRUE,
                                    pad         = 0.05,
                                    ncol        = NULL,
                                    combine     = TRUE,
                                    colors      = NULL,
                                    cells       = NULL) {
  .assert_seurat(obj)
  if (length(obj@images) == 0) {
    stop("`obj` has no images in obj@images -- SpatialFeaturePlotFixed() ",
         "needs tissue images.")
  }
  if (is.null(images)) {
    images <- names(obj@images)
  } else {
    missing_images <- setdiff(images, names(obj@images))
    if (length(missing_images)) {
      stop("Image(s) not found in obj@images: ",
           paste(missing_images, collapse = ", "), ". Available: ",
           paste(names(obj@images), collapse = ", "))
    }
  }
  if (!is.character(features) || length(features) == 0) {
    stop("`features` must be a non-empty character vector.")
  }

  if (!is.null(assay)) Seurat::DefaultAssay(obj) <- assay

  fetch_cells <- if (is.null(cells)) colnames(obj) else intersect(colnames(obj), cells)
  vals <- Seurat::FetchData(obj, vars = features, layer = layer, cells = fetch_cells)

  if (is.null(colors)) {
    colors <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
  }

  x <- y <- value <- NULL  # NSE

  panels <- list()
  for (feat in features) {
    feat_vals <- vals[[feat]]
    names(feat_vals) <- rownames(vals)
    feat_range <- range(feat_vals, na.rm = TRUE)

    for (img in images) {
      img_cells <- SeuratObject::Cells(obj[[img]])
      img_cells <- intersect(img_cells, names(feat_vals))
      if (!is.null(cells)) img_cells <- intersect(img_cells, cells)
      base <- .spatial_panel_base(obj, img, crop = crop, image.alpha = image.alpha,
                                  pad = pad, cells = img_cells)
      pc <- base$coords
      pc$value <- feat_vals[pc$cell]

      title <- if (length(features) > 1 && length(images) > 1) {
        paste0(feat, " - ", img)
      } else if (length(features) > 1) {
        feat
      } else if (length(images) > 1) {
        img
      } else {
        NULL
      }

      p <- base$plot +
        ggplot2::geom_point(
          data    = pc,
          mapping = ggplot2::aes(x = x, y = y, color = value),
          size    = pt.size, alpha = pt.alpha, stroke = pt.stroke
        ) +
        ggplot2::scale_color_gradientn(colors = colors, limits = feat_range,
                                       name = feat) +
        ggplot2::labs(title = title)

      panels[[length(panels) + 1]] <- p
    }
  }

  if (is.null(ncol) && length(features) > 1 && length(images) > 1) {
    ncol <- length(images)
  }

  .combine_spatial_plots(panels, ncol = ncol, combine = combine)
}
