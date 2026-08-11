#' Spatial dimension plot, without Seurat's GeomSpatial point-size bug
#'
#' A from-scratch \code{ggplot2} replacement for
#' \code{Seurat::SpatialDimPlot()}: plots the tissue image with spots
#' colored by a discrete grouping variable on top. Built with plain
#' \code{geom_point()} rather than Seurat's custom \code{GeomSpatial} layer,
#' because \code{pt.size.factor} has real, documented regressions in recent
#' Seurat releases where it silently has no effect on rendered spot size
#' (\url{https://github.com/satijalab/seurat/issues/9491},
#' \url{https://github.com/satijalab/seurat/issues/6179}). \code{pt.size}
#' here is an ordinary \code{geom_point(size = ...)} value and is guaranteed
#' to change spot size.
#'
#' @param obj A Seurat object with at least one image in \code{obj@images}.
#' @param group.by Metadata column to color spots by. Defaults to
#'   \code{Seurat::Idents(obj)}.
#' @param images Character vector of image names in \code{obj@images} to
#'   plot. Defaults to all of them (one panel each, like
#'   \code{Seurat::SpatialDimPlot()}'s default multi-sample behavior).
#' @param pt.size,pt.alpha,pt.stroke Passed to \code{geom_point()} as
#'   \code{size}, \code{alpha}, \code{stroke}.
#' @param image.alpha Dims the tissue image toward white; \code{1} = full
#'   opacity, \code{0} = fully white. Default \code{1}.
#' @param crop Logical; if \code{TRUE} (default), zooms to the bounding box
#'   of the plotted spots (plus \code{pad}). If \code{FALSE}, shows the full
#'   image canvas, matching \code{Seurat::SpatialDimPlot(crop = FALSE)}.
#' @param pad Fractional padding added around the spot bounding box when
#'   \code{crop = TRUE}. Default \code{0.05}.
#' @param ncol Number of columns when multiple images are plotted. Defaults
#'   to \code{min(length(images), 4)}.
#' @param combine Logical; if \code{TRUE} (default), returns a single
#'   combined \code{patchwork} object with a shared/collected legend across
#'   panels. If \code{FALSE}, returns a list of per-image \code{ggplot}
#'   objects instead.
#' @param colors Optional vector of colors, one per level of \code{group.by}
#'   (in level order). Defaults to a \code{grDevices::hcl()} rainbow, same
#'   as \code{ggplot2}'s default discrete palette approach.
#' @param cells Optional character vector of cell/spot barcodes to restrict
#'   plotting to. Defaults to all cells present in each image.
#' @return A combined \code{patchwork} object (default), or a list of
#'   \code{ggplot} objects if \code{combine = FALSE}.
#' @examples
#' \dontrun{
#' SpatialDimPlotFixed(visium, group.by = "seurat_clusters")
#' SpatialDimPlotFixed(visium, group.by = "rctd_dominant", pt.size = 2.5)
#' }
#' @importFrom ggplot2 aes geom_point labs scale_color_manual
#' @export
SpatialDimPlotFixed <- function(obj,
                                group.by    = NULL,
                                images      = NULL,
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
    stop("`obj` has no images in obj@images -- SpatialDimPlotFixed() needs ",
         "tissue images.")
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

  if (is.null(group.by)) {
    group_raw <- Seurat::Idents(obj)
    legend_title <- "Identity"
  } else {
    if (!group.by %in% colnames(obj@meta.data)) {
      stop("`group.by` column '", group.by, "' not found in obj@meta.data.")
    }
    group_raw <- obj@meta.data[[group.by]]
    legend_title <- group.by
  }
  if (!is.factor(group_raw)) group_raw <- factor(group_raw)
  all_levels <- levels(group_raw)
  group_vals <- as.character(group_raw)
  names(group_vals) <- rownames(obj@meta.data)

  if (is.null(colors)) {
    n <- length(all_levels)
    colors <- grDevices::hcl(h = seq(15, 375, length.out = n + 1),
                             c = 100, l = 65)[seq_len(n)]
  } else if (length(colors) < length(all_levels)) {
    stop("`colors` has ", length(colors), " value(s) but `group.by` has ",
         length(all_levels), " level(s).")
  }

  x <- y <- group <- NULL  # NSE

  plots <- lapply(images, function(img) {
    img_cells <- SeuratObject::Cells(obj[[img]])
    if (!is.null(cells)) img_cells <- intersect(img_cells, cells)
    base <- .spatial_panel_base(obj, img, crop = crop, image.alpha = image.alpha,
                                pad = pad, cells = img_cells)
    pc <- base$coords
    pc$group <- factor(group_vals[pc$cell], levels = all_levels)

    base$plot +
      ggplot2::geom_point(
        data    = pc,
        mapping = ggplot2::aes(x = x, y = y, color = group),
        size    = pt.size, alpha = pt.alpha, stroke = pt.stroke
      ) +
      ggplot2::scale_color_manual(values = colors, limits = all_levels,
                                  drop = FALSE, name = legend_title,
                                  na.value = "grey80") +
      ggplot2::labs(title = if (length(images) > 1) img else NULL)
  })

  .combine_spatial_plots(plots, ncol = ncol, combine = combine)
}
