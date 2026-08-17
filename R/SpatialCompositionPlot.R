#' Multi-cell-type composition pie plot for spatial spots
#'
#' Draws a pie glyph at each spot's tissue coordinate showing its full
#' per-cell-type mixture (e.g. RCTD \code{"full"}-mode
#' \code{rctd_<celltype>} proportions from \code{\link{RunRCTD}}), instead
#' of collapsing each spot to a single dominant color the way
#' \code{Seurat::SpatialDimPlot()} would on \code{rctd_dominant}. Requires
#' the \code{scatterpie} package (\code{Suggests} only).
#'
#' @details
#' Visium objects can have thousands of spots; drawing one pie per spot gets
#' visually and computationally unusable well before that. \code{n_spots_max}
#' caps how many spots are plotted (a random subsample beyond that, with a
#' warning). \code{pie_scale} controls glyph size directly rather than
#' trying to infer a good default from spot spacing.
#'
#' @param obj A Seurat object.
#' @param weight_cols Character vector of numeric metadata columns holding
#'   each spot's composition. Defaults to every \code{rctd_<celltype>}
#'   column \code{\link{RunRCTD}} writes (excluding \code{rctd_dominant},
#'   which is a label, not a weight). \code{scatterpie} normalizes each
#'   spot's slice sizes to that spot's own row total, so columns don't need
#'   to already sum to exactly 1.
#' @param image Name of the image in \code{obj@images} to pull spot
#'   coordinates from. Defaults to the first image.
#' @param pie_scale Pie glyph size, passed to
#'   \code{scatterpie::geom_scatterpie(pie_scale = ...)}. Default \code{0.4}.
#' @param n_spots_max Maximum number of spots to plot; if more are present,
#'   a random subsample of this size is drawn and a warning is emitted.
#'   Default \code{2000}.
#' @param donut Logical; if \code{TRUE}, draws a donut instead of a solid
#'   pie (\code{scatterpie::geom_scatterpie(donut_radius = 0.5)}). Default
#'   \code{FALSE}.
#' @param colors Optional vector of fill colors, one per \code{weight_cols}
#'   entry (in order). Defaults to \code{ggplot2}'s standard discrete
#'   palette.
#' @param seed Random seed used for the \code{n_spots_max} subsample, for
#'   reproducibility. Default \code{1}.
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' visium <- RunRCTD(visium, reference = ref, celltype_col = "cell_type")
#' SpatialCompositionPlot(visium)
#' SpatialCompositionPlot(visium, donut = TRUE, pie_scale = 0.6)
#' }
#' @importFrom ggplot2 aes coord_fixed ggplot labs scale_fill_manual
#' @export
SpatialCompositionPlot <- function(obj,
                                   weight_cols = NULL,
                                   image       = NULL,
                                   pie_scale   = 0.4,
                                   n_spots_max = 2000,
                                   donut       = FALSE,
                                   colors      = NULL,
                                   seed        = 1) {
  if (!requireNamespace("scatterpie", quietly = TRUE)) {
    stop("Package 'scatterpie' is required for SpatialCompositionPlot(). ",
         "Install with: install.packages('scatterpie')")
  }
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")

  if (is.null(weight_cols)) {
    weight_cols <- grep("^rctd_", colnames(obj@meta.data), value = TRUE)
    weight_cols <- setdiff(weight_cols, "rctd_dominant")
    if (length(weight_cols) == 0) {
      stop("No `rctd_<celltype>` columns found in obj@meta.data and ",
           "`weight_cols` wasn't supplied. Run RunRCTD() first, or pass ",
           "`weight_cols` explicitly.")
    }
  }
  missing_cols <- setdiff(weight_cols, colnames(obj@meta.data))
  if (length(missing_cols)) {
    stop("Column(s) not found in obj@meta.data: ",
         paste(missing_cols, collapse = ", "))
  }
  not_numeric <- weight_cols[!vapply(obj@meta.data[weight_cols], is.numeric,
                                     logical(1))]
  if (length(not_numeric)) {
    stop("`weight_cols` must be numeric columns; not numeric: ",
         paste(not_numeric, collapse = ", "))
  }

  if (is.null(image)) {
    if (length(obj@images) == 0) {
      stop("`obj` has no images in obj@images -- SpatialCompositionPlot() ",
           "needs tissue coordinates.")
    }
    image <- names(obj@images)[1]
  } else if (!image %in% names(obj@images)) {
    stop("Image '", image, "' not found in obj@images. Available: ",
         paste(names(obj@images), collapse = ", "))
  }

  coords <- as.data.frame(.get_fov_coords(obj, image))
  coords$cell <- rownames(coords)

  md <- obj@meta.data[coords$cell, weight_cols, drop = FALSE]
  df <- cbind(coords, md)
  df <- df[stats::complete.cases(df), ]

  if (nrow(df) == 0) {
    stop("No spots with complete coordinates and weight_cols data to plot.")
  }
  if (nrow(df) > n_spots_max) {
    warning(sprintf(
      "%d spots exceeds n_spots_max = %d -- randomly subsampling to %d ",
      nrow(df), n_spots_max, n_spots_max),
      "for plotting speed/legibility. Set n_spots_max higher to plot more, ",
      "or lower to speed things up further.", call. = FALSE)
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
    df <- df[sample(nrow(df), n_spots_max), ]
  }

  x <- y <- cell <- NULL  # silence R CMD check NSE notes

  p <- ggplot2::ggplot() +
    scatterpie::geom_scatterpie(
      mapping      = ggplot2::aes(x = x, y = y, group = cell),
      data         = df,
      cols         = weight_cols,
      pie_scale    = pie_scale,
      donut_radius = if (isTRUE(donut)) 0.5 else 0,
      color        = NA
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(fill = "Cell type", x = NULL, y = NULL) +
    Ol_Reliable()

  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }
  p
}
