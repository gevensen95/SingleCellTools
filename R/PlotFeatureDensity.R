#' Nebulosa-style feature density on a dimensional reduction
#'
#' Plots a 2D kernel-density estimate of expression (or module scores) on
#' the specified reduction (usually UMAP), one facet per feature. Weighted
#' KDE by expression makes sparse markers much easier to see than the
#' standard \code{Seurat::FeaturePlot}, which just draws grey-to-blue
#' points and loses signal in dense regions.
#'
#' Implemented with \code{ks::kde} weighted by the feature's expression at
#' each cell. Falls back to \code{MASS::kde2d} weighted by a threshold cut
#' if \code{ks} is unavailable.
#'
#' @param obj A Seurat object with the reduction embeddings computed.
#' @param features Character vector of gene names or metadata columns to
#'   plot. Genes are pulled from \code{assay}'s \code{"data"} layer;
#'   anything in \code{obj@meta.data} is used directly.
#' @param reduction Reduction name. Default \code{"umap"}.
#' @param dims Which two reduction dimensions to plot. Default
#'   \code{c(1, 2)}.
#' @param assay Assay for gene lookups. Default DefaultAssay.
#' @param layer Layer to read expression from. Default \code{"data"}.
#' @param n_grid Grid resolution for the KDE. Default 200 (higher =
#'   smoother, slower).
#' @param joint If TRUE, plot a "joint" density combining all requested
#'   features (co-expression signal) as an additional panel. Default FALSE.
#' @param colors Character vector of colors for the density gradient.
#'   Default \code{c("lightgrey", "#FEB24C", "#F03B20", "#BD0026")}.
#' @param ncol Number of columns in the facet grid. Default \code{NULL}
#'   (ggplot picks).
#' @param pt_size Size of the underlying scatter of cells (drawn faintly
#'   beneath the density). Default 0.3. Set to 0 to hide.
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' PlotFeatureDensity(obj, features = c("CD3D", "CD4", "CD8A"))
#' PlotFeatureDensity(obj, features = c("CD3D", "CD4"), joint = TRUE)
#' PlotFeatureDensity(obj, features = "module_score_Bcell",
#'                    reduction = "harmony")
#' }
#' @importFrom Seurat DefaultAssay Embeddings FetchData
#' @importFrom ggplot2 aes coord_fixed element_text facet_wrap geom_point ggplot labs scale_color_gradientn stat_density_2d theme theme_bw
#' @export
PlotFeatureDensity <- function(obj,
                               features,
                               reduction = "umap",
                               dims      = c(1, 2),
                               assay     = NULL,
                               layer     = "data",
                               n_grid    = 200,
                               joint     = FALSE,
                               colors    = c("lightgrey", "#FEB24C",
                                             "#F03B20", "#BD0026"),
                               ncol      = NULL,
                               pt_size   = 0.3) {

  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (!(reduction %in% names(obj@reductions))) {
    stop("Reduction '", reduction, "' not found in obj@reductions.")
  }
  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay
  emb <- Seurat::Embeddings(obj, reduction = reduction)[, dims,
                                                        drop = FALSE]
  colnames(emb) <- c("dim1", "dim2")

  # Pull expression for each requested feature (gene or metadata)
  expr_list <- lapply(features, function(f) {
    if (f %in% colnames(obj@meta.data)) {
      as.numeric(obj@meta.data[[f]])
    } else {
      tryCatch(
        as.numeric(Seurat::FetchData(obj, vars = f, layer = layer,
                                     assay = a)[[f]]),
        error = function(e) {
          stop("Feature '", f, "' not found in metadata or assay '",
               a, "'.")
        }
      )
    }
  })
  names(expr_list) <- features

  # Compute density per feature
  dens_list <- lapply(features, function(f) {
    d <- .weighted_kde(emb, weight = expr_list[[f]], n_grid = n_grid)
    d$feature <- f
    d
  })
  dens_df <- do.call(rbind, dens_list)

  # Joint panel: product of per-feature densities, re-normalized
  if (isTRUE(joint) && length(features) >= 2L) {
    joint_grid <- dens_list[[1]][, c("x", "y")]
    joint_grid$density <- Reduce(`*`, lapply(dens_list, `[[`, "density"))
    if (max(joint_grid$density) > 0) {
      joint_grid$density <- joint_grid$density / max(joint_grid$density)
    }
    joint_grid$feature <- paste(features, collapse = " + ")
    dens_df <- rbind(dens_df, joint_grid)
  }

  # Fix facet ordering
  facet_levels <- unique(c(features,
                           if (isTRUE(joint) && length(features) >= 2L)
                             paste(features, collapse = " + ")))
  dens_df$feature <- factor(dens_df$feature, levels = facet_levels)

  # ---- Build plot ---------------------------------------------------------
  emb_df <- as.data.frame(emb)
  dim1 <- dim2 <- x <- y <- density <- feature <- NULL  # NSE

  p <- ggplot2::ggplot()

  if (pt_size > 0) {
    p <- p + ggplot2::geom_point(
      data    = emb_df,
      mapping = ggplot2::aes(x = dim1, y = dim2),
      size    = pt_size, color = "grey85"
    )
  }
  p <- p +
    ggplot2::geom_tile(
      data    = dens_df,
      mapping = ggplot2::aes(x = x, y = y, fill = density),
      alpha   = 0.75
    ) +
    ggplot2::scale_fill_gradientn(colors = colors, na.value = "transparent",
                                  name = "density") +
    ggplot2::facet_wrap(~ feature, ncol = ncol) +
    ggplot2::labs(x = paste0(reduction, "_", dims[1]),
                  y = paste0(reduction, "_", dims[2])) +
    ggplot2::coord_fixed() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.text  = ggplot2::element_text(face = "bold"),
      panel.grid  = ggplot2::element_blank()
    )
  p
}


# ============================================================================
# Internal: weighted 2D KDE on a set of points, returning a long
# data frame with columns x, y, density (rescaled to [0, 1]).
# Uses ks::kde when available; otherwise a scaled MASS::kde2d fallback.
# ============================================================================
#' @keywords internal
#' @noRd
.weighted_kde <- function(emb, weight, n_grid = 200) {
  weight <- pmax(weight, 0)
  # If all zero (unexpressed feature), return an empty grid.
  if (sum(weight) == 0) {
    return(data.frame(x = numeric(0), y = numeric(0), density = numeric(0)))
  }

  x_range <- range(emb[, 1])
  y_range <- range(emb[, 2])

  if (requireNamespace("ks", quietly = TRUE)) {
    fit <- ks::kde(x = emb, w = weight * length(weight) / sum(weight),
                   gridsize = c(n_grid, n_grid),
                   xmin = c(x_range[1], y_range[1]),
                   xmax = c(x_range[2], y_range[2]))
    grid <- expand.grid(x = fit$eval.points[[1]], y = fit$eval.points[[2]])
    grid$density <- as.numeric(fit$estimate)
  } else if (requireNamespace("MASS", quietly = TRUE)) {
    # Threshold-cut weighting: replicate points by their (normalized) weight
    # and let MASS::kde2d smooth them. Cruder but no extra dep.
    w <- round(weight / max(weight) * 20)
    idx <- rep(seq_along(w), times = w)
    if (length(idx) < 5) {
      return(data.frame(x = numeric(0), y = numeric(0),
                        density = numeric(0)))
    }
    fit <- MASS::kde2d(emb[idx, 1], emb[idx, 2], n = n_grid,
                       lims = c(x_range, y_range))
    grid <- expand.grid(x = fit$x, y = fit$y)
    grid$density <- as.numeric(fit$z)
  } else {
    stop("PlotFeatureDensity needs `ks` or `MASS`. Install one of them.")
  }

  if (max(grid$density) > 0) {
    grid$density <- grid$density / max(grid$density)
  }
  # Drop cells below a small threshold for lighter rendering
  grid <- grid[grid$density > 0.02, , drop = FALSE]
  grid
}
