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
#'   When combined with \code{split.by}, the joint panel is computed
#'   separately within each split group.
#' @param split.by Optional metadata column name to split cells by (same
#'   idea as \code{Seurat::FeaturePlot(split.by = )}). Each level becomes
#'   its own column of panels, with \code{features} (and the \code{joint}
#'   panel, if requested) as rows -- a \code{facet_grid(feature ~
#'   split.by)} layout instead of the default \code{facet_wrap(~
#'   feature)}. The density grid for every panel is still evaluated over
#'   the FULL embedding's x/y range (not just that group's cells), so
#'   panels stay spatially aligned and comparable side by side; only the
#'   weighting (which cells contribute density) and the faint background
#'   scatter are restricted per group. Cells with \code{NA} in this column
#'   are dropped with a warning. Default \code{NULL} (no split).
#' @param colors Character vector of colors for the density gradient.
#'   Default \code{c("lightgrey", "#FEB24C", "#F03B20", "#BD0026")}.
#' @param ncol Number of columns in the facet grid. Ignored when
#'   \code{split.by} is set (that layout is a 2D grid driven by
#'   \code{feature} and \code{split.by}, not a wrapped single row).
#'   Default \code{NULL} (ggplot picks).
#' @param pt_size Size of the underlying scatter of cells (drawn faintly
#'   beneath the density). With \code{split.by}, each panel's scatter is
#'   restricted to that group's own cells. Default 0.3. Set to 0 to hide.
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' PlotFeatureDensity(obj, features = c("CD3D", "CD4", "CD8A"))
#' PlotFeatureDensity(obj, features = c("CD3D", "CD4"), joint = TRUE)
#' PlotFeatureDensity(obj, features = "module_score_Bcell",
#'                    reduction = "harmony")
#' PlotFeatureDensity(obj, features = c("CD3D", "CD4"), split.by = "genotype")
#' }
#' @importFrom Seurat DefaultAssay Embeddings FetchData
#' @importFrom ggplot2 aes coord_fixed element_text facet_grid facet_wrap geom_point ggplot labs scale_color_gradientn stat_density_2d theme theme_bw
#' @export
PlotFeatureDensity <- function(obj,
                               features,
                               reduction = "umap",
                               dims      = c(1, 2),
                               assay     = NULL,
                               layer     = "data",
                               n_grid    = 200,
                               joint     = FALSE,
                               split.by  = NULL,
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

  # ---- split.by: resolve groups, drop NA cells, reduce to one code path --
  # No split.by is treated as a single implicit group containing every
  # cell -- the density/joint logic below runs exactly once either way,
  # so the split and non-split cases can't silently drift apart.
  if (is.null(split.by)) {
    keep <- rep(TRUE, nrow(emb))
    group_f <- factor(rep("__all__", nrow(emb)))
  } else {
    if (!split.by %in% colnames(obj@meta.data)) {
      stop("`split.by` column '", split.by, "' not found in obj@meta.data.")
    }
    split_vals <- obj@meta.data[[split.by]]
    keep <- !is.na(split_vals)
    if (any(!keep)) {
      warning(sum(!keep), " cell(s) have NA in `split.by` ('", split.by,
              "') and were dropped.", call. = FALSE)
    }
    group_f <- if (is.factor(split_vals)) {
      droplevels(split_vals[keep])
    } else {
      factor(split_vals[keep])
    }
    if (!is.null(ncol)) {
      message("`ncol` is ignored when `split.by` is set ",
              "(layout uses facet_grid(feature ~ split.by)).")
    }
  }
  if (any(!keep)) {
    emb <- emb[keep, , drop = FALSE]
    expr_list <- lapply(expr_list, function(x) x[keep])
  }
  group_levels <- levels(group_f)

  # Grid bounds come from the FULL (post NA-drop) embedding, computed once,
  # so every panel -- whether split or not -- shares the same spatial
  # extent and is directly comparable. Each group is only ever restricted
  # in which cells WEIGHT the density, never in the grid it's evaluated on.
  x_range <- range(emb[, 1])
  y_range <- range(emb[, 2])

  .trim <- function(d) d[d$density > 0.02, , drop = FALSE]

  # Compute density (+ optional joint panel) once per group. With no
  # split.by this runs exactly once, over every cell -- identical to the
  # previous (pre-split.by) behavior.
  panel_list <- lapply(group_levels, function(g) {
    idx_g <- which(as.character(group_f) == g)
    emb_g <- emb[idx_g, , drop = FALSE]

    dens_full <- lapply(features, function(f) {
      .weighted_kde(emb_g, weight = expr_list[[f]][idx_g], n_grid = n_grid,
                    min_density = 0, x_range = x_range, y_range = y_range)
    })
    names(dens_full) <- features

    # The low-density rows are dropped per panel below, and separately
    # after the joint product -- never before it. Thresholding first is
    # what broke the joint panel: each feature keeps a different subset of
    # grid points, so the density vectors have different lengths, and
    # multiplying them made R recycle -- pairing feature A's density at
    # one grid point with feature B's at an unrelated one, and emitting
    # "longer object length is not a multiple of shorter object length".
    # The panel rendered, it was just wrong.
    dens_list <- lapply(features, function(f) {
      d <- .trim(dens_full[[f]])
      d$feature <- f
      d
    })
    dens_df_g <- do.call(rbind, dens_list)

    # Joint panel: product of per-feature densities, re-normalized
    if (isTRUE(joint) && length(features) >= 2L) {
      # Every feature is evaluated on the same lattice (same embedding,
      # same x/y range, same gridsize), so the untrimmed grids line up row
      # for row. Verify rather than assume: .weighted_kde() short-circuits
      # to a 0-row frame for a feature that is zero everywhere (within
      # this group), and multiplying through that would silently drop the
      # joint panel's meaning.
      ref <- dens_full[[1]]
      aligned <- nrow(ref) > 0 && all(vapply(dens_full, function(d) {
        nrow(d) == nrow(ref) &&
          isTRUE(all.equal(d$x, ref$x)) && isTRUE(all.equal(d$y, ref$y))
      }, logical(1)))

      if (!aligned) {
        warning("Skipping the joint density panel",
                if (!is.null(split.by)) paste0(" for ", split.by, " = '", g, "'") else "",
                ": the per-feature density grids don't align (usually ",
                "because one requested feature is zero in every cell",
                if (!is.null(split.by)) " of this group" else "",
                "). The individual panels are unaffected.", call. = FALSE)
      } else {
        joint_grid <- ref[, c("x", "y")]
        joint_grid$density <- Reduce(`*`, lapply(dens_full, `[[`, "density"))
        if (max(joint_grid$density) > 0) {
          joint_grid$density <- joint_grid$density / max(joint_grid$density)
        }
        joint_grid <- .trim(joint_grid)
        joint_grid$feature <- paste(features, collapse = " + ")
        dens_df_g <- rbind(dens_df_g, joint_grid)
      }
    }

    dens_df_g$split_group <- g
    dens_df_g
  })
  dens_df <- do.call(rbind, panel_list)

  # Fix facet ordering
  facet_levels <- unique(c(features,
                           if (isTRUE(joint) && length(features) >= 2L)
                             paste(features, collapse = " + ")))
  dens_df$feature <- factor(dens_df$feature, levels = facet_levels)

  # ---- Build plot ---------------------------------------------------------
  emb_df <- as.data.frame(emb)
  emb_df$split_group <- as.character(group_f)
  dim1 <- dim2 <- x <- y <- density <- feature <- split_group <- NULL  # NSE

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
    (if (is.null(split.by)) {
      ggplot2::facet_wrap(~ feature, ncol = ncol)
    } else {
      ggplot2::facet_grid(feature ~ split_group)
    }) +
    ggplot2::labs(x = paste0(reduction, "_", dims[1]),
                  y = paste0(reduction, "_", dims[2])) +
    ggplot2::coord_fixed() +
    Ol_Reliable() +
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
# x_range/y_range let a caller pin the grid to a fixed extent (e.g. the full
# embedding's range) rather than auto-deriving it from `emb` -- needed when
# `emb` is only a subset of cells (one split.by group) but the grid must
# still span the same area as every other group's panel for the panels to
# be visually comparable. NULL (the default) reproduces the old
# behavior: derive the range from whatever `emb` is passed.
# ============================================================================
#' @keywords internal
#' @noRd
.weighted_kde <- function(emb, weight, n_grid = 200, min_density = 0.02,
                          x_range = NULL, y_range = NULL) {
  weight <- pmax(weight, 0)
  # If all zero (unexpressed feature), return an empty grid.
  if (sum(weight) == 0) {
    return(data.frame(x = numeric(0), y = numeric(0), density = numeric(0)))
  }

  if (is.null(x_range)) x_range <- range(emb[, 1])
  if (is.null(y_range)) y_range <- range(emb[, 2])

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
  # Drop cells below a small threshold for lighter rendering. Callers that
  # need to combine grids across features MUST pass min_density = 0 and
  # trim afterwards -- trimming here first leaves each feature with a
  # different subset of grid points, so the grids no longer align.
  if (min_density > 0) {
    grid <- grid[grid$density > min_density, , drop = FALSE]
  }
  grid
}
