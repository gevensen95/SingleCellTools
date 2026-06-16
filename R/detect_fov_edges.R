#' Detect cells near the outer edge of a spatial FOV
#'
#' Two methods are supported, chosen by \code{method}:
#'
#' \describe{
#'   \item{\code{"bbox"} (default)}{Compute the axis-aligned bounding box
#'     around the cells in the FOV, then mark each cell whose minimum
#'     perpendicular distance to a box edge is below a threshold. The
#'     threshold scales with the median nearest-neighbor distance, so it
#'     adapts automatically to whatever coordinate units you're using
#'     (micron, pixel, etc.). This is the fast, deterministic option and
#'     is the right tool when the tissue (or imaged region) is roughly
#'     convex within each FOV — e.g. CosMx / Xenium tiles.}
#'   \item{\code{"angular"}}{For each cell, look at the angles to its
#'     \code{k} nearest neighbors and flag the cell if the largest gap
#'     exceeds \code{gap_threshold} AND its mean k-NN distance exceeds
#'     \code{density_factor * median(mean_nn_dist)}. Catches concave tissue
#'     edges, tear margins, and internal sparse regions that the bbox
#'     method misses. Slower and noisier.}
#' }
#'
#' Detection is run iteratively in both modes. Pass 1 identifies the outer
#' boundary; the next pass removes those cells, recomputes (the bounding
#' box or the angular statistics) on the survivors, and identifies the
#' new outer boundary, etc. The output metadata column records the
#' iteration in which each cell was caught:
#' \code{0} for interior, \code{1} for the outermost ring, \code{2} for
#' the next ring, ...
#'
#' @param obj A Seurat object with spatial coordinates.
#' @param fovs Character vector of FOV / image names. \code{NULL} (default)
#'   processes every image attached to the object.
#' @param method \code{"bbox"} (default) or \code{"angular"}. See Details.
#' @param bbox_factor (bbox method) Threshold multiplier on the median
#'   nearest-neighbor distance. A cell is flagged as an edge if its
#'   distance to the bounding box is less than
#'   \code{bbox_factor * median_nn_distance}. Default 2 — roughly "within
#'   two cell-widths of the box edge". Bump higher to catch a thicker
#'   ring per pass; drop lower to catch only the outermost ring.
#' @param k (angular method) Number of nearest neighbors per cell.
#'   Default 10.
#' @param gap_threshold (angular method) Minimum angular gap in radians.
#'   Default \code{2 * pi / 3} (120 degrees).
#' @param density_factor (angular method) Multiplier on the median mean-NN
#'   distance for the density check. Default 1. Set <= 0 to disable.
#' @param n_iterations Number of edge-peeling iterations. Default 2.
#' @param label_col Name of the metadata column to write. Default
#'   \code{"edge_layer"}. \code{0} = interior, \code{1} = outermost ring,
#'   \code{2} = next ring, ...
#' @return The Seurat object with a new integer metadata column.
#' @importFrom Seurat GetTissueCoordinates
#' @importFrom RANN nn2
#' @importFrom stats median
#' @export
detect_fov_edges <- function(obj,
                             fovs           = NULL,
                             method         = c("bbox", "angular"),
                             bbox_factor    = 1,
                             k              = 10,
                             gap_threshold  = 2 * pi / 3,
                             density_factor = 1,
                             n_iterations   = 2,
                             label_col      = "edge_layer") {

  method <- match.arg(method)
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (is.null(fovs)) fovs <- names(obj@images)
  if (!length(fovs)) stop("No FOVs in obj@images.")
  missing_fovs <- setdiff(fovs, names(obj@images))
  if (length(missing_fovs)) {
    stop("FOV(s) not in obj@images: ", paste(missing_fovs, collapse = ", "))
  }
  if (n_iterations < 1) stop("`n_iterations` must be >= 1.")

  all_layers <- setNames(integer(ncol(obj)), colnames(obj))

  for (fov in fovs) {
    coords <- Seurat::GetTissueCoordinates(obj[[fov]], which = "centroids")
    if ("cell" %in% colnames(coords)) {
      rownames(coords) <- coords$cell
      cells_in_fov <- intersect(coords$cell, colnames(obj))
      coords_mat   <- as.matrix(coords[cells_in_fov, c("x", "y")])
    } else {
      cells_in_fov <- intersect(rownames(coords), colnames(obj))
      coords_mat   <- as.matrix(coords[cells_in_fov, c("x", "y")])
    }
    if (length(cells_in_fov) < max(k + 1, 4)) {
      message(sprintf("  '%s': only %d cells, skipping",
                      fov, length(cells_in_fov)))
      next
    }

    # Median nearest-neighbor distance is the unit for both methods'
    # thresholds. Computed once across all cells in the FOV so the unit
    # is stable across iterations (the global density estimate, rather
    # than the active subset).
    nn_seed_k <- min(10, nrow(coords_mat) - 1)
    nn_seed   <- RANN::nn2(coords_mat, coords_mat,
                           k = nn_seed_k + 1)$nn.dists[, -1, drop = FALSE]
    median_nn <- stats::median(rowMeans(nn_seed))

    message(sprintf("--- '%s': %d cells, method=%s, %d iteration(s) ---",
                    fov, length(cells_in_fov), method, n_iterations))

    layer  <- integer(nrow(coords_mat))
    active <- seq_len(nrow(coords_mat))

    for (iter in seq_len(n_iterations)) {
      if (length(active) < max(k + 1, 4)) {
        message(sprintf("  Layer %d: too few cells remaining, stopping", iter))
        break
      }
      sub <- coords_mat[active, , drop = FALSE]

      if (method == "bbox") {
        # --- Bounding-box method ----------------------------------------
        xmin <- min(sub[, 1]); xmax <- max(sub[, 1])
        ymin <- min(sub[, 2]); ymax <- max(sub[, 2])
        dist_to_box <- pmin(
          sub[, 1] - xmin,
          xmax     - sub[, 1],
          sub[, 2] - ymin,
          ymax     - sub[, 2]
        )
        cutoff  <- bbox_factor * median_nn
        is_edge <- dist_to_box <= cutoff

      } else {
        # --- Angular-gap method (with density gate) ---------------------
        k_use   <- min(k, nrow(sub) - 1)
        nn_res  <- RANN::nn2(sub, sub, k = k_use + 1)
        nn      <- nn_res$nn.idx[, -1, drop = FALSE]
        nn_dist <- nn_res$nn.dists[, -1, drop = FALSE]
        mean_nn_dist <- rowMeans(nn_dist)
        density_cutoff <- if (density_factor > 0) {
          stats::median(mean_nn_dist) * density_factor
        } else -Inf
        is_edge <- vapply(seq_len(nrow(sub)), function(i) {
          neighbors <- sub[nn[i, ], , drop = FALSE]
          dx <- neighbors[, 1] - sub[i, 1]
          dy <- neighbors[, 2] - sub[i, 2]
          angles <- sort(atan2(dy, dx))
          gaps <- c(diff(angles),
                    (2 * pi + angles[1]) - angles[length(angles)])
          (max(gaps) > gap_threshold) &&
            (mean_nn_dist[i] > density_cutoff)
        }, logical(1))
      }

      layer[active[is_edge]] <- iter
      message(sprintf("  Layer %d: %d cells flagged", iter, sum(is_edge)))
      active <- active[!is_edge]
    }

    all_layers[cells_in_fov] <- layer
  }

  obj[[label_col]] <- all_layers[colnames(obj)]
  obj
}
