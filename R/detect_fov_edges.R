#' Detect cells near the outer edge of a spatial FOV
#'
#' For each cell, finds its \code{k} nearest neighbors, computes the
#' angular gap between consecutive neighbors (sorted by angle around the
#' cell), and marks the cell as an "edge" cell when \strong{both}:
#' \enumerate{
#'   \item the largest angular gap exceeds \code{gap_threshold} (the cell
#'     has a notable empty arc on one side), AND
#'   \item the cell's neighborhood is sparser than typical for the FOV
#'     (mean distance to its \code{k} nearest neighbors exceeds
#'     \code{density_factor} times the median across active cells).
#' }
#' The density requirement is what separates true outer-edge cells from
#' interior cells that happen to have one large angular gap because of
#' local density fluctuations.
#'
#' Detection is run iteratively. Pass 1 identifies the outer boundary;
#' on the next pass those cells are removed and the new outer boundary
#' is identified, etc. The output metadata column records the iteration
#' in which each cell was caught: \code{0} for interior, \code{1} for
#' the outermost ring, \code{2} for the next ring, ...
#'
#' Works for any spatial assay (Visium, Xenium, MERFISH, ...). For
#' Visium-specific edge detection that also considers tears and the
#' capture-area boundary, see \code{EdgeDetectionVisium}.
#'
#' @param obj A Seurat object with spatial coordinates.
#' @param fovs Character vector of FOV / image names. \code{NULL} (default)
#'   processes every image attached to the object.
#' @param k Number of nearest neighbors per cell. Default 10.
#' @param gap_threshold Minimum angular gap (in radians) for a cell to be
#'   considered an edge candidate. Default \code{2 * pi / 3} (120 degrees).
#' @param density_factor Multiplier on the median mean-neighbor distance.
#'   A cell only qualifies as an edge if its mean k-NN distance is greater
#'   than \code{density_factor * median(mean_nn_dist)}. Default 1 — set
#'   higher (e.g. 1.5, 2.0) to require sparser neighborhoods (fewer false
#'   positives, may miss some true edges); set to a value <= 0 to disable
#'   the density check entirely and use angular gap alone.
#' @param n_iterations Number of edge-peeling iterations. Default 2.
#' @param label_col Name of the metadata column to write. Default
#'   \code{"edge_layer"}. \code{0} = interior, \code{1} = outermost ring,
#'   \code{2} = next ring, ...
#' @return The Seurat object with a new integer metadata column
#'   \code{label_col}.
#' @importFrom Seurat GetTissueCoordinates
#' @importFrom RANN nn2
#' @importFrom stats median
#' @export
detect_fov_edges <- function(obj,
                             fovs           = NULL,
                             k              = 10,
                             gap_threshold  = 2 * pi / 3,
                             density_factor = 1,
                             n_iterations   = 2,
                             label_col      = "edge_layer") {

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
    if (length(cells_in_fov) < k + 1) {
      message(sprintf("  '%s': only %d cells, skipping",
                      fov, length(cells_in_fov)))
      next
    }

    message(sprintf("--- '%s': %d cells, %d edge iterations ---",
                    fov, length(cells_in_fov), n_iterations))

    layer  <- integer(nrow(coords_mat))
    active <- seq_len(nrow(coords_mat))

    for (iter in seq_len(n_iterations)) {
      if (length(active) < k + 1) {
        message(sprintf("  Layer %d: <%d cells remaining, stopping",
                        iter, k + 1))
        break
      }
      sub <- coords_mat[active, , drop = FALSE]
      k_use <- min(k, nrow(sub) - 1)

      nn_res  <- RANN::nn2(sub, sub, k = k_use + 1)
      nn      <- nn_res$nn.idx[, -1, drop = FALSE]
      nn_dist <- nn_res$nn.dists[, -1, drop = FALSE]

      # Per-cell mean distance to k neighbors — proxy for local density.
      # Edge cells should have larger mean NN distance than interior cells
      # because some of their "neighbors" are across an empty arc.
      mean_nn_dist <- rowMeans(nn_dist)
      density_cutoff <- if (density_factor > 0) {
        stats::median(mean_nn_dist) * density_factor
      } else {
        -Inf  # disable density gate; angular gap alone decides
      }

      is_edge <- vapply(seq_len(nrow(sub)), function(i) {
        neighbors <- sub[nn[i, ], , drop = FALSE]
        dx <- neighbors[, 1] - sub[i, 1]
        dy <- neighbors[, 2] - sub[i, 2]
        angles <- sort(atan2(dy, dx))
        gaps <- c(diff(angles),
                  (2 * pi + angles[1]) - angles[length(angles)])
        has_gap         <- max(gaps) > gap_threshold
        is_low_density  <- mean_nn_dist[i] > density_cutoff
        has_gap && is_low_density
      }, logical(1))

      layer[active[is_edge]] <- iter
      message(sprintf("  Layer %d: %d cells flagged", iter, sum(is_edge)))
      active <- active[!is_edge]
    }

    all_layers[cells_in_fov] <- layer
  }

  obj[[label_col]] <- all_layers[colnames(obj)]
  obj
}
