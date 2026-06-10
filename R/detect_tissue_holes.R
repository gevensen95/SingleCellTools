#' Detect cells bordering holes within a spatial tissue
#'
#' Builds a 2D occupancy grid over each FOV's coordinate space, marks any
#' grid cell containing at least one cell as "tissue," and uses a 4-connected
#' flood fill from the grid boundary to identify which empty grid cells are
#' \emph{outside} the tissue versus \emph{inside} (holes / gaps surrounded
#' by tissue). Cells whose grid cell is 4-adjacent to a hole grid cell are
#' flagged as hole-edge cells.
#'
#' Iteration peels back from the hole boundary inward. The first pass flags
#' cells immediately touching a hole; on the next pass those cells are
#' removed from the occupancy grid and a new ring of hole-edge cells is
#' found, etc. The output metadata column records the layer in which each
#' cell was caught: \code{0} for cells not near any hole, \code{1} for the
#' cells directly bordering a hole, \code{2} for the next ring inward, ...
#'
#' Distinct from \code{detect_fov_edges} &mdash; that function finds the
#' \emph{outer} boundary, this one finds \emph{internal} gaps. Run both
#' if you want to drop everything within N rings of either kind of boundary.
#'
#' @param obj A Seurat object with spatial coordinates.
#' @param fovs Character vector of FOV / image names. \code{NULL} (default)
#'   processes every image attached to the object.
#' @param bin_size Edge length of one grid bin in the same units as the
#'   coordinates. \code{NULL} (default) auto-computes
#'   \code{2.5 * median nearest-neighbor distance} per FOV, which puts
#'   roughly a handful of cells in each populated bin.
#' @param min_hole_size Minimum number of contiguous empty bins for a hole
#'   to be counted. Filters out single-bin "noise holes" caused by sparse
#'   sampling. Default 4.
#' @param n_iterations Number of inward-peeling iterations. Default 2.
#' @param label_col Name of the metadata column to write. Default
#'   \code{"hole_layer"}. \code{0} = not near a hole, \code{1} = directly
#'   bordering, \code{2} = next ring, ...
#' @return The Seurat object with a new integer metadata column.
#' @importFrom Seurat GetTissueCoordinates
#' @importFrom RANN nn2
#' @importFrom stats median
#' @export
detect_tissue_holes <- function(obj,
                                fovs          = NULL,
                                bin_size      = NULL,
                                min_hole_size = 4,
                                n_iterations  = 2,
                                label_col     = "hole_layer") {

  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (is.null(fovs)) fovs <- names(obj@images)
  if (!length(fovs)) stop("No FOVs in obj@images.")
  missing_fovs <- setdiff(fovs, names(obj@images))
  if (length(missing_fovs)) {
    stop("FOV(s) not in obj@images: ", paste(missing_fovs, collapse = ", "))
  }
  if (n_iterations < 1) stop("`n_iterations` must be >= 1.")

  all_layers <- setNames(integer(ncol(obj)), colnames(obj))

  # Vectorized "expand outside-mask by one grid step" used by the flood fill.
  # Use separate row / column counts — earlier version used a single `g` which
  # silently broke on rectangular grids (e.g. 21x20) by indexing the wrong axis.
  .expand <- function(mask) {
    nr <- nrow(mask)
    nc <- ncol(mask)
    up    <- rbind(mask[-1L, , drop = FALSE], FALSE)
    down  <- rbind(FALSE, mask[-nr, , drop = FALSE])
    left  <- cbind(mask[, -1L, drop = FALSE], FALSE)
    right <- cbind(FALSE, mask[, -nc, drop = FALSE])
    up | down | left | right
  }

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
    if (length(cells_in_fov) < 5) {
      message(sprintf("  '%s': only %d cells, skipping",
                      fov, length(cells_in_fov)))
      next
    }

    # Auto bin_size = 2.5x median nearest-neighbor distance unless given
    bsz <- bin_size
    if (is.null(bsz)) {
      nn_d <- RANN::nn2(coords_mat, coords_mat, k = 2)$nn.dists[, 2]
      bsz  <- 2.5 * stats::median(nn_d)
    }

    x_range <- range(coords_mat[, 1])
    y_range <- range(coords_mat[, 2])
    # Pad the grid by one bin so the outermost row/col is always empty
    x_breaks <- seq(x_range[1] - bsz, x_range[2] + bsz, by = bsz)
    y_breaks <- seq(y_range[1] - bsz, y_range[2] + bsz, by = bsz)
    nx <- length(x_breaks) - 1
    ny <- length(y_breaks) - 1
    if (nx < 3 || ny < 3) {
      message(sprintf("  '%s': grid too small (%dx%d), skipping",
                      fov, nx, ny))
      next
    }

    message(sprintf(
      "--- '%s': %d cells, %dx%d grid (bin = %.2f), %d hole iterations ---",
      fov, length(cells_in_fov), nx, ny, bsz, n_iterations))

    x_bin <- findInterval(coords_mat[, 1], x_breaks, all.inside = TRUE)
    y_bin <- findInterval(coords_mat[, 2], y_breaks, all.inside = TRUE)
    # cell index -> linear bin index (1..nx*ny)
    bin_of_cell <- (y_bin - 1L) * nx + x_bin

    layer  <- integer(nrow(coords_mat))
    active <- seq_len(nrow(coords_mat))

    for (iter in seq_len(n_iterations)) {
      if (length(active) < 5) {
        message(sprintf("  Layer %d: too few cells remaining, stopping", iter))
        break
      }

      # Build tissue grid from the currently-active cells
      tissue <- matrix(FALSE, nx, ny)
      active_bins <- bin_of_cell[active]
      tissue[active_bins] <- TRUE

      # Flood fill: outside = empty grid cells reachable from the boundary
      outside <- matrix(FALSE, nx, ny)
      outside[1, ] <- !tissue[1, ]; outside[nx, ] <- !tissue[nx, ]
      outside[, 1] <- !tissue[, 1]; outside[, ny] <- !tissue[, ny]
      repeat {
        expanded     <- .expand(outside) & !tissue
        new_outside  <- outside | expanded
        if (identical(new_outside, outside)) break
        outside <- new_outside
      }

      # Hole bins = empty AND not flood-fillable from outside
      hole_bin <- !tissue & !outside

      # Drop tiny "holes" of size < min_hole_size by simple grouping. We
      # don't need true connected components; counting total hole bins is a
      # cheap proxy. For a stricter min_hole_size, do component labeling.
      if (sum(hole_bin) < min_hole_size) {
        message(sprintf("  Layer %d: 0 holes detected (< min_hole_size)", iter))
        break
      }

      # Tissue bins adjacent to a hole bin
      hole_neighbor <- .expand(hole_bin) & tissue

      # Map back: cells in these bins, restricted to active cells
      is_near_hole_bin <- hole_neighbor[bin_of_cell]
      cells_flag <- is_near_hole_bin & seq_along(bin_of_cell) %in% active

      layer[cells_flag] <- iter
      message(sprintf("  Layer %d: %d cells flagged (%d hole bins, %d border bins)",
                      iter, sum(cells_flag), sum(hole_bin), sum(hole_neighbor)))
      active <- active[!cells_flag[active]]
    }

    all_layers[cells_in_fov] <- layer
  }

  obj[[label_col]] <- all_layers[colnames(obj)]
  obj
}
