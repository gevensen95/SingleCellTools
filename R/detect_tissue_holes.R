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
#' cells directly bordering a hole, \code{2} for the next ring, ...
#'
#' \strong{Marker-gene exclusion.} Set \code{exclude_gene} to ignore holes
#' that are biologically meaningful gaps (a vessel, a sinusoid) rather
#' than tears. The function computes the mean expression of
#' \code{exclude_gene} per tissue bin (a "bin expression vector"), then
#' decides which holes look like real anatomy:
#' \itemize{
#'   \item \strong{Data-adaptive thresholding (default).} The bin
#'     expression vector's quantile at \code{sensitivity} is used as the
#'     cutoff. Bins above that cutoff count as "high expression." A hole
#'     is excluded if the mean of its border-bin expression exceeds the
#'     cutoff. This adapts automatically to each FOV's overall
#'     expression level and is the recommended path.
#'   \item \strong{Absolute thresholding.} If
#'     \code{exclude_gene_threshold} is set to a number, that value is
#'     used as the cutoff directly (regardless of the distribution). Use
#'     this when you want a fixed biological cutoff (e.g. "any hole
#'     bordered by cells expressing Glul above log-normalized 1.0 is a
#'     central vein").
#' }
#'
#' Distinct from \code{detect_fov_edges} — that function finds the
#' \emph{outer} boundary, this one finds \emph{internal} gaps. Run both
#' if you want to drop everything within N rings of either kind of
#' boundary.
#'
#' @param obj A Seurat object with spatial coordinates.
#' @param fovs Character vector of FOV / image names. \code{NULL} (default)
#'   processes every image attached to the object.
#' @param bin_size Edge length of one grid bin in the same units as the
#'   coordinates. \code{NULL} (default) auto-computes
#'   \code{2.5 * median nearest-neighbor distance} per FOV.
#' @param min_hole_size Minimum number of contiguous empty bins (4-connected)
#'   for a hole to count. Default 4.
#' @param exclude_gene Optional gene symbol. When supplied, holes whose
#'   bordering tissue bins have a mean expression of this gene above the
#'   chosen cutoff are skipped (their border cells are not flagged).
#'   \code{NULL} (default) disables the filter. Examples: \code{"Glul"}
#'   for liver central veins, \code{"Pecam1"} / \code{"Cdh5"} for
#'   endothelial-lined vessels.
#' @param sensitivity Quantile of the per-bin mean-expression
#'   distribution to use as the cutoff (only when
#'   \code{exclude_gene_threshold = NULL}). Default 0.75 — bins above
#'   the 75th percentile are treated as high-expression. Higher values
#'   (e.g. 0.9) are stricter (fewer holes excluded); lower values
#'   (e.g. 0.5) are more permissive. Range 0--1.
#' @param exclude_gene_threshold Absolute expression cutoff. When
#'   \code{NULL} (default), the cutoff is computed from
#'   \code{sensitivity}; set to a numeric value to override the
#'   data-adaptive path and use a fixed threshold instead.
#' @param exclude_gene_assay Assay to read expression from. \code{NULL}
#'   (default) uses \code{DefaultAssay(obj)}.
#' @param exclude_gene_layer Layer to read from. Default \code{"data"}
#'   (log-normalized).
#' @param n_iterations Number of inward-peeling iterations. Default 2.
#' @param label_col Name of the metadata column to write. Default
#'   \code{"hole_layer"}.
#' @return The Seurat object with a new integer metadata column.
#' @importFrom Seurat GetTissueCoordinates FetchData DefaultAssay
#' @importFrom RANN nn2
#' @importFrom stats median quantile
#' @export
detect_tissue_holes <- function(obj,
                                fovs                    = NULL,
                                bin_size                = NULL,
                                min_hole_size           = 4,
                                exclude_gene            = NULL,
                                sensitivity             = 0.75,
                                exclude_gene_threshold  = NULL,
                                exclude_gene_assay      = NULL,
                                exclude_gene_layer      = "data",
                                n_iterations            = 2,
                                label_col               = "hole_layer") {

  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (is.null(fovs)) fovs <- names(obj@images)
  if (!length(fovs)) stop("No FOVs in obj@images.")
  missing_fovs <- setdiff(fovs, names(obj@images))
  if (length(missing_fovs)) {
    stop("FOV(s) not in obj@images: ", paste(missing_fovs, collapse = ", "))
  }
  if (n_iterations < 1) stop("`n_iterations` must be >= 1.")
  if (sensitivity < 0 || sensitivity > 1) {
    stop("`sensitivity` must be between 0 and 1.")
  }

  use_gene_filter <- !is.null(exclude_gene)
  use_sensitivity <- use_gene_filter && is.null(exclude_gene_threshold)

  if (use_gene_filter) {
    a <- if (is.null(exclude_gene_assay)) Seurat::DefaultAssay(obj)
                                          else exclude_gene_assay
    if (!exclude_gene %in% rownames(obj[[a]])) {
      stop("`exclude_gene` ('", exclude_gene, "') not found in assay '",
           a, "'.")
    }
    mode_msg <- if (use_sensitivity) {
      sprintf("data-adaptive cutoff (sensitivity=%.2f)", sensitivity)
    } else {
      sprintf("absolute cutoff (threshold=%.3f)", exclude_gene_threshold)
    }
    message(sprintf("  exclude_gene = '%s' (assay=%s, layer=%s, %s)",
                    exclude_gene, a, exclude_gene_layer, mode_msg))
  }

  all_layers <- setNames(integer(ncol(obj)), colnames(obj))

  # Vectorized "expand mask by one grid step" (4-connected). Separate
  # nrow/ncol so it works on rectangular grids.
  .expand <- function(mask) {
    nr <- nrow(mask); nc <- ncol(mask)
    up    <- rbind(mask[-1L, , drop = FALSE], FALSE)
    down  <- rbind(FALSE, mask[-nr, , drop = FALSE])
    left  <- cbind(mask[, -1L, drop = FALSE], FALSE)
    right <- cbind(FALSE, mask[, -nc, drop = FALSE])
    up | down | left | right
  }

  # 4-connected component labeling.
  .label_components <- function(mask) {
    nr <- nrow(mask); nc <- ncol(mask)
    labels <- matrix(0L, nr, nc)
    current <- 0L
    for (start in which(mask)) {
      if (labels[start] != 0L) next
      current <- current + 1L
      queue <- start
      while (length(queue) > 0) {
        idx <- queue[1]; queue <- queue[-1]
        if (labels[idx] != 0L || !mask[idx]) next
        labels[idx] <- current
        i <- ((idx - 1L) %% nr) + 1L
        j <- ((idx - 1L) %/% nr) + 1L
        if (i > 1L && mask[i - 1L, j] && labels[i - 1L, j] == 0L)
          queue <- c(queue, (j - 1L) * nr + (i - 1L))
        if (i < nr && mask[i + 1L, j] && labels[i + 1L, j] == 0L)
          queue <- c(queue, (j - 1L) * nr + (i + 1L))
        if (j > 1L && mask[i, j - 1L] && labels[i, j - 1L] == 0L)
          queue <- c(queue, (j - 2L) * nr + i)
        if (j < nc && mask[i, j + 1L] && labels[i, j + 1L] == 0L)
          queue <- c(queue, j * nr + i)
      }
    }
    list(labels = labels, n = current)
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

    # Fetch the exclude_gene expression vector once per FOV. We'll
    # re-aggregate it to per-bin means inside the iteration loop since
    # the set of active cells changes.
    fov_expr <- NULL
    if (use_gene_filter) {
      fov_expr <- tryCatch({
        Seurat::FetchData(
          obj, vars = exclude_gene, cells = cells_in_fov,
          layer = exclude_gene_layer)[[exclude_gene]]
      }, error = function(e) {
        warning("  FetchData failed for '", exclude_gene,
                "' in FOV '", fov, "': ", conditionMessage(e),
                "; skipping gene filter for this FOV.")
        NULL
      })
      if (!is.null(fov_expr)) names(fov_expr) <- cells_in_fov
    }

    bsz <- bin_size
    if (is.null(bsz)) {
      nn_d <- RANN::nn2(coords_mat, coords_mat, k = 2)$nn.dists[, 2]
      bsz  <- 2.5 * stats::median(nn_d)
    }

    x_range <- range(coords_mat[, 1])
    y_range <- range(coords_mat[, 2])
    x_breaks <- seq(x_range[1] - bsz, x_range[2] + bsz, by = bsz)
    y_breaks <- seq(y_range[1] - bsz, y_range[2] + bsz, by = bsz)
    nx <- length(x_breaks) - 1L
    ny <- length(y_breaks) - 1L
    if (nx < 3 || ny < 3) {
      message(sprintf("  '%s': grid too small (%dx%d), skipping",
                      fov, nx, ny))
      next
    }

    message(sprintf(
      "--- '%s': %d cells, %dx%d grid (bin=%.2f), %d iter ---",
      fov, length(cells_in_fov), nx, ny, bsz, n_iterations))

    x_bin <- findInterval(coords_mat[, 1], x_breaks, all.inside = TRUE)
    y_bin <- findInterval(coords_mat[, 2], y_breaks, all.inside = TRUE)
    bin_of_cell <- (y_bin - 1L) * nx + x_bin

    layer  <- integer(nrow(coords_mat))
    active <- seq_len(nrow(coords_mat))

    for (iter in seq_len(n_iterations)) {
      if (length(active) < 5) {
        message(sprintf("  Layer %d: too few cells remaining, stopping", iter))
        break
      }

      tissue <- matrix(FALSE, nx, ny)
      tissue[bin_of_cell[active]] <- TRUE

      # Flood fill outside-mask from the grid boundary
      outside <- matrix(FALSE, nx, ny)
      outside[1, ]  <- !tissue[1, ];  outside[nx, ] <- !tissue[nx, ]
      outside[, 1]  <- !tissue[, 1];  outside[, ny] <- !tissue[, ny]
      repeat {
        nxt <- outside | (.expand(outside) & !tissue)
        if (identical(nxt, outside)) break
        outside <- nxt
      }
      hole_bin <- !tissue & !outside

      if (!any(hole_bin)) {
        message(sprintf("  Layer %d: 0 hole bins detected", iter))
        break
      }

      # Per-bin mean expression of the marker gene, computed over the
      # currently-active cells. Used both to derive the data-adaptive
      # cutoff and to score each hole's border bins.
      bin_mean_expr <- NULL
      bin_threshold <- NULL
      if (use_gene_filter && !is.null(fov_expr)) {
        active_cells <- cells_in_fov[active]
        active_expr  <- fov_expr[active_cells]
        bin_mean_expr <- tapply(active_expr,
                                bin_of_cell[active],
                                mean, na.rm = TRUE)
        # tapply returns a named vector keyed on bin index (as string)
        if (use_sensitivity) {
          bin_threshold <- stats::quantile(
            bin_mean_expr, probs = sensitivity, na.rm = TRUE)
        } else {
          bin_threshold <- exclude_gene_threshold
        }
        message(sprintf(
          "  Layer %d: gene '%s' bin-mean range %.3f to %.3f; cutoff %.3f",
          iter, exclude_gene,
          min(bin_mean_expr, na.rm = TRUE),
          max(bin_mean_expr, na.rm = TRUE),
          bin_threshold))
      }

      # Label connected components and decide per-hole.
      comp <- .label_components(hole_bin)
      keep_hole_mask <- matrix(FALSE, nx, ny)
      n_total <- comp$n; n_small <- 0L; n_excl <- 0L; n_kept <- 0L

      for (cid in seq_len(comp$n)) {
        comp_mask <- comp$labels == cid
        if (sum(comp_mask) < min_hole_size) {
          n_small <- n_small + 1L
          next
        }

        # Marker-gene exclusion: is the mean of border-bin expression
        # above the cutoff?
        if (use_gene_filter && !is.null(bin_mean_expr)) {
          border_bins <- .expand(comp_mask) & tissue
          if (any(border_bins)) {
            border_bin_idx <- which(border_bins)
            be <- bin_mean_expr[as.character(border_bin_idx)]
            be <- be[!is.na(be)]
            if (length(be) > 0 && mean(be) > bin_threshold) {
              n_excl <- n_excl + 1L
              next
            }
          }
        }

        keep_hole_mask[comp_mask] <- TRUE
        n_kept <- n_kept + 1L
      }

      if (n_kept == 0L) {
        message(sprintf(
          "  Layer %d: %d hole(s); %d below min_hole_size, %d gene-excluded; nothing flagged",
          iter, n_total, n_small, n_excl))
        break
      }

      # Tissue bins adjacent to a kept hole
      hole_neighbor <- .expand(keep_hole_mask) & tissue
      is_near <- hole_neighbor[bin_of_cell]
      cells_flag <- is_near & seq_along(bin_of_cell) %in% active
      layer[cells_flag] <- iter

      message(sprintf(
        "  Layer %d: %d hole(s) (%d kept; %d small, %d gene-excluded); %d cells flagged",
        iter, n_total, n_kept, n_small, n_excl, sum(cells_flag)))
      active <- active[!cells_flag[active]]
    }

    all_layers[cells_in_fov] <- layer
  }

  obj[[label_col]] <- all_layers[colnames(obj)]
  obj
}
