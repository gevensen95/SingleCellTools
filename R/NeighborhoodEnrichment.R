#' Neighborhood enrichment between cell types across one or many FOVs
#'
#' For each pair of \code{group.by} labels, counts how often they appear as
#' k-nearest neighbors of each other and compares to a random permutation
#' null. Returns z-scores, empirical p-values, BH-adjusted q-values, and
#' the underlying observed / expected count matrices.
#'
#' Neighbors are computed \strong{within each FOV} (cells in different FOVs
#' are never considered neighbors), but observed and null counts are pooled
#' across FOVs before computing the final statistics. This gives a single
#' tissue-level enrichment estimate while respecting FOV boundaries.
#'
#' \strong{Niche designation (optional).} When \code{assign_niches = TRUE},
#' the function also computes a per-cell neighborhood composition vector
#' (the cell-type fractions among each cell's k nearest neighbors),
#' clusters cells across all FOVs into \code{n_niches} groups using the
#' chosen method, and returns the per-cell niche labels in \code{$niche}.
#' If \code{add_to_meta = TRUE} (the default), these labels are also written
#' to \code{obj@meta.data[[niche_col]]} on a copy of \code{obj}, which is
#' returned as \code{$obj} in the result list alongside the stats -- e.g.
#' \code{res <- NeighborhoodEnrichment(combined, ...); combined <- res$obj}.
#' Cells not covered by any requested FOV (or with an \code{NA}
#' \code{group.by} label) get \code{NA} in this column.
#'
#' \strong{P-value calculation.} Two-sided empirical p-value: for each pair
#' (i, j), let \eqn{p_{up} = (1 + \#\{perm_{ij} \ge obs_{ij}\}) /
#' (n\_perm + 1)} and \eqn{p_{down}} symmetrically. The reported p-value is
#' \eqn{2 \cdot \min(p_{up}, p_{down})}, capped at 1.
#'
#' @param obj A Seurat object with spatial coordinates.
#' @param group.by Metadata column with cell-type / cluster labels.
#' @param fovs Character vector of FOV / image names. \code{NULL} (default)
#'   uses every image attached to the object.
#' @param k Number of nearest neighbors per cell. Default 10.
#' @param n_perm Number of permutations for the null. Default 100.
#' @param seed Random seed. Also reset immediately before niche clustering
#'   (when \code{assign_niches = TRUE}), so the niche assignment is
#'   reproducible on its own and doesn't shift when \code{n_perm} changes.
#' @param assign_niches Logical; if TRUE (default), also cluster cells by
#'   their neighborhood composition into \code{n_niches} groups and return
#'   the per-cell niche labels. If \code{FALSE}, only the enrichment stats
#'   are computed/returned.
#' @param n_niches Number of niche clusters to produce when
#'   \code{assign_niches = TRUE}. Default 5.
#' @param niche_method Clustering method: \code{"kmeans"} (default, scales
#'   to large datasets) or \code{"hclust"} (Ward's linkage on Euclidean
#'   distance; quadratic memory, suitable for < ~10k cells).
#' @param niche_prefix String prepended to each numeric cluster id to form
#'   the niche label. Default \code{"niche_"} (e.g. \code{"niche_3"}).
#' @param add_to_meta Logical; if TRUE (default) and \code{assign_niches =
#'   TRUE}, attach a copy of \code{obj} with the per-cell niche labels
#'   written to \code{obj@meta.data[[niche_col]]} as \code{$obj} in the
#'   returned list (see Details). Has no effect if \code{assign_niches =
#'   FALSE}.
#' @param niche_col Name of the metadata column to write niche labels to
#'   when \code{add_to_meta = TRUE}. Default \code{"niche"}.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{\code{z}}{Matrix of z-scores: rows = focal cell type,
#'       cols = neighbor cell type.}
#'     \item{\code{p}}{Matrix of two-sided empirical p-values.}
#'     \item{\code{padj}}{BH-adjusted q-values.}
#'     \item{\code{observed}}{Pooled observed neighbor-count matrix.}
#'     \item{\code{expected}}{Pooled permutation-mean count matrix.}
#'     \item{\code{results}}{Long-form data frame with columns
#'       \code{focal}, \code{neighbor}, \code{observed}, \code{expected},
#'       \code{z}, \code{p}, \code{padj}, sorted by \code{padj}.}
#'     \item{\code{niche}}{(only if \code{assign_niches = TRUE}) named
#'       character vector of niche labels, one per cell across all FOVs.}
#'     \item{\code{composition}}{(only if \code{assign_niches = TRUE})
#'       cells x cell-type matrix of neighborhood composition fractions
#'       (rows sum to ~1).}
#'     \item{\code{obj}}{(only if \code{assign_niches = TRUE} and
#'       \code{add_to_meta = TRUE}) a copy of the input Seurat object with
#'       the per-cell niche labels written to
#'       \code{obj@meta.data[[niche_col]]} (see Details).}
#'   }
#' @importFrom Seurat GetTissueCoordinates
#' @importFrom RANN nn2
#' @importFrom stats sd p.adjust kmeans dist hclust cutree
#' @export
NeighborhoodEnrichment <- function(obj,
                                   group.by,
                                   fovs          = NULL,
                                   k             = 10,
                                   n_perm        = 100,
                                   seed          = 1,
                                   assign_niches = TRUE,
                                   n_niches      = 5,
                                   niche_method  = c("kmeans", "hclust"),
                                   niche_prefix  = "niche_",
                                   add_to_meta   = TRUE,
                                   niche_col     = "niche") {

  niche_method <- match.arg(niche_method)
  .assert_seurat(obj)
  if (is.null(fovs)) fovs <- names(obj@images)
  if (!length(fovs)) stop("No FOVs in obj@images.")
  missing_fovs <- setdiff(fovs, names(obj@images))
  if (length(missing_fovs)) {
    stop("FOV(s) not in obj@images: ", paste(missing_fovs, collapse = ", "))
  }
  if (!group.by %in% colnames(obj@meta.data)) {
    stop("group.by '", group.by, "' not in metadata.")
  }
  if (isTRUE(assign_niches) && n_niches < 2) {
    stop("`n_niches` must be >= 2.")
  }

  # ---- Collect coords + labels per FOV (single pass) ---------------------
  per_fov <- vector("list", length(fovs))
  names(per_fov) <- fovs
  for (fov in fovs) {
    coords_mat <- .get_fov_coords(obj, fov, cells = rownames(obj@meta.data))
    cells_in_fov <- rownames(coords_mat)
    labels <- as.character(obj@meta.data[cells_in_fov, group.by, drop = TRUE])
    keep <- !is.na(labels)
    per_fov[[fov]] <- list(
      coords = coords_mat[keep, , drop = FALSE],
      labels = labels[keep]
    )
  }
  all_types <- sort(unique(unlist(lapply(per_fov, function(x) x$labels))))
  n_types   <- length(all_types)
  if (!n_types) stop("No non-NA labels across the requested FOVs.")

  empty_mat <- function() {
    matrix(0, nrow = n_types, ncol = n_types,
           dimnames = list(all_types, all_types))
  }
  # Vectorized cross-tab rather than one table()/factor() call per cell:
  # flattening `neighbor_idx_mat` (n cells x k neighbors) column-major and
  # repeating `focal_lab` once per neighbor column lines each neighbor slot
  # up with its own cell's focal label, so a single 2-D table() over all
  # (cell, neighbor-slot) pairs reproduces exactly what the old per-cell
  # loop accumulated -- just without the per-cell R-level overhead, which
  # matters here since this runs once per permutation (n_perm + 1 times
  # per FOV).
  count_table <- function(focal_lab, neighbor_idx_mat, label_vec) {
    focal_rep     <- rep(focal_lab, times = ncol(neighbor_idx_mat))
    neighbor_flat <- label_vec[as.vector(neighbor_idx_mat)]
    tab <- table(factor(focal_rep,     levels = all_types),
                factor(neighbor_flat, levels = all_types))
    matrix(as.numeric(tab), nrow = n_types, ncol = n_types,
           dimnames = list(all_types, all_types))
  }

  # ---- Observed and per-permutation counts, summed across FOVs -----------
  # When assign_niches = TRUE we also build a per-cell neighborhood
  # composition matrix (cells x cell-type), since the k-NN result we
  # compute for the test is exactly the data the clustering needs.
  set.seed(seed)
  obs <- empty_mat()
  perms_array <- array(0,
                       dim      = c(n_types, n_types, n_perm),
                       dimnames = list(all_types, all_types, NULL))
  per_cell_comp_list <- if (assign_niches) vector("list", length(fovs)) else NULL
  if (assign_niches) names(per_cell_comp_list) <- fovs

  for (fov in fovs) {
    coords <- per_fov[[fov]]$coords
    labels <- per_fov[[fov]]$labels
    if (nrow(coords) < 2) {
      message(sprintf("  '%s': < 2 cells, skipping", fov))
      next
    }
    k_use <- min(k, nrow(coords) - 1)
    nn <- RANN::nn2(coords, coords, k = k_use + 1)$nn.idx[, -1, drop = FALSE]
    obs <- obs + count_table(labels, nn, labels)
    for (p in seq_len(n_perm)) {
      shuf <- sample(labels)
      perms_array[, , p] <- perms_array[, , p] + count_table(shuf, nn, shuf)
    }

    # Per-cell neighborhood composition for clustering. Vectorized the same
    # way as count_table() above: reshape the looked-up neighbor labels back
    # into an n-cell x k-neighbor matrix, then count matches per type with
    # rowSums() (looping over the -- typically small -- number of types
    # rather than over every cell).
    if (assign_niches) {
      neighbor_lab_mat <- matrix(labels[as.vector(nn)], nrow = nrow(nn))
      comp <- matrix(0, nrow = nrow(coords), ncol = n_types,
                     dimnames = list(rownames(coords), all_types))
      for (t in seq_len(n_types)) {
        comp[, t] <- rowSums(neighbor_lab_mat == all_types[t])
      }
      per_cell_comp_list[[fov]] <- comp
    }
  }

  # ---- Test statistics ---------------------------------------------------
  perm_mean <- apply(perms_array, c(1, 2), mean)
  perm_sd   <- apply(perms_array, c(1, 2), stats::sd)
  perm_sd_safe <- perm_sd
  perm_sd_safe[perm_sd_safe == 0] <- 1e-9
  z <- (obs - perm_mean) / perm_sd_safe

  obs_expanded <- array(rep(obs, n_perm),
                        dim = c(n_types, n_types, n_perm))
  ge <- apply(perms_array >= obs_expanded, c(1, 2), sum)
  le <- apply(perms_array <= obs_expanded, c(1, 2), sum)
  p_up   <- (1 + ge) / (n_perm + 1)
  p_down <- (1 + le) / (n_perm + 1)
  p_two  <- pmin(2 * pmin(p_up, p_down), 1)

  padj <- matrix(stats::p.adjust(as.vector(p_two), method = "BH"),
                 nrow = n_types, ncol = n_types,
                 dimnames = dimnames(obs))

  pairs <- expand.grid(
    focal    = all_types,
    neighbor = all_types,
    stringsAsFactors = FALSE
  )
  pairs$observed <- as.vector(obs)
  pairs$expected <- as.vector(perm_mean)
  pairs$z        <- as.vector(z)
  pairs$p        <- as.vector(p_two)
  pairs$padj     <- as.vector(padj)
  pairs <- pairs[order(pairs$padj, -abs(pairs$z)), ]
  rownames(pairs) <- NULL

  out <- list(
    z        = z,
    p        = p_two,
    padj     = padj,
    observed = obs,
    expected = perm_mean,
    results  = pairs
  )

  # ---- Niche assignment (optional) ---------------------------------------
  if (assign_niches) {
    all_comp <- do.call(rbind, per_cell_comp_list)
    # Normalize each row to fractions so the clustering is over
    # composition, not size.
    row_n <- rowSums(all_comp)
    row_n[row_n == 0] <- 1
    comp_frac <- all_comp / row_n

    # Reset the seed here rather than relying on whatever RNG state the
    # permutation loop above left behind -- otherwise the niche clustering
    # (a logically independent step) would silently change whenever n_perm
    # changes, even with the same `seed`.
    set.seed(seed)
    if (niche_method == "kmeans") {
      km <- stats::kmeans(comp_frac, centers = n_niches, nstart = 25,
                          iter.max = 50)
      niche_int <- km$cluster
    } else {
      hc <- stats::hclust(stats::dist(comp_frac), method = "ward.D2")
      niche_int <- stats::cutree(hc, k = n_niches)
    }
    niche_labels <- paste0(niche_prefix, niche_int)
    names(niche_labels) <- rownames(comp_frac)
    out$niche       <- niche_labels
    out$composition <- comp_frac

    # ---- Optionally attach obj (with niche labels) to the return list ----
    # Cells outside the requested FOV(s), or with an NA `group.by` label
    # (and therefore no neighborhood composition), get NA here.
    if (isTRUE(add_to_meta)) {
      niche_col_vals <- rep(NA_character_, nrow(obj@meta.data))
      names(niche_col_vals) <- rownames(obj@meta.data)
      common <- intersect(names(niche_labels), names(niche_col_vals))
      niche_col_vals[common] <- niche_labels[common]
      obj@meta.data[[niche_col]] <- unname(niche_col_vals)
      out$obj <- obj
    }
  }

  out
}
