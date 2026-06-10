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
#' to \code{obj@meta.data[[niche_col]]} \emph{on the object passed in by the
#' caller}: \code{obj} is updated in place in the calling environment as a
#' side effect, so after calling, e.g., \code{NeighborhoodEnrichment(combined,
#' ...)}, \code{combined} itself will have the new \code{niche_col} column,
#' separate from (and in addition to) the stats list returned by this
#' function. This requires \code{obj} to be passed as a plain variable name;
#' if it is not (e.g. it is the result of a function call or \code{$}/\code{[[}
#' expression), a warning is issued and the calling environment is left
#' untouched. Cells not covered by any requested FOV (or with an \code{NA}
#' \code{group.by} label) get \code{NA} in this column.
#'
#' \strong{P-value calculation.} Two-sided empirical p-value: for each pair
#' (i, j), let \eqn{p_{up} = (1 + \#\{perm_{ij} \ge obs_{ij}\}) /
#' (n\_perm + 1)} and \eqn{p_{down}} symmetrically. The reported p-value is
#' \eqn{2 \cdot \min(p_{up}, p_{down})}, capped at 1.
#'
#' @param obj A Seurat object with spatial coordinates.
#' @param fovs Character vector of FOV / image names. \code{NULL} (default)
#'   uses every image attached to the object.
#' @param group.by Metadata column with cell-type / cluster labels.
#' @param k Number of nearest neighbors per cell. Default 10.
#' @param n_perm Number of permutations for the null. Default 100.
#' @param seed Random seed.
#' @param assign_niches Logical; if TRUE, cluster cells by their
#'   neighborhood composition into \code{n_niches} groups and return the
#'   per-cell niche labels. Default \code{TRUE} (preserves the original
#'   return shape; only stats are computed).
#' @param n_niches Number of niche clusters to produce when
#'   \code{assign_niches = TRUE}. Default 5.
#' @param niche_method Clustering method: \code{"kmeans"} (default, scales
#'   to large datasets) or \code{"hclust"} (Ward's linkage on Euclidean
#'   distance; quadratic memory, suitable for < ~10k cells).
#' @param niche_prefix String prepended to each numeric cluster id to form
#'   the niche label. Default \code{"niche_"} (e.g. \code{"niche_3"}).
#' @param add_to_meta Logical; if TRUE (default) and \code{assign_niches =
#'   TRUE}, write the per-cell niche labels to
#'   \code{obj@meta.data[[niche_col]]} \emph{on the caller's copy of
#'   \code{obj}}, updating it in place in the calling environment (see
#'   Details). Has no effect if \code{assign_niches = FALSE}.
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
#'   }
#'   If \code{assign_niches = TRUE} and \code{add_to_meta = TRUE}, the
#'   variable passed in as \code{obj} is \emph{also} updated in place in the
#'   calling environment with a new \code{obj@meta.data[[niche_col]]} column
#'   (see Details) -- i.e. this function produces two separate, independently
#'   useful results: the returned stats list, and the modified Seurat object
#'   under its original variable name.
#' @importFrom Seurat GetTissueCoordinates
#' @importFrom RANN nn2
#' @importFrom stats sd p.adjust kmeans dist hclust cutree
#' @export
NeighborhoodEnrichment <- function(obj,
                                   fovs          = NULL,
                                   group.by,
                                   k             = 10,
                                   n_perm        = 100,
                                   seed          = 1,
                                   assign_niches = TRUE,
                                   n_niches      = 5,
                                   niche_method  = c("kmeans", "hclust"),
                                   niche_prefix  = "niche_",
                                   add_to_meta   = TRUE,
                                   niche_col     = "niche") {

  # Captured before any argument matching/promises change `obj` itself, so
  # that (if it's a plain variable name) we can write the updated object
  # back to that same variable in the caller's environment below.
  obj_expr <- substitute(obj)

  niche_method <- match.arg(niche_method)
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
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
    coords <- Seurat::GetTissueCoordinates(obj[[fov]], which = "centroids")
    if ("cell" %in% colnames(coords)) {
      rownames(coords) <- coords$cell
      coords_mat <- as.matrix(coords[, c("x", "y")])
    } else {
      coords_mat <- as.matrix(coords[, c("x", "y")])
    }
    cells_in_fov <- intersect(rownames(coords_mat), rownames(obj@meta.data))
    coords_mat <- coords_mat[cells_in_fov, , drop = FALSE]
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
  count_table <- function(focal_lab, neighbor_idx_mat, label_vec) {
    m <- empty_mat()
    for (i in seq_along(focal_lab)) {
      tab <- table(factor(label_vec[neighbor_idx_mat[i, ]], levels = all_types))
      m[focal_lab[i], ] <- m[focal_lab[i], ] + tab
    }
    m
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

    # Per-cell neighborhood composition for clustering
    if (assign_niches) {
      comp <- matrix(0, nrow = nrow(coords), ncol = n_types,
                     dimnames = list(rownames(coords), all_types))
      for (i in seq_len(nrow(nn))) {
        tab <- table(factor(labels[nn[i, ]], levels = all_types))
        comp[i, ] <- as.numeric(tab)
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

    # ---- Optionally write the niche labels back onto `obj` ---------------
    # Cells outside the requested FOV(s), or with an NA `group.by` label
    # (and therefore no neighborhood composition), get NA here.
    if (isTRUE(add_to_meta)) {
      niche_col_vals <- rep(NA_character_, nrow(obj@meta.data))
      names(niche_col_vals) <- rownames(obj@meta.data)
      common <- intersect(names(niche_labels), names(niche_col_vals))
      niche_col_vals[common] <- niche_labels[common]
      obj@meta.data[[niche_col]] <- unname(niche_col_vals)

      # Write the updated object back to the variable the caller passed in
      # for `obj`, so it ends up alongside `out` as a second, independent
      # object in the caller's environment. Only possible if `obj` was
      # passed as a plain variable name (not e.g. `my_list$obj` or the
      # result of a function call).
      if (is.symbol(obj_expr)) {
        assign(deparse(obj_expr), obj, envir = parent.frame())
      } else {
        warning("`add_to_meta = TRUE` but `obj` was not passed as a plain ",
                "variable name, so the updated object with the '", niche_col,
                "' column could not be written back to the calling ",
                "environment. Pass a variable (e.g. ",
                "`NeighborhoodEnrichment(my_obj, ...)`), or set ",
                "`add_to_meta = FALSE` and add the niche column yourself ",
                "from `$niche` in the returned list.")
      }
    }
  }

  out
}
