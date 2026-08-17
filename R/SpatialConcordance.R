#' Spatial concordance of a categorical label, per sample
#'
#' A per-sample QC/validity check for a categorical spatial assignment
#' (zonation call, spatial domain, niche label, etc.): for each cell/spot,
#' what fraction of its \code{k} nearest spatial neighbors share its own
#' \code{group.by} label, averaged across all cells and compared to a
#' permutation null (the same labels randomly reshuffled among the same
#' spatial neighbors). A real, spatially coherent domain assignment should
#' score well above the null; a noisy or over-fragmented one won't.
#'
#' @details
#' This answers a different question than \code{\link{NeighborhoodEnrichment}}.
#' \code{NeighborhoodEnrichment()} pools cells across every requested image
#' and asks \emph{which pairs} of labels are enriched near each other
#' (a focal x neighbor label z-score matrix). \code{SpatialConcordance()}
#' asks a single, coarser question \emph{per sample}: is this one label
#' column spatially coherent at all? That makes it a natural post-hoc check
#' after calling zones/domains/niches on each sample independently, before
#' trusting the labels downstream.
#'
#' Neighbors are computed separately within each image (a cell is never
#' considered a neighbor of a cell in a different image), matching
#' \code{\link{NeighborhoodEnrichment}}'s convention. The null distribution
#' for a given image is built by permuting that image's labels among its own
#' cells only, so a sample with a genuinely small number of distinct labels
#' won't be penalized (or credited) using another sample's label
#' distribution.
#'
#' \strong{P-value calculation.} Two-sided empirical p-value: let
#' \eqn{p_{up} = (1 + \#\{null \ge obs\}) / (n\_perm + 1)} and
#' \eqn{p_{down}} symmetrically; \eqn{p = 2 \cdot \min(p_{up}, p_{down})},
#' capped at 1 -- the same formula \code{\link{NeighborhoodEnrichment}}
#' uses. \code{padj} is the BH adjustment of \code{p} across samples.
#'
#' @param obj A Seurat object with spatial coordinates (Visium or
#'   FOV-based -- Xenium/CosMx/MERFISH).
#' @param group.by Metadata column with the categorical label to test
#'   (a zonation call, spatial domain, niche assignment, etc.).
#' @param images Character vector of image names in \code{obj@images} to
#'   test. Defaults to all of them (one row per image in the output).
#' @param k Number of nearest spatial neighbors per cell. Default \code{6}.
#'   Capped at \code{n_cells - 1} for a sample with fewer usable cells than
#'   that.
#' @param n_perm Number of label permutations for the null. Default
#'   \code{100}.
#' @param exclude Label value(s) to drop before computing concordance --
#'   e.g. an "Unclassified"/"Unknown" catch-all bucket that shouldn't count
#'   as spatially meaningful. Default \code{NA} (drops only missing labels;
#'   pass e.g. \code{c(NA, "Unclassified")} to also drop a named bucket).
#' @param seed Random seed for the permutations. Default \code{1}.
#' @return A data frame with one row per image: \code{sample}, \code{n_spots}
#'   (cells retained after \code{exclude} filtering), \code{k_used},
#'   \code{observed} (mean same-label neighbor fraction), \code{null_mean},
#'   \code{null_sd}, \code{z}, \code{p}, \code{padj}. Images with fewer than
#'   2 usable cells are skipped with a message.
#' @examples
#' \dontrun{
#' SpatialConcordance(visium, group.by = "spatial_domain")
#' SpatialConcordance(visium, group.by = "Zone_final",
#'                    exclude = c(NA, "Unclassified"), k = 6)
#' }
#' @importFrom stats sd p.adjust
#' @importFrom dplyr bind_rows
#' @export
SpatialConcordance <- function(obj,
                               group.by,
                               images  = NULL,
                               k       = 6,
                               n_perm  = 100,
                               exclude = NA,
                               seed    = 1) {
  .assert_seurat(obj)
  if (length(obj@images) == 0) {
    stop("`obj` has no images in obj@images -- SpatialConcordance() needs ",
         "spatial coordinates.")
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
  if (!group.by %in% colnames(obj@meta.data)) {
    stop("`group.by` column '", group.by, "' not found in obj@meta.data.")
  }
  if (!is.numeric(k) || length(k) != 1 || k < 1) {
    stop("`k` must be a single positive number.")
  }
  if (!is.numeric(n_perm) || length(n_perm) != 1 || n_perm < 1) {
    stop("`n_perm` must be a single positive number.")
  }

  # Vectorized same-label neighbor fraction: for each cell, what share of
  # its k nearest neighbors carry its own label. Mirrors the
  # neighbor_lab_mat idiom NeighborhoodEnrichment() uses for its per-cell
  # composition matrix, rather than a per-cell vapply loop.
  .concordance <- function(labels_vec, nn_idx) {
    neighbor_lab_mat <- matrix(labels_vec[as.vector(nn_idx)], nrow = nrow(nn_idx))
    same <- neighbor_lab_mat == labels_vec
    mean(rowMeans(same))
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  out <- vector("list", length(images))
  names(out) <- images

  for (img in images) {
    coords <- .get_fov_coords(obj, img, cells = rownames(obj@meta.data))
    cells_in_img <- rownames(coords)
    labels <- as.character(obj@meta.data[cells_in_img, group.by, drop = TRUE])

    keep <- rep(TRUE, length(labels))
    if (anyNA(exclude)) keep <- keep & !is.na(labels)
    non_na_exclude <- exclude[!is.na(exclude)]
    if (length(non_na_exclude)) keep <- keep & !(labels %in% non_na_exclude)

    coords <- coords[keep, , drop = FALSE]
    labels <- labels[keep]

    if (length(labels) < 2) {
      message(sprintf("  '%s': < 2 usable cells after `exclude` filtering, skipping", img))
      next
    }

    k_use <- min(k, length(labels) - 1)
    nn <- RANN::nn2(coords, coords, k = k_use + 1)$nn.idx[, -1, drop = FALSE]

    observed  <- .concordance(labels, nn)
    null_dist <- replicate(n_perm, .concordance(sample(labels), nn))
    null_mean <- mean(null_dist)
    null_sd   <- stats::sd(null_dist)
    null_sd_safe <- if (null_sd == 0) 1e-9 else null_sd
    z <- (observed - null_mean) / null_sd_safe

    p_up   <- (1 + sum(null_dist >= observed)) / (n_perm + 1)
    p_down <- (1 + sum(null_dist <= observed)) / (n_perm + 1)
    p_two  <- min(2 * min(p_up, p_down), 1)

    out[[img]] <- data.frame(
      sample    = img,
      n_spots   = length(labels),
      k_used    = k_use,
      observed  = observed,
      null_mean = null_mean,
      null_sd   = null_sd,
      z         = z,
      p         = p_two,
      stringsAsFactors = FALSE
    )
  }

  res <- dplyr::bind_rows(out)
  if (nrow(res)) res$padj <- stats::p.adjust(res$p, method = "BH")
  res
}
