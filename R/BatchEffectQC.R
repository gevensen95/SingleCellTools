#' Quantify batch mixing / integration quality
#'
#' Computes several standard metrics that quantify how well integration
#' has mixed batches while preserving biological structure. Meant to be
#' run \emph{before} and \emph{after} integration to decide whether the
#' integration helped.
#'
#' Metrics returned:
#' \describe{
#'   \item{\code{batch_asw}}{Average silhouette width using \code{batch_col}
#'     as the grouping. Lower is better mixing (closer to 0). Values much
#'     greater than 0 indicate cells cluster by batch.}
#'   \item{\code{celltype_asw}}{Average silhouette width using
#'     \code{celltype_col} as the grouping. Higher is better biology
#'     preservation.}
#'   \item{\code{knn_mixing}}{Mean fraction of each cell's k-NN neighbors
#'     that come from a different batch than the cell itself. \code{1 -
#'     expected_same_batch_frac} is 1.0 for perfect mixing. We report the
#'     ratio observed / expected so 1.0 = perfectly mixed, < 1 = batchy.}
#'   \item{\code{knn_purity}}{Mean fraction of each cell's k-NN neighbors
#'     that share its \code{celltype_col} label. 1.0 = perfect biology
#'     preservation.}
#' }
#'
#' Per-cluster breakdowns are returned in \code{per_cluster} when
#' \code{celltype_col} is supplied.
#'
#' @param obj A Seurat object with a reduction to evaluate.
#' @param reduction Reduction name (usually the pre- or post-integration
#'   embedding, e.g. \code{"pca"} or \code{"harmony"}). Default \code{"pca"}.
#' @param dims Number of dimensions to use. Default 20.
#' @param batch_col Metadata column identifying batches / samples. Required.
#' @param celltype_col Optional metadata column with cell-type labels.
#' @param k Neighbors for kNN mixing / purity. Default 30.
#' @param n_cells_max Subsample this many cells before computing silhouette
#'   (silhouette is O(n^2)). Default 5000. Set to \code{Inf} to disable.
#' @param seed RNG seed for the subsample. Default 42.
#' @return A list with \code{summary} (named numeric vector) and, if
#'   \code{celltype_col} was supplied, \code{per_cluster} (data frame).
#' @examples
#' \dontrun{
#' # Compare unintegrated vs Harmony
#' pre  <- BatchEffectQC(obj, reduction = "pca",
#'                       batch_col = "orig.ident",
#'                       celltype_col = "predicted_cell_type")
#' post <- BatchEffectQC(obj, reduction = "harmony",
#'                       batch_col = "orig.ident",
#'                       celltype_col = "predicted_cell_type")
#' rbind(pre$summary, post$summary)
#' }
#' @importFrom Seurat Embeddings
#' @importFrom RANN nn2
#' @export
BatchEffectQC <- function(obj,
                          reduction    = "pca",
                          dims         = 20,
                          batch_col    = NULL,
                          celltype_col = NULL,
                          k            = 30,
                          n_cells_max  = 5000,
                          seed         = 42) {

  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (is.null(batch_col)) stop("`batch_col` is required.")
  if (!batch_col %in% colnames(obj@meta.data)) {
    stop("`batch_col` '", batch_col, "' not found in metadata.")
  }
  if (!(reduction %in% names(obj@reductions))) {
    stop("Reduction '", reduction, "' not found in obj@reductions.")
  }

  emb <- Seurat::Embeddings(obj, reduction = reduction)
  dims <- min(dims, ncol(emb))
  emb <- emb[, seq_len(dims), drop = FALSE]

  batch <- as.character(obj@meta.data[[batch_col]])
  celltype <- if (!is.null(celltype_col)) {
    as.character(obj@meta.data[[celltype_col]])
  } else NULL

  # kNN-based mixing / purity (no subsample needed; scales with n * k)
  nn <- RANN::nn2(data = emb, query = emb, k = k + 1)$nn.idx[, -1, drop = FALSE]
  mixing <- rowMeans(matrix(batch[nn], ncol = k) != batch)
  # Expected mixing if labels were random: 1 - sum(p_b^2)
  freqs <- as.numeric(table(batch)) / length(batch)
  expected_mixing <- 1 - sum(freqs ^ 2)
  knn_mixing_ratio <- mean(mixing) / max(1e-8, expected_mixing)

  purity <- NA_real_
  per_cluster <- NULL
  if (!is.null(celltype)) {
    purity <- rowMeans(matrix(celltype[nn], ncol = k) == celltype)
    per_cluster <- data.frame(
      cluster    = unique(celltype),
      n_cells    = as.integer(table(celltype)[unique(celltype)]),
      knn_purity = vapply(unique(celltype),
                          function(cl) mean(purity[celltype == cl]),
                          numeric(1)),
      knn_mixing = vapply(unique(celltype),
                          function(cl) mean(mixing[celltype == cl]) /
                                       max(1e-8, expected_mixing),
                          numeric(1)),
      stringsAsFactors = FALSE
    )
    per_cluster <- per_cluster[order(per_cluster$cluster), ]
    rownames(per_cluster) <- NULL
  }

  # Silhouette (subsample for O(n^2))
  batch_asw <- NA_real_
  celltype_asw <- NA_real_
  if (requireNamespace("cluster", quietly = TRUE)) {
    set.seed(seed)
    idx <- if (nrow(emb) > n_cells_max) {
      sample.int(nrow(emb), n_cells_max)
    } else {
      seq_len(nrow(emb))
    }
    sub_emb   <- emb[idx, , drop = FALSE]
    sub_batch <- as.integer(factor(batch[idx]))
    d <- stats::dist(sub_emb)
    sil_b <- cluster::silhouette(sub_batch, d)
    batch_asw <- mean(sil_b[, "sil_width"])

    if (!is.null(celltype)) {
      sub_ct <- as.integer(factor(celltype[idx]))
      # Silhouette needs >= 2 groups
      if (length(unique(sub_ct)) >= 2L) {
        sil_c <- cluster::silhouette(sub_ct, d)
        celltype_asw <- mean(sil_c[, "sil_width"])
      }
    }
  } else {
    message("Package 'cluster' not installed; skipping silhouette metrics.")
  }

  summary <- c(
    batch_asw        = batch_asw,
    celltype_asw     = celltype_asw,
    knn_mixing       = knn_mixing_ratio,
    knn_purity       = if (!is.null(celltype)) mean(purity) else NA_real_,
    expected_mixing  = expected_mixing,
    n_cells          = nrow(emb),
    n_batches        = length(unique(batch)),
    n_celltypes      = if (!is.null(celltype)) length(unique(celltype))
                       else NA_real_
  )

  list(summary = summary, per_cluster = per_cluster)
}
