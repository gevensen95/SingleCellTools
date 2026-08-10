#' Summarize RCTD deconvolution proportions by cluster
#'
#' Companion to \code{\link{RunRCTD}}. \code{RunRCTD()} writes a per-spot
#' \code{rctd_dominant} column -- each spot's single highest-weight cell
#' type -- but collapsing every spot to one winner before looking at the
#' cluster level throws away exactly the minority signal RCTD's mixture
#' modeling was meant to preserve: average \code{rctd_dominant} (or
#' anything derived from it) across a cluster's spots and low-weight cell
#' types disappear before you ever see them.
#'
#' This instead averages the full \code{rctd_<celltype>} proportion
#' columns within each cluster first, then picks a per-cluster label from
#' that averaged composition -- the same threshold-based "best label, else
#' Unknown" logic \code{\link{AnnotateClusters}} uses internally, just
#' applied to RCTD's already mixture-aware proportions instead of
#' re-scoring raw expression. That sidesteps issues like an ambient/
#' dominant cell type's expression leaking into every cluster's marker
#' score (a real problem for \code{AnnotateClusters} on Visium data --
#' see its \code{visium_mode} docs) since RCTD's deconvolution has already
#' accounted for that when estimating proportions.
#'
#' @param obj A Seurat object with \code{rctd_<celltype>} columns in
#'   \code{obj@meta.data} (i.e. \code{\link{RunRCTD}} has already been run
#'   with \code{write_metadata = TRUE}).
#' @param cluster_col Metadata column holding cluster ids. Default
#'   \code{"seurat_clusters"}.
#' @param weight_cols Character vector of the proportion columns to
#'   aggregate. \code{NULL} (default) uses every \code{rctd_<celltype>}
#'   column, excluding the derived \code{rctd_dominant}/
#'   \code{rctd_max_weight}/\code{rctd_spot_class} columns
#'   \code{\link{RunRCTD}} also writes.
#' @param min_score Minimum average proportion below which a cluster is
#'   labeled \code{unassigned_label}. Default \code{NA} (disabled),
#'   matching \code{\link{AnnotateClusters}}' own default.
#' @param min_margin Minimum margin between the best and second-best
#'   average proportion; clusters below this margin are also labeled
#'   \code{unassigned_label}. Default \code{NA} (disabled).
#' @param unassigned_label Label applied to clusters that fail the
#'   thresholds. Default \code{"Unknown"}.
#' @return A list with elements \code{composition} (cluster x cell-type
#'   matrix of mean RCTD proportions) and \code{labels} (named character
#'   vector, cluster id -> assigned label). Nothing is written back to
#'   \code{obj} -- this returns values rather than mutating metadata, so
#'   attach the result to your object (or not) as you see fit.
#' @examples
#' \dontrun{
#' visium <- RunRCTD(visium, reference = ref, celltype_col = "cell_type")
#' res <- SummarizeRCTDByCluster(visium)
#' res$labels        # named vector: cluster -> most likely cell type
#' res$composition    # cluster x cell-type mean-proportion matrix
#'
#' # Flag ambiguous clusters instead of always returning a best guess
#' res <- SummarizeRCTDByCluster(visium, min_score = 0.3, min_margin = 0.1)
#' }
#' @export
SummarizeRCTDByCluster <- function(obj,
                                   cluster_col      = "seurat_clusters",
                                   weight_cols      = NULL,
                                   min_score        = NA,
                                   min_margin       = NA,
                                   unassigned_label = "Unknown") {

  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (!cluster_col %in% colnames(obj@meta.data)) {
    stop("Cluster column '", cluster_col, "' not found in metadata.")
  }

  if (is.null(weight_cols)) {
    weight_cols <- grep("^rctd_", colnames(obj@meta.data), value = TRUE)
    weight_cols <- setdiff(weight_cols,
                           c("rctd_dominant", "rctd_max_weight", "rctd_spot_class"))
    if (length(weight_cols) == 0) {
      stop("No `rctd_<celltype>` columns found in obj@meta.data and ",
           "`weight_cols` wasn't supplied. Run RunRCTD(write_metadata = TRUE) ",
           "first, or pass `weight_cols` explicitly.")
    }
  }
  missing_cols <- setdiff(weight_cols, colnames(obj@meta.data))
  if (length(missing_cols)) {
    stop("Column(s) not found in obj@meta.data: ",
         paste(missing_cols, collapse = ", "))
  }
  not_numeric <- weight_cols[!vapply(obj@meta.data[weight_cols], is.numeric,
                                     logical(1))]
  if (length(not_numeric)) {
    stop("`weight_cols` must be numeric columns; not numeric: ",
         paste(not_numeric, collapse = ", "))
  }

  clusters_vec <- as.character(obj@meta.data[[cluster_col]])
  w <- as.matrix(obj@meta.data[, weight_cols, drop = FALSE])

  # do.call(rbind, lapply(...)) rather than t(sapply(...)) -- same
  # single-weight-column simplification hazard as AnnotateClusters' matching
  # step (see the comment there): with only one weight column, sapply()
  # would collapse each cluster's colMeans() result down to a plain vector
  # before t() ever sees it, producing a 1-row matrix instead of one row per
  # cluster. rbind() of the per-cluster vectors is robust to this regardless
  # of column count. na.rm = TRUE so a spot with an NA proportion (e.g.
  # outside every image) doesn't wipe out an otherwise-informative cluster
  # mean.
  composition <- do.call(rbind, lapply(
    split(seq_len(nrow(w)), clusters_vec),
    function(idx) colMeans(w[idx, , drop = FALSE], na.rm = TRUE)
  ))
  rownames(composition) <- levels(factor(clusters_vec))
  cell_types <- sub("^rctd_", "", weight_cols)
  colnames(composition) <- cell_types

  labels <- .assign_with_unassigned(composition, cell_types, min_score,
                                    min_margin, unassigned_label)

  list(composition = composition, labels = labels)
}
