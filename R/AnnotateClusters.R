#' Annotate clusters with cell-type labels
#'
#' Two annotation modes:
#' \describe{
#'   \item{\code{method = "marker"} (default)}{Score each cell on each cell
#'     type's marker set with \code{UCell::AddModuleScore_UCell}, average the
#'     scores per cluster, and assign each cluster the cell type with the
#'     highest mean score. Requires \code{UCell}.}
#'   \item{\code{method = "singler"}}{Run SingleR against a reference
#'     SummarizedExperiment / Seurat object and propagate per-cell labels up
#'     to a per-cluster majority vote. Requires \code{SingleR}.}
#' }
#'
#' \strong{Unknown / unassigned clusters.} Two thresholds:
#' \itemize{
#'   \item \code{min_score}: clusters whose best score (marker mode) or
#'     majority-vote fraction (SingleR mode) is below this value get the
#'     \code{unassigned_label}.
#'   \item \code{min_margin}: clusters whose best score doesn't beat the
#'     second-best by at least this margin are flagged as ambiguous and
#'     also get the \code{unassigned_label}.
#' }
#'
#' By default both thresholds are \code{NA}, which disables the
#' "Unknown" checks entirely — every cluster gets the top-scoring label,
#' regardless of how weak or ambiguous that score is. To re-enable the
#' checks, pass an explicit numeric value, or pass \code{NULL} to fall
#' back to the mode-specific defaults:
#' \itemize{
#'   \item \strong{Marker mode}: \code{min_score = 0.1},
#'     \code{min_margin = 0.05} (UCell scores are bounded in \code{(0, 1)}
#'     and typically peak around 0.15-0.35 for matched signatures).
#'   \item \strong{SingleR mode}: \code{min_score = 0.5},
#'     \code{min_margin = 0.2} (require a majority vote that beats the
#'     runner-up by at least 20 percentage points).
#' }
#'
#' @param obj A Seurat object with clusters already computed.
#' @param method \code{"marker"} or \code{"singler"}.
#' @param markers Named list of marker gene vectors, one element per cell
#'   type. Required when \code{method = "marker"}.
#' @param reference A reference for SingleR (SummarizedExperiment or
#'   Seurat object). Required when \code{method = "singler"}.
#' @param ref_label_col For \code{method = "singler"}: the column in
#'   the reference's colData / metadata holding cell-type labels.
#' @param cluster_col Metadata column holding the cluster ids to annotate.
#'   Default \code{"seurat_clusters"}.
#' @param assay Assay to use. Default DefaultAssay.
#' @param new_col Name of the new metadata column. Default \code{"cell_type"}.
#' @param min_score Minimum score below which a cluster is labeled
#'   \code{unassigned_label}. Default \code{NA} (check disabled).
#'   \code{NULL} uses the mode-specific default; an explicit number
#'   overrides it.
#' @param min_margin Minimum margin between best and 2nd-best score / vote.
#'   Default \code{NA} (check disabled). \code{NULL} uses the
#'   mode-specific default; an explicit number overrides it.
#' @param unassigned_label Label applied to clusters that fail the
#'   \code{min_score} or \code{min_margin} test. Default \code{"Unknown"}.
#' @return The Seurat object with a new metadata column \code{new_col}.
#' @importFrom Seurat DefaultAssay
#' @export
AnnotateClusters <- function(obj,
                             method           = c("marker", "singler"),
                             markers          = NULL,
                             reference        = NULL,
                             ref_label_col    = "label.main",
                             cluster_col      = "seurat_clusters",
                             assay            = NULL,
                             new_col          = "cell_type",
                             min_score        = NA,
                             min_margin       = NA,
                             unassigned_label = "Unknown") {

  method <- match.arg(method)
  if (!cluster_col %in% colnames(obj@meta.data)) {
    stop("Cluster column '", cluster_col, "' not found in metadata.")
  }
  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay

  # ---- Resolve NULL thresholds to mode-specific defaults -----------------
  # (NA is now the default and disables the check; NULL opts back into the
  # mode-specific defaults below.)
  if (is.null(min_score)) {
    min_score <- if (method == "marker") 0.10 else 0.50
    message(sprintf("Using default min_score = %.2f for method '%s' ",
                    min_score, method),
            "(pass NA to disable)")
  }
  if (is.null(min_margin)) {
    min_margin <- if (method == "marker") 0.05 else 0.20
    message(sprintf("Using default min_margin = %.2f for method '%s' ",
                    min_margin, method),
            "(pass NA to disable)")
  }

  # Helper: apply unassigned rules to a score matrix (rows = clusters,
  # cols = candidate labels). Returns a named character vector of labels.
  .assign_with_unassigned <- function(score_mat, candidate_labels) {
    best_idx   <- apply(score_mat, 1, which.max)
    best_score <- apply(score_mat, 1, max)
    sorted_per_row <- apply(score_mat, 1, sort, decreasing = TRUE)
    second_score <- if (ncol(score_mat) >= 2) sorted_per_row[2, ] else 0
    margins <- best_score - second_score
    labs <- candidate_labels[best_idx]
    if (!is.na(min_score)) {
      labs[best_score < min_score] <- unassigned_label
    }
    if (!is.na(min_margin)) {
      labs[margins < min_margin] <- unassigned_label
    }
    setNames(labs, rownames(score_mat))
  }

  if (method == "marker") {
    if (!requireNamespace("UCell", quietly = TRUE)) {
      stop("'UCell' is required for marker-based annotation.")
    }
    if (is.null(markers) || !is.list(markers) || is.null(names(markers))) {
      stop("`markers` must be a named list (cell_type -> character vector of genes).")
    }

    # ---- Prefilter marker genes to those present in the object ------------
    # Genes absent from the assay (species/symbol mismatches, filtered-out
    # genes, etc.) silently get a near-zero UCell score for every cell type
    # that includes them, which can drag whole signatures below min_score.
    # Drop them up front and report what was removed.
    available_genes <- rownames(obj[[a]])
    filtered_markers <- lapply(markers, function(g) intersect(g, available_genes))
    n_orig    <- vapply(markers, length, integer(1))
    n_dropped <- n_orig - vapply(filtered_markers, length, integer(1))
    if (any(n_dropped > 0)) {
      for (ct in names(markers)[n_dropped > 0]) {
        dropped_genes <- setdiff(markers[[ct]], filtered_markers[[ct]])
        message(sprintf(
          "  '%s': dropping %d/%d marker gene(s) not found in assay '%s': %s",
          ct, n_dropped[ct], n_orig[ct], a,
          paste(dropped_genes, collapse = ", ")))
      }
    }
    empty <- names(filtered_markers)[vapply(filtered_markers, length, integer(1)) == 0]
    if (length(empty) > 0) {
      stop("The following cell type(s) have no marker genes present in assay '",
           a, "': ", paste(empty, collapse = ", "),
           ". Check gene symbol case/format and species (e.g. 'Cd3e' vs 'CD3E').")
    }
    markers <- filtered_markers

    message(sprintf("--- Scoring %d cell-type signatures with UCell ---",
                    length(markers)))
    obj <- UCell::AddModuleScore_UCell(obj, features = markers, assay = a,
                                       name = "_score")
    score_cols <- paste0(names(markers), "_score")
    if (!all(score_cols %in% colnames(obj@meta.data))) {
      stop("UCell did not produce expected columns; check your marker list names.")
    }

    # Average score per cluster, assign max
    clusters <- as.character(obj@meta.data[[cluster_col]])
    score_mat <- as.matrix(obj@meta.data[, score_cols, drop = FALSE])
    avg_per_cluster <- t(sapply(split(seq_len(nrow(score_mat)), clusters),
                                function(idx) colMeans(score_mat[idx, , drop = FALSE])))
    rownames(avg_per_cluster) <- levels(factor(clusters))
    colnames(avg_per_cluster) <- names(markers)

    lookup <- .assign_with_unassigned(avg_per_cluster, names(markers))
    obj[[new_col]] <- unname(lookup[clusters])
    return(obj)
  }

  # SingleR mode
  if (!requireNamespace("SingleR", quietly = TRUE)) {
    stop("'SingleR' is required for SingleR-based annotation.")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("'SummarizedExperiment' is required for SingleR-based annotation.")
  }
  if (is.null(reference)) stop("`reference` is required for SingleR mode.")

  message("--- Running SingleR ---")
  test_mat <- SeuratObject::LayerData(obj, assay = a, layer = "data")
  if (inherits(reference, "Seurat")) {
    ref_mat    <- SeuratObject::LayerData(reference,
                                          assay = SeuratObject::DefaultAssay(reference),
                                          layer = "data")
    ref_labels <- reference@meta.data[[ref_label_col]]
  } else {
    ref_mat    <- SummarizedExperiment::assay(reference, "logcounts")
    ref_labels <- SummarizedExperiment::colData(reference)[[ref_label_col]]
  }
  if (is.null(ref_labels)) {
    stop("Reference column '", ref_label_col, "' not found.")
  }

  # ---- Prefilter to genes shared between object and reference -------------
  # SingleR scores correlations over the intersection of gene sets anyway,
  # but doing it explicitly here lets us catch the "near-zero overlap"
  # failure mode early (species/symbol mismatch) instead of silently getting
  # poor correlations and lots of low-confidence/"Unknown" calls.
  common_genes <- intersect(rownames(test_mat), rownames(ref_mat))
  if (length(common_genes) == 0) {
    stop("No genes in common between object (", nrow(test_mat),
         " genes) and reference (", nrow(ref_mat), " genes). ",
         "Check gene symbol case/format and species (e.g. 'Cd3e' vs 'CD3E').")
  }
  if (length(common_genes) < nrow(test_mat) || length(common_genes) < nrow(ref_mat)) {
    message(sprintf(
      "  Using %d gene(s) shared between object (%d) and reference (%d)",
      length(common_genes), nrow(test_mat), nrow(ref_mat)))
  }
  if (length(common_genes) < 100) {
    warning("Only ", length(common_genes), " gene(s) overlap between object ",
            "and reference; SingleR results may be unreliable. Check gene ",
            "symbol case/format and species.")
  }
  test_mat <- test_mat[common_genes, , drop = FALSE]
  ref_mat  <- ref_mat[common_genes, , drop = FALSE]

  pred <- SingleR::SingleR(test = test_mat, ref = ref_mat, labels = ref_labels)

  # Per-cluster vote fractions; then apply min_score / min_margin on fractions
  per_cell_labels <- pred$labels
  names(per_cell_labels) <- colnames(test_mat)
  clusters <- as.character(obj@meta.data[[cluster_col]])
  uniq_labels <- sort(unique(per_cell_labels))
  vote_frac <- t(sapply(split(seq_along(clusters), clusters), function(idx) {
    tab <- table(factor(per_cell_labels[idx], levels = uniq_labels))
    as.numeric(tab) / max(1, sum(tab))
  }))
  rownames(vote_frac) <- levels(factor(clusters))
  colnames(vote_frac) <- uniq_labels

  lookup <- .assign_with_unassigned(vote_frac, uniq_labels)
  obj[[new_col]] <- unname(lookup[clusters])
  obj
}
