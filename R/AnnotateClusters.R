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
#' \strong{Non-specific marker filtering (marker mode).} In tissues with
#' one dominant cell type (e.g. hepatocytes in liver), highly-expressed
#' "markers" of that type often leak into every other cluster — enough
#' that the dominant signature beats genuinely correct but
#' lower-magnitude signatures in non-dominant clusters. Set
#' \code{filter_nonspecific = TRUE} (default) to compute per-cluster mean
#' expression of every marker gene before scoring, identify genes that
#' are "high" in too many clusters, and drop them from the marker lists.
#' A gene is "high" in a cluster if its mean expression there exceeds
#' \code{high_relative_to_max} times its maximum cluster-mean across all
#' clusters; a gene is "non-specific" overall if it is "high" in more
#' than \code{nonspecific_max_fraction} of the clusters. The dropped
#' genes are reported in a message.
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
#' "Unknown" checks entirely. To re-enable, pass an explicit numeric
#' value, or pass \code{NULL} to fall back to mode-specific defaults
#' (marker: \code{min_score = 0.1, min_margin = 0.05}; SingleR:
#' \code{0.5, 0.2}).
#'
#' @param obj A Seurat object with clusters already computed.
#' @param method \code{"marker"} or \code{"singler"}.
#' @param markers Named list of marker gene vectors, one element per cell
#'   type. Required when \code{method = "marker"}.
#' @param reference A reference for SingleR. Required when
#'   \code{method = "singler"}.
#' @param ref_label_col For SingleR: the column in the reference's
#'   colData / metadata holding cell-type labels.
#' @param cluster_col Metadata column holding the cluster ids to annotate.
#'   Default \code{"seurat_clusters"}.
#' @param assay Assay to use. Default DefaultAssay.
#' @param new_col Name of the new metadata column. Default
#'   \code{"predicted_cell_type"}.
#' @param filter_nonspecific Logical; if TRUE (default) drop marker genes
#'   that are broadly expressed across many clusters before scoring. Only
#'   used in marker mode.
#' @param nonspecific_max_fraction Fraction of clusters above which a
#'   marker is considered non-specific. Default 0.5 — a gene that is
#'   "high" in more than half of the clusters is dropped.
#' @param high_relative_to_max A gene counts as "high" in a cluster if
#'   its mean expression there is at least this fraction of the gene's
#'   maximum cluster-mean across all clusters. Default 0.5 — within 2x
#'   of the gene's max cluster mean. Lower values (e.g. 0.3) are more
#'   permissive about what counts as "high" and so will drop more
#'   markers; higher values (e.g. 0.7) are stricter.
#' @param min_score Minimum score below which a cluster is labeled
#'   \code{unassigned_label}. Default \code{NA} (disabled). \code{NULL}
#'   uses the mode-specific default.
#' @param min_margin Minimum margin between best and 2nd-best score / vote.
#'   Default \code{NA} (disabled). \code{NULL} uses the mode-specific default.
#' @param unassigned_label Label applied to clusters that fail the
#'   thresholds. Default \code{"Unknown"}.
#' @return The Seurat object with a new metadata column \code{new_col}.
#' @importFrom Seurat DefaultAssay FetchData
#' @export
AnnotateClusters <- function(obj,
                             method                   = c("marker", "singler"),
                             markers                  = NULL,
                             reference                = NULL,
                             ref_label_col            = "label.main",
                             cluster_col              = "seurat_clusters",
                             assay                    = NULL,
                             new_col                  = "predicted_cell_type",
                             filter_nonspecific       = TRUE,
                             nonspecific_max_fraction = 0.5,
                             high_relative_to_max     = 0.5,
                             min_score                = NA,
                             min_margin               = NA,
                             unassigned_label         = "Unknown") {

  method <- match.arg(method)
  if (!cluster_col %in% colnames(obj@meta.data)) {
    stop("Cluster column '", cluster_col, "' not found in metadata.")
  }
  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay

  # ---- Resolve NULL thresholds to mode-specific defaults -----------------
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

    # ---- Prefilter: drop genes not in the assay ---------------------------
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
      warning("Dropping cell type(s) with no marker genes present in assay '",
              a, "': ", paste(empty, collapse = ", "),
              ". Check gene symbol case/format and species (e.g. 'Cd3e' vs 'CD3E'). ",
              "Remaining cell types will still be scored; dropped types cannot be ",
              "assigned and clusters that would have matched them will fall to the ",
              "next-best label (or '", unassigned_label, "' if thresholds are set).",
              call. = FALSE)
      filtered_markers <- filtered_markers[setdiff(names(filtered_markers), empty)]
    }
    if (length(filtered_markers) == 0) {
      stop("No cell types have any marker genes present in assay '", a,
           "'. Cannot run marker-based annotation. Check gene symbol case/format ",
           "and species (e.g. 'Cd3e' vs 'CD3E').")
    }
    markers <- filtered_markers

    # ---- Filter non-specific (broadly-expressed) marker genes -------------
    # Pre-compute per-cluster mean expression of every marker gene. A gene
    # whose mean is comparable across many clusters isn't really a marker
    # — it's a broadly-expressed gene that will contaminate signature
    # scores. Drop those genes before they reach UCell.
    if (isTRUE(filter_nonspecific)) {
      all_marker_genes <- unique(unlist(markers))
      clusters_vec     <- as.character(obj@meta.data[[cluster_col]])

      expr_df <- tryCatch(
        Seurat::FetchData(obj, vars = all_marker_genes,
                          layer = "data"),
        error = function(e) NULL
      )

      if (is.null(expr_df)) {
        message("  filter_nonspecific: FetchData failed; skipping filter.")
      } else {
        # Cluster x gene matrix of mean expression
        cluster_means <- t(sapply(
          split(seq_along(clusters_vec), clusters_vec),
          function(idx) colMeans(expr_df[idx, , drop = FALSE])
        ))
        rownames(cluster_means) <- levels(factor(clusters_vec))

        n_clusters <- nrow(cluster_means)
        is_nonspecific <- vapply(colnames(cluster_means), function(g) {
          g_means <- cluster_means[, g]
          mx <- max(g_means)
          if (mx <= 0) return(TRUE)              # absent everywhere
          frac_high <- sum(g_means >= high_relative_to_max * mx) / n_clusters
          frac_high > nonspecific_max_fraction
        }, logical(1))

        nonspecific_genes <- names(is_nonspecific)[is_nonspecific]
        if (length(nonspecific_genes) > 0) {
          message(sprintf(
            "  Dropping %d non-specific marker gene(s) (high in > %.0f%% of %d clusters): %s",
            length(nonspecific_genes),
            100 * nonspecific_max_fraction,
            n_clusters,
            paste(nonspecific_genes, collapse = ", ")))

          # Tell the user which signatures lost which genes
          for (ct in names(markers)) {
            removed_here <- intersect(markers[[ct]], nonspecific_genes)
            if (length(removed_here)) {
              message(sprintf("    '%s' lost: %s", ct,
                              paste(removed_here, collapse = ", ")))
            }
          }
          markers <- lapply(markers,
                            function(g) setdiff(g, nonspecific_genes))
        } else {
          message("  filter_nonspecific: no non-specific markers detected.")
        }

        # Check no signature has been wiped out
        empty <- names(markers)[vapply(markers, length, integer(1)) == 0]
        if (length(empty) > 0) {
          stop("After non-specific filtering, signature(s) have no markers ",
               "left: ", paste(empty, collapse = ", "),
               ". Loosen the filter (e.g. nonspecific_max_fraction = 0.7) ",
               "or provide more cluster-specific markers.")
        }
      }
    }

    # ---- Score with UCell -------------------------------------------------
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

  # SingleR mode (unchanged) ------------------------------------------------
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

  common_genes <- intersect(rownames(test_mat), rownames(ref_mat))
  if (length(common_genes) == 0) {
    stop("No genes in common between object (", nrow(test_mat),
         " genes) and reference (", nrow(ref_mat), " genes). ",
         "Check gene symbol case/format and species.")
  }
  if (length(common_genes) < nrow(test_mat) || length(common_genes) < nrow(ref_mat)) {
    message(sprintf(
      "  Using %d gene(s) shared between object (%d) and reference (%d)",
      length(common_genes), nrow(test_mat), nrow(ref_mat)))
  }
  if (length(common_genes) < 100) {
    warning("Only ", length(common_genes), " gene(s) overlap between object ",
            "and reference; SingleR results may be unreliable.")
  }
  test_mat <- test_mat[common_genes, , drop = FALSE]
  ref_mat  <- ref_mat[common_genes, , drop = FALSE]

  pred <- SingleR::SingleR(test = test_mat, ref = ref_mat, labels = ref_labels)

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
