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
#' \strong{Tumor mode.} \code{tumor_mode = TRUE} adjusts defaults for
#' malignant samples: it disables \code{filter_nonspecific} (which is
#' counterproductive when tumor cells drive up per-gene maxima and make
#' real markers look non-specific) and tightens \code{high_relative_to_max}
#' if you do leave the filter on. It only changes \emph{defaults} — any
#' argument you pass explicitly still wins. Best practice for tumor data:
#' identify malignant cells first with InferCNV / CopyKAT / SCEVAN, then
#' annotate the immune/stromal compartment with this function and run
#' subtype scoring (proliferation, EMT, hypoxia, etc.) on the malignant
#' subset separately. Tumor-mode also issues a warning that labels are
#' approximate.
#'
#' \strong{Visium / spatial data.} This function is a winner-takes-all
#' classifier; Visium spots are 1-10 cells of mixed types, so per-spot
#' labels lose minority populations by construction. Use a dedicated
#' deconvolution tool -- \code{\link{RunRCTD}} in this package (wraps
#' \code{RCTD}/\code{spacexr}), or \code{cell2location} / \code{CARD} /
#' \code{SpotLight} -- to estimate per-spot cell-type \emph{proportions}
#' from a reference single-cell dataset; that's the recommended primary
#' annotation strategy for Visium, not this function.
#'
#' If you do still want a cluster-level descriptive label on Visium (e.g.
#' no suitable single-cell reference for RCTD, or alongside RCTD's
#' proportions rather than instead of them), \code{visium_mode = TRUE}
#' applies the fallback recipe this section used to just describe in
#' prose: disables \code{filter_nonspecific} by default (like
#' \code{tumor_mode}, only changes defaults -- explicit args still win)
#' and warns with the same reminders (bulk up marker sets to 30-50+
#' genes, consider \code{return_scores = "cluster"} to inspect minority
#' signal directly, prefer \code{\link{RunRCTD}}). Applies regardless of
#' \code{method} -- the mixed-spot caveat is equally true for
#' \code{"singler"} -- but only \code{filter_nonspecific} actually changes,
#' since marker mode is the only mode with that knob.
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
#' @param tumor_mode Logical; if TRUE, applies tumor-friendly defaults
#'   (\code{filter_nonspecific = FALSE}, \code{high_relative_to_max = 0.7})
#'   and issues a warning that labels are approximate. Defaults only —
#'   explicit args still win. Default FALSE.
#' @param visium_mode Logical; if TRUE, applies the Visium fallback recipe
#'   described above (\code{filter_nonspecific = FALSE} by default, marker
#'   mode only) and warns pointing at \code{\link{RunRCTD}} as the
#'   preferred alternative. Defaults only — explicit args still win.
#'   Default FALSE.
#' @param filter_nonspecific Logical; if TRUE drop marker genes that are
#'   broadly expressed across many clusters before scoring. Default TRUE
#'   (or FALSE if \code{tumor_mode = TRUE} and not set explicitly).
#'   Marker mode only.
#' @param nonspecific_max_fraction Fraction of clusters above which a
#'   marker is considered non-specific. Default 0.5.
#' @param high_relative_to_max A gene counts as "high" in a cluster if
#'   its mean expression there is at least this fraction of the gene's
#'   maximum cluster-mean across all clusters. Default 0.5 (0.7 in
#'   tumor mode unless overridden).
#' @param min_detection_frac Minimum fraction of cells/spots in which at
#'   least one of a signature's marker genes is detected (count > 0).
#'   Signatures below this threshold are dropped with a warning — useful
#'   on sparse Visium data or when a marker list is for the wrong tissue.
#'   Default \code{NA} (disabled). Marker mode only.
#' @param min_score Minimum score below which a cluster is labeled
#'   \code{unassigned_label}. Default \code{NA} (disabled). \code{NULL}
#'   uses the mode-specific default.
#' @param min_margin Minimum margin between best and 2nd-best score / vote.
#'   Default \code{NA} (disabled). \code{NULL} uses the mode-specific default.
#' @param unassigned_label Label applied to clusters that fail the
#'   thresholds. Default \code{"Unknown"}.
#' @param return_scores One of \code{"none"} (default), \code{"cluster"},
#'   or \code{"cell"}. \code{"none"} returns the annotated Seurat object
#'   (legacy behavior). \code{"cluster"} returns a list
#'   \code{list(obj, scores)} where \code{scores} is the cluster × label
#'   score matrix (UCell mean scores in marker mode, SingleR vote
#'   fractions in singler mode). \code{"cell"} returns a list
#'   \code{list(obj, scores)} where \code{scores} is the cell × label
#'   matrix (UCell per-cell scores in marker mode, SingleR per-cell
#'   scores / labels in singler mode).
#' @return Either the annotated Seurat object (\code{return_scores =
#'   "none"}), or a list with elements \code{obj} and \code{scores}.
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
                             tumor_mode               = FALSE,
                             visium_mode              = FALSE,
                             filter_nonspecific       = TRUE,
                             nonspecific_max_fraction = 0.5,
                             high_relative_to_max     = 0.5,
                             min_detection_frac       = NA,
                             min_score                = NA,
                             min_margin               = NA,
                             unassigned_label         = "Unknown",
                             return_scores            = c("none", "cluster",
                                                          "cell")) {

  method        <- match.arg(method)
  return_scores <- match.arg(return_scores)

  if (!cluster_col %in% colnames(obj@meta.data)) {
    stop("Cluster column '", cluster_col, "' not found in metadata.")
  }
  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay

  # ---- Tumor mode: apply softer defaults (only if user didn't override) ---
  # Uses missing() to detect arguments that weren't explicitly passed,
  # so a user who supplied filter_nonspecific = TRUE still gets it on
  # even in tumor mode.
  if (isTRUE(tumor_mode)) {
    if (missing(filter_nonspecific)) {
      filter_nonspecific <- FALSE
    }
    if (missing(high_relative_to_max)) {
      high_relative_to_max <- 0.7
    }
    warning(
      "tumor_mode = TRUE: marker-based labels in tumors are approximate. ",
      "Tumor cells often lose canonical markers and gain ",
      "tumor-program signatures (proliferation, EMT, hypoxia) not in ",
      "normal-tissue marker lists. Consider identifying malignant cells ",
      "first (InferCNV / CopyKAT / SCEVAN), then annotating the ",
      "immune/stromal compartment with this function and scoring tumor ",
      "subtypes separately.",
      call. = FALSE
    )
  }

  # ---- Visium mode: apply the documented Visium fallback recipe ----------
  # Same missing()-gated approach as tumor_mode above -- only changes
  # defaults, and only filter_nonspecific actually has anything to change
  # (marker mode only; SingleR mode has no equivalent knob). The warning
  # applies regardless of method, since the underlying "spot is a mixture
  # of cells" caveat is the same either way.
  if (isTRUE(visium_mode)) {
    if (method == "marker" && missing(filter_nonspecific)) {
      filter_nonspecific <- FALSE
    }
    warning(
      "visium_mode = TRUE: this is still a winner-takes-all classifier. ",
      "A Visium spot is a mixture of 1-10 cells of possibly different ",
      "types, so per-spot/per-cluster labels lose minority populations ",
      "by construction -- RunRCTD() (spacexr) is the recommended primary ",
      "annotation strategy for Visium and estimates per-spot cell-type ",
      "proportions instead. If you're using this function anyway (e.g. ",
      "no suitable single-cell reference for RCTD, or as a quick ",
      "descriptive label alongside RCTD's proportions), consider bulking ",
      "up marker sets to 30-50+ genes and return_scores = 'cluster' to ",
      "inspect minority signal directly.",
      call. = FALSE
    )
  }

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

  # NB: .assign_with_unassigned() itself is a top-level internal function
  # (see bottom of file) rather than a closure here, so it's independently
  # unit-testable with a synthetic score matrix.

  if (method == "marker") {
    if (!requireNamespace("UCell", quietly = TRUE)) {
      stop("'UCell' is required for marker-based annotation.")
    }
    if (is.null(markers) || !is.list(markers) || is.null(names(markers))) {
      stop("`markers` must be a named list (cell_type -> character vector of genes).")
    }
    not_char <- names(markers)[!vapply(markers, is.character, logical(1))]
    if (length(not_char) > 0) {
      stop("`markers` element(s) must be character vectors of gene names; ",
           "not character: ", paste(not_char, collapse = ", "))
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
              ". Check gene symbol case/format and species (e.g. 'Cd3e' vs 'CD3E').",
              call. = FALSE)
      filtered_markers <- filtered_markers[setdiff(names(filtered_markers), empty)]
    }
    if (length(filtered_markers) == 0) {
      stop("No cell types have any marker genes present in assay '", a,
           "'. Cannot run marker-based annotation.")
    }
    markers <- filtered_markers

    # ---- Detection filter: drop signatures with too-sparse detection ------
    # For each signature, count cells where >=1 marker gene has count > 0.
    # Useful for Visium / sparse datasets where a signature might score on
    # pure noise if hardly any cell detects any of its markers, and for
    # catching wrong-tissue marker lists early.
    if (!is.na(min_detection_frac)) {
      if (min_detection_frac < 0 || min_detection_frac > 1) {
        stop("`min_detection_frac` must be between 0 and 1 (or NA).")
      }
      n_cells <- ncol(obj[[a]])
      sig_frac <- vapply(markers, function(g) {
        expr <- tryCatch(
          Seurat::FetchData(obj, vars = g, layer = "counts"),
          error = function(e) NULL
        )
        if (is.null(expr) || ncol(expr) == 0) return(0)
        sum(rowSums(expr > 0) > 0) / n_cells
      }, numeric(1))
      drop_sig <- names(sig_frac)[sig_frac < min_detection_frac]
      if (length(drop_sig)) {
        message(sprintf(
          "  Dropping %d signature(s) below min_detection_frac = %.2f:",
          length(drop_sig), min_detection_frac))
        for (ct in drop_sig) {
          message(sprintf("    '%s': detected in %.1f%% of cells",
                          ct, 100 * sig_frac[ct]))
        }
        markers <- markers[setdiff(names(markers), drop_sig)]
      }
      if (length(markers) == 0) {
        stop("All signatures fell below `min_detection_frac`. Lower the ",
             "threshold or check that your marker list matches the tissue.")
      }
    }

    # ---- Filter non-specific (broadly-expressed) marker genes -------------
    if (isTRUE(filter_nonspecific)) {
      all_marker_genes <- unique(unlist(markers))
      clusters_vec     <- as.character(obj@meta.data[[cluster_col]])

      # genes (rows) x cells (cols) sparse matrix, subset to just the marker
      # genes -- kept sparse via GetAssayData()/Matrix::rowMeans() rather
      # than FetchData(), which returns a dense cells x genes data.frame for
      # the whole marker gene set (can be sizable: all cells x all marker
      # genes, densified, just to average it right back down per cluster).
      expr_mat <- tryCatch({
        genes_present <- intersect(all_marker_genes, rownames(obj[[a]]))
        Seurat::GetAssayData(obj, assay = a, layer = "data")[genes_present, , drop = FALSE]
      }, error = function(e) NULL)

      if (is.null(expr_mat)) {
        message("  filter_nonspecific: FetchData failed; skipping filter.")
      } else {
        cluster_ids <- levels(factor(clusters_vec))
        # sapply over clusters -> genes (rows) x clusters (cols); transpose
        # to match the clusters (rows) x genes (cols) orientation downstream
        # code expects (same shape FetchData()+colMeans() produced before).
        cluster_means <- t(sapply(cluster_ids, function(cl) {
          idx <- which(clusters_vec == cl)
          Matrix::rowMeans(expr_mat[, idx, drop = FALSE])
        }))
        rownames(cluster_means) <- cluster_ids

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
    colnames(score_mat) <- names(markers)  # strip "_score" suffix
    # do.call(rbind, lapply(...)) rather than t(sapply(...)): when there's
    # only one surviving cell type (e.g. every other signature got dropped
    # for having no marker genes in the assay), score_mat has a single
    # column, so each per-cluster colMeans() result is a length-1 vector.
    # sapply() then simplifies the whole result to a plain vector (one
    # entry per cluster, not per label) before t() ever sees it, producing
    # a 1-row matrix instead of one row per cluster -- rbind() of the
    # per-cluster vectors is robust to this regardless of column count.
    # na.rm = TRUE so a few NA-scored cells (e.g. a signature fully
    # undetected in those cells) don't wipe out an otherwise-informative
    # mean for the rest of the cluster.
    avg_per_cluster <- do.call(rbind, lapply(
      split(seq_len(nrow(score_mat)), clusters),
      function(idx) colMeans(score_mat[idx, , drop = FALSE], na.rm = TRUE)
    ))
    rownames(avg_per_cluster) <- levels(factor(clusters))
    colnames(avg_per_cluster) <- names(markers)

    lookup <- .assign_with_unassigned(avg_per_cluster, names(markers),
                                      min_score, min_margin, unassigned_label)
    obj[[new_col]] <- unname(lookup[clusters])

    return(.maybe_return_scores(obj, avg_per_cluster, score_mat,
                                return_scores))
  }

  # SingleR mode -----------------------------------------------------------
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
  # do.call(rbind, lapply(...)) rather than t(sapply(...)) -- see the
  # matching comment in marker mode above; this has the identical failure
  # mode when the reference contributes only a single label.
  vote_frac <- do.call(rbind, lapply(split(seq_along(clusters), clusters), function(idx) {
    tab <- table(factor(per_cell_labels[idx], levels = uniq_labels))
    as.numeric(tab) / max(1, sum(tab))
  }))
  rownames(vote_frac) <- levels(factor(clusters))
  colnames(vote_frac) <- uniq_labels

  lookup <- .assign_with_unassigned(vote_frac, uniq_labels,
                                    min_score, min_margin, unassigned_label)
  obj[[new_col]] <- unname(lookup[clusters])

  # For SingleR "cell" scores: prefer the full scores matrix from pred$scores
  # if available; otherwise fall back to the per-cell labels as a one-column
  # data frame.
  cell_scores <- if (!is.null(pred$scores)) {
    as.matrix(pred$scores)
  } else {
    matrix(per_cell_labels, ncol = 1,
           dimnames = list(names(per_cell_labels), "label"))
  }

  .maybe_return_scores(obj, vote_frac, cell_scores, return_scores)
}


# ============================================================================
# Internal helper: pick the best-scoring label per row, applying the
# min_score/min_margin "Unknown" thresholds. Shared by marker mode
# (rows = clusters, avg UCell score per label) and SingleR mode
# (rows = clusters, majority-vote fraction per label).
# ============================================================================

#' @keywords internal
#' @noRd
.assign_with_unassigned <- function(score_mat, candidate_labels,
                                    min_score, min_margin, unassigned_label) {
  # NA scores are treated as -Inf (never best/second-best) rather than
  # left to propagate -- max()/sort() default to na.rm = FALSE, so a
  # single NA cell-type score can otherwise turn into an NA cluster mean,
  # then an NA best_score/margin for the whole cluster, silently bypassing
  # the min_score/min_margin thresholds (and can make apply(..., sort,
  # ...) fail to simplify into a matrix at all). Most plausible on sparse
  # data -- e.g. Visium -- where a signature can go fully undetected
  # across a subset of cells.
  all_na <- apply(is.na(score_mat), 1, all)
  sm <- score_mat
  sm[is.na(sm)] <- -Inf
  best_idx   <- apply(sm, 1, which.max)
  best_score <- apply(sm, 1, max)
  sorted_per_row <- apply(sm, 1, sort, decreasing = TRUE)
  second_score <- if (ncol(sm) >= 2) sorted_per_row[2, ] else -Inf
  margins <- best_score - second_score
  labs <- candidate_labels[best_idx]
  if (!is.na(min_score)) {
    labs[best_score < min_score] <- unassigned_label
  }
  if (!is.na(min_margin)) {
    labs[margins < min_margin] <- unassigned_label
  }
  # A cluster with literally no non-NA score for any label has zero
  # signal to assign from -- always Unknown, regardless of whether
  # min_score/min_margin are enabled (which default to disabled).
  labs[all_na] <- unassigned_label
  setNames(labs, rownames(score_mat))
}


# ============================================================================
# Internal helper: build the return value based on `return_scores`.
# ============================================================================

#' @keywords internal
#' @noRd
.maybe_return_scores <- function(obj,
                                 cluster_scores,
                                 cell_scores,
                                 return_scores) {
  switch(
    return_scores,
    none    = obj,
    cluster = list(obj = obj, scores = cluster_scores),
    cell    = list(obj = obj, scores = cell_scores)
  )
}
