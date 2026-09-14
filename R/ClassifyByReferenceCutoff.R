#' Classify cells against a reference-derived quantile cutoff
#'
#' Generalizes a "score each cell against marker panels, then classify by
#' comparing to a fixed reference-derived cutoff" pattern originally written
#' by hand for liver zonation calls (pericentral / periportal / midlobular
#' from marker UCell scores), but the logic itself has nothing
#' zonation-specific about it: any small number of marker panels scored with
#' \code{UCell::AddModuleScore_UCell}, a reference/baseline population to
#' anchor the cutoff to, and an argmax-vs-cutoff decision rule with an
#' explicit "none of them" fallback, fits the same shape.
#'
#' \strong{The rule.} For each cell restricted to \code{annotation_col \%in\%
#' target} (e.g. \code{Final_Annotation == "Hepatocytes"}), compute each
#' marker panel's UCell score relative to that panel's reference cutoff. If
#' every panel's score is at or below its cutoff, the cell is
#' \code{"Unclassified"}; otherwise it's assigned the panel with the largest
#' margin above cutoff. Cells outside \code{target} keep their original
#' \code{annotation_col} label in the output column -- exactly the "zone
#' calls for hepatocytes, cell-type labels for everything else" pattern the
#' zonation scripts this generalizes always wanted.
#'
#' \strong{The cutoff itself} is the \code{cutoff_pct} quantile of each
#' marker panel's score, computed once over a reference/baseline subset of
#' \code{target} cells (\code{ref_subset}, e.g. \code{Age.group == "4 mo"}
#' or \code{Condition \%in\% baseline_conditions}) and then applied to
#' \emph{every} cell in \code{target} -- baseline and non-baseline alike.
#' Pass \code{group_by} (e.g. \code{"Sex"}, \code{"Genotype"}) to compute and
#' apply a separate cutoff per group instead of one pooled cutoff; leave it
#' \code{NULL} (default) when the cutoff shouldn't vary by any grouping
#' variable -- including when that variable doesn't exist in this dataset at
#' all (no need to synthesize a constant column just to satisfy this
#' function).
#'
#' @param obj A Seurat object.
#' @param markers Named list of marker gene vectors, one per panel, e.g.
#'   \code{list(pericentral = c("Glul", "Cyp2e1"), periportal = c("Ass1",
#'   "Pck1"))}. Genes not present in \code{obj} are dropped per-panel with a
#'   message (\code{verbose = TRUE}); a panel left with zero genes errors.
#' @param annotation_col Metadata column giving each cell's current
#'   cell-type / region label (e.g. \code{"Final_Annotation"}).
#' @param target Value(s) of \code{annotation_col} to classify (e.g.
#'   \code{"Hepatocytes"}). Cells with any other \code{annotation_col} value
#'   are left alone and keep their existing label in the output.
#' @param ref_subset Optional logical expression (NSE, not a string -- same
#'   convention as \code{\link{SubsetSpatial}}'s \code{subset}), evaluated
#'   against the metadata of \code{target} cells only, that narrows them
#'   down to the reference/baseline population the cutoff is computed from
#'   (e.g. \code{Age.group == "4 mo"}). \code{NULL} (default) uses every
#'   \code{target} cell as the reference.
#' @param group_by Optional metadata column to compute a separate cutoff per
#'   level of (e.g. \code{"Sex"}). \code{NULL} (default) pools all reference
#'   cells into a single cutoff. Every group value seen among \code{target}
#'   cells being classified must also appear among the reference cells, or
#'   this errors naming the missing value(s).
#' @param cutoff_pct Quantile used as the cutoff, computed on the reference
#'   set. Default \code{0.75}.
#' @param assay Assay UCell scores are computed from when
#'   \code{compute_scores = TRUE}. Default \code{DefaultAssay(obj)}.
#' @param score_suffix Suffix UCell appends to each panel name to form its
#'   score column (passed to \code{AddModuleScore_UCell(..., name =
#'   score_suffix)}); also the suffix stripped off when reporting the
#'   winning panel. Default \code{"_UCell"} (UCell's own default).
#' @param out_col Metadata column to write classifications to. Default
#'   \code{"Zone_final"}.
#' @param compute_scores If \code{TRUE} (default), runs
#'   \code{UCell::AddModuleScore_UCell} on \code{markers} first. Set
#'   \code{FALSE} to reuse score columns already present in \code{obj}
#'   (\code{paste0(names(markers), score_suffix)}) -- e.g. if you already
#'   scored with different UCell arguments than this function would use.
#' @param return_cutoffs If \code{TRUE}, returns \code{list(object, cutoffs)}
#'   instead of just the classified object -- useful for sanity-checking the
#'   cutoff values before trusting them. Default \code{FALSE}.
#' @param verbose Message dropped marker genes and print the reference
#'   cutoffs. Default \code{TRUE}.
#' @return \code{obj} with \code{out_col} added (or \code{list(object,
#'   cutoffs)} if \code{return_cutoffs = TRUE}).
#' @examples
#' \dontrun{
#' markers <- list(pericentral = c("Glul", "Cyp2e1", "Oat"),
#'                periportal  = c("Ass1", "Pck1", "Sds"),
#'                midlobular  = c("Hamp2", "Igfbp2", "Ccnd1"))
#'
#' # Single-sex dataset -- no group_by needed
#' cosmx <- ClassifyByReferenceCutoff(
#'   cosmx, markers,
#'   annotation_col = "Final_Annotation", target = "Hepatocytes",
#'   ref_subset = Age.group == "4 mo", cutoff_pct = 0.75
#' )
#'
#' # Sex-specific baseline cutoffs
#' visium <- ClassifyByReferenceCutoff(
#'   visium, markers,
#'   annotation_col = "Annotation", target = "Hepatocytes",
#'   ref_subset = Condition %in% c("Female 4mo WT", "Male 4mo WT"),
#'   group_by = "Sex", cutoff_pct = 0.6
#' )
#' }
#' @export
ClassifyByReferenceCutoff <- function(obj,
                                      markers,
                                      annotation_col,
                                      target,
                                      ref_subset     = NULL,
                                      group_by       = NULL,
                                      cutoff_pct     = 0.75,
                                      assay          = NULL,
                                      score_suffix   = "_UCell",
                                      out_col        = "Zone_final",
                                      compute_scores = TRUE,
                                      return_cutoffs = FALSE,
                                      verbose        = TRUE) {

  .assert_seurat(obj)
  if (!is.list(markers) || is.null(names(markers)) || any(names(markers) == "")) {
    stop("`markers` must be a named list of character vectors, e.g. ",
         "list(pericentral = c(...), periportal = c(...)).")
  }
  if (!annotation_col %in% colnames(obj@meta.data)) {
    stop("`annotation_col` ('", annotation_col, "') not found in obj@meta.data.")
  }
  if (!is.null(group_by) && !group_by %in% colnames(obj@meta.data)) {
    stop("`group_by` ('", group_by, "') not found in obj@meta.data.")
  }
  # Mirrors SubsetSpatial()'s own NSE convention/comment: enquo() captures a
  # literal NULL default the same way it captures a caller-supplied
  # expression, so quo_is_null() (not is.null(ref_subset)) is what actually
  # distinguishes "no ref_subset given" from a real expression.
  ref_subset_q <- rlang::enquo(ref_subset)

  ## ---- 1. Marker panels, filtered to genes present; UCell scores --------
  markers_filtered <- lapply(markers, function(g) intersect(g, rownames(obj)))
  empty <- names(markers_filtered)[lengths(markers_filtered) == 0]
  if (length(empty) > 0) {
    stop("Marker set(s) have no genes present in `obj`: ",
         paste(empty, collapse = ", "))
  }
  if (isTRUE(verbose)) {
    dropped <- mapply(function(g, kept) setdiff(g, kept), markers, markers_filtered,
                      SIMPLIFY = FALSE)
    dropped <- dropped[lengths(dropped) > 0]
    if (length(dropped) > 0) {
      message("Marker genes dropped (not found in obj):")
      for (nm in names(dropped)) {
        message("  ", nm, ": ", paste(dropped[[nm]], collapse = ", "))
      }
    }
  }

  zone_cols <- paste0(names(markers_filtered), score_suffix)

  if (isTRUE(compute_scores)) {
    a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay
    obj <- UCell::AddModuleScore_UCell(obj, features = markers_filtered, assay = a,
                                       name = score_suffix)
    if (!all(zone_cols %in% colnames(obj@meta.data))) {
      stop("UCell did not produce the expected score columns; check marker ",
           "names / `score_suffix`.")
    }
  } else if (!all(zone_cols %in% colnames(obj@meta.data))) {
    stop("compute_scores = FALSE but expected score column(s) missing: ",
         paste(setdiff(zone_cols, colnames(obj@meta.data)), collapse = ", "),
         ". Compute them first (e.g. UCell::AddModuleScore_UCell) or set ",
         "compute_scores = TRUE.")
  }

  ## ---- 2. Reference/baseline set: target cell type, optionally filtered -
  meta <- obj@meta.data
  target_idx <- which(as.character(meta[[annotation_col]]) %in% target)
  if (length(target_idx) == 0) {
    stop("No cells found with `", annotation_col, "` %in% c(",
         paste(shQuote(target), collapse = ", "), ").")
  }

  if (!rlang::quo_is_null(ref_subset_q)) {
    keep <- rlang::eval_tidy(ref_subset_q, data = meta[target_idx, , drop = FALSE])
    if (!is.logical(keep)) {
      stop("`ref_subset` must evaluate to a logical vector over the target ",
           "cells' metadata.")
    }
    ref_idx <- target_idx[which(keep)]
  } else {
    ref_idx <- target_idx
  }
  if (length(ref_idx) == 0) {
    stop("No cells left for the reference/baseline set after `ref_subset` -- ",
         "check that its columns/values exist in obj@meta.data.")
  }

  ref <- meta[ref_idx, zone_cols, drop = FALSE]

  ## ---- 3. Cutoffs: quantile(cutoff_pct) of the reference set, optionally
  ##         split by `group_by` (e.g. Sex, Genotype) ----------------------
  if (!is.null(group_by)) {
    ref_groups <- as.character(meta[[group_by]][ref_idx])
    cutoffs <- sapply(split(ref, ref_groups),
                      function(df) apply(df, 2, stats::quantile, probs = cutoff_pct))
  } else {
    cutoffs <- apply(ref, 2, stats::quantile, probs = cutoff_pct)
  }
  if (isTRUE(verbose)) {
    message(sprintf(
      "Reference cutoffs (%.3g quantile%s, n = %d cell%s):",
      cutoff_pct,
      if (!is.null(group_by)) paste0(", by ", group_by) else "",
      length(ref_idx), if (length(ref_idx) == 1) "" else "s"))
    print(cutoffs)
  }

  ## ---- 4. Classify every target-annotation cell: argmax vs. its cutoff --
  scores <- as.matrix(meta[target_idx, zone_cols, drop = FALSE])

  if (!is.null(group_by)) {
    grp_vals <- as.character(meta[[group_by]][target_idx])
    missing_groups <- setdiff(unique(grp_vals), colnames(cutoffs))
    if (length(missing_groups) > 0) {
      stop("`", group_by, "` value(s) among target cells not present in the ",
           "reference cutoffs -- check `", group_by, "` levels: ",
           paste(missing_groups, collapse = ", "))
    }
    rel <- scores - t(cutoffs[, grp_vals, drop = FALSE])
  } else {
    # sweep() (matched by name, not position) rather than a bare
    # `scores - cutoffs`: matrix-minus-vector recycles column-major (down
    # each column), not per-column, so a plain subtraction would silently
    # misalign whenever nrow(scores) != length(cutoffs).
    rel <- sweep(scores, 2, cutoffs[colnames(scores)], "-")
  }

  zone_call <- apply(rel, 1, function(x) {
    if (all(x <= 0)) return("Unclassified")
    sub(paste0(score_suffix, "$"), "", names(which.max(x)))
  })

  out <- as.character(meta[[annotation_col]])
  out[target_idx] <- zone_call
  obj[[out_col]] <- out

  if (isTRUE(return_cutoffs)) {
    return(list(object = obj, cutoffs = cutoffs))
  }
  obj
}
