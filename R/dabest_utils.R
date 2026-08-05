# Internal helpers shared by dabestr-backed "estimation plot" wrappers.
# Not exported. Used by both CompositionEstimationPlot() (proportions by
# cell type/cluster) and NicheCoExpressEstimationPlot() (co-expression
# scores by niche x gene pair) -- both wrap the same long-form shape, one
# row per (sample, group, condition), so the dabestr-facing logic lives
# here once rather than being duplicated in each.

#' @keywords internal
#' @noRd
.dabest_from_long <- function(df,
                              group_col,
                              y_col,
                              condition_col,
                              idx          = NULL,
                              effect       = c("mean_diff", "median_diff", "cohens_d",
                                               "hedges_g", "cliffs_delta", "cohens_h"),
                              group_levels = NULL) {
  effect <- match.arg(effect)

  if (!requireNamespace("dabestr", quietly = TRUE)) {
    stop("'dabestr' is required. Install with: install.packages('dabestr')")
  }
  if (!requireNamespace("rlang", quietly = TRUE)) {
    stop("'rlang' is required (should already be installed as a ",
         "SingleCellTools dependency).")
  }
  for (col in c(group_col, y_col, condition_col)) {
    if (!col %in% colnames(df)) stop("Column '", col, "' not found in `df`.")
  }

  cond_levels <- unique(as.character(df[[condition_col]]))

  if (is.null(idx)) {
    if (length(cond_levels) != 2) {
      stop("`idx` must be supplied explicitly when `condition_col` has ",
           length(cond_levels), " level(s) (found: ",
           paste(cond_levels, collapse = ", "),
           "). Pass idx = c(<reference>, <test>) to control which two ",
           "(and in what order) are compared.")
    }
    idx <- sort(cond_levels)
    message("`idx` not supplied -- comparing '", idx[1], "' (reference) vs '",
            idx[2], "' (test), alphabetical order. Pass idx = c(...) to ",
            "control this.")
  } else {
    if (length(idx) != 2) stop("`idx` must have exactly 2 elements.")
    missing_idx <- setdiff(idx, cond_levels)
    if (length(missing_idx)) {
      stop("`idx` value(s) not found in `df[[condition_col]]`: ",
           paste(missing_idx, collapse = ", "))
    }
  }

  all_groups <- unique(as.character(df[[group_col]]))
  if (is.null(group_levels)) group_levels <- all_groups
  missing_groups <- setdiff(group_levels, all_groups)
  if (length(missing_groups)) {
    stop("`group_levels` not found in `df[[group_col]]`: ",
         paste(missing_groups, collapse = ", "))
  }

  effect_fn <- switch(effect,
    mean_diff    = dabestr::mean_diff,
    median_diff  = dabestr::median_diff,
    cohens_d     = dabestr::cohens_d,
    hedges_g     = dabestr::hedges_g,
    cliffs_delta = dabestr::cliffs_delta,
    cohens_h     = dabestr::cohens_h
  )

  out <- vector("list", length(group_levels))
  names(out) <- group_levels

  for (g in group_levels) {
    sub <- df[as.character(df[[group_col]]) == g, , drop = FALSE]
    sub <- sub[as.character(sub[[condition_col]]) %in% idx, , drop = FALSE]
    sub[[condition_col]] <- factor(as.character(sub[[condition_col]]), levels = idx)
    sub[[y_col]] <- as.numeric(sub[[y_col]])

    # Both idx levels must have >= 1 row for THIS group. Checking nrow(sub)
    # > 0 isn't enough: a group can have plenty of rows but all from one
    # condition (e.g. a niche/cluster only present in "anterior" samples,
    # none in "posterior"). dabestr::load() itself would fail on that with
    # an opaque "<level> not present in x" error deep inside dabestr; catch
    # it here instead with a message that names the group and the missing
    # level(s), and skip this group rather than aborting the whole batch --
    # relevant when group_levels spans many groups (e.g. group_levels =
    # NULL / niches = NULL) and only some of them actually have both
    # conditions represented.
    present <- idx %in% as.character(sub[[condition_col]])
    if (!all(present)) {
      warning("Skipping group '", g, "': idx level(s) ",
              paste(idx[!present], collapse = ", "),
              " have 0 rows for this group after filtering `condition_col`. ",
              "This group isn't present in both conditions in the ",
              "underlying data -- not a bug, just nothing to compare here.")
      next
    }

    # dabestr::load() takes `x`/`y` as bare (tidy-eval) column references,
    # not strings -- rlang::inject() + sym() builds the call with the
    # actual column-name symbols substituted in before dispatch, so this
    # works regardless of how load() captures its own arguments internally.
    loaded <- rlang::inject(dabestr::load(
      data = sub,
      x    = !!rlang::sym(condition_col),
      y    = !!rlang::sym(y_col),
      idx  = idx
    ))

    out[[g]] <- effect_fn(loaded)
  }

  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) {
    stop("No requested group had both idx levels (", paste(idx, collapse = ", "),
         ") represented -- nothing to plot. Requested group(s): ",
         paste(group_levels, collapse = ", "), ".")
  }

  out
}
