#' Test for shifts in cell-type composition between conditions
#'
#' Wraps \code{speckle::propeller} — an empirical-Bayes-moderated ANOVA on
#' arcsin-square-root or logit-transformed proportions — for testing
#' whether cell-type composition differs between two or more conditions.
#' Propeller correctly handles the compositional nature of proportions
#' (each cell type competes with the others for probability mass) and
#' controls FDR across cell types.
#'
#' If \code{speckle} is unavailable, falls back to per-cluster
#' beta-regression (via \code{betareg}) as a simpler alternative, and if
#' that is also unavailable, to a per-cluster Wilcoxon test on the raw
#' proportions with BH FDR correction.
#'
#' @param obj A Seurat object.
#' @param cluster_col Metadata column holding cluster / cell-type ids.
#'   Default \code{"seurat_clusters"}. Ignored when \code{weight_cols} is
#'   supplied.
#' @param sample_col Metadata column identifying biological replicates
#'   (typically \code{"orig.ident"}).
#' @param condition_col Metadata column holding the condition to compare.
#' @param weight_cols Optional character vector of numeric metadata columns
#'   holding a continuous per-cell/per-spot composition -- e.g. the
#'   \code{rctd_<celltype>} proportion columns \code{\link{RunRCTD}} writes
#'   in \code{"full"} mode. When supplied, \code{cluster_col} is ignored and
#'   the test runs on the mean of each weight column within each sample,
#'   rather than on proportions derived from discrete cluster counts. This
#'   uses the full continuous mixture RCTD estimates instead of collapsing
#'   each spot to a single dominant cell type first, which throws away
#'   information "full" mode specifically exists to capture. Only
#'   \code{method = "betareg"} or \code{"wilcox"} support this --
#'   \code{"propeller"} needs a discrete per-cell cluster factor, not
#'   continuous weights, and errors if combined with \code{weight_cols}.
#' @param transform Transformation applied by propeller before ANOVA.
#'   One of \code{"asin"} (arcsin-sqrt, default; propeller's recommendation)
#'   or \code{"logit"}. Ignored when \code{weight_cols} is supplied.
#' @param method Backend: \code{"auto"} (default; try propeller, then
#'   betareg, then wilcox -- or, with \code{weight_cols}, betareg then
#'   wilcox), \code{"propeller"}, \code{"betareg"}, or \code{"wilcox"}.
#' @return A data frame with one row per cluster (or, with \code{weight_cols},
#'   one row per weight column) and columns \code{cluster},
#'   \code{mean_prop_<level>} for each condition level, \code{effect},
#'   \code{stat}, \code{pvalue}, \code{padj}, and \code{method}. When the
#'   \code{propeller} backend runs, an additional \code{stat_type} column
#'   ("T" or "F") records which test \code{speckle::propeller()} actually
#'   used: a t-test (\code{effect} = \code{PropRatio}) when
#'   \code{condition_col} has exactly 2 levels, or an ANOVA F-test (no
#'   single ratio, so \code{effect} is \code{NA}) when it has more than 2.
#' @examples
#' \dontrun{
#' res <- CompositionalTest(obj,
#'                          sample_col    = "orig.ident",
#'                          condition_col = "treatment")
#' subset(res, padj < 0.05)
#'
#' # Continuous RCTD "full"-mode weights instead of a discrete cluster call
#' celltypes <- grep("^rctd_", colnames(visium@meta.data), value = TRUE)
#' celltypes <- setdiff(celltypes, "rctd_dominant")
#' res_rctd <- CompositionalTest(visium,
#'                               weight_cols   = celltypes,
#'                               sample_col    = "orig.ident",
#'                               condition_col = "treatment",
#'                               method        = "betareg")
#' }
#' @export
CompositionalTest <- function(obj,
                              cluster_col   = "seurat_clusters",
                              sample_col    = "orig.ident",
                              condition_col = NULL,
                              weight_cols   = NULL,
                              transform     = c("asin", "logit"),
                              method        = c("auto", "propeller",
                                                "betareg", "wilcox")) {

  transform <- match.arg(transform)
  method    <- match.arg(method)
  if (is.null(condition_col)) stop("`condition_col` is required.")

  # ---- Continuous-weights mode (e.g. RCTD "full"-mode rctd_<celltype>) ----
  if (!is.null(weight_cols)) {
    if (method == "propeller") {
      stop("`method = 'propeller'` requires a discrete per-cell cluster ",
           "factor (speckle::propeller() doesn't accept continuous weights) ",
           "-- use `method = 'betareg'` or `'wilcox'` with `weight_cols`.")
    }
    missing_cols <- setdiff(c(weight_cols, sample_col, condition_col),
                            colnames(obj@meta.data))
    if (length(missing_cols)) {
      stop("Column(s) not found in obj@meta.data: ",
           paste(missing_cols, collapse = ", "))
    }
    md <- obj@meta.data
    not_numeric <- weight_cols[!vapply(md[weight_cols], is.numeric, logical(1))]
    if (length(not_numeric)) {
      stop("`weight_cols` must be numeric columns; not numeric: ",
           paste(not_numeric, collapse = ", "))
    }

    sample_vec <- as.character(md[[sample_col]])
    cond_vec   <- as.character(md[[condition_col]])
    sc <- unique(data.frame(sample = sample_vec, cond = cond_vec,
                            stringsAsFactors = FALSE))
    dupes <- sc$sample[duplicated(sc$sample)]
    if (length(dupes)) {
      stop("Samples have inconsistent condition values: ",
           paste(utils::head(dupes, 5), collapse = ", "))
    }

    n_per_sample <- table(sample_vec)
    eps <- 1 / (2 * max(n_per_sample))

    prop_df <- do.call(rbind, lapply(weight_cols, function(wc) {
      means <- tapply(md[[wc]], sample_vec, mean, na.rm = TRUE)
      data.frame(sample = names(means), cluster = wc,
                prop = as.numeric(means), stringsAsFactors = FALSE)
    }))
    prop_df$prop_bounded <- pmin(pmax(prop_df$prop, eps), 1 - eps)
    prop_df$cond <- unname(setNames(sc$cond, sc$sample)[prop_df$sample])

    backend <- method
    if (backend == "auto") {
      backend <- if (requireNamespace("betareg", quietly = TRUE)) "betareg" else "wilcox"
    }

    res <- switch(
      backend,
      betareg = .comp_test_betareg_continuous(prop_df),
      wilcox  = .comp_test_wilcox_continuous(prop_df)
    )
    res$method <- backend

    mp <- stats::aggregate(prop ~ cluster + cond, data = prop_df, FUN = mean)
    for (lv in unique(prop_df$cond)) {
      col <- paste0("mean_prop_", lv)
      res[[col]] <- mp$prop[match(paste(res$cluster, lv),
                                  paste(mp$cluster, mp$cond))]
    }
    res <- res[order(res$padj, na.last = TRUE), ]
    rownames(res) <- NULL
    return(res)
  }

  # ---- Discrete cluster_col mode (original behavior) ----------------------
  for (col in c(cluster_col, sample_col, condition_col)) {
    if (!col %in% colnames(obj@meta.data)) {
      stop("Column '", col, "' not found in obj@meta.data.")
    }
  }

  md <- obj@meta.data
  md_ss <- data.frame(
    cell    = rownames(md),
    cluster = as.character(md[[cluster_col]]),
    sample  = as.character(md[[sample_col]]),
    cond    = as.character(md[[condition_col]]),
    stringsAsFactors = FALSE
  )
  # Sanity: condition constant within sample
  sc <- unique(md_ss[, c("sample", "cond")])
  dupes <- sc$sample[duplicated(sc$sample)]
  if (length(dupes)) {
    stop("Samples have inconsistent condition values: ",
         paste(head(dupes, 5), collapse = ", "))
  }

  # Choose backend
  backend <- method
  if (backend == "auto") {
    backend <- if (requireNamespace("speckle", quietly = TRUE))     "propeller"
               else if (requireNamespace("betareg", quietly = TRUE)) "betareg"
               else                                                   "wilcox"
  }

  # Dispatch
  res <- switch(
    backend,
    propeller = .comp_test_propeller(md_ss, transform),
    betareg   = .comp_test_betareg(md_ss),
    wilcox    = .comp_test_wilcox(md_ss)
  )
  res$method <- backend

  # Add per-condition mean proportions for interpretability
  prop_df <- as.data.frame(table(sample = md_ss$sample,
                                 cluster = md_ss$cluster),
                           responseName = "n_cells",
                           stringsAsFactors = FALSE)
  tot <- stats::aggregate(n_cells ~ sample, data = prop_df, FUN = sum)
  names(tot)[2] <- "n_total"
  prop_df <- merge(prop_df, tot, by = "sample")
  prop_df$prop <- prop_df$n_cells / pmax(1, prop_df$n_total)
  prop_df$cond <- unname(setNames(sc$cond, sc$sample)[prop_df$sample])
  mean_props <- stats::aggregate(prop ~ cluster + cond,
                                 data = prop_df, FUN = mean)
  for (lv in unique(mean_props$cond)) {
    col <- paste0("mean_prop_", lv)
    res[[col]] <- mean_props$prop[
      match(paste(res$cluster, lv),
            paste(mean_props$cluster, mean_props$cond))
    ]
  }
  res <- res[order(res$padj, na.last = TRUE), ]
  rownames(res) <- NULL
  res
}


# ============================================================================
# Backend: propeller
# ============================================================================
#' @keywords internal
#' @noRd
.comp_test_propeller <- function(md_ss, transform) {
  if (!requireNamespace("speckle", quietly = TRUE)) {
    stop("'speckle' is required for propeller.")
  }
  res <- speckle::propeller(
    clusters  = md_ss$cluster,
    sample    = md_ss$sample,
    group     = md_ss$cond,
    transform = transform
  )
  # speckle::propeller() dispatches internally on how many levels `group`
  # has: exactly 2 -> propeller.ttest() (PropRatio + Tstatistic columns);
  # more than 2 -> propeller.anova() instead (no PropRatio -- there's no
  # single ratio across >2 groups -- and Fstatistic in place of Tstatistic).
  # Assuming Tstatistic always exists crashes data.frame() below in the
  # >2-group case (res$Tstatistic is NULL, not NA, so its 0-length can't be
  # recycled against the other columns) -- check for either column instead.
  stat_col <- if ("Tstatistic" %in% colnames(res)) {
    "Tstatistic"
  } else if ("Fstatistic" %in% colnames(res)) {
    "Fstatistic"
  } else {
    NA
  }
  data.frame(
    cluster   = rownames(res),
    effect    = if ("PropRatio" %in% colnames(res)) res$PropRatio else NA,
    stat      = if (!is.na(stat_col)) res[[stat_col]] else NA,
    stat_type = if (is.na(stat_col)) NA_character_ else sub("statistic$", "", stat_col),
    pvalue    = res$P.Value,
    padj      = res$FDR,
    stringsAsFactors = FALSE
  )
}


# ============================================================================
# Backend: per-cluster beta regression on the raw proportions
# ============================================================================
#' @keywords internal
#' @noRd
.comp_test_betareg <- function(md_ss) {
  if (!requireNamespace("betareg", quietly = TRUE)) {
    stop("'betareg' is required for betareg backend.")
  }
  prop_df <- as.data.frame(table(sample = md_ss$sample,
                                 cluster = md_ss$cluster),
                           responseName = "n", stringsAsFactors = FALSE)
  tot <- stats::aggregate(n ~ sample, data = prop_df, FUN = sum)
  names(tot)[2] <- "n_total"
  prop_df <- merge(prop_df, tot, by = "sample")
  prop_df$prop <- prop_df$n / pmax(1, prop_df$n_total)
  # betareg needs strictly (0, 1)
  eps <- 1 / (2 * max(prop_df$n_total))
  prop_df$prop_bounded <- pmin(pmax(prop_df$prop, eps), 1 - eps)
  cond_map <- unique(md_ss[, c("sample", "cond")])
  prop_df$cond <- unname(setNames(cond_map$cond,
                                  cond_map$sample)[prop_df$sample])

  clus <- unique(prop_df$cluster)
  res <- lapply(clus, function(cl) {
    d <- prop_df[prop_df$cluster == cl, ]
    if (length(unique(d$cond)) < 2 || nrow(d) < 4) {
      return(data.frame(cluster = cl, effect = NA, stat = NA,
                        pvalue = NA, padj = NA))
    }
    fit <- tryCatch(
      betareg::betareg(prop_bounded ~ cond, data = d),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(data.frame(cluster = cl, effect = NA, stat = NA,
                        pvalue = NA, padj = NA))
    }
    smry <- summary(fit)$coefficients$mean
    # Row 2 is the condition effect
    data.frame(
      cluster = cl,
      effect  = smry[2, "Estimate"],
      stat    = smry[2, "z value"],
      pvalue  = smry[2, "Pr(>|z|)"],
      padj    = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, res)
  out$padj <- stats::p.adjust(out$pvalue, method = "BH")
  out
}


# ============================================================================
# Backend: per-cluster Wilcoxon on raw proportions
# ============================================================================
#' @keywords internal
#' @noRd
.comp_test_wilcox <- function(md_ss) {
  prop_df <- as.data.frame(table(sample = md_ss$sample,
                                 cluster = md_ss$cluster),
                           responseName = "n", stringsAsFactors = FALSE)
  tot <- stats::aggregate(n ~ sample, data = prop_df, FUN = sum)
  names(tot)[2] <- "n_total"
  prop_df <- merge(prop_df, tot, by = "sample")
  prop_df$prop <- prop_df$n / pmax(1, prop_df$n_total)
  cond_map <- unique(md_ss[, c("sample", "cond")])
  prop_df$cond <- unname(setNames(cond_map$cond,
                                  cond_map$sample)[prop_df$sample])

  clus <- unique(prop_df$cluster)
  res <- lapply(clus, function(cl) {
    d <- prop_df[prop_df$cluster == cl, ]
    lv <- sort(unique(d$cond))
    if (length(lv) != 2) {
      return(data.frame(cluster = cl, effect = NA, stat = NA,
                        pvalue = NA, padj = NA))
    }
    a <- d$prop[d$cond == lv[1]]
    b <- d$prop[d$cond == lv[2]]
    if (length(a) < 2 || length(b) < 2) {
      return(data.frame(cluster = cl, effect = NA, stat = NA,
                        pvalue = NA, padj = NA))
    }
    w <- suppressWarnings(stats::wilcox.test(b, a))
    data.frame(
      cluster = cl,
      effect  = mean(b) - mean(a),
      stat    = unname(w$statistic),
      pvalue  = w$p.value,
      padj    = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, res)
  out$padj <- stats::p.adjust(out$pvalue, method = "BH")
  out
}


# ============================================================================
# Backend: per-weight-column beta regression on continuous means (weight_cols
# mode -- e.g. RCTD "full"-mode rctd_<celltype> proportions). Mirrors
# .comp_test_betareg() above but takes an already-built prop_df (sample,
# cluster, prop, prop_bounded, cond) rather than building one from discrete
# cluster counts via table().
# ============================================================================
#' @keywords internal
#' @noRd
.comp_test_betareg_continuous <- function(prop_df) {
  if (!requireNamespace("betareg", quietly = TRUE)) {
    stop("'betareg' is required for betareg backend.")
  }
  clus <- unique(prop_df$cluster)
  res <- lapply(clus, function(cl) {
    d <- prop_df[prop_df$cluster == cl, ]
    if (length(unique(d$cond)) < 2 || nrow(d) < 4) {
      return(data.frame(cluster = cl, effect = NA, stat = NA,
                        pvalue = NA, padj = NA))
    }
    fit <- tryCatch(
      betareg::betareg(prop_bounded ~ cond, data = d),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(data.frame(cluster = cl, effect = NA, stat = NA,
                        pvalue = NA, padj = NA))
    }
    smry <- summary(fit)$coefficients$mean
    data.frame(
      cluster = cl,
      effect  = smry[2, "Estimate"],
      stat    = smry[2, "z value"],
      pvalue  = smry[2, "Pr(>|z|)"],
      padj    = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, res)
  out$padj <- stats::p.adjust(out$pvalue, method = "BH")
  out
}


# ============================================================================
# Backend: per-weight-column Wilcoxon on continuous means (weight_cols mode).
# Mirrors .comp_test_wilcox() above but takes an already-built prop_df.
# ============================================================================
#' @keywords internal
#' @noRd
.comp_test_wilcox_continuous <- function(prop_df) {
  clus <- unique(prop_df$cluster)
  res <- lapply(clus, function(cl) {
    d <- prop_df[prop_df$cluster == cl, ]
    lv <- sort(unique(d$cond))
    if (length(lv) != 2) {
      return(data.frame(cluster = cl, effect = NA, stat = NA,
                        pvalue = NA, padj = NA))
    }
    a <- d$prop[d$cond == lv[1]]
    b <- d$prop[d$cond == lv[2]]
    if (length(a) < 2 || length(b) < 2) {
      return(data.frame(cluster = cl, effect = NA, stat = NA,
                        pvalue = NA, padj = NA))
    }
    w <- suppressWarnings(stats::wilcox.test(b, a))
    data.frame(
      cluster = cl,
      effect  = mean(b) - mean(a),
      stat    = unname(w$statistic),
      pvalue  = w$p.value,
      padj    = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, res)
  out$padj <- stats::p.adjust(out$pvalue, method = "BH")
  out
}
