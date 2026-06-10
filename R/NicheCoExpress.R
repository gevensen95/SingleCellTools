#' Differential Gene Co-expression Across Spatial Niches
#'
#' Computes per-sample gene-pair co-expression scores within each spatial
#' niche, then tests for differences between two conditions. Adapted in
#' spirit from \code{katlande/scCoExpress}.
#'
#' Co-expression is measured by the Manders Overlap Coefficient (MOC)
#' between a pair of genes, computed across cells in one
#' (sample x niche) subset. MOC is the cosine similarity of the two
#' genes' expression vectors:
#' \deqn{MOC(A, B) = \sum_k A_k B_k / \sqrt{\sum_k A_k^2 \sum_k B_k^2}}
#' The observed MOC is expressed as a ratio against an abundance-matched
#' background of random gene pairs so that single-cell / spatial sparsity
#' doesn't drive the signal. Each gene is L2-normalized once per subset,
#' which makes background sampling fast (one dot product per pair).
#'
#' Differences from scCoExpress:
#' \enumerate{
#'   \item Work is partitioned by an upstream spatial-niche label.
#'   \item Co-expression is computed only among a preset gene list,
#'     within each niche.
#'   \item It is computed per sample, so each niche has one score per
#'     sample.
#'   \item Per niche and per gene pair, the per-sample scores are
#'     compared between two conditions with a statistical test and
#'     multiple-testing correction.
#' }
#'
#' \strong{Cell-type composition controls} (optional, via \code{celltype_col}):
#' \itemize{
#'   \item \code{celltypes}: restrict the analysis to listed cell types.
#'   \item \code{max_celltype_frac}: cap any one cell type at this
#'     fraction by random downsampling.
#'   \item \code{center_celltype}: subtract per-cell-type means before
#'     computing MOC. Converts the metric from "log2 ratio vs background"
#'     to "z-score vs background" and removes between-type composition
#'     effects fully. The strongest available control.
#' }
#' When \code{celltype_col} is set, the output also gains a per-niche
#' \code{comp_diff} column (absolute difference in dominant-cell-type
#' fraction between conditions) and a \code{comp_flag} boolean, so
#' compositional confounds in the differential test are visible.
#'
#' @param seurat_obj A pre-processed Seurat object; cells must already be
#'   assigned to niches.
#' @param genes Either a character vector of genes (all pairwise
#'   combinations are scored) OR a 2-column data.frame / matrix of
#'   specific gene pairs (columns \code{geneA}, \code{geneB}).
#' @param niche_col meta.data column with the spatial niche label.
#' @param sample_col meta.data column identifying biological samples.
#' @param condition_col meta.data column with the condition / group.
#' @param conditions Length-2 character vector: \code{c(reference, test)}.
#'   Defaults to the two levels found.
#' @param niches Subset of niches to analyse. Default: all present.
#' @param assay,layer Assay and layer to read. Default \code{RNA} /
#'   \code{data}. If the assay has the requested layer split into
#'   per-sample pieces (Seurat v5, e.g. \code{data.1}, \code{data.2}, ...),
#'   they are joined first so the expression matrix covers every cell.
#' @param bg_n Number of background pairs sampled per target pair.
#'   Default 100.
#' @param bg_mode \code{"partition"} (bin genes by abundance, sample
#'   matches from the same bin) or \code{"local"} (sample matches from
#'   the genes with the closest abundance).
#' @param n_partitions Number of abundance bins for \code{"partition"}
#'   mode.
#' @param local_window Window size for \code{"local"} mode.
#' @param celltype_col Optional meta.data column with cell-type labels.
#'   Required for the composition controls below.
#' @param celltypes Restrict to these cell types only.
#' @param max_celltype_frac Fraction in \code{(0, 1\]}; downsample any
#'   cell type exceeding this fraction.
#' @param center_celltype If TRUE, mean-centre each gene within each
#'   cell type before computing MOC.
#' @param comp_flag_thresh Threshold above which the per-niche
#'   dominant-cell-type-fraction difference flips \code{comp_flag} to
#'   TRUE. Default 0.15.
#' @param min_cells Minimum cells per (sample x niche) to be scored.
#' @param min_samples Minimum samples per condition to run the differential
#'   test for a given niche x pair. Default 2 (the practical floor for both
#'   \code{"wilcox"} and \code{"t"}). With exactly 2 vs 2, the Wilcoxon test
#'   has only a few possible p-values (coarse resolution), and the t-test
#'   has 2 degrees of freedom -- both run, but treat results from very small
#'   groups cautiously. Below 2, \code{stats::sd()} / the tests themselves
#'   are undefined, so values \code{< 2} are treated as 2.
#' @param test \code{"wilcox"} (default) or \code{"t"}.
#' @param p_adjust_method Passed to \code{\link[stats]{p.adjust}}.
#' @param seed RNG seed for reproducible background sampling.
#' @param verbose Print progress.
#'
#' @return A list:
#'   \describe{
#'     \item{\code{per_sample}}{Long data.frame of per-sample x niche x
#'       pair co-expression scores.}
#'     \item{\code{stats}}{Per-niche x pair differential co-expression:
#'       means, \code{delta_log2}, statistic, p-value, BH-adjusted q-value.
#'       If \code{celltype_col} was set, also \code{comp_diff} and
#'       \code{comp_flag}.}
#'     \item{\code{composition}}{(only if \code{celltype_col} set)
#'       Per (sample x niche) cell-type fractions actually used, plus
#'       dominant type and its fraction.}
#'   }
#'
#' @examples
#' \dontrun{
#' panel <- c("Cd3e", "Cd8a", "Pdcd1", "Foxp3", "Ifng", "Gzmb")
#' res <- nicheCoExpress(
#'   seurat_obj    = so,
#'   genes         = panel,
#'   niche_col     = "niche",
#'   sample_col    = "sample",
#'   condition_col = "condition",
#'   conditions    = c("healthy", "tumor")
#' )
#' subset(res$stats, p_adj < 0.05 & delta_log2 > 0)
#' plotNicheCoExpress(res, type = "heatmap")
#' }
#'
#' @seealso \code{\link{plotNicheCoExpress}}
#' @importFrom SeuratObject LayerData Layers JoinLayers
#' @importFrom stats wilcox.test t.test p.adjust quantile sd
#' @importFrom utils combn head
#' @export
NicheCoExpress <- function(seurat_obj,
                           genes,
                           niche_col         = "niche",
                           sample_col        = "sample",
                           condition_col     = "condition",
                           conditions        = NULL,
                           niches            = NULL,
                           assay             = "RNA",
                           layer             = "data",
                           bg_n              = 100,
                           bg_mode           = c("partition", "local"),
                           n_partitions      = 25,
                           local_window      = 100,
                           celltype_col      = NULL,
                           celltypes         = NULL,
                           max_celltype_frac = NULL,
                           center_celltype   = FALSE,
                           comp_flag_thresh  = 0.15,
                           min_cells         = 20,
                           min_samples       = 2,
                           test              = c("wilcox", "t"),
                           p_adjust_method   = "BH",
                           seed              = 1,
                           verbose           = TRUE) {

  bg_mode <- match.arg(bg_mode)
  test    <- match.arg(test)
  if (!is.null(seed)) set.seed(seed)
  if (min_samples < 2) {
    warning("`min_samples` must be >= 2 for the differential test; using 2.")
    min_samples <- 2
  }

  md <- seurat_obj@meta.data
  for (col in c(niche_col, sample_col, condition_col)) {
    if (!col %in% colnames(md)) stop("meta.data has no column '", col, "'.")
  }
  use_ct <- !is.null(celltype_col)
  if (use_ct && !celltype_col %in% colnames(md)) {
    stop("meta.data has no column '", celltype_col, "'.")
  }
  if ((!is.null(celltypes) || !is.null(max_celltype_frac)) && !use_ct) {
    stop("celltypes / max_celltype_frac require celltype_col to be set.")
  }
  if (center_celltype && !use_ct) {
    stop("center_celltype = TRUE requires celltype_col to be set.")
  }
  score_type <- if (center_celltype) "zscore" else "log2ratio"

  # ---- Resolve gene pairs --------------------------------------------------
  if (is.null(dim(genes))) {
    genes <- unique(as.character(genes))
    if (length(genes) < 2) stop("`genes` must contain at least 2 genes.")
    cmb   <- utils::combn(genes, 2)
    pairs <- data.frame(geneA = cmb[1, ], geneB = cmb[2, ],
                        stringsAsFactors = FALSE)
  } else {
    pairs <- as.data.frame(genes, stringsAsFactors = FALSE)[, 1:2]
    colnames(pairs) <- c("geneA", "geneB")
  }
  if (nrow(pairs) == 0) stop("No gene pairs to score.")
  if (verbose) message("Scoring ", nrow(pairs), " gene pair(s).")

  # ---- Resolve conditions --------------------------------------------------
  cond_levels <- if (is.null(conditions)) {
    sort(unique(as.character(md[[condition_col]])))
  } else as.character(conditions)
  if (length(cond_levels) != 2) {
    stop("Exactly two conditions required; got ", length(cond_levels), ".")
  }
  if (verbose) {
    message("Comparing conditions: ", cond_levels[1],
            " (ref) vs ", cond_levels[2], " (test).")
  }

  # ---- Resolve niches ------------------------------------------------------
  if (is.null(niches)) niches <- sort(unique(as.character(md[[niche_col]])))

  # ---- Pull expression matrix ---------------------------------------------
  # If the requested `layer` is split into per-sample pieces (Seurat v5,
  # e.g. 'data.1', 'data.2', ... rather than a single 'data' layer),
  # LayerData() can return just one of those pieces -- a matrix covering
  # only a subset of `colnames(seurat_obj)`. Indexing that matrix below by
  # `cells` (drawn from the *full* meta.data) then fails with "subscript
  # out of bounds" for any cell outside that piece. Join the requested
  # layer's pieces first so `expr_all` covers every cell.
  obj_layers   <- SeuratObject::Layers(seurat_obj[[assay]])
  layer_pieces <- grep(paste0("^", layer, "($|\\.)"), obj_layers, value = TRUE)
  if (length(layer_pieces) > 1) {
    if (verbose) {
      message("Joining ", length(layer_pieces), " '", layer,
              "' layer(s) in assay '", assay, "'.")
    }
    seurat_obj <- tryCatch(
      SeuratObject::JoinLayers(seurat_obj, assay = assay, layers = layer_pieces),
      error = function(e) {
        warning("JoinLayers() failed (", conditionMessage(e),
                "); proceeding with unjoined layers -- some cells may be ",
                "excluded below.")
        seurat_obj
      }
    )
  }

  expr_all <- as.matrix(SeuratObject::LayerData(seurat_obj,
                                                assay = assay,
                                                layer = layer))
  wanted  <- unique(c(pairs$geneA, pairs$geneB))
  missing <- setdiff(wanted, rownames(expr_all))
  if (length(missing)) {
    warning(length(missing), " requested gene(s) not in assay: ",
            paste(utils::head(missing, 10), collapse = ", "),
            if (length(missing) > 10) ", ...")
  }

  # Cells present in meta.data but absent from `expr_all` (e.g. because
  # `assay`/`layer` doesn't have data for every cell, even after the join
  # above) are dropped up front with a warning, rather than failing later
  # when `expr_all` is subset by `cells`.
  expr_cells    <- colnames(expr_all)
  missing_cells <- setdiff(rownames(md), expr_cells)
  if (length(missing_cells)) {
    warning(length(missing_cells), " cell(s) in meta.data have no '", layer,
            "' data in assay '", assay, "' and will be excluded: ",
            paste(utils::head(missing_cells, 5), collapse = ", "),
            if (length(missing_cells) > 5) ", ...")
  }

  # ---- 1) Per (sample x niche) co-expression -------------------------------
  samples_md  <- as.character(md[[sample_col]])
  per_sample  <- list()
  composition <- list()

  for (nz in niches) {
    for (sp in unique(samples_md[as.character(md[[niche_col]]) == nz])) {
      cells <- rownames(md)[as.character(md[[niche_col]]) == nz &
                              samples_md == sp]
      cells <- intersect(cells, expr_cells)

      # Cell-type composition controls
      if (use_ct) {
        ct <- as.character(md[cells, celltype_col])
        names(ct) <- cells
        if (!is.null(celltypes)) ct <- ct[ct %in% celltypes]
        if (length(ct) > 0 && !is.null(max_celltype_frac)) {
          ct <- ct[.balance_cells(ct, max_celltype_frac)]
        }
        cells <- names(ct)
      }

      if (length(cells) < min_cells) {
        if (verbose) message("  skip niche=", nz, " sample=", sp,
                             " (", length(cells), " cells < min_cells)")
        next
      }
      cond_sp <- unique(as.character(md[cells, condition_col]))
      if (length(cond_sp) != 1) {
        warning("Sample ", sp, " spans >1 condition; using the first.")
        cond_sp <- cond_sp[1]
      }
      if (!cond_sp %in% cond_levels) next

      sub <- expr_all[, cells, drop = FALSE]
      ct_sub <- if (use_ct) {
        setNames(as.character(md[cells, celltype_col]), cells)
      } else NULL

      res <- .subset_coexpr(sub, pairs,
                            bg_n         = bg_n,
                            bg_mode      = bg_mode,
                            n_partitions = n_partitions,
                            local_window = local_window,
                            ct           = ct_sub,
                            center       = center_celltype)
      if (is.null(res)) next
      res$niche   <- nz
      res$sample  <- sp
      res$condition <- cond_sp
      res$n_cells <- length(cells)

      if (use_ct) {
        tab  <- sort(table(as.character(md[cells, celltype_col])),
                     decreasing = TRUE)
        frac <- tab / sum(tab)
        comp_row <- data.frame(
          niche             = nz,
          sample            = sp,
          condition         = cond_sp,
          n_cells           = length(cells),
          n_celltypes       = length(tab),
          dominant_celltype = names(frac)[1],
          dominant_frac     = as.numeric(frac[1]),
          stringsAsFactors  = FALSE
        )
        for (k in names(frac)) {
          comp_row[[paste0("frac_", k)]] <- as.numeric(frac[k])
        }
        composition[[paste(nz, sp, sep = "__")]] <- comp_row
        res$dominant_frac <- as.numeric(frac[1])
      }

      per_sample[[paste(nz, sp, sep = "__")]] <- res
    }
  }

  if (length(per_sample) == 0) stop("No sample x niche group passed filters.")
  per_sample <- do.call(rbind, per_sample)
  rownames(per_sample) <- NULL

  comp_tab <- if (use_ct && length(composition)) {
    all_cols <- unique(unlist(lapply(composition, colnames)))
    composition <- lapply(composition, function(x) {
      miss <- setdiff(all_cols, colnames(x))
      for (m in miss) x[[m]] <- 0
      x[, all_cols, drop = FALSE]
    })
    do.call(rbind, composition)
  } else NULL
  if (!is.null(comp_tab)) rownames(comp_tab) <- NULL

  # ---- 2) Differential co-expression --------------------------------------
  per_sample$pair <- paste(per_sample$geneA, per_sample$geneB, sep = "_")
  stats_list <- list()

  for (nz in unique(per_sample$niche)) {
    sub_nz <- per_sample[per_sample$niche == nz, ]
    for (pr in unique(sub_nz$pair)) {
      d  <- sub_nz[sub_nz$pair == pr & !is.na(sub_nz$coexpr), ]
      v1 <- d$coexpr[d$condition == cond_levels[1]]
      v2 <- d$coexpr[d$condition == cond_levels[2]]

      n1 <- length(v1); n2 <- length(v2)
      enough <- n1 >= min_samples && n2 >= min_samples

      pval <- NA_real_; stat <- NA_real_
      if (enough) {
        # With small (e.g. 2 vs 2) groups, wilcox.test() commonly warns
        # "cannot compute exact p-value with ties" and falls back to a
        # normal approximation; t.test() can warn about near-zero variance.
        # Both are expected at this sample size and not actionable, so they
        # are suppressed -- a failure still returns NULL via tryCatch and
        # leaves p-value/statistic as NA.
        tt <- tryCatch(
          suppressWarnings(
            if (test == "wilcox") stats::wilcox.test(v2, v1)
            else                  stats::t.test(v2, v1)
          ),
          error = function(e) NULL
        )
        if (!is.null(tt)) {
          pval <- tt$p.value
          stat <- unname(tt$statistic)
        }
      }

      stats_list[[paste(nz, pr, sep = "__")]] <- data.frame(
        niche      = nz,
        geneA      = d$geneA[1],
        geneB      = d$geneB[1],
        pair       = pr,
        n_ref      = n1,
        n_test     = n2,
        mean_ref   = if (n1 > 0) mean(v1) else NA_real_,
        mean_test  = if (n2 > 0) mean(v2) else NA_real_,
        delta_log2 = if (n1 > 0 && n2 > 0) mean(v2) - mean(v1) else NA_real_,
        statistic  = stat,
        p_value    = pval,
        stringsAsFactors = FALSE
      )
    }
  }

  stats_df <- do.call(rbind, stats_list)
  rownames(stats_df) <- NULL

  # Multiple-testing correction within niche
  stats_df$p_adj <- NA_real_
  for (nz in unique(stats_df$niche)) {
    idx <- stats_df$niche == nz & !is.na(stats_df$p_value)
    stats_df$p_adj[idx] <- stats::p.adjust(stats_df$p_value[idx],
                                           method = p_adjust_method)
  }

  # Composition confound indicator
  if (!is.null(comp_tab)) {
    comp_diff_by_niche <- vapply(unique(stats_df$niche), function(nz) {
      d  <- comp_tab[comp_tab$niche == nz, ]
      m1 <- mean(d$dominant_frac[d$condition == cond_levels[1]], na.rm = TRUE)
      m2 <- mean(d$dominant_frac[d$condition == cond_levels[2]], na.rm = TRUE)
      if (is.nan(m1) || is.nan(m2)) NA_real_ else abs(m2 - m1)
    }, numeric(1))
    names(comp_diff_by_niche) <- unique(stats_df$niche)
    stats_df$comp_diff <- comp_diff_by_niche[as.character(stats_df$niche)]
    stats_df$comp_flag <- !is.na(stats_df$comp_diff) &
                          stats_df$comp_diff > comp_flag_thresh
  }

  stats_df <- stats_df[order(stats_df$niche, stats_df$p_adj), ]

  attr(stats_df,   "conditions") <- cond_levels
  attr(stats_df,   "score_type") <- score_type
  attr(per_sample, "score_type") <- score_type
  out <- list(per_sample = per_sample, stats = stats_df)
  if (!is.null(comp_tab)) out$composition <- comp_tab
  out
}


#' Visualise the output of nicheCoExpress
#'
#' Two views:
#' \describe{
#'   \item{\code{type = "heatmap"} (default)}{Niche (rows) by gene-pair
#'     (cols) heatmap of \code{delta_log2} &mdash; the change in
#'     co-expression between the two conditions (test minus reference).
#'     Diverging palette centred at 0: red = gained co-expression in the
#'     test condition, blue = lost. Significance stars are drawn from
#'     \code{p_adj}.}
#'   \item{\code{type = "scores"}}{Per-sample co-expression scores for
#'     one or more gene pairs, split by condition and faceted by niche.
#'     Box + jitter so you can see the replicate-level spread the
#'     statistics are testing.}
#' }
#'
#' @param res The list returned by \code{\link{nicheCoExpress}}.
#' @param type \code{"heatmap"} or \code{"scores"}.
#' @param pairs (scores mode) character vector of pair IDs
#'   (\code{"geneA_geneB"}) to show; default = top \code{top_n} most
#'   significant.
#' @param niches Optional subset of niches to display.
#' @param top_n In heatmap: max pairs to show, ranked by smallest p_adj;
#'   in scores: number of top pairs if \code{pairs} is NULL.
#' @param sig_levels Named thresholds for significance stars.
#' @return A \code{ggplot} object.
#' @seealso \code{\link{nicheCoExpress}}
#' @importFrom ggplot2 ggplot aes geom_tile geom_text geom_hline geom_boxplot geom_jitter scale_fill_gradient2 facet_grid labs theme_minimal theme_bw theme element_text element_blank
#' @export
plotNicheCoExpress <- function(res,
                               type       = c("heatmap", "scores"),
                               pairs      = NULL,
                               niches     = NULL,
                               top_n      = 30,
                               sig_levels = c(`***` = 0.001,
                                              `**`  = 0.01,
                                              `*`   = 0.05)) {
  type <- match.arg(type)
  stats_df <- res$stats
  conds <- attr(stats_df, "conditions")
  if (is.null(conds)) conds <- c("reference", "test")
  score_type <- attr(stats_df, "score_type")
  if (is.null(score_type)) score_type <- "log2ratio"
  is_z <- identical(score_type, "zscore")
  fill_lab  <- if (is_z) expression(Delta ~ "co-expr (z)")
               else      expression(Delta * " log"[2] ~ "co-expr")
  ylab      <- if (is_z) "within-cell-type co-expr (z vs background)"
               else      expression("log"[2] ~ "(observed MOC / background MOC)")
  zero_note <- if (is_z) "0 = matched random within-type correlation"
               else      "0 = at background co-expression"

  star <- function(p) {
    vapply(p, function(pp) {
      if (is.na(pp)) return("")
      hit <- names(sig_levels)[pp <= sig_levels]
      if (length(hit)) hit[which.min(sig_levels[hit])] else ""
    }, character(1))
  }

  if (!is.null(niches)) stats_df <- stats_df[stats_df$niche %in% niches, ]
  if (nrow(stats_df) == 0) stop("No rows to plot after filtering.")

  if (type == "heatmap") {
    ord <- stats_df[!is.na(stats_df$p_adj), ]
    pair_rank  <- tapply(ord$p_adj, ord$pair, min, na.rm = TRUE)
    keep_pairs <- names(sort(pair_rank))[seq_len(min(top_n, length(pair_rank)))]
    if (!is.null(pairs)) keep_pairs <- intersect(pairs, unique(stats_df$pair))
    df <- stats_df[stats_df$pair %in% keep_pairs, ]
    df$label <- star(df$p_adj)

    pair_order <- names(sort(tapply(df$delta_log2, df$pair,
                                    function(x) mean(x, na.rm = TRUE))))
    df$pair  <- factor(df$pair,  levels = pair_order)
    df$niche <- factor(df$niche, levels = sort(unique(df$niche)))

    lim <- max(abs(df$delta_log2), na.rm = TRUE)
    pair <- niche <- delta_log2 <- label <- NULL  # NSE silencing
    ggplot2::ggplot(df, ggplot2::aes(x = pair, y = niche, fill = delta_log2)) +
      ggplot2::geom_tile(color = "grey90") +
      ggplot2::geom_text(ggplot2::aes(label = label), size = 4,
                         vjust = 0.78, color = "black") +
      ggplot2::scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0, limits = c(-lim, lim), name = fill_lab) +
      ggplot2::labs(
        x = "gene pair", y = "niche",
        title = paste0("Differential co-expression: ",
                       conds[2], " vs ", conds[1]),
        subtitle = "red = higher in test condition  |  * p_adj<.05  ** <.01  *** <.001"
      ) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid  = ggplot2::element_blank())

  } else {
    ps <- res$per_sample
    if (!is.null(niches)) ps <- ps[ps$niche %in% niches, ]
    if (is.null(pairs)) {
      ord <- stats_df[!is.na(stats_df$p_adj), ]
      pair_rank <- tapply(ord$p_adj, ord$pair, min, na.rm = TRUE)
      pairs <- names(sort(pair_rank))[seq_len(min(top_n, length(pair_rank)))]
    }
    df <- ps[ps$pair %in% pairs & !is.na(ps$coexpr), ]
    if (nrow(df) == 0) stop("No per-sample scores to plot after filtering.")
    df$condition <- factor(df$condition, levels = conds)
    df$pair      <- factor(df$pair,      levels = pairs)

    condition <- coexpr <- pair <- niche <- NULL  # NSE silencing
    ggplot2::ggplot(df, ggplot2::aes(x = condition, y = coexpr,
                                     fill = condition)) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.6) +
      ggplot2::geom_jitter(width = 0.15, height = 0, size = 1.8,
                           alpha = 0.8, ggplot2::aes(color = condition)) +
      ggplot2::facet_grid(pair ~ niche, switch = "y") +
      ggplot2::labs(
        x = NULL, y = ylab,
        title = "Per-sample co-expression by niche",
        subtitle = paste0("each point = one sample; ", zero_note)
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(legend.position   = "none",
                     strip.text.y.left = ggplot2::element_text(angle = 0))
  }
}


# ============================================================================
# Internal helpers
# ============================================================================

#' Observed + background MOC for a set of target gene pairs in one subset
#'
#' Internal helper for \code{nicheCoExpress}. Not exported.
#' @keywords internal
#' @noRd
.subset_coexpr <- function(expr, pairs, bg_n = 100,
                           bg_mode = c("partition", "local"),
                           n_partitions = 25, local_window = 100,
                           min_expr_genes = 5,
                           ct = NULL, center = FALSE) {
  bg_mode <- match.arg(bg_mode)

  keep <- rowSums(expr) > 0
  expr <- expr[keep, , drop = FALSE]
  if (nrow(expr) < min_expr_genes) return(NULL)

  mu <- rowMeans(expr)

  if (center) {
    if (is.null(ct)) stop("center = TRUE requires cell-type labels (ct).")
    ct   <- ct[colnames(expr)]
    cmat <- expr
    for (t in unique(ct)) {
      cols <- which(ct == t)
      if (length(cols) >= 2) {
        m <- rowMeans(cmat[, cols, drop = FALSE])
        cmat[, cols] <- cmat[, cols] - m
      } else {
        cmat[, cols] <- 0
      }
    }
    mat <- cmat
  } else {
    mat <- expr
  }

  l2 <- sqrt(rowSums(mat^2))
  l2[l2 == 0] <- 1
  norm_expr <- mat / l2

  moc <- function(gA, gB) sum(norm_expr[gA, ] * norm_expr[gB, ])

  if (bg_mode == "partition") {
    brks <- unique(stats::quantile(mu, probs = seq(0, 1, length.out = n_partitions + 1)))
    bins <- cut(mu, breaks = brks, include.lowest = TRUE, labels = FALSE)
    names(bins) <- names(mu)
  }

  candidates <- function(g) {
    if (bg_mode == "partition") {
      cand <- names(bins)[bins == bins[[g]]]
    } else {
      d <- abs(mu - mu[[g]])
      cand <- names(sort(d))[seq_len(min(local_window + 1, length(d)))]
    }
    setdiff(cand, g)
  }

  out <- vector("list", nrow(pairs))
  for (i in seq_len(nrow(pairs))) {
    gA <- pairs$geneA[i]; gB <- pairs$geneB[i]
    if (!gA %in% rownames(expr) || !gB %in% rownames(expr)) {
      out[[i]] <- data.frame(geneA = gA, geneB = gB,
                             moc_obs = NA_real_, moc_bg = NA_real_,
                             bg_sd = NA_real_, coexpr = NA_real_,
                             stringsAsFactors = FALSE)
      next
    }
    moc_obs <- moc(gA, gB)
    cA <- setdiff(candidates(gA), c(gA, gB))
    cB <- setdiff(candidates(gB), c(gA, gB))
    if (length(cA) == 0 || length(cB) == 0) {
      moc_bg <- NA_real_; bg_sd <- NA_real_
    } else {
      sA <- sample(cA, bg_n, replace = TRUE)
      sB <- sample(cB, bg_n, replace = TRUE)
      ok <- sA != sB
      sA <- sA[ok]; sB <- sB[ok]
      bg_vals <- vapply(seq_along(sA), function(j) moc(sA[j], sB[j]), numeric(1))
      moc_bg <- mean(bg_vals, na.rm = TRUE)
      bg_sd  <- stats::sd(bg_vals,  na.rm = TRUE)
    }
    coexpr <- if (center) {
      if (is.na(moc_bg) || is.na(bg_sd) || bg_sd == 0) NA_real_
      else (moc_obs - moc_bg) / bg_sd
    } else {
      if (is.na(moc_bg) || moc_bg <= 0) NA_real_
      else log2(moc_obs / moc_bg)
    }
    out[[i]] <- data.frame(geneA = gA, geneB = gB, moc_obs = moc_obs,
                           moc_bg = moc_bg, bg_sd = bg_sd, coexpr = coexpr,
                           stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}


#' Cap the dominant cell type within a set of cells by downsampling
#'
#' Internal helper. Returns the kept cell IDs.
#' @keywords internal
#' @noRd
.balance_cells <- function(ct, max_frac) {
  n <- length(ct)
  if (is.null(max_frac) || n == 0) return(names(ct))
  cap <- floor(max_frac * n)
  if (cap < 1) return(names(ct))
  keep <- character(0)
  for (lvl in unique(ct)) {
    ids <- names(ct)[ct == lvl]
    if (length(ids) > cap) ids <- sample(ids, cap)
    keep <- c(keep, ids)
  }
  keep
}
