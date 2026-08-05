# Tests for R/NicheCoExpress.R: NicheCoExpress(), plotNicheCoExpress(), and
# the internal helpers .subset_coexpr() / .balance_cells(). See
# helper-NicheCoExpress.R for the shared synthetic-object fixture.

# ============================================================================
# .balance_cells() -- internal cell-type-composition downsampler
# ============================================================================

test_that(".balance_cells satisfies max_frac at a fixed point (regression)", {
  # 80/20 split at max_frac = 0.5 used to cap once against the *original*
  # count (80 -> 50, but 50 of the remaining 70 cells is 71%, still over
  # the cap). The fix iterates the cap against the shrinking total; here
  # it should land at exactly 50/50.
  ct <- setNames(c(rep("A", 80), rep("B", 20)), paste0("c", 1:100))
  keep <- SingleCellTools:::.balance_cells(ct, max_frac = 0.5)
  tab <- table(ct[keep])
  fracs <- tab / sum(tab)
  expect_true(all(fracs <= 0.5 + 1e-9))
  expect_equal(unname(tab[["A"]]), unname(tab[["B"]]))
})

test_that(".balance_cells only caps the cell type(s) actually over the limit", {
  ct <- setNames(c(rep("A", 80), rep("B", 15), rep("C", 5)), paste0("c", 1:100))
  keep <- SingleCellTools:::.balance_cells(ct, max_frac = 0.5)
  tab <- table(ct[keep])
  expect_true(max(tab) / sum(tab) <= 0.5 + 1e-9)
  expect_equal(unname(tab[["B"]]), 15)
  expect_equal(unname(tab[["C"]]), 5)
})

test_that(".balance_cells is a no-op with a single cell type present", {
  ct <- setNames(rep("A", 50), paste0("c", 1:50))
  keep <- SingleCellTools:::.balance_cells(ct, max_frac = 0.1)
  expect_length(keep, 50)
})

test_that(".balance_cells warns and skips downsampling when max_frac is too restrictive", {
  ct <- setNames(c(rep("A", 34), rep("B", 33), rep("C", 33)), paste0("c", 1:100))
  expect_warning(
    keep <- SingleCellTools:::.balance_cells(ct, max_frac = 0.05),
    "too restrictive"
  )
  expect_length(keep, 100)
})

test_that(".balance_cells returns all cells unchanged when max_frac is NULL", {
  ct <- setNames(c(rep("A", 80), rep("B", 20)), paste0("c", 1:100))
  keep <- SingleCellTools:::.balance_cells(ct, max_frac = NULL)
  expect_setequal(keep, names(ct))
})

test_that(".balance_cells leaves an already-balanced set untouched", {
  ct <- setNames(c(rep("A", 10), rep("B", 10), rep("C", 10)), paste0("c", 1:30))
  keep <- SingleCellTools:::.balance_cells(ct, max_frac = 0.5)
  expect_length(keep, 30)
})


# ============================================================================
# .subset_coexpr() -- internal per-(sample x niche) MOC scorer
# ============================================================================

test_that(".subset_coexpr scores a genuinely co-expressed pair above an independent one", {
  set.seed(42)
  n <- 200
  # All genes (signal + filler) are drawn from the same abundance level so
  # partition-mode binning has plenty of candidates per bin; only the
  # *correlation* structure differs between G1/G2 (shared latent factor)
  # and the rest (independent draws).
  base <- rgamma(n, shape = 3, rate = 1)
  g1 <- pmax(base + rnorm(n, sd = 0.1), 0.01)
  g2 <- pmax(base + rnorm(n, sd = 0.1), 0.01)           # co-expressed with g1
  g3 <- pmax(rgamma(n, shape = 3, rate = 1), 0.01)       # independent
  filler <- t(replicate(30, pmax(rgamma(n, shape = 3, rate = 1), 0.01)))
  rownames(filler) <- paste0("F", seq_len(nrow(filler)))
  expr <- rbind(G1 = g1, G2 = g2, G3 = g3, filler)
  colnames(expr) <- paste0("cell", seq_len(n))

  pairs <- data.frame(geneA = c("G1", "G3"), geneB = c("G2", "G1"),
                      stringsAsFactors = FALSE)
  res <- SingleCellTools:::.subset_coexpr(expr, pairs, bg_n = 100,
                                          bg_mode = "partition",
                                          n_partitions = 4)
  g1g2 <- res[res$geneA == "G1" & res$geneB == "G2", ]
  g3g1 <- res[res$geneA == "G3" & res$geneB == "G1", ]
  expect_false(is.na(g1g2$coexpr))
  expect_false(is.na(g3g1$coexpr))
  expect_true(g1g2$coexpr > g3g1$coexpr)
  expect_true(g1g2$coexpr > 0)  # above background
})

test_that(".subset_coexpr returns NA scores for a gene absent from the expression matrix", {
  expr <- matrix(rgamma(6 * 50, 2), nrow = 6,
                dimnames = list(paste0("G", 1:6), paste0("c", 1:50)))
  pairs <- data.frame(geneA = "G1", geneB = "Ghost", stringsAsFactors = FALSE)
  res <- SingleCellTools:::.subset_coexpr(expr, pairs, bg_n = 20)
  expect_true(is.na(res$coexpr))
  expect_true(is.na(res$moc_obs))
})

test_that(".subset_coexpr returns NULL when too few expressed genes remain", {
  expr <- matrix(0, nrow = 3, ncol = 20,
                dimnames = list(paste0("G", 1:3), paste0("c", 1:20)))
  expr["G1", ] <- rgamma(20, 2)  # only 1 gene with nonzero expression
  pairs <- data.frame(geneA = "G1", geneB = "G2", stringsAsFactors = FALSE)
  res <- SingleCellTools:::.subset_coexpr(expr, pairs, bg_n = 20, min_expr_genes = 5)
  expect_null(res)
})

test_that(".subset_coexpr requires cell-type labels when center = TRUE", {
  expr <- matrix(rgamma(6 * 50, 2), nrow = 6,
                dimnames = list(paste0("G", 1:6), paste0("c", 1:50)))
  pairs <- data.frame(geneA = "G1", geneB = "G2", stringsAsFactors = FALSE)
  expect_error(
    SingleCellTools:::.subset_coexpr(expr, pairs, bg_n = 20, center = TRUE, ct = NULL),
    "cell-type"
  )
})

test_that(".subset_coexpr with center = TRUE returns a finite z-score-like value", {
  set.seed(7)
  # Only 6 genes total, so partition-mode binning (default n_partitions =
  # 25) would put most genes alone in their own bin with no candidates to
  # sample a background from. "local" mode ranks by abundance distance
  # instead and, with only 5 other genes to choose from, naturally uses
  # all of them as candidates regardless of window size.
  expr <- matrix(rgamma(6 * 60, 2), nrow = 6,
                dimnames = list(paste0("G", 1:6), paste0("c", 1:60)))
  ct <- setNames(rep(c("Tcell", "Myeloid"), each = 30), colnames(expr))
  pairs <- data.frame(geneA = "G1", geneB = "G2", stringsAsFactors = FALSE)
  res <- SingleCellTools:::.subset_coexpr(expr, pairs, bg_n = 30, bg_mode = "local",
                                          ct = ct, center = TRUE)
  expect_false(is.na(res$coexpr))
})


# ============================================================================
# NicheCoExpress() -- exported, end-to-end on a synthetic Seurat object
# ============================================================================

test_that("NicheCoExpress runs end-to-end and returns the documented structure", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1)
  res <- suppressWarnings(suppressMessages(
    NicheCoExpress(
      obj, genes = c("G1", "G2", "G3"),
      niche_col = "niche", sample_col = "sample", condition_col = "condition",
      conditions = c("healthy", "tumor"), min_cells = 5, verbose = FALSE,
      bg_mode = "local"  # only 6 genes in the fixture; see helper file note
    )
  ))
  expect_type(res, "list")
  expect_true(all(c("per_sample", "stats") %in% names(res)))
  expect_true(all(c("geneA", "geneB", "coexpr", "niche", "sample", "condition")
                 %in% colnames(res$per_sample)))
  expect_true(all(c("niche", "pair", "delta", "p_value", "p_adj")
                 %in% colnames(res$stats)))
  expect_equal(attr(res$stats, "conditions"), c("healthy", "tumor"))
})

test_that("NicheCoExpress requires at least 2 genes when given a gene vector", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1)
  expect_error(
    NicheCoExpress(obj, genes = "G1", conditions = c("healthy", "tumor"),
                   verbose = FALSE),
    "at least 2 genes"
  )
})

test_that("NicheCoExpress errors clearly on a missing meta.data column", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1)
  expect_error(
    NicheCoExpress(obj, genes = c("G1", "G2"), niche_col = "not_a_column",
                   verbose = FALSE),
    "not_a_column"
  )
})

test_that("NicheCoExpress requires exactly two conditions", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1)
  obj$condition3 <- obj$condition
  obj$condition3[1] <- "other"  # injects a 3rd distinct level
  expect_error(
    NicheCoExpress(obj, genes = c("G1", "G2"), condition_col = "condition3",
                   verbose = FALSE),
    "Exactly two conditions"
  )
})

test_that("NicheCoExpress drops self-pairs from a custom pair table (regression)", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1)
  pairs <- data.frame(geneA = c("G1", "G3"), geneB = c("G1", "G4"),
                      stringsAsFactors = FALSE)  # first row is a self-pair
  expect_warning(
    res <- suppressMessages(
      NicheCoExpress(obj, genes = pairs, conditions = c("healthy", "tumor"),
                     min_cells = 5, verbose = FALSE, bg_mode = "local")
    ),
    "self-pair"
  )
  expect_false(any(res$stats$geneA == res$stats$geneB))
  expect_true(all(res$per_sample$pair == "G3_G4"))
})

test_that("NicheCoExpress errors clearly when the requested layer doesn't exist (regression)", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1)
  expect_error(
    NicheCoExpress(obj, genes = c("G1", "G2"), layer = "nonexistent_layer",
                   conditions = c("healthy", "tumor"), verbose = FALSE),
    "not found"
  )
})

test_that("NicheCoExpress resolves a sample spanning >1 condition via majority vote (regression)", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1)
  # Corrupt a minority of S1's cells (normally all 'healthy') to 'tumor' so
  # S1 spans both conditions within niche N1, with 'healthy' still the
  # majority -- majority-vote resolution should keep it 'healthy'.
  md <- obj@meta.data
  s1_n1 <- rownames(md)[md$sample == "S1" & md$niche == "N1"]
  flip <- sample(s1_n1, 3)
  obj$condition[flip] <- "tumor"

  expect_warning(
    res <- suppressMessages(
      NicheCoExpress(obj, genes = c("G1", "G2"),
                     conditions = c("healthy", "tumor"),
                     min_cells = 5, verbose = FALSE, bg_mode = "local")
    ),
    "spans >1 condition"
  )
  row <- res$per_sample[res$per_sample$niche == "N1" & res$per_sample$sample == "S1", ]
  expect_equal(unique(row$condition), "healthy")
})

test_that("NicheCoExpress requires celltype_col when max_celltype_frac is set", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1)
  expect_error(
    NicheCoExpress(obj, genes = c("G1", "G2"), max_celltype_frac = 0.5,
                   conditions = c("healthy", "tumor"), verbose = FALSE),
    "celltype_col"
  )
})

test_that("NicheCoExpress center_celltype requires celltype_col", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1)
  expect_error(
    NicheCoExpress(obj, genes = c("G1", "G2"), center_celltype = TRUE,
                   conditions = c("healthy", "tumor"), verbose = FALSE),
    "celltype_col"
  )
})

test_that("NicheCoExpress omits composition output when celltype_col is not set", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1, with_celltypes = FALSE)
  res <- suppressWarnings(suppressMessages(
    NicheCoExpress(obj, genes = c("G1", "G2"),
                   conditions = c("healthy", "tumor"), min_cells = 5, verbose = FALSE,
                   bg_mode = "local")
  ))
  expect_false("composition" %in% names(res))
})

test_that("NicheCoExpress composition output respects max_celltype_frac (regression)", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 2, with_celltypes = FALSE)
  # Force an 80/20 imbalance per (niche x sample) group so the composition
  # cap fixed in .balance_cells() is actually exercised end-to-end.
  md <- obj@meta.data
  ct <- character(nrow(md))
  for (grp in unique(paste(md$niche, md$sample))) {
    idx <- which(paste(md$niche, md$sample) == grp)
    n <- length(idx)
    n_t <- round(n * 0.8)
    ct[idx] <- sample(c(rep("Tcell", n_t), rep("Myeloid", n - n_t)))
  }
  obj$celltype <- ct

  res <- suppressWarnings(suppressMessages(
    NicheCoExpress(
      obj, genes = c("G1", "G2"), celltype_col = "celltype",
      max_celltype_frac = 0.5, conditions = c("healthy", "tumor"),
      min_cells = 5, verbose = FALSE, bg_mode = "local"
    )
  ))
  expect_true("composition" %in% names(res))
  expect_true(all(res$composition$dominant_frac <= 0.5 + 1e-9))
})

test_that("NicheCoExpress warns and coerces min_samples < 2 to 2", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1)
  expect_warning(
    res <- suppressMessages(
      NicheCoExpress(obj, genes = c("G1", "G2"),
                     conditions = c("healthy", "tumor"),
                     min_cells = 5, min_samples = 1, verbose = FALSE,
                     bg_mode = "local")
    ),
    "min_samples"
  )
  expect_type(res, "list")
})

test_that("NicheCoExpress supports test = 't' as an alternative to Wilcoxon", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1)
  res <- suppressWarnings(suppressMessages(
    NicheCoExpress(obj, genes = c("G1", "G2"),
                   conditions = c("healthy", "tumor"),
                   min_cells = 5, verbose = FALSE, bg_mode = "local",
                   test = "t")
  ))
  expect_true(any(!is.na(res$stats$p_value)))
})

test_that("NicheCoExpress p_adjust_scope = 'niche' (default) corrects within each niche separately (regression)", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1)
  res <- suppressWarnings(suppressMessages(
    NicheCoExpress(obj, genes = c("G1", "G2", "G3"),
                   conditions = c("healthy", "tumor"),
                   min_cells = 5, verbose = FALSE, bg_mode = "local")
  ))
  expect_equal(attr(res$stats, "p_adjust_scope"), "niche")
  for (nz in unique(res$stats$niche)) {
    idx <- res$stats$niche == nz & !is.na(res$stats$p_value)
    expect_equal(res$stats$p_adj[idx],
                stats::p.adjust(res$stats$p_value[idx], method = "BH"))
  }
})

test_that("NicheCoExpress p_adjust_scope = 'global' corrects jointly across all niche x pair tests", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1)
  res <- suppressWarnings(suppressMessages(
    NicheCoExpress(obj, genes = c("G1", "G2", "G3"),
                   conditions = c("healthy", "tumor"),
                   min_cells = 5, verbose = FALSE, bg_mode = "local",
                   p_adjust_scope = "global")
  ))
  expect_equal(attr(res$stats, "p_adjust_scope"), "global")
  idx <- !is.na(res$stats$p_value)
  expect_equal(res$stats$p_adj[idx],
              stats::p.adjust(res$stats$p_value[idx], method = "BH"))
})

test_that("NicheCoExpress reports comp_diff/comp_flag on stats when celltype_col is set (regression)", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1, with_celltypes = TRUE)
  res <- suppressWarnings(suppressMessages(
    NicheCoExpress(obj, genes = c("G1", "G2"), celltype_col = "celltype",
                   conditions = c("healthy", "tumor"), min_cells = 5, verbose = FALSE,
                   bg_mode = "local", comp_flag_thresh = 0.15)
  ))
  expect_true(all(c("comp_diff", "comp_flag") %in% colnames(res$stats)))
  expect_type(res$stats$comp_flag, "logical")
  present <- !is.na(res$stats$comp_diff)
  expect_equal(res$stats$comp_flag[present],
              res$stats$comp_diff[present] > 0.15)
})

test_that("NicheCoExpress sets score_type to log2ratio by default and zscore when center_celltype = TRUE", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1, with_celltypes = TRUE)

  res_log2 <- suppressWarnings(suppressMessages(
    NicheCoExpress(obj, genes = c("G1", "G2"),
                   conditions = c("healthy", "tumor"), min_cells = 5, verbose = FALSE,
                   bg_mode = "local")
  ))
  expect_equal(attr(res_log2$stats, "score_type"), "log2ratio")
  expect_equal(attr(res_log2$per_sample, "score_type"), "log2ratio")

  res_z <- suppressWarnings(suppressMessages(
    NicheCoExpress(obj, genes = c("G1", "G2"), celltype_col = "celltype",
                   center_celltype = TRUE,
                   conditions = c("healthy", "tumor"), min_cells = 5, verbose = FALSE,
                   bg_mode = "local")
  ))
  expect_equal(attr(res_z$stats, "score_type"), "zscore")
  expect_equal(attr(res_z$per_sample, "score_type"), "zscore")
})


# ============================================================================
# plotNicheCoExpress() -- visualisation
# ============================================================================

test_that("plotNicheCoExpress returns a ggplot for both heatmap and scores modes", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1)
  res <- suppressWarnings(suppressMessages(
    NicheCoExpress(obj, genes = c("G1", "G2", "G3"),
                   conditions = c("healthy", "tumor"), min_cells = 5, verbose = FALSE,
                   bg_mode = "local")
  ))
  p1 <- plotNicheCoExpress(res, type = "heatmap")
  p2 <- plotNicheCoExpress(res, type = "scores")
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("plotNicheCoExpress errors when the niche filter leaves nothing to plot", {
  .skip_if_no_seurat()
  obj <- .make_coexpr_object(seed = 1)
  res <- suppressWarnings(suppressMessages(
    NicheCoExpress(obj, genes = c("G1", "G2"),
                   conditions = c("healthy", "tumor"), min_cells = 5, verbose = FALSE,
                   bg_mode = "local")
  ))
  expect_error(
    plotNicheCoExpress(res, niches = "NotARealNiche"),
    "No rows"
  )
})


# ============================================================================
# NicheCoExpressEstimationPlot() -- dabestr-backed effect-size plots
# ============================================================================

test_that("NicheCoExpressEstimationPlot returns a single plot for one niche x pair", {
  .skip_if_no_seurat()
  .skip_if_missing("dabestr")
  obj <- .make_coexpr_object(seed = 1)
  res <- suppressWarnings(suppressMessages(
    NicheCoExpress(obj, genes = c("G1", "G2"),
                   conditions = c("healthy", "tumor"), min_cells = 5, verbose = FALSE,
                   bg_mode = "local")
  ))
  p <- NicheCoExpressEstimationPlot(res, niches = "N1", pairs = "G1_G2")
  expect_false(is.null(p))
  expect_false(identical(names(p), "N1 | G1_G2"))
})

test_that("NicheCoExpressEstimationPlot returns a named list for multiple niche x pair combinations", {
  .skip_if_no_seurat()
  .skip_if_missing("dabestr")
  obj <- .make_coexpr_object(seed = 1)
  res <- suppressWarnings(suppressMessages(
    NicheCoExpress(obj, genes = c("G1", "G2"),
                   conditions = c("healthy", "tumor"), min_cells = 5, verbose = FALSE,
                   bg_mode = "local")
  ))
  plots <- NicheCoExpressEstimationPlot(res)
  expect_true(is.list(plots))
  expect_setequal(names(plots), c("N1 | G1_G2", "N2 | G1_G2"))
})

test_that("NicheCoExpressEstimationPlot defaults idx to attr(res$stats, 'conditions'), not alphabetical", {
  .skip_if_no_seurat()
  .skip_if_missing("dabestr")
  obj <- .make_coexpr_object(seed = 1)
  # Conditions requested in reverse-alphabetical order on purpose: if the
  # idx default silently fell back to alphabetical sorting instead of
  # reusing NicheCoExpress()'s own resolved order, this would be the case
  # that catches it (see R/NicheCoExpressEstimationPlot.R Details).
  res <- suppressWarnings(suppressMessages(
    NicheCoExpress(obj, genes = c("G1", "G2"),
                   conditions = c("tumor", "healthy"), min_cells = 5, verbose = FALSE,
                   bg_mode = "local")
  ))
  expect_equal(attr(res$stats, "conditions"), c("tumor", "healthy"))

  msgs <- testthat::capture_messages(
    p <- NicheCoExpressEstimationPlot(res, niches = "N1", pairs = "G1_G2")
  )
  expect_false(any(grepl("not supplied", msgs)))
  expect_false(is.null(p))
})

test_that("NicheCoExpressEstimationPlot errors when niches/pairs filter leaves nothing", {
  .skip_if_no_seurat()
  .skip_if_missing("dabestr")
  obj <- .make_coexpr_object(seed = 1)
  res <- suppressWarnings(suppressMessages(
    NicheCoExpress(obj, genes = c("G1", "G2"),
                   conditions = c("healthy", "tumor"), min_cells = 5, verbose = FALSE,
                   bg_mode = "local")
  ))
  expect_error(
    NicheCoExpressEstimationPlot(res, niches = "NotARealNiche"),
    "No per-sample scores"
  )
})

test_that("NicheCoExpressEstimationPlot errors on a malformed res", {
  expect_error(
    NicheCoExpressEstimationPlot(list(per_sample = data.frame())),
    "NicheCoExpress"
  )
  expect_error(
    NicheCoExpressEstimationPlot(data.frame(x = 1)),
    "NicheCoExpress"
  )
})
