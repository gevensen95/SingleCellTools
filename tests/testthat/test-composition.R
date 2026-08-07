# Tests for CellComposition(), CompositionAnalysis(), CompositionBarplot(),
# and CompositionalTest(). All use the shared .make_small_seurat() fixture,
# which already carries seurat_clusters / sample / condition / celltype
# metadata columns.

# ============================================================================
# CellComposition()
# ============================================================================

test_that("CellComposition style = 'none' returns per-sample proportions summing to 1", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 100, n_clusters = 3, n_samples = 4)
  df <- CellComposition(obj, cluster_col = "seurat_clusters",
                        sample_col = "sample", style = "none")
  expect_true(all(c("sample", "cluster", "n_cells", "n_sample_total", "prop")
                 %in% colnames(df)))
  totals <- stats::aggregate(prop ~ sample, data = df, FUN = sum)
  expect_equal(totals$prop, rep(1, nrow(totals)), tolerance = 1e-8)
})

test_that("CellComposition normalize = 'cluster' sums to 1 per cluster instead", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 100, n_clusters = 3, n_samples = 4)
  df <- CellComposition(obj, cluster_col = "seurat_clusters", sample_col = "sample",
                        style = "none", normalize = "cluster")
  totals <- stats::aggregate(prop ~ cluster, data = df, FUN = sum)
  expect_equal(totals$prop, rep(1, nrow(totals)), tolerance = 1e-8)
})

test_that("CellComposition style = 'stack' returns a list(df, plot)", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 60)
  out <- CellComposition(obj, sample_col = "sample", style = "stack")
  expect_named(out, c("df", "plot"))
  expect_s3_class(out$plot, "ggplot")
})

test_that("CellComposition style = 'box'/'line' require group_col", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 60)
  expect_error(CellComposition(obj, sample_col = "sample", style = "box"),
              "group_col")
  expect_error(CellComposition(obj, sample_col = "sample", style = "line"),
              "group_col")
})

test_that("CellComposition style = 'box' works with group_col supplied", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 60)
  out <- CellComposition(obj, sample_col = "sample", group_col = "condition",
                         style = "box")
  expect_s3_class(out$plot, "ggplot")
  expect_true("group" %in% colnames(out$df))
})

test_that("CellComposition errors on missing columns", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  expect_error(CellComposition(obj, cluster_col = "nope", sample_col = "sample"),
              "nope")
  expect_error(CellComposition(obj, sample_col = "sample", group_col = "nope",
                               style = "none"),
              "nope")
})

test_that("CellComposition errors on non-Seurat input", {
  expect_error(CellComposition(list(1, 2), sample_col = "sample"), "Seurat object")
})


# ============================================================================
# CompositionAnalysis()
# ============================================================================

test_that("CompositionAnalysis returns counts/proportions tibbles with expected columns", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 100)
  res <- CompositionAnalysis(obj, group_col = "seurat_clusters", sample_col = "sample")
  expect_true(all(c("sample", "group", "n_cells", "prop") %in% colnames(res$counts)))
  expect_null(res$test)
  by_sample <- stats::aggregate(prop ~ sample, data = res$counts, FUN = sum)
  expect_equal(by_sample$prop, rep(1, nrow(by_sample)), tolerance = 1e-8)
})

test_that("CompositionAnalysis attaches a condition column when condition_col is set", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 100)
  res <- CompositionAnalysis(obj, group_col = "seurat_clusters",
                             sample_col = "sample", condition_col = "condition")
  expect_true("condition" %in% colnames(res$counts))
})

test_that("CompositionAnalysis runs a chisq test when requested (with a pseudoreplication warning)", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 150, n_samples = 6)
  expect_warning(
    res <- CompositionAnalysis(obj, group_col = "seurat_clusters", sample_col = "sample",
                               condition_col = "condition", test = "chisq"),
    "CompositionalTest"
  )
  expect_s3_class(res$test, "htest")
})

test_that("CompositionAnalysis requires condition_col when a test is requested", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 60)
  expect_error(
    CompositionAnalysis(obj, group_col = "seurat_clusters",
                        sample_col = "sample", test = "chisq"),
    "condition_col"
  )
})

test_that("CompositionAnalysis errors on missing columns", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 30)
  expect_error(CompositionAnalysis(obj, group_col = "nope", sample_col = "sample"), "nope")
})


# ============================================================================
# .dabest_from_long() -- internal helper shared with a future
# NicheCoExpress-facing estimation-plot wrapper. Tested directly here
# (not just indirectly through CompositionEstimationPlot()) since it owns
# all the idx/group_levels validation and the dabestr NSE bridging.
# ============================================================================

test_that(".dabest_from_long errors on a missing column", {
  .skip_if_missing("dabestr")
  df <- data.frame(sample = c("S1", "S2"), group = c("g1", "g1"),
                   prop = c(0.1, 0.2), condition = c("A", "B"),
                   stringsAsFactors = FALSE)
  expect_error(.dabest_from_long(df, group_col = "group", y_col = "prop",
                                 condition_col = "nope"),
              "nope")
})

test_that(".dabest_from_long errors when idx doesn't have exactly 2 elements", {
  .skip_if_missing("dabestr")
  df <- data.frame(sample = paste0("S", 1:4), group = "g1",
                   prop = stats::runif(4), condition = c("A", "A", "B", "B"),
                   stringsAsFactors = FALSE)
  expect_error(.dabest_from_long(df, "group", "prop", "condition", idx = "A"),
              "exactly 2")
})

test_that(".dabest_from_long errors when idx values aren't present in condition_col", {
  .skip_if_missing("dabestr")
  df <- data.frame(sample = paste0("S", 1:4), group = "g1",
                   prop = stats::runif(4), condition = c("A", "A", "B", "B"),
                   stringsAsFactors = FALSE)
  expect_error(.dabest_from_long(df, "group", "prop", "condition", idx = c("A", "Z")),
              "Z")
})

test_that(".dabest_from_long errors when group_levels aren't present in group_col", {
  .skip_if_missing("dabestr")
  df <- data.frame(sample = paste0("S", 1:4), group = "g1",
                   prop = stats::runif(4), condition = c("A", "A", "B", "B"),
                   stringsAsFactors = FALSE)
  expect_error(.dabest_from_long(df, "group", "prop", "condition",
                                 idx = c("A", "B"), group_levels = "nope"),
              "nope")
})

test_that(".dabest_from_long errors without idx when condition_col has more than 2 levels", {
  .skip_if_missing("dabestr")
  df <- data.frame(sample = paste0("S", 1:6), group = "g1",
                   prop = stats::runif(6), condition = rep(c("A", "B", "C"), 2),
                   stringsAsFactors = FALSE)
  expect_error(.dabest_from_long(df, "group", "prop", "condition"), "idx")
})

test_that(".dabest_from_long messages and auto-picks alphabetical idx when not supplied", {
  .skip_if_missing("dabestr")
  df <- data.frame(sample = paste0("S", 1:6), group = "g1",
                   prop = stats::runif(6), condition = rep(c("B", "A"), 3),
                   stringsAsFactors = FALSE)
  expect_message(
    out <- .dabest_from_long(df, "group", "prop", "condition"),
    "A.*reference.*B.*test"
  )
  expect_true("g1" %in% names(out))
})

test_that(".dabest_from_long returns dabestr effect-size objects, named by group", {
  .skip_if_missing("dabestr")
  set.seed(1)
  df <- data.frame(
    sample    = rep(paste0("S", 1:6), 2),
    group     = rep(c("g1", "g2"), each = 6),
    prop      = stats::runif(12),
    condition = rep(rep(c("A", "B"), each = 3), 2),
    stringsAsFactors = FALSE
  )
  out <- .dabest_from_long(df, "group", "prop", "condition", idx = c("A", "B"))
  expect_setequal(names(out), c("g1", "g2"))
})

test_that(".dabest_from_long skips (with a warning) a group present in only one condition, keeping the rest", {
  .skip_if_missing("dabestr")
  set.seed(1)
  df <- data.frame(
    sample    = c(paste0("S", 1:6), paste0("S", 7:9)),
    group     = c(rep("g1", 6), rep("g2", 3)),
    prop      = stats::runif(9),
    # g1: both conditions present (3 vs 3). g2: "A" only -- no "B" rows at all.
    condition = c(rep(c("A", "B"), each = 3), rep("A", 3)),
    stringsAsFactors = FALSE
  )
  expect_warning(
    out <- .dabest_from_long(df, "group", "prop", "condition", idx = c("A", "B")),
    "g2.*B"
  )
  expect_setequal(names(out), "g1")
})

test_that(".dabest_from_long errors when every requested group is missing a condition", {
  .skip_if_missing("dabestr")
  # g1 has only condition "A", g2 has only condition "B" -- so both idx
  # levels exist SOMEWHERE in the data (passing the top-level "does idx
  # exist at all" check), but no individual group has both conditions
  # represented, so every group gets skipped and the final "no group had
  # both levels" stop() is what should actually fire.
  df <- data.frame(
    sample    = paste0("S", 1:6),
    group     = c(rep("g1", 3), rep("g2", 3)),
    prop      = stats::runif(6),
    condition = c(rep("A", 3), rep("B", 3)),
    stringsAsFactors = FALSE
  )
  expect_warning(
    expect_error(
      .dabest_from_long(df, "group", "prop", "condition", idx = c("A", "B")),
      "No requested group"
    ),
    "g1.*B"
  )
})


# ============================================================================
# CompositionEstimationPlot()
# ============================================================================

test_that("CompositionEstimationPlot returns a single plot for one group_levels entry", {
  .skip_if_missing("Seurat", "SeuratObject", "dabestr")
  obj <- .make_small_seurat(seed = 1, n_cells = 150, n_samples = 6, n_clusters = 3)
  res <- CompositionAnalysis(obj, group_col = "seurat_clusters", sample_col = "sample",
                             condition_col = "condition")
  one_level <- levels(obj$seurat_clusters)[1]
  p <- CompositionEstimationPlot(res, group_levels = one_level, idx = c("A", "B"))
  expect_false(is.null(p))
  # single group_levels entry -> the plot itself, not a list keyed by group
  expect_false(identical(names(p), one_level))
})

test_that("CompositionEstimationPlot returns a named list for multiple group levels", {
  .skip_if_missing("Seurat", "SeuratObject", "dabestr")
  obj <- .make_small_seurat(seed = 1, n_cells = 150, n_samples = 6, n_clusters = 3)
  res <- CompositionAnalysis(obj, group_col = "seurat_clusters", sample_col = "sample",
                             condition_col = "condition")
  plots <- CompositionEstimationPlot(res, idx = c("A", "B"))
  expect_true(is.list(plots))
  expect_setequal(names(plots), as.character(levels(obj$seurat_clusters)))
})

test_that("CompositionEstimationPlot errors when comp has no condition column", {
  .skip_if_missing("Seurat", "SeuratObject", "dabestr")
  obj <- .make_small_seurat(seed = 1, n_cells = 60)
  res <- CompositionAnalysis(obj, group_col = "seurat_clusters", sample_col = "sample")
  expect_error(CompositionEstimationPlot(res), "condition_col")
})

test_that("CompositionEstimationPlot errors when idx is required but not supplied", {
  .skip_if_missing("dabestr")
  fake_comp <- list(proportions = data.frame(
    sample    = rep(c("S1", "S2", "S3"), each = 2),
    group     = rep(c("g1", "g2"), 3),
    prop      = runif(6),
    condition = rep(c("A", "B", "C"), each = 2),
    stringsAsFactors = FALSE
  ))
  expect_error(CompositionEstimationPlot(fake_comp), "idx")
})

test_that("CompositionEstimationPlot errors when comp isn't a CompositionAnalysis()-shaped list", {
  expect_error(CompositionEstimationPlot(list(counts = data.frame())), "CompositionAnalysis")
  expect_error(CompositionEstimationPlot(data.frame(x = 1)), "CompositionAnalysis")
})

test_that("CompositionEstimationPlot rejects an unknown effect value", {
  # match.arg() fails before `comp`/dabestr are ever touched, so this needs
  # no Seurat/dabestr fixtures at all.
  expect_error(CompositionEstimationPlot(list(proportions = data.frame()), effect = "nope"))
})

test_that("CompositionEstimationPlot passes through a custom effect argument", {
  .skip_if_missing("Seurat", "SeuratObject", "dabestr")
  obj <- .make_small_seurat(seed = 1, n_cells = 150, n_samples = 6, n_clusters = 3)
  res <- CompositionAnalysis(obj, group_col = "seurat_clusters", sample_col = "sample",
                             condition_col = "condition")
  one_level <- levels(obj$seurat_clusters)[1]
  p <- CompositionEstimationPlot(res, group_levels = one_level, idx = c("A", "B"),
                                 effect = "cohens_h")
  expect_false(is.null(p))
})


# ============================================================================
# CompositionBarplot()
# ============================================================================

test_that("CompositionBarplot works directly from a Seurat object", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 60)
  p <- CompositionBarplot(obj, group_col = "seurat_clusters", sample_col = "sample")
  expect_s3_class(p, "ggplot")
})

test_that("CompositionBarplot requires group_col/sample_col for Seurat input", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 30)
  expect_error(CompositionBarplot(obj), "group_col.*sample_col|sample_col.*group_col")
})

test_that("CompositionBarplot accepts a pre-computed proportions tibble", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 60)
  res <- CompositionAnalysis(obj, group_col = "seurat_clusters", sample_col = "sample")
  p <- CompositionBarplot(res$proportions, style = "grouped", y = "n_cells")
  expect_s3_class(p, "ggplot")
})

test_that("CompositionBarplot errors when a manual data frame is missing required columns", {
  expect_error(
    CompositionBarplot(data.frame(sample = "a", group = "b")),
    "missing columns"
  )
})


# ============================================================================
# CompositionalTest()
# ============================================================================

test_that("CompositionalTest (wilcox backend) returns one row per cluster with padj", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 200, n_samples = 6, n_clusters = 3)
  res <- CompositionalTest(obj, cluster_col = "seurat_clusters",
                           sample_col = "sample", condition_col = "condition",
                           method = "wilcox")
  expect_true(all(c("cluster", "effect", "stat", "pvalue", "padj", "method")
                 %in% colnames(res)))
  expect_equal(unique(res$method), "wilcox")
  expect_equal(nrow(res), length(unique(obj$seurat_clusters)))
})

test_that("CompositionalTest requires condition_col", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 60)
  expect_error(CompositionalTest(obj, sample_col = "sample"), "condition_col")
})

test_that("CompositionalTest errors when a sample has inconsistent condition values", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 60, n_samples = 2)
  md <- obj@meta.data
  s1_cells <- rownames(md)[md$sample == "S1"]
  obj$condition[s1_cells[1]] <- setdiff(unique(obj$condition), obj$condition[s1_cells[1]])[1]
  expect_error(
    CompositionalTest(obj, sample_col = "sample", condition_col = "condition"),
    "inconsistent condition"
  )
})

test_that("CompositionalTest (propeller backend) runs when speckle is installed", {
  .skip_if_missing("Seurat", "SeuratObject", "speckle")
  obj <- .make_small_seurat(seed = 1, n_cells = 200, n_samples = 6, n_clusters = 3)
  res <- CompositionalTest(obj, sample_col = "sample", condition_col = "condition",
                           method = "propeller")
  expect_equal(unique(res$method), "propeller")
  expect_true("padj" %in% colnames(res))
})

test_that("CompositionalTest (betareg backend) runs when betareg is installed", {
  .skip_if_missing("Seurat", "SeuratObject", "betareg")
  obj <- .make_small_seurat(seed = 1, n_cells = 200, n_samples = 6, n_clusters = 3)
  res <- CompositionalTest(obj, sample_col = "sample", condition_col = "condition",
                           method = "betareg")
  expect_equal(unique(res$method), "betareg")
})


# ============================================================================
# CompositionalTest() -- weight_cols (continuous) mode, e.g. RCTD "full"-mode
# rctd_<celltype> proportions. Synthetic numeric metadata columns stand in
# for real RCTD output -- weight_cols mode only cares that the columns are
# numeric, not where they came from.
# ============================================================================

test_that("CompositionalTest weight_cols mode (wilcox backend) returns one row per weight column", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 200, n_samples = 6, n_clusters = 3)
  set.seed(1)
  obj$w_typeA <- stats::runif(ncol(obj))
  obj$w_typeB <- stats::runif(ncol(obj))
  res <- CompositionalTest(obj, weight_cols = c("w_typeA", "w_typeB"),
                           sample_col = "sample", condition_col = "condition",
                           method = "wilcox")
  expect_true(all(c("cluster", "effect", "stat", "pvalue", "padj", "method")
                 %in% colnames(res)))
  expect_equal(unique(res$method), "wilcox")
  expect_setequal(res$cluster, c("w_typeA", "w_typeB"))
})

test_that("CompositionalTest weight_cols mode ignores cluster_col", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 100, n_samples = 4, n_clusters = 3)
  obj$w_typeA <- stats::runif(ncol(obj))
  res <- CompositionalTest(obj, cluster_col = "nonexistent_column",
                           weight_cols = "w_typeA",
                           sample_col = "sample", condition_col = "condition",
                           method = "wilcox")
  expect_equal(res$cluster, "w_typeA")
})

test_that("CompositionalTest weight_cols + method = 'propeller' errors clearly", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 60, n_samples = 4)
  obj$w_typeA <- stats::runif(ncol(obj))
  expect_error(
    CompositionalTest(obj, weight_cols = "w_typeA", sample_col = "sample",
                      condition_col = "condition", method = "propeller"),
    "propeller"
  )
})

test_that("CompositionalTest weight_cols errors on a non-numeric column", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 60, n_samples = 4)
  expect_error(
    CompositionalTest(obj, weight_cols = "seurat_clusters", sample_col = "sample",
                      condition_col = "condition", method = "wilcox"),
    "numeric"
  )
})

test_that("CompositionalTest weight_cols errors on a missing column", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 60, n_samples = 4)
  expect_error(
    CompositionalTest(obj, weight_cols = "nope", sample_col = "sample",
                      condition_col = "condition", method = "wilcox"),
    "nope"
  )
})

test_that("CompositionalTest weight_cols mode (betareg backend) runs when betareg is installed", {
  .skip_if_missing("Seurat", "SeuratObject", "betareg")
  obj <- .make_small_seurat(seed = 1, n_cells = 200, n_samples = 6, n_clusters = 3)
  obj$w_typeA <- stats::runif(ncol(obj), min = 0.1, max = 0.9)
  obj$w_typeB <- stats::runif(ncol(obj), min = 0.1, max = 0.9)
  res <- CompositionalTest(obj, weight_cols = c("w_typeA", "w_typeB"),
                           sample_col = "sample", condition_col = "condition",
                           method = "betareg")
  expect_equal(unique(res$method), "betareg")
  expect_setequal(res$cluster, c("w_typeA", "w_typeB"))
})
