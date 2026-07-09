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

test_that("CompositionAnalysis runs a chisq test when requested", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 150, n_samples = 6)
  res <- CompositionAnalysis(obj, group_col = "seurat_clusters", sample_col = "sample",
                             condition_col = "condition", test = "chisq")
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
