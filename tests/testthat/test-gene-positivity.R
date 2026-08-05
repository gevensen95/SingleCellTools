# Tests for GenePositivityAnalysis() and GenePositivityEstimationPlot() --
# the gene-positivity counterparts to CompositionAnalysis() /
# CompositionEstimationPlot() (see test-composition.R). Both reuse the
# .dabest_from_long() internal helper, whose own validation logic is already
# covered directly in test-composition.R -- not duplicated here.

# ============================================================================
# GenePositivityAnalysis()
# ============================================================================

test_that("GenePositivityAnalysis returns counts/proportions tibbles with expected columns", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 100)
  obj <- AddGenePositivity(obj, genes = c("Gene1", "Gene2"))
  res <- GenePositivityAnalysis(obj, genes = c("Gene1", "Gene2"), sample_col = "sample")
  expect_true(all(c("sample", "gene", "n_pos", "n_total", "prop_pos")
                 %in% colnames(res$counts)))
  expect_false("group" %in% colnames(res$counts))
  expect_null(res$test)
  expect_true(all(res$counts$prop_pos >= 0 & res$counts$prop_pos <= 1))
  expect_identical(res$counts, res$proportions)
})

test_that("GenePositivityAnalysis errors when a positivity column is missing", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 60)
  obj <- AddGenePositivity(obj, genes = "Gene1")
  expect_error(
    GenePositivityAnalysis(obj, genes = c("Gene1", "Gene2"), sample_col = "sample"),
    "Gene2_pos"
  )
})

test_that("GenePositivityAnalysis attaches a condition column when condition_col is set", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 100)
  obj <- AddGenePositivity(obj, genes = "Gene1")
  res <- GenePositivityAnalysis(obj, genes = "Gene1", sample_col = "sample",
                                condition_col = "condition")
  expect_true("condition" %in% colnames(res$counts))
})

test_that("GenePositivityAnalysis stratifies by group_col when supplied", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 100)
  obj <- AddGenePositivity(obj, genes = "Gene1")
  res <- GenePositivityAnalysis(obj, genes = "Gene1", sample_col = "sample",
                                group_col = "celltype")
  expect_true("group" %in% colnames(res$counts))
  expect_setequal(res$counts$group, unique(obj$celltype))
})

test_that("GenePositivityAnalysis runs a chisq test per gene when requested (with a pseudoreplication warning)", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 150, n_samples = 6)
  obj <- AddGenePositivity(obj, genes = c("Gene1", "Gene2"))
  expect_warning(
    res <- GenePositivityAnalysis(obj, genes = c("Gene1", "Gene2"), sample_col = "sample",
                                  condition_col = "condition", test = "chisq"),
    "pseudoreplication"
  )
  expect_setequal(names(res$test), c("Gene1", "Gene2"))
  expect_s3_class(res$test[["Gene1"]], "htest")
})

test_that("GenePositivityAnalysis names tests '<gene> | <group>' when group_col is used", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 150, n_samples = 6)
  obj <- AddGenePositivity(obj, genes = "Gene1")
  expect_warning(
    res <- GenePositivityAnalysis(obj, genes = "Gene1", sample_col = "sample",
                                  condition_col = "condition", group_col = "celltype",
                                  test = "fisher"),
    "pseudoreplication"
  )
  expected_names <- paste("Gene1", unique(as.character(obj$celltype)), sep = " | ")
  expect_setequal(names(res$test), expected_names)
})

test_that("GenePositivityAnalysis requires condition_col when a test is requested", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 60)
  obj <- AddGenePositivity(obj, genes = "Gene1")
  expect_error(
    GenePositivityAnalysis(obj, genes = "Gene1", sample_col = "sample", test = "chisq"),
    "condition_col"
  )
})

test_that("GenePositivityAnalysis errors on missing sample_col/group_col columns", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 30)
  obj <- AddGenePositivity(obj, genes = "Gene1")
  expect_error(
    GenePositivityAnalysis(obj, genes = "Gene1", sample_col = "nope"), "nope"
  )
  expect_error(
    GenePositivityAnalysis(obj, genes = "Gene1", sample_col = "sample", group_col = "nope"),
    "nope"
  )
})

test_that("GenePositivityAnalysis errors on non-Seurat input", {
  expect_error(GenePositivityAnalysis(list(1, 2), genes = "Gene1", sample_col = "sample"),
              "Seurat object")
})


# ============================================================================
# GenePositivityEstimationPlot()
# ============================================================================

test_that("GenePositivityEstimationPlot returns a single plot for one gene", {
  .skip_if_missing("Seurat", "SeuratObject", "dabestr")
  obj <- .make_small_seurat(seed = 1, n_cells = 150, n_samples = 6)
  obj <- AddGenePositivity(obj, genes = c("Gene1", "Gene2"))
  res <- GenePositivityAnalysis(obj, genes = c("Gene1", "Gene2"), sample_col = "sample",
                                condition_col = "condition")
  p <- GenePositivityEstimationPlot(res, genes = "Gene1", idx = c("A", "B"))
  expect_false(is.null(p))
  expect_false(identical(names(p), "Gene1"))
})

test_that("GenePositivityEstimationPlot returns a named list for multiple genes", {
  .skip_if_missing("Seurat", "SeuratObject", "dabestr")
  obj <- .make_small_seurat(seed = 1, n_cells = 150, n_samples = 6)
  obj <- AddGenePositivity(obj, genes = c("Gene1", "Gene2"))
  res <- GenePositivityAnalysis(obj, genes = c("Gene1", "Gene2"), sample_col = "sample",
                                condition_col = "condition")
  plots <- GenePositivityEstimationPlot(res, idx = c("A", "B"))
  expect_true(is.list(plots))
  expect_setequal(names(plots), c("Gene1", "Gene2"))
})

test_that("GenePositivityEstimationPlot uses '<gene> | <group>' naming when group_col was used", {
  .skip_if_missing("Seurat", "SeuratObject", "dabestr")
  obj <- .make_small_seurat(seed = 1, n_cells = 150, n_samples = 6)
  obj <- AddGenePositivity(obj, genes = "Gene1")
  res <- GenePositivityAnalysis(obj, genes = "Gene1", sample_col = "sample",
                                condition_col = "condition", group_col = "celltype")
  plots <- GenePositivityEstimationPlot(res, idx = c("A", "B"))
  expected_names <- paste("Gene1", unique(as.character(obj$celltype)), sep = " | ")
  expect_setequal(names(plots), expected_names)
})

test_that("GenePositivityEstimationPlot errors when gpa has no condition column", {
  .skip_if_missing("Seurat", "SeuratObject", "dabestr")
  obj <- .make_small_seurat(seed = 1, n_cells = 60)
  obj <- AddGenePositivity(obj, genes = "Gene1")
  res <- GenePositivityAnalysis(obj, genes = "Gene1", sample_col = "sample")
  expect_error(GenePositivityEstimationPlot(res), "condition_col")
})

test_that("GenePositivityEstimationPlot errors when group_levels is supplied without group_col", {
  .skip_if_missing("Seurat", "SeuratObject", "dabestr")
  obj <- .make_small_seurat(seed = 1, n_cells = 60)
  obj <- AddGenePositivity(obj, genes = "Gene1")
  res <- GenePositivityAnalysis(obj, genes = "Gene1", sample_col = "sample",
                                condition_col = "condition")
  expect_error(GenePositivityEstimationPlot(res, group_levels = "Type1"), "group_col")
})

test_that("GenePositivityEstimationPlot errors when gpa isn't a GenePositivityAnalysis()-shaped list", {
  expect_error(GenePositivityEstimationPlot(list(counts = data.frame())),
              "GenePositivityAnalysis")
  expect_error(GenePositivityEstimationPlot(data.frame(x = 1)), "GenePositivityAnalysis")
})

test_that("GenePositivityEstimationPlot rejects an unknown effect value", {
  expect_error(GenePositivityEstimationPlot(list(proportions = data.frame()), effect = "nope"))
})

test_that("GenePositivityEstimationPlot passes through a custom effect argument", {
  .skip_if_missing("Seurat", "SeuratObject", "dabestr")
  obj <- .make_small_seurat(seed = 1, n_cells = 150, n_samples = 6)
  obj <- AddGenePositivity(obj, genes = "Gene1")
  res <- GenePositivityAnalysis(obj, genes = "Gene1", sample_col = "sample",
                                condition_col = "condition")
  p <- GenePositivityEstimationPlot(res, genes = "Gene1", idx = c("A", "B"),
                                    effect = "cohens_h")
  expect_false(is.null(p))
})
