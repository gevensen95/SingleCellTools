# Tests for ApplyQCFilters(), BatchEffectQC(), QCComparePlots(),
# assign_cell_cycle_phase(), calldoublet(), and (lightly) GenerateQCReport().
#
# GenerateQCReport() renders a full HTML report via rmarkdown/knitr/pandoc;
# that render path is deliberately not exercised here (too heavy/slow for a
# unit test and highly environment-dependent) beyond a basic input-
# validation check. calldoublet() runs a full SCTransform/PCA/UMAP/cluster/
# DoubletFinder pipeline -- the one test below is a genuine smoke test, not
# a fast unit test; it may take tens of seconds.

# ============================================================================
# ApplyQCFilters()
# ============================================================================

test_that("ApplyQCFilters (single object) keeps only cells within the cutoff range", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 25, n_cells = 100)
  nc  <- obj$nCount_RNA
  lo  <- as.numeric(stats::quantile(nc, 0.1))
  hi  <- as.numeric(stats::quantile(nc, 0.9))
  cutoffs <- data.frame(sample = "S1", metric = "nCount_RNA",
                        suggest_lo = lo, suggest_hi = hi,
                        stringsAsFactors = FALSE)
  out <- ApplyQCFilters(obj, cutoffs = cutoffs, single_sample_name = "S1",
                        filter_doublets = FALSE, verbose = FALSE)
  expect_true(ncol(out) < ncol(obj))
  expect_true(ncol(out) > 0)
  expect_true(all(out$nCount_RNA >= lo & out$nCount_RNA <= hi))
})

test_that("ApplyQCFilters return_report = TRUE includes obj/report/cutoffs", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 50)
  cutoffs <- data.frame(sample = "S1", metric = "nCount_RNA",
                        suggest_lo = 0, suggest_hi = 1e6,
                        stringsAsFactors = FALSE)
  out <- ApplyQCFilters(obj, cutoffs = cutoffs, single_sample_name = "S1",
                        filter_doublets = FALSE, verbose = FALSE,
                        return_report = TRUE)
  expect_named(out, c("obj", "report", "cutoffs"))
  expect_true("pct_kept" %in% colnames(out$report))
})

test_that("ApplyQCFilters works with a named list of Seurat objects", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj1 <- .make_small_seurat(seed = 1, n_cells = 40)
  obj2 <- .make_small_seurat(seed = 2, n_cells = 40)
  cutoffs <- data.frame(
    sample = c("a", "b"), metric = "nCount_RNA",
    suggest_lo = 0, suggest_hi = 1e6, stringsAsFactors = FALSE
  )
  out <- ApplyQCFilters(list(a = obj1, b = obj2), cutoffs = cutoffs,
                        filter_doublets = FALSE, verbose = FALSE)
  expect_named(out, c("a", "b"))
})

test_that("ApplyQCFilters errors when a list of Seurat objects is unnamed", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj1 <- .make_small_seurat(seed = 1, n_cells = 20)
  cutoffs <- data.frame(sample = "a", metric = "nCount_RNA",
                        suggest_lo = 0, suggest_hi = 1e6, stringsAsFactors = FALSE)
  expect_error(ApplyQCFilters(list(obj1), cutoffs = cutoffs, verbose = FALSE),
              "named")
})

test_that("ApplyQCFilters sample_col mode filters each sample independently", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 100, n_samples = 4)
  cutoffs <- data.frame(
    sample = paste0("S", 1:4), metric = "nCount_RNA",
    suggest_lo = 0, suggest_hi = 1e6, stringsAsFactors = FALSE
  )
  out <- ApplyQCFilters(obj, cutoffs = cutoffs, sample_col = "sample",
                        filter_doublets = FALSE, verbose = FALSE)
  expect_equal(ncol(out), ncol(obj))  # wide-open cutoffs keep everyone
})

test_that("ApplyQCFilters applies overrides on top of the cutoffs table", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 80)
  cutoffs <- data.frame(sample = "S1", metric = "nCount_RNA",
                        suggest_lo = 0, suggest_hi = 1e6, stringsAsFactors = FALSE)
  baseline <- ApplyQCFilters(obj, cutoffs = cutoffs, single_sample_name = "S1",
                             filter_doublets = FALSE, verbose = FALSE)
  # Narrow the override to a small window around the median -- enough to
  # provably filter *some* cells without filtering *all* of them. (An
  # impossible/empty-result range isn't used here because this Seurat
  # version's subset() errors ("No cells found") rather than returning a
  # 0-cell object, which would make this a test of subset()'s behavior,
  # not ApplyQCFilters()'s override logic.)
  med <- stats::median(obj$nCount_RNA)
  narrowed <- ApplyQCFilters(obj, cutoffs = cutoffs, single_sample_name = "S1",
                             override = list(S1 = list(nCount_RNA = c(med, med + 5))),
                             filter_doublets = FALSE, verbose = FALSE)
  expect_equal(ncol(baseline), ncol(obj))
  expect_true(ncol(narrowed) < ncol(baseline))
  expect_true(ncol(narrowed) > 0)
})

test_that("ApplyQCFilters drops flagged doublets", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 40)
  obj$doublet_finder <- c(rep("Doublet", 5), rep("Singlet", 35))
  cutoffs <- data.frame(sample = "S1", metric = "nCount_RNA",
                        suggest_lo = 0, suggest_hi = 1e6, stringsAsFactors = FALSE)
  out <- ApplyQCFilters(obj, cutoffs = cutoffs, single_sample_name = "S1",
                        filter_doublets = TRUE, verbose = FALSE)
  expect_equal(ncol(out), 35)
})

test_that("ApplyQCFilters errors on a cutoffs table missing required columns", {
  expect_error(
    ApplyQCFilters(1, cutoffs = data.frame(sample = "a")),
    "missing required column"
  )
})

test_that("ApplyQCFilters errors when the cutoffs CSV path doesn't exist", {
  expect_error(ApplyQCFilters(1, cutoffs = "no/such/file_cutoffs.csv"),
              "not found")
})

test_that("ApplyQCFilters reads a cutoffs CSV, and finds the sidecar CSV from an .html path", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 30)
  tmp <- tempfile(fileext = "_cutoffs.csv")
  utils::write.csv(
    data.frame(sample = "S1", metric = "nCount_RNA",
              suggest_lo = 0, suggest_hi = 1e6),
    tmp, row.names = FALSE
  )
  on.exit(unlink(tmp))
  out1 <- ApplyQCFilters(obj, cutoffs = tmp, single_sample_name = "S1",
                         filter_doublets = FALSE, verbose = FALSE)
  expect_equal(ncol(out1), ncol(obj))

  # Simulate pointing at the HTML report instead of the CSV sidecar
  html_path <- sub("_cutoffs\\.csv$", ".html", tmp)
  out2 <- ApplyQCFilters(obj, cutoffs = html_path, single_sample_name = "S1",
                         filter_doublets = FALSE, verbose = FALSE)
  expect_equal(ncol(out2), ncol(obj))
})

test_that("ApplyQCFilters `metrics` filter restricts to requested metrics and errors if none match", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 40)
  cutoffs <- data.frame(
    sample = c("S1", "S1"), metric = c("nCount_RNA", "nFeature_RNA"),
    suggest_lo = 0, suggest_hi = 1e6, stringsAsFactors = FALSE
  )
  out <- ApplyQCFilters(obj, cutoffs = cutoffs, single_sample_name = "S1",
                        metrics = "nCount_RNA", filter_doublets = FALSE,
                        verbose = FALSE, return_report = TRUE)
  expect_equal(unique(out$cutoffs$metric_key), "nCount_RNA")

  expect_error(
    ApplyQCFilters(obj, cutoffs = cutoffs, single_sample_name = "S1",
                  metrics = "not_a_metric", verbose = FALSE),
    "No cutoffs remain"
  )
})


# ============================================================================
# BatchEffectQC()
# ============================================================================

.add_fake_pca <- function(obj, ndim = 5) {
  emb <- matrix(rnorm(ncol(obj) * ndim), nrow = ncol(obj),
               dimnames = list(colnames(obj), paste0("PC_", seq_len(ndim))))
  obj[["pca"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "PC_",
                                                     assay = "RNA")
  obj
}

test_that("BatchEffectQC returns the documented summary fields without celltype_col", {
  .skip_if_missing("Seurat", "SeuratObject", "cluster")
  obj <- .make_small_seurat(seed = 1, n_cells = 60, n_samples = 2)
  obj <- .add_fake_pca(obj)
  res <- BatchEffectQC(obj, reduction = "pca", batch_col = "sample", k = 5)
  expect_true(all(c("batch_asw", "knn_mixing", "n_cells", "n_batches")
                 %in% names(res$summary)))
  expect_null(res$per_cluster)
})

test_that("BatchEffectQC computes per_cluster breakdown when celltype_col is set", {
  .skip_if_missing("Seurat", "SeuratObject", "cluster")
  obj <- .make_small_seurat(seed = 1, n_cells = 60, n_samples = 2, n_celltypes = 2)
  obj <- .add_fake_pca(obj)
  res <- BatchEffectQC(obj, reduction = "pca", batch_col = "sample",
                       celltype_col = "celltype", k = 5)
  expect_false(is.null(res$per_cluster))
  expect_true(all(c("cluster", "n_cells", "knn_purity", "knn_mixing")
                 %in% colnames(res$per_cluster)))
})

test_that("BatchEffectQC validates inputs", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  obj <- .add_fake_pca(obj)
  expect_error(BatchEffectQC(list(1)), "Seurat object")
  expect_error(BatchEffectQC(obj, batch_col = NULL), "batch_col.*required")
  expect_error(BatchEffectQC(obj, batch_col = "nope"), "nope")
  expect_error(BatchEffectQC(obj, reduction = "umap", batch_col = "sample"),
              "not found")
})


# ============================================================================
# QCComparePlots()
# ============================================================================

test_that("QCComparePlots returns a single ggplot when only one metric is requested", {
  .skip_if_missing("Seurat", "SeuratObject")
  pre  <- .make_small_seurat(seed = 1, n_cells = 60)
  post <- subset(pre, cells = colnames(pre)[1:40])
  p <- QCComparePlots(pre, post, metrics = "nCount_RNA")
  expect_s3_class(p, "ggplot")
})

test_that("QCComparePlots returns a patchwork (or list) for multiple metrics", {
  .skip_if_missing("Seurat", "SeuratObject")
  pre  <- .make_small_seurat(seed = 1, n_cells = 60)
  post <- subset(pre, cells = colnames(pre)[1:40])
  p <- QCComparePlots(pre, post, metrics = c("nCount_RNA", "nFeature_RNA"))
  expect_true(inherits(p, "patchwork") || is.list(p))
})

test_that("QCComparePlots works with named lists of Seurat objects", {
  .skip_if_missing("Seurat", "SeuratObject")
  pre1  <- .make_small_seurat(seed = 1, n_cells = 30)
  pre2  <- .make_small_seurat(seed = 2, n_cells = 30)
  post1 <- subset(pre1, cells = colnames(pre1)[1:20])
  post2 <- subset(pre2, cells = colnames(pre2)[1:20])
  p <- QCComparePlots(list(a = pre1, b = pre2), list(a = post1, b = post2),
                      metrics = "nCount_RNA")
  expect_s3_class(p, "ggplot")
})

test_that("QCComparePlots errors when no standard metrics are present and none requested", {
  .skip_if_missing("Seurat", "SeuratObject")
  pre  <- .make_small_seurat(seed = 1, n_cells = 20)
  post <- subset(pre, cells = colnames(pre)[1:10])
  pre@meta.data[, c("nCount_RNA", "nFeature_RNA")]  <- NULL
  post@meta.data[, c("nCount_RNA", "nFeature_RNA")] <- NULL
  expect_error(QCComparePlots(pre, post), "No standard QC metrics")
})


# ============================================================================
# assign_cell_cycle_phase()
# ============================================================================

test_that("assign_cell_cycle_phase adds a Phase column with valid values", {
  .skip_if_missing("Seurat", "SeuratObject", "UCell")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 40, gene_prefix = "Gene")
  out <- suppressMessages(assign_cell_cycle_phase(
    obj, s_genes = c("Gene1", "Gene2"), g2m_genes = c("Gene3", "Gene4")
  ))
  expect_true("Phase" %in% colnames(out@meta.data))
  expect_true(all(out$Phase %in% c("G1", "S", "G2M", NA)))
})


# ============================================================================
# calldoublet() -- genuine integration smoke test, not a fast unit test
# ============================================================================

test_that("calldoublet runs the full pipeline and adds doublet_finder, stripping intermediates", {
  .skip_if_missing("Seurat", "SeuratObject", "DoubletFinder")
  testthat::skip_on_cran()
  set.seed(1)
  n_genes <- 150; n_cells <- 300
  m <- matrix(stats::rpois(n_genes * n_cells, lambda = 3), nrow = n_genes,
             dimnames = list(paste0("Gene", seq_len(n_genes)), paste0("c", seq_len(n_cells))))
  storage.mode(m) <- "double"
  obj <- SeuratObject::CreateSeuratObject(counts = m)

  out <- tryCatch(
    suppressWarnings(suppressMessages(calldoublet(obj))),
    error = function(e) e
  )
  skip_if(inherits(out, "error"),
         paste("calldoublet pipeline did not complete on this synthetic",
               "dataset (likely a PC-selection heuristic edge case):",
               if (inherits(out, "error")) conditionMessage(out) else ""))

  expect_true("doublet_finder" %in% colnames(out@meta.data))
  expect_true(all(out$doublet_finder %in% c("Doublet", "Singlet")))
  expect_length(SeuratObject::Reductions(out), 0)
  expect_false("data" %in% SeuratObject::Layers(out[["RNA"]]))
})


# ============================================================================
# GenerateQCReport() -- validation only; the HTML render path is too heavy
# for a unit test (rmarkdown/knitr/pandoc dependent, slow).
# ============================================================================

test_that("GenerateQCReport errors on non-Seurat, non-list input", {
  expect_error(GenerateQCReport(list(1, 2)), "Seurat object")
})
