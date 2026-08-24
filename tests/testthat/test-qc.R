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

test_that("ApplyQCFilters skips a metric whose suggest_lo == suggest_hi", {
  # Regression: a degenerate [lo, hi] kept only cells whose value was
  # exactly that number. It shows up when a metric is constant across a
  # sample -- so the samples needing filtering least got gutted.
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 60)
  obj$flat_metric <- 5                      # constant -> lo == hi
  cutoffs <- data.frame(
    sample     = c("S1", "S1"),
    metric     = c("flat_metric", "nCount_RNA"),
    suggest_lo = c(5, 0),
    suggest_hi = c(5, 1e6),
    stringsAsFactors = FALSE
  )
  out <- ApplyQCFilters(obj, cutoffs = cutoffs, single_sample_name = "S1",
                        filter_doublets = FALSE, verbose = FALSE,
                        return_report = TRUE)
  # Degenerate metric contributed nothing; the wide-open one kept everything.
  expect_equal(ncol(out$obj), ncol(obj))
  expect_false("flat_metric" %in% out$report$metric)
  expect_true("nCount_RNA" %in% out$report$metric)
})

test_that("ApplyQCFilters skips metrics with missing or inverted cutoffs", {
  # An NA bound used to poison the mask: `v >= NA` is NA, so `keep` carried
  # NAs, sum(keep) was NA, and cells[keep] produced NA cell names.
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 60)

  na_cut <- data.frame(sample = "S1", metric = "nCount_RNA",
                       suggest_lo = NA_real_, suggest_hi = 1e6,
                       stringsAsFactors = FALSE)
  out_na <- ApplyQCFilters(obj, cutoffs = na_cut, single_sample_name = "S1",
                           filter_doublets = FALSE, verbose = FALSE)
  expect_equal(ncol(out_na), ncol(obj))

  # Inverted bounds would keep zero cells; skip, but warn -- the table is wrong.
  inv_cut <- data.frame(sample = "S1", metric = "nCount_RNA",
                        suggest_lo = 1e6, suggest_hi = 0,
                        stringsAsFactors = FALSE)
  expect_warning(
    out_inv <- ApplyQCFilters(obj, cutoffs = inv_cut, single_sample_name = "S1",
                              filter_doublets = FALSE, verbose = FALSE),
    "inverted"
  )
  expect_equal(ncol(out_inv), ncol(obj))
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

test_that("QCComparePlots tolerates per-sample metadata column names (e.g. pANN_*)", {
  # Regression: DoubletFinder names its pANN column
  # "pANN_<pN>_<pK>_<nExp>", and pK/nExp are fit per sample -- so a list of
  # objects has matching column COUNTS but one differently-named column
  # each. That made rbind.data.frame() fail with "names do not match
  # previous names", which named neither the column nor the sample.
  .skip_if_missing("Seurat", "SeuratObject")
  mk <- function(seed, pk, nexp) {
    o <- .make_small_seurat(seed = seed, n_cells = 30)
    o@meta.data[[sprintf("pANN_0.25_%s_%d", pk, nexp)]] <-
      stats::runif(ncol(o))
    o
  }
  pre1 <- mk(1, 0.005, 22);  pre2 <- mk(2, 0.24, 1069)
  post1 <- subset(pre1, cells = colnames(pre1)[1:20])
  post2 <- subset(pre2, cells = colnames(pre2)[1:20])

  expect_warning(
    p <- QCComparePlots(list(a = pre1, b = pre2), list(a = post1, b = post2),
                        metrics = "nCount_RNA"),
    "not present in every"
  )
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
  # Three cell populations, each with its own 30-gene marker module, plus
  # per-cell library size factors. Both parts are load-bearing:
  #
  #  * Structure (vs flat Poisson noise) is needed at all, or the PC stdev
  #    curve is flat, calldoublet()'s PC-selection heuristics find nothing,
  #    and the whole test skips -- which is how it behaved originally,
  #    meaning none of the assertions below had ever actually run.
  #
  #  * But the separation must stay MODERATE (4x enrichment over background,
  #    not 30x) and cells must vary in depth. With sharply separated,
  #    equal-depth populations, paramSweep()'s large-pK grid points give
  #    every real cell an identical fraction of artificial-doublet
  #    neighbours -- pANN becomes constant, its kernel density estimate has
  #    zero bandwidth, and summarizeSweep() dies inside KernSmooth::bkde()
  #    with "missing value where TRUE/FALSE needed". These values were
  #    checked against all 159 (pN, pK) grid points paramSweep visits here;
  #    a 30x-enrichment version degenerated at 8 of them.
  n_genes <- 200; n_cells <- 300; n_grp <- 3
  grp    <- rep(seq_len(n_grp), length.out = n_cells)
  size   <- stats::rlnorm(n_cells, 0, 0.35)
  lambda <- matrix(1, nrow = n_genes, ncol = n_cells)
  for (g in seq_len(n_grp)) {
    lambda[((g - 1) * 30 + 1):(g * 30), grp == g] <- 4
  }
  lambda <- sweep(lambda, 2, size, "*")
  m <- matrix(stats::rpois(n_genes * n_cells, lambda), nrow = n_genes,
              dimnames = list(paste0("Gene", seq_len(n_genes)),
                              paste0("c", seq_len(n_cells))))
  storage.mode(m) <- "double"
  obj <- SeuratObject::CreateSeuratObject(counts = methods::as(m, "CsparseMatrix"))

  out <- tryCatch(
    suppressWarnings(suppressMessages(calldoublet(obj))),
    error = function(e) e
  )
  # Skip ONLY on the PC-selection edge case, and fail on anything else. A
  # blanket skip_if(inherits(out, "error")) turns every real regression in
  # this function into a silent green run.
  if (inherits(out, "error")) {
    skip_if(grepl("significant PCs", conditionMessage(out)),
            paste("PC-selection heuristic found nothing on this fixture:",
                  conditionMessage(out)))
    stop(out)
  }

  expect_true("doublet_finder" %in% colnames(out@meta.data))
  expect_true(all(out$doublet_finder %in% c("Doublet", "Singlet")))
  expect_length(SeuratObject::Reductions(out), 0)

  # Assert on the WHOLE layer set, not just absence of "data". The old
  # `obj[["RNA"]][["data"]] <- NULL` deletion was a silent no-op on Assay5,
  # and a narrow expect_false() would still miss a stray scale.data or a v5
  # split variant (data.1, scale.data.2, ...).
  expect_equal(SeuratObject::Layers(out[["RNA"]]), "counts")

  # Both DoubletFinder columns must be renamed to stable names. Any surviving
  # run-parameterized name (pANN_<pN>_<pK>_<nExp>, DF.classifications_...)
  # differs per sample and breaks every downstream function that stacks
  # metadata across a sample list.
  expect_true("doublet_pANN" %in% colnames(out@meta.data))
  expect_length(grep("^pANN_", colnames(out@meta.data)), 0)
  expect_length(grep("^DF\\.classifications", colnames(out@meta.data)), 0)

  # nExp must come from doublet_rate, not from the fitted pK. With the
  # default 0.075 and homotypic adjustment, the called rate should land in a
  # plausible band -- not the 0.4%-18% spread that keying nExp to pK produced.
  expect_lt(mean(out$doublet_finder == "Doublet"), 0.15)
})

test_that("calldoublet rejects an implausible doublet_rate before doing any work", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  expect_error(calldoublet(obj, doublet_rate = 1.5), "must be a single number")
  expect_error(calldoublet(obj, doublet_rate = -0.1), "must be a single number")
  expect_error(calldoublet(obj, doublet_rate = c(0.05, 0.1)), "must be a single number")
  expect_error(calldoublet(obj, doublet_rate = NA), "must be a single number")
})


# ============================================================================
# GenerateQCReport() -- validation only; the HTML render path is too heavy
# for a unit test (rmarkdown/knitr/pandoc dependent, slow).
# ============================================================================

test_that("GenerateQCReport errors on non-Seurat, non-list input", {
  expect_error(GenerateQCReport(list(1, 2)), "Seurat object")
})
