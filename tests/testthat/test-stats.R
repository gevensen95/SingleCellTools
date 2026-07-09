# Tests for CompareMarkers() (pure data.frame, no Seurat needed),
# CellSuiteSummary(), and PseudobulkDE().

# ============================================================================
# CompareMarkers() -- pure data.frame comparison, no Seurat object needed
# ============================================================================

.make_de_pair <- function() {
  # 8 genes shared between a/b (deliberately covering every category),
  # plus 2 genes only in `a` and 1 only in `b` to exercise the inner-join
  # behavior. `a` and `b` intentionally use different column-naming
  # conventions to exercise .detect_col()'s auto-detection.
  a <- data.frame(
    gene       = c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10"),
    avg_log2FC = c(2, -2, 2, 2, 0.1, 0.1, -0.2, 0.05, 3, -3),
    p_val_adj  = c(0.001, 0.001, 0.001, 0.001, 0.9, 0.8, 0.7, 0.95, 0.001, 0.001),
    stringsAsFactors = FALSE
  )
  b <- data.frame(
    Gene  = c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G11"),
    logFC = c(1.5, -1.8, -2, 0.1, 2, 0.2, -0.1, 0.05, 4),
    FDR   = c(0.001, 0.01, 0.001, 0.9, 0.001, 0.85, 0.7, 0.9, 0.001),
    stringsAsFactors = FALSE
  )
  list(a = a, b = b)
}

test_that("CompareMarkers auto-detects column names and merges on shared genes only", {
  de <- .make_de_pair()
  out <- CompareMarkers(de$a, de$b, labels = c("a", "b"), plot = FALSE)
  expect_equal(sort(out$merged$gene), paste0("G", 1:8))
  expect_true(all(c("log2FC_a", "log2FC_b", "padj_a", "padj_b", "category")
                 %in% colnames(out$merged)))
})

test_that("CompareMarkers categorizes shared_up / shared_down / opposite_sign / only_a / only_b correctly", {
  de <- .make_de_pair()
  out <- CompareMarkers(de$a, de$b, labels = c("a", "b"), plot = FALSE)
  cat_of <- function(g) as.character(out$merged$category[out$merged$gene == g])

  expect_equal(cat_of("G1"), "shared up")
  expect_equal(cat_of("G2"), "shared down")
  expect_equal(cat_of("G3"), "opposite sign")
  expect_equal(cat_of("G4"), "only a")
  expect_equal(cat_of("G5"), "only b")
  expect_equal(cat_of("G6"), "n.s. in both")

  expect_equal(unname(out$overlap["shared_up"]),     1)
  expect_equal(unname(out$overlap["shared_down"]),   1)
  expect_equal(unname(out$overlap["opposite_sign"]), 1)
  expect_equal(unname(out$overlap["only_a"]),        1)
  expect_equal(unname(out$overlap["only_b"]),        1)
})

test_that("CompareMarkers returns a Fisher test on the significance overlap", {
  de <- .make_de_pair()
  out <- CompareMarkers(de$a, de$b, plot = FALSE)
  expect_s3_class(out$fisher, "htest")
})

test_that("CompareMarkers plot = FALSE omits the plot; plot = TRUE returns a ggplot", {
  de <- .make_de_pair()
  out_noplot <- CompareMarkers(de$a, de$b, plot = FALSE)
  expect_false("plot" %in% names(out_noplot))
  out_plot <- CompareMarkers(de$a, de$b, plot = TRUE)
  expect_s3_class(out_plot$plot, "ggplot")
})

test_that("CompareMarkers falls back to rownames when no gene column is found", {
  de <- .make_de_pair()
  a2 <- de$a; rownames(a2) <- a2$gene; a2$gene <- NULL
  b2 <- de$b; rownames(b2) <- b2$Gene; b2$Gene <- NULL
  out <- CompareMarkers(a2, b2, plot = FALSE)
  expect_true(all(paste0("G", 1:8) %in% out$merged$gene))
})

test_that("CompareMarkers validates inputs", {
  de <- .make_de_pair()
  expect_error(CompareMarkers(list(1), de$b), "data frames")
  expect_error(CompareMarkers(de$a, de$b, labels = "only_one"), "length-2")
})

test_that("CompareMarkers errors when a fold-change/p-value column can't be auto-detected", {
  bad <- data.frame(gene = "G1", weird_col = 1, stringsAsFactors = FALSE)
  de <- .make_de_pair()
  expect_error(CompareMarkers(bad, de$b, plot = FALSE), "auto-detect")
})


# ============================================================================
# CellSuiteSummary()
# ============================================================================

test_that("CellSuiteSummary reports basic object stats with the right class", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 15, n_cells = 40, n_clusters = 3)
  s <- CellSuiteSummary(obj, cluster_col = "seurat_clusters", verbose = FALSE)
  expect_s3_class(s, "CellSuiteSummary")
  expect_equal(s$n_cells, 40)
  expect_equal(s$n_genes, 15)
  expect_equal(s$default_assay, "RNA")
  expect_equal(nrow(s$cluster_counts), length(unique(obj$seurat_clusters)))
  expect_true(all(s$qc_summary$metric %in% c("nCount_RNA", "nFeature_RNA")))
})

test_that("CellSuiteSummary includes sample_counts when sample_col is supplied", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 40, n_samples = 4)
  s <- CellSuiteSummary(obj, sample_col = "sample", verbose = FALSE)
  expect_false(is.null(s$sample_counts))
  expect_equal(sum(s$sample_counts$n_cells), 40)
})

test_that("CellSuiteSummary print method runs without erroring", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 30, n_samples = 2)
  s <- CellSuiteSummary(obj, sample_col = "sample", verbose = FALSE)
  expect_no_error(out <- capture.output(print(s)))
  expect_true(any(grepl("Seurat object summary", out)))
})

test_that("CellSuiteSummary errors on non-Seurat input", {
  expect_error(CellSuiteSummary(list(1, 2)), "Seurat object")
})

test_that("CellSuiteSummary with top_markers > 0 runs FindAllMarkers without erroring, even with no real signal", {
  .skip_if_missing("Seurat", "SeuratObject")
  # Synthetic counts are pure noise, so it's plausible FindAllMarkers finds
  # zero genes passing padj < 0.05 for some/all clusters -- this exercises
  # that edge case (rbind of empty-but-well-formed per-cluster frames)
  # rather than assuming there will be real markers.
  obj <- .make_small_seurat(seed = 3, n_genes = 20, n_cells = 60, n_clusters = 2)
  expect_no_error(
    s <- suppressWarnings(suppressMessages(
      CellSuiteSummary(obj, cluster_col = "seurat_clusters",
                       top_markers = 3, verbose = FALSE)
    ))
  )
  expect_true(is.null(s$top_markers) || is.data.frame(s$top_markers))
})


# ============================================================================
# PseudobulkDE()
# ============================================================================

.make_pseudobulk_obj <- function(seed = 1, n_genes = 40, n_per_sample = 30,
                                 n_samples_per_cond = 4, n_clusters = 2) {
  # Built by pre-allocating and filling by explicit column index, rather
  # than do.call(rbind/cbind, <named list>) -- that pattern is fragile:
  # rbind.data.frame on a *named* list of data.frames can rewrite row
  # names when reconciling the list's own tags against each piece's
  # row.names, which then silently desyncs from a separately-built
  # expression matrix's colnames (see helper-NicheCoExpress.R for the
  # same issue root-caused in detail). Filling pre-allocated structures by
  # index avoids it entirely.
  set.seed(seed)
  samples  <- paste0("D", seq_len(n_samples_per_cond * 2))
  cond_map <- setNames(rep(c("treated", "control"), each = n_samples_per_cond), samples)

  n_total <- length(samples) * n_per_sample
  all_ids <- character(n_total)
  sample_v <- character(n_total)
  condition_v <- character(n_total)
  cluster_v <- integer(n_total)
  counts <- matrix(0, nrow = n_genes, ncol = n_total,
                   dimnames = list(paste0("Gene", seq_len(n_genes)), NULL))

  for (i in seq_along(samples)) {
    sp <- samples[i]
    n  <- n_per_sample
    cols <- ((i - 1) * n + 1):(i * n)
    ids  <- paste0(sp, "_c", seq_len(n))

    all_ids[cols]     <- ids
    sample_v[cols]    <- sp
    condition_v[cols] <- cond_map[[sp]]
    cluster_v[cols]   <- sample(seq_len(n_clusters) - 1, n, replace = TRUE)

    lambda <- if (cond_map[[sp]] == "treated") 6 else 3
    counts[, cols] <- stats::rpois(n_genes * n, lambda = lambda)
  }
  colnames(counts) <- all_ids
  storage.mode(counts) <- "double"

  meta <- data.frame(
    sample           = sample_v,
    condition        = condition_v,
    seurat_clusters  = factor(cluster_v),
    row.names        = all_ids,
    stringsAsFactors = FALSE
  )

  obj <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)
  SeuratObject::LayerData(obj, assay = "RNA", layer = "data") <- log1p(counts)
  obj
}

test_that("PseudobulkDE single-group mode returns results + normalized_counts", {
  .skip_if_missing("Seurat", "SeuratObject", "DESeq2")
  obj <- .make_pseudobulk_obj(seed = 1)
  res <- suppressMessages(PseudobulkDE(
    obj, sample_col = "sample", condition_col = "condition",
    ident_1 = "treated", ident_2 = "control", verbose = FALSE
  ))
  expect_named(res, c("results", "normalized_counts"))
  expect_true(all(c("gene", "baseMean", "log2FC", "pvalue", "padj")
                 %in% colnames(res$results)))
  expect_true(nrow(res$results) > 0)
})

test_that("PseudobulkDE multi-cluster mode returns a named list per cluster", {
  .skip_if_missing("Seurat", "SeuratObject", "DESeq2")
  obj <- .make_pseudobulk_obj(seed = 2, n_clusters = 2)
  res <- suppressMessages(PseudobulkDE(
    obj, sample_col = "sample", condition_col = "condition",
    ident_1 = "treated", ident_2 = "control",
    cluster_col = "seurat_clusters", verbose = FALSE
  ))
  expect_setequal(names(res), c("0", "1"))
  expect_true(all(vapply(res, function(x) all(c("results", "normalized_counts") %in% names(x)),
                        logical(1))))
})

test_that("PseudobulkDE group_by/group_value subsets before pseudobulking", {
  .skip_if_missing("Seurat", "SeuratObject", "DESeq2")
  obj <- .make_pseudobulk_obj(seed = 3, n_clusters = 2)
  res <- suppressMessages(PseudobulkDE(
    obj, sample_col = "sample", condition_col = "condition",
    ident_1 = "treated", ident_2 = "control",
    group_by = "seurat_clusters", group_value = "0", verbose = FALSE
  ))
  expect_true(nrow(res$results) > 0)
})

test_that("PseudobulkDE requires ident_1/ident_2 or a full contrast", {
  .skip_if_missing("Seurat", "SeuratObject", "DESeq2")
  obj <- .make_pseudobulk_obj(seed = 1)
  expect_error(
    PseudobulkDE(obj, sample_col = "sample", condition_col = "condition", verbose = FALSE),
    "ident_1"
  )
})

test_that("PseudobulkDE errors on missing metadata columns", {
  .skip_if_missing("Seurat", "SeuratObject", "DESeq2")
  obj <- .make_pseudobulk_obj(seed = 1)
  expect_error(
    PseudobulkDE(obj, sample_col = "nope", condition_col = "condition",
                ident_1 = "treated", ident_2 = "control", verbose = FALSE),
    "nope"
  )
})

test_that("PseudobulkDE errors when a condition has too few pseudobulk samples", {
  .skip_if_missing("Seurat", "SeuratObject", "DESeq2")
  obj <- .make_pseudobulk_obj(seed = 1, n_samples_per_cond = 4)
  expect_error(
    PseudobulkDE(obj, sample_col = "sample", condition_col = "condition",
                ident_1 = "treated", ident_2 = "control",
                min_cells_per_sample = 10000, verbose = FALSE),
    "need >="
  )
})
