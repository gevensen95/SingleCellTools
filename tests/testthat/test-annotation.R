# Tests for AnnotateClusters() (marker mode; SingleR mode is validation-only
# since it needs an external reference package) and AnnotateWithReference()
# (scmap backend end-to-end; celltypist/scanvi are validation-only since
# they require a configured Python environment with celltypist/scvi-tools
# installed, which isn't something a unit test should assume).

# ============================================================================
# AnnotateClusters() -- marker mode
# ============================================================================

.make_annotation_obj <- function(seed = 1, n_cells = 80) {
  set.seed(seed)
  cluster_v <- rep(c(0, 1), each = n_cells / 2)
  genes <- c("MarkerA1", "MarkerA2", "MarkerB1", "MarkerB2", paste0("Filler", 1:16))
  cells <- paste0("c", seq_len(n_cells))
  counts <- matrix(stats::rpois(length(genes) * n_cells, lambda = 2),
                   nrow = length(genes), dimnames = list(genes, cells))
  idx0 <- which(cluster_v == 0)
  idx1 <- which(cluster_v == 1)
  # Boost the "correct" markers within each cluster's cells so scoring has a
  # real, unambiguous signal to recover.
  counts["MarkerA1", idx0] <- counts["MarkerA1", idx0] + stats::rpois(length(idx0), 15)
  counts["MarkerA2", idx0] <- counts["MarkerA2", idx0] + stats::rpois(length(idx0), 15)
  counts["MarkerB1", idx1] <- counts["MarkerB1", idx1] + stats::rpois(length(idx1), 15)
  counts["MarkerB2", idx1] <- counts["MarkerB2", idx1] + stats::rpois(length(idx1), 15)
  storage.mode(counts) <- "double"

  meta <- data.frame(seurat_clusters = factor(cluster_v),
                     row.names = cells, stringsAsFactors = FALSE)
  obj <- SeuratObject::CreateSeuratObject(counts = methods::as(counts, "CsparseMatrix"), meta.data = meta)
  SeuratObject::LayerData(obj, assay = "RNA", layer = "data") <- log1p(counts)
  obj
}

test_that("AnnotateClusters (marker mode) recovers the correct cell type per cluster", {
  .skip_if_missing("Seurat", "SeuratObject", "UCell")
  obj <- .make_annotation_obj(seed = 1)
  markers <- list(TypeA = c("MarkerA1", "MarkerA2"), TypeB = c("MarkerB1", "MarkerB2"))
  out <- suppressMessages(AnnotateClusters(
    obj, method = "marker", markers = markers, cluster_col = "seurat_clusters",
    filter_nonspecific = FALSE
  ))
  expect_equal(unique(out$predicted_cell_type[out$seurat_clusters == "0"]), "TypeA")
  expect_equal(unique(out$predicted_cell_type[out$seurat_clusters == "1"]), "TypeB")
})

test_that("AnnotateClusters return_scores = 'cluster'/'cell' return the documented structures", {
  .skip_if_missing("Seurat", "SeuratObject", "UCell")
  obj <- .make_annotation_obj(seed = 1)
  markers <- list(TypeA = c("MarkerA1", "MarkerA2"), TypeB = c("MarkerB1", "MarkerB2"))

  out_cluster <- suppressMessages(AnnotateClusters(
    obj, markers = markers, cluster_col = "seurat_clusters",
    filter_nonspecific = FALSE, return_scores = "cluster"
  ))
  expect_named(out_cluster, c("obj", "scores"))
  expect_equal(sort(rownames(out_cluster$scores)), c("0", "1"))
  expect_setequal(colnames(out_cluster$scores), c("TypeA", "TypeB"))

  out_cell <- suppressMessages(AnnotateClusters(
    obj, markers = markers, cluster_col = "seurat_clusters",
    filter_nonspecific = FALSE, return_scores = "cell"
  ))
  expect_equal(nrow(out_cell$scores), ncol(obj))
})

test_that("AnnotateClusters min_score threshold labels low-confidence clusters Unknown", {
  .skip_if_missing("Seurat", "SeuratObject", "UCell")
  obj <- .make_annotation_obj(seed = 1)
  markers <- list(TypeA = c("MarkerA1", "MarkerA2"), TypeB = c("MarkerB1", "MarkerB2"))
  # A hardcoded "impossibly high" threshold isn't reliable here: with only
  # 20 genes total and markers boosted this aggressively, UCell scores for
  # the matching cluster can genuinely approach the theoretical max, so a
  # fixed 0.999 cutoff might not actually exceed every cluster's score.
  # Instead, probe the real scores first (thresholds disabled), then set
  # min_score just above the observed max -- guaranteed to fail every
  # cluster regardless of how extreme this fixture's separation is.
  probe <- suppressMessages(AnnotateClusters(
    obj, markers = markers, cluster_col = "seurat_clusters",
    filter_nonspecific = FALSE, min_score = NA, min_margin = NA,
    return_scores = "cluster"
  ))
  above_max <- max(probe$scores) + 0.01
  out <- suppressMessages(AnnotateClusters(
    obj, markers = markers, cluster_col = "seurat_clusters",
    filter_nonspecific = FALSE, min_score = above_max, min_margin = NA
  ))
  expect_true(all(out$predicted_cell_type == "Unknown"))
})

test_that("AnnotateClusters tumor_mode changes defaults and warns", {
  .skip_if_missing("Seurat", "SeuratObject", "UCell")
  obj <- .make_annotation_obj(seed = 1)
  markers <- list(TypeA = c("MarkerA1", "MarkerA2"), TypeB = c("MarkerB1", "MarkerB2"))
  expect_warning(
    out <- suppressMessages(AnnotateClusters(
      obj, markers = markers, cluster_col = "seurat_clusters", tumor_mode = TRUE
    )),
    "tumor_mode"
  )
  expect_true("predicted_cell_type" %in% colnames(out@meta.data))
})

test_that("AnnotateClusters drops cell types whose markers aren't in the assay, errors if none remain", {
  .skip_if_missing("Seurat", "SeuratObject", "UCell")
  obj <- .make_annotation_obj(seed = 1)
  markers <- list(
    TypeA  = c("MarkerA1", "MarkerA2"),
    Ghost  = c("NotAGene1", "NotAGene2")
  )
  expect_warning(
    out <- suppressMessages(AnnotateClusters(
      obj, markers = markers, cluster_col = "seurat_clusters", filter_nonspecific = FALSE
    )),
    "Ghost"
  )
  expect_false("Ghost" %in% out$predicted_cell_type)

  expect_error(
    suppressMessages(AnnotateClusters(
      obj, markers = list(Ghost = c("NotAGene1")), cluster_col = "seurat_clusters"
    )),
    "No cell types"
  )
})

test_that("AnnotateClusters validates its inputs", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_annotation_obj(seed = 1)
  expect_error(
    AnnotateClusters(obj, markers = list(a = "MarkerA1"), cluster_col = "nope"),
    "nope"
  )
  expect_error(
    AnnotateClusters(obj, markers = NULL, cluster_col = "seurat_clusters"),
    "markers"
  )
})

test_that("AnnotateClusters rejects a markers element that isn't a character vector", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_annotation_obj(seed = 1)
  expect_error(
    AnnotateClusters(obj, markers = list(TypeA = 1:3, TypeB = "MarkerB1"),
                     cluster_col = "seurat_clusters"),
    "TypeA"
  )
})

test_that("AnnotateClusters visium_mode disables filter_nonspecific by default and warns", {
  .skip_if_missing("Seurat", "SeuratObject", "UCell")
  obj <- .make_annotation_obj(seed = 1)
  markers <- list(TypeA = c("MarkerA1", "MarkerA2"), TypeB = c("MarkerB1", "MarkerB2"))
  expect_warning(
    out <- suppressMessages(AnnotateClusters(
      obj, markers = markers, cluster_col = "seurat_clusters", visium_mode = TRUE
    )),
    "visium_mode"
  )
  expect_true("predicted_cell_type" %in% colnames(out@meta.data))
})

test_that("AnnotateClusters visium_mode still warns for method = 'singler' even though it has no knob to change", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_annotation_obj(seed = 1)
  # No SingleR/reference needed -- the visium_mode warning fires before
  # method-specific work starts, so this errors on the *next* check
  # (missing `reference`), not on the warning itself.
  expect_warning(
    tryCatch(
      AnnotateClusters(obj, method = "singler", cluster_col = "seurat_clusters",
                       visium_mode = TRUE),
      error = function(e) NULL
    ),
    "visium_mode"
  )
})

test_that("AnnotateClusters visium_mode respects an explicit filter_nonspecific", {
  .skip_if_missing("Seurat", "SeuratObject", "UCell")
  obj <- .make_annotation_obj(seed = 1)
  markers <- list(TypeA = c("MarkerA1", "MarkerA2"), TypeB = c("MarkerB1", "MarkerB2"))
  # filter_nonspecific = TRUE passed explicitly -- visium_mode's default
  # (FALSE) must not override it.
  expect_warning(
    out <- suppressMessages(AnnotateClusters(
      obj, markers = markers, cluster_col = "seurat_clusters",
      visium_mode = TRUE, filter_nonspecific = TRUE
    )),
    "visium_mode"
  )
  expect_true("predicted_cell_type" %in% colnames(out@meta.data))
})


# ============================================================================
# .assign_with_unassigned() -- internal per-row best-label picker shared by
# marker mode and SingleR mode
# ============================================================================

test_that(".assign_with_unassigned picks the max-scoring label per row", {
  score_mat <- matrix(c(0.8, 0.2, 0.1, 0.9), nrow = 2, byrow = TRUE,
                      dimnames = list(c("c0", "c1"), c("TypeA", "TypeB")))
  out <- SingleCellTools:::.assign_with_unassigned(
    score_mat, c("TypeA", "TypeB"), min_score = NA, min_margin = NA,
    unassigned_label = "Unknown"
  )
  expect_equal(unname(out["c0"]), "TypeA")
  expect_equal(unname(out["c1"]), "TypeB")
})

test_that(".assign_with_unassigned treats a fully-NA row as Unknown regardless of thresholds", {
  score_mat <- matrix(c(NA_real_, NA_real_, 0.1, 0.9), nrow = 2, byrow = TRUE,
                      dimnames = list(c("c0", "c1"), c("TypeA", "TypeB")))
  out <- SingleCellTools:::.assign_with_unassigned(
    score_mat, c("TypeA", "TypeB"), min_score = NA, min_margin = NA,
    unassigned_label = "Unknown"
  )
  expect_equal(unname(out["c0"]), "Unknown")
  expect_equal(unname(out["c1"]), "TypeB")
})

test_that(".assign_with_unassigned ignores a partial NA rather than propagating it to the whole row", {
  # TypeB's NA must not stop TypeA's real 0.8 score from winning c0.
  score_mat <- matrix(c(0.8, NA_real_, 0.1, 0.9), nrow = 2, byrow = TRUE,
                      dimnames = list(c("c0", "c1"), c("TypeA", "TypeB")))
  out <- SingleCellTools:::.assign_with_unassigned(
    score_mat, c("TypeA", "TypeB"), min_score = NA, min_margin = NA,
    unassigned_label = "Unknown"
  )
  expect_equal(unname(out["c0"]), "TypeA")
})

test_that(".assign_with_unassigned still applies min_score/min_margin on non-NA rows", {
  score_mat <- matrix(c(0.11, 0.10, 0.9, 0.1), nrow = 2, byrow = TRUE,
                      dimnames = list(c("low_margin", "confident"), c("TypeA", "TypeB")))
  out <- SingleCellTools:::.assign_with_unassigned(
    score_mat, c("TypeA", "TypeB"), min_score = 0.05, min_margin = 0.05,
    unassigned_label = "Unknown"
  )
  expect_equal(unname(out["low_margin"]), "Unknown")   # margin 0.01 < 0.05
  expect_equal(unname(out["confident"]), "TypeA")
})


# ============================================================================
# AnnotateWithReference() -- scmap backend end-to-end; celltypist/scanvi are
# validation-only (require an external, unlikely-to-be-configured Python env)
# ============================================================================

.make_ref_query_pair <- function(seed = 1, n_genes = 60, n_cells = 60) {
  set.seed(seed)
  genes <- paste0("Gene", seq_len(n_genes))
  make_obj <- function(label_prefix) {
    cells <- paste0(label_prefix, seq_len(n_cells))
    ct <- rep(c("Tcell", "Bcell"), each = n_cells / 2)
    m <- matrix(stats::rpois(n_genes * n_cells, lambda = 3), nrow = n_genes,
               dimnames = list(genes, cells))
    # A couple of genes strongly mark each type, so projection has a real
    # signal to recover.
    is_t <- ct == "Tcell"
    m["Gene1", is_t]  <- m["Gene1", is_t]  + stats::rpois(sum(is_t), 20)
    m["Gene2", is_t]  <- m["Gene2", is_t]  + stats::rpois(sum(is_t), 20)
    m["Gene3", !is_t] <- m["Gene3", !is_t] + stats::rpois(sum(!is_t), 20)
    m["Gene4", !is_t] <- m["Gene4", !is_t] + stats::rpois(sum(!is_t), 20)
    storage.mode(m) <- "double"
    meta <- data.frame(cell_type = ct, row.names = cells, stringsAsFactors = FALSE)
    obj <- SeuratObject::CreateSeuratObject(counts = methods::as(m, "CsparseMatrix"), meta.data = meta)
    SeuratObject::LayerData(obj, assay = "RNA", layer = "data") <- log1p(m)
    obj
  }
  list(reference = make_obj("ref"), query = make_obj("qry"))
}

test_that("AnnotateWithReference (scmap-cluster) projects labels onto the query", {
  .skip_if_missing("Seurat", "SeuratObject", "scmap", "SingleCellExperiment",
                   "SummarizedExperiment")
  pair <- .make_ref_query_pair(seed = 1)
  out <- suppressWarnings(suppressMessages(AnnotateWithReference(
    pair$query, method = "scmap", reference = pair$reference,
    ref_label_col = "cell_type", scmap_method = "cluster"
  )))
  expect_true("predicted_cell_type" %in% colnames(out@meta.data))
  expect_true("predicted_cell_type_score" %in% colnames(out@meta.data))
  expect_true(all(out$predicted_cell_type %in% c("Tcell", "Bcell", "unassigned")))
})

test_that("AnnotateWithReference min_score relabels low-confidence cells as Unknown", {
  .skip_if_missing("Seurat", "SeuratObject", "scmap", "SingleCellExperiment",
                   "SummarizedExperiment")
  pair <- .make_ref_query_pair(seed = 2)
  out <- suppressWarnings(suppressMessages(AnnotateWithReference(
    pair$query, method = "scmap", reference = pair$reference,
    ref_label_col = "cell_type", scmap_method = "cluster",
    min_score = 1.01  # above scmap's [0,1] similarity range -> everyone fails
  )))
  expect_true(all(out$predicted_cell_type == "Unknown"))
})

test_that("AnnotateWithReference errors on a non-Seurat object and an invalid method", {
  expect_error(AnnotateWithReference(list(1)), "Seurat object")
  expect_error(
    AnnotateWithReference(1, method = "not_a_method"),
    "should be one of"
  )
})

test_that("AnnotateWithReference (scmap) requires reference and ref_label_col", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_annotation_obj(seed = 1)
  expect_error(
    AnnotateWithReference(obj, method = "scmap"),
    "reference"
  )
})
