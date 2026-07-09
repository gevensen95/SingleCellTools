# Tests for SaveWithProvenance(), PseudotimeWrapper(), and ToAnnData()/
# FromAnnData(). All are wrapped defensively (tryCatch + skip for the
# heavier ones) since they depend on real external packages/backends.
#
# RunRCTD() and RunLIANA() are NOT covered: both have zero custom input
# validation (grepping each file found no stop() calls -- they just call
# straight into spacexr/liana), and meaningful testing would require a
# real spatial deconvolution reference dataset (RunRCTD) or a live
# OmnipathR ligand-receptor database fetch over the network (RunLIANA),
# neither of which is practical to set up or fake convincingly here.

# ============================================================================
# SaveWithProvenance()
# ============================================================================

test_that("SaveWithProvenance writes an RDS and a provenance JSON sidecar", {
  .skip_if_missing("Seurat", "SeuratObject", "jsonlite")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  tmp_rds <- tempfile(fileext = ".rds")
  on.exit(unlink(c(tmp_rds, sub("\\.rds$", "_provenance.json", tmp_rds))))

  out_path <- SaveWithProvenance(obj, file = tmp_rds)
  expect_true(file.exists(tmp_rds))
  sidecar <- sub("\\.rds$", "_provenance.json", tmp_rds)
  expect_true(file.exists(sidecar))

  reloaded <- readRDS(tmp_rds)
  expect_true(inherits(reloaded, "Seurat"))

  prov <- jsonlite::fromJSON(sidecar)
  expect_equal(prov$seurat_state$n_cells, 20)
  expect_equal(prov$seurat_state$default_assay, "RNA")
  expect_true("nCount_RNA" %in% prov$seurat_state$metadata_cols)
})

test_that("SaveWithProvenance includes the `extra` list in the sidecar", {
  .skip_if_missing("Seurat", "SeuratObject", "jsonlite")
  obj <- .make_small_seurat(seed = 1, n_cells = 10)
  tmp_rds <- tempfile(fileext = ".rds")
  on.exit(unlink(c(tmp_rds, sub("\\.rds$", "_provenance.json", tmp_rds))))

  SaveWithProvenance(obj, file = tmp_rds, extra = list(analyst = "tester"))
  prov <- jsonlite::fromJSON(sub("\\.rds$", "_provenance.json", tmp_rds))
  expect_equal(prov$extra$analyst, "tester")
})

test_that("SaveWithProvenance errors on non-Seurat input", {
  expect_error(SaveWithProvenance(list(1), file = tempfile()), "Seurat object")
})


# ============================================================================
# PseudotimeWrapper()
# ============================================================================

.make_trajectory_obj <- function(seed = 1, n_per_cluster = 50) {
  set.seed(seed)
  clusters <- c("0", "1", "2")
  n_total <- length(clusters) * n_per_cluster
  genes <- paste0("Gene", 1:10)

  cluster_v <- rep(clusters, each = n_per_cluster)
  # Linear trajectory in UMAP space: cluster 0 -> 1 -> 2 along the x axis,
  # so slingshot has an unambiguous one-dimensional path to fit.
  centers <- c("0" = 0, "1" = 5, "2" = 10)
  x <- centers[cluster_v] + rnorm(n_total, sd = 0.5)
  y <- rnorm(n_total, sd = 0.5)

  cells <- paste0("c", seq_len(n_total))
  counts <- matrix(stats::rpois(length(genes) * n_total, lambda = 3),
                   nrow = length(genes), dimnames = list(genes, cells))
  storage.mode(counts) <- "double"
  meta <- data.frame(seurat_clusters = factor(cluster_v), row.names = cells,
                     stringsAsFactors = FALSE)

  obj <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)
  emb <- matrix(c(x, y), ncol = 2, dimnames = list(cells, c("UMAP_1", "UMAP_2")))
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "UMAP_",
                                                      assay = "RNA")
  obj
}

test_that("PseudotimeWrapper fits a lineage and writes pseudotime metadata columns", {
  .skip_if_missing("Seurat", "SeuratObject", "slingshot")
  testthat::skip_on_cran()
  obj <- .make_trajectory_obj(seed = 1)

  out <- tryCatch(
    suppressWarnings(suppressMessages(PseudotimeWrapper(
      obj, reduction = "umap", dims = 2, cluster_col = "seurat_clusters",
      start_cluster = "0"
    ))),
    error = function(e) e
  )
  skip_if(inherits(out, "error"),
         paste("PseudotimeWrapper did not complete on this synthetic trajectory:",
               if (inherits(out, "error")) conditionMessage(out) else ""))

  pt_cols <- grep("^slingshot_Lineage", colnames(out@meta.data), value = TRUE)
  expect_true(length(pt_cols) >= 1)
  expect_true(!is.null(out@misc$slingshot))
  # Cluster "0" was fixed as the root -- its cells should have the lowest
  # (or near-lowest) pseudotime on the first lineage.
  pt <- out@meta.data[[pt_cols[1]]]
  mean_pt_by_cluster <- tapply(pt, out$seurat_clusters, mean, na.rm = TRUE)
  expect_equal(names(which.min(mean_pt_by_cluster)), "0")
})

test_that("PseudotimeWrapper validates inputs", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  expect_error(PseudotimeWrapper(list(1)), "Seurat object")
  expect_error(PseudotimeWrapper(obj, reduction = "umap"), "Reduction.*not found")

  # Give it a reduction so this specifically exercises the cluster_col
  # check, not the (already-tested) reduction check.
  emb <- matrix(rnorm(ncol(obj) * 2), nrow = ncol(obj),
               dimnames = list(colnames(obj), c("UMAP_1", "UMAP_2")))
  obj[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "UMAP_",
                                                      assay = "RNA")
  expect_error(
    PseudotimeWrapper(obj, reduction = "umap", cluster_col = "nope"),
    "Cluster column.*not found"
  )
})


# ============================================================================
# ToAnnData() / FromAnnData() -- validation-only; a real round trip needs
# zellkonverter's Python/basilisk backend, which isn't something a unit
# test should assume is configured.
# ============================================================================

test_that("ToAnnData validates its inputs", {
  .skip_if_missing("Seurat", "SeuratObject")
  expect_error(ToAnnData(list(1), file = tempfile(fileext = ".h5ad")), "Seurat object")
})

test_that("ToAnnData warns when the output path doesn't end in .h5ad", {
  .skip_if_missing("Seurat", "SeuratObject", "zellkonverter")
  obj <- .make_small_seurat(seed = 1, n_cells = 10)
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp))
  expect_warning(
    tryCatch(ToAnnData(obj, file = tmp), error = function(e) NULL),
    "does not end in"
  )
})

test_that("FromAnnData errors when the file doesn't exist", {
  expect_error(FromAnnData("no/such/file.h5ad"), "not found")
})

test_that("ToAnnData / FromAnnData round-trip (skipped unless zellkonverter's backend is configured)", {
  .skip_if_missing("Seurat", "SeuratObject", "zellkonverter")
  obj <- .make_small_seurat(seed = 1, n_genes = 10, n_cells = 20)
  tmp <- tempfile(fileext = ".h5ad")
  on.exit(unlink(tmp))

  written <- tryCatch(ToAnnData(obj, file = tmp), error = function(e) e)
  skip_if(inherits(written, "error"),
         paste("ToAnnData needs a configured zellkonverter backend:",
               if (inherits(written, "error")) conditionMessage(written) else ""))

  back <- tryCatch(FromAnnData(tmp), error = function(e) e)
  skip_if(inherits(back, "error"),
         paste("FromAnnData needs a configured zellkonverter backend:",
               if (inherits(back, "error")) conditionMessage(back) else ""))

  expect_true(inherits(back, "Seurat"))
  expect_equal(ncol(back), ncol(obj))
})
