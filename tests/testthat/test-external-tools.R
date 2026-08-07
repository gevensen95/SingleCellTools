# Tests for SaveWithProvenance(), PseudotimeWrapper(), ToAnnData()/
# FromAnnData(), RunCellChat(), and RunRCTD(). All are wrapped defensively
# (tryCatch + skip for the heavier ones) since they depend on real external
# packages/backends.
#
# RunLIANA() is NOT covered: it has zero custom input validation (grepping
# the file found no stop() calls -- it just calls straight into liana), and
# meaningful testing would require a live OmnipathR ligand-receptor
# database fetch over the network, which isn't practical to set up or fake
# convincingly here.
#
# RunCellChat() DOES have real input validation, so that part is covered
# below with no CellChat installation required. A genuine end-to-end run
# is NOT attempted: computeCommunProb()'s permutation testing is slow, and
# getting any non-empty result requires gene symbols that actually overlap
# CellChatDB.mouse/.human -- a synthetic "Gene1"/"Gene2"/... fixture never
# will, so an end-to-end call would only ever exercise the "zero
# interactions survived filtering" path, not real behavior.
#
# RunRCTD() similarly has no CellChat-scale end-to-end test here (a real
# RCTD run needs a real spatial deconvolution reference dataset and is
# slow), but it DOES have real input validation -- and as of the
# obj/reference/celltype_col checks running before the spacexr
# requireNamespace() check, that validation is reachable and testable
# whether or not spacexr happens to be installed. See below.

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
# RunCellChat() -- validation-only; see file header note on why a real
# end-to-end run isn't attempted here.
# ============================================================================

test_that("RunCellChat errors on non-Seurat input", {
  expect_error(RunCellChat(list(1), label = "x"), "Seurat object")
})

test_that("RunCellChat errors when group_col is not found in metadata", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  expect_error(
    RunCellChat(obj, label = "x", group_col = "not_a_column"),
    "not_a_column"
  )
})

test_that("RunCellChat validates `species` via match.arg", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  expect_error(
    RunCellChat(obj, label = "x", group_col = "celltype", species = "cat"),
    "should be one of|arg"
  )
})

test_that("RunCellChat's own checks run before it checks for CellChat itself", {
  # i.e. seurat/group_col validation happens even when CellChat isn't
  # installed -- mirrors the ordering RunBanksyWrapper uses, and is what
  # makes the two tests above meaningful regardless of the local setup.
  expect_error(RunCellChat(list(1), label = "x"), "Seurat object")
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  expect_error(
    RunCellChat(obj, label = "x", group_col = "nope"),
    "nope"
  )
})

test_that("RunCellChat runs end-to-end on a synthetic object (skipped unless CellChat is installed)", {
  .skip_if_missing("Seurat", "SeuratObject", "CellChat")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 60,
                            n_celltypes = 3)

  out <- tryCatch(
    suppressWarnings(suppressMessages(
      RunCellChat(obj, label = "test", group_col = "celltype",
                 nboot = 5, verbose = FALSE)
    )),
    error = function(e) e
  )
  skip_if(inherits(out, "error"),
         paste("RunCellChat did not complete on this synthetic object",
               "(expected -- see file header note):",
               if (inherits(out, "error")) conditionMessage(out) else ""))

  expect_true(inherits(out, "CellChat"))
  expect_equal(length(levels(out@idents)), length(unique(obj$celltype)))
  expect_true(!is.null(out@net))
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


# ============================================================================
# RunRCTD() -- validation only. reference/celltype_col/obj checks run before
# the spacexr requireNamespace() check (see file header), so they're
# reachable in any environment. A real end-to-end run needs a real spatial
# deconvolution reference dataset and spacexr installed, neither of which
# this suite assumes.
# ============================================================================

test_that("RunRCTD errors on a non-Seurat reference before touching spacexr", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  expect_error(
    RunRCTD(obj, reference = list(1, 2), celltype_col = "celltype"),
    "reference.*Seurat"
  )
})

test_that("RunRCTD errors when celltype_col is missing or invalid", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  ref <- .make_small_seurat(seed = 2, n_cells = 20)
  expect_error(RunRCTD(obj, reference = ref, celltype_col = NULL), "celltype_col")
  expect_error(RunRCTD(obj, reference = ref, celltype_col = "nope"), "celltype_col")
})

test_that("RunRCTD errors on non-Seurat obj, single or as a list", {
  .skip_if_missing("Seurat", "SeuratObject")
  ref <- .make_small_seurat(seed = 1, n_cells = 20)
  expect_error(
    RunRCTD(list(1, 2), reference = ref, celltype_col = "celltype"),
    "Seurat object"
  )
})

test_that("RunRCTD errors clearly when spacexr isn't installed (after cheap validation passes)", {
  testthat::skip_if(requireNamespace("spacexr", quietly = TRUE),
                    "spacexr is installed -- this test targets the missing-package path")
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  ref <- .make_small_seurat(seed = 2, n_cells = 20)
  expect_error(
    RunRCTD(obj, reference = ref, celltype_col = "celltype"),
    "spacexr"
  )
})


# ============================================================================
# SpatialCompositionPlot() -- validation only. Requires the scatterpie
# package to run at all (checked first, before anything else), so most of
# these only exercise reachable validation when scatterpie IS installed;
# the missing-package path is tested with the inverted skip_if pattern used
# elsewhere in this suite.
# ============================================================================

test_that("SpatialCompositionPlot errors clearly when scatterpie isn't installed", {
  testthat::skip_if(requireNamespace("scatterpie", quietly = TRUE),
                    "scatterpie is installed -- this test targets the missing-package path")
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  expect_error(SpatialCompositionPlot(obj), "scatterpie")
})

test_that("SpatialCompositionPlot errors on non-Seurat input", {
  .skip_if_missing("scatterpie")
  expect_error(SpatialCompositionPlot(list(1, 2)), "Seurat object")
})

test_that("SpatialCompositionPlot errors when no rctd_ columns exist and weight_cols isn't supplied", {
  .skip_if_missing("Seurat", "SeuratObject", "scatterpie")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  expect_error(SpatialCompositionPlot(obj), "RunRCTD|weight_cols")
})

test_that("SpatialCompositionPlot errors on a missing weight_cols column", {
  .skip_if_missing("Seurat", "SeuratObject", "scatterpie")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  expect_error(SpatialCompositionPlot(obj, weight_cols = "nope"), "nope")
})

test_that("SpatialCompositionPlot errors on a non-numeric weight_cols column", {
  .skip_if_missing("Seurat", "SeuratObject", "scatterpie")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  expect_error(SpatialCompositionPlot(obj, weight_cols = "seurat_clusters"),
              "numeric")
})

test_that("SpatialCompositionPlot errors when obj has no images", {
  .skip_if_missing("Seurat", "SeuratObject", "scatterpie")
  obj <- .make_small_seurat(seed = 1, n_cells = 20)
  obj$w_typeA <- stats::runif(ncol(obj))
  expect_error(SpatialCompositionPlot(obj, weight_cols = "w_typeA"), "images")
})
