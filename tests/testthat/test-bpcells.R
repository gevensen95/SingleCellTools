# Tests for ConvertToBPCells() and the on_disk/bpcells_dir params it added
# to CreateRNAObjects()/CreateVisiumObjects()/LoadXenium2().
#
# BPCells is Suggests-only and isn't installed in this environment (or most
# CI), so these tests are written to be environment-adaptive rather than
# unconditionally skipped:
#   - Tests that need BPCells to actually run a conversion are gated with
#     skip_if_not_installed("BPCells") -- they simply won't run here.
#   - The "errors clearly when BPCells is missing" tests are gated the
#     other way (skip_if(requireNamespace("BPCells", quietly = TRUE), ...))
#     so they only run in an environment where BPCells is genuinely absent
#     -- which is this one, so these are real, meaningful assertions right
#     now, not just placeholders.
# The three loader-level "errors immediately" tests below rely on the
# BPCells check running before any real directory I/O (see the up-front
# requireNamespace() gate placed alongside each function's existing
# workers/arrow checks) -- same reasoning already used for the
# file_type/use_quantile/outs-type checks documented at the top of
# test-object-construction.R.

test_that("ConvertToBPCells errors clearly when BPCells isn't installed", {
  testthat::skip_if(requireNamespace("BPCells", quietly = TRUE),
                    "BPCells is installed; this test only applies when it's absent.")
  obj <- .make_small_seurat()
  expect_error(ConvertToBPCells(obj), "BPCells")
})

test_that("CreateRNAObjects errors immediately on on_disk = TRUE without BPCells", {
  testthat::skip_if(requireNamespace("BPCells", quietly = TRUE),
                    "BPCells is installed; this test only applies when it's absent.")
  expect_error(
    CreateRNAObjects(data_dirs = "/definitely/does/not/exist", on_disk = TRUE),
    "BPCells"
  )
})

test_that("CreateVisiumObjects errors immediately on on_disk = TRUE without BPCells", {
  testthat::skip_if(requireNamespace("BPCells", quietly = TRUE),
                    "BPCells is installed; this test only applies when it's absent.")
  expect_error(
    CreateVisiumObjects(data_dirs = "/definitely/does/not/exist", on_disk = TRUE),
    "BPCells"
  )
})

test_that("LoadXenium2 errors immediately on on_disk = TRUE without BPCells", {
  testthat::skip_if(requireNamespace("BPCells", quietly = TRUE),
                    "BPCells is installed; this test only applies when it's absent.")
  # outs = "matrix" only, so this can't instead hit the (also Suggests-only)
  # arrow gate first regardless of whether arrow happens to be installed.
  expect_error(
    LoadXenium2(data_dir = "/definitely/does/not/exist", sample_name = "test",
               outs = "matrix", on_disk = TRUE),
    "BPCells"
  )
})

test_that("ConvertToBPCells errors clearly on an assay that doesn't exist", {
  testthat::skip_if_not_installed("BPCells")
  obj <- .make_small_seurat()
  expect_error(
    ConvertToBPCells(obj, assay = "nonexistent", path = tempfile("bpcells_")),
    "not found"
  )
})

test_that("ConvertToBPCells errors clearly on a ChromatinAssay instead of attempting a blind Assay5 coercion", {
  testthat::skip_if_not_installed("BPCells")
  testthat::skip_if_not_installed("Signac")

  counts <- matrix(
    stats::rpois(20, lambda = 2), nrow = 5, ncol = 4,
    dimnames = list(
      c("chr1-100-200", "chr1-300-400", "chr1-500-600",
       "chr1-700-800", "chr1-900-1000"),
      paste0("cell", 1:4)
    )
  )
  chrom_assay <- Signac::CreateChromatinAssay(
    counts = methods::as(counts, "CsparseMatrix"), sep = c("-", "-")
  )
  obj <- SeuratObject::CreateSeuratObject(counts = chrom_assay, assay = "ATAC")

  expect_error(
    ConvertToBPCells(obj, assay = "ATAC", path = tempfile("bpcells_")),
    "ChromatinAssay"
  )
})

test_that("ConvertToBPCells round-trips counts to an on-disk BPCells matrix", {
  testthat::skip_if_not_installed("BPCells")
  obj  <- .make_small_seurat()
  path <- tempfile("bpcells_")

  out <- ConvertToBPCells(obj, path = path)

  expect_true(dir.exists(file.path(path, "counts")))
  expect_true("counts" %in% SeuratObject::Layers(out[["RNA"]]))
  expect_equal(
    as.matrix(SeuratObject::LayerData(out, layer = "counts")),
    as.matrix(SeuratObject::LayerData(obj,  layer = "counts")),
    ignore_attr = TRUE
  )
})

test_that("ConvertToBPCells warns and skips a layer that isn't present", {
  testthat::skip_if_not_installed("BPCells")
  obj  <- .make_small_seurat()
  path <- tempfile("bpcells_")
  expect_warning(
    ConvertToBPCells(obj, layers = c("counts", "scale.data"), path = path),
    "scale.data"
  )
})

test_that("ConvertToBPCells refuses to silently overwrite an existing on-disk matrix", {
  testthat::skip_if_not_installed("BPCells")
  obj  <- .make_small_seurat()
  path <- tempfile("bpcells_")
  ConvertToBPCells(obj, path = path)

  expect_error(ConvertToBPCells(obj, path = path), "overwrite")
  expect_no_error(ConvertToBPCells(obj, path = path, overwrite = TRUE))
})

test_that("ConvertToBPCells accepts a named list of objects and preserves names", {
  testthat::skip_if_not_installed("BPCells")
  objs <- list(a = .make_small_seurat(seed = 1), b = .make_small_seurat(seed = 2))
  path <- tempfile("bpcells_")

  out <- ConvertToBPCells(objs, path = path)

  expect_named(out, c("a", "b"))
  expect_true(dir.exists(file.path(path, "a", "counts")))
  expect_true(dir.exists(file.path(path, "b", "counts")))
})
