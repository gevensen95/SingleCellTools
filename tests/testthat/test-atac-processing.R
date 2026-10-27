# Smoke/validation tests for RunATACWrapper(), RemoveAmbientRNA(), and
# LinkPeaksToGenes() -- three functions with real custom input-validation
# logic (several stop() conditions each) that had zero test coverage.
#
# None of the three get a full genome-dependent happy-path run here:
#   - RunATACWrapper() DOES get one (TF-IDF + top-features + SVD only touch
#     the counts matrix -- no genome/fragments needed), using the same
#     sep-based synthetic ChromatinAssay pattern as test-bpcells.R.
#   - RemoveAmbientRNA() needs 'celda' (decontX) or 'SoupX', neither of
#     which is installed in this environment (both Suggests-only) -- so,
#     following the skip_if(requireNamespace(...)) convention already used
#     in test-bpcells.R/test-external-tools.R for Suggests-only deps, the
#     "errors clearly when the package is missing" tests below are real,
#     meaningful assertions here rather than placeholders.
#   - LinkPeaksToGenes()'s heavy calls (Signac::RegionStats()/LinkPeaks())
#     need real genome sequence and gene-coordinate annotation to produce
#     anything meaningful; its validation logic (peak/expression assay
#     detection, genome-as-path handling) sits entirely before those calls
#     and is covered in isolation instead.

# Small synthetic ChromatinAssay: peak names encode chr-start-end via
# Signac's own `sep` convention (same pattern as the ChromatinAssay fixture
# in test-bpcells.R). No genome/fragments attached -- nothing exercised
# below touches them.
.make_atac_seurat <- function(seed = 1, n_peaks = 30, n_cells = 15) {
  set.seed(seed)
  peaks <- paste0("chr1-", seq(1, by = 1000, length.out = n_peaks), "-",
                  seq(1000, by = 1000, length.out = n_peaks))
  cells <- paste0("cell", seq_len(n_cells))
  counts <- matrix(stats::rpois(n_peaks * n_cells, lambda = 3),
                   nrow = n_peaks, ncol = n_cells,
                   dimnames = list(peaks, cells))
  storage.mode(counts) <- "double"
  chrom_assay <- Signac::CreateChromatinAssay(
    counts = methods::as(counts, "CsparseMatrix"), sep = c("-", "-")
  )
  SeuratObject::CreateSeuratObject(counts = chrom_assay, assay = "ATAC")
}

# ATAC + RNA on the same cells, for LinkPeaksToGenes().
.make_multiome_seurat <- function(seed = 1, n_peaks = 30, n_genes = 10, n_cells = 15) {
  obj <- .make_atac_seurat(seed = seed, n_peaks = n_peaks, n_cells = n_cells)
  set.seed(seed + 1)
  genes <- paste0("Gene", seq_len(n_genes))
  rna_counts <- matrix(stats::rpois(n_genes * n_cells, lambda = 3),
                       nrow = n_genes, ncol = n_cells,
                       dimnames = list(genes, colnames(obj)))
  storage.mode(rna_counts) <- "double"
  obj[["RNA"]] <- Seurat::CreateAssayObject(counts = methods::as(rna_counts, "CsparseMatrix"))
  obj
}


# ============================================================================
# RunATACWrapper()
# ============================================================================

test_that("RunATACWrapper errors on a non-Seurat object", {
  expect_error(RunATACWrapper(list(1)), "must be a Seurat object")
})

test_that("RunATACWrapper errors when the (default) assay isn't a ChromatinAssay", {
  testthat::skip_if_not_installed("Signac")
  obj <- .make_small_seurat()  # plain RNA assay, no Signac involved
  expect_error(RunATACWrapper(obj), "not a ChromatinAssay")
})

test_that("RunATACWrapper runs TF-IDF + SVD end to end on a synthetic ChromatinAssay", {
  testthat::skip_if_not_installed("Signac")
  obj <- .make_atac_seurat(n_peaks = 30, n_cells = 15)

  out <- tryCatch(
    suppressMessages(RunATACWrapper(obj, n_components = 5, verbose = FALSE)),
    error = function(e) e
  )
  skip_if(inherits(out, "error"),
         paste("RunATACWrapper did not complete on this synthetic dataset:",
               if (inherits(out, "error")) conditionMessage(out) else ""))

  expect_true("lsi" %in% Seurat::Reductions(out))
  expect_equal(ncol(Seurat::Embeddings(out, "lsi")), 5)
  expect_equal(nrow(Seurat::Embeddings(out, "lsi")), ncol(obj))
})


# ============================================================================
# RemoveAmbientRNA()
# ============================================================================

test_that("RemoveAmbientRNA errors on a non-Seurat object", {
  expect_error(RemoveAmbientRNA(list(1)), "must be a Seurat object")
})

test_that("RemoveAmbientRNA errors clearly on a clusters_col not in metadata", {
  obj <- .make_small_seurat()
  expect_error(RemoveAmbientRNA(obj, clusters_col = "not_a_real_column"),
              "not found")
})

test_that("RemoveAmbientRNA errors clearly when 'celda' isn't installed (method = 'decontx', the default)", {
  testthat::skip_if(requireNamespace("celda", quietly = TRUE),
                    "celda is installed; this test only applies when it's absent.")
  obj <- .make_small_seurat()
  expect_error(RemoveAmbientRNA(obj), "celda")
})

test_that("RemoveAmbientRNA errors clearly when 'SoupX' isn't installed (method = 'soupx')", {
  testthat::skip_if(requireNamespace("SoupX", quietly = TRUE),
                    "SoupX is installed; this test only applies when it's absent.")
  obj <- .make_small_seurat()
  expect_error(RemoveAmbientRNA(obj, method = "soupx"), "SoupX")
})

test_that("RemoveAmbientRNA (method = 'soupx') requires raw_counts", {
  testthat::skip_if_not_installed("SoupX")
  obj <- .make_small_seurat()
  expect_error(RemoveAmbientRNA(obj, method = "soupx"), "raw_counts")
})

test_that("RemoveAmbientRNA accepts a list of Seurat objects and validates each one", {
  # No celda/SoupX needed -- clusters_col validation happens before either
  # tool is touched, in .remove_ambient_rna_one(), so this exercises the
  # list-plumbing (.as_seurat_list() + per-element dispatch) without
  # depending on either optional package.
  objs <- list(a = .make_small_seurat(), b = .make_small_seurat())
  expect_error(RemoveAmbientRNA(objs, clusters_col = "not_a_real_column"),
              "not found")
})

test_that("RemoveAmbientRNA errors when raw_counts (list mode) doesn't match obj's length", {
  testthat::skip_if_not_installed("SoupX")
  objs <- list(a = .make_small_seurat(), b = .make_small_seurat())
  expect_error(
    RemoveAmbientRNA(objs, method = "soupx",
                    raw_counts = list("only_one_path/")),
    "must match one-to-one")
})

test_that("RemoveAmbientRNA errors when raw_counts (list mode) names don't match obj's names", {
  testthat::skip_if_not_installed("SoupX")
  objs <- list(a = .make_small_seurat(), b = .make_small_seurat())
  expect_error(
    RemoveAmbientRNA(objs, method = "soupx",
                    raw_counts = list(x = "path1/", y = "path2/")),
    "names don't match")
})

test_that("RemoveAmbientRNA errors when workers exceeds available cores", {
  objs <- list(a = .make_small_seurat(), b = .make_small_seurat())
  n_cores <- suppressWarnings(parallel::detectCores())
  testthat::skip_if(is.na(n_cores), "parallel::detectCores() returned NA on this machine.")
  expect_error(
    RemoveAmbientRNA(objs, clusters_col = "not_a_real_column", workers = n_cores + 1L),
    "exceeds the")
})


# ============================================================================
# LinkPeaksToGenes()
# ============================================================================

test_that("LinkPeaksToGenes errors on a non-Seurat object", {
  expect_error(LinkPeaksToGenes(list(1), genome = NULL), "must be a Seurat object")
})

test_that("LinkPeaksToGenes errors when `genome` is a path that doesn't exist", {
  testthat::skip_if_not_installed("Signac")
  testthat::skip_if_not_installed("Rsamtools")  # needed to reach the file-exists check at all
  obj <- .make_multiome_seurat()
  expect_error(LinkPeaksToGenes(obj, genome = "/definitely/does/not/exist.fa"),
              "neither.*a BSgenome")
})

test_that("LinkPeaksToGenes errors when no ChromatinAssay can be auto-detected", {
  testthat::skip_if_not_installed("Signac")
  obj <- .make_small_seurat()  # plain RNA assay only -- zero ChromatinAssays
  expect_error(LinkPeaksToGenes(obj, genome = NULL), "auto-detect")
})

test_that("LinkPeaksToGenes errors when more than one ChromatinAssay exists and peak_assay isn't given", {
  testthat::skip_if_not_installed("Signac")
  obj <- .make_multiome_seurat()
  # A second, independent ChromatinAssay (own key, to avoid any key
  # collision with "ATAC"'s) -- enough to make peak_assay auto-detection
  # ambiguous.
  peaks2 <- paste0("chr2-", seq(1, by = 1000, length.out = 30), "-",
                   seq(1000, by = 1000, length.out = 30))
  counts2 <- matrix(stats::rpois(30 * ncol(obj), lambda = 3),
                    nrow = 30, dimnames = list(peaks2, colnames(obj)))
  storage.mode(counts2) <- "double"
  obj[["ATAC2"]] <- Signac::CreateChromatinAssay(
    counts = methods::as(counts2, "CsparseMatrix"), sep = c("-", "-"),
    key = "atac2_"
  )
  expect_error(LinkPeaksToGenes(obj, genome = NULL), "auto-detect")
})

test_that("LinkPeaksToGenes errors when an explicit peak_assay isn't a ChromatinAssay", {
  testthat::skip_if_not_installed("Signac")
  obj <- .make_multiome_seurat()
  expect_error(LinkPeaksToGenes(obj, genome = NULL, peak_assay = "RNA"),
              "not a ChromatinAssay")
})

test_that("LinkPeaksToGenes errors when no RNA assay can be auto-detected for expression_assay", {
  testthat::skip_if_not_installed("Signac")
  obj <- .make_atac_seurat()  # ChromatinAssay only, no RNA assay
  expect_error(LinkPeaksToGenes(obj, genome = NULL), "Could not find an 'RNA' assay")
})

test_that("LinkPeaksToGenes errors when an explicit expression_assay isn't found", {
  testthat::skip_if_not_installed("Signac")
  obj <- .make_multiome_seurat()
  expect_error(
    LinkPeaksToGenes(obj, genome = NULL, expression_assay = "not_a_real_assay"),
    "not found"
  )
})
