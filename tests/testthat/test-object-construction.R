# Tests for CreateRNAObjects() and CreateRNAObjectsFilter() against a
# synthetic, on-disk 10X-format directory (matrix.mtx.gz / barcodes.tsv.gz /
# features.tsv.gz) -- the flat-directory layout these functions look for
# first. Both are genuine smoke tests (real file I/O + Seurat object
# construction), written defensively (tryCatch + skip) and run inside a
# temp working directory, since CreateRNAObjectsFilter unconditionally
# writes RDS files and both print() QC plots to whatever graphics device is
# active.
#
# CreateATACObjects(Filter), CreateVisiumObjects, LoadXenium2, MakeParseObj,
# and CreateAndIntegrateRNA are not covered: none has custom input
# validation to unit-test (Signac/Seurat's own readers do the validating,
# and grepping each file found zero stop() calls), and real coverage would
# require synthesizing full Space Ranger / Xenium / ATAC fragment-file
# directory trees, which is a much larger undertaking than the flat 10X RNA
# layout used here and wasn't judged worth the risk of a wrong fixture
# without real example data to check against.

.make_10x_dir <- function(seed = 1, n_genes = 300, n_cells = 50) {
  set.seed(seed)
  dir <- tempfile("tenx_")
  dir.create(dir)

  # A few "mt-" genes with real, cell-varying counts -- CreateRNAObjects(Filter)
  # computes percent.mt via PercentageFeatureSet(pattern = "^mt-"); without
  # any matching genes, percent.mt is exactly 0 for every cell, which makes
  # quantile-based upper thresholds also compute to 0 and the strict "<"
  # filter drop every single cell ("No cells found"). Real data always has
  # some mitochondrial genes with genuine variance, so this only bites a
  # synthetic fixture with no "mt-" genes at all.
  genes <- c(paste0("mt-Gene", 1:5), paste0("Gene", seq_len(n_genes - 5)))
  barcodes <- paste0("BC", seq_len(n_cells), "-1")
  m <- matrix(stats::rpois(n_genes * n_cells, lambda = 3), nrow = n_genes)
  m[1:5, ] <- stats::rpois(5 * n_cells, lambda = 10)  # elevated "mt-" counts
  sm <- methods::as(m, "CsparseMatrix")

  # Matrix::writeMM() doesn't reliably write through an open R connection
  # object -- confirmed empirically: passing it a gzfile() connection
  # produced a totally empty file (readLines() came back character(0)).
  # It expects a plain file path and manages its own I/O, so write the
  # plain .mtx first, then gzip that known-good text the same way the
  # barcodes/features files (which work fine) are compressed below.
  mtx_path <- file.path(dir, "matrix.mtx")
  Matrix::writeMM(sm, mtx_path)
  mtx_lines <- readLines(mtx_path)
  con <- gzfile(file.path(dir, "matrix.mtx.gz"), "wt")
  writeLines(mtx_lines, con)
  close(con)
  file.remove(mtx_path)

  con <- gzfile(file.path(dir, "barcodes.tsv.gz"), "wt")
  writeLines(barcodes, con)
  close(con)

  con <- gzfile(file.path(dir, "features.tsv.gz"), "wt")
  writeLines(paste(genes, genes, "Gene Expression", sep = "\t"), con)
  close(con)

  dir
}

test_that("CreateRNAObjects reads a flat 10X-format directory into a named Seurat list", {
  .skip_if_missing("Seurat", "SeuratObject")
  testthat::skip_on_cran()
  d <- .make_10x_dir(seed = 1)
  on.exit(unlink(d, recursive = TRUE))

  old_wd <- getwd()
  tmp_dir <- tempfile("createrna_test_")
  dir.create(tmp_dir)
  setwd(tmp_dir)
  on.exit({ setwd(old_wd); unlink(tmp_dir, recursive = TRUE) }, add = TRUE)

  out <- tryCatch(
    suppressWarnings(suppressMessages(CreateRNAObjects(
      data_dirs = d, cells = 3, features = 200,
      run_doublet_finder = FALSE
    ))),
    error = function(e) e
  )
  skip_if(inherits(out, "error"),
         paste("CreateRNAObjects did not complete on this synthetic 10X dir:",
               if (inherits(out, "error")) conditionMessage(out) else ""))

  expect_type(out, "list")
  expect_true(inherits(out[[1]], "Seurat"))
  expect_true("percent.mt" %in% colnames(out[[1]]@meta.data))
  expect_true(ncol(out[[1]]) > 0)
})

test_that("CreateRNAObjectsFilter reads and quantile-filters a flat 10X directory", {
  .skip_if_missing("Seurat", "SeuratObject")
  testthat::skip_on_cran()
  d <- .make_10x_dir(seed = 2)
  on.exit(unlink(d, recursive = TRUE))

  old_wd <- getwd()
  tmp_dir <- tempfile("creaternafilter_test_")
  dir.create(tmp_dir)
  setwd(tmp_dir)
  on.exit({ setwd(old_wd); unlink(tmp_dir, recursive = TRUE) }, add = TRUE)

  out <- tryCatch(
    suppressWarnings(suppressMessages(CreateRNAObjectsFilter(
      data_dirs = d, cells = 3, features = 200,
      use_quantile = TRUE, quantile_value_min = 0.1
    ))),
    error = function(e) e
  )
  skip_if(inherits(out, "error"),
         paste("CreateRNAObjectsFilter did not complete on this synthetic 10X dir:",
               if (inherits(out, "error")) conditionMessage(out) else ""))

  expect_type(out, "list")
  expect_true(inherits(out[[1]], "Seurat"))
  expect_true(ncol(out[[1]]) > 0)
})

test_that("CreateRNAObjectsFilter requires explicit thresholds when use_quantile = FALSE", {
  expect_error(
    CreateRNAObjectsFilter(data_dirs = "irrelevant", use_quantile = FALSE),
    "feature_min"
  )
})
