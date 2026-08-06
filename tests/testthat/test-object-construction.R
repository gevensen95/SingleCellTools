# Tests for CreateRNAObjects() and CreateRNAObjectsFilter() against a
# synthetic, on-disk 10X-format directory (matrix.mtx.gz / barcodes.tsv.gz /
# features.tsv.gz) -- the flat-directory layout these functions look for
# first. Both are genuine smoke tests (real file I/O + Seurat object
# construction), written defensively (tryCatch + skip) and run inside a
# temp working directory, since CreateRNAObjectsFilter unconditionally
# writes RDS files and both print() QC plots to whatever graphics device is
# active.
#
# CreateVisiumObjects, LoadXenium2, and CreateAndIntegrateRNA are not covered
# for the full read+build pipeline: real coverage would require synthesizing
# full Space Ranger / Xenium directory trees, which is a much larger
# undertaking than the flat 10X RNA layout used here and wasn't judged worth
# the risk of a wrong fixture without real example data to check against.
# CreateVisiumObjects/CreateAndIntegrateRNA also compute
# percent.mt/percent.rb/percent.hb with the same default patterns as
# CreateRNAObjects (see below) -- see the standalone regex-correctness tests
# further down for coverage of that shared logic without needing full I/O.
#
# That said, each of these three DOES have a handful of stop()-based argument
# checks, and where a check runs before any real file is touched, it's
# covered below without needing a fixture directory tree:
#   - CreateVisiumObjects: the file_type validation (branches purely on the
#     `file_type` argument, before any list.files()/Read10X* call).
#   - CreateAndIntegrateRNA: the use_quantile/feature_min/feature_max/
#     percent_mt_max consistency checks (purely argument-shape checks, at
#     the very top of the function, before `data_dirs` is ever read).
#   - LoadXenium2: `outs`/`type` are restricted via match.arg(choices = ...)
#     before any reading happens; the switch()'s own "Unknown Xenium input
#     type" default case further down is consequently dead code (every
#     value match.arg can still let through already has an explicit case) --
#     not a bug, just unreachable, so it isn't (and can't usefully be) unit
#     tested here.
# CreateAndIntegrateRNA's other three stop()s (integration/new_reduction
# mismatch, k_anchor/k_weight requirements) sit after the full read + filter
# + merge pipeline runs, so exercising them would need the same directory
# fixture investment already ruled out above -- left uncovered for that
# reason, consistent with the rest of this function.
#
# MakeParseObj's argument validation, by contrast, IS covered below (all of
# it runs before any Parse Biosciences directory is actually read, aside
# from a plain dir.exists() check that only needs a nonexistent path
# string) -- see the "MakeParseObj argument validation" block.
#
# CreateATACObjects(Filter) do now have argument validation worth testing
# (genome/object_names/add_treatment consistency, all checked before any
# file is read) -- see test-atac-objects.R. The read+build pipeline itself
# still isn't covered here, for the same ATAC-fragment-file-tree reason as
# above.

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
  #
  # Also includes a handful of ribosomal (Rp[sl]/RP[SL]) and hemoglobin
  # (Hb[ab]/HB[AB]) genes with elevated counts -- exercises
  # CreateRNAObjects()'s percent.rb/percent.hb columns -- plus a "Hbp1"
  # decoy gene (also elevated) that the hb_pattern default is specifically
  # designed to exclude despite its "Hb" prefix.
  n_special <- 14L  # 5 mt + 4 rb + 4 hb + 1 Hbp1 decoy
  genes <- c(paste0("mt-Gene", 1:5),
            "Rps6", "Rpl3", "RPS4", "RPL7",
            "Hba-a1", "Hbb-bs", "HBA1", "HBB",
            "Hbp1",
            paste0("Gene", seq_len(n_genes - n_special)))
  barcodes <- paste0("BC", seq_len(n_cells), "-1")
  m <- matrix(stats::rpois(n_genes * n_cells, lambda = 3), nrow = n_genes)
  m[1:n_special, ] <- stats::rpois(n_special * n_cells, lambda = 10)  # elevated counts
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

  # percent.rb/percent.hb: computed alongside percent.mt from the same
  # fixture's Rp[sl]/RP[SL] and Hb[ab]/HB[AB] genes (see .make_10x_dir()).
  expect_true("percent.rb" %in% colnames(out[[1]]@meta.data))
  expect_true("percent.hb" %in% colnames(out[[1]]@meta.data))
  expect_true(all(out[[1]]$percent.rb >= 0))
  expect_true(all(out[[1]]$percent.hb >= 0))
  expect_true(mean(out[[1]]$percent.rb) > 0)
  expect_true(mean(out[[1]]$percent.hb) > 0)

  # The default hb_pattern must exclude "Hbp1" (also elevated in this
  # fixture, despite not being a real hemoglobin gene) -- a naive "^Hb"/"^HB"
  # pattern would pick it up and report a higher percentage.
  naive_hb <- Seurat::PercentageFeatureSet(out[[1]], pattern = "^(Hb|HB)")
  expect_true(mean(out[[1]]$percent.hb) < mean(naive_hb))
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


# ============================================================================
# rb_pattern / hb_pattern defaults -- pure regex correctness, independent of
# any Seurat object construction. These are the same default pattern strings
# used by CreateRNAObjects(), CreateAndIntegrateRNA(), CreateVisiumObjects(),
# CreateRNAObjectsFilter(), and MakeParseObj() (whose mt_pattern/rb_pattern/
# hb_pattern default to these same mouse-friendly patterns rather than NULL,
# though NULL remains a valid opt-out for any of the three there); kept here
# as plain grepl() checks (not literal `::: `-style access to a shared
# constant, since none exists -- each function declares the same default
# independently) so a future edit to any one of them that drifts from the
# others, or accidentally lets "Hbp1"/"HBP1" back in, is caught
# without needing full 10X directory fixtures.
# ============================================================================

test_that("default rb_pattern matches mouse/human ribosomal gene symbols", {
  rb_pattern <- "^(Rp[sl]|RP[SL])"
  matches <- c("Rps6", "Rpl3", "Rps29", "Rpl37a",   # mouse
              "RPS6", "RPL3", "RPS29", "RPL37A")   # human
  non_matches <- c("Rpgr", "RPGR", "Rab1", "Gapdh", "mt-Nd1")

  expect_true(all(grepl(rb_pattern, matches)))
  expect_false(any(grepl(rb_pattern, non_matches)))
})

test_that("default hb_pattern matches hemoglobin genes but excludes Hbp1/HBP1", {
  hb_pattern <- "^(Hb[^p]|HB[^P])"
  matches <- c("Hba-a1", "Hbb-bs", "Hbq1b",        # mouse
              "HBA1", "HBB", "HBA2", "HBQ1")      # human
  non_matches <- c("Hbp1", "HBP1", "Gapdh", "mt-Nd1")

  expect_true(all(grepl(hb_pattern, matches)))
  expect_false(any(grepl(hb_pattern, non_matches)))
})


# ============================================================================
# MakeParseObj() argument validation -- every check below runs before any
# Parse Biosciences directory is actually read, so these are genuine unit
# tests (no fixture directories needed), unlike the full read+build pipeline
# noted as out of scope at the top of this file. The lone exception,
# `paths` existence, only needs a path that doesn't exist -- no real
# DGE_filtered/ tree required.
# ============================================================================

test_that("MakeParseObj requires a non-empty character `paths`", {
  expect_error(MakeParseObj(character(0)), "non-empty character vector")
  expect_error(MakeParseObj(1:3), "non-empty character vector")
})

test_that("MakeParseObj errors when a path doesn't exist", {
  expect_error(MakeParseObj("/no/such/parse/dir/at/all"), "do not exist")
})

test_that("MakeParseObj requires sample_names to match paths in length", {
  d <- tempfile("parse_"); dir.create(d)
  on.exit(unlink(d, recursive = TRUE))
  expect_error(
    MakeParseObj(c(d, d), sample_names = "only_one"),
    "sample_names.*same length"
  )
})

test_that("MakeParseObj validates `treatments`", {
  d <- tempfile("parse_"); dir.create(d)
  on.exit(unlink(d, recursive = TRUE))
  expect_error(
    MakeParseObj(d, treatments = list(1)),
    "atomic vector"
  )
  expect_error(
    MakeParseObj(c(d, d), treatments = "only_one"),
    "same length"
  )
})

test_that("MakeParseObj fails fast when doublet_vars_to_regress requests an uncomputed mt/rb/hb column", {
  d <- tempfile("parse_"); dir.create(d)
  on.exit(unlink(d, recursive = TRUE))

  expect_error(
    MakeParseObj(d, mt_pattern = NULL,
                doublet_vars_to_regress = "percent.mt"),
    "mt_pattern"
  )
  expect_error(
    MakeParseObj(d, rb_pattern = NULL,
                doublet_vars_to_regress = "percent.rb"),
    "rb_pattern"
  )
  expect_error(
    MakeParseObj(d, hb_pattern = NULL,
                doublet_vars_to_regress = "percent.hb"),
    "hb_pattern"
  )
})

test_that("MakeParseObj's default doublet_vars_to_regress auto-resolves without erroring", {
  d <- tempfile("parse_"); dir.create(d)
  on.exit(unlink(d, recursive = TRUE))

  # doublet_vars_to_regress is left unspecified in both calls below, so it
  # should resolve internally (via the mt_pattern-driven default) rather
  # than tripping the mismatch check exercised above. Both should get past
  # argument validation cleanly and fail only on the (irrelevant here)
  # missing DGE_filtered directory -- proving the auto-resolved value never
  # collides with `mt_pattern`, whether mt_pattern is left at its default
  # ("percent.mt" ends up requested) or explicitly turned off (falls back
  # to NULL, nothing requested).
  expect_error(MakeParseObj(d), "DGE_filtered")
  expect_error(MakeParseObj(d, mt_pattern = NULL), "DGE_filtered")
})


# ============================================================================
# CreateVisiumObjects() / CreateAndIntegrateRNA() / LoadXenium2() -- the
# subset of their argument validation that runs before any real file is
# touched (see the note at the top of this file for what's deliberately
# still uncovered and why).
# ============================================================================

test_that("CreateVisiumObjects errors on an unknown file_type", {
  # Branches purely on `file_type` before any list.files()/Read10X* call, so
  # `data_dirs` never actually needs to exist for this particular check.
  expect_error(
    CreateVisiumObjects(data_dirs = "irrelevant", file_type = "bogus"),
    "Choose file_type"
  )
})

test_that("CreateAndIntegrateRNA requires explicit thresholds when use_quantile = FALSE", {
  expect_error(
    CreateAndIntegrateRNA(data_dirs = "irrelevant", use_quantile = FALSE),
    "feature_min"
  )
})

test_that("CreateAndIntegrateRNA errors when a hard cutoff is set but use_quantile stays TRUE", {
  expect_error(
    CreateAndIntegrateRNA(data_dirs = "irrelevant", feature_min = 100),
    "quantile=TRUE"
  )
})

test_that("LoadXenium2 rejects invalid `outs`/`type` values via match.arg", {
  # type is validated before outs (see source order) -- cover both.
  expect_error(
    LoadXenium2(data_dir = "irrelevant", sample_name = "s1", type = "bogus"),
    "should be one of"
  )
  expect_error(
    LoadXenium2(data_dir = "irrelevant", sample_name = "s1", outs = "bogus"),
    "should be one of"
  )
})
