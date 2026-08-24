# Tests for CombineParseRounds(), CombineCellRangerRounds(), and the shared
# format-agnostic core they both call (combine-rounds-utils.R).
#
# .reindex() and .combine_one_core() are plain package-internal functions
# (not nested inside either exported function), so they're unit-tested
# directly with small in-memory matrices below -- this is the most reliable
# way to pin down the actual union/reindex/sum/metadata-reconciliation
# arithmetic, independent of any file-format detail.
#
# CombineParseRounds()'s and CombineCellRangerRounds()'s own argument
# validation is covered without any fixture, following the same pattern used
# for MakeParseObj() in test-object-construction.R. The full read -> combine
# -> write round trip for each platform is covered with a small synthetic
# on-disk fixture (Parse's DGE_filtered layout / a CellRanger triplet,
# respectively) as a genuine smoke test, written defensively (tryCatch +
# skip_if), matching this suite's existing precedent for real file I/O.
#
# .read_round() / .read_round_cellranger() / .combine_one() / .write_one()
# are closures defined *inside* each exported function (not top-level
# package functions), so they can't be unit-tested in isolation -- they're
# only exercised indirectly through the full-function smoke tests below.


# ============================================================================
# .reindex() -- pure matrix remapping, no combining involved yet
# ============================================================================

test_that(".reindex remaps rows onto a superset of keys with zero-fill", {
  m <- Matrix::Matrix(c(1, 2, 3, 4), nrow = 2,
                      dimnames = list(c("g1", "g2"), c("c1", "c2")),
                      sparse = TRUE)
  out <- .reindex(m, old_keys = c("g1", "g2"),
                  new_keys = c("g0", "g1", "g2", "g3"), margin = "rows")

  expect_equal(rownames(out), c("g0", "g1", "g2", "g3"))
  expect_equal(colnames(out), colnames(m))
  expect_equal(as.numeric(out["g1", ]), as.numeric(m["g1", ]))
  expect_equal(as.numeric(out["g2", ]), as.numeric(m["g2", ]))
  expect_true(all(out["g0", ] == 0))
  expect_true(all(out["g3", ] == 0))
})

test_that(".reindex remaps columns onto a superset of keys with zero-fill", {
  m <- Matrix::Matrix(c(1, 2, 3, 4), nrow = 2,
                      dimnames = list(c("g1", "g2"), c("c1", "c2")),
                      sparse = TRUE)
  out <- .reindex(m, old_keys = colnames(m),
                  new_keys = c("c0", "c1", "c2", "c3"), margin = "cols")

  expect_equal(colnames(out), c("c0", "c1", "c2", "c3"))
  expect_equal(rownames(out), rownames(m))
  expect_equal(as.numeric(out[, "c1"]), as.numeric(m[, "c1"]))
  expect_equal(as.numeric(out[, "c2"]), as.numeric(m[, "c2"]))
  expect_true(all(out[, "c0"] == 0))
  expect_true(all(out[, "c3"] == 0))
})


# ============================================================================
# .combine_one_core() -- the format-agnostic union/reindex/sum/metadata core,
# exercised directly with small hand-verifiable in-memory rounds. Both rounds
# below deliberately overlap on genes {g2, g3} and barcodes {bcB, bcC}, and
# each has one gene and one barcode the other round doesn't -- exactly the
# "resequenced, not re-prepped" shape these functions exist for.
# ============================================================================

# `counts_mat` is supplied genes x cells (rows named by `genes`, in the same
# order as columns named by `barcodes`).
.make_core_round <- function(genes, barcodes, counts_mat, extra_cell_cols = list()) {
  m <- methods::as(counts_mat, "CsparseMatrix")
  dimnames(m) <- list(genes, barcodes)
  cells_df <- data.frame(cell_id = barcodes, stringsAsFactors = FALSE)
  for (nm in names(extra_cell_cols)) cells_df[[nm]] <- extra_cell_cols[[nm]]
  list(
    counts = m,
    genes  = data.frame(gene_id = genes, stringsAsFactors = FALSE),
    cells  = cells_df,
    path   = "unit-test"
  )
}

test_that(".combine_one_core sums shared barcodes, passes through round-only barcodes, unions genes, and reconciles metadata", {
  r1 <- .make_core_round(
    genes      = c("g1", "g2", "g3"),
    barcodes   = c("bcA", "bcB", "bcC"),
    counts_mat = matrix(c(1, 0, 5,   2, 1, 0,   0, 3, 0), nrow = 3, byrow = TRUE),
    extra_cell_cols = list(species = c("human", "human", "human"))
  )
  r2 <- .make_core_round(
    genes      = c("g2", "g3", "g4"),
    barcodes   = c("bcB", "bcC", "bcD"),
    counts_mat = matrix(c(10, 0, 1,   0, 2, 0,   1, 1, 1), nrow = 3, byrow = TRUE),
    extra_cell_cols = list(species = c("human", "mouse", "mouse"))
  )

  out <- suppressWarnings(.combine_one_core(
    r1, r2, "SampleA",
    gene_key_col = "gene_id", cell_id_col = "cell_id",
    metadata_conflict = "warn"
  ))

  expect_equal(sort(colnames(out$counts)), c("bcA", "bcB", "bcC", "bcD"))
  expect_equal(sort(rownames(out$counts)), c("g1", "g2", "g3", "g4"))

  m <- as.matrix(out$counts)
  # bcA: round1-only barcode -- passed through unchanged, g4 (round2-only
  # gene) zero-filled.
  expect_equal(unname(m["g1", "bcA"]), 1)
  expect_equal(unname(m["g2", "bcA"]), 2)
  expect_equal(unname(m["g4", "bcA"]), 0)
  # bcB/bcC: shared barcodes -- overlapping genes summed across rounds;
  # round-unique genes contribute only from the round that has them.
  expect_equal(unname(m["g2", "bcB"]), 1 + 10)
  expect_equal(unname(m["g3", "bcB"]), 3 + 0)
  expect_equal(unname(m["g1", "bcB"]), 0)      # round2 has no g1
  expect_equal(unname(m["g4", "bcB"]), 0 + 1)  # round1 has no g4
  expect_equal(unname(m["g1", "bcC"]), 5 + 0)
  expect_equal(unname(m["g3", "bcC"]), 0 + 2)
  # bcD: round2-only barcode.
  expect_equal(unname(m["g2", "bcD"]), 1)
  expect_equal(unname(m["g1", "bcD"]), 0)

  expect_equal(out$summary$n_cells_round1, 3)
  expect_equal(out$summary$n_cells_round2, 3)
  expect_equal(out$summary$n_cells_shared, 2)
  expect_equal(out$summary$n_cells_round1_only, 1)
  expect_equal(out$summary$n_cells_round2_only, 1)
  expect_equal(out$summary$n_cells_combined, 4)
  expect_equal(out$summary$n_genes_round1, 3)
  expect_equal(out$summary$n_genes_round2, 3)
  expect_equal(out$summary$n_genes_combined, 4)

  # gene_count/tscp_count are recomputed from the combined matrix, not
  # passed through from either round.
  expect_equal(out$cells$tscp_count[out$cells$cell_id == "bcB"], sum(m[, "bcB"]))
  expect_equal(out$cells$gene_count[out$cells$cell_id == "bcB"], sum(m[, "bcB"] > 0))

  # Non-count metadata: round-only barcodes pass through their one round's
  # value; a shared barcode with agreeing rounds keeps that value; a shared
  # barcode with disagreeing rounds keeps round 1's (see warning test below).
  by_id <- setNames(out$cells$species, out$cells$cell_id)
  expect_equal(unname(by_id["bcA"]), "human")
  expect_equal(unname(by_id["bcD"]), "mouse")
  expect_equal(unname(by_id["bcB"]), "human")
  expect_equal(unname(by_id["bcC"]), "human")
})

test_that(".combine_one_core warns on a metadata disagreement and keeps round 1's value", {
  r1 <- .make_core_round(genes = "g1", barcodes = "bcX",
                         counts_mat = matrix(1, nrow = 1),
                         extra_cell_cols = list(species = "human"))
  r2 <- .make_core_round(genes = "g1", barcodes = "bcX",
                         counts_mat = matrix(1, nrow = 1),
                         extra_cell_cols = list(species = "mouse"))

  expect_warning(
    out <- .combine_one_core(r1, r2, "S", gene_key_col = "gene_id",
                             cell_id_col = "cell_id", metadata_conflict = "warn"),
    "disagrees between rounds"
  )
  expect_equal(out$cells$species, "human")
})

test_that(".combine_one_core errors on a metadata disagreement when metadata_conflict = 'error'", {
  r1 <- .make_core_round(genes = "g1", barcodes = "bcX",
                         counts_mat = matrix(1, nrow = 1),
                         extra_cell_cols = list(species = "human"))
  r2 <- .make_core_round(genes = "g1", barcodes = "bcX",
                         counts_mat = matrix(1, nrow = 1),
                         extra_cell_cols = list(species = "mouse"))

  expect_error(
    .combine_one_core(r1, r2, "S", gene_key_col = "gene_id",
                      cell_id_col = "cell_id", metadata_conflict = "error"),
    "disagrees between rounds"
  )
})

test_that(".combine_one_core keeps a metadata column that exists in round 2 only", {
  # Regression: selecting round 2's value on `!is.na(idx1)` alone (barcode
  # present in round 1) rather than "round 1 has the column AND the barcode"
  # wrote NA over round 2's value on every SHARED barcode, for any column
  # round 1 didn't carry -- silent data loss.
  r1 <- .make_core_round("g1", c("bcX", "bcY"), matrix(c(1, 2), nrow = 1))
  r2 <- .make_core_round("g1", c("bcX", "bcZ"), matrix(c(3, 4), nrow = 1),
                         extra_cell_cols = list(cell_class = c("neuron", "glia")))

  out <- .combine_one_core(r1, r2, "S", gene_key_col = "gene_id",
                           cell_id_col = "cell_id", metadata_conflict = "warn")
  by_id <- setNames(out$cells$cell_class, out$cells$cell_id)
  expect_equal(unname(by_id["bcX"]), "neuron")  # shared barcode, round-2-only column
  expect_equal(unname(by_id["bcZ"]), "glia")    # round-2-only barcode
  expect_true(is.na(by_id[["bcY"]]))            # round-1-only barcode, no source
})

test_that(".combine_one_core drops Seurat-style per-assay count columns", {
  # Regression: "nCount"/"nFeature" were matched as exact names, so Seurat's
  # actual "nCount_RNA"/"nFeature_RNA" survived as stale round-1 values
  # contradicting the recomputed tscp_count/gene_count next to them.
  r1 <- .make_core_round("g1", "bcX", matrix(5, nrow = 1),
                         extra_cell_cols = list(nCount_RNA = 5, nFeature_RNA = 1))
  r2 <- .make_core_round("g1", "bcX", matrix(7, nrow = 1),
                         extra_cell_cols = list(nCount_RNA = 7, nFeature_RNA = 1))

  out <- .combine_one_core(r1, r2, "S", gene_key_col = "gene_id",
                           cell_id_col = "cell_id", metadata_conflict = "warn")
  expect_false("nCount_RNA" %in% colnames(out$cells))
  expect_false("nFeature_RNA" %in% colnames(out$cells))
  expect_equal(out$cells$tscp_count, 12)
})

test_that(".combine_one_core handles factor metadata columns with different levels", {
  # Regression: `v1 != v2` on two factors with different level sets threw
  # "level sets of factors are different" and aborted the whole combine.
  r1 <- .make_core_round("g1", c("bcX", "bcY"), matrix(c(1, 2), nrow = 1),
                         extra_cell_cols = list(grp = factor(c("a", "b"))))
  r2 <- .make_core_round("g1", "bcX", matrix(3, nrow = 1),
                         extra_cell_cols = list(grp = factor("a")))

  out <- .combine_one_core(r1, r2, "S", gene_key_col = "gene_id",
                           cell_id_col = "cell_id", metadata_conflict = "warn")
  expect_type(out$cells$grp, "character")
  expect_equal(out$cells$grp[out$cells$cell_id == "bcX"], "a")
  expect_equal(out$cells$grp[out$cells$cell_id == "bcY"], "b")
})

test_that(".combine_one_core fills round-2-only genes by column NAME, not position", {
  # Regression: `genes_out[need_r2, ] <- r2$genes[...]` assigns positionally,
  # so a differing gene-table column ORDER silently transposed values between
  # columns, and a differing column COUNT silently recycled them.
  mk <- function(g, b, v, genes_df) {
    m <- methods::as(matrix(v, nrow = length(g)), "CsparseMatrix")
    dimnames(m) <- list(g, b)
    list(counts = m, genes = genes_df,
         cells = data.frame(cell_id = b, stringsAsFactors = FALSE), path = "x")
  }
  r1 <- mk("g1", "bcX", 1,
           data.frame(gene_id = "g1", gene_name = "ALPHA", stringsAsFactors = FALSE))
  # Same two columns, reversed order.
  r2 <- mk("g2", "bcX", 1,
           data.frame(gene_name = "BETA", gene_id = "g2", stringsAsFactors = FALSE))

  out <- .combine_one_core(r1, r2, "S", gene_key_col = "gene_id",
                           cell_id_col = "cell_id", metadata_conflict = "warn")
  by_id <- setNames(out$genes$gene_name, out$genes$gene_id)
  expect_equal(unname(by_id["g1"]), "ALPHA")
  expect_equal(unname(by_id["g2"]), "BETA")

  # Round 2 missing a column entirely -> NA for its genes, not a recycled value.
  r2_narrow <- mk("g2", "bcX", 1,
                  data.frame(gene_id = "g2", stringsAsFactors = FALSE))
  out2 <- .combine_one_core(r1, r2_narrow, "S", gene_key_col = "gene_id",
                            cell_id_col = "cell_id", metadata_conflict = "warn")
  expect_true(is.na(out2$genes$gene_name[out2$genes$gene_id == "g2"]))
})

test_that(".reindex errors when a key is missing from the target key set", {
  m <- Matrix::Matrix(c(1, 2), nrow = 2,
                      dimnames = list(c("g1", "g2"), "c1"), sparse = TRUE)
  expect_error(
    .reindex(m, old_keys = c("g1", "g2"), new_keys = "g1", margin = "rows"),
    "absent from"
  )
})

test_that(".combine_one_core warns when a gene key has duplicates within a round", {
  r1 <- .make_core_round(genes = c("g1", "g1"), barcodes = "bcX",
                         counts_mat = matrix(c(1, 2), nrow = 2))
  r2 <- .make_core_round(genes = "g1", barcodes = "bcX",
                         counts_mat = matrix(1, nrow = 1))

  expect_warning(
    .combine_one_core(r1, r2, "S", gene_key_col = "gene_id",
                      cell_id_col = "cell_id", metadata_conflict = "warn"),
    "duplicate values"
  )
})


# ============================================================================
# CombineParseRounds() / CombineCellRangerRounds() argument validation --
# all of it runs before any round is actually read, same pattern as
# MakeParseObj()'s validation tests in test-object-construction.R.
# ============================================================================

test_that("CombineParseRounds validates round1_paths/round2_paths/sample_names/output_dir", {
  d1 <- tempfile("round1_"); dir.create(d1)
  d2 <- tempfile("round2_"); dir.create(d2)
  on.exit(unlink(c(d1, d2), recursive = TRUE))

  expect_error(CombineParseRounds(character(0), character(0), output_dir = tempfile()),
              "non-empty character vector")
  expect_error(CombineParseRounds(c(d1, d2), d1, output_dir = tempfile()),
              "same length")
  expect_error(CombineParseRounds("/no/such/round1", d1, output_dir = tempfile()),
              "do not exist")
  expect_error(CombineParseRounds(d1, "/no/such/round2", output_dir = tempfile()),
              "do not exist")
  expect_error(
    CombineParseRounds(c(d1, d2), c(d1, d2), sample_names = "only_one",
                      output_dir = tempfile()),
    "same length"
  )
  expect_error(
    CombineParseRounds(c(d1, d2), c(d1, d2), sample_names = c("S", "S"),
                      output_dir = tempfile()),
    "unique"
  )
  expect_error(CombineParseRounds(d1, d2), "single path")
})

test_that("CombineCellRangerRounds validates round1_paths/round2_paths/sample_names/output_dir", {
  d1 <- tempfile("round1_"); dir.create(d1)
  d2 <- tempfile("round2_"); dir.create(d2)
  on.exit(unlink(c(d1, d2), recursive = TRUE))

  expect_error(CombineCellRangerRounds(character(0), character(0), output_dir = tempfile()),
              "non-empty character vector")
  expect_error(CombineCellRangerRounds(c(d1, d2), d1, output_dir = tempfile()),
              "same length")
  expect_error(CombineCellRangerRounds("/no/such/round1", d1, output_dir = tempfile()),
              "do not exist")
  expect_error(CombineCellRangerRounds(d1, "/no/such/round2", output_dir = tempfile()),
              "do not exist")
  expect_error(
    CombineCellRangerRounds(c(d1, d2), c(d1, d2), sample_names = "only_one",
                            output_dir = tempfile()),
    "same length"
  )
  expect_error(
    CombineCellRangerRounds(c(d1, d2), c(d1, d2), sample_names = c("S", "S"),
                            output_dir = tempfile()),
    "unique"
  )
  expect_error(CombineCellRangerRounds(d1, d2), "single path")
})


# ============================================================================
# Full read -> combine -> write round trip, Parse format. Genuine smoke test
# (real file I/O), written defensively per this suite's existing precedent
# for such tests (see .make_10x_dir()'s callers in test-object-construction.R).
# ============================================================================

# counts_mat supplied genes x cells; Parse's on-disk orientation is cells x
# genes, so this transposes before writing (mirroring what
# CombineParseRounds()'s own writer does in reverse on the way out).
.make_parse_round_dir <- function(genes, barcodes, counts_mat, species) {
  dir <- tempfile("parse_round_")
  dge_dir <- file.path(dir, "DGE_filtered")
  dir.create(dge_dir, recursive = TRUE)

  m <- methods::as(t(counts_mat), "CsparseMatrix")
  Matrix::writeMM(m, file.path(dge_dir, "count_matrix.mtx"))

  utils::write.csv(
    data.frame(gene_id = genes, gene_name = genes, stringsAsFactors = FALSE),
    file.path(dge_dir, "all_genes.csv"), row.names = FALSE
  )
  utils::write.csv(
    data.frame(bc_wells = barcodes, species = species, stringsAsFactors = FALSE),
    file.path(dge_dir, "cell_metadata.csv"), row.names = FALSE
  )
  dir
}

test_that("CombineParseRounds combines two rounds end-to-end and writes readable output", {
  testthat::skip_on_cran()
  r1_dir <- .make_parse_round_dir(
    genes = c("g1", "g2", "g3"), barcodes = c("bcA", "bcB", "bcC"),
    counts_mat = matrix(c(1, 0, 5,   2, 1, 0,   0, 3, 0), nrow = 3, byrow = TRUE),
    species = c("human", "human", "human")
  )
  r2_dir <- .make_parse_round_dir(
    genes = c("g2", "g3", "g4"), barcodes = c("bcB", "bcC", "bcD"),
    counts_mat = matrix(c(10, 0, 1,   0, 2, 0,   1, 1, 1), nrow = 3, byrow = TRUE),
    species = c("human", "mouse", "mouse")
  )
  out_dir <- tempfile("parse_combined_")
  on.exit(unlink(c(r1_dir, r2_dir, out_dir), recursive = TRUE))

  # NB: called directly, NOT wrapped in tryCatch + skip_if. The only inputs
  # here are a temp dir and base/Matrix I/O -- there is no optional dependency
  # or environment quirk that could legitimately make this un-runnable, so
  # anything that goes wrong is a real defect and must fail rather than skip.
  result <- suppressWarnings(CombineParseRounds(
    round1_paths = r1_dir, round2_paths = r2_dir,
    sample_names = "SampleA", output_dir = out_dir
  ))

  expect_equal(result$summary$n_cells_shared, 2)
  expect_equal(result$summary$n_cells_combined, 4)
  expect_equal(result$summary$n_genes_combined, 4)

  dge_out <- file.path(result$paths[1], "DGE_filtered")
  expect_true(file.exists(file.path(dge_out, "count_matrix.mtx")))
  genes_out <- utils::read.csv(file.path(dge_out, "all_genes.csv"))
  cells_out <- utils::read.csv(file.path(dge_out, "cell_metadata.csv"))
  expect_equal(nrow(genes_out), 4)
  expect_equal(nrow(cells_out), 4)
  expect_true(all(c("g1", "g2", "g3", "g4") %in% genes_out$gene_id))
  expect_true(all(c("bcA", "bcB", "bcC", "bcD") %in% cells_out$bc_wells))

  # Round-trip the counts too, not just the table shapes -- checking only
  # nrow() would pass even if the summed values written to disk were wrong.
  # Parse's on-disk orientation is cells x genes, so transpose back.
  back <- Matrix::t(Matrix::readMM(file.path(dge_out, "count_matrix.mtx")))
  dimnames(back) <- list(genes_out$gene_id, cells_out$bc_wells)
  expect_equal(back["g2", "bcB"], 1 + 10)   # shared gene, shared barcode: summed
  expect_equal(back["g4", "bcB"], 0 + 1)    # round-2-only gene on a shared barcode
  expect_equal(back["g1", "bcA"], 1)        # round-1-only barcode: unchanged
  expect_equal(back["g1", "bcD"], 0)        # round-1-only gene on a round-2-only barcode
})


# ============================================================================
# Full read -> combine -> write round trip, CellRanger format. Uncompressed
# files are fine here -- .read_10x_triplet() (reused by
# CombineCellRangerRounds()) matches barcodes/features/matrix names whether
# or not they're gzipped, unlike the gzip-only fixture used for
# CreateRNAObjects() in test-object-construction.R.
# ============================================================================

# counts_mat supplied genes x cells -- native 10X orientation, no transpose
# needed (unlike the Parse fixture above).
.make_cellranger_round_dir <- function(genes, barcodes, counts_mat) {
  dir <- tempfile("cr_round_")
  dir.create(dir)
  sm <- methods::as(counts_mat, "CsparseMatrix")
  Matrix::writeMM(sm, file.path(dir, "matrix.mtx"))
  writeLines(barcodes, file.path(dir, "barcodes.tsv"))
  writeLines(paste(genes, genes, "Gene Expression", sep = "\t"),
            file.path(dir, "features.tsv"))
  dir
}

test_that("CombineCellRangerRounds combines two rounds end-to-end, stripping GEM-well suffixes by default", {
  testthat::skip_on_cran()
  r1_dir <- .make_cellranger_round_dir(
    genes = c("g1", "g2", "g3"), barcodes = c("bcA-1", "bcB-1", "bcC-1"),
    counts_mat = matrix(c(1, 0, 5,   2, 1, 0,   0, 3, 0), nrow = 3, byrow = TRUE)
  )
  # Same physical cells (bcB/bcC) resequenced under a different GEM-well
  # label ("-2") -- this is exactly the case strip_barcode_suffix exists for.
  r2_dir <- .make_cellranger_round_dir(
    genes = c("g2", "g3", "g4"), barcodes = c("bcB-2", "bcC-2", "bcD-2"),
    counts_mat = matrix(c(10, 0, 1,   0, 2, 0,   1, 1, 1), nrow = 3, byrow = TRUE)
  )
  out_dir <- tempfile("cr_combined_")
  on.exit(unlink(c(r1_dir, r2_dir, out_dir), recursive = TRUE))

  result <- CombineCellRangerRounds(
    round1_paths = r1_dir, round2_paths = r2_dir,
    sample_names = "SampleA", output_dir = out_dir
  )

  expect_equal(result$summary$n_cells_shared, 2)
  expect_equal(result$summary$n_cells_combined, 4)
  expect_equal(result$summary$n_genes_combined, 4)

  mtx_out <- file.path(result$paths[1], "filtered_feature_bc_matrix")
  expect_true(file.exists(file.path(mtx_out, "matrix.mtx")))
  barcodes_out <- readLines(file.path(mtx_out, "barcodes.tsv"))
  expect_true(all(c("bcA", "bcB", "bcC", "bcD") %in% barcodes_out))
  expect_false(any(grepl("-[0-9]+$", barcodes_out)))  # suffixes stripped in output
})

test_that("CombineCellRangerRounds treats barcodes as entirely distinct across rounds when strip_barcode_suffix = FALSE", {
  testthat::skip_on_cran()
  r1_dir <- .make_cellranger_round_dir(
    genes = c("g1", "g2"), barcodes = c("bcB-1", "bcC-1"),
    counts_mat = matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
  )
  r2_dir <- .make_cellranger_round_dir(
    genes = c("g1", "g2"), barcodes = c("bcB-2", "bcC-2"),
    counts_mat = matrix(c(5, 6, 7, 8), nrow = 2, byrow = TRUE)
  )
  out_dir <- tempfile("cr_combined_nostrip_")
  on.exit(unlink(c(r1_dir, r2_dir, out_dir), recursive = TRUE))

  result <- CombineCellRangerRounds(
    round1_paths = r1_dir, round2_paths = r2_dir,
    sample_names = "SampleA", output_dir = out_dir,
    strip_barcode_suffix = FALSE
  )

  # With suffixes preserved, "bcB-1"/"bcB-2" etc never match across rounds --
  # every barcode looks round-unique, since the suffix reflects the GEM well
  # a round was processed under, not anything about the physical cell.
  expect_equal(result$summary$n_cells_shared, 0)
  expect_equal(result$summary$n_cells_round1_only, 2)
  expect_equal(result$summary$n_cells_round2_only, 2)
  expect_equal(result$summary$n_cells_combined, 4)
})

test_that("CombineCellRangerRounds errors when suffix-stripping creates duplicate barcodes within one round", {
  testthat::skip_on_cran()
  # "bcX-1" and "bcX-2" both present in the SAME round -- collide once
  # stripped, which strip_barcode_suffix = TRUE (the default) must catch
  # rather than silently letting them get merged.
  r1_dir <- .make_cellranger_round_dir(
    genes = "g1", barcodes = c("bcX-1", "bcX-2"),
    counts_mat = matrix(c(1, 2), nrow = 1)
  )
  r2_dir <- .make_cellranger_round_dir(
    genes = "g1", barcodes = "bcX-1",
    counts_mat = matrix(1, nrow = 1)
  )
  out_dir <- tempfile("cr_combined_dupe_")
  on.exit(unlink(c(r1_dir, r2_dir, out_dir), recursive = TRUE))

  expect_error(
    CombineCellRangerRounds(
      round1_paths = r1_dir, round2_paths = r2_dir,
      sample_names = "SampleA", output_dir = out_dir
    ),
    "duplicate cell IDs"
  )
})
