# Tests for the small, dependency-light utility functions: get_all_children
# (GO.db), parse_polygons, combine_metadata, detect_gene_id_type /
# check_gene_ids_across_objects, check_duplicate_genes,
# strip_workflow_artifacts, the .cluster_mat internal helper (shared by
# MarkerPlot/MarkerPctPlot), and a .onAttach smoke test.
#
# Note: DetectGeneIDs.R and DetectGenes.R both define `detect_gene_id_type`
# and `check_gene_ids_across_objects`. Files load alphabetically, so
# DetectGenes.R's versions (which add a `verbose` argument) are the ones
# that actually end up active; DetectGeneIDs.R's definitions are silently
# shadowed dead code. Tests below exercise the effective (verbose-capable)
# behavior.

# ============================================================================
# get_all_children() -- GO term recursion
# ============================================================================

test_that("get_all_children returns NULL for a term with no children", {
  testthat::skip_if_not_installed("GO.db")
  # A syntactically valid but non-existent GO id: mget(..., ifnotfound = NA)
  # resolves to NA for it, and the function returns NULL when every direct
  # child lookup is NA -- true regardless of the installed GO.db content/
  # version, so this doesn't depend on knowing real GO hierarchy details.
  result <- get_all_children("GO:9999999")
  expect_null(result)
})

test_that("get_all_children verbose = TRUE emits a message, verbose = FALSE is silent", {
  testthat::skip_if_not_installed("GO.db")
  expect_message(get_all_children("GO:9999999", verbose = TRUE), "Collecting all child GO terms")
  expect_silent(get_all_children("GO:9999999", verbose = FALSE))
})


# ============================================================================
# parse_polygons() -- WKT-style POLYGON string parsing
# ============================================================================

test_that("parse_polygons parses POLYGON strings into per-row x/y data frames", {
  df <- data.frame(
    EntityID = c("cellA", "cellB"),
    geometry = c("POLYGON ((1 2, 3 4, 5 6))",
                "POLYGON ((0 0, 1 0, 1 1, 0 1))"),
    stringsAsFactors = FALSE
  )
  result <- parse_polygons(df)
  expect_type(result, "list")
  expect_named(result, c("cellA", "cellB"))
  expect_equal(result$cellA$x, c(1, 3, 5))
  expect_equal(result$cellA$y, c(2, 4, 6))
  expect_equal(nrow(result$cellB), 4)
  expect_equal(result$cellB$x, c(0, 1, 1, 0))
})

test_that("parse_polygons leaves the list unnamed when EntityID is absent", {
  df <- data.frame(geometry = "POLYGON ((1 2, 3 4, 5 6))", stringsAsFactors = FALSE)
  result <- parse_polygons(df)
  expect_length(result, 1)
  expect_null(names(result))
})

test_that("parse_polygons errors clearly on malformed geometry instead of silently producing garbage", {
  df <- data.frame(
    EntityID = c("cellA", "cellB"),
    geometry = c("POLYGON ((1 2, 3 4, 5 6))", "not a polygon at all"),
    stringsAsFactors = FALSE
  )
  expect_error(parse_polygons(df), "not matching the expected")
})


# ============================================================================
# .cluster_mat() -- internal hierarchical-clustering helper
# ============================================================================

test_that(".cluster_mat validates the clustering method", {
  mat <- matrix(rnorm(20), nrow = 5)
  expect_error(SingleCellTools:::.cluster_mat(mat, method = "not_a_method"),
              "clustering method")
})

test_that(".cluster_mat validates the distance argument", {
  mat <- matrix(rnorm(20), nrow = 5)
  expect_error(SingleCellTools:::.cluster_mat(mat, distance = "not_a_distance"),
              "distance")
})

test_that(".cluster_mat returns an hclust object using correlation distance by default", {
  set.seed(1)
  mat <- matrix(rnorm(50), nrow = 5)
  hc <- SingleCellTools:::.cluster_mat(mat)
  expect_s3_class(hc, "hclust")
  expect_equal(length(hc$order), nrow(mat))
})

test_that(".cluster_mat tolerates a zero-variance row (undefined correlation -> 0)", {
  mat <- rbind(rnorm(20), rnorm(20), rep(5, 20))  # row 3 is constant
  # stats::cor() expectedly warns "the standard deviation is zero" for the
  # constant row; the function handles this deliberately by coercing the
  # resulting NA correlations to 0 before clustering, so the warning is
  # expected noise here, not a real problem.
  expect_no_error(hc <- suppressWarnings(SingleCellTools:::.cluster_mat(mat)))
  expect_s3_class(hc, "hclust")
})

test_that(".cluster_mat accepts a precomputed dist object", {
  set.seed(2)
  mat <- matrix(rnorm(30), nrow = 6)
  d <- stats::dist(mat)
  hc <- SingleCellTools:::.cluster_mat(mat, distance = d)
  expect_s3_class(hc, "hclust")
})


# ============================================================================
# combine_metadata() -- stack meta.data across a list of Seurat objects
# ============================================================================

.make_group_obj <- function(n, seed) {
  set.seed(seed)
  m <- matrix(stats::rpois(5 * n, 3), nrow = 5,
             dimnames = list(paste0("G", 1:5), paste0("c", seed, "_", seq_len(n))))
  storage.mode(m) <- "double"
  SeuratObject::CreateSeuratObject(
    counts = m,
    meta.data = data.frame(group = sample(c("x", "y"), n, replace = TRUE),
                           row.names = colnames(m), stringsAsFactors = FALSE)
  )
}

test_that("combine_metadata stacks meta.data with a sample-id and cell-id column", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj1 <- .make_group_obj(10, 1)
  obj2 <- .make_group_obj(15, 2)
  result <- combine_metadata(list(a = obj1, b = obj2))
  expect_equal(nrow(result), 25)
  expect_true(all(c("sample", "cell_id", "group") %in% colnames(result)))
  expect_setequal(unique(result$sample), c("a", "b"))
})

test_that("combine_metadata accepts a single Seurat object (not a list)", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_group_obj(10, 1)
  result <- combine_metadata(obj)
  expect_equal(nrow(result), 10)
})

test_that("combine_metadata auto-names an unnamed list as obj_1, obj_2, ...", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj1 <- .make_group_obj(5, 1)
  obj2 <- .make_group_obj(5, 2)
  result <- combine_metadata(list(obj1, obj2))
  expect_setequal(unique(result$sample), c("obj_1", "obj_2"))
})

test_that("combine_metadata tolerates objects with different metadata columns", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj1 <- .make_group_obj(5, 1)
  obj2 <- .make_group_obj(5, 2)
  obj2$extra_col <- "x"
  result <- combine_metadata(list(a = obj1, b = obj2))
  expect_true("extra_col" %in% colnames(result))
  expect_true(all(is.na(result$extra_col[result$sample == "a"])))
})


# ============================================================================
# detect_gene_id_type() / check_gene_ids_across_objects()
# ============================================================================

test_that("detect_gene_id_type recognizes gene symbols", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 10, n_cells = 5, gene_prefix = "Actb")
  res <- detect_gene_id_type(obj, verbose = FALSE)
  expect_equal(res$guess, "symbol")
  expect_length(res$examples, 5)
})

test_that("detect_gene_id_type recognizes Ensembl human gene IDs", {
  .skip_if_missing("Seurat", "SeuratObject")
  genes <- sprintf("ENSG%011d", seq_len(10))
  m <- matrix(stats::rpois(10 * 5, 3), nrow = 10, dimnames = list(genes, paste0("c", 1:5)))
  storage.mode(m) <- "double"
  obj <- SeuratObject::CreateSeuratObject(counts = methods::as(m, "CsparseMatrix"))
  res <- detect_gene_id_type(obj, verbose = FALSE)
  expect_equal(res$guess, "ensembl")
})

test_that("detect_gene_id_type recognizes Entrez ids", {
  .skip_if_missing("Seurat", "SeuratObject")
  genes <- as.character(1000:1009)
  m <- matrix(stats::rpois(10 * 5, 3), nrow = 10, dimnames = list(genes, paste0("c", 1:5)))
  storage.mode(m) <- "double"
  obj <- SeuratObject::CreateSeuratObject(counts = methods::as(m, "CsparseMatrix"))
  res <- detect_gene_id_type(obj, verbose = FALSE)
  expect_equal(res$guess, "entrez")
})

test_that("detect_gene_id_type verbose = TRUE emits progress messages", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 5, n_cells = 5, gene_prefix = "Actb")
  expect_message(detect_gene_id_type(obj, verbose = TRUE), "Detecting gene ID type")
})

test_that("check_gene_ids_across_objects summarizes multiple objects with mismatched ID types", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj_sym <- .make_small_seurat(seed = 1, n_genes = 10, n_cells = 5, gene_prefix = "Actb")
  genes_ens <- sprintf("ENSG%011d", seq_len(10))
  m <- matrix(stats::rpois(10 * 5, 3), nrow = 10, dimnames = list(genes_ens, paste0("c", 1:5)))
  storage.mode(m) <- "double"
  obj_ens <- SeuratObject::CreateSeuratObject(counts = methods::as(m, "CsparseMatrix"))

  res <- check_gene_ids_across_objects(list(rna = obj_sym, ref = obj_ens), verbose = FALSE)
  expect_equal(nrow(res$summary), 2)
  expect_setequal(res$summary$guess, c("symbol", "ensembl"))
  expect_setequal(res$summary$object, c("rna", "ref"))
})


# ============================================================================
# check_duplicate_genes()
# ============================================================================

test_that("check_duplicate_genes reports zero duplicates for a clean object", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 10, n_cells = 5)
  out <- capture.output(result <- check_duplicate_genes(obj))
  expect_length(result[[1]], 0)
})

test_that("check_duplicate_genes accepts a single Seurat object (not a list)", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 5, n_cells = 5)
  out <- capture.output(result <- check_duplicate_genes(obj))
  expect_named(result, "obj_1")
})

test_that("check_duplicate_genes reports duplicated rownames when present", {
  .skip_if_missing("Seurat", "SeuratObject")
  genes <- c("G1", "G2", "G2", "G3")
  m <- matrix(stats::rpois(4 * 5, 3), nrow = 4, dimnames = list(genes, paste0("c", 1:5)))
  storage.mode(m) <- "double"
  # This Seurat/SeuratObject version validates rownames at construction
  # time and throws a hard error on duplicates ("Duplicate rownames not
  # allowed") rather than silently de-duplicating them -- so there's no way
  # to get a real Seurat object with duplicate rownames via the normal
  # constructor here. Skip gracefully rather than erroring the test suite;
  # check_duplicate_genes() is still useful for objects that pick up
  # duplicates via other paths (merges, subsets, hand-built assays).
  obj <- tryCatch(SeuratObject::CreateSeuratObject(counts = methods::as(m, "CsparseMatrix")),
                 error = function(e) NULL)
  testthat::skip_if(
    is.null(obj) || !any(duplicated(rownames(obj[["RNA"]]))),
    "This Seurat/SeuratObject version rejects or de-duplicates duplicate rownames on object creation"
  )
  out <- capture.output(result <- check_duplicate_genes(obj))
  expect_true("G2" %in% result[[1]])
})


# ============================================================================
# strip_workflow_artifacts()
# ============================================================================

test_that("strip_workflow_artifacts drops the data layer, variable features, and reductions", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 30)
  SeuratObject::VariableFeatures(obj) <- rownames(obj)[1:5]
  emb <- matrix(rnorm(30 * 2), nrow = 30, dimnames = list(colnames(obj), c("PC_1", "PC_2")))
  obj[["pca"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "PC_", assay = "RNA")

  out <- strip_workflow_artifacts(obj, assay = "RNA")
  expect_length(SeuratObject::Reductions(out), 0)
  expect_length(SeuratObject::VariableFeatures(out), 0)
  expect_false("data" %in% SeuratObject::Layers(out[["RNA"]]))
})

test_that("strip_workflow_artifacts keep_reductions preserves only the named reduction(s)", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_genes = 20, n_cells = 30)
  emb1 <- matrix(rnorm(60), nrow = 30, dimnames = list(colnames(obj), c("PC_1", "PC_2")))
  emb2 <- matrix(rnorm(60), nrow = 30, dimnames = list(colnames(obj), c("harmony_1", "harmony_2")))
  obj[["pca"]]     <- SeuratObject::CreateDimReducObject(embeddings = emb1, key = "PC_", assay = "RNA")
  obj[["harmony"]] <- SeuratObject::CreateDimReducObject(embeddings = emb2, key = "harmony_", assay = "RNA")

  out <- strip_workflow_artifacts(obj, keep_reductions = "harmony")
  expect_setequal(SeuratObject::Reductions(out), "harmony")
})

test_that("strip_workflow_artifacts works on a list of Seurat objects and preserves the list shape", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj1 <- .make_small_seurat(seed = 1, n_cells = 10)
  obj2 <- .make_small_seurat(seed = 2, n_cells = 10)
  out <- strip_workflow_artifacts(list(a = obj1, b = obj2))
  expect_type(out, "list")
  expect_named(out, c("a", "b"))
  expect_true(all(vapply(out, inherits, logical(1), "Seurat")))
})

test_that("strip_workflow_artifacts errors on non-Seurat input", {
  expect_error(strip_workflow_artifacts(list(1, 2, 3)), "Seurat object")
})


# ============================================================================
# .onAttach() -- package startup hook (smoke test only)
# ============================================================================

test_that(".onAttach can be invoked directly without erroring", {
  # .onAttach normally runs once via library(SingleCellTools). This is a
  # smoke test only -- it doesn't assert which dependencies get attached,
  # since that depends on what's installed on the machine running the tests.
  expect_no_error(
    suppressPackageStartupMessages(
      SingleCellTools:::.onAttach(libname = NULL, pkgname = "SingleCellTools")
    )
  )
})


# ============================================================================
# .resolve_workers() -- shared workers-vs-cores validation helper, used by
# CreateRNAObjects/CreateVisiumObjects/CreateAndIntegrateRNA/
# CreateATACObjects/CreateATACObjectsFilter/MakeParseObj, all of which
# default `workers` to their sample count. parallel::detectCores() is mocked
# throughout so these are deterministic regardless of the machine actually
# running the tests.
# ============================================================================

test_that(".resolve_workers returns workers unchanged when within available cores", {
  testthat::local_mocked_bindings(detectCores = function(...) 8L, .package = "parallel")
  expect_equal(.resolve_workers(4, n_samples = 4, was_default = TRUE), 4)
  expect_equal(.resolve_workers(8, n_samples = 8, was_default = TRUE), 8)
})

test_that(".resolve_workers errors and names a usable core count when workers exceeds detectCores()", {
  testthat::local_mocked_bindings(detectCores = function(...) 4L, .package = "parallel")
  expect_error(.resolve_workers(12, n_samples = 12, was_default = TRUE),
              "workers = 12 exceeds the 4 core")
  expect_error(.resolve_workers(12, n_samples = 12, was_default = TRUE),
              "workers = 4")
})

test_that(".resolve_workers's error message reflects whether workers was defaulted or explicit", {
  testthat::local_mocked_bindings(detectCores = function(...) 2L, .package = "parallel")
  expect_error(.resolve_workers(6, n_samples = 6, was_default = TRUE),
              "defaulted to the sample count, 6")

  err <- tryCatch(.resolve_workers(6, n_samples = 6, was_default = FALSE),
                  error = function(e) conditionMessage(e))
  expect_false(grepl("defaulted to the sample count", err, fixed = TRUE))
})

test_that(".resolve_workers skips validation when detectCores() returns NA", {
  # Mirrors the same detectCores()-returns-NA edge case RunRCTD.R guards
  # against (minimal/HPC shells missing wc/nproc) -- can't validate against
  # an unknown core count, so this should pass `workers` through untouched
  # rather than erroring.
  testthat::local_mocked_bindings(detectCores = function(...) NA_integer_, .package = "parallel")
  expect_equal(.resolve_workers(999, n_samples = 999, was_default = TRUE), 999)
})


# ============================================================================
# workers default = sample count, across the six multi-sample loader
# functions -- checked at the formals level rather than by actually running
# each function, since most need real 10x/Visium/ATAC/Parse directories.
# ============================================================================

test_that("workers defaults to the sample-count expression in all six loader functions", {
  expect_equal(deparse(formals(CreateRNAObjects)$workers), "length(data_dirs)")
  expect_equal(deparse(formals(CreateVisiumObjects)$workers), "length(data_dirs)")
  expect_equal(deparse(formals(CreateAndIntegrateRNA)$workers), "length(data_dirs)")
  expect_equal(deparse(formals(CreateATACObjects)$workers), "length(data_dirs)")
  expect_equal(deparse(formals(CreateATACObjectsFilter)$workers), "length(data_dirs)")
  expect_equal(deparse(formals(MakeParseObj)$workers), "length(paths)")
})
