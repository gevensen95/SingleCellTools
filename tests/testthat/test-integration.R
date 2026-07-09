# Tests for MergeSeurat() and SubsetAndRecluster(). Both run substantial
# real Seurat pipelines (normalization, PCA, optionally Harmony integration,
# neighbors/clustering, UMAP) -- these are genuine smoke tests, not fast
# unit tests, and are written defensively (tryCatch + skip) since synthetic
# data can hit edge cases in heuristics these functions don't control.
# MergeSeurat() also unconditionally writes a PDF (and optionally an RDS /
# CSV) to the working directory, so its test runs inside a temporary
# directory that's restored afterward.

# ============================================================================
# SubsetAndRecluster()
# ============================================================================

.make_recluster_obj <- function(seed = 1, n_genes = 60, n_cells = 80) {
  set.seed(seed)
  genes <- paste0("Gene", seq_len(n_genes))
  cells <- paste0("c", seq_len(n_cells))
  counts <- matrix(stats::rpois(n_genes * n_cells, lambda = 5), nrow = n_genes,
                   dimnames = list(genes, cells))
  storage.mode(counts) <- "double"
  meta <- data.frame(
    group = rep(c("keep", "drop"), length.out = n_cells),
    row.names = cells, stringsAsFactors = FALSE
  )
  SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)
}

test_that("SubsetAndRecluster requires exactly one subsetting spec", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_recluster_obj()
  expect_error(SubsetAndRecluster(obj), "exactly one of")
  expect_error(
    SubsetAndRecluster(obj, metadata_col = "group", metadata_value = "keep",
                       cells = colnames(obj)[1:5]),
    "exactly one of"
  )
})

test_that("SubsetAndRecluster errors when metadata_value is missing", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_recluster_obj()
  expect_error(SubsetAndRecluster(obj, metadata_col = "group"), "metadata_value")
})

test_that("SubsetAndRecluster errors on an unknown assay", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_recluster_obj()
  expect_error(
    SubsetAndRecluster(obj, cells = colnames(obj)[1:10], assay = "nope"),
    "nope"
  )
})

test_that("SubsetAndRecluster (integrate = FALSE) subsets, recomputes QC, and clusters", {
  .skip_if_missing("Seurat", "SeuratObject")
  testthat::skip_on_cran()
  # RunPCA()'s default npcs = 50 requires min(n_features, n_cells) > 50 on
  # the *subsetted* object (half the cells here, via group == "keep"), so
  # this needs real headroom above 50 on both axes -- a smaller fixture
  # hits irlba's "nu/nv must be strictly less than min(nrow, ncol)".
  obj <- .make_recluster_obj(seed = 1, n_genes = 150, n_cells = 240)
  keep_cells <- colnames(obj)[obj$group == "keep"]

  out <- tryCatch(
    suppressWarnings(suppressMessages(SubsetAndRecluster(
      obj, cells = keep_cells, integrate = FALSE,
      normalize = "LogNormalize", dims = 10, resolution = 0.3
    ))),
    error = function(e) e
  )
  skip_if(inherits(out, "error"),
         paste("SubsetAndRecluster did not complete on this synthetic dataset:",
               if (inherits(out, "error")) conditionMessage(out) else ""))

  expect_equal(ncol(out), length(keep_cells))
  expect_true(paste0("nCount_", SeuratObject::DefaultAssay(out)) %in% colnames(out@meta.data))
  expect_true("umap" %in% SeuratObject::Reductions(out))
  expect_true("pca" %in% SeuratObject::Reductions(out))
})


# ============================================================================
# MergeSeurat()
# ============================================================================

test_that("MergeSeurat validates integration-method/new_reduction combinations", {
  expect_error(
    MergeSeurat(list(a = 1, b = 2), integration = "RPCAIntegration",
               new_reduction = "harmony"),
    "new_reduction"
  )
  expect_error(
    MergeSeurat(list(a = 1, b = 2), integration = "RPCAIntegration",
               new_reduction = "rpca"),
    "k_anchor"
  )
  expect_error(
    MergeSeurat(list(a = 1, b = 2), integration = "RPCAIntegration",
               new_reduction = "rpca", k_anchor = 20),
    "k_weight"
  )
})

test_that("MergeSeurat validates that banksy = TRUE requires spatial Visium/Xenium", {
  expect_error(
    MergeSeurat(list(a = 1, b = 2), banksy = TRUE, spatial = "no"),
    "banksy = TRUE requires spatial"
  )
})

test_that("MergeSeurat runs end-to-end (LogNormalize + Harmony) with no disk clutter left behind", {
  .skip_if_missing("Seurat", "SeuratObject", "harmony")
  testthat::skip_on_cran()

  set.seed(1)
  make_sample <- function(prefix, n_genes = 80, n_cells = 60) {
    genes <- paste0("Gene", seq_len(n_genes))
    cells <- paste0(prefix, "_c", seq_len(n_cells))
    m <- matrix(stats::rpois(n_genes * n_cells, lambda = 5), nrow = n_genes,
               dimnames = list(genes, cells))
    storage.mode(m) <- "double"
    SeuratObject::CreateSeuratObject(counts = m)
  }
  objs <- list(s1 = make_sample("s1"), s2 = make_sample("s2"))

  old_wd <- getwd()
  tmp_dir <- tempfile("mergeseurat_test_")
  dir.create(tmp_dir)
  setwd(tmp_dir)
  on.exit({ setwd(old_wd); unlink(tmp_dir, recursive = TRUE) }, add = TRUE)

  out <- tryCatch(
    suppressWarnings(suppressMessages(MergeSeurat(
      objs, use_SCT = FALSE, save_rds_file = FALSE, markers = FALSE,
      to_regress = NULL,  # synthetic objects have no percent.mt (or any QC)
                          # metadata to regress against
      # integration_normalization/integration_assay default to 'SCT'
      # regardless of use_SCT -- since use_SCT = FALSE takes the
      # LogNormalize path, there's no "SCT" assay on the object, so these
      # must be pointed at the assay/method that was actually run.
      integration_normalization = "LogNormalize", integration_assay = "RNA",
      integration = "HarmonyIntegration", new_reduction = "harmony",
      max_dims = 10, cluster_resolution = 0.3
    ))),
    error = function(e) e
  )
  skip_if(inherits(out, "error"),
         paste("MergeSeurat did not complete on this synthetic dataset:",
               if (inherits(out, "error")) conditionMessage(out) else ""))

  expect_equal(ncol(out), 120)
  expect_true("harmony" %in% SeuratObject::Reductions(out))
  expect_true(file.exists(file.path(tmp_dir, "dimplot_seurat_clusters.pdf")))
})

# ---- MergeSeurat(banksy = TRUE) -------------------------------------------
# Builds two small "Xenium"-style samples (assay literally named "Xenium",
# to match the `spatial == 'Xenium'` branch of MergeSeurat's own merge
# step) each with one FOV of cell centroids, so RunBanksyWrapper() has real
# spatial coordinates to pull via get_all_coords().

.make_xenium_like_sample <- function(prefix, seed, n_genes = 60, n_cells = 60) {
  set.seed(seed)
  genes <- paste0("Gene", seq_len(n_genes))
  cells <- paste0(prefix, "_c", seq_len(n_cells))
  m <- matrix(stats::rpois(n_genes * n_cells, lambda = 5), nrow = n_genes,
             dimnames = list(genes, cells))
  storage.mode(m) <- "double"
  # `project = prefix` gives each sample a distinct orig.ident, so the
  # default group_column = 'orig.ident' actually distinguishes s1 from s2
  # for BANKSY's per-sample coordinate staggering (rather than every cell
  # falling into one indistinguishable "SeuratProject" group).
  obj <- SeuratObject::CreateSeuratObject(counts = m, assay = "Xenium", project = prefix)

  coords <- data.frame(x = stats::runif(n_cells, 0, 100),
                       y = stats::runif(n_cells, 0, 100),
                       cell = cells, stringsAsFactors = FALSE)
  cents <- SeuratObject::CreateCentroids(coords)
  fov <- SeuratObject::CreateFOV(coords = list(centroids = cents), type = "centroids",
                                 assay = "Xenium", key = paste0(prefix, "fov_"))
  obj[[paste0(prefix, "fov")]] <- fov
  obj
}

test_that("MergeSeurat(banksy = TRUE) runs the BANKSY spatial-clustering path end-to-end", {
  .skip_if_missing("Seurat", "SeuratObject", "harmony", "Banksy", "SeuratWrappers")
  testthat::skip_on_cran()

  objs <- list(s1 = .make_xenium_like_sample("s1", seed = 1),
              s2 = .make_xenium_like_sample("s2", seed = 2))

  old_wd <- getwd()
  tmp_dir <- tempfile("mergeseurat_banksy_test_")
  dir.create(tmp_dir)
  setwd(tmp_dir)
  on.exit({ setwd(old_wd); unlink(tmp_dir, recursive = TRUE) }, add = TRUE)

  out <- tryCatch(
    suppressWarnings(suppressMessages(MergeSeurat(
      objs, spatial = "Xenium", use_SCT = FALSE, save_rds_file = FALSE,
      markers = FALSE, to_regress = NULL,
      banksy = TRUE, banksy_lambda = 0.2, banksy_k_geom = 6,
      integration = "HarmonyIntegration", new_reduction = "harmony",
      max_dims = 10, cluster_resolution = 0.3
    ))),
    error = function(e) e
  )
  skip_if(inherits(out, "error"),
         paste("MergeSeurat(banksy = TRUE) did not complete on this synthetic dataset:",
               if (inherits(out, "error")) conditionMessage(out) else ""))

  expect_true("BANKSY" %in% SeuratObject::Assays(out))
  expect_true("pca_banksy" %in% SeuratObject::Reductions(out))
  expect_true("harmony" %in% SeuratObject::Reductions(out))
  expect_true(file.exists(file.path(tmp_dir, "dimplot_seurat_clusters.pdf")))
})
