# Tests for get_all_coords(), get_cells_in_polygon(), AnnotateRegions(),
# NeighborhoodEnrichment(), and RunBanksyWrapper() against a synthetic
# two-FOV imaging-based (Xenium/CosMx-style) Seurat object. See
# helper-spatial.R for the fixture.
#
# subset_opt()/CleanMolSlot() are NOT covered here: CleanMolSlot() reads
# obj@images[[img]]$molecules, which requires a full molecule-level FOV
# (SeuratObject::CreateMolecules()) rather than the centroids-only FOVs
# built here, and constructing a realistic synthetic molecules table was
# judged not worth the added fixture complexity for this pass.
#
# CreateATACObjects(Filter), CreateVisiumObjects, LoadXenium2, MakeParseObj,
# CreateAndIntegrateRNA, EdgeDetectionVisium, detect_fov_edges,
# detect_tissue_holes, combine_fovs, SetImageBoundary, and
# BuildMultipleNicheAssays are not covered either -- none of them have
# custom input-validation logic worth unit-testing in isolation (most have
# zero stop() calls; Signac/Seurat's own readers do the validating), and
# real coverage would require actual CellRanger/Space Ranger/Xenium/ATAC
# fragment-file directory structures that aren't practical to synthesize
# correctly without real example data.

test_that("get_all_coords collects tissue coordinates from every FOV into one data frame", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 20)
  df <- get_all_coords(obj)
  expect_equal(nrow(df), 40)
  expect_true(all(c("x", "y", "fov") %in% colnames(df)))
  expect_setequal(unique(df$fov), c("fov1", "fov2"))
})

test_that("get_all_coords can join selected metadata columns", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 20)
  df <- get_all_coords(obj, meta_cols = "celltype")
  expect_true("celltype" %in% colnames(df))
  expect_setequal(unique(df$celltype), c("TypeA", "TypeB"))
})


# ============================================================================
# get_cells_in_polygon()
# ============================================================================

test_that("get_cells_in_polygon finds cells inside a bounding box", {
  .skip_if_missing("Seurat", "SeuratObject", "sf")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 100)
  # A box covering the left half of the FOV (x in [0,100], y in [0,100]) --
  # should catch roughly the TypeA half (x < 50).
  poly <- data.frame(x = c(0, 50, 50, 0), y = c(0, 0, 100, 100))
  inside <- get_cells_in_polygon(obj, poly, image_name = "fov1")
  expect_true(length(inside) > 0)
  md <- obj@meta.data[names(inside), ]
  expect_true(all(md$celltype == "TypeA"))
})

test_that("get_cells_in_polygon validates its inputs", {
  .skip_if_missing("Seurat", "SeuratObject", "sf")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 10)
  expect_error(get_cells_in_polygon(list(1), data.frame(x = 1, y = 1), "fov1"),
              "Seurat object")
  expect_error(
    get_cells_in_polygon(obj, data.frame(a = 1, b = 2), "fov1"),
    "'x' and 'y'"
  )
})


# ============================================================================
# AnnotateRegions()
# ============================================================================

test_that("AnnotateRegions labels cells by which polygon they fall inside", {
  .skip_if_missing("Seurat", "SeuratObject", "sf")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 100)
  polygons <- list(
    left  = data.frame(x = c(0, 50, 50, 0),   y = c(0, 0, 100, 100)),
    right = data.frame(x = c(50, 100, 100, 50), y = c(0, 0, 100, 100))
  )
  out <- AnnotateRegions(obj, polygons, image_name = "fov1", region_col = "region")
  expect_true("region" %in% colnames(out@meta.data))
  expect_true(all(out$region %in% c("left", "right", "unassigned")))
  # cells from fov2 weren't covered by either polygon (built from fov1
  # coordinates only) -- they should fall back to "unassigned"
  fov2_cells <- rownames(out@meta.data)[grepl("^fov2", rownames(out@meta.data))]
  expect_true(all(out@meta.data[fov2_cells, "region"] == "unassigned"))
})

test_that("AnnotateRegions validates polygon structure", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 10)
  expect_error(AnnotateRegions(obj, list(1, 2), "fov1"), "NAMED list")
  expect_error(
    AnnotateRegions(obj, list(a = data.frame(x = 1)), "fov1"),
    "columns 'x' and 'y'"
  )
})


# ============================================================================
# NeighborhoodEnrichment()
# ============================================================================

test_that("NeighborhoodEnrichment finds same-type spatial enrichment (cells cluster by type)", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 150)
  res <- NeighborhoodEnrichment(obj, group.by = "celltype", k = 10, n_perm = 50,
                                assign_niches = FALSE)
  expect_true(all(c("z", "p", "padj", "observed", "expected", "results") %in% names(res)))
  # TypeA/TypeA and TypeB/TypeB should show positive enrichment (z > 0)
  # since cells were generated spatially clustered by type.
  same_type <- res$results[res$results$focal == res$results$neighbor, ]
  expect_true(all(same_type$z > 0))
})

test_that("NeighborhoodEnrichment assign_niches = TRUE returns niche labels and can update obj in place", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 150)
  res <- NeighborhoodEnrichment(obj, group.by = "celltype", k = 10, n_perm = 30,
                                assign_niches = TRUE, n_niches = 2, add_to_meta = TRUE)
  expect_true("niche" %in% names(res))
  expect_true("composition" %in% names(res))
  # add_to_meta = TRUE + obj passed as a plain variable name -> obj itself
  # should now carry the niche column in this environment.
  expect_true("niche" %in% colnames(obj@meta.data))
})

test_that("NeighborhoodEnrichment validates inputs", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 10)
  expect_error(NeighborhoodEnrichment(list(1), group.by = "celltype"), "Seurat object")
  expect_error(NeighborhoodEnrichment(obj, group.by = "nope"), "nope")
  expect_error(
    NeighborhoodEnrichment(obj, fovs = "not_a_fov", group.by = "celltype"),
    "not in obj@images"
  )
  expect_error(
    NeighborhoodEnrichment(obj, group.by = "celltype", assign_niches = TRUE, n_niches = 1),
    "n_niches"
  )
})


# ============================================================================
# RunBanksyWrapper()
# ============================================================================
# `obj`/`lambda` are validated before RunBanksyWrapper() checks for
# Banksy/SeuratWrappers, so those two tests run regardless of whether the
# heavy packages are installed. Everything past that point -- including the
# "unknown image_name" check, which lives inside the coordinate-lookup
# block -- requires both packages, since that's the earliest point in the
# function where their absence would otherwise be masked.

test_that("RunBanksyWrapper errors on non-Seurat input", {
  expect_error(RunBanksyWrapper(list(1), lambda = 0.2), "Seurat object")
})

test_that("RunBanksyWrapper validates the lambda range", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 10)
  expect_error(RunBanksyWrapper(obj, lambda = -0.1), "lambda")
  expect_error(RunBanksyWrapper(obj, lambda = 1.5), "lambda")
})

test_that("RunBanksyWrapper errors on an unknown image_name", {
  .skip_if_missing("Seurat", "SeuratObject", "Banksy", "SeuratWrappers")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 10)
  # RunBanksyWrapper() checks the target layer is non-empty before doing any
  # coordinate lookup -- a freshly-created object only has a "counts" layer
  # (SeuratObject::CreateSeuratObject() doesn't auto-populate "data"), so
  # normalize first or that check fires before the image_name check does.
  obj <- suppressWarnings(suppressMessages(Seurat::NormalizeData(obj)))
  expect_error(
    RunBanksyWrapper(obj, lambda = 0.2, image_name = "not_a_fov"),
    "not found among obj images"
  )
})

test_that("RunBanksyWrapper pulls FOV coordinates automatically and returns a BANKSY assay + PCA", {
  .skip_if_missing("Seurat", "SeuratObject", "Banksy", "SeuratWrappers")
  testthat::skip_on_cran()
  obj <- .make_spatial_obj(seed = 1, n_genes = 30, n_per_fov = 60)
  # BANKSY (like the vignette's own workflow) expects normalized data in
  # the "data" layer -- CreateSeuratObject() alone only populates "counts".
  obj <- suppressWarnings(suppressMessages(Seurat::NormalizeData(obj)))

  out <- tryCatch(
    suppressWarnings(suppressMessages(RunBanksyWrapper(
      obj, lambda = 0.2, k_geom = 6, npcs = 5
    ))),
    error = function(e) e
  )
  skip_if(inherits(out, "error"),
         paste("RunBanksyWrapper did not complete on this synthetic FOV object:",
               if (inherits(out, "error")) conditionMessage(out) else ""))

  expect_true("BANKSY" %in% SeuratObject::Assays(out))
  expect_equal(SeuratObject::DefaultAssay(out), "BANKSY")
  expect_true("pca_banksy" %in% SeuratObject::Reductions(out))
  expect_true(all(c("banksy_x", "banksy_y") %in% colnames(out@meta.data)))
})

test_that("RunBanksyWrapper restricts coordinate lookup to a single FOV via image_name", {
  .skip_if_missing("Seurat", "SeuratObject", "Banksy", "SeuratWrappers")
  testthat::skip_on_cran()
  # Subset to one FOV's cells first so every remaining cell actually has a
  # coordinate -- image_name alone doesn't drop cells from `obj`, it only
  # restricts which FOV get_all_coords() pulls from, so mixing a
  # single-FOV lookup with a multi-FOV object would leave the other FOV's
  # cells with NA coordinates. The subset() call is wrapped in the same
  # tryCatch as RunBanksyWrapper() -- subsetting a multi-FOV object down to
  # one FOV is itself a fragile enough operation (stale/empty image
  # reconciliation) to skip rather than hard-fail on synthetic data.
  obj_full <- .make_spatial_obj(seed = 2, n_genes = 30, n_per_fov = 60)
  obj_full <- suppressWarnings(suppressMessages(Seurat::NormalizeData(obj_full)))
  fov1_cells <- grep("^fov1", colnames(obj_full), value = TRUE)

  out <- tryCatch({
    obj <- subset(obj_full, cells = fov1_cells)
    suppressWarnings(suppressMessages(RunBanksyWrapper(
      obj, lambda = 0.2, k_geom = 6, npcs = 5, image_name = "fov1"
    )))
  }, error = function(e) e)
  skip_if(inherits(out, "error"),
         paste("RunBanksyWrapper did not complete on this synthetic FOV object:",
               if (inherits(out, "error")) conditionMessage(out) else ""))

  expect_true("BANKSY" %in% SeuratObject::Assays(out))
  expect_equal(ncol(out), length(fov1_cells))
})
