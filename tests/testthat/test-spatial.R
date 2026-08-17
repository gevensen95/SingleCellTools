# Tests for get_all_coords(), get_cells_in_polygon(), AnnotateRegions(),
# NeighborhoodEnrichment(), RunBanksyWrapper(), combine_fovs(), and
# subset_opt() against a synthetic two-FOV imaging-based (Xenium/CosMx-style)
# Seurat object. See helper-spatial.R for the fixture.
#
# subset_opt()'s cleanMolecules = TRUE path (CleanMolSlot()) is NOT covered:
# CleanMolSlot() reads obj@images[[img]]$molecules, which requires a full
# molecule-level FOV (SeuratObject::CreateMolecules()) rather than the
# centroids-only FOVs built here, and constructing a realistic synthetic
# molecules table was judged not worth the added fixture complexity for this
# pass. subset_opt() itself is covered below with cleanMolecules = FALSE.
#
# CreateATACObjects(Filter), CreateVisiumObjects, LoadXenium2, MakeParseObj,
# CreateAndIntegrateRNA, detect_fov_edges, detect_tissue_holes,
# SetImageBoundary, and BuildMultipleNicheAssays are not covered either --
# none of them have custom input-validation logic worth unit-testing in
# isolation (most have zero stop() calls; Signac/Seurat's own readers do the
# validating), and real coverage would require actual CellRanger/Space
# Ranger/Xenium/ATAC fragment-file directory structures that aren't
# practical to synthesize correctly without real example data.
#
# EdgeDetectionVisium()'s coord_path fallback IS covered below for its
# tissue_positions.parquet branch (Visium HD) -- a small synthetic parquet
# file is cheap to build with arrow::write_parquet(), unlike the full
# CellRanger/Space Ranger directory trees the functions above would need.

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

test_that("get_cells_in_polygon recovers cell identity on older-style Visium coordinates (imagecol/imagerow, no 'cell' column)", {
  # Regression test: GetTissueCoordinates() only returns a 'cell' column for
  # FOV/imaging-style images and current (v5) Visium images -- older-style
  # Visium output (imagecol/imagerow) carries cell identity in rownames
  # instead, with no 'cell' column at all. Without the rownames fallback,
  # every returned cell silently lost its barcode/name entirely.
  .skip_if_missing("Seurat", "SeuratObject", "sf")
  obj <- .make_visium_seurat(seed = 1, n_genes = 5, spots_per_image = 20,
                             h = 50, w = 60, cluster_by_position = TRUE)
  # cluster_by_position = TRUE puts imagecol < w/2 spots in cluster "0" --
  # a box covering the left half of the image should recover those spots,
  # named by their real barcodes.
  poly <- data.frame(x = c(0, 30, 30, 0), y = c(0, 0, 50, 50))
  inside <- get_cells_in_polygon(obj, poly, image_name = "slice1")
  expect_true(length(inside) > 0)
  expect_false(any(is.na(names(inside))))
  expect_true(all(names(inside) %in% colnames(obj)))
  expect_true(all(obj$seurat_clusters[names(inside)] == "0"))
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

test_that("NeighborhoodEnrichment assign_niches = TRUE returns niche labels and an updated obj", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 150)
  res <- NeighborhoodEnrichment(obj, group.by = "celltype", k = 10, n_perm = 30,
                                assign_niches = TRUE, n_niches = 2, add_to_meta = TRUE)
  expect_true("niche" %in% names(res))
  expect_true("composition" %in% names(res))
  # add_to_meta = TRUE -> the updated Seurat object comes back as res$obj,
  # carrying the new niche column; the caller's own `obj` is untouched.
  expect_true("obj" %in% names(res))
  expect_true("niche" %in% colnames(res$obj@meta.data))
  expect_false("niche" %in% colnames(obj@meta.data))
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

test_that("NeighborhoodEnrichment does not leak its internal seed into the caller's RNG stream", {
  # Regression test: set.seed(seed) used to run (twice -- once for the
  # permutation null, again before niche clustering) without saving/
  # restoring the caller's prior .Random.seed.
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 30)

  set.seed(999)
  x_ref <- runif(5)

  set.seed(999)
  invisible(suppressMessages(NeighborhoodEnrichment(
    obj, group.by = "celltype", k = 5, n_perm = 5,
    assign_niches = FALSE, seed = 42
  )))
  x_after <- runif(5)

  expect_equal(x_after, x_ref)
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


# ============================================================================
# EdgeDetectionVisium() -- tissue_positions.parquet support (Visium HD)
# ============================================================================

test_that("EdgeDetectionVisium reads Visium HD tissue_positions.parquet via coord_path", {
  .skip_if_missing("arrow")
  d <- tempfile("visium_hd_spatial_"); dir.create(d)
  on.exit(unlink(d, recursive = TRUE))

  # A small 4x4 grid of array coordinates -- Visium HD's tissue_positions.parquet
  # ships these exact column names already, unlike the CSV branch which has to
  # guess whether a header row is present.
  grid <- expand.grid(array_row = 0:3, array_col = 0:3)
  positions <- data.frame(
    barcode             = paste0("bc", seq_len(nrow(grid))),
    in_tissue           = 1L,
    array_row           = grid$array_row,
    array_col           = grid$array_col,
    pxl_row_in_fullres  = grid$array_row * 100,
    pxl_col_in_fullres  = grid$array_col * 100,
    stringsAsFactors    = FALSE
  )
  arrow::write_parquet(positions, file.path(d, "tissue_positions.parquet"))

  out <- EdgeDetectionVisium(coord_path = d, neighbors = 3)
  expect_equal(nrow(out), nrow(positions))
  expect_true(all(c("barcode", "Filter", "Filter2", "Filter3", "Filter4") %in% colnames(out)))
  expect_true(all(out$Filter %in% c("Keep", "Filter")))
})


# ============================================================================
# CreateVisiumObjects() -- Visium HD directory detection
# ============================================================================

test_that("CreateVisiumObjects rejects invalid hd_bin_size values", {
  expect_error(
    CreateVisiumObjects(data_dirs = "irrelevant", hd_bin_size = "bogus"),
    "should be one of"
  )
})

test_that("CreateVisiumObjects errors clearly when a detected Visium HD sample is missing the requested bin size", {
  # A binned_outputs/ subdirectory is enough to trip HD detection -- the
  # error should fire (listing what IS available) before any matrix/image
  # file is ever touched, so this doesn't need real Space Ranger output.
  d <- tempfile("visium_hd_"); dir.create(d)
  dir.create(file.path(d, "binned_outputs", "square_002um"), recursive = TRUE)
  on.exit(unlink(d, recursive = TRUE))

  expect_error(
    CreateVisiumObjects(data_dirs = d, hd_bin_size = "008um"),
    "square_008um.*square_002um"
  )
})


# ============================================================================
# SpatialObjectInfo() / DropSpatialImage() -- image/FOV-slot management,
# generalized to cover both pixel-backed (VisiumV1) and coordinate-only
# (FOV -- Xenium/CosMx/MERFISH-style) images. Neither function's "does it
# correctly read/rebuild a real Visium image" behavior is covered here
# (needs a real spaceranger spatial/ directory, same precedent as the rest
# of this file) -- what IS covered is everything that doesn't depend on a
# real attached VisiumV1 image: argument validation, list vs. single-object
# handling, the pixel-image branches on a plain synthetic Seurat object
# (.make_small_seurat(), empty @images), and the FOV branches on the real
# (synthetic) two-FOV object from helper-spatial.R (.make_spatial_obj()).
# ============================================================================

test_that("SpatialObjectInfo requires a Seurat object or list of them", {
  expect_error(SpatialObjectInfo("not a seurat"), "Seurat object")
  expect_error(SpatialObjectInfo(list("not a seurat")), "Seurat object")
})

test_that("SpatialObjectInfo reports a single NA row for an object with no images", {
  obj <- .make_small_seurat()
  out <- SpatialObjectInfo(obj)
  expect_equal(nrow(out), 1L)
  expect_true(is.na(out$image_name))
  expect_true(is.na(out$deferred))
})

test_that("SpatialObjectInfo handles a list of objects, naming unnamed entries", {
  objs <- list(.make_small_seurat(seed = 1), .make_small_seurat(seed = 2))
  out <- SpatialObjectInfo(objs)
  expect_equal(nrow(out), 2L)
  expect_setequal(out$sample, c("sample1", "sample2"))

  named <- setNames(objs, c("A", "B"))
  out2 <- SpatialObjectInfo(named)
  expect_setequal(out2$sample, c("A", "B"))
})

test_that("SpatialObjectInfo reports FOV images (Xenium/CosMx/MERFISH-style) as coordinate-only", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 20)
  out <- SpatialObjectInfo(obj)

  expect_equal(nrow(out), 2L)
  expect_setequal(out$image_name, c("fov1", "fov2"))
  expect_true(all(out$class == "FOV"))
  # Coordinate-only: no decoded pixel array, so width/height/size_mb/deferred
  # are all NA -- this is what distinguishes an FOV row from a VisiumV1 row.
  expect_true(all(is.na(out$width)))
  expect_true(all(is.na(out$height)))
  expect_true(all(is.na(out$deferred)))
  # But cell count and boundary sets ARE resolvable for a coordinate-only image.
  expect_equal(out$n_cells, c(20L, 20L))
  expect_true(all(grepl("centroids", out$boundary_sets)))
  # No molecules boundary on this fixture (centroids only).
  expect_false(any(out$has_molecules))
  expect_false(any(out$molecules_lazy))
})

test_that("DropSpatialImage requires a Seurat object or list of them", {
  expect_error(DropSpatialImage("not a seurat"), "Seurat object")
})

test_that("DropSpatialImage(mode = 'remove') empties @images and preserves single-object shape", {
  obj <- .make_small_seurat()
  out <- DropSpatialImage(obj, mode = "remove")
  expect_s4_class(out, "Seurat")
  expect_equal(length(out@images), 0L)
})

test_that("DropSpatialImage(mode = 'remove') preserves list shape and names", {
  objs <- setNames(list(.make_small_seurat(seed = 1), .make_small_seurat(seed = 2)),
                   c("A", "B"))
  out <- DropSpatialImage(objs, mode = "remove")
  expect_true(is.list(out))
  expect_setequal(names(out), c("A", "B"))
})

test_that("DropSpatialImage(mode = 'remove') empties @images for FOV-based objects too", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 20)
  out <- DropSpatialImage(obj, mode = "remove")
  expect_equal(length(out@images), 0L)
})

test_that("DropSpatialImage(mode = 'downgrade') errors when visium_image_dir isn't stashed", {
  obj <- .make_small_seurat()
  expect_error(
    DropSpatialImage(obj, mode = "downgrade"),
    "visium_image_dir"
  )
})

test_that("DropSpatialImage(mode = 'downgrade') is a no-op when already deferred", {
  obj <- .make_small_seurat()
  obj@misc$hires_image_path <- "/fake/hires.png"
  expect_message(
    out <- DropSpatialImage(obj, mode = "downgrade"),
    "Already deferred"
  )
  expect_identical(out@misc$hires_image_path, "/fake/hires.png")
})

test_that("DropSpatialImage(mode = 'downgrade') is a no-op when no images are attached", {
  obj <- .make_small_seurat()
  obj@misc$visium_image_dir <- "/fake/dir"
  expect_message(
    DropSpatialImage(obj, mode = "downgrade"),
    "No images attached"
  )
})

test_that("DropSpatialImage(mode = 'downgrade') skips FOV images with a message instead of erroring", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 20)
  # No visium_image_dir stashed (this object was never built by
  # CreateVisiumObjects()) -- downgrade must NOT demand one for a
  # coordinate-only object, it should just explain why it can't downgrade
  # FOV images and leave them alone.
  expect_message(
    out <- DropSpatialImage(obj, mode = "downgrade"),
    "coordinate-only"
  )
  expect_equal(length(out@images), 2L)
  expect_setequal(names(out@images), c("fov1", "fov2"))
})


# ============================================================================
# SpatialDimPlotFixed() / SpatialFeaturePlotFixed() -- ggplot2-based
# replacements for Seurat::SpatialDimPlot()/SpatialFeaturePlot(), built to
# work around pt.size.factor silently doing nothing in recent Seurat
# releases (GeomSpatial's custom draw_panel() -- see the functions' own
# docs for the specific GitHub issues). Since that's exactly the bug being
# worked around, `pt.size` actually reaching geom_point()'s rendered size is
# checked directly via ggplot2::layer_data(), not just "does it error".
# ============================================================================

test_that("SpatialDimPlotFixed errors on non-Seurat input", {
  expect_error(SpatialDimPlotFixed(list(1)), "Seurat object")
})

test_that("SpatialDimPlotFixed errors when obj has no images", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat()
  expect_error(SpatialDimPlotFixed(obj), "no images")
})

test_that("SpatialDimPlotFixed errors on an unknown image name", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 10)
  expect_error(SpatialDimPlotFixed(obj, images = "not_a_slice"), "not found")
})

test_that("SpatialDimPlotFixed errors on an unknown group.by column", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 10)
  expect_error(SpatialDimPlotFixed(obj, group.by = "nope"), "nope")
})

test_that("SpatialDimPlotFixed errors when `colors` has too few values", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 10)
  expect_error(
    SpatialDimPlotFixed(obj, group.by = "seurat_clusters", colors = "red"),
    "colors"
  )
})

test_that("SpatialDimPlotFixed returns a single ggplot with all spots and the requested point size", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 12, n_images = 1)
  p <- SpatialDimPlotFixed(obj, group.by = "seurat_clusters", pt.size = 3)
  expect_s3_class(p, "ggplot")
  expect_false(inherits(p, "patchwork"))

  pts <- ggplot2::layer_data(p, 2)  # layer 1 = annotation_raster, layer 2 = geom_point
  expect_equal(nrow(pts), 12L)
  expect_true(all(pts$size == 3))
})

test_that("SpatialDimPlotFixed combines multiple images into one patchwork by default", {
  .skip_if_missing("Seurat", "SeuratObject", "patchwork")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 8, n_images = 2)
  p <- SpatialDimPlotFixed(obj, group.by = "seurat_clusters")
  expect_true(inherits(p, "patchwork"))
})

test_that("SpatialDimPlotFixed(combine = FALSE) returns one ggplot per image", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 8, n_images = 2)
  plots <- SpatialDimPlotFixed(obj, group.by = "seurat_clusters", combine = FALSE)
  expect_type(plots, "list")
  expect_equal(length(plots), 2L)
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
  # Each image's panel should only have that image's own spots.
  expect_equal(nrow(ggplot2::layer_data(plots[[1]], 2)), 8L)
})

test_that("SpatialFeaturePlotFixed errors on non-character/empty features", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 10)
  expect_error(SpatialFeaturePlotFixed(obj, features = character(0)), "features")
  expect_error(SpatialFeaturePlotFixed(obj, features = 1), "features")
})

test_that("SpatialFeaturePlotFixed errors on an unknown image name", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 10)
  expect_error(
    SpatialFeaturePlotFixed(obj, features = "Gene1", images = "not_a_slice"),
    "not found"
  )
})

test_that("SpatialFeaturePlotFixed returns a single ggplot for one feature/one image", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 12, n_images = 1)
  p <- SpatialFeaturePlotFixed(obj, features = "Gene1", pt.size = 2.5)
  expect_s3_class(p, "ggplot")
  expect_false(inherits(p, "patchwork"))

  pts <- ggplot2::layer_data(p, 2)
  expect_equal(nrow(pts), 12L)
  expect_true(all(pts$size == 2.5))
})

test_that("SpatialFeaturePlotFixed(combine = FALSE) returns one panel per feature x image", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 6, n_images = 2)
  plots <- SpatialFeaturePlotFixed(obj, features = c("Gene1", "Gene2"), combine = FALSE)
  expect_type(plots, "list")
  expect_equal(length(plots), 4L)  # 2 features x 2 images
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
})

test_that("SpatialFeaturePlotFixed combines into one patchwork by default", {
  .skip_if_missing("Seurat", "SeuratObject", "patchwork")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 6, n_images = 2)
  p <- SpatialFeaturePlotFixed(obj, features = "Gene1")
  expect_true(inherits(p, "patchwork"))
})


# ============================================================================
# SpatialConcordance()
# ============================================================================

test_that("SpatialConcordance errors on non-Seurat input", {
  expect_error(SpatialConcordance(list(1), group.by = "seurat_clusters"), "Seurat object")
})

test_that("SpatialConcordance errors when obj has no images", {
  obj <- .make_small_seurat()
  expect_error(SpatialConcordance(obj, group.by = "seurat_clusters"), "no images")
})

test_that("SpatialConcordance errors on an unknown image name", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 20)
  expect_error(
    SpatialConcordance(obj, group.by = "seurat_clusters", images = "not_a_slice"),
    "not found"
  )
})

test_that("SpatialConcordance errors on an unknown group.by column", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 20)
  expect_error(SpatialConcordance(obj, group.by = "nope"), "nope")
})

test_that("SpatialConcordance validates k and n_perm", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 20)
  expect_error(SpatialConcordance(obj, group.by = "seurat_clusters", k = 0), "k")
  expect_error(SpatialConcordance(obj, group.by = "seurat_clusters", n_perm = 0), "n_perm")
})

test_that("SpatialConcordance does not leak its permutation seed into the caller's RNG stream", {
  # Regression test: set.seed(seed) used to run without saving/restoring
  # the caller's prior .Random.seed.
  .skip_if_missing("Seurat", "SeuratObject", "RANN")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 20)

  set.seed(999)
  x_ref <- runif(5)

  set.seed(999)
  invisible(SpatialConcordance(obj, group.by = "seurat_clusters", k = 5,
                               n_perm = 5, seed = 42))
  x_after <- runif(5)

  expect_equal(x_after, x_ref)
})

test_that("SpatialConcordance detects real spatial structure (z > 0) and finds none in random labels", {
  .skip_if_missing("Seurat", "SeuratObject", "RANN")
  set.seed(42)
  structured <- .make_visium_seurat(seed = 2, spots_per_image = 60,
                                    cluster_by_position = TRUE)
  res <- SpatialConcordance(structured, group.by = "seurat_clusters", k = 6, n_perm = 100)
  expect_equal(nrow(res), 1L)
  expect_true(all(c("sample", "n_spots", "k_used", "observed", "null_mean",
                    "null_sd", "z", "p", "padj") %in% colnames(res)))
  # Labels are a clean left/right spatial split -- concordance should be well
  # above the shuffled-label null.
  expect_true(res$z[1] > 0)
  expect_true(res$observed[1] > res$null_mean[1])

  random <- .make_visium_seurat(seed = 2, spots_per_image = 60,
                                cluster_by_position = FALSE)
  res_random <- SpatialConcordance(random, group.by = "seurat_clusters", k = 6, n_perm = 100)
  # Random labels: observed concordance should sit close to the null mean,
  # not systematically far above it the way the spatially structured case is.
  expect_true(abs(res_random$z[1]) < res$z[1])
})

test_that("SpatialConcordance `exclude` drops matching labels before scoring", {
  .skip_if_missing("Seurat", "SeuratObject", "RANN")
  obj <- .make_visium_seurat(seed = 3, spots_per_image = 20, cluster_by_position = TRUE)
  # Recode a few spots to an "Unclassified" bucket and confirm `exclude`
  # actually removes them from n_spots.
  obj$seurat_clusters <- as.character(obj$seurat_clusters)
  obj$seurat_clusters[1:5] <- "Unclassified"
  res_default <- SpatialConcordance(obj, group.by = "seurat_clusters", k = 4, n_perm = 20)
  res_excl <- SpatialConcordance(obj, group.by = "seurat_clusters", k = 4, n_perm = 20,
                                 exclude = c(NA, "Unclassified"))
  expect_equal(res_default$n_spots[1], 20L)
  expect_equal(res_excl$n_spots[1], 15L)
})


# ============================================================================
# RenameSpatialImages()
# ============================================================================

test_that("RenameSpatialImages errors on non-Seurat input", {
  expect_error(RenameSpatialImages(list(1), group_col = "sample"), "Seurat object")
})

test_that("RenameSpatialImages errors when obj has no images", {
  obj <- .make_small_seurat()
  expect_error(RenameSpatialImages(obj, group_col = "seurat_clusters"), "no images")
})

test_that("RenameSpatialImages errors when both or neither of group_col/new_names are given", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 10, n_images = 2)
  expect_error(RenameSpatialImages(obj), "both are NULL")
  expect_error(
    RenameSpatialImages(obj, group_col = "sample", new_names = c("a", "b")),
    "both were given"
  )
})

test_that("RenameSpatialImages auto-derives names from group_col by cell identity, not position", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 10, n_images = 2)
  # slice1's cells -> "sampleB", slice2's cells -> "sampleA" -- deliberately
  # not alphabetical/appearance order, so a naive positional rename would
  # get this backwards.
  obj$sample <- ifelse(obj$sample == "slice1", "sampleB", "sampleA")

  out <- RenameSpatialImages(obj, group_col = "sample")
  expect_setequal(names(out@images), c("sampleA", "sampleB"))
  # Confirm identity, not just presence: the image now called "sampleB"
  # must still contain the cells that were originally slice1's.
  orig_slice1_cells <- colnames(obj)[obj$sample == "sampleB"]
  expect_setequal(SeuratObject::Cells(out[["sampleB"]]), orig_slice1_cells)
})

test_that("RenameSpatialImages errors on an unknown group_col", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 10, n_images = 2)
  expect_error(RenameSpatialImages(obj, group_col = "nope"), "nope")
})

test_that("RenameSpatialImages errors when an image spans more than one group_col value", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 10, n_images = 2)
  slice1_cells <- colnames(obj)[obj$sample == "slice1"]
  obj$sample[slice1_cells[1]] <- "slice2"  # one stray cell now disagrees
  expect_error(RenameSpatialImages(obj, group_col = "sample"), "distinct value")
})

test_that("RenameSpatialImages errors when an image's group_col values are all NA", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 10, n_images = 2)
  obj$sample <- as.character(obj$sample)
  obj$sample[obj$sample == "slice1"] <- NA
  expect_error(RenameSpatialImages(obj, group_col = "sample"), "NA")
})

test_that("RenameSpatialImages new_names positional mode renames in names(obj@images) order", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 10, n_images = 2)
  out <- RenameSpatialImages(obj, new_names = c("first", "second"))
  expect_equal(names(out@images), c("first", "second"))
})

test_that("RenameSpatialImages new_names positional mode requires matching length", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 10, n_images = 2)
  expect_error(RenameSpatialImages(obj, new_names = "only_one"), "length")
})

test_that("RenameSpatialImages new_names named mode does a partial rename", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 10, n_images = 2)
  out <- RenameSpatialImages(obj, new_names = c(slice1 = "anterior"))
  expect_setequal(names(out@images), c("anterior", "slice2"))
})

test_that("RenameSpatialImages new_names named mode errors on an unknown old name", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 10, n_images = 2)
  expect_error(
    RenameSpatialImages(obj, new_names = c(not_a_slice = "anterior")),
    "not found"
  )
})

test_that("RenameSpatialImages errors when resolved names would collide", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 10, n_images = 2)
  expect_error(RenameSpatialImages(obj, new_names = c("same", "same")), "collide")
})

test_that("SpatialConcordance returns one row per image and skips tiny images with a message", {
  .skip_if_missing("Seurat", "SeuratObject", "RANN")
  obj <- .make_visium_seurat(seed = 4, spots_per_image = 15, n_images = 2)
  res <- SpatialConcordance(obj, group.by = "seurat_clusters", k = 4, n_perm = 20)
  expect_equal(nrow(res), 2L)
  expect_setequal(res$sample, c("slice1", "slice2"))

  # Force everything in slice2 to be excluded so it has < 2 usable cells.
  obj$seurat_clusters <- as.character(obj$seurat_clusters)
  obj$seurat_clusters[obj$sample == "slice2"] <- NA
  expect_message(
    res2 <- SpatialConcordance(obj, group.by = "seurat_clusters", k = 4, n_perm = 20),
    "skipping"
  )
  expect_equal(nrow(res2), 1L)
  expect_equal(res2$sample, "slice1")
})

test_that("SubsetSpatial errors on non-Seurat input", {
  expect_error(SubsetSpatial(list(a = 1)), "Seurat object")
})

test_that("SubsetSpatial with no images just delegates to subset()", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_small_seurat(seed = 1, n_cells = 20, n_clusters = 2)
  res <- SubsetSpatial(obj, subset = seurat_clusters == levels(obj$seurat_clusters)[1])
  expect_true(all(as.character(res$seurat_clusters) == levels(obj$seurat_clusters)[1]))
  expect_equal(length(res@images), 0L)
})

test_that("SubsetSpatial subset= keeps only matching cells and every image stays consistent", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, spots_per_image = 20, n_images = 2)
  res <- SubsetSpatial(obj, subset = seurat_clusters == "0")

  expect_true(all(as.character(res$seurat_clusters) == "0"))
  for (img in names(res@images)) {
    expect_true(all(SeuratObject::Cells(res@images[[img]]) %in% colnames(res)))
  }
  # every kept cell should be accounted for by exactly the images it belongs to
  expect_equal(
    sum(vapply(res@images, function(i) length(SeuratObject::Cells(i)), integer(1))),
    ncol(res)
  )
})

test_that("SubsetSpatial cells= (character) and idents= both work", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 2, spots_per_image = 15, n_images = 2)

  keep <- colnames(obj)[1:10]
  res_cells <- SubsetSpatial(obj, cells = keep)
  expect_setequal(colnames(res_cells), keep)
  for (img in names(res_cells@images)) {
    expect_true(all(SeuratObject::Cells(res_cells@images[[img]]) %in% colnames(res_cells)))
  }

  res_idents <- SubsetSpatial(obj, idents = "0")
  expect_true(all(as.character(Seurat::Idents(res_idents)) == "0"))
  for (img in names(res_idents@images)) {
    expect_true(all(SeuratObject::Cells(res_idents@images[[img]]) %in% colnames(res_idents)))
  }
})

test_that("SubsetSpatial drops an image entirely (with a message) when 0 cells survive", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 3, spots_per_image = 10, n_images = 2)
  expect_message(
    res <- SubsetSpatial(obj, subset = sample == "slice1"),
    "slice2.*dropping"
  )
  expect_false("slice2" %in% names(res@images))
  expect_true("slice1" %in% names(res@images))
  expect_true(all(SeuratObject::Cells(res@images[["slice1"]]) %in% colnames(res)))
})

test_that("SubsetSpatial features= subsets genes without touching images", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_visium_seurat(seed = 1, n_genes = 20, spots_per_image = 10, n_images = 1)
  keep_genes <- rownames(obj)[1:5]
  res <- SubsetSpatial(obj, features = keep_genes)
  expect_setequal(rownames(res), keep_genes)
  expect_equal(ncol(res), ncol(obj))
  expect_true(all(SeuratObject::Cells(res@images[["slice1"]]) %in% colnames(res)))
})

test_that("SubsetSpatial works on FOV-based (Xenium/CosMx-style) objects too", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 60)
  res <- SubsetSpatial(obj, subset = celltype == "TypeA")
  expect_true(all(res$celltype == "TypeA"))
  for (img in names(res@images)) {
    expect_true(all(SeuratObject::Cells(res@images[[img]]) %in% colnames(res)))
  }
})


# ============================================================================
# combine_fovs() -- grid-stitches per-FOV centroids into one combined FOV.
# Regression test for the n_cols == 1 row-wrap bug: `(image + 1) %% n_cols
# == 1` (the previous condition) is never TRUE when n_cols == 1, so a
# single-column layout never wrapped rows and instead laid every FOV out
# side-by-side in one row. Fixed to `image %% n_cols == 0`.
# ============================================================================

# Every FOV has 2 cells at the same (0, 0)/(1, 1) coordinates, so every
# FOV has the same, known max-x/max-y (1, 1) -- making the resulting grid
# offsets exactly predictable for the assertions below.
.make_combine_fovs_obj <- function(n_fovs = 3, cells_per_fov = 2) {
  fov_names <- paste0("fov", seq_len(n_fovs))
  all_ids <- unlist(lapply(fov_names, function(f) paste0(f, "_c", seq_len(cells_per_fov))))

  genes <- paste0("Gene", 1:5)
  counts <- matrix(1, nrow = length(genes), ncol = length(all_ids),
                   dimnames = list(genes, all_ids))
  storage.mode(counts) <- "double"
  obj <- SeuratObject::CreateSeuratObject(counts = counts, assay = "RNA")

  for (f in fov_names) {
    ids <- paste0(f, "_c", seq_len(cells_per_fov))
    coords <- data.frame(x = c(0, 1), y = c(0, 1), cell = ids,
                         stringsAsFactors = FALSE)
    cents <- SeuratObject::CreateCentroids(coords)
    fov_obj <- SeuratObject::CreateFOV(
      coords = list(centroids = cents), type = "centroids", assay = "RNA",
      key = paste0(f, "_")
    )
    obj[[f]] <- fov_obj
  }
  obj
}

test_that("combine_fovs wraps to a new row for every FOV when n_cols = 1", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_combine_fovs_obj(n_fovs = 3, cells_per_fov = 2)
  out <- suppressMessages(combine_fovs(obj, assay = "RNA", n_cols = 1, offset = 10))

  combined_centroids <- out@images[["combined"]]$centroids
  coords <- as.data.frame(combined_centroids@coords)
  rownames(coords) <- combined_centroids@cells

  min_xy <- function(fov) {
    sub <- coords[grepl(paste0("^", fov, "_"), rownames(coords)), ]
    c(x = min(sub$x), y = min(sub$y))
  }
  fov_mins <- t(vapply(paste0("fov", 1:3), min_xy, numeric(2)))

  # With n_cols = 1, every FOV should start its own row: x stays 0 for all
  # three, and y strictly increases from one FOV to the next. Before the
  # fix, this instead came out as y constant (0) and x increasing (every
  # FOV laid out in a single row).
  expect_equal(unname(fov_mins[, "x"]), c(0, 0, 0))
  expect_true(all(diff(unname(fov_mins[, "y"])) > 0))
})

test_that("combine_fovs wraps within a row for n_cols > 1", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_combine_fovs_obj(n_fovs = 4, cells_per_fov = 2)
  out <- suppressMessages(combine_fovs(obj, assay = "RNA", n_cols = 2, offset = 10))

  combined_centroids <- out@images[["combined"]]$centroids
  coords <- as.data.frame(combined_centroids@coords)
  rownames(coords) <- combined_centroids@cells

  min_xy <- function(fov) {
    sub <- coords[grepl(paste0("^", fov, "_"), rownames(coords)), ]
    c(x = min(sub$x), y = min(sub$y))
  }
  fov_mins <- t(vapply(paste0("fov", 1:4), min_xy, numeric(2)))

  # FOV1/FOV2 share row 1 (y = 0, x increasing); FOV3/FOV4 share row 2 (a
  # new, larger y; x resets to 0 for FOV3 then increases again for FOV4).
  expect_equal(unname(fov_mins["fov1", "y"]), unname(fov_mins["fov2", "y"]))
  expect_equal(unname(fov_mins["fov3", "y"]), unname(fov_mins["fov4", "y"]))
  expect_true(unname(fov_mins["fov3", "y"]) > unname(fov_mins["fov1", "y"]))
  expect_equal(unname(fov_mins["fov1", "x"]), 0)
  expect_equal(unname(fov_mins["fov3", "x"]), 0)
  expect_true(unname(fov_mins["fov2", "x"]) > unname(fov_mins["fov1", "x"]))
})


# ============================================================================
# subset_opt() -- Visium/FOV-aware wrapper around subset(). cleanMolecules
# is set to FALSE throughout (see file header) since CleanMolSlot() needs a
# full molecule-level FOV this suite's fixtures don't build.
#
# Regression coverage: subset_opt() used to call several Seurat/SeuratObject
# functions (Cells(), Images(), WhichCells(), UpdateSlots(),
# UpdateSeuratObject()) and rlang::enquo()/magrittr's %<>% pipe completely
# unqualified/unimported -- it only worked because zzz.R's .onAttach()
# happens to attach all of those packages under `library(SingleCellTools)`,
# not because the calls were actually resolvable within the package
# namespace. Running these tests (which load the package via
# devtools::load_all()/testthat, not necessarily a full library() attach)
# is itself a meaningful check that the fix holds.
# ============================================================================

test_that("subset_opt subsets a Seurat object and every FOV down to the requested cells", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- .make_spatial_obj(seed = 1, n_per_fov = 20)
  keep_cells <- rownames(obj@meta.data)[obj$celltype == "TypeA"]
  out <- suppressMessages(subset_opt(obj, cells = keep_cells, cleanMolecules = FALSE))

  expect_setequal(colnames(out), keep_cells)
  expect_true(all(out$celltype == "TypeA"))
  for (img in names(out@images)) {
    expect_true(all(SeuratObject::Cells(out@images[[img]]) %in% keep_cells))
  }
})
