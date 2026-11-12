# Tests for detect_fov_edges() / detect_tissue_holes(): list input (shape
# preserved, results identical to per-object calls, fovs matched per
# object) and VisiumV1 coordinates in detect_tissue_holes(). Fixtures:
# .make_spatial_obj() (two FOV images, fov1/fov2) and .make_visium_seurat()
# (VisiumV1 images) from helper-spatial.R, plus a grid with a hole below.

# 30 x 30 grid of cells with an empty 10 x 10 square in the middle, as a
# single centroid FOV named `fov`.
.make_holey_obj <- function(fov = "grid") {
  g <- expand.grid(x = 0:29, y = 0:29)
  g <- g[!(g$x >= 10 & g$x <= 19 & g$y >= 10 & g$y <= 19), ]
  g$cell <- paste0(fov, "_c", seq_len(nrow(g)))
  counts <- matrix(1, nrow = 5, ncol = nrow(g),
                   dimnames = list(paste0("Gene", 1:5), g$cell))
  obj <- SeuratObject::CreateSeuratObject(counts = methods::as(counts, "CsparseMatrix"),
                                          assay = "RNA")
  obj[[fov]] <- SeuratObject::CreateFOV(
    coords = list(centroids = SeuratObject::CreateCentroids(g[, c("x", "y", "cell")])),
    type = "centroids", assay = "RNA", key = paste0(fov, "_"))
  obj
}

# ---- detect_fov_edges() ------------------------------------------------------

test_that("detect_fov_edges returns a named list for list input, matching per-object results", {
  .skip_if_missing("Seurat", "SeuratObject")
  o1 <- .make_spatial_obj(seed = 1); o2 <- .make_spatial_obj(seed = 2)
  res <- suppressMessages(detect_fov_edges(list(a = o1, b = o2), n_iterations = 2))

  expect_type(res, "list")
  expect_named(res, c("a", "b"))
  expect_true(all(vapply(res, inherits, logical(1), "Seurat")))
  one <- suppressMessages(detect_fov_edges(o2, n_iterations = 2))
  expect_equal(res$b$edge_layer, one$edge_layer)
  expect_true(any(res$a$edge_layer > 0))
})

test_that("detect_fov_edges keeps an unnamed list unnamed and a single object single", {
  .skip_if_missing("Seurat", "SeuratObject")
  o1 <- .make_spatial_obj(seed = 1)
  res <- suppressMessages(detect_fov_edges(list(o1, o1)))
  expect_length(res, 2)
  expect_null(names(res))
  expect_s4_class(suppressMessages(detect_fov_edges(o1)), "Seurat")
})

test_that("detect_fov_edges matches `fovs` per object in a list", {
  .skip_if_missing("Seurat", "SeuratObject")
  o1  <- .make_spatial_obj(seed = 1)                       # images fov1, fov2
  vis <- .make_visium_seurat(spots_per_image = 40)         # image slice1
  res <- suppressMessages(detect_fov_edges(list(xen = o1, vis = vis), fovs = "fov1"))

  in_fov2 <- grepl("^fov2_", colnames(res$xen))
  expect_true(all(res$xen$edge_layer[in_fov2] == 0))       # not processed
  expect_true(any(res$xen$edge_layer[!in_fov2] > 0))
  expect_true(all(res$vis$edge_layer == 0))                # no requested FOVs
  expect_error(suppressMessages(detect_fov_edges(list(xen = o1, vis = vis), fovs = "nope")),
               "not in any object")
})

test_that("detect_fov_edges validates input", {
  .skip_if_missing("Seurat", "SeuratObject")
  expect_error(detect_fov_edges("not a seurat object"), "Seurat")
  expect_error(detect_fov_edges(list(.make_spatial_obj(), "x")), "Seurat")
  expect_error(detect_fov_edges(.make_spatial_obj(), fovs = "nope"), "not in obj@images")
})

# ---- detect_tissue_holes() ---------------------------------------------------

test_that("detect_tissue_holes flags the cells bordering an internal hole", {
  .skip_if_missing("Seurat", "SeuratObject")
  obj <- suppressMessages(detect_tissue_holes(.make_holey_obj(), n_iterations = 1))
  xy  <- SeuratObject::GetTissueCoordinates(obj[["grid"]])
  rownames(xy) <- xy$cell
  flagged <- colnames(obj)[obj$hole_layer == 1]
  expect_gt(length(flagged), 0)
  # flagged cells sit next to the hole, not at the outer border
  fx <- xy[flagged, "x"]; fy <- xy[flagged, "y"]
  expect_true(all(fx >= 5 & fx <= 24 & fy >= 5 & fy <= 24))
})

test_that("detect_tissue_holes works on VisiumV1 images (imagerow/imagecol)", {
  # Regression: it used to read x/y columns directly and failed with
  # "undefined columns selected" on VisiumV1 images.
  .skip_if_missing("Seurat", "SeuratObject")
  vis <- .make_visium_seurat(spots_per_image = 40, n_images = 2)
  out <- suppressMessages(detect_tissue_holes(vis, n_iterations = 1))
  expect_true("hole_layer" %in% colnames(out@meta.data))
  expect_true(all(out$hole_layer %in% 0:1))
})

test_that("detect_tissue_holes accepts a mixed list and matches per-object results", {
  .skip_if_missing("Seurat", "SeuratObject")
  holey <- .make_holey_obj()
  vis   <- .make_visium_seurat(spots_per_image = 40)
  res   <- suppressMessages(detect_tissue_holes(list(holey = holey, vis = vis), n_iterations = 1))

  expect_named(res, c("holey", "vis"))
  one <- suppressMessages(detect_tissue_holes(holey, n_iterations = 1))
  expect_equal(res$holey$hole_layer, one$hole_layer)
  expect_true("hole_layer" %in% colnames(res$vis@meta.data))
})

test_that("detect_tissue_holes names the object when exclude_gene is missing from it", {
  .skip_if_missing("Seurat", "SeuratObject")
  a <- .make_spatial_obj(n_genes = 20)       # has Gene20
  b <- .make_spatial_obj(n_genes = 10)       # does not
  expect_error(suppressWarnings(suppressMessages(
    detect_tissue_holes(list(a = a, b = b), exclude_gene = "Gene20"))), "\\[b\\]")
})
