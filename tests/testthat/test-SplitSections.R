# Tests for SplitSections() against a synthetic object with two FOVs:
#   fovA: three well-separated sections along y + a 10-cell debris cluster
#   fovB: one section + two sections 2 units apart that DBSCAN joins into one
#         piece (a narrow gap, like serial sections laid down touching)
# `truth` records which simulated section each cell came from, and a numeric
# `fov` column stands in for the CosMx-style column that must not be touched.

.make_sections_obj <- function(seed = 1) {
  set.seed(seed)
  disk <- function(n, cx, cy, r, truth) {
    th <- stats::runif(n, 0, 2 * pi); rr <- r * sqrt(stats::runif(n))
    data.frame(x = cx + rr * cos(th), y = cy + rr * sin(th), truth = truth,
               stringsAsFactors = FALSE)
  }
  layout <- list(
    fovA = rbind(disk(150, 0,   0, 10, "A1"),
                 disk(150, 0, 100, 10, "A2"),
                 disk(150, 0, 200, 10, "A3"),
                 disk( 10, 0, 300,  1, "debris")),
    fovB = rbind(disk(150, 0,   0, 10, "B1"),
                 disk(150, 0, 100,  8, "B2low"),
                 disk(150, 0, 118,  8, "B2high"))
  )

  pieces <- lapply(seq_along(layout), function(i) {
    d <- layout[[i]]
    d$cell <- paste0(names(layout)[i], "_c", seq_len(nrow(d)))
    d$fov_num <- i
    d
  })
  all_cells <- do.call(rbind, pieces)
  n <- nrow(all_cells)

  counts <- matrix(stats::rpois(10 * n, 2), nrow = 10,
                   dimnames = list(paste0("Gene", 1:10), all_cells$cell))
  storage.mode(counts) <- "double"
  meta <- data.frame(truth = all_cells$truth, fov = all_cells$fov_num,
                     row.names = all_cells$cell, stringsAsFactors = FALSE)
  obj <- SeuratObject::CreateSeuratObject(counts = methods::as(counts, "CsparseMatrix"),
                                          meta.data = meta, assay = "RNA")
  for (i in seq_along(layout)) {
    f <- names(layout)[i]
    cents <- SeuratObject::CreateCentroids(pieces[[i]][, c("x", "y", "cell")])
    obj[[f]] <- SeuratObject::CreateFOV(coords = list(centroids = cents), type = "centroids",
                                        assay = "RNA", key = paste0(f, "_"))
  }
  obj
}

test_that("SplitSections finds sections, orders them by y, and drops debris", {
  .skip_if_missing("Seurat", "SeuratObject", "dbscan")
  obj0 <- .make_sections_obj()
  obj  <- suppressMessages(SplitSections(obj0, min_cells = 50))
  md   <- obj@meta.data

  expect_true(all(c("section_fov", "section_auto", "section_id") %in% colnames(md)))
  expect_equal(unname(obj$fov), unname(obj0$fov))          # existing `fov` column untouched

  a <- md[md$section_fov == "fovA", ]
  expect_equal(sort(unique(stats::na.omit(a$section_auto))), c("S01", "S02", "S03"))
  expect_true(all(a$section_auto[a$truth == "A1"] == "S01"))
  expect_true(all(a$section_auto[a$truth == "A3"] == "S03"))
  expect_true(all(is.na(a$section_auto[a$truth == "debris"])))
  expect_equal(unique(a$section_id[a$truth == "A2"]), "fovA_S02")

  b <- md[md$section_fov == "fovB", ]
  expect_equal(sort(unique(stats::na.omit(b$section_auto))), c("S01", "S02"))  # touching pair = 1 piece
})

test_that("SplitSections applies manual rename, merge and drop; unlisted FOVs become NA", {
  .skip_if_missing("Seurat", "SeuratObject", "dbscan")
  fixes <- data.frame(fov = "fovA", from = c("S01", "S02", "S03"), to = c("sec1", "sec1", NA))
  obj <- suppressMessages(SplitSections(.make_sections_obj(), manual = fixes, min_cells = 50))
  md  <- obj@meta.data
  a   <- md[md$section_fov == "fovA", ]

  expect_true(all(a$section_id[a$truth %in% c("A1", "A2")] == "fovA_sec1"))
  expect_true(all(is.na(a$section_id[a$truth == "A3"])))
  expect_equal(unique(a$section_auto[a$truth == "A3"]), "S03")   # automatic label kept
  expect_true(all(is.na(md$section_id[md$section_fov == "fovB"])))
  expect_equal(obj@misc$SplitSections$manual$fov, rep("fovA", 3))
})

test_that("SplitSections splits touching sections with 'a;b', labelled left to right", {
  .skip_if_missing("Seurat", "SeuratObject", "dbscan")
  fixes <- data.frame(fov = "fovB", from = c("S01", "S02"), to = c("sec3", "b;a"))
  obj <- suppressMessages(SplitSections(.make_sections_obj(), manual = fixes, min_cells = 50))
  b   <- obj@meta.data[obj$section_fov == "fovB", ]

  expect_true(all(b$section_id[b$truth == "B1"] == "fovB_sec3"))
  expect_true(all(b$section_id[b$truth == "B2low"]  == "fovB_b"))   # lower y -> first name
  expect_true(all(b$section_id[b$truth == "B2high"] == "fovB_a"))
})

test_that(".ss_split_valley cuts at the gap when neighbouring sections differ in size", {
  # A k-means midpoint cut lands inside the wider section here (~53); the
  # gap search should cut in the 60-70 gap.
  set.seed(3)
  a <- stats::runif(300, 0, 60)
  b <- stats::runif(150, 70, 100)
  lab <- .ss_split_valley(c(a, b), k = 2, bw = 1)
  expect_true(all(lab[seq_along(a)] == 1L))
  expect_true(all(lab[length(a) + seq_along(b)] == 2L))
  expect_equal(.ss_split_valley(a, k = 1, bw = 1), rep(1L, length(a)))
})

test_that(".ss_split_valley finds the gap even when both k-means centres fall on one side of it", {
  # 300 cells over 0-60 next to 100 cells over 70-80: the best 2-means split
  # puts both centres left of the gap (~19 and ~62), so searching only
  # between the centres would miss it.
  set.seed(4)
  a <- stats::runif(300, 0, 60)
  b <- stats::runif(100, 70, 80)
  lab <- .ss_split_valley(c(a, b), k = 2, bw = 1)
  expect_true(all(lab == rep(1:2, c(300, 100))))
})

test_that(".ss_split_valley handles three uneven sections", {
  set.seed(5)
  y <- c(stats::runif(400, 0, 80), stats::runif(120, 88, 100), stats::runif(250, 110, 150))
  lab <- .ss_split_valley(y, k = 3, bw = 1)
  expect_true(all(lab == rep(1:3, c(400, 120, 250))))
})

test_that(".ss_split_valley falls back to k-means midpoints when there is no gap", {
  set.seed(6)
  y <- stats::runif(2000, 0, 100)
  lab <- .ss_split_valley(y, k = 2, bw = 1)
  expect_gt(mean(lab == 1L), 0.4)
  expect_lt(mean(lab == 1L), 0.6)
})

test_that("SplitSections reads a fixes CSV; blank `to` drops the piece; extra columns are ignored", {
  .skip_if_missing("Seurat", "SeuratObject", "dbscan")
  f <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(fov = "fovA", from = c("S01", "S02"), to = c("x", ""),
                              comment = c("keep", "drop")), f, row.names = FALSE)
  obj <- suppressMessages(SplitSections(.make_sections_obj(), manual = f, min_cells = 50))
  a   <- obj@meta.data[obj$section_fov == "fovA", ]

  expect_true(all(a$section_id[a$truth == "A1"] == "fovA_x"))
  expect_true(all(is.na(a$section_id[a$truth == "A2"])))
  expect_true(all(a$section_id[a$truth == "A3"] == "fovA_S03"))  # listed FOV, unlisted piece
})

test_that("SplitSections validates the fixes table", {
  .skip_if_missing("Seurat", "SeuratObject", "dbscan")
  obj0 <- .make_sections_obj()
  expect_error(SplitSections(obj0, manual = data.frame(fov = "fovA", from = "S01"), min_cells = 50),
               "missing column")
  expect_warning(suppressMessages(SplitSections(
    obj0, manual = data.frame(fov = c("fovA", "nope"), from = "S01", to = "s"), min_cells = 50)),
    "not in object")
  expect_error(SplitSections("not a seurat object"))
})

test_that("SplitSections does not leak its seed into the caller's RNG stream", {
  .skip_if_missing("Seurat", "SeuratObject", "dbscan")
  obj0 <- .make_sections_obj()
  set.seed(999); x_ref <- stats::runif(5)
  set.seed(999)
  invisible(suppressMessages(SplitSections(obj0, min_cells = 50, seed = 42)))
  expect_equal(stats::runif(5), x_ref)
})

test_that("SplitSections works on older VisiumV1 images (imagerow/imagecol coordinates)", {
  .skip_if_missing("Seurat", "SeuratObject", "dbscan")
  vis <- .make_visium_seurat(spots_per_image = 40, n_images = 2)
  out <- suppressMessages(SplitSections(vis, min_cells = 5))
  expect_setequal(unique(out$section_fov), names(vis@images))
})

test_that("SplitSections writes check plots, with a _final plot only for FOVs in `manual`", {
  .skip_if_missing("Seurat", "SeuratObject", "dbscan")
  skip_if_not(capabilities("png"))
  d <- tempfile("ss_plots_")
  invisible(suppressMessages(SplitSections(
    .make_sections_obj(), min_cells = 50, plot_dir = d,
    manual = data.frame(fov = "fovA", from = "S01", to = "s1"))))
  expect_true(all(file.exists(file.path(d, c("fovA.png", "fovA_pieces.png", "fovA_final.png",
                                              "fovB.png", "fovB_pieces.png")))))
  expect_false(file.exists(file.path(d, "fovB_final.png")))
})
