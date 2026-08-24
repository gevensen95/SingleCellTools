# Shared fixture for spatial/FOV-based functions: get_all_coords(),
# get_cells_in_polygon(), AnnotateRegions(), NeighborhoodEnrichment(),
# subset_opt()/CleanMolSlot(). Builds a small imaging-based (Xenium/CosMx
# style) Seurat object with two FOVs, cell centroids, and cell-type labels
# clustered in space so neighborhood-enrichment tests have real spatial
# structure to recover.
#
# meta.data is built by pre-allocating a single data.frame and filling it
# by explicit row index, rather than do.call(rbind, <named list>) -- that
# pattern is fragile (rbind.data.frame on a *named* list can rewrite row
# names when reconciling the list's own tags against each piece's own
# row.names, desyncing from the separately-built counts matrix's colnames;
# see helper-NicheCoExpress.R for the same issue root-caused in detail).
# Filling a pre-allocated data.frame by index avoids it entirely.

.make_spatial_obj <- function(seed = 1, n_genes = 20, n_per_fov = 100) {
  set.seed(seed)
  fov_names <- c("fov1", "fov2")
  n_total <- length(fov_names) * n_per_fov

  all_ids  <- character(n_total)
  celltype_v <- character(n_total)
  coords_list <- list()

  for (i in seq_along(fov_names)) {
    fov <- fov_names[i]
    rows <- ((i - 1) * n_per_fov + 1):(i * n_per_fov)
    ids <- paste0(fov, "_c", seq_len(n_per_fov))

    # Two spatial blocks within the FOV: x < 50 -> "TypeA", x >= 50 -> "TypeB"
    # so neighbors of the same type genuinely cluster together in space.
    x <- runif(n_per_fov, 0, 100)
    y <- runif(n_per_fov, 0, 100)
    ct <- ifelse(x < 50, "TypeA", "TypeB")

    all_ids[rows]    <- ids
    celltype_v[rows] <- ct
    coords_list[[fov]] <- data.frame(x = x, y = y, cell = ids,
                                     stringsAsFactors = FALSE)
  }

  counts <- matrix(stats::rpois(n_genes * n_total, lambda = 3), nrow = n_genes,
                   dimnames = list(paste0("Gene", seq_len(n_genes)), all_ids))
  storage.mode(counts) <- "double"

  meta <- data.frame(celltype = celltype_v, row.names = all_ids,
                     stringsAsFactors = FALSE)

  obj <- SeuratObject::CreateSeuratObject(counts = methods::as(counts, "CsparseMatrix"), meta.data = meta, assay = "RNA")

  for (fov in fov_names) {
    cents <- SeuratObject::CreateCentroids(coords_list[[fov]])
    fov_obj <- SeuratObject::CreateFOV(
      coords = list(centroids = cents), type = "centroids", assay = "RNA",
      key = paste0(fov, "_")  # distinct key per FOV; default derives from the
                              # assay and collides across FOVs otherwise
    )
    obj[[fov]] <- fov_obj
  }
  obj
}

# Fixture for SpatialDimPlotFixed()/SpatialFeaturePlotFixed()/
# SpatialConcordance(): a small Seurat object with one or more real VisiumV1
# images (pixel array + scale.factors + imagerow/imagecol coordinates), built
# directly with methods::new() rather than via Seurat::Read10X_Image() --
# there's no real Space Ranger spatial/ directory to read, and constructing
# the S4 object by hand is the standard pattern this file already uses for
# FOV objects above. Each image gets its own non-overlapping slice of spots
# (image_names x spots_per_image), mirroring a real multi-sample merged
# Visium object, so multi-panel tests have per-image spot counts that
# actually differ from the object-wide total.
#
# cluster_by_position = TRUE assigns seurat_clusters by imagecol (left half
# of the image = cluster "0", right half = cluster "1") instead of randomly
# -- SpatialConcordance() needs a fixture with real spatial structure to
# confirm it can actually recover above-null concordance, the same way
# .make_spatial_obj()'s x<50 TypeA/TypeB split serves NeighborhoodEnrichment.
.make_visium_seurat <- function(seed = 1, n_genes = 20, spots_per_image = 15,
                                h = 50, w = 60, n_images = 1,
                                cluster_by_position = FALSE) {
  set.seed(seed)
  genes <- paste0("Gene", seq_len(n_genes))
  image_names <- paste0("slice", seq_len(n_images))
  n_spots <- spots_per_image * n_images
  spots <- paste0("spot", seq_len(n_spots))

  counts <- matrix(stats::rpois(n_genes * n_spots, lambda = 3), nrow = n_genes,
                   dimnames = list(genes, spots))
  storage.mode(counts) <- "double"

  # `scalefactors` isn't an S4 class -- it's an S3 class defined in Seurat
  # (not SeuratObject) via `structure(list(...), class = "scalefactors")`
  # and registered with `setOldClass()` so S4 slots (VisiumV1@scale.factors)
  # can hold it. `setOldClass()` is what produces the "virtual class" error
  # if you try `methods::new("scalefactors", ...)` on it -- build the S3
  # object directly instead, exactly like Seurat's own (possibly
  # unexported) `scalefactors()` constructor does internally.
  sf <- structure(
    list(spot = 1, fiducial = 1, hires = 1, lowres = 0.5),
    class = "scalefactors"
  )

  # Build each image's coordinates up front so cluster_by_position can derive
  # labels from them before the Seurat object (which needs the full meta.data
  # up front) is created.
  coords_list <- vector("list", n_images)
  cluster_v   <- character(n_spots)
  for (i in seq_len(n_images)) {
    img_name <- image_names[i]
    idx <- ((i - 1) * spots_per_image + 1):(i * spots_per_image)
    these_spots <- spots[idx]
    imagecol <- stats::runif(spots_per_image, 5, w - 5)
    imagerow <- stats::runif(spots_per_image, 5, h - 5)
    coords_list[[i]] <- data.frame(
      tissue    = 1L,
      row       = sample(0:20, spots_per_image, replace = TRUE),
      col       = sample(0:20, spots_per_image, replace = TRUE),
      imagerow  = imagerow,
      imagecol  = imagecol,
      row.names = these_spots
    )
    cluster_v[idx] <- if (isTRUE(cluster_by_position)) {
      ifelse(imagecol < w / 2, "0", "1")
    } else {
      sample(c("0", "1"), spots_per_image, replace = TRUE)
    }
  }

  meta <- data.frame(
    seurat_clusters = factor(cluster_v),
    sample           = rep(image_names, each = spots_per_image),
    row.names        = spots,
    stringsAsFactors = FALSE
  )

  obj <- SeuratObject::CreateSeuratObject(counts = methods::as(counts, "CsparseMatrix"), meta.data = meta, assay = "Spatial")
  SeuratObject::LayerData(obj, assay = "Spatial", layer = "data") <- log1p(counts)
  Seurat::Idents(obj) <- meta$seurat_clusters

  for (i in seq_len(n_images)) {
    img_name <- image_names[i]
    img_arr <- array(stats::runif(h * w * 3), dim = c(h, w, 3))
    coords <- coords_list[[i]]
    vis_img <- methods::new(
      "VisiumV1",
      assay         = "Spatial",
      key           = paste0(img_name, "_"),
      image         = img_arr,
      scale.factors = sf,
      coordinates   = coords,
      spot.radius   = 0.008
    )
    obj@images[[img_name]] <- vis_img
  }
  obj
}
