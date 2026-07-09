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

  obj <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta, assay = "RNA")

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
