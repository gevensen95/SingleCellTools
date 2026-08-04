# Dummy spatial single-cell data: 2 Xenium/CosMx-style samples, 100 genes
# each, cell centroids + FOV. Not part of the package -- throwaway dev/test
# script, not sourced on package load.
#
# Built for two things:
#   1. NeighborhoodEnrichment() / NicheCoExpress()  -- cell types are
#      spatially clustered (not randomly scattered), so same-type spatial
#      enrichment is real, not noise.
#   2. MergeSeurat(banksy = TRUE) / RunBanksyWrapper() -- assay is named
#      "Xenium" and each sample carries its own FOV, matching the
#      spatial = "Xenium" path already used elsewhere in this package
#      (see tests/testthat/test-integration.R's .make_xenium_like_sample).

library(Seurat)

make_dummy_spatial_sample <- function(prefix, seed, n_genes = 100, n_cells = 250) {
  set.seed(seed)
  genes <- paste0("Gene", seq_len(n_genes))
  cells <- paste0(prefix, "_c", seq_len(n_cells))

  # Two spatial blocks (x < 50 -> "TypeA", x >= 50 -> "TypeB") so cells of
  # the same type genuinely cluster together in space -- gives
  # NeighborhoodEnrichment() real same-type enrichment to detect, rather
  # than fully random coordinates.
  x  <- stats::runif(n_cells, 0, 100)
  y  <- stats::runif(n_cells, 0, 100)
  ct <- ifelse(x < 50, "TypeA", "TypeB")

  # The first 10 genes are mild "TypeA markers" (higher mean in TypeA cells),
  # so cell-type-aware analyses have a real expression signal to find, not
  # just a spatial/metadata label with nothing behind it.
  lambda_mat <- matrix(3, nrow = n_genes, ncol = n_cells,
                       dimnames = list(genes, cells))
  marker_genes <- genes[1:10]
  lambda_mat[marker_genes, ct == "TypeA"] <- 8

  counts <- matrix(stats::rpois(n_genes * n_cells, lambda = as.vector(lambda_mat)),
                   nrow = n_genes, dimnames = list(genes, cells))
  storage.mode(counts) <- "double"

  meta <- data.frame(celltype = ct, row.names = cells, stringsAsFactors = FALSE)

  # assay = "Xenium" + project = prefix (-> orig.ident) matches what
  # MergeSeurat(spatial = "Xenium") and RunBanksyWrapper() expect.
  obj <- CreateSeuratObject(counts = counts, meta.data = meta,
                            assay = "Xenium", project = prefix)

  coords <- data.frame(x = x, y = y, cell = cells, stringsAsFactors = FALSE)
  cents  <- CreateCentroids(coords)
  fov    <- CreateFOV(coords = list(centroids = cents), type = "centroids",
                      assay = "Xenium", key = paste0(prefix, "fov_"))
  obj[[paste0(prefix, "fov")]] <- fov
  obj
}

sample1 <- make_dummy_spatial_sample("sample1", seed = 1)
sample2 <- make_dummy_spatial_sample("sample2", seed = 2)
spatial_list <- list(sample1 = sample1, sample2 = sample2)

# ---- 1. NeighborhoodEnrichment() on a single sample ------------------------
# (Operates on one Seurat object at a time -- run per sample, or merge first
# if you want enrichment computed across both together.)
ne1 <- NeighborhoodEnrichment(sample1, group.by = "celltype", k = 10,
                              n_perm = 50, assign_niches = TRUE, n_niches = 2)
str(ne1, max.level = 1)

# ---- 2. MergeSeurat(banksy = TRUE) across both samples ---------------------
merged <- MergeSeurat(
  spatial_list,
  spatial                   = "Xenium",
  use_SCT                   = FALSE,
  to_regress                = NULL,
  banksy                    = TRUE,
  banksy_lambda              = 0.2,
  banksy_k_geom              = 10,
  integration                = "HarmonyIntegration",
  new_reduction               = "harmony",
  max_dims                    = 10,
  markers                     = FALSE,
  save_rds_file               = FALSE
)

# ---- 3. NicheCoExpress() -- needs more than 2 samples to be meaningful -----
# NicheCoExpress() compares two *conditions*, each needing >= 2 samples
# (min_samples) for its statistical test. With only 2 samples total here,
# putting one sample in each condition gives n=1 per condition, which isn't
# enough for a real comparison. If you want to exercise NicheCoExpress
# properly, either:
#   (a) treat sample1/sample2 as two replicates of the SAME condition (no
#       real "differential" test, just confirms the plumbing runs), or
#   (b) ask for 4 dummy samples (2 per condition) instead of 2.
