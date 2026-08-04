# Dummy Seurat objects for exercising MergeSeurat()'s HarmonyIntegration
# path (the harmony::RunHarmony() direct-call fix). Not part of the package
# -- this is a throwaway dev/test script, not sourced on package load.

library(Seurat)
set.seed(42)

make_dummy_sample <- function(prefix, n_genes = 200, n_cells = 300,
                               batch_shift = 1) {
  genes <- paste0("Gene", seq_len(n_genes))
  cells <- paste0(prefix, "_c", seq_len(n_cells))

  # batch_shift scales every gene's mean count for this sample -- a crude
  # stand-in for a real batch effect, so the two samples visibly separate
  # on plain PCA but should mix back together on the "harmony" reduction.
  counts <- matrix(
    stats::rpois(n_genes * n_cells, lambda = 3 * batch_shift),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(genes, cells)
  )
  storage.mode(counts) <- "double"

  # project = prefix sets orig.ident, which is MergeSeurat's default
  # group_column -- i.e. what Harmony batch-corrects on.
  CreateSeuratObject(counts = counts, project = prefix)
}

obj1 <- make_dummy_sample("sample1", batch_shift = 1)
obj2 <- make_dummy_sample("sample2", batch_shift = 2.5)

seurat_list <- list(sample1 = obj1, sample2 = obj2)

# ---- Try MergeSeurat() with the Harmony fix --------------------------------
# use_SCT = FALSE / to_regress = NULL because this synthetic data has no
# percent.mt or other real QC metadata to regress against.
merged <- MergeSeurat(
  seurat_list,
  use_SCT                   = FALSE,
  to_regress                = NULL,
  integration                = "HarmonyIntegration",
  integration_normalization  = "LogNormalize",
  integration_assay          = "RNA",
  new_reduction              = "harmony",
  max_dims                   = 10,
  markers                    = FALSE,
  save_rds_file              = FALSE
)

# ---- Sanity check: harmony should mix samples1 & 2; raw pca should not ----
DimPlot(merged, reduction = "pca",     group.by = "orig.ident") +
  ggplot2::ggtitle("Before integration (pca)")
DimPlot(merged, reduction = "harmony", group.by = "orig.ident") +
  ggplot2::ggtitle("After Harmony integration")
