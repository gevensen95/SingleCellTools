# General-purpose synthetic Seurat object fixture shared across test files
# for functions that just need *some* reasonably-shaped Seurat RNA object
# with counts/data layers and a few common metadata columns -- not specific
# to any one analysis function. Functions with more specific structural
# needs (spatial FOVs, ATAC fragments, niche co-expression layout, etc.)
# get their own dedicated helper files (see helper-NicheCoExpress.R).

.skip_if_missing <- function(...) {
  pkgs <- c(...)
  for (p in pkgs) testthat::skip_if_not_installed(p)
}

# n_genes genes x n_cells cells, `n_clusters` seurat_clusters levels, one
# `sample` column with `n_samples` levels, one `condition` column derived
# from sample (2 levels, alternating), and one `celltype` column with
# `n_celltypes` levels. Counts are small non-negative integers (Poisson);
# `data` layer is a simple log1p transform -- good enough for structural /
# validation tests, not meant to represent realistic biology.
.make_small_seurat <- function(seed = 1, n_genes = 30, n_cells = 80,
                               n_clusters = 3, n_samples = 4,
                               n_celltypes = 2, gene_prefix = "Gene") {
  set.seed(seed)
  genes <- paste0(gene_prefix, seq_len(n_genes))
  cells <- paste0("cell", seq_len(n_cells))

  counts <- matrix(
    stats::rpois(n_genes * n_cells, lambda = 3),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(genes, cells)
  )
  storage.mode(counts) <- "double"

  samples  <- paste0("S", seq_len(n_samples))
  sample_v <- sample(samples, n_cells, replace = TRUE)
  cond_map <- setNames(rep(c("A", "B"), length.out = n_samples), samples)

  meta <- data.frame(
    seurat_clusters  = factor(sample(seq_len(n_clusters) - 1, n_cells, replace = TRUE)),
    sample           = sample_v,
    condition        = unname(cond_map[sample_v]),
    celltype         = sample(paste0("Type", seq_len(n_celltypes)), n_cells, replace = TRUE),
    row.names        = cells,
    stringsAsFactors = FALSE
  )

  obj <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)
  SeuratObject::LayerData(obj, assay = "RNA", layer = "data") <- log1p(counts)
  Seurat::Idents(obj) <- meta$seurat_clusters
  obj
}
