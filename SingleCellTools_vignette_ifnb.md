# SingleCellTools Vignette: PBMC Integration with the ifnb Dataset

**Package:** `SingleCellTools`  
**Data:** `ifnb` — 13,999 human PBMCs (control vs. IFN-β stimulated), from `SeuratData`  
**Goal:** Walk through a complete multi-sample scRNA-seq workflow using only built-in
Seurat datasets — no raw CellRanger output required.

---

## Table of Contents

1. [Setup](#1-setup)
2. [Load and Inspect the Data](#2-load-and-inspect-the-data)
3. [Split into a Sample List](#3-split-into-a-sample-list)
4. [Doublet Detection — `calldoublet()`](#4-doublet-detection--calldoublet)
5. [Gene ID Check — `detect_gene_id_type()` and `check_gene_ids_across_objects()`](#5-gene-id-check)
6. [Merge and Integrate — `MergeSeurat()`](#6-merge-and-integrate--mergeseurat)
7. [Annotate Cell Types — `MarkerPlot()`](#7-annotate-cell-types--markerplot)
8. [Flag Gene-Positive Cells — `AddGenePositivity()`](#8-flag-gene-positive-cells--addgenepositivity)
9. [Cell-Cycle Scoring — `assign_cell_cycle_phase()`](#9-cell-cycle-scoring--assign_cell_cycle_phase)
10. [Cell-Type Composition](#10-cell-type-composition)
11. [Session Info](#11-session-info)

---

## 1. Setup

### Install required packages

```r
# SingleCellTools
remotes::install_github("gevensen95/SingleCellTools")

# DoubletFinder (GitHub-only)
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")

# Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("glmGamPoi", "UCell"))

# SeuratData — provides the ifnb dataset
remotes::install_github("satijalab/seurat-data")
```

### Install and load the ifnb dataset

```r
library(SeuratData)
InstallData("ifnb")   # downloads ~45 MB once; cached for future sessions
```

### Load libraries

```r
library(SingleCellTools)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(dplyr)
```

---

## 2. Load and Inspect the Data

```r
data("ifnb")

# Basic overview
ifnb
#> An object of class Seurat
#> 14053 features across 13999 samples within 1 assay
#> Active assay: RNA (14053 features, 0 variable features)
#> 2 layers present: counts, data

# The key metadata column is 'stim': CTRL or STIM
table(ifnb$stim)
#> CTRL  STIM
#> 6548  7451

# Also has donor IDs
table(ifnb$donor)
```

The `ifnb` object contains two biological conditions — resting PBMCs (`CTRL`) and PBMCs
treated with interferon-beta (`STIM`) — collected from the same donors. This makes it
an ideal test case for batch-corrected integration: the biology changes, but we want
to align the donors properly.

Before going further, make sure the counts layer is accessible. The `ifnb` object from
`SeuratData` ships as a v3 Seurat object; join the layers so it behaves like a standard
v5 object:

```r
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])
```

---

## 3. Split into a Sample List

`SingleCellTools` integration functions expect a **named list of Seurat objects**, one per
sample. The cleanest approach here is to split by donor within each condition, giving
us one object per (donor, condition) pair — exactly what you'd have coming out of
`CreateRNAObjects()` with real CellRanger output.

```r
# Add a combined sample ID column
ifnb$sample_id <- paste(ifnb$stim, ifnb$donor, sep = "_")
table(ifnb$sample_id)

# Split into a list — one object per sample
sample_list <- SplitObject(ifnb, split.by = "sample_id")

# How many cells per sample?
sapply(sample_list, ncol)
```

We also need `percent.mt` in metadata — `calldoublet()` uses it by default as a
regression covariate, and `MergeSeurat()` will regress it out during normalization.
Human mitochondrial genes start with `MT-` (uppercase).

```r
sample_list <- lapply(sample_list, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
})
```

---

## 4. Doublet Detection — `calldoublet()`

`calldoublet()` wraps the full DoubletFinder workflow — normalize, PCA, auto-detect
significant PCs, pK sweep, doublet calling — into a single function. It strips
intermediate reductions and normalized layers from the returned object, keeping
it lightweight.

> **Note:** DoubletFinder is designed for droplet-based data where each GEM is one cell.
> Running it on a pre-split object that already has cell-type labels is still valid — it
> finds technical doublets within the sample.

```r
# Run calldoublet() on every sample
# This may take a few minutes per sample
sample_list <- setNames(
  lapply(seq_along(sample_list), function(i) {
    message(sprintf("[%d/%d] Calling doublets: %s",
                    i, length(sample_list), names(sample_list)[i]))
    calldoublet(
      obj                = sample_list[[i]],
      samplenameIndex    = i,
      normalization      = "LogNormalize",
      vars.to.regress    = "percent.mt",
      cluster_resolution = 0.1
    )
  }),
  names(sample_list)
)

# Inspect doublet calls across samples
lapply(sample_list, function(obj) table(obj$doublet_finder))
```

Review the doublet fractions — a healthy droplet experiment typically yields 1–8%
doublets. If a sample looks unusually high, check the cell count (very high cell
loading increases multiplet rate) or QC metrics.

Filter doublets before integration:

```r
sample_list <- lapply(sample_list, function(obj) {
  n_before <- ncol(obj)
  obj      <- subset(obj, doublet_finder == "Singlet")
  message(sprintf("  Dropped %d doublets; %d singlets remaining",
                  n_before - ncol(obj), ncol(obj)))
  obj
})
```

Apply a basic QC filter while we're here:

```r
sample_list <- lapply(sample_list, function(obj) {
  subset(obj,
    nFeature_RNA > 200  &
    nFeature_RNA < 5000 &
    percent.mt   < 20
  )
})

# Confirm cell counts
sapply(sample_list, ncol)
```

---

## 5. Gene ID Check

Before merging samples from different sources it is good practice to confirm that all
objects use the same gene identifier format. For `ifnb` this is trivially true since
it all came from one object, but in real projects samples often come from different
pipelines or genome annotations.

```r
# Check the identifier type for a single object
detect_gene_id_type(sample_list[[1]])
#> [1] "HGNC symbol"

# Check consistency across the whole list
check_gene_ids_across_objects(sample_list)
#> All objects use: HGNC symbol
```

If this ever returns a mix (e.g., some objects have Ensembl IDs and others have gene
symbols), you will lose almost all genes at the merge step. Fix identifier format before
proceeding.

---

## 6. Merge and Integrate — `MergeSeurat()`

`MergeSeurat()` handles normalization, PCA, integration, clustering, UMAP, and
(optionally) marker detection in one call. We use Harmony integration here, which is
fast and works well for this dataset.

```r
integrated <- MergeSeurat(
  seurat_objects = sample_list,

  # Regression during SCTransform
  to_regress     = "percent.mt",
  use_SCT        = TRUE,

  # Dimensionality
  max_dims       = 20,

  # Harmony integration
  integration    = "HarmonyIntegration",
  new_reduction  = "harmony",

  # Clustering
  cluster_resolution = 0.5,

  # Save outputs
  save_rds_file = TRUE,
  file_name     = "ifnb",

  # Run FindAllMarkers and save marker table + dot plot
  markers       = TRUE
)
```

When this finishes, your working directory will contain:

- `ifnb_merged_seurat_objects.rds` — the integrated object
- `dimplot_seurat_clusters.pdf` — UMAP colored by Seurat cluster
- `markers_all.csv` — `FindAllMarkers` results
- `marker_plot.pdf` — top-10 marker dot plot per cluster

### Quick sanity checks

```r
# UMAP colored by cluster
DimPlot(integrated, label = TRUE, repel = TRUE)

# Verify that CTRL and STIM cells are interleaved (good integration)
DimPlot(integrated, group.by = "stim", cols = c("steelblue", "tomato"))

# Cells per cluster
table(Idents(integrated))
```

If CTRL and STIM cells form separate blobs by condition in the UMAP, the integration
did not fully correct for the stimulation effect. In that case, try increasing
`max_dims`, adjusting `cluster_resolution`, or switching to `RPCAIntegration` with
`k_anchor = 20, k_weight = 100`.

---

## 7. Annotate Cell Types — `MarkerPlot()`

Now we assign cell-type labels to clusters using canonical PBMC markers. `MarkerPlot()`
builds an annotated dot plot, groups the genes by cell type, and clusters the
identities by expression correlation so similar populations appear next to each other.

```r
# Canonical human PBMC marker panel
pbmc_markers <- data.frame(
  Gene = c(
    # T cells
    "CD3D",  "CD3E",  "CD8A",   "CD8B",
    # CD4 T cells
    "CD4",   "IL7R",  "CCR7",
    # NK cells
    "GNLY",  "NKG7",  "KLRD1",
    # B cells
    "MS4A1", "CD79A", "CD79B",
    # Monocytes (classical)
    "CD14",  "LYZ",   "CST3",
    # Monocytes (non-classical)
    "FCGR3A","MS4A7",
    # Dendritic cells
    "FCER1A","CST3",
    # Platelets
    "PPBP"
  ),
  CellType = c(
    rep("T cell", 4),
    rep("CD4 T", 3),
    rep("NK", 3),
    rep("B cell", 3),
    rep("CD14 Mono", 3),
    rep("CD16 Mono", 2),
    rep("DC", 2),
    "Platelet"
  )
)

# Ensure we're plotting against the RNA assay
DefaultAssay(integrated) <- "RNA"

p <- MarkerPlot(
  obj              = integrated,
  genes            = pbmc_markers,
  assay            = "RNA",
  cluster          = TRUE,       # cluster identities by correlation
  show.annotations = TRUE,
  maxsize          = 5,
  label.fontsize   = 3
)
print(p)
ggsave("markerplot_ifnb.pdf", p, width = 14, height = 9)
```

`MarkerPlot()` will automatically drop any gene absent from the assay or with zero
expression, and report those drops as messages — no manual filtering needed.

### Assign cell-type labels

After reviewing the dot plot, rename clusters:

```r
# Adjust these mappings based on what you see in the plot
new_labels <- c(
  "0"  = "CD14 Mono",
  "1"  = "CD4 T naive",
  "2"  = "CD4 T memory",
  "3"  = "CD14 Mono",
  "4"  = "B cell",
  "5"  = "NK",
  "6"  = "CD8 T",
  "7"  = "CD16 Mono",
  "8"  = "T cell",
  "9"  = "DC",
  "10" = "B cell",
  "11" = "Platelet"
)

integrated <- RenameIdents(integrated, new_labels)
integrated$cell_type <- Idents(integrated)

# Labeled UMAP
DimPlot(integrated, label = TRUE, repel = TRUE) +
  ggtitle("PBMC cell types — ifnb dataset") +
  theme(legend.position = "right")
ggsave("umap_celltypes_ifnb.pdf", width = 10, height = 7)
```

---

## 8. Flag Gene-Positive Cells — `AddGenePositivity()`

`AddGenePositivity()` adds a logical metadata column for each gene indicating whether
a cell expresses it above a threshold. This is useful for gating, subsetting, or
computing co-expression fractions.

```r
# Flag expression of key lineage markers
integrated <- AddGenePositivity(
  seurat_objects = integrated,
  genes          = c("CD3D", "CD14", "MS4A1", "GNLY", "CD8A"),
  layer          = "counts",
  threshold      = 0,       # any count > 0 is positive
  suffix         = "_pos"
)

# New logical columns in metadata
head(integrated@meta.data[, c("CD3D_pos", "CD14_pos", "MS4A1_pos",
                               "GNLY_pos", "CD8A_pos")])

# Fraction positive per cell type
positivity_summary <- integrated@meta.data %>%
  group_by(cell_type) %>%
  summarise(
    pct_CD3D  = mean(CD3D_pos)  * 100,
    pct_CD14  = mean(CD14_pos)  * 100,
    pct_MS4A1 = mean(MS4A1_pos) * 100,
    pct_GNLY  = mean(GNLY_pos)  * 100,
    pct_CD8A  = mean(CD8A_pos)  * 100
  )
print(positivity_summary)
```

You can also use positivity flags for subsetting:

```r
# Isolate CD14+ monocytes for downstream re-analysis
mono <- subset(integrated, CD14_pos == TRUE)
ncol(mono)
```

---

## 9. Cell-Cycle Scoring — `assign_cell_cycle_phase()`

`assign_cell_cycle_phase()` uses UCell module scoring to assign S, G2M, or G1 phase to
each cell. It is more robust than Seurat's `CellCycleScoring` on datasets with lower
sequencing depth because UCell normalizes for library size internally.

```r
integrated <- assign_cell_cycle_phase(integrated)

# New metadata columns: S.Score, G2M.Score, Phase
table(integrated$Phase)
#>   G1  G2M    S
#> 9821 1043 1872  (approximate)

# Visualize phase on UMAP
DimPlot(integrated, group.by = "Phase",
        cols = c("G1" = "grey80", "S" = "steelblue", "G2M" = "tomato"))
ggsave("umap_cellcycle_ifnb.pdf", width = 8, height = 6)

# Phase breakdown per cell type
table(integrated$cell_type, integrated$Phase)
```

If cycling cells cluster separately from their resting counterparts and you want to
remove this confound, re-run `MergeSeurat()` with cell-cycle scores added to
`to_regress`:

```r
integrated <- MergeSeurat(
  seurat_objects = sample_list,
  to_regress     = c("percent.mt", "S.Score", "G2M.Score"),
  use_SCT        = TRUE,
  integration    = "HarmonyIntegration",
  new_reduction  = "harmony",
  cluster_resolution = 0.5
)
```

---

## 10. Cell-Type Composition

A natural follow-up question after integration is whether the proportions of each cell
type differ between CTRL and STIM. The `CompositionAnalysis()` function (if available
in your version) computes per-sample counts and proportions:

```r
comp <- CompositionAnalysis(
  integrated,
  group_col    = "cell_type",
  sample_col   = "sample_id",
  condition_col = "stim"
)

# Inspect the proportion table
head(comp$proportions)

# Stacked bar plot
CompositionBarplot(comp, fill = "cell_type", facet = "stim")
ggsave("composition_barplot_ifnb.pdf", width = 10, height = 5)
```

For a manual alternative using base Seurat + dplyr:

```r
comp_manual <- integrated@meta.data %>%
  count(stim, cell_type) %>%
  group_by(stim) %>%
  mutate(proportion = n / sum(n) * 100) %>%
  ungroup()

ggplot(comp_manual, aes(x = stim, y = proportion, fill = cell_type)) +
  geom_col() +
  labs(x = "Condition", y = "% of cells", fill = "Cell type",
       title = "PBMC composition: CTRL vs. IFN-β") +
  theme_classic()
ggsave("composition_manual_ifnb.pdf", width = 6, height = 5)
```

The IFN-β stimulated condition typically shows a shift in monocyte proportions and
upregulation of ISGs (interferon-stimulated genes) across all cell types — a well-
characterized response that makes `ifnb` a good sanity check that your integration
pipeline is working correctly.

---

## 11. Session Info

```r
sessionInfo()
```

Key packages used in this vignette:

| Package | Role |
|---|---|
| `SingleCellTools` | Doublet calling, integration, annotation, plotting |
| `Seurat` | Core data structure and analysis functions |
| `SeuratData` | `ifnb` dataset |
| `DoubletFinder` | Underlying doublet detection engine |
| `UCell` | Module scoring for cell-cycle phase assignment |
| `harmony` | Batch correction (called via `IntegrateLayers`) |
| `dplyr` / `ggplot2` | Data wrangling and plotting |

---

*For questions or bug reports, open an issue at
[github.com/gevensen95/SingleCellTools/issues](https://github.com/gevensen95/SingleCellTools/issues).*
