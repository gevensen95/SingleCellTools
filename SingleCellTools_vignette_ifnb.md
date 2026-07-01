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
5. [Gene ID Check](#5-gene-id-check)
6. [QC Report and Filter — `GenerateQCReport()` / `ApplyQCFilters()` / `QCComparePlots()`](#6-qc-report-and-filter)
7. [Merge and Integrate — `MergeSeurat()`](#7-merge-and-integrate--mergeseurat)
8. [Integration Quality — `BatchEffectQC()`](#8-integration-quality--batcheffectqc)
9. [Annotate Cell Types — Three Approaches](#9-annotate-cell-types--three-approaches)
   - 9.1 [Marker Dot Plots — `MarkerPlot()` / `MarkerPctPlot()`](#91-marker-dot-plots--markerplot--markerpctplot)
   - 9.2 [Cluster-level Marker Scoring — `AnnotateClusters()`](#92-cluster-level-marker-scoring--annotateclusters)
   - 9.3 [Reference-based — `AnnotateWithReference()`](#93-reference-based--annotatewithreference)
10. [Flag Gene-Positive Cells — `AddGenePositivity()` / `PlotGenePositivity()`](#10-flag-gene-positive-cells--addgenepositivity--plotgenepositivity)
11. [Cell-Cycle Scoring — `assign_cell_cycle_phase()`](#11-cell-cycle-scoring--assign_cell_cycle_phase)
12. [Cell-Type Composition — `CellComposition()` / `CompositionalTest()`](#12-cell-type-composition--cellcomposition--compositionaltest)
13. [Differential Expression — `PseudobulkDE()` / `PlotVolcano()`](#13-differential-expression--pseudobulkde--plotvolcano)
14. [Cell-Cell Communication — `RunLIANA()`](#14-cell-cell-communication--runliana)
15. [Save with Provenance — `SaveWithProvenance()`](#15-save-with-provenance--savewithprovenance)
16. [Session Info](#16-session-info)

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

## 6. QC Report and Filter

Rather than eyeballing violin plots and picking cutoffs by hand, generate a full HTML QC report with recommended per-sample cutoffs and apply them programmatically.

```r
# Generate HTML report + machine-readable sidecar CSV
GenerateQCReport(
  sample_list,
  output_file    = "qc/ifnb_qc.html",
  metadata_cols  = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  mad_multiplier = 3,
  doublet_col    = "doublet_finder"
)
# writes qc/ifnb_qc.html AND qc/ifnb_qc_cutoffs.csv
```

Apply the recommended cutoffs — doublet filtering happens by default:

```r
filtered <- ApplyQCFilters(
  sample_list,
  cutoffs         = "qc/ifnb_qc_cutoffs.csv",
  filter_doublets = TRUE,
  return_report   = TRUE
)
sample_list <- filtered$obj
filtered$report                # per-sample per-metric retention stats
```

Confirm at a glance that the filter did what you intended:

```r
QCComparePlots(
  pre     = sample_list_pre,   # save this before overwriting
  post    = sample_list,
  metrics = c("nCount_RNA", "nFeature_RNA", "percent.mt")
)
```

This replaces the manual `subset(obj, nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)` step from earlier versions of this vignette.

---

## 7. Merge and Integrate — `MergeSeurat()`

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

## 8. Integration Quality — `BatchEffectQC()`

`BatchEffectQC()` quantifies how well the integration mixed CTRL and STIM (or across donors) versus preserving biology. Run it on both the pre-integration `pca` and post-integration `harmony` reductions to see the effect numerically:

```r
pre  <- BatchEffectQC(integrated, reduction = "pca",
                      batch_col    = "stim",
                      celltype_col = NULL)     # no labels yet
post <- BatchEffectQC(integrated, reduction = "harmony",
                      batch_col    = "stim",
                      celltype_col = NULL)

rbind(pre$summary, post$summary)
#          batch_asw knn_mixing expected_mixing ...
# [pca]         0.28       0.61            0.94
# [harmony]     0.03       0.93            0.94
```

Ideally `knn_mixing` ratio approaches 1 (perfectly random neighborhoods across conditions), and `batch_asw` drops toward 0 (cells no longer cluster by condition). Once you have cell-type labels, re-run with `celltype_col` set to check that biological structure was preserved (`knn_purity` and `celltype_asw` should stay high).

---

## 9. Annotate Cell Types — Three Approaches

`SingleCellTools` v2.4 offers three complementary paths for cell-type annotation. For PBMC-style well-studied tissues, Azimuth is usually the fastest and most accurate; marker-based approaches are more transparent and don't require internet access.

### 9.1 Marker Dot Plots — `MarkerPlot()` / `MarkerPctPlot()`

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

### 9.2 Cluster-level Marker Scoring — `AnnotateClusters()`

Skip the manual `RenameIdents` step by having `AnnotateClusters` score each cluster on every marker set with UCell and pick the winner:

```r
pbmc_marker_list <- list(
  T_cell     = c("CD3D", "CD3E", "CD8A", "CD8B"),
  CD4_T      = c("CD4", "IL7R", "CCR7"),
  NK         = c("GNLY", "NKG7", "KLRD1"),
  B_cell     = c("MS4A1", "CD79A", "CD79B"),
  CD14_Mono  = c("CD14", "LYZ", "CST3"),
  CD16_Mono  = c("FCGR3A", "MS4A7"),
  DC         = c("FCER1A"),
  Platelet   = c("PPBP")
)

integrated <- AnnotateClusters(
  integrated,
  method              = "marker",
  markers             = pbmc_marker_list,
  cluster_col         = "seurat_clusters",
  new_col             = "predicted_cell_type",
  filter_nonspecific  = TRUE,
  min_score           = 0.1,
  min_margin          = 0.05
)

table(integrated$predicted_cell_type)

# Inspect the full per-cluster score matrix
scored <- AnnotateClusters(integrated, markers = pbmc_marker_list,
                           method = "marker",
                           return_scores = "cluster")
scored$scores        # cluster x cell-type UCell mean scores
```

### 9.3 Reference-based — `AnnotateWithReference()`

For PBMCs specifically, Azimuth is essentially the gold standard: it projects onto a large, curated reference and provides both coarse (`l1`) and fine (`l2`) labels with confidence scores.

```r
integrated <- AnnotateWithReference(
  integrated,
  reference          = "pbmcref",
  annotation_levels  = c("l1", "l2"),
  min_score          = 0.5,        # Unknown for low-confidence calls
  unassigned_label   = "Unknown"
)

table(integrated$predicted.celltype.l1)
table(integrated$predicted.celltype.l2)

# Compare to your marker-based labels
table(marker = integrated$predicted_cell_type,
      azimuth = integrated$predicted.celltype.l1)
```

For the rest of the vignette we'll use `cell_type` — pick whichever of the three above you trust most and assign:

```r
integrated$cell_type <- integrated$predicted.celltype.l2   # or predicted_cell_type
Idents(integrated) <- integrated$cell_type
```

---

## 10. Flag Gene-Positive Cells — `AddGenePositivity()` / `PlotGenePositivity()`

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

Visualize positivity across cell types:

```r
# Percent-positive per cluster, one bar per gene
PlotGenePositivity(integrated,
                   c("CD3D", "CD14", "MS4A1", "GNLY", "CD8A"))

# Heatmap for a longer panel
PlotGenePositivity(integrated,
                   c("CD3D", "CD4", "CD8A", "IL7R", "CD14", "LYZ",
                     "FCGR3A", "MS4A1", "GNLY", "NKG7"),
                   style = "heatmap", max_pct = 90)

# Co-expression combinations for a T/B binary
PlotGenePositivity(integrated, c("CD3D", "MS4A1"), style = "combo")
```

---

## 11. Cell-Cycle Scoring — `assign_cell_cycle_phase()`

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

## 12. Cell-Type Composition — `CellComposition()` / `CompositionalTest()`

A natural follow-up question after integration is whether the proportions of each cell type differ between CTRL and STIM.

```r
# Compute proportions + plot
comp <- CellComposition(
  integrated,
  cluster_col = "cell_type",
  sample_col  = "sample_id",
  group_col   = "stim",
  style       = "box"
)
comp$plot
ggsave("composition_ifnb.pdf", comp$plot, width = 10, height = 5)

# Underlying tidy data
head(comp$df)     # sample, cluster, n_cells, prop, group
```

Test statistical significance of the composition shift using propeller (empirical-Bayes-moderated ANOVA on arcsin-sqrt transformed proportions):

```r
comp_test <- CompositionalTest(
  integrated,
  cluster_col   = "cell_type",
  sample_col    = "sample_id",
  condition_col = "stim",
  method        = "auto"        # propeller if available, else betareg, else wilcox
)
subset(comp_test, padj < 0.05)
```

The IFN-β stimulated condition typically shows a shift in monocyte proportions and upregulation of ISGs across all cell types — a well-characterized response that makes `ifnb` a good sanity check that your integration and composition-testing pipelines are working correctly.

---

## 13. Differential Expression — `PseudobulkDE()` / `PlotVolcano()`

The correct way to identify STIM-vs-CTRL DE genes is per-cell-type pseudobulk DE, not `FindMarkers`. `PseudobulkDE()` aggregates per (sample, cell_type) and runs DESeq2:

```r
# Every cell type at once
de_all <- PseudobulkDE(
  integrated,
  sample_col    = "sample_id",
  condition_col = "stim",
  ident_1       = "STIM",
  ident_2       = "CTRL",
  cluster_col   = "cell_type",
  min_cells_per_sample      = 10,
  min_samples_per_condition = 2
)

# Inspect one cell type
head(de_all$CD14_Mono$results)
```

Volcano plot for one cell type:

```r
PlotVolcano(de_all$CD14_Mono$results,
            fc_threshold = 1,
            p_threshold  = 0.05,
            top_n        = 20,
            label_genes  = c("ISG15", "IFI6", "MX1", "IFIT1", "IFIT3"))
```

Compare the DE profiles across two cell types (do they share the same ISGs?):

```r
cmp <- CompareMarkers(
  de_all$CD14_Mono$results,
  de_all$CD4_T$results,
  labels = c("CD14_Mono", "CD4_T")
)
cmp$overlap
cmp$plot                   # log2FC-vs-log2FC scatter, colored by category
```

---

## 14. Cell-Cell Communication — `RunLIANA()`

Compare ligand-receptor interactions between conditions by running LIANA on each stim/ctrl split, then intersecting the top interactions:

```r
# CTRL
lr_ctrl <- RunLIANA(
  subset(integrated, stim == "CTRL"),
  idents_col = "cell_type",
  method     = "consensus",
  min_cells  = 10
)

# STIM
lr_stim <- RunLIANA(
  subset(integrated, stim == "STIM"),
  idents_col = "cell_type",
  method     = "consensus",
  min_cells  = 10
)

# Top interactions unique to STIM
top_ctrl <- head(lr_ctrl$ligand.complex_receptor.complex, 100)
top_stim <- head(lr_stim$ligand.complex_receptor.complex, 100)
setdiff(top_stim, top_ctrl)   # candidate stim-specific interactions
```

For a differential ligand-receptor pipeline that formally tests condition effects, LIANA's own `liana_deconstruct` or CellChat's differential analysis are the recommended paths — `RunLIANA()` is a starting point.

---

## 15. Save with Provenance — `SaveWithProvenance()`

Save the annotated object with a JSON provenance sidecar so downstream readers can inspect what analysis it went through without loading the `.rds`:

```r
SaveWithProvenance(
  integrated,
  file    = "results/ifnb_integrated.rds",
  git_dir = getwd(),
  extra   = list(project = "ifnb_demo",
                 analyst = "K. Evensen",
                 date    = format(Sys.Date()))
)
# Writes results/ifnb_integrated.rds
#    AND results/ifnb_integrated_provenance.json
```

For a quick summary of the object's state after all this work:

```r
CellSuiteSummary(integrated,
                 cluster_col = "cell_type",
                 sample_col  = "sample_id",
                 top_markers = 5)
```

---

## 16. Session Info

```r
sessionInfo()
```

Key packages used in this vignette:

| Package | Role |
|---|---|
| `SingleCellTools` | Doublet calling, QC report/filter, integration, annotation, DE, composition, LIANA, provenance |
| `Seurat` | Core data structure and analysis functions |
| `SeuratData` | `ifnb` dataset |
| `DoubletFinder` | Underlying doublet detection engine |
| `UCell` | Module scoring for annotation and cell-cycle phase |
| `harmony` | Batch correction (called via `IntegrateLayers`) |
| `DESeq2` | Pseudobulk differential expression |
| `speckle` | Propeller composition test |
| `Azimuth` | Reference-based annotation |
| `liana` | Ligand-receptor consensus scoring |
| `patchwork` | QCComparePlots grid layout |
| `ks` | 2D KDE for `PlotFeatureDensity` |
| `jsonlite` | Provenance sidecar |
| `dplyr` / `ggplot2` | Data wrangling and plotting |

---

*For questions or bug reports, open an issue at
[github.com/gevensen95/SingleCellTools/issues](https://github.com/gevensen95/SingleCellTools/issues).*
