SingleCellTools Vignette: PBMC Integration with the ifnb Dataset
================

**Package:** `SingleCellTools`  
**Data:** `ifnb` — 13,999 human PBMCs (control vs. IFN-β stimulated),
from `SeuratData`  
**Goal:** Walk through a complete multi-sample scRNA-seq workflow using
only built-in Seurat datasets — no raw CellRanger output required.

------------------------------------------------------------------------

## Table of Contents

1.  [Setup](#1-setup)
2.  [Load and Inspect the Data](#2-load-and-inspect-the-data)
3.  [Split into a Sample List](#3-split-into-a-sample-list)
4.  [Doublet Detection —
    `calldoublet()`](#4-doublet-detection--calldoublet)
5.  [Gene ID Check](#5-gene-id-check)
6.  [QC Report and Filter — `GenerateQCReport()` / `ApplyQCFilters()` /
    `QCComparePlots()`](#6-qc-report-and-filter)
7.  [Merge and Integrate —
    `MergeSeurat()`](#7-merge-and-integrate--mergeseurat)
8.  [Integration Quality —
    `BatchEffectQC()`](#8-integration-quality--batcheffectqc)
9.  [Annotate Cell Types — Three
    Approaches](#9-annotate-cell-types--three-approaches)
    - 9.1 [Marker Dot Plots — `MarkerPlot()` /
      `MarkerPctPlot()`](#91-marker-dot-plots--markerplot--markerpctplot)
    - 9.2 [Cluster-level Marker Scoring —
      `AnnotateClusters()`](#92-cluster-level-marker-scoring--annotateclusters)
    - 9.3 [Reference-based —
      `AnnotateWithReference()`](#93-reference-based--annotatewithreference)
10. [Flag Gene-Positive Cells — `AddGenePositivity()` /
    `PlotGenePositivity()`](#10-flag-gene-positive-cells--addgenepositivity--plotgenepositivity)
    - 10.1 [`GenePositivityAnalysis()` and
      `GenePositivityEstimationPlot()`](#101-genepositivityanalysis-and-genepositivityestimationplot)
11. [Cell-Cycle Scoring —
    `assign_cell_cycle_phase()`](#11-cell-cycle-scoring--assign_cell_cycle_phase)
12. [Cell-Type Composition — `CellComposition()` /
    `CompositionalTest()`](#12-cell-type-composition--cellcomposition--compositionaltest)
    - 12.1 [`CompositionAnalysis()` and
      `CompositionEstimationPlot()`](#121-compositionanalysis-and-compositionestimationplot)
13. [Differential Expression — `PseudobulkDE()` /
    `PlotVolcano()`](#13-differential-expression--pseudobulkde--plotvolcano)
14. [Cell-Cell Communication —
    `RunLIANA()`](#14-cell-cell-communication--runliana)
    - 14.1 [`RunCellChat()`](#141-runcellchat)
15. [Save with Provenance —
    `SaveWithProvenance()`](#15-save-with-provenance--savewithprovenance)
16. [Session Info](#16-session-info)

------------------------------------------------------------------------

## 1. Setup

### Install required packages

``` r
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

``` r
library(SeuratData)
InstallData("ifnb")   # downloads ~45 MB once; cached for future sessions
```

### Load libraries

``` r
library(SingleCellTools)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(dplyr)
```

------------------------------------------------------------------------

## 2. Load and Inspect the Data

``` r
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

The `ifnb` object contains two biological conditions — resting PBMCs
(`CTRL`) and PBMCs treated with interferon-beta (`STIM`) — collected
from the same donors. This makes it an ideal test case for
batch-corrected integration: the biology changes, but we want to align
the donors properly.

Before going further, make sure the counts layer is accessible. The
`ifnb` object from `SeuratData` ships as a v3 Seurat object; join the
layers so it behaves like a standard v5 object:

``` r
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])
```

> **A note on memory for real datasets.** At ~14K cells, `ifnb`
> comfortably fits in memory and none of this matters here. For a real
> dataset large enough that it doesn’t – hundreds of thousands of cells,
> many samples, or both – `ConvertToBPCells()` moves any existing Seurat
> object’s counts layer to an on-disk `BPCells` matrix, in place, so the
> rest of a normal Seurat pipeline keeps working against it unmodified.
> It works retroactively on any object regardless of how it was built
> (this `ifnb` object included, just to show the call):
>
> ``` r
> # Illustrative only -- ifnb is small enough that this isn't actually needed.
> ifnb_on_disk <- ConvertToBPCells(ifnb, path = "bpcells_ifnb")
> ```
>
> `CreateRNAObjects()` also has this built in via `on_disk = TRUE`, for
> when you’re reading from raw CellRanger output rather than a
> pre-packaged `SeuratData` object. See
> [`SingleCellTools_vignette_bpcells.md`](SingleCellTools_vignette_bpcells.md).

------------------------------------------------------------------------

## 3. Split into a Sample List

`SingleCellTools` integration functions expect a **named list of Seurat
objects**, one per sample. The cleanest approach here is to split by
donor within each condition, giving us one object per (donor, condition)
pair — exactly what you’d have coming out of `CreateRNAObjects()` with
real CellRanger output.

``` r
# Add a combined sample ID column
ifnb$sample_id <- paste(ifnb$stim, ifnb$donor, sep = "_")
table(ifnb$sample_id)

# Split into a list — one object per sample
sample_list <- SplitObject(ifnb, split.by = "sample_id")

# How many cells per sample?
sapply(sample_list, ncol)
```

We also need `percent.mt` in metadata — `calldoublet()` uses it by
default as a regression covariate, and `MergeSeurat()` will regress it
out during normalization. Human mitochondrial genes start with `MT-`
(uppercase).

``` r
sample_list <- lapply(sample_list, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
})
```

------------------------------------------------------------------------

## 4. Doublet Detection — `calldoublet()`

`calldoublet()` wraps the full DoubletFinder workflow — normalize, PCA,
auto-detect significant PCs, pK sweep, doublet calling — into a single
function. It strips intermediate reductions and normalized layers from
the returned object, keeping it lightweight.

> **Note:** DoubletFinder is designed for droplet-based data where each
> GEM is one cell. Running it on a pre-split object that already has
> cell-type labels is still valid — it finds technical doublets within
> the sample.

> **Speed:** most of that per-sample time is DoubletFinder’s pK
> parameter sweep, which reprocesses a real+artificial-doublet matrix
> from scratch for each of 6 fixed pN values. `pk_sweep_max_cells`
> (default `4000`) estimates pK from a random subsample instead of every
> cell, and `sweep_cores` parallelizes the 6 pN passes internally — see
> [Section 3.3 of the main
> vignette](SingleCellTools_vignette.html#33-doublet-detection--calldoublet)
> for the tradeoffs. Both are shown below.

``` r
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
      cluster_resolution = 0.1,
      pk_sweep_max_cells = 4000,
      sweep_cores        = 1
    )
  }),
  names(sample_list)
)

# Inspect doublet calls across samples
lapply(sample_list, function(obj) table(obj$doublet_finder))
```

Review the doublet fractions — a healthy droplet experiment typically
yields 1–8% doublets. If a sample looks unusually high, check the cell
count (very high cell loading increases multiplet rate) or QC metrics.

Filter doublets before integration:

``` r
sample_list <- lapply(sample_list, function(obj) {
  n_before <- ncol(obj)
  obj      <- subset(obj, doublet_finder == "Singlet")
  message(sprintf("  Dropped %d doublets; %d singlets remaining",
                  n_before - ncol(obj), ncol(obj)))
  obj
})
```

Apply a basic QC filter while we’re here:

``` r
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

------------------------------------------------------------------------

## 5. Gene ID Check

Before merging samples from different sources it is good practice to
confirm that all objects use the same gene identifier format. For `ifnb`
this is trivially true since it all came from one object, but in real
projects samples often come from different pipelines or genome
annotations.

``` r
# Check the identifier type for a single object
detect_gene_id_type(sample_list[[1]])
#> [1] "HGNC symbol"

# Check consistency across the whole list
check_gene_ids_across_objects(sample_list)
#> All objects use: HGNC symbol
```

If this ever returns a mix (e.g., some objects have Ensembl IDs and
others have gene symbols), you will lose almost all genes at the merge
step. Fix identifier format before proceeding.

------------------------------------------------------------------------

## 6. QC Report and Filter

Rather than eyeballing violin plots and picking cutoffs by hand,
generate a full HTML QC report with recommended per-sample cutoffs and
apply them programmatically.

``` r
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

``` r
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

``` r
QCComparePlots(
  pre     = sample_list_pre,   # save this before overwriting
  post    = sample_list,
  metrics = c("nCount_RNA", "nFeature_RNA", "percent.mt")
)
```

This replaces the manual
`subset(obj, nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)`
step from earlier versions of this vignette.

------------------------------------------------------------------------

## 7. Merge and Integrate — `MergeSeurat()`

`MergeSeurat()` handles normalization, PCA, integration, clustering,
UMAP, and (optionally) marker detection in one call. We use Harmony
integration here, which is fast and works well for this dataset.

``` r
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

``` r
# UMAP colored by cluster
DimPlot(integrated, label = TRUE, repel = TRUE)

# Verify that CTRL and STIM cells are interleaved (good integration)
DimPlot(integrated, group.by = "stim", cols = c("steelblue", "tomato"))

# Cells per cluster
table(Idents(integrated))
```

If CTRL and STIM cells form separate blobs by condition in the UMAP, the
integration did not fully correct for the stimulation effect. In that
case, try increasing `max_dims`, adjusting `cluster_resolution`, or
switching to `RPCAIntegration` with `k_anchor = 20, k_weight = 100`.

------------------------------------------------------------------------

## 8. Integration Quality — `BatchEffectQC()`

`BatchEffectQC()` quantifies how well the integration mixed CTRL and
STIM (or across donors) versus preserving biology. Run it on both the
pre-integration `pca` and post-integration `harmony` reductions to see
the effect numerically:

``` r
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

Ideally `knn_mixing` ratio approaches 1 (perfectly random neighborhoods
across conditions), and `batch_asw` drops toward 0 (cells no longer
cluster by condition). Once you have cell-type labels, re-run with
`celltype_col` set to check that biological structure was preserved
(`knn_purity` and `celltype_asw` should stay high).

------------------------------------------------------------------------

## 9. Annotate Cell Types — Three Approaches

`SingleCellTools` v2.4 offers three complementary paths for cell-type
annotation. For PBMC-style well-studied tissues, Azimuth is usually the
fastest and most accurate; marker-based approaches are more transparent
and don’t require internet access.

### 9.1 Marker Dot Plots — `MarkerPlot()` / `MarkerPctPlot()`

Now we assign cell-type labels to clusters using canonical PBMC markers.
`MarkerPlot()` builds an annotated dot plot, groups the genes by cell
type, and clusters the identities by expression correlation so similar
populations appear next to each other.

``` r
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

`MarkerPlot()` will automatically drop any gene absent from the assay or
with zero expression, and report those drops as messages — no manual
filtering needed.

**Auto-sizing for large panels.** The example above sets
`label.fontsize` explicitly. Leave it (and `axis_text_size`) at their
`NULL` default instead and `MarkerPlot()` scales both, plus a suggested
figure size, down as the gene count grows — handy for marker panels much
bigger than this one. `save_path` saves directly at that computed size:

``` r
MarkerPlot(integrated, pbmc_markers, save_path = "markerplot_ifnb.pdf")
```

### Assign cell-type labels

After reviewing the dot plot, rename clusters:

``` r
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

------------------------------------------------------------------------

### 9.2 Cluster-level Marker Scoring — `AnnotateClusters()`

Skip the manual `RenameIdents` step by having `AnnotateClusters` score
each cluster on every marker set with UCell and pick the winner:

``` r
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

For PBMCs, the CellTypist immune models give consistent, granular labels
without needing to supply a reference dataset. `Immune_All_Low.pkl` is
the low-resolution PBMC-friendly model; `Immune_All_High.pkl` is the
fine-resolution variant.

``` r
integrated <- AnnotateWithReference(
  integrated,
  method          = "celltypist",
  model           = "Immune_All_Low.pkl",
  majority_voting = TRUE,
  min_score       = 0.5
)

table(integrated$predicted_cell_type)

# Compare to your marker-based labels
table(marker    = integrated$predicted_cell_type,
      annotate  = integrated$domain)      # from 9.2 AnnotateClusters
```

Alternative backends work with the same interface — supply your own
labeled reference:

``` r
# scANVI: semi-supervised, best when reference and query come from
# different technologies
AnnotateWithReference(integrated, method = "scanvi",
                      reference = pbmc_ref, ref_label_col = "cell_type",
                      batch_col = "sample_id")

# scmap: R-native, no Python required
AnnotateWithReference(integrated, method = "scmap",
                      reference = pbmc_ref, ref_label_col = "cell_type",
                      scmap_method = "cluster")
```

For the rest of the vignette we’ll use `cell_type` — pick whichever of
the three above you trust most and assign:

``` r
integrated$cell_type <- integrated$predicted_cell_type   # from CellTypist
Idents(integrated) <- integrated$cell_type
```

------------------------------------------------------------------------

## 10. Flag Gene-Positive Cells — `AddGenePositivity()` / `PlotGenePositivity()`

`AddGenePositivity()` adds a logical metadata column for each gene
indicating whether a cell expresses it above a threshold. This is useful
for gating, subsetting, or computing co-expression fractions.

``` r
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

``` r
# Isolate CD14+ monocytes for downstream re-analysis
mono <- subset(integrated, CD14_pos == TRUE)
ncol(mono)
```

Visualize positivity across cell types:

``` r
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

### 10.1 `GenePositivityAnalysis()` and `GenePositivityEstimationPlot()`

The positivity summary above is a per-cell-type snapshot; it doesn’t say
whether a gene’s positivity rate actually shifted between CTRL and STIM.
`GenePositivityAnalysis()` computes per-sample positivity rates —
stratified by `cell_type` here — plus an optional chi-square/Fisher test
across `stim`, the same `sample_id`/`stim` design
`CompositionAnalysis()` (12.1) uses – and the same caveat: the test
pools cells across donors within each condition, so it’s not a valid
replicate-level test (it emits a `warning()` saying so). There’s no
propeller-style sample-level test for gene positivity here, so
`GenePositivityEstimationPlot()` below is the result to trust for an
effect that actually reflects the donor-level replicates:

``` r
gpa <- GenePositivityAnalysis(
  integrated,
  genes         = c("CD14", "GNLY"),
  sample_col    = "sample_id",
  condition_col = "stim",
  group_col     = "cell_type",
  test          = "chisq"
)
gpa$proportions               # sample, gene, group, n_pos, n_total, prop_pos, condition
gpa$test[["CD14 | CD14 Mono"]]   # chi-square result for that gene x cell-type combination (pooled-cell caveat above)
```

As with `CompositionAnalysis()`/`CompositionEstimationPlot()`, the test
above answers “is there a difference”; `GenePositivityEstimationPlot()`
answers “by how much, with what uncertainty” using per-sample positivity
rates and a bootstrap 95% CI:

``` r
# CD14 positivity in CD14+ monocytes, CTRL vs STIM
GenePositivityEstimationPlot(gpa, genes = "CD14", group_levels = "CD14 Mono",
                             idx = c("CTRL", "STIM"))

# Every gene x cell-type combination present, as a named list of plots
plots <- GenePositivityEstimationPlot(gpa, idx = c("CTRL", "STIM"))
plots[["CD14 | CD14 Mono"]]

# Cohen's h -- the effect size designed specifically for comparing two proportions
GenePositivityEstimationPlot(gpa, genes = "CD14", group_levels = "CD14 Mono",
                             idx = c("CTRL", "STIM"), effect = "cohens_h")
```

#### How to read a `dabestr` estimation plot

This is the first estimation plot in this vignette, so it’s worth
spelling out what’s actually on it – the same layout applies to
`CompositionEstimationPlot()` (12.1) later on too
([`dabestr`](https://acclab.github.io/dabestr/) implements “estimation
statistics”: Ho et al. 2019, *Moving beyond P values: data analysis with
estimation graphics*, Nature Methods).

Each plot has two stacked panels sharing an x-axis (the two `idx`
conditions, here `CTRL`/`STIM`):

- **Top panel — raw data.** Every individual value feeding the
  comparison is plotted as a swarm along its condition’s column. That’s
  one point per *donor/sample* (`sample_id`), not per cell – the whole
  point of these functions is to show the per-sample values a bare
  p-value collapses into a single number.
- **Bottom panel — effect size.** The chosen `effect` (`mean_diff` by
  default) between STIM and CTRL, drawn as a single point with a
  vertical bar for its 95% bootstrap confidence interval (5000
  resamples, by default). A dashed horizontal line marks the reference
  condition’s (CTRL’s) value, so the effect-size point/CI can be read
  directly against it.

Because the CI comes from resampling the actual donors rather than a
parametric formula, it doesn’t assume normality – appropriate here since
`n` is the number of donors (often single digits), not the (much larger,
and non-independent) number of cells. There’s no p-value threshold to
eyeball; instead look at whether the CI includes zero (no effect) and
how wide it is – a wide CI with few donors is telling you the estimate
is uncertain, which a bare significant/non-significant p-value would
hide.

`effect` options, all available on every `*EstimationPlot()` function in
this package:

| `effect` | What it measures |
|----|----|
| `mean_diff` (default) | Test mean − reference mean, in the original units (e.g. proportion points). |
| `median_diff` | Same, using medians – less sensitive to one outlier donor. |
| `cohens_d` | Standardized mean difference (pooled-SD units); comparable across differently-scaled measurements. |
| `hedges_g` | `cohens_d` with a small-sample bias correction – prefer this over `cohens_d` with few donors per condition. |
| `cliffs_delta` | Non-parametric, rank-based (akin to a standardized Mann-Whitney effect) – robust to outliers and non-normal spread, no distributional assumptions. |
| `cohens_h` | Designed specifically for comparing two proportions – often the more principled choice for `prop_pos`/`prop` columns (bounded 0-1) rather than a raw mean difference. |

------------------------------------------------------------------------

## 11. Cell-Cycle Scoring — `assign_cell_cycle_phase()`

`assign_cell_cycle_phase()` uses UCell module scoring to assign S, G2M,
or G1 phase to each cell. It is more robust than Seurat’s
`CellCycleScoring` on datasets with lower sequencing depth because UCell
normalizes for library size internally.

``` r
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

If cycling cells cluster separately from their resting counterparts and
you want to remove this confound, re-run `MergeSeurat()` with cell-cycle
scores added to `to_regress`:

``` r
integrated <- MergeSeurat(
  seurat_objects = sample_list,
  to_regress     = c("percent.mt", "S.Score", "G2M.Score"),
  use_SCT        = TRUE,
  integration    = "HarmonyIntegration",
  new_reduction  = "harmony",
  cluster_resolution = 0.5
)
```

------------------------------------------------------------------------

## 12. Cell-Type Composition — `CellComposition()` / `CompositionalTest()`

A natural follow-up question after integration is whether the
proportions of each cell type differ between CTRL and STIM.

``` r
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

Test statistical significance of the composition shift using propeller
(empirical-Bayes-moderated ANOVA on arcsin-sqrt transformed
proportions):

``` r
comp_test <- CompositionalTest(
  integrated,
  cluster_col   = "cell_type",
  sample_col    = "sample_id",
  condition_col = "stim",
  method        = "auto"        # propeller if available, else betareg, else wilcox
)
subset(comp_test, padj < 0.05)
```

The IFN-β stimulated condition typically shows a shift in monocyte
proportions and upregulation of ISGs across all cell types — a
well-characterized response that makes `ifnb` a good sanity check that
your integration and composition-testing pipelines are working
correctly.

### 12.1 `CompositionAnalysis()` and `CompositionEstimationPlot()`

`CompositionalTest()` above answers “does composition differ between
CTRL and STIM” with a p-value per cell type. `CompositionAnalysis()`
computes the same long-format counts/proportions but is the function
whose output `CompositionEstimationPlot()` consumes to instead show *how
big* that shift is, with a bootstrap 95% confidence interval, rather
than just a p-value – see 10.1 for how to read that plot (raw per-donor
swarm + effect size/CI panels) and what each `effect` option means.
Since `ifnb` has real CTRL/STIM replicates (multiple donors per
condition via `sample_id`), this is a genuine estimation-plot demo, not
a syntax stub. Note that `comp2$test` below pools cells across donors
rather than testing on the per-donor proportions `CompositionalTest()`
above uses, so it emits a `warning()` and is a much rougher check than
the propeller/wilcox result already computed above — included here only
because `test = "chisq"` is a `CompositionAnalysis()` option, not as the
test to actually report:

``` r
comp2 <- CompositionAnalysis(
  integrated,
  group_col     = "cell_type",
  sample_col    = "sample_id",
  condition_col = "stim",
  test          = "chisq"
)
comp2$test   # chi-square result across all cell types at once (pooled-cell caveat above)

# Monocytes are the cell type most affected by IFN-β stimulation --
# CompositionEstimationPlot() shows the effect size directly instead of a p-value
CompositionEstimationPlot(comp2, group_levels = "CD14 Mono",
                          idx = c("CTRL", "STIM"))

# Every cell type present, as a named list of plots
plots <- CompositionEstimationPlot(comp2, idx = c("CTRL", "STIM"))
plots[["CD14 Mono"]]

# Cohen's h -- the effect size designed specifically for comparing two proportions
CompositionEstimationPlot(comp2, idx = c("CTRL", "STIM"),
                          group_levels = "CD14 Mono", effect = "cohens_h")
```

------------------------------------------------------------------------

## 13. Differential Expression — `PseudobulkDE()` / `PlotVolcano()`

The correct way to identify STIM-vs-CTRL DE genes is per-cell-type
pseudobulk DE, not `FindMarkers`. `PseudobulkDE()` aggregates per
(sample, cell_type) and runs DESeq2:

``` r
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

``` r
PlotVolcano(de_all$CD14_Mono$results,
            fc_threshold = 1,
            p_threshold  = 0.05,
            top_n        = 20,
            label_genes  = c("ISG15", "IFI6", "MX1", "IFIT1", "IFIT3"))
```

Compare the DE profiles across two cell types (do they share the same
ISGs?):

``` r
cmp <- CompareMarkers(
  de_all$CD14_Mono$results,
  de_all$CD4_T$results,
  labels = c("CD14_Mono", "CD4_T")
)
cmp$overlap
cmp$plot                   # log2FC-vs-log2FC scatter, colored by category
```

------------------------------------------------------------------------

## 14. Cell-Cell Communication — `RunLIANA()`

Compare ligand-receptor interactions between conditions by running LIANA
on each stim/ctrl split, then intersecting the top interactions:

``` r
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

For a differential ligand-receptor pipeline that formally tests
condition effects, LIANA’s own `liana_deconstruct` or CellChat’s
differential analysis (14.1 below) are the recommended paths —
`RunLIANA()` is a starting point.

### 14.1 `RunCellChat()`

`RunCellChat()` wraps the full `CellChat` pipeline (`createCellChat`
through `aggregateNet`) into one call, for one Seurat subset at a time —
the same CTRL/STIM-split pattern as `RunLIANA()` above, but giving you
`CellChat`’s own signaling-pathway-level output and comparison plots
instead of a flat interaction table. `ifnb` is human, so pass
`species = "human"`; `group_col` points at the same `cell_type` labels
used throughout this vignette (`RunCellChat()`’s default `group_col` is
`"signaling_group"`, so it must be overridden here).

``` r
cc_ctrl <- RunCellChat(
  subset(integrated, stim == "CTRL"),
  label     = "CTRL",
  group_col = "cell_type",
  species   = "human"
)

cc_stim <- RunCellChat(
  subset(integrated, stim == "STIM"),
  label     = "STIM",
  group_col = "cell_type",
  species   = "human"
)
```

`cc_ctrl`/`cc_stim` are ordinary `CellChat` objects, so any `CellChat`
visualization works directly on them. Merging the two lets `CellChat`’s
own comparison functions show what actually shifted between conditions:

``` r
merged <- CellChat::mergeCellChat(list(CTRL = cc_ctrl, STIM = cc_stim),
                                  add.names = c("CTRL", "STIM"))

# Total interactions/strength, CTRL vs STIM
CellChat::compareInteractions(merged, show.legend = FALSE, group = c(1, 2))

# Which signaling pathways gained/lost strength in STIM
CellChat::rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)

# Bubble plot of specific source -> target interactions, both conditions side by side
CellChat::netVisual_bubble(merged, sources.use = "CD14 Mono", targets.use = "CD4 T",
                           comparison = c(1, 2))
```

`nboot` (bootstrap iterations for `computeCommunProb()`’s significance
test, default 25) and `min_cells` (minimum group size to keep an
interaction, default 10) are the two knobs most worth raising for a real
analysis — the low defaults keep a first pass fast, not to be used as
final results:

``` r
cc_stim_full <- RunCellChat(
  subset(integrated, stim == "STIM"),
  label     = "STIM",
  group_col = "cell_type",
  species   = "human",
  nboot     = 100
)
```

------------------------------------------------------------------------

## 15. Save with Provenance — `SaveWithProvenance()`

Save the annotated object with a JSON provenance sidecar so downstream
readers can inspect what analysis it went through without loading the
`.rds`:

``` r
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

For a quick summary of the object’s state after all this work:

``` r
CellSuiteSummary(integrated,
                 cluster_col = "cell_type",
                 sample_col  = "sample_id",
                 top_markers = 5)
```

------------------------------------------------------------------------

## 16. Session Info

``` r
sessionInfo()
```

Key packages used in this vignette:

| Package | Role |
|----|----|
| `SingleCellTools` | Doublet calling, QC report/filter, integration, annotation, DE, composition, LIANA, provenance |
| `Seurat` | Core data structure and analysis functions |
| `SeuratData` | `ifnb` dataset |
| `DoubletFinder` | Underlying doublet detection engine |
| `UCell` | Module scoring for annotation and cell-cycle phase |
| `harmony` | Batch correction (called via `IntegrateLayers`) |
| `DESeq2` | Pseudobulk differential expression |
| `speckle` | Propeller composition test |
| `CellTypist` (Python, via `reticulate`) | Default reference-based annotation |
| `scmap` | R-native alternative reference-based annotation |
| `liana` | Ligand-receptor consensus scoring |
| `patchwork` | QCComparePlots grid layout |
| `ks` | 2D KDE for `PlotFeatureDensity` |
| `jsonlite` | Provenance sidecar |
| `dplyr` / `ggplot2` | Data wrangling and plotting |

------------------------------------------------------------------------

*For questions or bug reports, open an issue at
[github.com/gevensen95/SingleCellTools/issues](https://github.com/gevensen95/SingleCellTools/issues).*
