# SingleCellTools Vignette & Tutorial

**Package:** `SingleCellTools`  
**Author:** K. Garrett Evensen, PhD — Bioinformatics Analyst II, The Salk Institute  
**Version:** 2.4  
**Source:** [github.com/gevensen95/SingleCellTools](https://github.com/gevensen95/SingleCellTools)

---

## Table of Contents

1. [Overview](#1-overview)
2. [Installation](#2-installation)
3. [Core scRNA-seq Workflow](#3-core-scrna-seq-workflow)
   - 3.1 [Loading Data — `CreateRNAObjects()`](#31-loading-data--creaternаobjects)
   - 3.2 [Interactive QC Filtering — `CreateRNAObjectsFilter()`](#32-interactive-qc-filtering--creatернаobjectsfilter)
   - 3.3 [Doublet Detection — `calldoublet()`](#33-doublet-detection--calldoublet)
   - 3.4 [Merging and Integration — `MergeSeurat()`](#34-merging-and-integration--mergeseurat)
   - 3.5 [Sub-clustering — `SubsetAndRecluster()`](#35-sub-clustering--subsetandrecluster)
   - 3.6 [Cell Annotation — `MarkerPlot()` / `MarkerPctPlot()`](#36-cell-annotation--markerplot--markerpctplot)
   - 3.7 [Marker-based Annotation — `AnnotateClusters()`](#37-marker-based-annotation--annotateclusters)
   - 3.8 [Reference-based Annotation — `AnnotateWithReference()`](#38-reference-based-annotation--annotatewithreference)
   - 3.9 [Batch Marker Discovery — `FindMarkersList()`](#39-batch-marker-discovery--findmarkerslist)
   - 3.10 [Gene Positivity — `AddGenePositivity()` / `PlotGenePositivity()`](#310-gene-positivity--addgenepositivity--plotgenepositivity)
4. [One-Shot Pipeline — `CreateAndIntegrateRNA()`](#4-one-shot-pipeline--createandintegrаterna)
5. [QC Reporting and Filtering](#5-qc-reporting-and-filtering)
   - 5.1 [`GenerateQCReport()`](#51-generateqcreport)
   - 5.2 [`ApplyQCFilters()`](#52-applyqcfilters)
   - 5.3 [`QCComparePlots()`](#53-qccompareplots)
6. [Differential Expression](#6-differential-expression)
   - 6.1 [`PseudobulkDE()`](#61-pseudobulkde)
   - 6.2 [`PlotVolcano()`](#62-plotvolcano)
   - 6.3 [`CompareMarkers()`](#63-comparemarkers)
7. [Composition Analysis](#7-composition-analysis)
   - 7.1 [`CellComposition()`](#71-cellcomposition)
   - 7.2 [`CompositionalTest()`](#72-compositionaltest)
8. [Cell-Cell Communication — `RunLIANA()`](#8-cell-cell-communication--runliana)
9. [Trajectory Analysis — `PseudotimeWrapper()`](#9-trajectory-analysis--pseudotimewrapper)
10. [Integration Quality — `BatchEffectQC()`](#10-integration-quality--batcheffectqc)
11. [Visualization — `PlotFeatureDensity()`](#11-visualization--plotfeaturedensity)
12. [Spatial Workflows](#12-spatial-workflows)
    - 12.1 [Visium — `CreateVisiumObjects()` + `EdgeDetectionVisium()`](#121-visium--createvisiumobjects--edgedetectionvisium)
    - 12.2 [Xenium — `LoadXenium2()`](#122-xenium--loadxenium2)
    - 12.3 [Parse Biosciences — `MakeParseObj()`](#123-parse-biosciences--makeparseobj)
    - 12.4 [FOV Edge / Tissue Hole Detection](#124-fov-edge--tissue-hole-detection)
    - 12.5 [Neighborhood Enrichment — `NeighborhoodEnrichment()`](#125-neighborhood-enrichment--neighborhoodenrichment)
    - 12.6 [Niche Co-expression — `NicheCoExpress()`](#126-niche-co-expression--nicheco express)
    - 12.7 [Visium Deconvolution — `RunRCTD()`](#127-visium-deconvolution--runrctd)
13. [scATAC-seq — `CreateATACObjects()`](#13-scatac-seq--createatacobjects)
14. [scanpy Interoperability — `ToAnnData()` / `FromAnnData()`](#14-scanpy-interoperability--toanndata--fromanndata)
15. [Reproducibility](#15-reproducibility)
    - 15.1 [`SaveWithProvenance()`](#151-savewithprovenance)
    - 15.2 [`CellSuiteSummary()`](#152-cellsuitesummary)
16. [Utility Functions](#16-utility-functions)
    - 16.1 [Gene ID Diagnostics](#161-gene-id-diagnostics)
    - 16.2 [Cell-Cycle Scoring](#162-cell-cycle-scoring)
    - 16.3 [Spatial Niche Assays](#163-spatial-niche-assays)
    - 16.4 [Spatial Polygon Tools](#164-spatial-polygon-tools)
    - 16.5 [GO Term Children](#165-go-term-children)
17. [Tips and Common Pitfalls](#17-tips-and-common-pitfalls)
18. [Session Info](#18-session-info)

---

## 1. Overview

`SingleCellTools` is an opinionated set of R wrapper functions that reduce the boilerplate of common single-cell and spatial analysis tasks. Rather than replacing Seurat or Signac, it wraps them — consolidating the five-to-ten function calls you write for every project into single, well-parametrized entry points.

The package covers:

- **Loading**: 10x Genomics, Parse Biosciences, Visium, Xenium, and scATAC-seq into Seurat objects
- **QC**: doublet calling (DoubletFinder), automated HTML QC reports with recommended cutoffs, one-shot filter application, pre/post filter comparison plots
- **Integration**: Harmony / RPCA / CCA / JointPCA with normalization, PCA, clustering, and UMAP in one call; quantitative batch-effect QC
- **Annotation**: marker-list scoring (UCell), SingleR reference correlation, Azimuth reference projection, cluster-level majority voting with `min_score` / `min_margin` thresholds; tumor-mode and Visium-safe defaults
- **DE and downstream**: pseudobulk DE (DESeq2), volcano plots, marker-set overlap comparison, cell-type composition testing (propeller / beta regression / Wilcoxon), ligand-receptor inference (LIANA)
- **Trajectory**: pseudotime inference (Slingshot)
- **Spatial**: edge / hole detection, angular-gap and bounding-box outlier calling, neighborhood enrichment, niche co-expression, RCTD deconvolution, spatial polygon ROI tools
- **Interop**: bidirectional Seurat ↔ AnnData conversion for scanpy compatibility
- **Reproducibility**: `.rds` + JSON provenance sidecar, one-command project summary card
- **Visualization**: annotated dot plots with automatic filtering, percent-positive heatmaps, gene-positivity overlays, 2D KDE feature density (Nebulosa-style)

---

## 2. Installation

### Core package

```r
# install.packages("remotes")
remotes::install_github("gevensen95/SingleCellTools")
```

### DoubletFinder (GitHub-only dependency)

```r
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
```

### Bioconductor dependencies

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  # Core scRNA-seq / spatial
  "EnsDb.Mmusculus.v79",  # mouse gene annotations
  "glmGamPoi",            # fast SCTransform fitting
  "GO.db",                # GO term lookups
  "UCell",                # module scoring
  "Signac",               # scATAC-seq
  # Differential expression and composition
  "DESeq2",
  "SummarizedExperiment",
  "S4Vectors",
  "edgeR",
  "speckle",              # propeller composition test
  "limma",
  # Trajectory
  "slingshot",
  # Reference-based annotation
  "SingleR",              # optional; used by AnnotateClusters(method="singler")
  # scanpy interop
  "zellkonverter"
))
```

### GitHub-only dependencies

```r
# Ligand-receptor inference
remotes::install_github("saezlab/liana")

# scmap for the R-native branch of AnnotateWithReference
BiocManager::install("scmap")

# CellTypist (Python default backend of AnnotateWithReference)
reticulate::py_install("celltypist")
# Or: pip install celltypist

# scANVI (optional Python backend of AnnotateWithReference)
reticulate::py_install("scvi-tools")
# Or: pip install scvi-tools

# Visium deconvolution
remotes::install_github("dmcable/spacexr")
```

### CRAN dependencies added in v2.4

```r
install.packages(c(
  "patchwork", "ks", "MASS", "cluster",
  "jsonlite", "betareg"
))
```

### Optional dependencies

```r
install.packages("sf")     # required only for get_cells_in_polygon()
install.packages("anndata"); reticulate::install_python()  # for FromAnnData Python reader
```

### Load the package

```r
library(SingleCellTools)
```

On load, `SingleCellTools` auto-attaches all of its `Imports` so functions like `%>%`, `ggplot()`, `filter()`, etc. are usable without namespace prefixes. Missing optional dependencies are reported as a single startup message and their functions become inert.

---

## 3. Core scRNA-seq Workflow

The typical pipeline (v2.4) looks like this:

```
CellRanger outputs
      ↓
CreateRNAObjects()             ← load, compute %mt, call doublets
      ↓
GenerateQCReport()             ← recommended cutoffs (HTML + CSV sidecar)
      ↓
ApplyQCFilters()               ← per-metric + doublet drop
      ↓
QCComparePlots()               ← verify effect of filtering
      ↓
MergeSeurat()                  ← normalize, integrate, cluster, UMAP
      ↓
BatchEffectQC()                ← integration diagnostics (optional)
      ↓
AnnotateClusters()             ← marker or SingleR labels
   or AnnotateWithReference()  ← Azimuth reference projection
      ↓
SubsetAndRecluster()           ← per-cell-type sub-clustering (optional)
      ↓
PseudobulkDE()                 ← condition contrasts
      ↓
PlotVolcano() / PlotFeatureDensity() / CellComposition() / RunLIANA()
      ↓
SaveWithProvenance()           ← .rds + JSON sidecar
```

---

### 3.1 Loading Data — `CreateRNAObjects()`

`CreateRNAObjects()` accepts a vector of directory paths, one per sample. Each directory can be:

- the raw CellRanger output folder (contains `outs/filtered_feature_bc_matrix/`)
- a `filtered_feature_bc_matrix/` folder directly
- a folder with loose matrix files (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`)
- a folder containing a filtered `.h5` file

The function auto-detects the layout so you do not need to adapt your paths.

```r
library(SingleCellTools)

# Paths to per-sample CellRanger output folders
samples <- c(
  "data/vehicle_rep1",
  "data/vehicle_rep2",
  "data/drugA_rep1",
  "data/drugA_rep2"
)

seurat_list <- CreateRNAObjects(
  data_dirs    = samples,

  # Basic QC thresholds applied at object creation
  cells        = 3,          # a gene must appear in ≥3 cells to be kept
  features     = 200,        # a cell must express ≥200 genes to be kept

  # Mitochondrial read percentage
  mt_pattern   = "^mt-",     # mouse mitochondrial genes; use "^MT-" for human

  # Optional: add a Treatment column to metadata
  treatment    = c("Vehicle", "Vehicle", "DrugA", "DrugA"),

  # Optional: give the list elements readable names
  object_names = c("Veh1", "Veh2", "Drug1", "Drug2"),

  # Doublet detection (runs calldoublet() on every sample)
  run_doublet_finder        = TRUE,
  doublet_normalization     = "LogNormalize",   # or "SCT"
  doublet_vars_to_regress   = "percent.mt",
  doublet_cluster_resolution = 0.1,

  # Set to TRUE to immediately drop doublets from each object;
  # FALSE (default) keeps them so you can inspect the labels first
  filter_doublets = FALSE
)
```

After running, each element of `seurat_list` is a Seurat object with:

- standard count data in the `RNA` assay
- `percent.mt` in metadata
- `Treatment` in metadata (if supplied)
- `doublet_finder` ("Singlet"/"Doublet") in metadata (if `run_doublet_finder = TRUE`)
- QC boxplots printed to the active graphics device

#### Inspecting QC before filtering

```r
# Violin plots per sample — decide your nFeature and percent.mt cutoffs here
library(Seurat)
library(patchwork)

VlnPlot(seurat_list[["Veh1"]], features = c("nFeature_RNA", "percent.mt"))
```

#### Filtering based on QC

```r
seurat_list <- lapply(seurat_list, function(obj) {
  subset(obj,
    nFeature_RNA > 500  &
    nFeature_RNA < 6000 &
    percent.mt   < 20   &
    doublet_finder == "Singlet"
  )
})
```

**Key parameters**

| Parameter | Default | Notes |
|---|---|---|
| `cells` | `3` | Minimum cells a gene must appear in |
| `features` | `200` | Minimum genes per cell |
| `mt_pattern` | `"^mt-"` | Mouse pattern; use `"^MT-"` for human |
| `treatment` | `NULL` | Vector, one value per sample |
| `run_doublet_finder` | `TRUE` | Calls `calldoublet()` under the hood |
| `filter_doublets` | `FALSE` | Set `TRUE` to auto-drop doublets |

---

### 3.2 Interactive QC Filtering — `CreateRNAObjectsFilter()`

If you prefer to choose cutoffs interactively by eye at the console, use the `*Filter` variant. It shows QC plots and prompts you to enter thresholds:

```r
seurat_list <- CreateRNAObjectsFilter(
  data_dirs  = samples,
  mt_pattern = "^mt-",
  treatment  = c("Vehicle", "Vehicle", "DrugA", "DrugA")
)
```

The function will pause after each QC plot and ask for:

- Maximum `percent.mt`
- Minimum and maximum `nFeature_RNA`

Use this when you have a new tissue type or data modality and don't yet know reasonable cutoffs.

---

### 3.3 Doublet Detection — `calldoublet()`

`calldoublet()` wraps the full DoubletFinder workflow into a single call. It is invoked automatically by `CreateRNAObjects()`, but you can also run it directly on any Seurat object.

```r
# Run on a single object
obj_with_doublets <- calldoublet(
  obj                = seurat_list[["Veh1"]],
  samplenameIndex    = 1,              # unused internally; kept for compatibility
  normalization      = "LogNormalize", # or "SCT"
  vars.to.regress    = "percent.mt",
  cluster_resolution = 0.1
)

# Check the result
table(obj_with_doublets$doublet_finder)
#> Doublet Singlet
#>     312    5841
```

**What happens under the hood:**

1. Normalize → PCA
2. Auto-detect significant PCs using the cumulative variance method (>90% variance or elbow)
3. UMAP + clustering at `cluster_resolution`
4. pK parameter sweep (no ground truth) — picks the pK at maximum BCmvn
5. Estimate homotypic doublet proportion from cluster sizes
6. Run `doubletFinder()`
7. **Clean up** — strips normalized layers, SCT assay, PCA, and UMAP reductions from the returned object to keep it lightweight

The returned object carries only the original counts plus `doublet_finder` metadata, ready to be merged.

---

### 3.4 Merging and Integration — `MergeSeurat()`

`MergeSeurat()` takes the filtered list and handles everything through UMAP in one call.

```r
integrated <- MergeSeurat(
  seurat_objects = seurat_list,

  # Variables to regress during normalization
  to_regress = "percent.mt",

  # Normalization method
  use_SCT = TRUE,           # FALSE uses LogNormalize / FindVariableFeatures / ScaleData

  # Dimensionality
  max_dims       = 20,      # number of PCs to use (ignored if use_elbow_plot = TRUE)
  use_elbow_plot = FALSE,   # set TRUE to pick dims interactively at the console

  # Integration
  integration     = "HarmonyIntegration",  # see table below for alternatives
  new_reduction   = "harmony",

  # Clustering
  cluster_resolution = 0.3,

  # Save outputs
  save_rds_file = TRUE,
  file_name     = "lung_experiment",

  # Automatically run FindAllMarkers after clustering
  markers = TRUE,

  # Restrict to genes shared across all samples before merging
  common_genes_only = FALSE
)
```

On completion, the working directory will contain:

- `lung_experiment_merged_seurat_objects.rds` — the integrated Seurat object
- `dimplot_seurat_clusters.pdf` — UMAP colored by cluster
- `markers_all.csv` — full `FindAllMarkers` results
- `marker_plot.pdf` — dot plot of top 10 markers per cluster

#### Integration methods

| `integration` argument | `new_reduction` to set | Notes |
|---|---|---|
| `"HarmonyIntegration"` | `"harmony"` | Default; fast; good for most datasets |
| `"RPCAIntegration"` | `"integrated.rpca"` | Requires `k_anchor` and `k_weight` |
| `"CCAIntegration"` | `"integrated.cca"` | Requires `k_anchor` and `k_weight` |
| `"JointPCAIntegration"` | `"integrated.jpca"` | Requires `k_anchor` and `k_weight` |

```r
# Example with RPCA
integrated <- MergeSeurat(
  seurat_objects  = seurat_list,
  integration     = "RPCAIntegration",
  new_reduction   = "integrated.rpca",
  k_anchor        = 20,
  k_weight        = 100
)
```

#### Spatial data

For spatial assays pass `spatial = "Visium"` or `spatial = "Xenium"`. The function will use the appropriate assay name internally and keep images in sync.

```r
integrated_visium <- MergeSeurat(
  seurat_objects = visium_list,
  spatial        = "Visium",
  use_SCT        = TRUE,
  integration    = "HarmonyIntegration",
  new_reduction  = "harmony"
)
```

---

### 3.5 Sub-clustering — `SubsetAndRecluster()`

Once you have a first-pass clustering and cell-type labels, it's common to want a finer look at one or a few types (e.g. the T-cell compartment). `SubsetAndRecluster()` chains subsetting → cleanup → optional re-normalization → PCA → optional integration → UMAP → clustering, with sensible defaults for Seurat v5's split-layer conventions.

```r
tcell <- SubsetAndRecluster(
  integrated,
  metadata_col   = "cell_type",
  metadata_value = c("CD4 T", "CD8 T", "Treg"),
  assay          = "RNA",
  normalize      = "auto",            # detect missing data/scale.data layers
  integrate      = TRUE,
  integration_method = "HarmonyIntegration",
  dims           = 20,
  resolution     = 0.3
)
```

Key features:

- `normalize = "auto"` inspects the working assay's layers and runs `NormalizeData` + `FindVariableFeatures` + `ScaleData` (or `SCTransform` when `normalization_method = "SCT"`) only if `data` / `scale.data` are missing. Prevents the common `'arg' should be "counts"` error when passing a raw-count subset to PCA.
- Cell-empty and sparsely-detected-gene pruning after subsetting, so downstream SCT doesn't fail with "subscript out of bounds".
- Best-effort `JoinLayers` on `counts` layers so QC metadata (`nCount_<assay>`, `nFeature_<assay>`) is recomputed correctly for the retained cells.

You can also subset by `idents` or by an explicit cell-barcode vector:

```r
# By cluster id
SubsetAndRecluster(integrated, idents = c("0", "3", "7"))

# By explicit cells
SubsetAndRecluster(integrated, cells = my_barcodes)
```

---

### 3.6 Cell Annotation — `MarkerPlot()` / `MarkerPctPlot()`

`MarkerPlot()` builds an annotated dot plot from a two-column data frame of genes and their cell-type labels. Genes are automatically grouped, sorted, and labeled along the right edge of the plot. Identities are optionally clustered by correlation so similar cell types appear adjacent.

```r
# Define a marker panel as a data frame
markers <- data.frame(
  Gene = c(
    "Sftpc", "Sftpb", "Abca3",        # AT2
    "Ager",  "Hopx",  "Pdpn",         # AT1
    "Trp63", "Krt5",  "Krt14",        # Basal
    "Scgb1a1", "Scgb3a2",             # Club
    "Pecam1", "Cdh5", "Kdr",          # Endothelial
    "Col1a1", "Postn", "Acta2"        # Fibroblast
  ),
  CellType = c(
    rep("AT2", 3),
    rep("AT1", 3),
    rep("Basal", 3),
    rep("Club", 2),
    rep("Endothelial", 3),
    rep("Fibroblast", 3)
  )
)

# Set identities to the annotated column before plotting
Idents(integrated) <- "seurat_clusters"

p <- MarkerPlot(
  obj              = integrated,
  genes            = markers,
  assay            = "RNA",
  cluster          = TRUE,      # cluster identities by expression correlation
  show.annotations = TRUE,      # draw cell-type labels along right edge
  maxsize          = 5,         # maximum dot size
  label.fontsize   = 3,         # annotation label size
  margin_factor    = 0.5        # increase if labels are clipped
)
print(p)

# Save
ggsave("markerplot_lung.pdf", p, width = 12, height = 8)
```

**Automatic filtering.** Before rendering, `MarkerPlot()` silently drops:

1. Genes not present in the assay's feature set
2. Genes with zero expression across all cells
3. Genes with no variation across identities (these would render as a uniform grey row)

Each dropped gene is reported by a message so you always know what was excluded.

**Customizing identity order.** If you have already annotated clusters and want a specific order, set the identity factor levels before calling:

```r
Idents(integrated) <- factor(
  integrated$cell_type,
  levels = c("AT2", "AT1", "Club", "Basal", "Endothelial", "Fibroblast")
)
MarkerPlot(integrated, markers, cluster = FALSE)  # disable correlation clustering
```

---

#### Percent-positive companion — `MarkerPctPlot()`

`MarkerPctPlot()` isolates the percent-positive signal (percent of cells in each cluster expressing each gene) as either a heatmap tile or a dot plot. Where `MarkerPlot` encodes avg expression + percent positive together, `MarkerPctPlot` shows only the "how many cells detect this gene" signal, which is more readable when the avg-expression scaling would compress important variation.

```r
MarkerPctPlot2(
  integrated, cellID,
  style        = "tile",         # or "dot"
  colors       = c("white", "firebrick"),
  max_pct      = 80,             # spread color range over informative 0-80%
  label_offset = 1.2             # right-margin cell-type labels
)
```

Both `MarkerPlot` and `MarkerPctPlot` accept the same two-column `genes` data frame and share the internal `.pull_dotplot_data` / `.cluster_mat` helpers (in `dotplot-helpers.R`).

---

### 3.7 Marker-based Annotation — `AnnotateClusters()`

`AnnotateClusters()` transfers cell-type labels to each cluster using either UCell module scoring (`method = "marker"`) or SingleR reference correlation (`method = "singler"`). Cluster labels are assigned by averaging per-cell scores per cluster and picking the winner, with configurable `min_score` and `min_margin` thresholds for "Unknown" fallback.

```r
markers_list <- list(
  T_cell = c("CD3D", "CD3E", "CD8A", "CD4"),
  B_cell = c("MS4A1", "CD79A", "CD19"),
  NK     = c("NKG7", "GNLY", "KLRD1"),
  Mono   = c("CD14", "LYZ", "FCGR3A")
)

integrated <- AnnotateClusters(
  integrated,
  method              = "marker",
  markers             = markers_list,
  cluster_col         = "seurat_clusters",
  new_col             = "predicted_cell_type",
  filter_nonspecific  = TRUE,      # drop broadly-expressed "markers"
  min_score           = 0.1,
  min_margin          = 0.05
)
table(integrated$predicted_cell_type)
```

**Tumor-friendly defaults.** In tumor samples, canonical markers get compressed toward a proliferating, EMT-mixed dedifferentiated state. Pass `tumor_mode = TRUE` to disable `filter_nonspecific` and tighten `high_relative_to_max`, and to trigger a warning that labels are approximate.

```r
obj <- AnnotateClusters(obj, markers = markers_list, tumor_mode = TRUE)
```

**Visium-friendly usage.** Each Visium spot mixes multiple cell types; the winner-takes-all classifier hides minority signal. Use `return_scores = "cluster"` to inspect the full per-cluster score matrix, `min_detection_frac` to drop signatures whose markers barely register, and `filter_nonspecific = FALSE`:

```r
res <- AnnotateClusters(
  visium, markers = markers_list,
  filter_nonspecific  = FALSE,
  min_detection_frac  = 0.05,
  return_scores       = "cluster"
)
res$obj               # annotated
res$scores            # cluster x cell-type matrix; heatmap this
```

---

### 3.8 Reference-based Annotation — `AnnotateWithReference()`

Reference-based annotation with a choice of three backends. All write the predicted label to `new_col` and (where available) a confidence score to `<new_col>_score`.

**Default: CellTypist.** A pre-trained logistic-regression classifier with a large model zoo (immune, gut, lung, kidney, tumor, developmental, and more). No user-supplied reference needed. Requires Python via `reticulate` + `pip install celltypist`.

```r
integrated <- AnnotateWithReference(
  integrated,
  method          = "celltypist",
  model           = "Immune_All_Low.pkl",   # or Immune_All_High, Lung, etc.
  majority_voting = TRUE,                    # smoothing via over-clustering
  min_score       = 0.5
)
table(integrated$predicted_cell_type)
```

**scANVI** — semi-supervised VAE reference mapping via `scvi-tools`. Best when reference and query come from different technologies or protocols. Requires a labeled reference Seurat / `.h5ad`.

```r
integrated <- AnnotateWithReference(
  integrated,
  method        = "scanvi",
  reference     = pbmc_ref_seurat,
  ref_label_col = "cell_type",
  batch_col     = "orig.ident",
  n_latent      = 30
)
```

**scmap** — R-native (Bioconductor), no Python required. Two variants: `scmap-cluster` (fast, cluster-level projection) and `scmap-cell` (product-quantization nearest neighbor). Requires a labeled reference.

```r
integrated <- AnnotateWithReference(
  integrated,
  method        = "scmap",
  reference     = pbmc_ref_seurat,
  ref_label_col = "cell_type",
  scmap_method  = "cluster",     # or "cell"
  threshold     = 0.5,
  n_features    = 500
)
```

Complements `AnnotateClusters()` — use CellTypist / scANVI / scmap for the coarse pass, then optionally refine sub-populations with marker-based scoring on a `SubsetAndRecluster` output.

---

### 3.9 Batch Marker Discovery — `FindMarkersList()`

Runs `FindAllMarkers` across a list of Seurat objects, sets `Idents()` to a common cell-type / cluster column, and returns a named list of filtered marker tables (padj + log2FC thresholds applied). Useful across a cohort of samples processed separately.

```r
markers_by_sample <- FindMarkersList(
  seurat_list,
  idents_col       = "cell_ontology_class",
  padj_threshold   = 0.05,
  log2fc_threshold = 1
)

# Combine across samples
all_markers <- do.call(rbind, markers_by_sample)
```

---

### 3.10 Gene Positivity — `AddGenePositivity()` / `PlotGenePositivity()`

`AddGenePositivity()` adds a logical metadata column for each gene, flagging cells where expression exceeds a threshold. It accepts a single Seurat object or a list and returns the same shape.

```r
# Single object
integrated <- AddGenePositivity(
  seurat_objects = integrated,
  genes          = c("Sftpc", "Ager", "Trp63"),
  layer          = "counts",   # use "data" for log-normalized layer
  threshold      = 0,          # cells with counts > 0 are positive
  suffix         = "_pos"
)

# New metadata columns: Sftpc_pos, Ager_pos, Trp63_pos (logical)
head(integrated@meta.data[, c("Sftpc_pos", "Ager_pos", "Trp63_pos")])
```

For a list of objects (e.g., before merging), the function restricts to genes present in **every** object, ensuring a consistent metadata schema:

```r
seurat_list <- AddGenePositivity(
  seurat_objects = seurat_list,
  genes          = c("Sftpc", "Ager", "Trp63"),
  layer          = "counts"
)
```

You can then use positivity columns for subsetting or as covariates:

```r
# Subset to AT2 cells (Sftpc-positive)
at2 <- subset(integrated, Sftpc_pos == TRUE)

# Fraction positive per cluster
tapply(integrated$Sftpc_pos, Idents(integrated), mean)
```

`PlotGenePositivity()` summarizes those columns visually — three style modes:

```r
# Percent-positive bar chart per cluster, one bar per gene
PlotGenePositivity(integrated, c("Sftpc", "Ager", "Trp63"))

# Many genes → tile heatmap
PlotGenePositivity(integrated, big_gene_vec, style = "heatmap")

# Co-expression combinations (CD3D+/CD4+, CD3D+/CD4-, ...)
PlotGenePositivity(integrated, c("CD3D", "CD4"), style = "combo")

# For a list of samples, one facet per sample
PlotGenePositivity(seurat_list, c("Sftpc", "Ager"))
```

---

## 4. One-Shot Pipeline — `CreateAndIntegrateRNA()`

For straightforward experiments where you trust the defaults, `CreateAndIntegrateRNA()` chains loading → QC → merging → integration in a single function call.

```r
integrated <- CreateAndIntegrateRNA(
  data_dirs          = samples,
  mt_pattern         = "^mt-",
  treatment          = c("Vehicle", "Vehicle", "DrugA", "DrugA"),
  nFeature_min       = 500,
  nFeature_max       = 6000,
  percent_mt_max     = 20,
  run_doublet_finder = TRUE,
  filter_doublets    = TRUE,
  integration        = "HarmonyIntegration",
  cluster_resolution = 0.3,
  markers            = TRUE
)
```

This is convenient for quick exploratory runs, but `CreateRNAObjects()` + `MergeSeurat()` gives you more control over intermediate QC decisions.

---

## 5. QC Reporting and Filtering

Three functions form a self-contained QC pipeline: `GenerateQCReport()` picks cutoffs, `ApplyQCFilters()` applies them, and `QCComparePlots()` shows what changed.

### 5.1 `GenerateQCReport()`

Produces a self-contained HTML report of per-sample QC — per-metric distributions, log-skew-aware transformation, spatial QC when applicable, doublet call summaries — plus a table of **recommended filter cutoffs** built from `median ± mad_multiplier × MAD` per (sample, metric).

```r
GenerateQCReport(
  sample_list,
  output_file   = "qc/qc_report.html",
  metadata_cols = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  mad_multiplier = 3,          # 3 MAD = "outlier" threshold
  doublet_col   = "doublet_finder"
)
```

Alongside the HTML, it writes a machine-readable sidecar (`qc/qc_report_cutoffs.csv`) that `ApplyQCFilters()` picks up automatically.

### 5.2 `ApplyQCFilters()`

Reads the sidecar CSV (or accepts a data frame or the HTML path directly) and subsets each Seurat object to cells that pass every metric's `[suggest_lo, suggest_hi]` range. Doublet filtering is included by default via `filter_doublets = TRUE`.

```r
sample_list_filtered <- ApplyQCFilters(
  sample_list,
  cutoffs         = "qc/qc_report_cutoffs.csv",
  filter_doublets = TRUE,
  doublet_col     = "doublet_finder",
  doublet_value   = "Doublet",
  return_report   = TRUE          # for a per-sample retention table
)
sample_list_filtered$report        # sample x metric x pct_kept
```

Handy variants:

```r
# Only some metrics — leave doublet score / percent.hb alone
ApplyQCFilters(sample_list, cutoffs = "qc/qc_report_cutoffs.csv",
               metrics = c("nCount_RNA", "nFeature_RNA", "percent.mt"))

# Override a single sample's mito cutoff without editing the CSV
ApplyQCFilters(sample_list, cutoffs = "qc/qc_report_cutoffs.csv",
               override = list(sample1 = list(percent.mt = c(0, 15))))
```

### 5.3 `QCComparePlots()`

Pre/post violin overlays of QC metrics, one panel per metric, faceted by sample, with `n_before → n_after (pct%)` annotations. Confirms at a glance that the filtering did what you intended.

```r
QCComparePlots(
  pre     = sample_list,
  post    = sample_list_filtered$obj,
  metrics = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  log_y   = c("nCount_RNA", "nFeature_RNA")
)
```

---

## 6. Differential Expression

### 6.1 `PseudobulkDE()`

The statistically correct way to compare conditions across donors: aggregate per-cell counts into per-sample-per-cluster pseudobulk samples, then run DESeq2. Cell-level DE (Wilcoxon, MAST) treats each cell as independent and inflates type-I error.

```r
# Single cluster / cell type
res <- PseudobulkDE(
  integrated,
  sample_col    = "orig.ident",
  condition_col = "treatment",
  ident_1       = "DrugA",
  ident_2       = "Vehicle",
  group_by      = "cell_type",
  group_value   = "T cell"
)
head(res$results)                # DESeq2 results table
head(res$normalized_counts)      # size-factor-normalized pseudobulk counts
```

For multi-cluster analysis, pass `cluster_col` to loop over every cluster:

```r
res_all <- PseudobulkDE(
  integrated,
  sample_col    = "orig.ident",
  condition_col = "treatment",
  ident_1       = "DrugA",
  ident_2       = "Vehicle",
  cluster_col   = "cell_type"
)
res_all$T_cell$results            # per-cluster access
de_long <- do.call(rbind, lapply(res_all, `[[`, "results"))
```

Custom design formulas work for batch adjustment, three-level factors, and interaction models:

```r
res <- PseudobulkDE(
  integrated,
  sample_col    = "orig.ident",
  condition_col = "treatment",
  cluster_col   = "cell_type",
  design        = ~ batch + treatment,
  contrast      = c("treatment", "DrugA", "Vehicle")
)
```

Genes are filtered with `edgeR::filterByExpr` when available (recommended), with a row-sum fallback otherwise.

### 6.2 `PlotVolcano()`

Volcano plot with auto-detection of common column-name conventions. Accepts FindMarkers-style (`avg_log2FC`, `p_val_adj`) or PseudobulkDE-style (`log2FC`, `padj`) inputs.

```r
PlotVolcano(res$results,
            fc_threshold = 1,     # |log2FC| threshold
            p_threshold  = 0.05,
            top_n        = 20,    # label top-N by combined criterion
            label_genes  = c("MYC", "TP53"))
```

### 6.3 `CompareMarkers()`

Merge two DE result tables on gene and classify each as `shared_up / shared_down / opposite_sign / only_a / only_b / n.s. in both`. Includes a Fisher's exact test on the shared-vs-unique contingency and a log2FC-vs-log2FC scatter.

```r
res_drug    <- PseudobulkDE(obj, sample_col = "orig.ident",
                            condition_col = "treatment",
                            ident_1 = "drug", ident_2 = "vehicle")
res_disease <- PseudobulkDE(obj, sample_col = "orig.ident",
                            condition_col = "status",
                            ident_1 = "disease", ident_2 = "healthy")

cmp <- CompareMarkers(
  res_drug$results, res_disease$results,
  labels = c("drug_vs_veh", "disease_vs_healthy")
)
cmp$overlap                    # counts per category
cmp$fisher                     # Fisher's exact test
cmp$plot                       # log2FC-vs-log2FC scatter
head(cmp$merged)
```

---

## 7. Composition Analysis

### 7.1 `CellComposition()`

Per-sample per-cluster proportions in a tidy data frame + ggplot. Three plot styles: stacked bars per sample, boxplots per group per cluster, or line plots for time-series-style comparisons.

```r
comp <- CellComposition(
  integrated,
  cluster_col = "cell_type",
  sample_col  = "orig.ident",
  group_col   = "treatment",
  style       = "box"
)
comp$plot
head(comp$df)                  # per-sample per-cluster proportions
```

### 7.2 `CompositionalTest()`

Statistical test for shifts in composition between conditions. Prefers `speckle::propeller` (arcsin-sqrt or logit-transformed ANOVA — the current best practice), falls back to per-cluster beta regression (`betareg`), falls back further to per-cluster Wilcoxon with BH correction.

```r
comp_test <- CompositionalTest(
  integrated,
  cluster_col   = "cell_type",
  sample_col    = "orig.ident",
  condition_col = "treatment",
  transform     = "asin",
  method        = "auto"
)
subset(comp_test, padj < 0.05)
```

Each row has `cluster`, per-level `mean_prop_<level>`, `effect`, `stat`, `pvalue`, `padj`, `method`.

---

## 8. Cell-Cell Communication — `RunLIANA()`

Ligand-receptor inference via `liana::liana_wrap` + `liana::liana_aggregate`. Default runs LIANA's consensus set (Connectome, logFC, NATMI, SCA, CellPhoneDB) and aggregates their rankings — more robust than any single method.

```r
Idents(integrated) <- integrated$cell_type
lr <- RunLIANA(
  integrated,
  idents_col = "cell_type",
  method     = "consensus",
  resource   = "Consensus",
  min_cells  = 10
)
head(lr)

# Restrict to interactions FROM T cells
lr_t <- RunLIANA(integrated, idents_col = "cell_type",
                 source_cells = "T cell")

# Mouse via ortholog translation
lr_mouse <- RunLIANA(integrated, idents_col = "cell_type",
                     use_ortho = TRUE)
```

---

## 9. Trajectory Analysis — `PseudotimeWrapper()`

Slingshot wrapper for pseudotime inference. Fits lineage curves through the specified reduction using cluster labels as anchors, writes one metadata column per detected lineage.

```r
integrated <- PseudotimeWrapper(
  integrated,
  reduction     = "umap",
  dims          = 2,
  cluster_col   = "seurat_clusters",
  start_cluster = "3",             # optional root
  prefix        = "slingshot"
)
# New columns: slingshot_Lineage1, slingshot_Lineage2, ...
FeaturePlot(integrated, features = "slingshot_Lineage1")

# The SlingshotDataSet itself, for downstream analysis
sds <- integrated@misc$slingshot
```

---

## 10. Integration Quality — `BatchEffectQC()`

Quantifies how well integration has mixed batches while preserving biological structure. Meant for pre-vs-post integration comparison.

```r
pre  <- BatchEffectQC(integrated, reduction = "pca",
                      batch_col    = "orig.ident",
                      celltype_col = "cell_type")
post <- BatchEffectQC(integrated, reduction = "harmony",
                      batch_col    = "orig.ident",
                      celltype_col = "cell_type")

rbind(pre$summary, post$summary)
#      batch_asw celltype_asw knn_mixing knn_purity ...
# pre       0.32          0.28       0.61       0.72
# post      0.04          0.35       0.94       0.78
```

Interpretation:

- `batch_asw` closer to 0 = better mixing.
- `celltype_asw` higher = better biology preservation.
- `knn_mixing` closer to 1 = perfectly random neighborhoods across batches (perfect mixing).
- `knn_purity` closer to 1 = each cell's neighbors share its label (biology preserved).

Per-cluster diagnostics are returned in `$per_cluster` when `celltype_col` is set.

---

## 11. Visualization — `PlotFeatureDensity()`

Nebulosa-style weighted 2D KDE on any reduction. Much more readable than `FeaturePlot` for sparse markers or module scores.

```r
PlotFeatureDensity(
  integrated,
  features  = c("CD3D", "CD4", "CD8A"),
  reduction = "umap"
)

# Joint co-expression panel
PlotFeatureDensity(integrated, features = c("CD3D", "CD4"), joint = TRUE)

# From a metadata column (e.g. a module score)
PlotFeatureDensity(integrated, features = "module_score_Bcell")
```

Uses `ks::kde` when available, falls back to `MASS::kde2d`.

---

## 12. Spatial Workflows

### 12.1 Visium — `CreateVisiumObjects()` + `EdgeDetectionVisium()`

#### Loading

```r
visium_dirs <- c(
  "data/visium/sample1/outs",
  "data/visium/sample2/outs"
)

visium_list <- CreateVisiumObjects(
  data_dirs    = visium_dirs,
  object_names = c("S1", "S2")
)
```

#### Edge detection

Spots at the capture-area boundary, tissue edge, and tissue tears systematically have abnormal UMI counts and should be removed before analysis. `EdgeDetectionVisium()` runs four iterations of nearest-neighbor filtering to identify these spots.

```r
# Run on the first sample
edge_df <- EdgeDetectionVisium(
  coord_path = "data/visium/sample1/outs/spatial",  # folder containing tissue_positions_list.csv
  seurat.obj = visium_list[["S1"]],   # ensures barcode order matches
  search     = "radius",
  neighbors  = 7
)

# edge_df has columns: barcode, Filter, Filter2, Filter3, Filter4
# Each Filter column represents one iteration; Filter4 is the most conservative.
head(edge_df[, c("barcode", "Filter", "Filter4")])
```

Choosing which iteration to use is a judgment call — inspect the spatial distribution:

```r
library(Seurat)

# Add the filter column to metadata
visium_list[["S1"]] <- AddMetaData(
  visium_list[["S1"]],
  metadata = setNames(edge_df$Filter4, edge_df$barcode),
  col.name = "edge_filter"
)

# Visualize
SpatialDimPlot(visium_list[["S1"]], group.by = "edge_filter")

# Remove edge spots
visium_list[["S1"]] <- subset(visium_list[["S1"]], edge_filter == "Keep")
```

> **Tip:** Always run `EdgeDetectionVisium()` before `MergeSeurat()`. Edge spots inflate variance and distort integration.

#### Integrating Visium samples

```r
integrated_visium <- MergeSeurat(
  seurat_objects     = visium_list,
  spatial            = "Visium",
  use_SCT            = TRUE,
  to_regress         = "percent.mt",
  integration        = "HarmonyIntegration",
  new_reduction      = "harmony",
  cluster_resolution = 0.4,
  markers            = FALSE
)
```

---

### 12.2 Xenium — `LoadXenium2()`

```r
xenium_list <- lapply(
  c("data/xenium/sample1", "data/xenium/sample2"),
  LoadXenium2
)
names(xenium_list) <- c("X1", "X2")

# Merge and integrate
integrated_xenium <- MergeSeurat(
  seurat_objects = xenium_list,
  spatial        = "Xenium",
  integration    = "HarmonyIntegration",
  new_reduction  = "harmony",
  use_SCT        = TRUE
)
```

---

### 12.3 Parse Biosciences — `MakeParseObj()`

Parse pipeline output uses a `DGE_filtered/` folder layout. `MakeParseObj()` handles the sub-library structure and optional sample splitting.

```r
parse_list <- MakeParseObj(
  data_dirs    = c("data/parse/run1/DGE_filtered",
                   "data/parse/run2/DGE_filtered"),
  mt_pattern   = "^mt-",
  object_names = c("Run1", "Run2")
)
```

The returned list is a standard list of Seurat objects that can be passed directly to `MergeSeurat()`.

---

### 12.4 FOV Edge / Tissue Hole Detection

Two complementary functions for finding cells near boundaries in single-cell spatial data (Xenium, CosMx, etc.). Complements `EdgeDetectionVisium()` which handles Visium's hex-grid geometry.

`detect_fov_edges()` — flags cells near the outer edge of each FOV. Two methods:

```r
# Bounding-box method (default): fast, deterministic, good for roughly convex FOVs
obj <- detect_fov_edges(obj,
                        method       = "bbox",
                        bbox_factor  = 2,        # ~2 cell-widths of the box edge
                        n_iterations = 2,
                        label_col    = "edge_layer")

# Angular-gap method: catches concave edges / tears the bbox misses
obj <- detect_fov_edges(obj,
                        method         = "angular",
                        k              = 10,
                        gap_threshold  = 2 * pi / 3,
                        density_factor = 1)
```

Output: an integer metadata column (`0 = interior, 1 = outer ring, 2 = next ring, ...`).

`detect_tissue_holes()` — flags cells bordering internal gaps / tears using a 2D occupancy grid + 4-connected flood fill.

```r
obj <- detect_tissue_holes2(obj,
                            bin_size     = NULL,      # auto = 2.5 * median NN dist
                            min_hole_size = 4,
                            n_iterations  = 2,
                            label_col     = "hole_layer")
```

**Marker-gene exclusion** — biologically meaningful gaps (e.g. liver central veins expressing Glul) can be skipped:

```r
obj <- detect_tissue_holes2(obj,
                            exclude_gene       = "Glul",
                            sensitivity        = 0.75,   # data-adaptive quantile
                            exclude_gene_layer = "data")
```

Together, these give you a clean pre-analysis filter:

```r
obj <- obj[, obj$edge_layer == 0 & obj$hole_layer == 0]
```

---

### 12.5 Neighborhood Enrichment — `NeighborhoodEnrichment()`

Tests pairwise cell-type co-localization in space by comparing the observed frequency of each (source, target) k-NN pair against a permutation null. Optionally clusters cells into "niches" based on their neighborhood composition.

```r
enrich <- NeighborhoodEnrichment(
  obj,
  celltype_col = "cell_type",
  k            = 10,
  n_perm       = 200,
  assign_niches = TRUE,     # k-means clustering of neighborhood vectors
  k_niches     = 6
)
enrich$results               # source x target enrichment z-scores
obj$niche <- enrich$obj$niche
```

### 12.6 Niche Co-expression — `NicheCoExpress()`

Differential co-expression across niches / regions using the Manders overlap coefficient, useful for asking whether two ligand-receptor genes are co-detected more often in specific tissue zones than expected.

```r
co <- NicheCoExpress(
  obj,
  gene_a       = "Vegfa",
  gene_b       = "Kdr",
  niche_col    = "niche",
  layer        = "counts"
)
co$per_niche                 # Manders coefficient per niche + p-value
co$plot                      # heatmap
```

### 12.7 Visium Deconvolution — `RunRCTD()`

For Visium, per-spot cell-type calls are fundamentally the wrong output — each spot is 1-10 cells of mixed types. RCTD deconvolves each spot as a mixture using a reference single-cell dataset.

```r
visium <- RunRCTD(
  visium,
  reference    = pbmc_ref,
  celltype_col = "cell_type",
  mode         = "full",         # or "doublet" / "multi"
  max_cells_per_ref_celltype = 10000,
  n_cores      = 8
)

# Per-cell-type proportion columns
colnames(visium@misc$rctd_weights)
SpatialFeaturePlot(visium, features = c("rctd_T_cell", "rctd_B_cell"))

# Winner-takes-all label per spot (for quick visualization)
SpatialDimPlot(visium, group.by = "rctd_dominant")
```

`RunRCTD()` should be the primary annotation strategy for Visium; use `AnnotateClusters()` only as a quick sanity check.

---

## 13. scATAC-seq — `CreateATACObjects()`

```r
atac_list <- CreateATACObjects(
  data_dirs    = c("data/atac/sample1", "data/atac/sample2"),
  object_names = c("ATAC1", "ATAC2"),
  genome       = "mm10"          # or "hg38"
)
```

For interactive cutoff selection on ATAC QC metrics (fragment count, TSS enrichment, nucleosome signal), use:

```r
atac_list <- CreateATACObjectsFilter(
  data_dirs = c("data/atac/sample1", "data/atac/sample2"),
  genome    = "mm10"
)
```

---

## 14. scanpy Interoperability — `ToAnnData()` / `FromAnnData()`

Bidirectional conversion between Seurat objects and AnnData `.h5ad` files. Useful for hand-off to Python-based tools (scanpy, scVI, cellrank, CellPhoneDB) or for reading Python-produced data back into R.

```r
# Export to AnnData (writes .h5ad on disk)
ToAnnData(integrated, file = "results/pbmc.h5ad", assay = "RNA", layer = "data")

# Read AnnData back to a Seurat object
obj <- FromAnnData(file = "data/tabula_sapiens/Lung.h5ad",
                   reader = "python",   # "python" (via reticulate) or "R" (via zellkonverter)
                   assay = "RNA")
```

The `reader = "python"` path uses `reticulate` + `anndata` (through `basilisk` when available). The `reader = "R"` path uses `zellkonverter::readH5AD`. A tryCatch fallback also handles the case where `as.Seurat` fails on Seurat v5 by building the object manually from the underlying SCE assays, colData, and reducedDims.

Batch conversion of all `.h5ad` files in a directory:

```r
files <- list.files("data/tabula_sapiens", pattern = "\\.h5ad$",
                    full.names = TRUE)
objs  <- lapply(files, FromAnnData)
names(objs) <- sub("\\.h5ad$", "", basename(files))
```

---

## 15. Reproducibility

### 15.1 `SaveWithProvenance()`

Saves a Seurat object as `.rds` and, next to it, a `<name>_provenance.json` sidecar recording session info, package versions, `DefaultAssay`, layers per assay, reductions present, metadata column names, cell/gene counts, and — optionally — the git SHA and dirty flag of the calling script directory. Future-you (or a reviewer) can inspect the sidecar without loading the RDS.

```r
SaveWithProvenance(
  integrated,
  file    = "results/pbmc_annotated.rds",
  git_dir = getwd(),
  extra   = list(analyst = "K. Evensen", project = "study42")
)
# Writes results/pbmc_annotated.rds AND
#        results/pbmc_annotated_provenance.json
```

### 15.2 `CellSuiteSummary()`

One-command overview: cell/gene counts, assays, reductions, per-cluster and per-sample counts, standard QC metric summaries, and optionally top markers per cluster. Suitable for a README, a paper method section, or a hand-off document.

```r
s <- CellSuiteSummary(integrated,
                      cluster_col = "cell_type",
                      sample_col  = "orig.ident",
                      top_markers = 5)
print(s)
```

Sample output:

```
Seurat object summary
---------------------------------------------
  cells:         12345
  genes:         21050
  default assay: RNA
  all assays:    RNA, SCT
  reductions:    pca, harmony, umap
  clusters:      12 (col: cell_type)

Cluster sizes
 cluster n_cells  pct
 T cell     4210 34.1
 Mono       3055 24.7
 ...

QC metric summaries
        metric median   q25   q75 min   max
    nCount_RNA   2500  1120  4800 500 22300
  nFeature_RNA    950   540  1500 250  4820
    percent.mt   3.10  1.80  5.60   0  19.9
```

---

## 16. Utility Functions

### 16.1 Gene ID Diagnostics

Before merging objects from different sources, verify that gene identifiers are in the same format. Mixed symbol/Ensembl objects will silently lose most genes at merge.

```r
# Check a single object
detect_gene_id_type(seurat_list[["Veh1"]])
#> [1] "MGI symbol"   # or "HGNC symbol", "Ensembl", "RefSeq", "Entrez"

# Check across an entire list — catches mismatches before they cost you
check_gene_ids_across_objects(seurat_list)
```

`DetectGenes()` provides bulk detection and quality flags on a specific feature set:

```r
DetectGenes(integrated, genes = c("Sftpc", "Ager", "ENSMUG00000001234"))
```

---

### 16.2 Cell-Cycle Scoring

`assign_cell_cycle_phase()` uses UCell module scoring (rather than Seurat's `CellCycleScoring`) and is therefore more robust on datasets with lower sequencing depth.

```r
# Requires mouse S and G2M gene lists
# The function uses built-in gene sets from Seurat by default
integrated <- assign_cell_cycle_phase(integrated)

# New metadata columns: S.Score, G2M.Score, Phase
table(integrated$Phase)
#>  G1  G2M    S
#> 4021  812 1023

# Regress out cell cycle in a downstream normalization
integrated <- SCTransform(integrated, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))
```

---

### 16.3 Spatial Niche Assays

`BuildMultipleNicheAssays()` constructs a spatial neighborhood ("niche") assay across a list of objects, then clusters niches with mini-batch k-means. This is useful for identifying recurring spatial patterns across multiple tissue sections.

```r
visium_list <- BuildMultipleNicheAssays(
  seurat_list     = visium_list,
  neighbors       = 6,     # spots to include in each neighborhood
  k               = 8,     # number of niche clusters
  assay           = "SCT"
)

# A new "niche" assay and "niche_cluster" metadata column are added to each object
SpatialDimPlot(visium_list[["S1"]], group.by = "niche_cluster")
```

---

### 16.4 Spatial Polygon Tools

These functions support region-of-interest (ROI) analysis in Xenium and other spatial modalities.

```r
# Pull tissue coordinates across all images in an object
coords <- get_all_coords(integrated_xenium)

# Parse polygon definitions (e.g., from a JSON or CSV exported from Xenium Explorer)
polys <- parse_polygons("my_rois.json")

# Identify which cells fall inside a polygon
cells_in_roi <- get_cells_in_polygon(
  coords  = coords,
  polygon = polys[["ROI_1"]]
)

# Subset the object to those cells
roi_obj <- subset_opt(integrated_xenium, cells = cells_in_roi)
```

`subset_opt()` is a drop-in replacement for `subset()` that keeps spatial FOVs and images synchronized with the retained cells — preventing the "stale image" errors that `subset()` can produce on Xenium and Visium objects.

```r
# Always prefer subset_opt() over subset() for spatial objects
clean <- subset_opt(integrated_xenium, seurat_clusters %in% c("0", "1", "3"))
```

---

### 16.5 GO Term Children

`get_all_children()` recursively collects all descendant GO terms for a given parent term. Useful for building gene sets from broad ontology categories.

```r
# Get all children of "response to oxidative stress" (GO:0006979)
children <- get_all_children("GO:0006979")
length(children)
#> [1] 47

# Use with a GO annotation database to collect relevant genes
library(org.Mm.eg.db)
genes_in_terms <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = children,
  columns = c("SYMBOL"),
  keytype = "GO"
)
```

---

## 17. Tips and Common Pitfalls

**Read the `?docs`.** Every function has documented parameters beyond what is shown here. Defaults are sensible but the knobs are worth knowing about — especially `to_regress`, `doublet_normalization`, and `cluster_resolution`.

**Use `subset_opt()` instead of `subset()` for spatial data.** Seurat's built-in `subset()` can leave stale image metadata attached to spatial objects after subsetting. `subset_opt()` keeps FOVs and images consistent.

**Run `EdgeDetectionVisium()` before merging.** Edge and tear spots are almost always the lowest-quality cells in a Visium experiment. Removing them before integration prevents them from distorting the low-dimensional embedding.

**Check gene IDs before merging samples from different pipelines.** If one pipeline produced Ensembl IDs and another produced gene symbols, the merge will silently keep only genes whose names match literally — which is zero. Use `check_gene_ids_across_objects()` to catch this before it wastes an integration run.

**`common_genes_only = TRUE` in `MergeSeurat()`.** When integrating across different panel-based platforms or across datasets that used different genome annotations, set `common_genes_only = TRUE`. The function reports how many genes are retained and drops cells that end up with zero counts on the shared gene set.

**Interactive variants for new data types.** `CreateRNAObjectsFilter()` and `CreateATACObjectsFilter()` are worth using the first time you analyze a new tissue or library prep. The interactive prompts force you to look at the distribution before committing to cutoffs.

**For co-expression analysis**, this package deliberately doesn't cover full functionality. See [katlande/scCoExpress](https://github.com/katlande/scCoExpress). `NicheCoExpress()` handles the specific case of niche-differential Manders overlap.

**Use `PseudobulkDE()` — not `FindMarkers()` — for condition contrasts.** Per-cell tests dramatically inflate significance when cells from the same donor aren't treated as replicates. `PseudobulkDE()` aggregates per (sample, cluster) and runs DESeq2 on the pseudobulk matrix, which is the current best practice.

**Use `RunRCTD()` — not `AnnotateClusters()` — for Visium.** Visium spots are cell-type mixtures by construction; a winner-takes-all classifier hides minority signal. Deconvolution is the right primary output.

**In tumor samples, pass `tumor_mode = TRUE` to `AnnotateClusters()`.** The non-specific-gene filter is actively counterproductive when tumor cells drive up per-gene maxima; `tumor_mode` disables it and tightens `high_relative_to_max`. Better yet: identify malignant cells first (InferCNV / CopyKAT / SCEVAN), annotate the immune/stromal compartment normally, and score tumor programs separately.

**Save with provenance.** `SaveWithProvenance()` writes a JSON sidecar next to your `.rds` recording package versions, layers, reductions, metadata columns, and (optionally) the git SHA. Small overhead, invaluable when the dataset outlives the analysis notebook.

**Run `QCComparePlots()` after `ApplyQCFilters()`.** Confirms at a glance that the filter did what the report suggested — sometimes an outlier sample dominates the MAD estimate and the cutoffs are too tight.

---

## 18. Session Info

```r
sessionInfo()
```

Expected key packages and versions:

| Package | Minimum version |
|---|---|
| R | ≥ 4.1 |
| Seurat | ≥ 4.9.9 |
| SeuratObject | ≥ 4.9.9 |
| Signac | ≥ 1.9.0 |
| DoubletFinder | latest GitHub |
| EnsDb.Mmusculus.v79 | ≥ 2.99.0 |
| UCell | ≥ 2.4.0 |
| DESeq2 | ≥ 1.42 |
| speckle | ≥ 1.2 |
| slingshot | ≥ 2.10 |
| Azimuth | ≥ 0.5 (GitHub) |
| spacexr | ≥ 2.2 (GitHub) |
| liana | ≥ 0.1 (GitHub) |
| ks | ≥ 1.14 |
| dplyr | ≥ 1.1.2 |
| ggplot2 | ≥ 3.4.4 |
| patchwork | ≥ 1.1.2 |

---

*Vignette generated for SingleCellTools v2.4. For bug reports and feature requests, open an issue at [github.com/gevensen95/SingleCellTools/issues](https://github.com/gevensen95/SingleCellTools/issues).*
