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
   - 3.5 [Cell Annotation — `MarkerPlot()`](#35-cell-annotation--markerplot)
   - 3.6 [Gene Positivity — `AddGenePositivity()`](#36-gene-positivity--addgenepositivity)
4. [One-Shot Pipeline — `CreateAndIntegrateRNA()`](#4-one-shot-pipeline--createandintegrаterna)
5. [Spatial Workflows](#5-spatial-workflows)
   - 5.1 [Visium — `CreateVisiumObjects()` + `EdgeDetectionVisium()`](#51-visium--createvisiumobjects--edgedetectionvisium)
   - 5.2 [Xenium — `LoadXenium2()`](#52-xenium--loadxenium2)
   - 5.3 [Parse Biosciences — `MakeParseObj()`](#53-parse-biosciences--makeparseobj)
6. [scATAC-seq — `CreateATACObjects()`](#6-scatac-seq--createatacobjects)
7. [Utility Functions](#7-utility-functions)
   - 7.1 [Gene ID Diagnostics](#71-gene-id-diagnostics)
   - 7.2 [Cell-Cycle Scoring](#72-cell-cycle-scoring)
   - 7.3 [Spatial Niche Analysis](#73-spatial-niche-analysis)
   - 7.4 [Spatial Polygon Tools](#74-spatial-polygon-tools)
   - 7.5 [GO Term Children](#75-go-term-children)
8. [Tips and Common Pitfalls](#8-tips-and-common-pitfalls)
9. [Session Info](#9-session-info)

---

## 1. Overview

`SingleCellTools` is an opinionated set of R wrapper functions that reduce the boilerplate of common single-cell analysis tasks. Rather than replacing Seurat or Signac, it wraps them — consolidating the five-to-ten function calls you write for every project into single, well-parametrized entry points.

The package covers:

- Reading **10x Genomics**, **Parse Biosciences**, **Visium**, **Xenium**, and **scATAC-seq** outputs into Seurat objects
- **Doublet calling** with DoubletFinder, with automatic pK optimization
- **Merging and integration** (Harmony, RPCA, CCA, or JointPCA) with normalization, PCA, clustering, and UMAP in one call
- **Edge-spot detection** for Visium data
- Annotated **dot plots** with automatic gene filtering and identity clustering
- Metadata helpers including gene positivity flags and cell-cycle phase scoring

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
  "EnsDb.Mmusculus.v79",  # mouse gene annotations
  "glmGamPoi",            # fast SCTransform fitting
  "GO.db",                # GO term lookups
  "UCell",                # module scoring
  "Signac"                # scATAC-seq
))
```

### Optional dependencies

```r
install.packages("sf")   # required only for get_cells_in_polygon()
```

### Load the package

```r
library(SingleCellTools)
```

---

## 3. Core scRNA-seq Workflow

The typical pipeline looks like this:

```
CellRanger outputs
      ↓
CreateRNAObjects()        ← load, compute %mt, call doublets
      ↓
[manual subset/filter]    ← subset() based on QC plots
      ↓
MergeSeurat()             ← normalize, integrate, cluster, UMAP
      ↓
MarkerPlot()              ← annotate cell types
      ↓
AddGenePositivity()       ← flag gene-positive cells
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

### 3.5 Cell Annotation — `MarkerPlot()`

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

### 3.6 Gene Positivity — `AddGenePositivity()`

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

## 5. Spatial Workflows

### 5.1 Visium — `CreateVisiumObjects()` + `EdgeDetectionVisium()`

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

### 5.2 Xenium — `LoadXenium2()`

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

### 5.3 Parse Biosciences — `MakeParseObj()`

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

## 6. scATAC-seq — `CreateATACObjects()`

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

## 7. Utility Functions

### 7.1 Gene ID Diagnostics

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

### 7.2 Cell-Cycle Scoring

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

### 7.3 Spatial Niche Analysis

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

### 7.4 Spatial Polygon Tools

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

### 7.5 GO Term Children

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

## 8. Tips and Common Pitfalls

**Read the `?docs`.** Every function has documented parameters beyond what is shown here. Defaults are sensible but the knobs are worth knowing about — especially `to_regress`, `doublet_normalization`, and `cluster_resolution`.

**Use `subset_opt()` instead of `subset()` for spatial data.** Seurat's built-in `subset()` can leave stale image metadata attached to spatial objects after subsetting. `subset_opt()` keeps FOVs and images consistent.

**Run `EdgeDetectionVisium()` before merging.** Edge and tear spots are almost always the lowest-quality cells in a Visium experiment. Removing them before integration prevents them from distorting the low-dimensional embedding.

**Check gene IDs before merging samples from different pipelines.** If one pipeline produced Ensembl IDs and another produced gene symbols, the merge will silently keep only genes whose names match literally — which is zero. Use `check_gene_ids_across_objects()` to catch this before it wastes an integration run.

**`common_genes_only = TRUE` in `MergeSeurat()`.** When integrating across different panel-based platforms or across datasets that used different genome annotations, set `common_genes_only = TRUE`. The function reports how many genes are retained and drops cells that end up with zero counts on the shared gene set.

**Interactive variants for new data types.** `CreateRNAObjectsFilter()` and `CreateATACObjectsFilter()` are worth using the first time you analyze a new tissue or library prep. The interactive prompts force you to look at the distribution before committing to cutoffs.

**For co-expression analysis**, this package deliberately doesn't cover it. See [katlande/scCoExpress](https://github.com/katlande/scCoExpress).

---

## 9. Session Info

```r
sessionInfo()
```

Expected key packages and versions:

| Package | Minimum version |
|---|---|
| R | ≥ 2.10 |
| Seurat | ≥ 4.9.9 |
| SeuratObject | ≥ 4.9.9 |
| Signac | ≥ 1.9.0 |
| DoubletFinder | latest GitHub |
| EnsDb.Mmusculus.v79 | ≥ 2.99.0 |
| UCell | ≥ 2.4.0 |
| dplyr | ≥ 1.1.2 |
| ggplot2 | ≥ 3.4.4 |

---

*Vignette generated for SingleCellTools v2.4. For bug reports and feature requests, open an issue at [github.com/gevensen95/SingleCellTools/issues](https://github.com/gevensen95/SingleCellTools/issues).*
