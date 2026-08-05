# SingleCellTools Vignette: Spatial Transcriptomics with Mouse Brain Visium Data

**Package:** `SingleCellTools`  
**Data:** `stxBrain` — 10x Genomics Visium mouse brain sections, available via `SeuratData`. We use all four serial sections it ships (`anterior1`, `anterior2`, `posterior1`, `posterior2`) so that anterior-vs-posterior brain region has two sections each, which is enough to run a real `NicheCoExpress()` comparison later in the vignette.  
**Goal:** Walk through a complete Visium workflow — edge detection, integration, annotation, niche analysis — using freely available public data.

---

## Table of Contents

1. [Setup](#1-setup)
2. [Load the Data](#2-load-the-data)
3. [QC and Percent Mitochondrial Reads](#3-qc-and-percent-mitochondrial-reads)
4. [Edge & Hole Detection — `detect_fov_edges()` / `detect_tissue_holes()`](#4-edge--hole-detection--detect_fov_edges--detect_tissue_holes)
5. [Merge and Integrate — `MergeSeurat()`](#5-merge-and-integrate--mergeseurat)
6. [Annotate Spatial Domains](#6-annotate-spatial-domains)
   - 6.1 [Marker Dot Plots — `MarkerPlot()` / `MarkerPctPlot()`](#61-marker-dot-plots--markerplot--markerpctplot)
   - 6.2 [Cluster-level Marker Scoring — `AnnotateClusters()`](#62-cluster-level-marker-scoring--annotateclusters)
7. [Visium Deconvolution — `RunRCTD()`](#7-visium-deconvolution--runrctd)
8. [Feature Density on the Tissue — `PlotFeatureDensity()`](#8-feature-density-on-the-tissue--plotfeaturedensity)
9. [Gene Positivity — `AddGenePositivity()` / `PlotGenePositivity()`](#9-gene-positivity--addgenepositivity--plotgenepositivity)
    - 9.1 [`GenePositivityAnalysis()` and `GenePositivityEstimationPlot()`](#91-genepositivityanalysis-and-genepositivityestimationplot)
10. [Spatial Niche Analysis — `BuildMultipleNicheAssays()`](#10-spatial-niche-analysis--buildmultiplenicheassays)
    - 10.1 [An Alternative: Spatial Domain Segmentation via `RunBanksyWrapper()`](#101-an-alternative-spatial-domain-segmentation-via-runbanksywrapper)
11. [Neighborhood Enrichment — `NeighborhoodEnrichment()`](#11-neighborhood-enrichment--neighborhoodenrichment)
12. [Niche Co-expression — `NicheCoExpress()`](#12-niche-co-expression--nicheco express)
    - 12.1 [Estimation Plot — `NicheCoExpressEstimationPlot()`](#121-estimation-plot--nichecoexpressestimationplot)
13. [Single-cell Spatial (Xenium / CosMx) — `detect_fov_edges()` / `detect_tissue_holes()`](#13-single-cell-spatial-xenium--cosmx--detect_fov_edges--detect_tissue_holes)
14. [Subsetting Spatial Objects — `subset_opt()`](#14-subsetting-spatial-objects--subset_opt)
15. [Tips Specific to Spatial Data](#15-tips-specific-to-spatial-data)
16. [Session Info](#16-session-info)

---

## 1. Setup

### Install required packages

```r
# SingleCellTools
remotes::install_github("gevensen95/SingleCellTools")

# Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("glmGamPoi", "UCell"))

# SeuratData — provides the stxBrain dataset
remotes::install_github("satijalab/seurat-data")

# ClusterR is needed for BuildMultipleNicheAssays
install.packages("ClusterR")
```

### Install and cache the dataset

```r
library(SeuratData)
InstallData("stxBrain")   # downloads ~200 MB once; cached for future sessions
```

### Load libraries

```r
library(SingleCellTools)
library(Seurat)
library(SeuratData)
library(ggplot2)
library(dplyr)
library(patchwork)
```

---

## 2. Load the Data

`stxBrain` ships four serial coronal mouse brain Visium sections: two adjacent
`anterior` sections and two adjacent `posterior` sections. Earlier drafts of this
vignette only loaded one section per region (`anterior1`/`posterior1`), which is
plenty for integration but leaves no way to demonstrate a real `NicheCoExpress()`
comparison later, since that function needs at least 2 samples *per condition*.
Loading all four sections and treating anterior-vs-posterior as the comparison of
interest gives us exactly that: 2 samples per side, using real public data rather
than anything synthetic.

```r
# Load all four sections
brain_ant1 <- LoadData("stxBrain", type = "anterior1")
brain_ant2 <- LoadData("stxBrain", type = "anterior2")
brain_pos1 <- LoadData("stxBrain", type = "posterior1")
brain_pos2 <- LoadData("stxBrain", type = "posterior2")

brain_list <- list(
  anterior1  = brain_ant1,
  anterior2  = brain_ant2,
  posterior1 = brain_pos1,
  posterior2 = brain_pos2
)

# Quick overview
lapply(brain_list, function(x) c(features = nrow(x), spots = ncol(x)))
#> $anterior1
#> features    spots
#>    31053     2696
#> $anterior2
#> features    spots
#>    31053     2825
#> $posterior1
#> features    spots
#>    31053     3353
#> $posterior2
#> features    spots
#>    31053     3289
```

> **A note on Visium resolution.** Each spot in a Visium capture area is ~55 µm in
> diameter — in brain tissue this typically captures the expression of 2–10 cells
> simultaneously. Clustering on Visium data therefore identifies **spatial domains**
> (brain regions with a coherent transcriptional signature) rather than individual cell
> types. Deconvolution methods (e.g., RCTD, SPOTlight) can estimate cell-type
> composition per spot, but are outside the scope of this vignette.

> **A note on what "region" means here.** Treating anterior-vs-posterior as a
> `NicheCoExpress()` condition later in this vignette is a real, honest use of public
> data — but it's an anatomical/spatial grouping, not an experimental treatment. It
> answers "does this gene pair co-express differently between anterior and posterior
> cortex," which is a legitimate biological question, but don't copy this pattern and
> assume it validates for treatment-vs-control designs. For that, swap in your own
> `sample_col`/`condition_col` from a dataset with real biological replicates per
> treatment arm.

Visualize the raw tissue sections to confirm the data loaded correctly (shown for
`anterior1`; the same call works for any of the four):

```r
SpatialFeaturePlot(brain_list$anterior1, features = "nCount_Spatial") +
  ggtitle("Anterior 1 — total UMI per spot")
```

---

## 3. QC and Percent Mitochondrial Reads

With four sections it's cleaner to loop rather than repeat each call four times:

```r
# Mouse mitochondrial genes: "^mt-" (lowercase)
brain_list <- lapply(brain_list, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj
})

# QC distributions -- shown for anterior1; inspect each section the same way
VlnPlot(brain_list$anterior1,
        features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
        pt.size  = 0.1) &
  theme(axis.text.x = element_blank())
```

Spatial scatter to see if low-quality spots are spatially clustered (often they are,
at tissue edges or tears):

```r
SpatialFeaturePlot(brain_list$anterior1, features = "percent.mt") +
  scale_fill_gradientn(colors = c("grey90", "red3")) +
  ggtitle("Percent mitochondrial reads — anterior 1")
```

Apply filters to every section. Visium thresholds are typically looser than
single-cell because spots have higher total counts:

```r
brain_list <- lapply(brain_list, function(obj) {
  subset(obj,
    nFeature_Spatial > 200  &
    nFeature_Spatial < 8000 &
    percent.mt < 25
  )
})

sapply(brain_list, ncol)
#>  anterior1  anterior2 posterior1 posterior2
#>       2611       2739       3244       3172
```

---

## 4. Edge & Hole Detection — `detect_fov_edges()` / `detect_tissue_holes()`

Spots at the capture-area boundary have systematically lower UMI counts and higher
noise; spots bordering internal tissue tears or holes have the same problem.
`detect_fov_edges()` flags the former (outer boundary), `detect_tissue_holes()`
flags the latter (internal gaps) — both work on any spatial modality since they
pull coordinates via `Seurat::GetTissueCoordinates()` rather than assuming a
hex-grid layout.

> **`detect_fov_edges()` vs. `EdgeDetectionVisium()`.** SingleCellTools also ships
> `EdgeDetectionVisium()`, which flags a spot as an edge if it has too few
> neighbors within a radius — a good fit for Visium's fixed hex-grid spacing
> (every interior spot has exactly 6 neighbors). `detect_fov_edges()`'s default
> `method = "bbox"`, used below, instead flags spots close to the axis-aligned
> bounding box of the tissue — fast and simple, but it assumes the tissue is
> roughly convex within each section. A coronal brain section is convex enough
> for this to work well; a tissue with a very irregular or crescent-shaped
> outline may need `method = "angular"` (also available on `detect_fov_edges()`)
> or `EdgeDetectionVisium()`'s neighbor-count approach instead. Neither method is
> a strict upgrade over the other — pick based on tissue shape.

```r
brain_list <- lapply(brain_list, function(obj) {
  obj <- detect_fov_edges(obj,
                          method       = "bbox",
                          bbox_factor  = 2,
                          n_iterations = 2,
                          label_col    = "edge_layer")
  obj <- detect_tissue_holes(obj,
                             min_hole_size = 4,
                             n_iterations  = 2,
                             label_col     = "hole_layer")
  obj
})

# 0 = interior, 1 = outermost ring/layer, 2 = next ring/layer, ... (anterior1 shown)
table(brain_list$anterior1$edge_layer)
table(brain_list$anterior1$hole_layer)
```

Visualize both layers before committing to a cutoff:

```r
SpatialDimPlot(brain_list$anterior1, group.by = "edge_layer") +
  ggtitle("Edge layers — anterior 1")

SpatialDimPlot(brain_list$anterior1, group.by = "hole_layer") +
  ggtitle("Hole layers — anterior 1")
```

> **Marker-gene exclusion for biologically meaningful gaps.** If a "hole" is real
> anatomy rather than a tear (a ventricle, a large vessel, a central vein), pass
> `exclude_gene` to `detect_tissue_holes()` so bins bordered by high expression of
> that gene are skipped rather than flagged — see Section 13 for a worked example
> with `Glul` in liver tissue. Not needed for this mouse brain dataset, but worth
> knowing about before you discard real anatomical structures in your own tissue.

If both layers look right, keep only interior spots (`layer == 0` in both —
mirroring exactly how Section 13 combines the same two columns for Xenium data):

```r
brain_list <- lapply(brain_list, function(obj) {
  subset(obj, edge_layer == 0 & hole_layer == 0)
})

sapply(brain_list, ncol)
#>  anterior1  anterior2 posterior1 posterior2
#>       2558       2685       3181       3115
```

> **Choosing `n_iterations`, `bbox_factor`, and `min_hole_size`.** `n_iterations = 2`
> peels back two rings/layers on each function; raise it for noisier tissue edges.
> `bbox_factor` (edge detection) and `min_hole_size` (hole detection, in grid bins)
> both trade off aggressiveness vs. over-trimming — validate by comparing UMI
> distributions of kept vs. flagged spots before committing, the same way you would
> with `EdgeDetectionVisium()`'s `Filter`/`Filter4` columns.

---

## 5. Merge and Integrate — `MergeSeurat()`

With four clean, filtered sections, we are ready to merge and integrate. Pass
`spatial = "Visium"` so `MergeSeurat()` handles the Spatial assay correctly and
keeps images in sync. `MergeSeurat()` takes the named list directly.

```r
integrated <- MergeSeurat(
  seurat_objects = brain_list,

  # Normalization
  use_SCT    = TRUE,
  to_regress = "percent.mt",

  # Dimensionality
  max_dims = 30,

  # Integration
  spatial        = "Visium",
  integration    = "HarmonyIntegration",
  new_reduction  = "harmony",

  # Clustering — slightly higher resolution for spatial data
  cluster_resolution = 0.5,

  # Output
  save_rds_file = TRUE,
  file_name     = "brain_visium",
  markers       = TRUE
)

# The list names become orig.ident; derive the anterior/posterior region label
# from that. This is the "condition" NicheCoExpress() will compare in Section 12 --
# see the note in Section 2 on what that comparison does and doesn't demonstrate.
integrated$region <- ifelse(grepl("^anterior", integrated$orig.ident),
                            "anterior", "posterior")
table(integrated$orig.ident, integrated$region)
```

After the run your working directory contains:
- `brain_visium_merged_seurat_objects.rds`
- `dimplot_seurat_clusters.pdf` — UMAP colored by cluster
- `markers_all.csv` and `marker_plot.pdf`

### Inspect the integration

```r
# UMAP by section — good integration means all four sections overlap
DimPlot(integrated, group.by = "orig.ident",
        cols = c("anterior1"  = "steelblue", "anterior2"  = "dodgerblue",
                 "posterior1" = "tomato",     "posterior2" = "firebrick")) +
  ggtitle("By section — should be interleaved")

# UMAP by region — the coarser anterior/posterior grouping used later
DimPlot(integrated, group.by = "region",
        cols = c("anterior" = "steelblue", "posterior" = "tomato")) +
  ggtitle("By region")

DimPlot(integrated, label = TRUE, repel = TRUE) +
  ggtitle("Seurat clusters")

# Spatial view of clusters on each section
SpatialDimPlot(integrated, label = FALSE) +
  ggtitle("Clusters mapped back to tissue")
```

---

## 6. Annotate Spatial Domains

### 6.1 Marker Dot Plots — `MarkerPlot()` / `MarkerPctPlot()`

We now assign neuroanatomical labels to clusters using canonical mouse brain cell-type
and region markers. `MarkerPlot()` groups genes by annotation, clusters identities by
expression correlation, and automatically drops genes with zero or uniform expression.

```r
brain_markers <- data.frame(
  Gene = c(
    # Neurons (pan)
    "Snap25", "Syt1",   "Rbfox3",
    # Excitatory neurons
    "Slc17a7", "Camk2a", "Grin2a",
    # Inhibitory neurons
    "Gad1",   "Gad2",   "Slc32a1",
    # Astrocytes
    "Gfap",   "Aqp4",   "Aldh1l1",
    # Oligodendrocytes
    "Mbp",    "Mog",    "Plp1",
    # OPCs
    "Pdgfra", "Cspg4",
    # Microglia
    "Cx3cr1", "P2ry12", "Tmem119",
    # Endothelial
    "Cldn5",  "Pecam1",
    # Choroid plexus
    "Ttr",    "Folr1"
  ),
  CellType = c(
    rep("Neuron",           3),
    rep("Excitatory",       3),
    rep("Inhibitory",       3),
    rep("Astrocyte",        3),
    rep("Oligodendrocyte",  3),
    rep("OPC",              2),
    rep("Microglia",        3),
    rep("Endothelial",      2),
    rep("Choroid plexus",   2)
  )
)

DefaultAssay(integrated) <- "SCT"

p <- MarkerPlot(
  obj              = integrated,
  genes            = brain_markers,
  assay            = "SCT",
  cluster          = TRUE,
  show.annotations = TRUE,
  maxsize          = 5,
  label.fontsize   = 3,
  margin_factor    = 0.6   # increase if right-edge labels are clipped
)
print(p)
ggsave("markerplot_brain.pdf", p, width = 14, height = 10)
```

### Assign neuroanatomical labels

After reviewing the dot plot, relabel clusters with their anatomical identity.
Your cluster–label mapping will differ from the example below depending on resolution
and the specific integration run — use your marker plot as the guide.

```r
# Example mapping — adjust to match what you observe
region_labels <- c(
  "0"  = "Oligodendrocyte",
  "1"  = "Excitatory neuron",
  "2"  = "Excitatory neuron",
  "3"  = "Inhibitory neuron",
  "4"  = "Astrocyte",
  "5"  = "Excitatory neuron",
  "6"  = "Microglia",
  "7"  = "OPC",
  "8"  = "Endothelial",
  "9"  = "Choroid plexus",
  "10" = "Inhibitory neuron"
)

integrated <- RenameIdents(integrated, region_labels)
integrated$spatial_domain <- Idents(integrated)

# Visualize labeled domains on tissue
SpatialDimPlot(integrated, group.by = "spatial_domain", label = FALSE) +
  ggtitle("Annotated spatial domains")
ggsave("spatial_domains_brain.pdf", width = 14, height = 7)
```

---

### 6.2 Cluster-level Marker Scoring — `AnnotateClusters()`

For Visium, use `AnnotateClusters()` with Visium-appropriate defaults — the winner-takes-all limitation applies here (each spot mixes cell types), so `return_scores = "cluster"` is often more informative than the winner label alone.

```r
brain_markers <- list(
  Neuron_Ex   = c("Slc17a7", "Camk2a"),      # excitatory
  Neuron_In   = c("Gad1", "Gad2"),           # inhibitory
  Astrocyte   = c("Gfap", "Aqp4"),
  Oligo       = c("Mbp", "Mog", "Plp1"),
  Endothelial = c("Pecam1", "Cldn5"),
  Microglia   = c("Cx3cr1", "P2ry12")
)

# Visium-safe defaults: no non-specific filter, low min_detection_frac
res <- AnnotateClusters(
  integrated,
  markers             = brain_markers,
  cluster_col         = "seurat_clusters",
  filter_nonspecific  = FALSE,
  min_detection_frac  = 0.05,
  return_scores       = "cluster",
  new_col             = "domain"
)
res$scores       # cluster x cell-type UCell score matrix
integrated <- res$obj

# Winner label is available too
table(integrated$domain)
```

**Deconvolution is preferred for Visium** — see section 7 below. Use marker scoring only as a sanity check.

---

## 7. Visium Deconvolution — `RunRCTD()`

Each Visium spot covers ~55 μm and contains 1-10 cells of mixed types. Winner-takes-all classifiers (`MarkerPlot`, `AnnotateClusters`) lose minority populations by construction. `RunRCTD()` deconvolves each spot as a mixture using a reference single-cell dataset.

```r
# You need a single-cell reference Seurat object with cell-type labels.
# For mouse brain, the Allen Brain Institute reference or Zeisel et al. work.
# Placeholder here:
brain_ref <- readRDS("data/allen_mouse_cortex.rds")

integrated <- RunRCTD(
  integrated,
  reference    = brain_ref,
  celltype_col = "cell_type",
  mode         = "full",              # or "doublet" for cleaner but binary calls
  max_cells_per_ref_celltype = 5000,
  n_cores      = 8
)

# Per-cell-type proportion columns on each spot
head(colnames(integrated@misc$rctd_weights))
# [1] "L2/3 IT" "L5 IT"   "L6 CT"   "Oligo"   ...

# Overlay proportions spatially
SpatialFeaturePlot(integrated,
                   features = c("rctd_Oligo", "rctd_Astro", "rctd_L2.3.IT"))

# Dominant type per spot for quick visualization
SpatialDimPlot(integrated, group.by = "rctd_dominant")
```

For laminar cortex or clearly zoned tissue, RCTD proportions map cleanly onto the anatomical structure — a strong sanity check that everything is working.

---

## 8. Feature Density on the Tissue — `PlotFeatureDensity()`

`PlotFeatureDensity()` gives a much cleaner view of sparse markers than `SpatialFeaturePlot()` when combined with the UMAP or a spatial 2D coordinate reduction. On a UMAP:

`PlotFeatureDensity()` isn't limited to genes — any numeric metadata column works,
including a module score. `AddModuleScore()` writes a column named `<name>1`; copy
it to a fixed, readable name so the call below doesn't depend on that numbering:

```r
integrated <- AddModuleScore(
  integrated,
  features = list(c("Mbp", "Mog", "Plp1")),   # oligodendrocyte / myelination markers
  name     = "OligoScore"
)
integrated$module_score_oligo <- integrated$OligoScore1

PlotFeatureDensity(
  integrated,
  features  = c("Mbp", "Gad1", "Slc17a7"),
  reduction = "umap"
)

# On the spatial coordinates directly (works when 'spatial' is registered
# as a reduction; otherwise use SpatialFeaturePlot for the tissue view).
PlotFeatureDensity(integrated,
                   features  = "module_score_oligo",
                   reduction = "umap")
```

For side-by-side comparison of a gene's density and its ligand's density:

```r
PlotFeatureDensity(integrated,
                   features = c("Vegfa", "Kdr"),
                   joint    = TRUE)
```

---

## 9. Gene Positivity — `AddGenePositivity()` / `PlotGenePositivity()`

Flag spots where a gene is detectably expressed (counts > 0). Useful for quick spatial
co-expression checks and for creating binary masks before niche analysis.

```r
integrated <- AddGenePositivity(
  seurat_objects = integrated,
  genes          = c("Snap25", "Gfap", "Mbp", "Cx3cr1", "Cldn5"),
  layer          = "counts",
  threshold      = 0,
  suffix         = "_pos"
)

# Fraction of spots positive per spatial domain
integrated@meta.data %>%
  group_by(spatial_domain) %>%
  summarise(
    pct_Snap25  = round(mean(Snap25_pos)  * 100, 1),
    pct_Gfap    = round(mean(Gfap_pos)    * 100, 1),
    pct_Mbp     = round(mean(Mbp_pos)     * 100, 1),
    pct_Cx3cr1  = round(mean(Cx3cr1_pos)  * 100, 1)
  ) %>%
  arrange(desc(pct_Snap25))
```

Spatial map of a positivity flag:

```r
SpatialFeaturePlot(integrated, features = "Gfap_pos") +
  scale_fill_manual(values = c("FALSE" = "grey90", "TRUE" = "firebrick")) +
  ggtitle("Gfap-positive spots")
```

Percent-positive summary across spatial domains via `PlotGenePositivity()`:

```r
PlotGenePositivity(integrated,
                   c("Snap25", "Gfap", "Mbp", "Cx3cr1"),
                   group.by = "spatial_domain",
                   style    = "heatmap",
                   max_pct  = 90)

# Or the co-expression view for two markers
PlotGenePositivity(integrated,
                   c("Snap25", "Gfap"),
                   group.by = "spatial_domain",
                   style    = "combo")
```

### 9.1 `GenePositivityAnalysis()` and `GenePositivityEstimationPlot()`

The summary above is pooled across all spots regardless of section. `GenePositivityAnalysis()`
computes positivity rates per **sample** (section) instead — stratified by `spatial_domain`
here — plus an optional chi-square/Fisher test across `region`, the same
`orig.ident`/`region` anterior-vs-posterior design `NicheCoExpress()` (Section 12) uses.
With only 4 sections total (2 per region) that test pools spots across sections rather
than testing on the 2-vs-2 sample-level rates, so it emits a `warning()` and shouldn't
be read as a real significance result here -- `GenePositivityEstimationPlot()` below,
which works directly off the per-sample rates, is the honest view of this comparison:

```r
gpa <- GenePositivityAnalysis(
  integrated,
  genes         = c("Gfap", "Mbp"),
  sample_col    = "orig.ident",
  condition_col = "region",
  group_col     = "spatial_domain",
  test          = "chisq"
)
gpa$proportions                    # sample, gene, group, n_pos, n_total, prop_pos, condition
gpa$test[["Mbp | Oligodendrocyte"]]   # chi-square result for that gene x domain combination (pooled-spot caveat above)
```

As with `CompositionAnalysis()`/`CompositionEstimationPlot()` in the other two
vignettes, the test above answers "is there a difference"; `GenePositivityEstimationPlot()`
answers "by how much, with what uncertainty" using per-sample positivity rates and a
bootstrap 95% CI -- the same caveat from Section 2 applies here too: "anterior" vs
"posterior" is an anatomical grouping (2 sections each), not a treatment arm.

```r
# Mbp positivity in the oligodendrocyte domain, anterior vs posterior
GenePositivityEstimationPlot(gpa, genes = "Mbp", group_levels = "Oligodendrocyte",
                             idx = c("anterior", "posterior"))

# Every gene x domain combination present, as a named list of plots
plots <- GenePositivityEstimationPlot(gpa, idx = c("anterior", "posterior"))
plots[["Mbp | Oligodendrocyte"]]

# Cohen's h -- the effect size designed specifically for comparing two proportions
GenePositivityEstimationPlot(gpa, genes = "Mbp", group_levels = "Oligodendrocyte",
                             idx = c("anterior", "posterior"), effect = "cohens_h")
```

#### How to read a `dabestr` estimation plot

This is the first estimation plot in this vignette, so it's worth spelling out what's
actually on it -- the same layout applies to `NicheCoExpressEstimationPlot()` (12.1)
later on too ([`dabestr`](https://acclab.github.io/dabestr/) implements "estimation
statistics": Ho et al. 2019, *Moving beyond P values: data analysis with estimation
graphics*, Nature Methods).

Each plot has two stacked panels sharing an x-axis (the two `idx` conditions, here
`anterior`/`posterior`):

- **Top panel — raw data.** Every individual value feeding the comparison is plotted
  as a swarm along its condition's column. That's one point per *section*
  (`orig.ident`), not per spot -- with only 4 sections total (2 per region) you'll
  literally see all 4 points, which is the honest picture this comparison can
  support, rather than a p-value implying more precision than 2-vs-2 actually gives.
- **Bottom panel — effect size.** The chosen `effect` (`mean_diff` by default) between
  posterior and anterior, drawn as a single point with a vertical bar for its 95%
  bootstrap confidence interval (5000 resamples, by default). A dashed horizontal
  line marks the reference condition's (anterior's) value, so the effect-size
  point/CI can be read directly against it.

Because the CI comes from resampling the actual sections rather than a parametric
formula, it doesn't assume normality -- but with only 2 sections per condition here,
resampling has very little to work with, so expect a wide CI; that width is itself
useful information (this comparison can detect large effects only, not subtle ones).
There's no p-value threshold to eyeball; instead look at whether the CI includes zero
(no effect) and how wide it is.

`effect` options, all available on every `*EstimationPlot()` function in this package:

| `effect` | What it measures |
|---|---|
| `mean_diff` (default) | Test mean − reference mean, in the original units (e.g. proportion points). |
| `median_diff` | Same, using medians -- less sensitive to one outlier section. |
| `cohens_d` | Standardized mean difference (pooled-SD units); comparable across differently-scaled measurements. |
| `hedges_g` | `cohens_d` with a small-sample bias correction -- prefer this over `cohens_d` with only 2 sections per condition. |
| `cliffs_delta` | Non-parametric, rank-based (akin to a standardized Mann-Whitney effect) -- robust to outliers and non-normal spread, no distributional assumptions. |
| `cohens_h` | Designed specifically for comparing two proportions -- often the more principled choice for `prop_pos`/`prop` columns (bounded 0-1) rather than a raw mean difference. |

---

## 10. Spatial Niche Analysis — `BuildMultipleNicheAssays()`

A spatial niche captures not just what a spot expresses, but what its *neighbors*
express. `BuildMultipleNicheAssays()` builds a neighborhood-composition assay across
all sections and then clusters spots by niche profile using mini-batch k-means. This
reveals recurring spatial patterns that transcend individual sections.

The function requires a `group.by` metadata column that classifies spots — here we
use the `spatial_domain` labels assigned in Section 6. It also needs the list of FOV
(image) names that correspond to each object, and `type = "visium"` to use the correct
coordinate accessor.

```r
# Split the integrated object back into per-section objects for niche building.
# BuildMultipleNicheAssays() works on a list, not the merged object.
brain_list_annotated <- SplitObject(integrated, split.by = "orig.ident")

# Confirm FOV names match what's in each object's images slot
lapply(brain_list_annotated, function(obj) names(obj@images))
#> $anterior1
#> [1] "anterior1"
#> $anterior2
#> [1] "anterior2"
#> $posterior1
#> [1] "posterior1"
#> $posterior2
#> [1] "posterior2"

niche_list <- BuildMultipleNicheAssays(
  list.object  = brain_list_annotated,
  list.fov     = as.list(names(brain_list_annotated)),  # must match names(obj@images)
  group.by     = "spatial_domain",
  assay        = "niche",
  cluster.name = "niches",
  neighbors.k  = 6,              # Visium spots have ~6 direct neighbors in the hexagonal grid
  niches.k.range = 4:15,         # test k = 4 through 15; inspect stability
  batch_size   = 20,
  num_init     = 20,
  type         = "visium"        # tells the function to use Visium coordinate accessors
)
```

`niche_list` now has four elements (`anterior1`, `anterior2`, `posterior1`,
`posterior2`). The examples below use `anterior1`/`posterior1` for illustration —
the same calls work for the other two sections.

`BuildMultipleNicheAssays()` adds:
- A `niche` assay to each object (features = cell-type labels, values = neighborhood composition)
- `niches.kmeans_4` through `niches.kmeans_15` metadata columns (one per k tested)

### Choosing k

Examine a few values of k spatially to find the one that produces interpretable,
anatomically consistent niche patterns:

```r
# Visualize different k values on the anterior section
p4  <- SpatialDimPlot(niche_list[["anterior1"]],
                      group.by = "niches.kmeans_4",  label = FALSE) + ggtitle("k = 4")
p8  <- SpatialDimPlot(niche_list[["anterior1"]],
                      group.by = "niches.kmeans_8",  label = FALSE) + ggtitle("k = 8")
p12 <- SpatialDimPlot(niche_list[["anterior1"]],
                      group.by = "niches.kmeans_12", label = FALSE) + ggtitle("k = 12")

p4 | p8 | p12
ggsave("niche_k_comparison_anterior.pdf", width = 18, height = 6)
```

Pick the k where niches map cleanly onto known anatomical boundaries (cortical layers,
white matter, hippocampus, etc.) without fragmented, salt-and-pepper patterns. For the
mouse brain, k = 6–10 often works well.

```r
# Once you've chosen k, add a clean "best_k" column for downstream use
best_k <- 8
for (i in seq_along(niche_list)) {
  niche_list[[i]]$best_niche <- niche_list[[i]][[paste0("niches.kmeans_", best_k)]]
}

# Compare niches across sections
SpatialDimPlot(niche_list[["anterior1"]],  group.by = "best_niche") +
  ggtitle("Anterior niches") |
SpatialDimPlot(niche_list[["posterior1"]], group.by = "best_niche") +
  ggtitle("Posterior niches")
ggsave("best_niches_both_sections.pdf", width = 14, height = 7)
```

### Niche composition heatmap

Which spatial domains make up each niche?

```r
library(tidyr)

# Use the anterior section as an example
niche_meta <- niche_list[["anterior1"]]@meta.data

comp <- niche_meta %>%
  count(best_niche, spatial_domain) %>%
  group_by(best_niche) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

ggplot(comp, aes(x = best_niche, y = spatial_domain, fill = pct)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("white", "steelblue", "navy"),
                       name   = "% of niche") +
  labs(x = "Niche", y = "Spatial domain",
       title = "Spatial domain composition per niche — anterior section") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("niche_composition_heatmap.pdf", width = 8, height = 6)
```

### 10.1 An Alternative: Spatial Domain Segmentation via `RunBanksyWrapper()`

`BuildMultipleNicheAssays()` above finds niches from neighborhood *composition* (what
cell-type labels surround each spot). `RunBanksyWrapper()` takes a different approach:
it wraps `SeuratWrappers::RunBanksy()`/`Banksy::computeBanksy()` to build a joint
expression + spatial-neighbor feature space directly from gene expression, then lets
you cluster on that instead. `lambda` controls the trade-off: low values (~0.2) favor
ordinary cell-typing, high values (~0.8) favor spatial domain segmentation — since
that's the goal here, use a high `lambda`.

This runs on `integrated` from Section 5 (the merged, SCTransform-normalized, harmony-
integrated object across all 4 sections) rather than the per-section `niche_list`.
Because BANKSY needs spatial coordinates and `integrated` has 4 sections' worth of
spots that would otherwise spatially overlap (each section's pixel coordinates start
near the origin independently), pass `group = "orig.ident"` so coordinates are
staggered per section before computing spatial neighbors -- skipping this would let
BANKSY treat spots from *different* sections as if they were physically adjacent.

```r
integrated <- RunBanksyWrapper(
  integrated,
  lambda   = 0.8,          # high lambda -- spatial domain segmentation, not cell-typing
  k_geom   = 15,
  group    = "orig.ident",  # keep the 4 sections' coordinates from overlapping
  assay_name = "BANKSY"
)

# RunBanksyWrapper() already ran PCA on the BANKSY assay (reduction "pca_banksy")
integrated <- FindNeighbors(integrated, reduction = "pca_banksy", dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5, cluster.name = "banksy_domain")

# Compare to the marker-based spatial_domain labels from Section 6
SpatialDimPlot(integrated, group.by = "banksy_domain") + ggtitle("BANKSY domains")
table(integrated$banksy_domain, integrated$spatial_domain)
```

> **Caveat inherited from `RunBanksy()`:** don't call `ScaleData()` on the `BANKSY`
> assay afterward -- it already holds the lambda-weighted BANKSY matrix in
> `scale.data`, and re-scaling would undo the effect of `lambda`.

Neither method is strictly better: `BuildMultipleNicheAssays()`'s niches are defined
by which cell types surround a spot (interpretable directly in terms of the
annotations from Section 6), while BANKSY domains are defined by a blend of a spot's
own expression and its neighbors' expression (can pick up transcriptional gradients
composition-based niches miss, at the cost of being less directly interpretable).

---

## 11. Neighborhood Enrichment — `NeighborhoodEnrichment()`

Tests pairwise cell-type co-localization by comparing observed k-NN pair frequencies to a permutation null. Reports enrichment z-scores per (source, target) pair, and optionally clusters cells into "niches" based on their neighborhood composition.

```r
enrich <- NeighborhoodEnrichment(
  integrated,
  group.by      = "domain",       # or rctd_dominant, or predicted labels
  k             = 10,
  n_perm        = 200,
  assign_niches = TRUE,
  n_niches      = 6
)
enrich$z                            # source x target z-score matrix
integrated <- enrich$obj            # now carries a 'niche' metadata column

# Visualize niches on the tissue
SpatialDimPlot(integrated, group.by = "niche")

# Enrichment heatmap
library(pheatmap)
pheatmap::pheatmap(enrich$z, cluster_rows = FALSE, cluster_cols = FALSE)
```

Complementary to `BuildMultipleNicheAssays()` — this one focuses on statistical enrichment of specific cell-type pairings; the other builds a full neighborhood assay usable in downstream Seurat workflows.

---

## 12. Niche Co-expression — `NicheCoExpress()`

Differential co-expression of a gene pair across niches using the Manders overlap coefficient, comparing two *conditions* (each needing at least 2 samples for the statistical test). Answers "is this ligand-receptor pair co-detected more often in the vascular niche in condition B than condition A?"

This is where loading all four `stxBrain` sections back in Section 2 pays off:
`sample_col = "orig.ident"` gives four samples (`anterior1`, `anterior2`, `posterior1`,
`posterior2`), and `condition_col = "region"` (set in Section 5) splits them 2-vs-2 —
exactly the minimum `NicheCoExpress()` needs to run a real test rather than a syntax
stub. As flagged back in Section 2: "region" here means anterior-vs-posterior brain
anatomy, not a treatment condition — read the result as "does this pair co-express
differently between these two brain regions," not as a stand-in for a drug/disease
comparison.

```r
DefaultAssay(integrated) <- 'RNA'
integrated <- JoinLayers(integrated)
integrated <- NormalizeData(integrated)

co <- NicheCoExpress(
  seurat_obj    = integrated,
  genes         = c("Vegfa", "Kdr"),   # or a 2-column data.frame of specific pairs
  niche_col     = "niche",             # set by NeighborhoodEnrichment() in Section 11
  sample_col    = "orig.ident",        # anterior1, anterior2, posterior1, posterior2
  condition_col = "region",            # anterior (n=2) vs posterior (n=2)
  layer         = "data"
)
co$stats               # per niche x pair: delta, statistic, p, p_adj
co$per_sample           # long-form per-sample x niche x pair co-expression scores

plotNicheCoExpress(co, type = "heatmap")   # niche x pair heatmap of delta, significance stars from p_adj
plotNicheCoExpress(co, type = "scores")    # per-sample score plots for top/selected pairs
```

### 12.1 Estimation Plot — `NicheCoExpressEstimationPlot()`

`NicheCoExpress()`'s Wilcoxon/t-test in `co$stats` answers "is there a difference."
`NicheCoExpressEstimationPlot()` answers "how big, with what uncertainty" for a
specific (niche, gene-pair) combination, using `dabestr` to show a bootstrap 95%
confidence interval on the effect size alongside the raw per-sample values -- see
9.1 for how to read that plot (raw per-sample swarm + effect size/CI panels) and
what each `effect` option means.

By default it reuses `attr(co$stats, "conditions")` for `idx`, so the reference/test
order matches `co$stats$delta` automatically — you don't need to re-specify which
region is which.

`niches` must match the actual values in `co$per_sample$niche`, not the bare cluster
number. `NeighborhoodEnrichment()` in Section 11 was called with its default
`niche_prefix = "niche_"`, so the 6 clusters requested there (`n_niches = 6`) are
labeled `"niche_1"` through `"niche_6"`, not `"1"` through `"6"` — check
`unique(co$per_sample$niche)` if you're not sure what's actually present:

```r
unique(co$per_sample$niche)

# One niche x pair combination
NicheCoExpressEstimationPlot(co, niches = "niche_1", pairs = "Vegfa_Kdr")

# Every niche x pair combination present, as a named list of plots
plots <- NicheCoExpressEstimationPlot(co)
plots[["niche_1 | Vegfa_Kdr"]]

# A different effect size -- e.g. Cliff's delta instead of the mean-difference default
NicheCoExpressEstimationPlot(co, niches = "niche_1", pairs = "Vegfa_Kdr", effect = "cliffs_delta")
```

---

## 13. Single-cell Spatial (Xenium / CosMx) — `detect_fov_edges()` / `detect_tissue_holes()`

`detect_fov_edges()` and `detect_tissue_holes()` aren't specific to imaging-based
data — Section 4 above already used both on Visium. They're especially handy for
Xenium/CosMx, though, since those platforms have no hex-grid structure at all (cells
sit at arbitrary micron coordinates), so a Visium-specific tool like
`EdgeDetectionVisium()` isn't an option there in the first place.

```r
# Cells near the FOV outer boundary (bbox method is fast and stable)
xenium <- detect_fov_edges(xenium,
                           method       = "bbox",
                           bbox_factor  = 2,
                           n_iterations = 2,
                           label_col    = "edge_layer")

# Cells bordering internal gaps / tears
xenium <- detect_tissue_holes(xenium,
                              min_hole_size = 4,
                              n_iterations  = 2,
                              label_col     = "hole_layer")
```

**Marker-gene exclusion** — biologically meaningful gaps (e.g. liver central veins expressing `Glul`, vessels expressing `Pecam1`) can be preserved:

```r
xenium <- detect_tissue_holes(xenium,
                              exclude_gene       = "Glul",
                              sensitivity        = 0.75,
                              exclude_gene_layer = "data")
```

Combine both filters:

```r
xenium <- xenium[, xenium$edge_layer == 0 & xenium$hole_layer == 0]
```

---

## 14. Subsetting Spatial Objects — `subset_opt()`

Always use `subset_opt()` instead of `subset()` for spatial Seurat objects. Seurat's
built-in `subset()` can leave stale image metadata attached after cell removal, which
causes downstream errors. `subset_opt()` keeps FOVs and images synchronized.

```r
# Isolate oligodendrocyte spots from one section
oligo_ant <- subset_opt(
  niche_list[["anterior1"]],
  subset = spatial_domain == "Oligodendrocyte"
)

cat("Oligodendrocyte spots:", ncol(oligo_ant), "\n")

# Visualize just oligodendrocytes on tissue
SpatialDimPlot(oligo_ant) + ggtitle("Oligodendrocyte spots — anterior")

# Subset by niche
white_matter_niche <- subset_opt(
  niche_list[["anterior1"]],
  subset = best_niche %in% c("1", "3")   # adjust niche IDs to match your run
)
```

`subset_opt()` also accepts an explicit `cells` vector, which is useful after
polygon-based selection:

```r
# Keep only cells in a specific niche AND expressing Mbp
mbp_oligo <- subset_opt(
  niche_list[["anterior1"]],
  cells = WhichCells(niche_list[["anterior1"]],
                     expression = spatial_domain == "Oligodendrocyte" & Mbp_pos == TRUE)
)
```

---

## 15. Tips Specific to Spatial Data

**Run `detect_fov_edges()`/`detect_tissue_holes()` before `MergeSeurat()`.** Edge and
hole-bordering spots are almost always the lowest-quality cells in the dataset.
Including them biases normalization and distorts the integration. (`EdgeDetectionVisium()`
is also available for this purpose — see the comparison note in Section 4 for when
to prefer one over the other.)

**`spatial = "Visium"` in `MergeSeurat()`.** Without this flag, the function won't
handle the `Spatial` assay correctly and image slots may be dropped or misaligned
after merging.

**Regress `percent.mt` even in Visium.** Although Visium spots contain multiple cells,
high mitochondrial read fractions still indicate damaged or necrotic tissue regions
and should be regressed.

**SCTransform is preferred over LogNormalize for Visium.** Visium spots have highly
variable total counts driven by spot cellularity. SCTransform's variance-stabilization
handles this better than library-size normalization.

**`neighbors.k = 6` for Visium niche assays.** Visium spots are arrayed in a
hexagonal grid with exactly 6 direct neighbors. Using `neighbors.k = 6` therefore
captures the immediate physical neighborhood without reaching into distant spots.

**Deconvolution for cell-type resolution.** Visium spots capture multiple cells.
The bundled `RunRCTD()` wraps `spacexr::RCTD` and writes per-cell-type proportion
columns onto every spot — this is the recommended primary annotation strategy for
Visium. `AnnotateClusters()` on Visium data works as a sanity check but loses
minority populations by construction (winner-takes-all).

**Visium-safe defaults for `AnnotateClusters()`.** If you do use marker scoring on
Visium, pass `filter_nonspecific = FALSE`, `min_detection_frac = 0.05`, and
`return_scores = "cluster"` so you can inspect the minority signal in the score
matrix.

**`check_gene_ids_across_objects()` before merging sections from different runs.**
Mouse gene symbols are consistent across standard 10x Genomics pipelines, but if you
ever mix CellRanger versions or genome builds, identifiers can drift.

```r
check_gene_ids_across_objects(brain_list)
#> All objects use: MGI symbol
```

---

## 16. Session Info

```r
sessionInfo()
```

Key packages used in this vignette:

| Package | Role |
|---|---|
| `SingleCellTools` | Edge/hole detection, integration, annotation, RCTD, niche + enrichment analysis |
| `Seurat` | Core data structure and spatial visualization |
| `SeuratData` | `stxBrain` Visium mouse brain dataset |
| `harmony` | Batch correction across sections (via `IntegrateLayers`) |
| `ClusterR` | Mini-batch k-means for niche clustering |
| `spacexr` | RCTD Visium deconvolution |
| `RANN` | Nearest-neighbor searches (`detect_fov_edges`, `detect_tissue_holes`, `NeighborhoodEnrichment`) |
| `ks` | 2D KDE for `PlotFeatureDensity` |
| `UCell` | Module scoring for `AnnotateClusters` |
| `dplyr` / `ggplot2` / `patchwork` | Data wrangling and plotting |

---

*For questions or bug reports, open an issue at
[github.com/gevensen95/SingleCellTools/issues](https://github.com/gevensen95/SingleCellTools/issues).*
