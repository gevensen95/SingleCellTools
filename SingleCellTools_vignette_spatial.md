# SingleCellTools Vignette: Spatial Transcriptomics with Mouse Brain Visium Data

**Package:** `SingleCellTools`  
**Data:** `stxBrain` — 10x Genomics Visium mouse brain sections (anterior + posterior), available via `SeuratData`  
**Goal:** Walk through a complete Visium workflow — edge detection, integration, annotation, niche analysis — using freely available public data.

---

## Table of Contents

1. [Setup](#1-setup)
2. [Load the Data](#2-load-the-data)
3. [QC and Percent Mitochondrial Reads](#3-qc-and-percent-mitochondrial-reads)
4. [Edge Detection — `EdgeDetectionVisium()`](#4-edge-detection--edgedetectionvisium)
5. [Merge and Integrate — `MergeSeurat()`](#5-merge-and-integrate--mergeseurat)
6. [Annotate Spatial Domains — `MarkerPlot()`](#6-annotate-spatial-domains--markerplot)
7. [Gene Positivity — `AddGenePositivity()`](#7-gene-positivity--addgenepositivity)
8. [Spatial Niche Analysis — `BuildMultipleNicheAssays()`](#8-spatial-niche-analysis--buildmultiplenicheassays)
9. [Subsetting Spatial Objects — `subset_opt()`](#9-subsetting-spatial-objects--subset_opt)
10. [Tips Specific to Visium Data](#10-tips-specific-to-visium-data)
11. [Session Info](#11-session-info)

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

The `stxBrain` dataset contains two coronal mouse brain Visium sections captured from
adjacent positions along the anterior–posterior axis. Having two sections from the same
tissue type makes this an ideal test case for multi-sample spatial integration.

```r
# Load both sections
brain_ant <- LoadData("stxBrain", type = "anterior1")
brain_pos <- LoadData("stxBrain", type = "posterior1")

# Quick overview
brain_ant
#> An object of class Seurat
#> 31053 features across 2696 samples within 1 assay
#> Active assay: Spatial (31053 features, 0 variable features)
#>  1 spatial image present: anterior1

brain_pos
#> An object of class Seurat
#> 31053 features across 3353 samples within 1 assay
#> Active assay: Spatial (31053 features, 0 variable features)
#>  1 spatial image present: posterior1
```

> **A note on Visium resolution.** Each spot in a Visium capture area is ~55 µm in
> diameter — in brain tissue this typically captures the expression of 2–10 cells
> simultaneously. Clustering on Visium data therefore identifies **spatial domains**
> (brain regions with a coherent transcriptional signature) rather than individual cell
> types. Deconvolution methods (e.g., RCTD, SPOTlight) can estimate cell-type
> composition per spot, but are outside the scope of this vignette.

Visualize the raw tissue sections to confirm the data loaded correctly:

```r
SpatialFeaturePlot(brain_ant, features = "nCount_Spatial") +
  ggtitle("Anterior — total UMI per spot")

SpatialFeaturePlot(brain_pos, features = "nCount_Spatial") +
  ggtitle("Posterior — total UMI per spot")
```

---

## 3. QC and Percent Mitochondrial Reads

```r
# Mouse mitochondrial genes: "^mt-" (lowercase)
brain_ant[["percent.mt"]] <- PercentageFeatureSet(brain_ant, pattern = "^mt-")
brain_pos[["percent.mt"]] <- PercentageFeatureSet(brain_pos, pattern = "^mt-")

# QC distributions
VlnPlot(brain_ant,
        features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
        pt.size  = 0.1) &
  theme(axis.text.x = element_blank())

VlnPlot(brain_pos,
        features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
        pt.size  = 0.1) &
  theme(axis.text.x = element_blank())
```

Spatial scatter to see if low-quality spots are spatially clustered (often they are,
at tissue edges or tears):

```r
SpatialFeaturePlot(brain_ant, features = "percent.mt") +
  scale_fill_gradientn(colors = c("grey90", "red3")) +
  ggtitle("Percent mitochondrial reads — anterior")
```

Apply filters. Visium thresholds are typically looser than single-cell because spots
have higher total counts:

```r
brain_ant <- subset(brain_ant,
  nFeature_Spatial > 200  &
  nFeature_Spatial < 8000 &
  percent.mt < 25
)

brain_pos <- subset(brain_pos,
  nFeature_Spatial > 200  &
  nFeature_Spatial < 8000 &
  percent.mt < 25
)

cat("Anterior spots remaining:", ncol(brain_ant), "\n")
cat("Posterior spots remaining:", ncol(brain_pos), "\n")
```

---

## 4. Edge Detection — `EdgeDetectionVisium()`

Spots at the capture-area boundary, tissue edge, and tissue tears have systematically
lower UMI counts and higher noise. `EdgeDetectionVisium()` runs four iterative rounds
of nearest-neighbor filtering to flag these spots.

The function expects a path to a directory containing the Visium `tissue_positions_list.csv`
file. When loading data via `SeuratData`, we reconstruct this file from the coordinates
already stored in the Seurat object's image slot.

```r
# Helper: extract Visium spot coordinates from a Seurat object
# and write them in the format EdgeDetectionVisium() expects.
write_visium_coords <- function(obj, image_name, out_dir) {
  # The @coordinates data frame has columns: tissue, row, col, imagerow, imagecol
  coords_raw <- obj@images[[image_name]]@coordinates

  edge_input <- data.frame(
    barcode            = rownames(coords_raw),
    in_tissue          = coords_raw$tissue,
    array_row          = coords_raw$row,
    array_col          = coords_raw$col,
    pxl_row_in_fullres = coords_raw$imagerow,
    pxl_col_in_fullres = coords_raw$imagecol
  )

  # Write WITHOUT a header row — EdgeDetectionVisium reads with header = FALSE
  write.table(
    edge_input,
    file      = file.path(out_dir, "tissue_positions_list.csv"),
    sep       = ",",
    row.names = FALSE,
    col.names = FALSE,
    quote     = FALSE
  )
  invisible(out_dir)
}

# Write coordinates for each section
ant_dir <- file.path(tempdir(), "anterior1")
pos_dir <- file.path(tempdir(), "posterior1")
dir.create(ant_dir, showWarnings = FALSE)
dir.create(pos_dir, showWarnings = FALSE)

write_visium_coords(brain_ant, "anterior1", ant_dir)
write_visium_coords(brain_pos, "posterior1", pos_dir)
```

> **Working with real CellRanger output?** Skip the helper above and point
> `coord_path` directly to your `outs/spatial/` folder, which already contains
> `tissue_positions_list.csv`.

Now run edge detection. Four filter iterations are returned as columns
`Filter`, `Filter2`, `Filter3`, `Filter4` — each one is progressively more
aggressive at peeling back the boundary.

```r
edge_ant <- EdgeDetectionVisium(
  coord_path = ant_dir,
  seurat.obj = brain_ant,
  search     = "radius",
  neighbors  = 7
)

edge_pos <- EdgeDetectionVisium(
  coord_path = pos_dir,
  seurat.obj = brain_pos,
  search     = "radius",
  neighbors  = 7
)

# Inspect how many spots are flagged at each iteration
table(edge_ant$Filter)   # iteration 1 (least aggressive)
table(edge_ant$Filter4)  # iteration 4 (most aggressive)
```

Add the filter results to metadata and visualize before committing to a cutoff:

```r
# Add Filter4 (outermost 4 rings removed) as metadata
brain_ant <- AddMetaData(
  brain_ant,
  metadata = setNames(edge_ant$Filter4, edge_ant$barcode),
  col.name = "edge_filter"
)
brain_pos <- AddMetaData(
  brain_pos,
  metadata = setNames(edge_pos$Filter4, edge_pos$barcode),
  col.name = "edge_filter"
)

# Visualize — spots to remove shown in red
SpatialDimPlot(brain_ant, group.by = "edge_filter",
               cols = c("Keep" = "grey80", "Filter" = "red3")) +
  ggtitle("Edge detection — anterior (Filter4)")

SpatialDimPlot(brain_pos, group.by = "edge_filter",
               cols = c("Keep" = "grey80", "Filter" = "red3")) +
  ggtitle("Edge detection — posterior (Filter4)")
```

If the red spots align with the visible tissue boundary and tears, apply the filter:

```r
brain_ant <- subset(brain_ant, edge_filter == "Keep")
brain_pos <- subset(brain_pos, edge_filter == "Keep")

cat("Anterior spots after edge removal:", ncol(brain_ant), "\n")
cat("Posterior spots after edge removal:", ncol(brain_pos), "\n")
```

> **How conservative to be?** `Filter` (1 iteration) removes only the outermost ring;
> `Filter4` removes the outer four rings. Start with `Filter2` or `Filter3` for most
> experiments — it captures tissue-edge effects without discarding too much of the
> boundary. Validate by comparing UMI distributions of "Keep" vs. "Filter" spots.

---

## 5. Merge and Integrate — `MergeSeurat()`

With two clean, filtered sections, we are ready to merge and integrate. Pass
`spatial = "Visium"` so `MergeSeurat()` handles the Spatial assay correctly and
keeps images in sync.

```r
brain_list <- list(anterior1 = brain_ant, posterior1 = brain_pos)

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
```

After the run your working directory contains:
- `brain_visium_merged_seurat_objects.rds`
- `dimplot_seurat_clusters.pdf` — UMAP colored by cluster
- `markers_all.csv` and `marker_plot.pdf`

### Inspect the integration

```r
# UMAP — good integration means anterior and posterior spots overlap
DimPlot(integrated, group.by = "orig.ident",
        cols = c("anterior1" = "steelblue", "posterior1" = "tomato")) +
  ggtitle("By section — should be interleaved")

DimPlot(integrated, label = TRUE, repel = TRUE) +
  ggtitle("Seurat clusters")

# Spatial view of clusters on each section
SpatialDimPlot(integrated, label = FALSE) +
  ggtitle("Clusters mapped back to tissue")
```

---

## 6. Annotate Spatial Domains — `MarkerPlot()`

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

## 7. Gene Positivity — `AddGenePositivity()`

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

---

## 8. Spatial Niche Analysis — `BuildMultipleNicheAssays()`

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
#> $posterior1
#> [1] "posterior1"

niche_list <- BuildMultipleNicheAssays(
  list.object  = brain_list_annotated,
  list.fov     = list("anterior1", "posterior1"),  # must match names(obj@images)
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

---

## 9. Subsetting Spatial Objects — `subset_opt()`

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

## 10. Tips Specific to Visium Data

**Run `EdgeDetectionVisium()` before `MergeSeurat()`.** Edge spots are almost always
the lowest-quality cells in the dataset. Including them biases normalization and
distorts the integration.

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
If you need per-cell-type maps rather than spatial domain maps, consider deconvolution
tools such as RCTD (`spacexr`), SPOTlight, or cell2location — these can be run
after `MergeSeurat()` on the integrated object.

**`check_gene_ids_across_objects()` before merging sections from different runs.**
Mouse gene symbols are consistent across standard 10x Genomics pipelines, but if you
ever mix CellRanger versions or genome builds, identifiers can drift.

```r
check_gene_ids_across_objects(brain_list)
#> All objects use: MGI symbol
```

---

## 11. Session Info

```r
sessionInfo()
```

Key packages used in this vignette:

| Package | Role |
|---|---|
| `SingleCellTools` | Edge detection, integration, annotation, niche analysis |
| `Seurat` | Core data structure and spatial visualization |
| `SeuratData` | `stxBrain` Visium mouse brain dataset |
| `harmony` | Batch correction across sections (via `IntegrateLayers`) |
| `ClusterR` | Mini-batch k-means for niche clustering |
| `RANN` | Nearest-neighbor searches in `EdgeDetectionVisium` |
| `dplyr` / `ggplot2` / `patchwork` | Data wrangling and plotting |

---

*For questions or bug reports, open an issue at
[github.com/gevensen95/SingleCellTools/issues](https://github.com/gevensen95/SingleCellTools/issues).*
