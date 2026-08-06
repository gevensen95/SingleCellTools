# Out-of-core spatial data for SingleCellTools — design doc

Status: draft, not yet implemented. Nothing in this doc has been coded.

## 1. Problem

BPCells solves one problem well: sparse count matrices that don't fit in RAM,
via a disk-backed format that Seurat v5 treats as just another layer. Spatial
data has *three* separate large-object problems, and only one of them is the
same shape as BPCells':

| Object | Example scale | Same problem as BPCells? |
|---|---|---|
| Count matrix (spots/bins x genes) | Visium HD: millions of bins | Yes — sparse matrix, BPCells applies as-is |
| Morphology image (H&E, IF) | Visium HD full-res / Xenium: multi-gigapixel | No — dense raster, not a matrix |
| Molecule table (per-transcript XY) | Xenium/MERFISH: 10⁸–10⁹ rows | No — point cloud, not a matrix |

This doc scopes the two problems BPCells doesn't already cover: images and
molecule tables.

## 2. Current state in this package (audited from source, not assumed)

**`CreateVisiumObjects()`** (`R/CreateVisiumObjects.R`, line ~91) calls
`Seurat::Read10X_Image()`, which reads the tissue image fully into memory as
an array on `obj@images$image`. For standard Visium this is the spaceranger
`hires`/`lowres` PNG (~2000 px, fine). For **Visium HD**, spaceranger also
ships a full-resolution image that can be tens of thousands of pixels per
side — this is the concrete case that breaks today, not a hypothetical one.

**`LoadXenium2()`** (`R/LoadXenium2.R`, line ~64-68), the `microns` branch:

```r
transcripts <- arrow::read_parquet(file.path(data_dir, "transcripts.parquet"))
transcripts <- subset(transcripts, qv >= mols.qv.threshold)
```

This reads the *entire* transcript table into an in-memory R data frame
before filtering. For a whole-slide Xenium run this is 10⁸+ rows and is
already the single most memory-hungry line in the package. This is the
clearest, already-existing target for a fix — independent of Visium HD.

Also found while auditing this file: `arrow`, `data.table`, and `R.utils`
are used (`arrow::read_parquet` unconditionally; `data.table`/`R.utils` via
a `requireNamespace()` fallback) but **none of the three are declared in
`DESCRIPTION`**. `arrow` isn't even guarded — if it's missing, the user gets
R's raw "there is no package called 'arrow'" instead of a clear message.
This is a pre-existing bug, not part of the original ask, but fixing the
molecule-table path means touching this code anyway, so it should be folded
into Phase 2 below rather than filed as a separate task.

**`GenerateQCReport()`**'s spatial section (line ~446-474) only calls
`Seurat::GetTissueCoordinates(image, which = "centroids")` to scatter-plot
QC metrics at spot/cell positions — it never touches the raster pixels.
That means it needs **no changes** as long as any new lazy image class
implements the `GetTissueCoordinates` generic, which Seurat's
`SpatialImage` virtual class already requires of every image subclass.

**`EdgeDetectionVisium()`** is coordinate-only (`GetTissueCoordinates` or a
CSV fallback) and never touches pixel data — out of scope entirely.

Net: the image problem is currently latent (breaks specifically on Visium
HD full-res images, not standard Visium), while the molecule-table problem
is live today in `LoadXenium2()` for any full Xenium run.

## 3. Subsystem A — lazy images (terra-backed)

**Goal:** stop loading full-resolution spatial images into an in-memory
array. Read only the tile/resolution actually needed for a given plot or
computation.

**Why terra:** it's the standard, maintained R interface to GDAL, and GDAL's
tiled/pyramidal raster format (Cloud-Optimized GeoTIFF) is the same solved
problem satellite imagery has used for years — `terra::rast()` memory-maps
a file and only reads the requested window/resolution. Reinventing this
would be strictly worse than depending on it. (Confirmed via user decision:
add `terra` as a new dependency.)

**Design choice: implement a `SpatialImage` subclass, not a side-channel.**
Seurat's image slot (`obj@images$<name>`) expects an object inheriting from
the S4 virtual class `Seurat::SpatialImage`, with generics `GetImage()`,
`GetTissueCoordinates()`, `dim()`, `Radius()`, etc. This is Seurat's actual
extension point (it's how `VisiumV1`/`VisiumV2`/`STARmap` etc. are
implemented) — a new lazy image class should be a proper subclass, not
data stashed in `@misc`, so it plugs into everything downstream
(`SpatialFeaturePlot`, this package's own plotting, `GenerateQCReport`) for
free.

Proposed:

```r
setClass("LazyVisiumImage", contains = "SpatialImage", representation(
  raster       = "ANY",    # a terra::SpatRaster handle (memory-mapped, not loaded)
  coordinates  = "data.frame",  # spot/cell positions in pixel space (small, eager)
  scale.factors = "list",
  spot.radius  = "numeric"
))

setMethod("GetImage", "LazyVisiumImage", function(object, mode = "raster", ...) {
  # Reads only the requested window/resolution via terra, on demand
})
setMethod("GetTissueCoordinates", "LazyVisiumImage", function(object, ...) {
  object@coordinates  # already in memory -- this is what GenerateQCReport uses today
})
```

New constructor, additive (does not touch existing behavior):

```r
CreateVisiumObjects(
  ...,
  image_backend = c("eager", "lazy")   # default "eager" = today's behavior, unchanged
)
```

`image_backend = "lazy"` builds a `LazyVisiumImage` from the full-res TIFF
via `terra::rast()` instead of `Read10X_Image()`. `"eager"` stays the
default so every existing script and vignette keeps working unmodified.

**Open question (needs your input before implementation):** does terra
require the image as a georeferenced raster file (GeoTIFF/COG), or can it
read spaceranger's plain PNG/JPEG tissue images directly? spaceranger's
`hires_image.png` isn't georeferenced. Likely answer: terra can read plain
PNG/TIFF as an unreferenced raster (pixel-coordinate CRS), which is
actually fine here since Visium's own coordinate system is already pixel
space — but this needs a spike against a real Visium HD output directory
before committing to the API above.

## 4. Subsystem B — lazy molecule tables (arrow-backed)

**Goal:** stop loading the full Xenium/MERFISH transcript table before
filtering it. This is the fix for the concrete bug found in section 2.

**Design:** replace the eager `arrow::read_parquet()` + `subset()` with
`arrow::open_dataset()`, which returns a lazy query object — filters
(`qv >= threshold`, gene subset, bounding-box) get pushed down and only the
matching rows are ever materialized into R memory. This is directly
analogous to how `Signac::CreateFragmentObject()` already handles ATAC
fragment files elsewhere in this package's dependency tree (tabix-indexed,
queried on demand, never fully loaded) — same pattern, different format.

```r
LoadXenium2(
  ...,
  microns_lazy = FALSE   # default FALSE = today's eager behavior, unchanged
)
```

When `microns_lazy = TRUE`, the `microns` branch returns an
`arrow::open_dataset()` connection (or a thin wrapper around one) instead
of a materialized data frame. The practical consequence: `CreateFOV()`
needs *some* molecule positions to build the object Seurat expects, so the
realistic contract is:

- Seurat's own `FOV`/`Molecules` object still gets a bounded, in-memory set
  (e.g., filtered by `qv` — already small relative to the raw table, or an
  explicit gene subset the user cares about for plotting).
- The lazy `arrow` connection is additionally attached (e.g.
  `obj@misc$molecules_lazy`, since there's no equivalent Seurat generic to
  hook into here the way `SpatialImage` exists for images) for later
  windowed/gene-subset queries without re-reading the full parquet file.

This is a real compromise, not a clean symmetrical solution like Subsystem
A — Seurat's `Molecules`/`FOV` classes assume in-memory positions, and
fighting that is out of scope. Flagging this explicitly rather than
glossing over it.

## 5. Dependency changes

- **`terra`**: new, `Imports` (needed unconditionally once `image_backend =
  "lazy"` is used) — or `Suggests` + `requireNamespace()` guard, matching
  this package's existing pattern for heavy optional deps (`sf`, `harmony`,
  `CellChat`, etc.). Recommend `Suggests`, since most existing standard
  Visium users will never touch the lazy path.
- **`arrow`**: currently used but **undeclared** — add to `DESCRIPTION`
  regardless of this feature (pre-existing bug). `Suggests` + guard, since
  `LoadXenium2()`'s `matrix`/`centroids`/`segmentations` branches don't
  need it.
- **`data.table`, `R.utils`**: same — already guarded with
  `requireNamespace()` in code but undeclared in `DESCRIPTION`. Add to
  `Suggests`.

## 6. Testing strategy

Matches this package's existing precedent (see `CreateRNAObjects()`'s
`workers` parameter, also not fixture-tested): argument validation
(`image_backend`/`microns_lazy` value checks, `match.arg` behavior) is
cheap to unit test without real data and should be. The actual lazy-read
behavior needs a real Visium HD or Xenium output directory to verify
against — out of scope for `testthat`, called out explicitly rather than
skipped silently.

## 7. Phasing

1. **Phase 0 (small, standalone, no new deps):** fix the undeclared
   `arrow`/`data.table`/`R.utils` dependency bug in `DESCRIPTION`. Useful
   regardless of whether the rest of this proceeds.
2. **Phase 1 (molecule tables):** `LoadXenium2()`'s `microns_lazy` option,
   since it fixes a live problem in already-shipped code and needs no new
   heavy system dependency (arrow's already in use).
3. **Phase 2 (images):** `LazyVisiumImage` class + `image_backend` option
   in `CreateVisiumObjects()`. Needs the terra spike (section 3's open
   question) resolved first.
4. **Phase 3:** extend the same `SpatialImage` subclass approach to Xenium
   morphology images, if Phase 2's design holds up in practice.

## 8. Explicitly out of scope for now

- Visium HD bin-matrix handling itself — already solvable today via
  BPCells + Seurat v5 layers, no new code needed here.
- Rewriting `CreateVisiumObjects()`/`LoadXenium2()`'s existing eager paths
  — both stay as the default; everything above is additive and opt-in.
- A unified "SpatialData-equivalent" storage format. The scverse
  `SpatialData` project (Zarr-backed, Python) already exists for that; this
  package targets Seurat-object interop, not a new storage standard.

## 9. Decisions still needed from you before writing code

- Phase order confirmed as 0 → 1 → 2 → 3 above, or reprioritize?
- `terra` as `Suggests` (guarded) vs `Imports` (always installed)?
- For Phase 2: does someone have a real Visium HD output directory to spike
  against, or should the terra-can-read-plain-PNG question be resolved via
  documentation/research first?
