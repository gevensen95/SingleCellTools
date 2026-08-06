# Out-of-core spatial data for SingleCellTools — design doc

Status: Phase 0, Phase 1, and Phase 2 (revised) implemented (see section 7).
Phase 3 still design-only. Section 3 was rewritten after spiking against a
real Visium HD directory — the original terra-based design was replaced by
something much simpler once the spike showed the images involved aren't
actually gigapixel-scale. The Visium HD directory-layout gap the spike
surfaced (section 9, originally) is now also fixed — see section 7's
"Visium HD directory support" entry. The image-slot management tools
(section 7, item 6) were generalized from Visium-only to also cover
Xenium/CosMx/MERFISH FOV objects, renamed `VisiumImageInfo()` →
`SpatialObjectInfo()` and `DropVisiumImage()` → `DropSpatialImage()`
accordingly (safe to rename outright — neither had been committed yet).

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
an array on `obj@images$image`. Originally assumed this breaks on Visium
HD's full-resolution image specifically — **spiked against a real Visium HD
mouse-brain output directory (section 3) and that assumption didn't hold.**
The actual problem turned out to be aggregate memory across many samples,
not any single image's resolution. See section 3 for what changed.

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

## 3. Subsystem A — deferred images (revised after a real-data spike)

**Original plan (superseded):** a terra-backed `SpatialImage` S4 subclass
implementing tiled/pyramidal windowed reads, on the assumption that Visium
HD full-resolution images are gigapixel-scale (like whole-slide-scanner
output). Spiked against a real Visium HD mouse-brain output directory to
confirm before building it. The assumption didn't hold, so this section was
rewritten rather than implemented as originally planned.

**What the spike found.** `spatial/tissue_hires_image.png` in the real
directory is 6000×5925 px (35 MP, ~103 MB decoded to a raw RGB array) — not
gigapixel. Checked inside `binned_outputs.tar.gz` too: each bin-size
subdirectory (`square_002um`/`square_008um`/`square_016um`) re-ships the
exact same images at the exact same resolution, not a bigger one. That's
every image spaceranger's own Visium HD output bundle ships. A genuinely
gigapixel image would only enter the picture if someone supplies their own
raw whole-slide-scanner source scan separately — `CreateVisiumObjects()`
doesn't ingest that today regardless of backend, eager or otherwise, so
it's not this function's problem to solve.

35 MP / ~103 MB is comfortably within what `Read10X_Image()` already
handles today for *one* sample. The couldn't-verify terra/GDAL question
(can it read a plain, non-georeferenced PNG?) became moot — not because it
was answered, but because the premise it was answering didn't apply.

**The real problem, once re-examined:** `CreateVisiumObjects()` takes a
*list* of `data_dirs` and eagerly decodes and holds a full image for every
sample before returning anything. At ~103 MB/sample, 10 samples is ~1 GB
just for images; 30-50 samples (a realistic Visium HD cohort size) is
3-5 GB — and it's wasted for the common case: `GenerateQCReport()`'s
spatial section (see section 2) never touches image pixels, only
`GetTissueCoordinates()`, yet today's code decodes every sample's image
regardless.

**Revised design — no new dependency needed.** Visium HD images aren't big
enough to need tiled/windowed reads; they're small enough to just not
*all* be decoded at once. `CreateVisiumObjects(..., image_backend =
"deferred")` attaches `tissue_lowres_image.png` (~1 MB decoded) instead of
`tissue_hires_image.png` to every sample — same
`Seurat::Read10X_Image(image.dir, image.name = ...)` call already used
today, just pointed at the small file via its existing `image.name`
parameter. Coordinates/scale-factors/spot-radius are unaffected (Seurat
computes those from the CSV/JSON files, not the image bytes), so
`EdgeDetectionVisium()` and `GenerateQCReport()` need zero changes. The
hires path gets stashed at `obj@misc$hires_image_path`; a new
`GetHiresVisiumImage(obj, attach = FALSE)` decodes it on demand for
whichever specific sample(s) you actually need full detail on, mirroring
the `QueryXeniumMolecules()` accessor pattern from Phase 1.

This deliberately avoids the S4-subclass route (a real `SpatialImage`
subclass implementing `GetImage()`/`dim()`/etc. correctly) since that has
no way to be verified in this environment (no R interpreter available to
test S4 method dispatch against the installed Seurat version) and getting
slot/generic details subtly wrong is a bad failure mode — code that loads
without error but silently renders at the wrong scale. Reusing
`Read10X_Image()` unmodified, just pointed at a different file, carries
none of that risk: every number it produces (coordinates, scale factors,
spot radius) comes from Seurat's own tested code path, unchanged.

**Known gap found during the spike — fixed as a follow-up (see section 7,
"Visium HD directory support").** `CreateVisiumObjects()`'s
directory-detection logic (`list.dirs(...)` for a folder named `spatial`)
did not work against a real Visium HD output layout. Regular Visium keeps
`tissue_positions.csv` directly inside `spatial/`; Visium HD nests it at
`binned_outputs/square_{002,008,016}um/spatial/tissue_positions.parquet`
instead — a `.parquet` file, not `.csv`, in a different subtree entirely.
Closed via a new `hd_bin_size` parameter, auto-detection of
`binned_outputs/`, and parquet support in `EdgeDetectionVisium()`'s
`coord_path` fallback.

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

- **`terra`**: not needed. Dropped after the section 3 spike showed Visium
  HD images aren't gigapixel-scale, so the tiled/windowed-read capability
  terra exists for was never actually required. No GDAL/system dependency
  added by this work.
- **`arrow`**: currently used but **undeclared** — added to `DESCRIPTION`
  regardless of this feature (pre-existing bug, Phase 0). `Suggests` +
  guard, since `LoadXenium2()`'s `matrix`/`centroids`/`segmentations`
  branches don't need it.
- **`data.table`, `R.utils`**: same — already guarded with
  `requireNamespace()` in code but were undeclared in `DESCRIPTION`. Added
  to `Suggests` (Phase 0).
- **`png`**: added to `Suggests` (Phase 2) — `GetHiresVisiumImage()` calls
  `png::readPNG()` directly. Guaranteed installed transitively via Seurat
  (which uses it internally in `Read10X_Image()`), but declaring it
  explicitly is correct practice since this package now calls it directly
  too.

## 6. Testing strategy

Matches this package's existing precedent (see `CreateRNAObjects()`'s
`workers` parameter, also not fixture-tested): argument validation
(`image_backend`/`microns_lazy` value checks, `match.arg` behavior) is
cheap to unit test without real data and should be. The actual lazy-read
behavior needs a real Visium HD or Xenium output directory to verify
against — out of scope for `testthat`, called out explicitly rather than
skipped silently.

## 7. Phasing

1. **Phase 0 (small, standalone, no new deps) — DONE.** `arrow`,
   `data.table`, `R.utils` declared in `DESCRIPTION` `Suggests`;
   `LoadXenium2()`'s `microns` branch now fails fast with a clear message
   if `arrow` isn't installed instead of an unconditional
   `arrow::read_parquet()` call.
2. **Phase 1 (molecule tables) — DONE.** `LoadXenium2(..., microns_lazy =
   TRUE)` now reads `transcripts.parquet` via `arrow::open_dataset()` with
   the QV filter/column selection pushed down to Arrow's query engine
   instead of materializing the full table, and attaches the unfiltered
   dataset connection at `obj@misc$molecules_lazy`. New
   `QueryXeniumMolecules(obj, genes=, x_range=, y_range=, qv_threshold=)`
   accessor queries that connection for a gene/region/QV subset without
   re-reading the file. `microns_lazy` defaults to `FALSE` — existing
   behavior is unchanged unless opted into. Argument-validation tests
   added; the actual lazy-read behavior against a real
   `transcripts.parquet` is untested here (needs real Xenium data),
   matching this file's existing precedent.
3. **Phase 2 (images) — DONE, redesigned after the spike.** Spiked against
   a real Visium HD directory (section 3); the terra/gigapixel premise
   didn't hold, so the terra-backed `SpatialImage` subclass was replaced
   with a much simpler fix: `CreateVisiumObjects(..., image_backend =
   "deferred")` attaches `tissue_lowres_image.png` instead of
   `tissue_hires_image.png` to every sample (same `Read10X_Image()` call,
   different `image.name`), stashing the hires path at
   `obj@misc$hires_image_path`. New `GetHiresVisiumImage(obj, attach =
   FALSE)` decodes the hires image on demand. No new dependency. Default
   stays `"eager"` — unchanged behavior unless opted into.
   Argument-validation tests added; found (but didn't fix at the time) a
   separate, real gap where `CreateVisiumObjects()` couldn't parse a
   genuine Visium HD directory layout at all (see section 3's "known
   gap") — that gap is now closed, see below.
4. **Visium HD directory support — DONE (follow-up to the "known gap"
   above).** `CreateVisiumObjects()` now auto-detects a Visium HD sample
   per-directory (presence of a `binned_outputs/` subdirectory) and
   redirects matrix/image/tissue-position reads into
   `binned_outputs/square_<hd_bin_size>um/`, via a new `hd_bin_size`
   parameter (`"008um"` default, `"002um"`/`"016um"` also available).
   Requires the archives pre-extracted, same convention as regular
   Visium. `EdgeDetectionVisium()`'s `coord_path` fallback now also reads
   `tissue_positions.parquet` (Visium HD's format) in addition to the
   existing CSV support, via `arrow::read_parquet()` (guarded,
   `Suggests`). Also fixed two pre-existing bugs found while tracing this
   code path, both independent of HD and affecting every
   `CreateVisiumObjects()` call: (1) `EdgeDetectionVisium(path_seurat,
   obj)` was called with arguments in the wrong order relative to its
   actual signature (`seurat.obj` first, `coord_path` second) — silently
   worked around today only because `coord_path`'s file-based fallback
   error path happened to still run, but would break outright as soon as
   the mismatch mattered; now called with named arguments. (2)
   `path_seurat` was built from `names(seurat_objects[i])`
   (`basename(data_dirs)`/`object_names`) instead of the actual
   `data_dirs[[i]]` path, which only worked by coincidence when a sample's
   directory happened to equal its own basename relative to the working
   directory; now resolved from `data_dirs[[i]]` directly. Tests added for
   the new bin-size detection/error path and the parquet-reading path
   (synthetic parquet fixture, not real Visium HD data).
5. **Phase 3 — not started, and no longer obviously needed.** Was
   "extend the terra `SpatialImage` subclass to Xenium images" — since
   Phase 2 dropped terra entirely, this needs re-scoping if Xenium
   morphology images turn out to have the same "aggregate multi-sample
   memory" shape as Visium HD did. Not investigated yet.
6. **Image-slot management tools — DONE (follow-up to Phase 2), then
   generalized to Xenium/CosMx/MERFISH.**
   `CreateVisiumObjects()` stashes `obj@misc$visium_image_dir`
   unconditionally (not just under `image_backend = "deferred"`). Two
   functions audit/free spatial memory:

   `SpatialObjectInfo(obj)` (originally `VisiumImageInfo()` — renamed
   during generalization, see below) reports per-sample image/FOV
   class, resolution, cell count, and memory size across a Seurat object
   or list, so you can see where memory is going without guessing.

   `DropSpatialImage(obj, mode = c("remove","downgrade"))` (originally
   `DropVisiumImage()`) frees it retroactively -- `"remove"` drops
   `@images` entirely regardless of platform (breaks anything needing
   `GetTissueCoordinates()` afterward, e.g. re-running
   `EdgeDetectionVisium()`/`GenerateQCReport()`'s spatial section or any
   FOV-based function like `NeighborhoodEnrichment()`, so only appropriate
   once spatial-coordinate work is done); `"downgrade"` rebuilds
   pixel-backed (Visium) images at lowres, keeping coordinates/QC working.
   Argument-validation and branch-logic tests added
   (no-op/already-deferred/no-images/missing-`visium_image_dir` paths);
   the actual pixel-image rebuild in `"downgrade"` mode needs a real
   Visium directory and is untested here, same precedent as everything
   else touching real spaceranger output.

   **Generalization to Xenium/CosMx/MERFISH (this pass).** Both functions
   originally only handled `VisiumV1`-class pixel images. Reworked to
   branch by *capability* rather than class name: whether
   `Seurat::GetImage(img, mode = "raw")` returns a real array (pixel-backed,
   e.g. `VisiumV1`) or not (coordinate-only, e.g. `FOV`). This makes them
   work on any `SpatialImage` subclass, not just the two seen in this
   codebase's own loaders --
   including Xenium (`LoadXenium2()`), CosMx (`Seurat::LoadNanostring()`),
   and MERFISH (built via `Seurat::CreateFOV()`, no dedicated loader in
   Seurat itself), since all of these produce standard `FOV` objects under
   Seurat's own object model. No new CosMx/MERFISH-specific loading code
   was written -- the existing FOV-based analysis functions in this package
   (`NeighborhoodEnrichment()`, `detect_fov_edges()`, `SetImageBoundary()`,
   etc.) already operate generically on any `FOV`, confirmed by their own
   source (`combine_fovs.R`/`SetImageBoundary.R` use `Boundaries()`,
   `img$centroids`, `img$molecules` unqualified, relying on `Seurat`/
   `SeuratObject` being in `Depends`); the only gap was the two
   Visium-specific memory-management helpers, now closed.

   `SpatialObjectInfo()` on an `FOV` image reports `n_cells` (via
   `Cells()`, falling back to the `centroids` boundary's cell count),
   `boundary_sets` (from `Boundaries()`), and best-effort
   `has_molecules`/`n_molecule_features`/`n_molecules` (via `img$molecules`,
   a named-by-gene list of point sets, each with an `@coords` matrix --
   the same accessor pattern `combine_fovs.R` already uses). A
   `molecules_lazy` column reports whether `obj@misc$molecules_lazy` is
   stashed (the Phase 1 Xenium arrow handle), independent of
   `has_molecules` -- the molecules already attached to an FOV are the
   QV-filtered subset either way, in-memory or not; `microns_lazy` only
   changes *how* that filtering happened during loading (query pushdown vs.
   full read), not the steady-state size of the attached object.

   `DropSpatialImage(mode = "downgrade")` on an FOV image: there's no
   lowres-equivalent concept for a coordinate-only image, so instead of
   attempting in-place S4 slot surgery to strip an FOV's `molecules` (rejected
   as too risky to do unverified -- no R interpreter here to confirm a
   stripped-down `Molecules`/`FOV` object stays valid), FOV images are left
   untouched with a message pointing at `mode = "remove"` (blunt, drops the
   whole image) or the existing Phase 1 `QueryXeniumMolecules()` /
   `microns_lazy` mechanism (the actual answer to "Xenium molecules use too
   much memory", solved separately since the QV-filtered subset attached to
   the FOV is already small relative to the raw transcript table).
   `mode = "remove"` itself needed no changes -- `o@images <- list()` was
   already platform-agnostic.

   Tested against `.make_spatial_obj()` (`tests/testthat/helper-spatial.R`),
   a real (synthetic) two-FOV `centroids`-only Seurat object already used
   by this package's other FOV-based function tests -- not a mock. No real
   CosMx or MERFISH example data exists in this repo (only the Visium HD
   directory under `Visium_HD_test/`), so the `n_molecule_features`/
   `n_molecules` best-effort path (which needs an FOV with a `molecules`
   boundary attached, not just `centroids`) is exercised structurally
   (falls through cleanly to `has_molecules = FALSE`/`NA` on the
   `centroids`-only fixture) but not against a real molecules-bearing FOV.

## 8. Explicitly out of scope for now

- Visium HD bin-matrix handling itself — already solvable today via
  BPCells + Seurat v5 layers, no new code needed here.
- Rewriting `CreateVisiumObjects()`/`LoadXenium2()`'s existing eager paths
  — both stay as the default; everything above is additive and opt-in.
- A unified "SpatialData-equivalent" storage format. The scverse
  `SpatialData` project (Zarr-backed, Python) already exists for that; this
  package targets Seurat-object interop, not a new storage standard.

## 9. Open items

Phases 0-2, Visium HD directory support, and the Xenium/CosMx/MERFISH
generalization of the image-management tools are done (section 7).
Remaining, not yet started:

- **Phase 3** (Xenium image memory): not investigated. Unclear whether
  Xenium morphology images have the same "many small images add up"
  shape Visium HD turned out to have, or a different one. Note this is
  distinct from the section-7-item-6 generalization work: that made the
  *existing* audit/drop tools platform-aware; it did not add any new
  ability to load or manage a Xenium/CosMx morphology image file (e.g.
  `morphology_focus.ome.tif`), which `LoadXenium2()` still never reads at
  all. If that turns out to matter, it's still open.
- **CosMx/MERFISH loaders**: this package still has no CosMx- or
  MERFISH-specific *loading* function (only `LoadXenium2()` for Xenium).
  The image-management generalization deliberately didn't add one --
  `SpatialObjectInfo()`/`DropSpatialImage()` work on any standard `FOV`
  object regardless of how it was built (Seurat's own `LoadNanostring()`
  for CosMx, or a hand-built `CreateFOV()` for MERFISH/anything else), so
  this was judged lower-risk than writing new, unverifiable
  platform-specific parsing code (no example data, no R interpreter to
  test against, same caution applied throughout this doc). Worth
  revisiting if there's an actual CosMx/MERFISH dataset to build and test
  a loader against.
- **`n_molecule_features`/`n_molecules` in `SpatialObjectInfo()`**: only
  exercised against a `centroids`-only synthetic FOV
  (`.make_spatial_obj()`), not a real FOV with a `molecules` boundary
  attached (the `img$molecules` named-list-of-`@coords` access pattern is
  copied from `combine_fovs.R`'s existing, presumably-working code, not
  independently verified here).
- **`segmented_outputs.tar.gz`** (Visium HD's cell-segmentation-derived
  output, separate from `binned_outputs.tar.gz`): not investigated at all.
  Unclear whether it fits this package's existing `CreateVisiumObjects()`
  model (spot/bin grid) or needs its own function closer to
  `LoadXenium2()`'s per-cell shape.
- **Untested against real data**: the Visium HD directory-navigation logic
  (section 7) was verified against the real directory structure (via `tar
  -tzf` listing) but not run end-to-end, since `binned_outputs.tar.gz`
  (8.5GB) wasn't extracted and there's no R interpreter in this
  environment to execute it. Worth a real smoke test before relying on it
  for actual analysis.
