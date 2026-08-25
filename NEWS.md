# SingleCellTools 2.8

Six commits since v2.7; 42 files, +4,040 / −250.

This release is mostly correctness work. Several long-standing bugs were
found by running tests that had been silently skipping for months — four of
them were only reachable once a blanket `skip_if()` was narrowed. One
affects the validity of previously generated results and needs action.

---

## ⚠️ Action required: re-run doublet calling

`calldoublet()` computed the expected doublet count as:

```r
nExp.poi <- round(optimal.pk * nrow(obj@meta.data))
```

`pK` is DoubletFinder's pANN **neighborhood size**, chosen from the data by
`find.pK()`. It is not a doublet rate and carries no information about how
many multiplets a run produced. Using it there made each sample's assumed
doublet rate an artifact of where its own BCmvn curve happened to peak.

Measured across one 8-sample experiment, the called-doublet rate ranged from
**0.4% to 18.4%**, rank-ordered by fitted pK rather than by anything
biological:

| | PN2 | PN1 | PBS2 | Alt1 | Alt2 | Alt3 | PN3 | PBS1 |
|---|---|---|---|---|---|---|---|---|
| fitted pK | .005 | .01 | .01 | .02 | .05 | .1 | .24 | .24 |
| % called doublet | 0.4 | 0.8 | 0.6 | 1.1 | 3.7 | 7.1 | 16.2 | 18.4 |

**Any object built with v2.7 or earlier carries doublet calls made this way**,
as does anything filtered on them. Re-run `calldoublet()` (or
`CreateRNAObjects(run_doublet_finder = TRUE)`) with an explicit
`doublet_rate`.

`doublet_rate` (default `0.075`) is now an argument: the assumed doublet
*formation* rate given how many cells were loaded. Scale it to your recovery
— for 10X Chromium roughly 0.8% per 1,000 cells recovered — rather than
accepting the default:

```r
objs <- CreateRNAObjects(dirs, doublet_rate = 0.03)   # ~4,000 cells recovered

objs <- Map(function(o, r) calldoublet(o, doublet_rate = r),
            objs, c(Alt1 = 0.035, Alt2 = 0.05, PBS1 = 0.04))
```

---

## Breaking changes

- **`calldoublet()` gained `doublet_rate`** and no longer derives `nExp` from
  `pK`. See above.
- **DoubletFinder's `pANN_<pN>_<pK>_<nExp>` column is renamed to
  `doublet_pANN`.** The old name embedded per-sample fitted values, so every
  object in a list carried a differently-named column — which broke any
  function that stacks metadata across samples. Code referencing the old
  names must be updated. The classification column is now selected by
  `grep("^DF.classifications")` rather than by position.
- **`assign_cell_cycle_phase()`'s G1 boundary is now `<=`, not `<`.** Phase
  calls will change on sparse data. See Bug fixes.
- **`ClusterR` moved from Imports to Suggests.** Only
  `BuildMultipleNicheAssays()` used it; that function now checks for it and
  fails fast with a clear message. Users who need niche clustering must
  install it explicitly.

## Installation

- **`DESCRIPTION` gained `Remotes:` entries for `DoubletFinder` and
  `liana`.** Both are `Imports` and neither is on CRAN or Bioconductor, so
  dependency resolution could not find them — breaking
  `devtools::install_github()`, `R CMD check`, and `build_vignettes()` for
  anyone who didn't already have them installed. `Remotes:` previously
  covered only the GitHub packages in `Suggests`.

## New features

- **`CombineParseRounds()` / `CombineCellRangerRounds()`** — combine two
  sequencing rounds of the *same* library by cell barcode, over a shared
  union/reindex/sum core. Barcodes in both rounds have counts summed; genes
  are unioned. Output is written in each platform's native format, ready to
  feed straight back into `MakeParseObj()` / `CreateRNAObjects()`.

  Note the documented caveat: summing two independently UMI-deduplicated
  matrices double-counts any molecule sequenced in both rounds. If you still
  have the FASTQs, `cellranger count --fastqs=/round1,/round2` (or
  concatenating and re-running `split-pipe`) is exact and preferable.

- **`TopMarkerPlot()`** — one call from a clustered object to a labeled marker
  panel: `FindAllMarkers()`, top *n* per cluster, then `MarkerPlot()`, so
  styling and auto-sizing match every other marker figure. Pass an existing
  markers table to skip the slow step, or `group_by` a metadata column to
  label by cell type. Genes topping two clusters are assigned to their
  best-scoring one; numeric cluster IDs are zero-padded past 9 so
  `MarkerPlot()`'s alphabetical group ordering matches numeric order.

- **`assign_cell_cycle_phase()` now works with no gene lists.** Defaults to
  `Seurat::cc.genes.updated.2019`, and matches symbols to the object's
  features case-insensitively — so human (`MCM5`) and mouse (`Mcm5`) both
  work with no `species` argument and no ortholog table to maintain. Verified
  against a 24,593-row human/mouse ortholog map: all 97 cell-cycle genes
  agree exactly, with no one-to-many ambiguity. Also gained an `assay`
  argument and a phase-count summary.

- **`MergeSeurat()` gained `marker_n`** (default 10) and now draws
  `marker_plot.pdf` through `TopMarkerPlot()`, so it gets right-edge cluster
  labels and the same content-aware sizing instead of a fixed 10×10 inch page.

- **`QCComparePlots()` gained `counts_size`** for the retention annotations.

- **`CreateRNAObjects()` gained `doublet_rate`**, `doublet_pk_sweep_max_cells`
  and `doublet_sweep_cores` passthroughs, plus an oversubscription guard that
  errors when `workers * doublet_sweep_cores` exceeds the core count.

## Bug fixes

### `calldoublet()`

- **Layer stripping was a silent no-op on Seurat v5.**
  `obj[["RNA"]][["data"]] <- NULL` neither errors nor removes anything on an
  `Assay5`, so every returned object still carried its `data` and
  `scale.data` layers despite the documentation promising otherwise. Now uses
  `LayerData()<-` and matches v5 split variants (`data.1`, `scale.data.2`).
  On a small test object this is a 56% size reduction; on real samples those
  dense layers dwarf the sparse counts they sit beside.
- **`min(co1, co2)` returned `NA` whenever either PC-selection heuristic came
  up empty**, aborting the function even when the other had found a perfectly
  usable cutoff. Now falls back to whichever succeeded, erroring only when
  both fail.
- **A degenerate pANN distribution surfaced as
  `missing value where TRUE/FALSE needed`** from deep inside
  `KernSmooth::bkde`, naming neither DoubletFinder nor the sample. Now
  translated into an explanation with a suggested fix.

### `QCComparePlots()`

- **Failed with `names do not match previous names`** when stacking a list of
  objects whose metadata carried per-sample column names (most often the
  `pANN_*` column above). List elements are now stacked on their common
  columns, with one warning naming what was dropped.
- **Retention annotations were clipped at the panel edge.** The label is text
  drawn in absolute units, so the y scale never expanded to fit it. Fixed
  with scale expansion sized to the label's line count plus `clip = "off"`;
  labels are also stacked over three lines so they stop colliding at 8+
  samples.

### `ApplyQCFilters()`

- **A `(sample, metric)` row with `suggest_lo == suggest_hi` kept only cells
  matching that exact value.** This arises whenever a metric is constant or
  near-constant across a sample, making the suggested quantiles coincide — so
  the samples needing filtering least were the ones being gutted. Degenerate
  rows are now skipped for that sample only, leaving other metrics in force.
- **A missing bound poisoned the mask.** `v >= NA` is `NA`, so `keep` carried
  `NA`s, `sum(keep)` became `NA`, and `cells[keep]` produced `NA` cell names
  that went on to `subset()`. Now skipped.
- **Inverted bounds (`lo > hi`)** would keep zero cells; now skipped with a
  warning, since it means the cutoffs table itself is malformed.

### `assign_cell_cycle_phase()`

- **No cell was ever called G1 on sparse data.** UCell returns exactly `0`
  for a cell expressing none of a signature's genes; on shallow data that is
  a large share of cells, often enough that the median — the default cutoff —
  is itself `0`. With the strict `score < threshold` test, `0 < 0` is never
  true, so every all-zero cell fell through both comparisons to `NA`. On a
  60%-zero sample that meant **0 G1 cells and 60% NA**, where the answer
  should be 60% G1. The boundary is now `<=`.
- **Re-running scored off stale columns.** The score-column `grep` matched
  both the previous run's columns and the fresh ones, and silently used the
  first. Prior columns are now cleared, making re-runs idempotent.
- Added a `requireNamespace("UCell")` guard.

### `PlotFeatureDensity()`

- **The joint panel was wrong, not just noisy.** Per-feature density grids are
  trimmed to different lengths before being multiplied, so R recycled and
  paired unrelated grid points. The product is now computed on the untrimmed
  grid and thresholded afterwards.

### `PseudobulkDE()`

- **`sample_col = "orig.ident"` was silently overwritten.**
  `AggregateExpression(return.seurat = TRUE)` rebuilds via
  `CreateSeuratObject()`, which always repopulates `orig.ident` with the new
  object's colnames. Since `orig.ident` is Seurat's own default sample column,
  this was the single most likely value a caller would pass. Aggregation now
  goes through guaranteed-safe temporary columns.
- **One-vs-rest designs produced duplicate colnames** whenever a sample
  contributed to both condition arms — previously masked by the `orig.ident`
  bug, which made colnames unique for the wrong reason. Colliding entries are
  disambiguated as `<sample>_<condition>`; `coldata$sample` still carries the
  undecorated label.

### `MergeSeurat()`

- `marker_plot.pdf` was a fixed 10×10 inches regardless of marker count — up
  to 100 gene labels in 10 inches. Now sized to content via `TopMarkerPlot()`,
  and skipped cleanly (rather than aborting the whole function at its last
  step) when no gene survives filtering.

### combine-rounds core

Four bugs in the shared core, all found by tests written before the code was
ever executed:

- Metadata columns present only in round 2 were overwritten with `NA` on every
  shared barcode.
- Seurat's `nCount_RNA` / `nFeature_RNA` survived as stale values, because the
  count-derived blocklist matched `"nCount"` exactly rather than as a prefix.
- Factor columns whose level sets differed between rounds aborted the combine.
- Round-2-only gene rows were filled by position, silently transposing values
  when the two rounds' gene tables had columns in a different order.

## Documentation

- **New QC vignette** — `vignette("SingleCellTools_QC_vignette")`. Part one is
  single-cell (report → edit cutoffs → apply → verify, doublets, batch
  effects, cell cycle); part two is spatial (object inspection, edge and hole
  detection, spatial concordance, image-safe subsetting).
- **The main vignette is now a real package vignette** — converted from a
  root-level `.md` to `vignettes/SingleCellTools_vignette.Rmd`, discoverable
  via `vignette()` / `browseVignettes()`. Added coverage for
  `PlotQCMetrics()`, `GetHiresVisiumImage()`, `combine_fovs()` and
  `SetImageBoundary()`, which had no worked example anywhere.
- **README dependency table rebuilt** against DESCRIPTION's actual
  Imports/Suggests split — DESeq2, SummarizedExperiment, edgeR, speckle,
  liana and OmnipathR are hard Imports, not the optional Suggests previously
  documented. New rows for `TopMarkerPlot()`, `CombineParseRounds()`,
  `CombineCellRangerRounds()`, and updated entries for `MergeSeurat()`,
  `calldoublet()` and `assign_cell_cycle_phase()`.

## Testing

- **Narrowed the blanket `skip_if(inherits(result, "error"))` guards.** These
  turned any regression into a silent green run; four of the bugs above were
  hiding behind one of them. Skips are now conditioned on the specific
  environmental cause they were meant to cover.
- **`calldoublet()`'s smoke test now actually runs.** Its fixture was uniform
  Poisson noise, which produces a flat PC stdev curve, so the PC-selection
  heuristics found nothing and the test skipped every time. It now uses a
  fixture with real per-cluster structure, verified across all 159 (pN, pK)
  grid points `paramSweep()` visits.
- **Test fixtures build sparse matrices** rather than handing dense ones to
  `CreateSeuratObject()`. Nearly every test in the suite builds a fixture, so
  the resulting `Data is of class matrix. Coercing to dgCMatrix.` warnings
  accounted for the large majority of a full run's output — enough to bury
  the real ones.
- Regression tests added for every fix above.
