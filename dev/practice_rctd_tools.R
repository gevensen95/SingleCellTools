#!/usr/bin/env Rscript
# =============================================================================
# Practice script: RunRCTD() multi-sample + QC fields, CompositionalTest()
# weight_cols mode, and SpatialCompositionPlot().
#
# Walks through all three RCTD-related additions on the stxBrain public
# Visium dataset (4 sections via SeuratData), paired with a single-cell
# reference you supply. See the CONFIG block below.
#
# STAGES
#   1. Load the 4 stxBrain sections (kept UNMERGED on purpose, to exercise
#      RunRCTD()'s new list-input path).
#   2. Load/build a reference with cell-type labels.
#   3. Run RunRCTD() on the list of 4 samples at once (mode = "doublet"),
#      confirming the reference is built once and reused per sample.
#   4. Inspect the new QC columns (rctd_spot_class, rctd_max_weight),
#      then re-run one sample under mode = "full" to see rctd_spot_class
#      correctly absent there.
#   5. Merge the 4 samples and run RunRCTD() again on the merged object,
#      exercising the multi-image coordinate-gathering fix.
#   6. Compare CompositionalTest() in discrete (rctd_dominant) vs.
#      continuous (weight_cols) mode.
#   7. Visualize with SpatialCompositionPlot() vs. plain SpatialDimPlot().
#
# REQUIREMENTS
#   SingleCellTools, Seurat, SeuratData (for stxBrain), spacexr (for
#   RunRCTD itself), scatterpie (for SpatialCompositionPlot), betareg
#   (for the continuous CompositionalTest example). The script checks for
#   each and skips/warns rather than hard-failing where it reasonably can.
#
# USAGE
#   1. Edit the CONFIG block below -- at minimum, point `reference_path` at
#      a single-cell reference RDS with a cell-type metadata column, or
#      set `use_pbmc3k_fallback = TRUE` to just exercise the code mechanics
#      without a biologically matched reference (see CONFIG comments).
#   2. Rscript dev/practice_rctd_tools.R
#      (or source() interactively from the repo root)
# =============================================================================

## ---- CONFIG -- edit before running -----------------------------------------

# Path to an RDS of a single-cell Seurat object with a cell-type metadata
# column. DO NOT use the full raw Allen Institute atlas here -- it's
# enormous (order of a million+ cells) and impractical to load locally.
# Seurat's own spatial vignette pairs stxBrain with a much smaller,
# pre-subsetted "allen_cortex.rds" they host specifically for this
# tutorial (order of 10-15k cortex cells, already has a cell-type column)
# -- that's the file to use here, not the full atlas:
#   https://www.dropbox.com/s/cuowvm4vrf65pvq/allen_cortex.rds?dl=1
# Download it once, save locally, and point this at that path. If even
# that feels big for your machine, max_cells_per_ref_celltype below
# downsamples it further after loading -- but the load itself is what
# was hurting you with the full atlas, and this file avoids that.
reference_path <- "~/Downloads/allen_cortex.rds"

# Column in the reference's metadata holding cell-type labels. In the
# allen_cortex.rds tutorial file above, Seurat's own vignette code uses
# "subclass" for this -- if you're using a different reference, check
# colnames(readRDS(reference_path)@meta.data) and adjust.
reference_celltype_col <- "class"

# If TRUE, ignores reference_path and uses SeuratData's pbmc3k as a stand-in
# reference instead -- purely to exercise RunRCTD()'s code path end to end
# when you don't have a real brain reference on hand yet. The result will
# NOT be biologically meaningful (PBMC genes barely overlap brain genes),
# but every function call below still runs the same way.
use_pbmc3k_fallback <- FALSE

# RCTD tuning (see the spatial vignette's Section 7.1 for the tradeoffs)
max_cells_per_ref_celltype <- 5000
n_cores                    <- 4

## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(SingleCellTools)
  library(Seurat)
})

if (!requireNamespace("spacexr", quietly = TRUE)) {
  stop("This script needs 'spacexr'. Install with: ",
       "remotes::install_github('dmcable/spacexr')")
}
if (!requireNamespace("SeuratData", quietly = TRUE)) {
  stop("This script needs 'SeuratData'. Install with: ",
       "remotes::install_github('satijalab/seurat-data')")
}

## ---- Stage 1: load the 4 stxBrain sections, unmerged -----------------------
message("=== Stage 1: loading stxBrain sections (unmerged) ===")

if (!"stxbrain.SeuratData" %in% rownames(utils::installed.packages())) {
  SeuratData::InstallData("stxBrain")
}

section_names <- c("anterior1", "anterior2", "posterior1", "posterior2")
visium_list <- setNames(
  lapply(section_names, function(s) SeuratData::LoadData("stxBrain", type = s)),
  section_names
)
message(sprintf("Loaded %d sections: %s", length(visium_list),
                paste(names(visium_list), collapse = ", ")))

## ---- Stage 2: reference ------------------------------------------------------
message("=== Stage 2: loading reference ===")

if (isTRUE(use_pbmc3k_fallback)) {
  message("use_pbmc3k_fallback = TRUE -- using pbmc3k as a mechanics-only ",
         "stand-in reference (NOT biologically matched to brain tissue).")
  if (!"pbmc3k.SeuratData" %in% rownames(utils::installed.packages())) {
    SeuratData::InstallData("pbmc3k")
  }
  ref <- SeuratData::LoadData("pbmc3k", type = "final")
  ref <- Seurat::UpdateSeuratObject(ref)
  reference_celltype_col <- "seurat_annotations"
  if (!reference_celltype_col %in% colnames(ref@meta.data)) {
    stop("pbmc3k fallback: expected a 'seurat_annotations' column and didn't ",
         "find one -- your installed pbmc3k.SeuratData version may differ. ",
         "Set reference_celltype_col to whatever cell-type column it has.")
  }
  ref <- ref[, !is.na(ref@meta.data[[reference_celltype_col]])]
} else {
  if (!file.exists(reference_path)) {
    stop("reference_path ('", reference_path, "') doesn't exist. Set it to a ",
         "real reference RDS in the CONFIG block, or set ",
         "use_pbmc3k_fallback = TRUE to practice the mechanics without one.")
  }
  ref <- readRDS(reference_path)
  if (!inherits(ref, "Seurat")) {
    stop("reference_path must point to an RDS containing a Seurat object.")
  }
  if (!reference_celltype_col %in% colnames(ref@meta.data)) {
    stop("Column '", reference_celltype_col, "' not found in the reference's ",
         "metadata. Set reference_celltype_col to the correct column name.")
  }
}
message(sprintf("Reference: %d cells, %d cell types (column '%s').",
                ncol(ref), length(unique(ref@meta.data[[reference_celltype_col]])),
                reference_celltype_col))

## ---- Stage 3: RunRCTD() on the list of 4 samples at once -------------------
message("=== Stage 3: RunRCTD() on all 4 sections (list input, mode = 'doublet') ===")
message("Watch for ONE 'Building RCTD reference' message followed by 4 tagged ",
       "'Building query puck'/'Running RCTD' blocks -- confirms the reference ",
       "is built once and reused, not rebuilt per sample.")

visium_list <- RunRCTD(
  visium_list,
  reference                  = ref,
  celltype_col                = reference_celltype_col,
  mode                        = "doublet",
  max_cells_per_ref_celltype  = max_cells_per_ref_celltype,
  n_cores                     = n_cores
)

## ---- Stage 4: inspect the new QC columns ------------------------------------
message("=== Stage 4: inspecting rctd_spot_class / rctd_max_weight ===")

for (s in names(visium_list)) {
  cat(sprintf("\n--- %s ---\n", s))
  print(table(visium_list[[s]]$rctd_spot_class))
  print(summary(visium_list[[s]]$rctd_max_weight))
}

message("Re-running anterior1 under mode = 'full' to confirm rctd_spot_class ",
       "is absent there (no spacexr-computed equivalent in full mode) while ",
       "rctd_max_weight is still present (it's our own heuristic, every mode).")
anterior1_full <- RunRCTD(
  visium_list[["anterior1"]],
  reference                  = ref,
  celltype_col                = reference_celltype_col,
  mode                        = "full",
  max_cells_per_ref_celltype  = max_cells_per_ref_celltype,
  n_cores                     = n_cores
)
stopifnot(
  !"rctd_spot_class" %in% colnames(anterior1_full@meta.data),
  "rctd_max_weight" %in% colnames(anterior1_full@meta.data)
)
message("Confirmed: rctd_spot_class absent, rctd_max_weight present under mode = 'full'.")

## ---- Stage 5: merge, then RunRCTD() again on the merged object -------------
message("=== Stage 5: merging all 4 sections, then RunRCTD() on the merged object ===")
message("This exercises the multi-image coordinate fix -- RunRCTD() gathers ",
       "tissue coordinates from every image in obj@images, not just the first, ",
       "so every section's spots should be included below.")

integrated <- merge(visium_list[[1]], visium_list[-1], add.cell.ids = names(visium_list))
integrated[["Spatial"]] <- SeuratObject::JoinLayers(integrated[["Spatial"]])
integrated$region <- ifelse(grepl("^anterior", integrated$orig.ident),
                            "anterior", "posterior")

integrated <- RunRCTD(
  integrated,
  reference                  = ref,
  celltype_col                = reference_celltype_col,
  mode                        = "full",
  max_cells_per_ref_celltype  = max_cells_per_ref_celltype,
  n_cores                     = n_cores
)
message(sprintf("Merged object: %d spots across %d sections, %d rctd_weights rows.",
                ncol(integrated), length(unique(integrated$orig.ident)),
                nrow(integrated@misc$rctd_weights)))
stopifnot(nrow(integrated@misc$rctd_weights) == ncol(integrated))
message("Confirmed: every spot across all 4 sections got an RCTD call.")

## ---- Stage 6: discrete vs. continuous CompositionalTest() ------------------
message("=== Stage 6: CompositionalTest() -- discrete (rctd_dominant) vs. continuous (weight_cols) ===")

res_discrete <- CompositionalTest(
  integrated,
  cluster_col   = "rctd_dominant",
  sample_col    = "orig.ident",
  condition_col = "region"
)
cat("\n--- Discrete (rctd_dominant) ---\n")
print(res_discrete[, c("cluster", "effect", "pvalue", "padj", "method")])

celltypes <- grep("^rctd_", colnames(integrated@meta.data), value = TRUE)
celltypes <- setdiff(celltypes, c("rctd_dominant", "rctd_max_weight", "rctd_spot_class"))

res_continuous <- tryCatch(
  CompositionalTest(
    integrated,
    weight_cols   = celltypes,
    sample_col    = "orig.ident",
    condition_col = "region",
    method        = "betareg"
  ),
  error = function(e) {
    message("betareg backend unavailable (", conditionMessage(e), "); ",
           "falling back to method = 'wilcox'.")
    CompositionalTest(
      integrated,
      weight_cols   = celltypes,
      sample_col    = "orig.ident",
      condition_col = "region",
      method        = "wilcox"
    )
  }
)
cat("\n--- Continuous (weight_cols) ---\n")
print(res_continuous[, c("cluster", "effect", "pvalue", "padj", "method")])

sig_discrete   <- res_discrete$cluster[!is.na(res_discrete$padj) & res_discrete$padj < 0.05]
sig_continuous <- res_continuous$cluster[!is.na(res_continuous$padj) & res_continuous$padj < 0.05]
message(sprintf(
  "Significant (padj < 0.05) -- discrete: %s | continuous: %s",
  if (length(sig_discrete)) paste(sig_discrete, collapse = ", ") else "(none)",
  if (length(sig_continuous)) paste(sig_continuous, collapse = ", ") else "(none)"
))

## ---- Stage 7: visualize ------------------------------------------------------
message("=== Stage 7: SpatialCompositionPlot() vs. SpatialDimPlot() ===")

if (requireNamespace("scatterpie", quietly = TRUE)) {
  p_dominant <- Seurat::SpatialDimPlot(integrated, group.by = "rctd_dominant",
                                       images = "anterior1")
  p_pie      <- SpatialCompositionPlot(subset(integrated, orig.ident == "anterior1"),
                                       n_spots_max = 500)

  outdir <- if (nzchar(Sys.getenv("RSTUDIO"))) getwd() else tempdir()
  pdf_path <- file.path(outdir, "rctd_practice_plots.pdf")
  grDevices::pdf(pdf_path, width = 8, height = 6)
  print(p_dominant)
  print(p_pie)
  grDevices::dev.off()
  message(sprintf(
    "Wrote %s -- compare the dominant-label plot (one color per spot) ",
    pdf_path),
    "against the pie plot (full mixture per spot). Spots that look ",
    "confidently single-type in the first plot but show a visible second ",
    "color in the pie plot are exactly what rctd_max_weight is meant to flag ",
    "-- cross-check low rctd_max_weight spots against the smallest dominant ",
    "pie slice.")
} else {
  message("'scatterpie' not installed -- skipping SpatialCompositionPlot(). ",
         "Install with install.packages('scatterpie') to run this stage.")
}

message("--- Practice script complete ---")
