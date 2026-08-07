#!/usr/bin/env Rscript
# Worker for dev/benchmark_bpcells.R -- runs ONE scenario (one data type x
# one on_disk setting) in its own fresh R process and writes timing/memory
# checkpoints to disk. Not normally run by hand, but it's a plain script and
# can be -- useful for debugging a single scenario in isolation:
#
#   config <- list(data_type = "rna", on_disk = TRUE,
#                  bpcells_dir = "bpcells_debug", data_dirs = c("data/s1", "data/s2"),
#                  hd_bin_size = "008um", run_downstream = TRUE,
#                  include_doublet_finder = FALSE, n_variable_features = 2000,
#                  out_rds = "debug_result.rds")
#   saveRDS(config, "debug_config.rds")
#   # Rscript dev/benchmark_bpcells_worker.R debug_config.rds

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript benchmark_bpcells_worker.R <config.rds>")
config <- readRDS(args[[1]])

suppressPackageStartupMessages({
  library(SingleCellTools)
  library(Seurat)
})

# CreateRNAObjects()/CreateVisiumObjects() both print() a QC ggplot as a side
# effect, which forces R to open a real graphics device. Under non-interactive
# Rscript with no display attached, that lazily opens the platform default
# (quartz on macOS), which can fail to close/flush cleanly at process exit --
# surfacing as "Error: write failed" / non-zero exit status AFTER all the
# actual work (including saveRDS() of this scenario's results) already
# succeeded. Open a null PDF device up front so that plotting has somewhere
# real and simple to go instead.
grDevices::pdf(NULL)

# Process resident memory in MB. /proc/self/status (Linux) reports actual
# physical memory the process holds -- including whatever BPCells' compiled
# backend allocates outside R's own allocator, which gc()-based accounting
# alone would miss. Falls back to summing gc()'s reported "(Mb)" columns
# (R-managed memory only) on platforms without /proc (e.g. macOS).
get_rss_mb <- function() {
  status_file <- "/proc/self/status"
  if (file.exists(status_file)) {
    lines <- readLines(status_file, warn = FALSE)
    vmrss_line <- grep("^VmRSS:", lines, value = TRUE)
    if (length(vmrss_line) == 1) {
      kb <- as.numeric(gsub("[^0-9]", "", vmrss_line))
      if (!is.na(kb)) return(kb / 1024)
    }
  }
  gc_info <- gc(FALSE)
  sum(gc_info[, 2])
}

checkpoints <- list()
t0     <- Sys.time()
t_prev <- t0

record <- function(stage) {
  now <- Sys.time()
  checkpoints[[stage]] <<- list(
    cumulative_sec = as.numeric(difftime(now, t0, units = "secs")),
    stage_sec      = as.numeric(difftime(now, t_prev, units = "secs")),
    rss_mb         = get_rss_mb()
  )
  t_prev <<- now
}

message(sprintf("[worker] %s / on_disk = %s / %d sample(s)",
                config$data_type, config$on_disk, length(config$data_dirs)))

# Clean up a stale on-disk matrix from a previous run of this same scenario
# so re-running the benchmark doesn't hit ConvertToBPCells()'s overwrite
# guard.
if (isTRUE(config$on_disk) && dir.exists(config$bpcells_dir)) {
  unlink(config$bpcells_dir, recursive = TRUE)
}

if (identical(config$data_type, "rna")) {
  objs <- CreateRNAObjects(
    data_dirs           = config$data_dirs,
    # Doublet calling is slow and adds a lot of run-to-run variance that has
    # nothing to do with on_disk's actual cost -- off by default so the
    # comparison isolates read + (optional) BPCells conversion time. Set
    # include_doublet_finder = TRUE in the orchestrator's CONFIG block to
    # benchmark the fully realistic pipeline instead.
    run_doublet_finder  = isTRUE(config$include_doublet_finder),
    on_disk             = config$on_disk,
    bpcells_dir         = config$bpcells_dir
  )
} else if (identical(config$data_type, "visium")) {
  objs <- CreateVisiumObjects(
    data_dirs   = config$data_dirs,
    hd_bin_size = config$hd_bin_size,
    on_disk     = config$on_disk,
    bpcells_dir = config$bpcells_dir
  )
} else {
  stop("Unknown data_type: ", config$data_type)
}
record("construction")
# NB: on_disk = TRUE converts as the LAST step of construction (after QC/
# doublet calling), so "construction" will always take a bit longer under
# on_disk = TRUE than on_disk = FALSE -- that's the expected, real cost of
# writing to disk. The payoff shows up in rss_mb here and at every stage
# after, not in this stage's timing.

# Merge all samples into one object before downstream analysis (if there's
# more than one). This is deliberately a plain merge() -- the same call
# MergeSeurat() uses internally -- not the full MergeSeurat() pipeline,
# which bundles merging together with SCTransform/Harmony/clustering/UMAP
# and writes RDS/PDF/CSV files as a side effect. That would conflate merge
# cost with a lot of unrelated cost and isn't what on_disk actually affects.
# Merging is the interesting stage for on_disk specifically: merging
# on-disk BPCells-backed layers should stay lazy/cheap, while merging
# in-memory matrices concatenates them in RAM.
if (length(objs) > 1) {
  obj <- suppressWarnings(merge(objs[[1]], objs[-1], add.cell.ids = names(objs)))
  record("merge")
} else {
  obj <- objs[[1]]
}

if (isTRUE(config$run_downstream)) {
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  record("NormalizeData")

  obj <- Seurat::FindVariableFeatures(obj, nfeatures = config$n_variable_features, verbose = FALSE)
  record("FindVariableFeatures")

  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  record("ScaleData")

  obj <- Seurat::RunPCA(obj, npcs = 30, verbose = FALSE)
  record("RunPCA")
}

disk_mb <- NA_real_
if (isTRUE(config$on_disk) && dir.exists(config$bpcells_dir)) {
  fs <- list.files(config$bpcells_dir, recursive = TRUE, full.names = TRUE)
  if (length(fs) > 0) disk_mb <- sum(file.size(fs)) / 1024^2
}

saveRDS(
  list(
    data_type   = config$data_type,
    on_disk     = config$on_disk,
    n_samples   = length(config$data_dirs),
    checkpoints = checkpoints,
    disk_mb     = disk_mb
  ),
  config$out_rds
)

message("[worker] done")
