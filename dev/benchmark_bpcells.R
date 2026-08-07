#!/usr/bin/env Rscript
# =============================================================================
# Benchmark: on_disk = FALSE vs on_disk = TRUE (BPCells) for
# CreateRNAObjects() and CreateVisiumObjects().
#
# WHY SEPARATE R PROCESSES PER SCENARIO
#   Comparing memory between on_disk settings only means something if each
#   scenario starts from a clean process -- running multiple scenarios in one
#   R session lets earlier allocations (freed but not necessarily returned to
#   the OS by R's allocator) inflate later measurements. This script runs
#   each scenario (data type x on_disk setting) as its own `Rscript`
#   subprocess via benchmark_bpcells_worker.R, then compares results
#   afterward. That worker file needs to sit alongside this one.
#
# WHAT GETS MEASURED
#   Per scenario: elapsed time and process RSS (resident memory) checkpointed
#   after object construction, after merging all samples into one object
#   (plain merge() -- see benchmark_bpcells_worker.R for why not
#   MergeSeurat()) if more than one sample was provided, and, optionally,
#   after each of NormalizeData/FindVariableFeatures/ScaleData/RunPCA on the
#   merged object. For on_disk = TRUE scenarios, also the resulting on-disk
#   footprint in MB.
#
#   The merge stage is where on_disk should matter most with multiple
#   samples: merging on-disk BPCells-backed layers should stay lazy/cheap,
#   while merging in-memory matrices concatenates them in RAM.
#
#   Expect construction itself to take slightly LONGER under on_disk = TRUE
#   (the BPCells conversion is a real cost, applied as the last construction
#   step) -- the payoff is in the RSS numbers at that stage and every stage
#   after, not construction time itself. See
#   SingleCellTools_vignette_bpcells.md for the full explanation of that
#   design choice.
#
# USAGE
#   1. Edit the CONFIG block below -- at minimum, set rna_dirs and/or
#      visium_dirs to real data directories.
#   2. Rscript dev/benchmark_bpcells.R
#      (or source() interactively, from an R session started in the repo
#      root or with `outdir`/paths set to full paths)
#
# REQUIREMENTS
#   SingleCellTools installed (obviously), and BPCells installed if you want
#   the on_disk = TRUE scenarios to actually run (they'll fail cleanly with a
#   warning and get skipped otherwise -- the on_disk = FALSE scenarios still
#   run fine without it).
#
# OUTPUT (written to `outdir`)
#   - bpcells_benchmark_results.csv -- one row per (scenario, stage), with
#     cumulative and per-stage elapsed time, process RSS, and (for
#     on_disk = TRUE) the on-disk footprint in MB.
#   - bpcells_benchmark_comparison.csv -- pivoted head-to-head: on_disk =
#     TRUE vs FALSE side by side per (data_type, stage), with a time ratio
#     and RSS MB/percent saved -- the direct "was it faster, was it lighter"
#     answer, printed to the console too.
#   - bpcells_benchmark_plot.pdf -- time and RSS comparison bar charts, if
#     ggplot2 is available (it's a hard SingleCellTools dependency, so it
#     will be).
# =============================================================================

## ---- CONFIG -- edit before running -----------------------------------------
rna_dirs    <- list.dirs('~/Downloads/data/rna', recursive = F)   # e.g. list.dirs("data/rna", recursive = FALSE)
visium_dirs <- list.dirs('~/Downloads/data/visium', recursive = F)   # e.g. list.dirs("data/visium_hd", recursive = FALSE)
hd_bin_size <- "008um"        # only matters if visium_dirs are HD samples: "002um"/"008um"/"016um"

bpcells_root <- file.path(tempdir(), "bpcells_benchmark")

run_downstream          <- TRUE   # also time Normalize/FindVariableFeatures/ScaleData/RunPCA
include_doublet_finder  <- FALSE  # RNA only; slow + noisy, off by default -- see worker script comments
n_variable_features     <- 2000

outdir <- getwd()
## -----------------------------------------------------------------------------

if (length(rna_dirs) == 0 && length(visium_dirs) == 0) {
  stop("Set rna_dirs and/or visium_dirs in the CONFIG block before running ",
       "(both are empty character(0) by default).")
}

# Standard Rscript idiom for finding the currently-running script's own
# directory, so this can locate benchmark_bpcells_worker.R next to it
# regardless of the working directory it's launched from. Falls back to
# "dev" for source()'d/interactive use, where --file= isn't set.
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  "dev"
}
script_dir  <- get_script_dir()
worker_path <- file.path(script_dir, "benchmark_bpcells_worker.R")
if (!file.exists(worker_path)) {
  stop("Could not find benchmark_bpcells_worker.R next to this script ",
      "(looked in '", script_dir, "'). Make sure both files stay together, ",
      "or edit `worker_path` above manually.")
}

run_scenario <- function(data_type, on_disk, data_dirs) {
  config <- list(
    data_type              = data_type,
    on_disk                = on_disk,
    bpcells_dir            = file.path(bpcells_root, paste0(data_type, "_", if (on_disk) "ondisk" else "inmemory")),
    data_dirs              = data_dirs,
    hd_bin_size            = hd_bin_size,
    run_downstream         = run_downstream,
    include_doublet_finder = include_doublet_finder,
    n_variable_features    = n_variable_features,
    out_rds                = tempfile(fileext = ".rds")
  )
  config_rds <- tempfile(fileext = ".rds")
  saveRDS(config, config_rds)

  message(sprintf("=== %s / on_disk = %s (%d sample(s)) ===", data_type, on_disk, length(data_dirs)))
  status <- system2("Rscript", args = c(shQuote(worker_path), shQuote(config_rds)))

  if (!identical(status, 0L) || !file.exists(config$out_rds)) {
    warning(sprintf(
      "Scenario %s / on_disk = %s failed (subprocess exit status %s). ", data_type, on_disk, status),
      "If on_disk = TRUE, confirm BPCells is installed for the R that `Rscript` resolves to ",
      "(can differ from your interactive session's R if you have multiple R installs). Skipping.",
      call. = FALSE)
    return(NULL)
  }
  readRDS(config$out_rds)
}

results <- list()
if (length(rna_dirs) > 0) {
  results[["rna_inmemory"]] <- run_scenario("rna", FALSE, rna_dirs)
  results[["rna_ondisk"]]   <- run_scenario("rna", TRUE,  rna_dirs)
}
if (length(visium_dirs) > 0) {
  results[["visium_inmemory"]] <- run_scenario("visium", FALSE, visium_dirs)
  results[["visium_ondisk"]]   <- run_scenario("visium", TRUE,  visium_dirs)
}
results <- Filter(Negate(is.null), results)

if (length(results) == 0) {
  stop("Every scenario failed -- see warnings above. Nothing to summarize.")
}

## ---- Summarize --------------------------------------------------------------
summary_rows <- lapply(names(results), function(nm) {
  r <- results[[nm]]
  do.call(rbind, lapply(names(r$checkpoints), function(stage) {
    cp <- r$checkpoints[[stage]]
    data.frame(
      scenario       = nm,
      data_type      = r$data_type,
      on_disk        = r$on_disk,
      n_samples      = r$n_samples,
      stage          = stage,
      cumulative_sec = cp$cumulative_sec,
      stage_sec      = cp$stage_sec,
      rss_mb         = cp$rss_mb,
      disk_mb        = r$disk_mb,
      stringsAsFactors = FALSE
    )
  }))
})
summary_df <- do.call(rbind, summary_rows)

print(summary_df)
csv_path <- file.path(outdir, "bpcells_benchmark_results.csv")
write.csv(summary_df, csv_path, row.names = FALSE)
message(sprintf("Wrote %s", csv_path))

## ---- Head-to-head comparison (speed AND memory, on_disk vs in-memory) -------
# The long-format summary_df above already has both timing (cumulative_sec/
# stage_sec) and memory (rss_mb) per stage -- this section just pivots that
# into a direct in-memory-vs-on-disk comparison per (data_type, stage), which
# is what you actually want to eyeball: is on_disk faster/slower, and by how
# much memory does it save.
inmem <- summary_df[!summary_df$on_disk, c("data_type", "stage", "cumulative_sec", "rss_mb", "n_samples")]
ondsk <- summary_df[summary_df$on_disk,  c("data_type", "stage", "cumulative_sec", "rss_mb", "disk_mb")]
names(inmem)[names(inmem) %in% c("cumulative_sec", "rss_mb")] <- c("sec_inmemory", "rss_mb_inmemory")
names(ondsk)[names(ondsk) %in% c("cumulative_sec", "rss_mb")] <- c("sec_ondisk", "rss_mb_ondisk")

comparison_df <- merge(inmem, ondsk, by = c("data_type", "stage"))
if (nrow(comparison_df) > 0) {
  comparison_df$time_ratio_ondisk_vs_inmemory <- comparison_df$sec_ondisk / comparison_df$sec_inmemory
  comparison_df$rss_mb_saved                  <- comparison_df$rss_mb_inmemory - comparison_df$rss_mb_ondisk
  comparison_df$rss_pct_saved                 <- 100 * comparison_df$rss_mb_saved / comparison_df$rss_mb_inmemory
  comparison_df <- comparison_df[, c("data_type", "stage", "n_samples",
                                     "sec_inmemory", "sec_ondisk", "time_ratio_ondisk_vs_inmemory",
                                     "rss_mb_inmemory", "rss_mb_ondisk", "rss_mb_saved", "rss_pct_saved",
                                     "disk_mb")]

  cat("\n--- Head-to-head: on_disk vs in-memory (time_ratio > 1 means on_disk was slower; rss_pct_saved > 0 means on_disk used less memory) ---\n")
  print(comparison_df)

  comparison_csv_path <- file.path(outdir, "bpcells_benchmark_comparison.csv")
  write.csv(comparison_df, comparison_csv_path, row.names = FALSE)
  message(sprintf("Wrote %s", comparison_csv_path))
} else {
  message("Skipping head-to-head comparison table: need both on_disk = FALSE and ",
         "on_disk = TRUE results for at least one data_type to compare (you may have ",
         "only run one side, or one side's scenario failed -- see warnings above).")
}

## ---- Optional plot -----------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE) && nrow(summary_df) > 0) {
  summary_df$stage <- factor(summary_df$stage, levels = unique(summary_df$stage))

  p_time <- ggplot2::ggplot(summary_df, ggplot2::aes(stage, cumulative_sec, fill = scenario)) +
    ggplot2::geom_col(position = ggplot2::position_dodge()) +
    ggplot2::labs(title = "Cumulative elapsed time by stage", y = "Seconds", x = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p_mem <- ggplot2::ggplot(summary_df, ggplot2::aes(stage, rss_mb, fill = scenario)) +
    ggplot2::geom_col(position = ggplot2::position_dodge()) +
    ggplot2::labs(title = "Process RSS by stage", y = "RSS (MB)", x = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  pdf_path <- file.path(outdir, "bpcells_benchmark_plot.pdf")
  grDevices::pdf(pdf_path, width = 9, height = 5)
  print(p_time)
  print(p_mem)
  grDevices::dev.off()
  message(sprintf("Wrote %s", pdf_path))
}

message("--- Benchmark complete ---")
