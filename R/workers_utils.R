# Internal helper shared across the multi-sample loader functions that
# support a `workers` argument for future.apply-based parallelization across
# samples: CreateRNAObjects, CreateVisiumObjects, CreateAndIntegrateRNA,
# CreateATACObjects, CreateATACObjectsFilter, MakeParseObj.
#
# Each of those functions defaults `workers` directly in its own signature
# to the sample count (e.g. `workers = length(data_dirs)`), rather than a
# fixed `1` -- so by default every sample is read/built in parallel, one
# worker per sample. This helper validates that default (or an explicit
# user-supplied value) against the cores actually available, and errors
# with a concrete suggested value rather than silently letting
# future::plan(multisession, workers = ...) oversubscribe the machine.

# Picks the future backend for the multi-worker plan() calls in
# CreateRNAObjects/CreateVisiumObjects/CreateAndIntegrateRNA/
# CreateATACObjects/CreateATACObjectsFilter/MakeParseObj.
#
# future::multisession workers are separate R *processes* talking over
# local sockets -- every argument going into a worker and every result
# coming back has to be serialized/deserialized. For the payloads these
# functions pass around (Seurat/Signac objects, which carry a @commands
# history slot that grows with every processing step recorded on the
# object) that serialization cost can dominate and even make more workers
# *slower* than fewer -- confirmed empirically: an 8-worker multisession
# run on 8 samples ran for two weeks without finishing, while the same
# work with future::multicore (forked workers, sharing memory via
# copy-on-write, no serialization needed) finished in under 9 minutes.
#
# future::supportsMulticore() is the right gate rather than a bare
# Sys.info()[["sysname"]] != "Windows" check: fork-based workers aren't
# available on Windows at all, and future itself disables multicore
# inside RStudio by default (forking a process that's running RStudio's
# own GUI event loop is a known crash risk) -- supportsMulticore()
# already encodes both exceptions, so this helper just defers to it and
# falls back to multisession whenever it returns FALSE.
#' @keywords internal
#' @noRd
.future_backend <- function(verbose = TRUE) {
  use_multicore <- isTRUE(future::supportsMulticore())
  if (isTRUE(verbose)) {
    # Cheap and easy to check after the fact via system.time()/elapsed-time
    # comparisons (as happened while diagnosing the two-week multisession
    # run above) -- but there's no reason a caller should ever have to
    # infer which backend got picked from timing alone when this can just
    # say so up front.
    message(sprintf("Using future backend: %s%s",
                    if (use_multicore) "multicore" else "multisession",
                    if (!use_multicore) " (multicore unavailable: Windows, or running inside RStudio)" else ""))
  }
  if (use_multicore) future::multicore else future::multisession
}

# Swaps in the future plan .future_backend() picks for `workers` parallel
# workers, returning a cleanup closure the caller must register itself via
# `on.exit(cleanup(), add = TRUE)` -- on.exit() is scoped to the frame that
# calls it, so registering it *inside* this helper would restore the old
# plan the moment this helper returns, not when the caller (which still has
# its whole per-sample loop left to run) eventually returns.
#
# Shared by the six workers-taking loader functions (CreateRNAObjects,
# CreateATACObjects, CreateATACObjectsFilter, CreateVisiumObjects,
# CreateAndIntegrateRNA, MakeParseObj) -- previously an identical
# requireNamespace()+plan()+on.exit() stanza, copy-pasted six times (the
# multicore-vs-multisession fix above needed six manual edits as a result of
# that duplication). CreateRNAObjects additionally clamps BLAS threads
# around its own call to this helper -- that part isn't shared, since it's
# the only one of the six with an empirically-confirmed need for it (see its
# own comment) -- so this helper only covers the plan swap itself, not any
# BLAS handling.
#
# Only call this when workers > 1.
#' @keywords internal
#' @noRd
.setup_future_plan <- function(workers) {
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Package 'future.apply' is required for workers > 1. ",
         "install.packages('future.apply')")
  }
  old_plan <- future::plan(.future_backend(), workers = workers)
  function() future::plan(old_plan)
}

#' @keywords internal
#' @noRd
.resolve_workers <- function(workers, n_samples, was_default) {
  # parallel::detectCores() can return NA on minimal/HPC container shells
  # missing the `wc`/`nproc` tools it shells out to (the same edge case
  # RunRCTD.R guards around spacexr's internal detectCores() calls) -- can't
  # safely validate `workers` against an unknown core count, so skip the
  # check rather than either erroring or silently trusting an unverifiable
  # `workers`.
  n_cores <- suppressWarnings(parallel::detectCores())
  if (is.na(n_cores)) {
    return(workers)
  }

  if (workers > n_cores) {
    # `was_default` (the caller's `missing(workers)`) distinguishes "you
    # explicitly asked for more workers than there are cores" from "this
    # defaulted to the sample count and that happens to exceed cores" --
    # the suggested fix is the same either way, but the reason shouldn't
    # claim a default was used when the caller actually set it themselves.
    reason <- if (isTRUE(was_default)) {
      sprintf(" (defaulted to the sample count, %d, since `workers` wasn't set explicitly)",
              n_samples)
    } else {
      ""
    }
    stop(sprintf(
      "workers = %d exceeds the %d core(s) available on this machine%s. Pass workers = %d (or lower) explicitly.",
      workers, n_cores, reason, n_cores))
  }

  workers
}
