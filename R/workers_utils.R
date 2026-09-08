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
# CreateRNAObjects/CreateATACObjects/CreateATACObjectsFilter.
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
.future_backend <- function() {
  if (isTRUE(future::supportsMulticore())) future::multicore else future::multisession
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
