#' Batch-load multiple 10X Genomics Xenium samples
#'
#' Thin multi-sample wrapper around \code{\link{LoadXenium2}}, bringing it
#' in line with this package's other batch loaders
#' (\code{\link{CreateRNAObjects}}, \code{\link{CreateATACObjects}},
#' \code{\link{CreateVisiumObjects}}, \code{\link{CreateAndIntegrateRNA}},
#' \code{\link{MakeParseObj}}) -- all of which take a vector of sample
#' directories and a \code{workers} argument for reading them in parallel,
#' rather than requiring a hand-rolled loop around a single-sample loader.
#'
#' @param data_dirs Character vector of Xenium output directories, one per
#'   sample.
#' @param sample_names Character vector of sample names, same length and
#'   order as \code{data_dirs}. \code{NULL} (default) uses
#'   \code{basename(data_dirs)}. Must be unique.
#' @param outs,type,mols.qv.threshold,microns_lazy Passed straight through
#'   to \code{\link{LoadXenium2}} for every sample -- see its docs.
#' @param on_disk Logical; if \code{TRUE}, move each returned object's
#'   Xenium counts layer to an on-disk BPCells matrix via
#'   \code{\link{LoadXenium2}}'s own \code{on_disk} handling. Default
#'   \code{FALSE}.
#' @param bpcells_dir Base directory for on-disk matrices when
#'   \code{on_disk = TRUE}; each sample gets its own
#'   \code{file.path(bpcells_dir, sample_name)} subdirectory so samples
#'   never collide. \code{NULL} (default) uses \code{\link{LoadXenium2}}'s
#'   own per-sample default.
#' @param workers Number of parallel workers to use (via
#'   \code{future.apply}) for reading each sample -- fully independent
#'   across samples. Defaults to \code{length(data_dirs)} (one worker per
#'   sample); errors up front if that (or an explicit value) exceeds
#'   \code{parallel::detectCores()}, naming the number of cores actually
#'   available. Pass \code{workers = 1} to run sequentially instead.
#'   \code{workers > 1} spins up that many parallel workers via
#'   \code{future::plan()} -- forked processes (\code{future::multicore})
#'   on Unix-likes outside RStudio, or background R sessions
#'   (\code{future::multisession}) on Windows / in RStudio, where forking
#'   isn't available -- restored on exit.
#' @return A named list of Seurat objects, one per sample, named by
#'   \code{sample_names}.
#' @examples
#' \dontrun{
#' xenium_dirs <- list.dirs("xenium_runs", recursive = FALSE)
#' xenium_objs <- CreateXeniumObjects(xenium_dirs, workers = 4)
#' }
#' @export
CreateXeniumObjects <- function(data_dirs,
                                sample_names      = NULL,
                                outs              = c("matrix", "microns"),
                                type              = c("centroids", "segmentations"),
                                mols.qv.threshold = 20,
                                microns_lazy      = FALSE,
                                on_disk           = FALSE,
                                bpcells_dir       = NULL,
                                workers           = length(data_dirs)) {

  if (!is.character(data_dirs) || length(data_dirs) == 0) {
    stop("`data_dirs` must be a non-empty character vector of directories.")
  }
  if (is.null(sample_names)) sample_names <- basename(data_dirs)
  if (length(sample_names) != length(data_dirs)) {
    stop("`sample_names` must be the same length as `data_dirs` (", length(data_dirs),
         "), got ", length(sample_names), ".")
  }
  if (anyDuplicated(sample_names)) {
    stop("`sample_names` must be unique: ",
         paste(unique(sample_names[duplicated(sample_names)]), collapse = ", "))
  }

  n_samples <- length(data_dirs)
  workers <- .resolve_workers(workers, n_samples = n_samples,
                             was_default = missing(workers))

  if (workers > 1) {
    # See workers_utils.R -- shared by all six (now seven) workers-taking loaders.
    cleanup <- .setup_future_plan(workers)
    on.exit(cleanup(), add = TRUE)
  }

  .load_one <- function(data_dir, sample_name) {
    LoadXenium2(
      data_dir          = data_dir,
      sample_name       = sample_name,
      outs              = outs,
      type              = type,
      mols.qv.threshold = mols.qv.threshold,
      microns_lazy      = microns_lazy,
      on_disk           = on_disk,
      bpcells_dir       = if (is.null(bpcells_dir)) NULL
                          else file.path(bpcells_dir, sample_name)
    )
  }

  message(sprintf(
    '--- Loading %d Xenium sample(s)%s ---',
    n_samples, if (workers > 1) sprintf(', %d parallel workers', workers) else ''))

  objs <- if (workers > 1) {
    future.apply::future_mapply(.load_one, data_dirs, sample_names,
                                SIMPLIFY = FALSE, future.seed = TRUE)
  } else {
    mapply(.load_one, data_dirs, sample_names, SIMPLIFY = FALSE)
  }
  names(objs) <- sample_names
  objs
}
