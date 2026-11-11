#' Remove ambient RNA contamination from a Seurat object's counts
#'
#' Wraps \code{decontX} (\code{celda} package) or \code{SoupX} to estimate
#' and subtract ambient/background RNA contamination, matching this
#' package's existing "known upstream tool, one call" wrappers
#' (\code{\link{RunLIANA}}, \code{\link{RunCellChat}},
#' \code{\link{RunRCTD}}, \code{\link{RunBanksyWrapper}}). Particularly
#' relevant for tissues where one dominant, fragile cell type contributes
#' outsized ambient signal -- hepatocytes in liver are a textbook case.
#'
#' \strong{Which method.} \code{method = "decontx"} (default) estimates
#' contamination directly from the filtered counts matrix already in
#' \code{obj} -- no extra input needed, at the cost of being a fully
#' data-driven estimate. \code{method = "soupx"} additionally needs the
#' raw, \emph{unfiltered} droplet matrix (every barcode SoupX should treat
#' as potential background, not just the cells Cell Ranger/your pipeline
#' kept) via \code{raw_counts}. Unlike decontX, SoupX's own automatic
#' contamination estimation (\code{SoupX::autoEstCont()}, used whenever
#' \code{contamination_fraction} is \code{NULL}) \emph{requires} cluster
#' labels -- confirmed by its own runtime error ("Clustering information
#' must be supplied, run setClusters first.") when none are set, not just
#' a "less accurate" fallback -- so \code{method = "soupx"} needs either
#' \code{clusters_col} or an explicit \code{contamination_fraction}; this
#' is checked up front with a clear error rather than left to fail deep
#' inside SoupX partway through a (possibly long, possibly parallel)
#' run.
#'
#' \strong{This overwrites the assay's counts layer by default.} Pass
#' \code{new_assay} to instead write the corrected counts into a new assay
#' and leave the original untouched.
#'
#' \strong{Multiple samples at once.} \code{obj} may instead be a
#' (optionally named) list of Seurat objects, same as
#' \code{\link{RunRCTD}}/\code{\link{RunCARD}}/\code{\link{RunSPOTlight}}
#' -- each is processed independently and the result comes back in the
#' same shape as \code{obj} (a bare object in, a bare object out; a list
#' in, a list out). This is the one lever for parallelizing
#' \code{RemoveAmbientRNA} itself: there's no per-object loop to
#' parallelize internally (each call to decontX/SoupX is already one shot
#' against a whole counts matrix), so parallelism here means farming
#' independent *objects* out to \code{workers} processes, not splitting
#' the work within a single object. Errors propagate immediately (no
#' per-object \code{tryCatch}), matching this package's other
#' \code{workers}-taking loaders.
#'
#' @param obj A Seurat object, or a (optionally named) list of them (see
#'   Details).
#' @param method \code{"decontx"} (default) or \code{"soupx"}. See Details.
#' @param assay Assay to correct. Default \code{DefaultAssay(obj)} (per
#'   object, when \code{obj} is a list).
#' @param clusters_col Metadata column of per-cell cluster/cell-type labels.
#'   Passed to decontX's \code{z} argument, or to SoupX's
#'   \code{setClusters()}. \code{NULL} (default): decontX clusters
#'   internally. For \code{method = "soupx"}, \code{NULL} is only valid
#'   together with an explicit \code{contamination_fraction} -- SoupX's
#'   automatic estimation has no unclustered mode (see Details). Applied
#'   by the same column name to every object when \code{obj} is a list.
#' @param raw_counts \code{method = "soupx"} only (required there): the
#'   raw, unfiltered droplet x gene counts matrix, or a path to a directory
#'   \code{Seurat::Read10X()} can read. If that path resolves to multiple
#'   feature-type matrices (e.g. a \code{cellranger multi} raw matrix
#'   bundling Gene Expression with Multiplexing Capture/Antibody Capture
#'   rows), the \code{"Gene Expression"} matrix is used automatically.
#'   Whatever the source, it's then subset (and reordered) to the same
#'   genes as \code{obj}'s filtered counts before being handed to
#'   \code{SoupX::SoupChannel()}, which requires the two to match exactly
#'   and does not align them itself -- errors clearly if any gene present
#'   in the filtered counts is missing from \code{raw_counts} entirely.
#'   When \code{obj} is a list, this can be a single value (used for every
#'   object -- unusual, since \code{soupx} normally needs one raw matrix
#'   per sample) or a list the same length as \code{obj}, one entry per
#'   object -- matched by name when both \code{obj} and this list are
#'   fully named, otherwise matched positionally.
#' @param contamination_fraction \code{method = "soupx"} only: skip
#'   \code{autoEstCont()} and use this fixed contamination fraction for
#'   every cell instead. \code{NULL} (default) estimates it automatically
#'   -- which requires \code{clusters_col} to also be set (see Details);
#'   supplying \code{contamination_fraction} directly is the way to use
#'   \code{method = "soupx"} \emph{without} cluster labels. Accepts a
#'   single value or a per-object list, exactly like \code{raw_counts}
#'   above, when \code{obj} is a list.
#' @param new_assay If set, write corrected counts to a new assay of this
#'   name instead of overwriting \code{assay}'s counts layer in place.
#'   \code{NULL} (default) overwrites in place.
#' @param round_to_integer If \code{TRUE} (default), round corrected counts
#'   to the nearest integer -- both tools return non-integer "corrected"
#'   values, and most downstream steps (DoubletFinder, DE) expect integer
#'   counts.
#' @param verbose Message progress and the mean estimated contamination
#'   fraction. Default \code{TRUE}.
#' @param workers Number of parallel workers to use (via
#'   \code{future.apply}) when \code{obj} is a list -- the one step this
#'   parallelizes (see Details); ignored for a single Seurat object.
#'   Default \code{1} (sequential), matching
#'   \code{\link{NeighborhoodEnrichment}}'s \code{workers} -- parallelism
#'   across objects is opt-in rather than automatic, since decontX/SoupX are
#'   themselves often memory/CPU-heavy, and a per-object default would risk
#'   oversubscribing a constrained machine (CI runners included) with no way
#'   to opt back out. Errors up front if an explicit value exceeds
#'   \code{parallel::detectCores()}, naming the number of cores actually
#'   available. \code{workers > 1} spins up that many parallel
#'   workers via \code{future::plan()} -- forked processes
#'   (\code{future::multicore}) on Unix-likes outside RStudio, or
#'   background R sessions (\code{future::multisession}) on Windows / in
#'   RStudio, where forking isn't available -- torn down again before this
#'   function returns.
#' @return \code{obj} (or a list of them, matching the shape of \code{obj}
#'   -- see Details), with corrected counts (in place or in
#'   \code{new_assay}) and a new \code{<method>_contamination} metadata
#'   column recording each cell's estimated contamination fraction.
#' @examples
#' \dontrun{
#' # decontX -- no extra input needed
#' obj <- RemoveAmbientRNA(obj, clusters_col = "seurat_clusters")
#'
#' # SoupX -- needs the raw, unfiltered matrix too
#' obj <- RemoveAmbientRNA(obj, method = "soupx",
#'                         raw_counts = "raw_feature_bc_matrix/",
#'                         clusters_col = "seurat_clusters")
#'
#' # A list of objects, decontX, in parallel (one worker per object)
#' obj_list <- RemoveAmbientRNA(obj_list, clusters_col = "seurat_clusters",
#'                             workers = length(obj_list))
#'
#' # A list of objects, SoupX, one raw matrix path per sample (named to
#' # match obj_list's own names)
#' obj_list <- RemoveAmbientRNA(
#'   obj_list, method = "soupx",
#'   raw_counts = list(sample1 = "data/sample1/raw_feature_bc_matrix/",
#'                     sample2 = "data/sample2/raw_feature_bc_matrix/"),
#'   clusters_col = "seurat_clusters")
#' }
#' @importFrom Seurat DefaultAssay GetAssayData SetAssayData CreateAssayObject Read10X
#' @export
RemoveAmbientRNA <- function(obj,
                             method                  = c("decontx", "soupx"),
                             assay                   = NULL,
                             clusters_col            = NULL,
                             raw_counts              = NULL,
                             contamination_fraction  = NULL,
                             new_assay               = NULL,
                             round_to_integer        = TRUE,
                             verbose                 = TRUE,
                             workers                 = 1) {

  method <- match.arg(method)

  parsed     <- .as_seurat_list(obj)
  objs       <- parsed$objs
  was_single <- parsed$was_single
  orig_names <- names(objs)
  n          <- length(objs)

  if (!was_single) {
    workers <- .resolve_workers(workers, n_samples = n, was_default = missing(workers))
  }

  # A single (non-list) raw_counts/contamination_fraction is recycled to
  # every object; a *list* must match objs one-to-one -- silently
  # recycling a list element-by-element (R's usual recycling rule) would
  # risk quietly pairing the wrong sample's raw matrix/fraction with
  # another sample, which for method = "soupx" (an incorrect background
  # estimate baked silently into "corrected" counts) is exactly the kind
  # of mistake this should catch rather than allow.
  .expand_per_sample <- function(x, arg_name) {
    if (is.null(x) || !is.list(x)) {
      return(rep(list(x), n))
    }
    if (length(x) != n) {
      stop(sprintf(
        "`%s` is a list of length %d but `obj` has %d element%s -- these must match one-to-one.",
        arg_name, length(x), n, if (n == 1) "" else "s"))
    }
    objs_named <- !is.null(orig_names) && all(nzchar(orig_names))
    x_named    <- !is.null(names(x))   && all(nzchar(names(x)))
    if (objs_named && x_named) {
      if (!setequal(names(x), orig_names)) {
        stop(sprintf("`%s`'s names don't match `obj`'s names.", arg_name))
      }
      x <- x[orig_names]
    }
    x
  }

  raw_counts_list <- .expand_per_sample(raw_counts, "raw_counts")
  contamination_fraction_list <- .expand_per_sample(contamination_fraction, "contamination_fraction")

  if (isTRUE(verbose) && !was_single) {
    message(sprintf(
      '--- Running RemoveAmbientRNA (method = "%s") on %d objects%s ---',
      method, n, if (workers > 1) sprintf(', %d parallel workers', workers) else ''))
  }

  # Every element is already guaranteed to be a Seurat object here --
  # .as_seurat_list() validated that above.
  .process_one <- function(i) {
    .remove_ambient_rna_one(objs[[i]], method = method, assay = assay,
                            clusters_col = clusters_col,
                            raw_counts = raw_counts_list[[i]],
                            contamination_fraction = contamination_fraction_list[[i]],
                            new_assay = new_assay,
                            round_to_integer = round_to_integer,
                            verbose = verbose)
  }

  if (!was_single && workers > 1) {
    cleanup <- .setup_future_plan(workers)
    on.exit(cleanup(), add = TRUE)
    objs <- future.apply::future_lapply(seq_len(n), .process_one, future.seed = TRUE)
  } else {
    objs <- lapply(seq_len(n), .process_one)
  }
  names(objs) <- orig_names

  if (was_single) return(objs[[1]])
  objs
}

# Does the actual decontX/SoupX work for exactly one Seurat object --
# unchanged from RemoveAmbientRNA's own body before it gained list/workers
# support above, just split out so the exported function can call it once
# per list element (sequentially, or in parallel via
# future.apply::future_lapply).
#' @keywords internal
#' @noRd
.remove_ambient_rna_one <- function(obj, method, assay, clusters_col,
                                    raw_counts, contamination_fraction,
                                    new_assay, round_to_integer, verbose) {

  if (!is.null(clusters_col) && !clusters_col %in% colnames(obj@meta.data)) {
    stop("`clusters_col` ('", clusters_col, "') not found in obj@meta.data.")
  }
  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay
  counts <- Seurat::GetAssayData(obj, assay = a, layer = "counts")

  if (method == "decontx") {
    if (!requireNamespace("celda", quietly = TRUE)) {
      stop("'celda' is required for method = 'decontx'. Install with ",
           "BiocManager::install('celda').")
    }
    z <- if (!is.null(clusters_col)) as.character(obj@meta.data[[clusters_col]]) else NULL
    if (isTRUE(verbose)) {
      message(sprintf("--- Running decontX (assay = '%s'%s) ---", a,
                      if (!is.null(z)) ", using provided clusters" else ""))
    }
    res <- celda::decontX(as.matrix(counts), z = z)
    corrected <- res$decontXcounts
    contamination <- stats::setNames(res$contamination, colnames(counts))

  } else {
    if (!requireNamespace("SoupX", quietly = TRUE)) {
      stop("'SoupX' is required for method = 'soupx'. Install with ",
           "install.packages('SoupX').")
    }
    if (is.null(clusters_col) && is.null(contamination_fraction)) {
      # SoupX::autoEstCont() (run below whenever contamination_fraction is
      # NULL) has no unclustered fallback -- confirmed by its own runtime
      # error ("Clustering information must be supplied, run setClusters
      # first."), not just reduced accuracy, despite what this function's
      # own docs used to say before that was confirmed against real SoupX
      # behavior. Caught here, before SoupChannel()/Read10X() even run, so
      # a list-mode call with workers > 1 fails immediately for every
      # object instead of spinning up a whole parallel plan only to have
      # every worker fail deep inside autoEstCont() with a much less
      # legible error.
      stop("method = 'soupx' needs either `clusters_col` (for SoupX's ",
           "per-cluster automatic contamination estimation) or an ",
           "explicit `contamination_fraction` -- without clusters, ",
           "SoupX::autoEstCont() errors internally with 'Clustering ",
           "information must be supplied, run setClusters first.'; there ",
           "is no unclustered automatic-estimation fallback.")
    }
    if (is.null(raw_counts)) {
      stop("`raw_counts` (the raw, unfiltered droplet matrix, or a path ",
           "Seurat::Read10X() can read) is required for method = 'soupx'.")
    }
    tod <- if (is.character(raw_counts)) {
      Seurat::Read10X(data.dir = raw_counts)
    } else {
      raw_counts
    }

    # Seurat::Read10X() returns a *list* of matrices, one per feature type,
    # whenever the directory's features.tsv.gz lists more than one type --
    # common for cellranger multi's combined raw_feature_bc_matrix, which
    # often bundles Gene Expression together with Multiplexing Capture
    # (CMO) and/or Antibody Capture rows even when the per-sample filtered
    # matrix (what `counts`/toc came from) only ever had Gene Expression.
    # Passing that whole list straight to SoupX::SoupChannel() either
    # errors opaquely or (if nrow() on the list happens not to error)
    # compares the wrong thing entirely, so pull out "Gene Expression"
    # specifically here, where the problem can be named clearly.
    if (is.list(tod)) {
      if (!"Gene Expression" %in% names(tod)) {
        stop("`raw_counts` resolved to multiple feature-type matrices (",
             paste(names(tod), collapse = ", "), ") but none is named ",
             "'Gene Expression' -- pass the specific matrix instead, e.g. ",
             "Seurat::Read10X(data.dir = ...)[['Gene Expression']].")
      }
      tod <- tod[["Gene Expression"]]
    }

    # SoupX::SoupChannel() requires tod and toc to have identical genes in
    # identical order -- it does not align them by name itself. tod (the
    # raw/unfiltered matrix) commonly has *more* genes than toc (`counts`,
    # already gene-filtered by e.g. CreateSeuratObject(min.cells = ...), or
    # missing probes CellRanger multi drops from the per-sample filtered
    # matrix but keeps in the pooled raw one) -- confirmed a real-world
    # mismatch, not just a theoretical one (see
    # github.com/constantAmateur/SoupX/issues/148). Subsetting tod to toc's
    # genes, in toc's order, is the fix recommended there; done by name
    # (not a positional truncation) so a real gene-order difference doesn't
    # silently miscorrect every count instead of erroring.
    missing_genes <- setdiff(rownames(counts), rownames(tod))
    if (length(missing_genes) > 0) {
      stop(sprintf(
        paste0("%d of %d gene(s) in the filtered counts (assay = '%s') are ",
               "not present in `raw_counts` -- cannot align the raw and ",
               "filtered matrices for SoupX. First missing: %s."),
        length(missing_genes), nrow(counts), a,
        paste(utils::head(missing_genes, 5), collapse = ", ")))
    }
    tod <- tod[rownames(counts), , drop = FALSE]

    if (isTRUE(verbose)) message(sprintf("--- Running SoupX (assay = '%s') ---", a))
    sc <- SoupX::SoupChannel(tod, counts)
    if (!is.null(clusters_col)) {
      sc <- SoupX::setClusters(
        sc, stats::setNames(as.character(obj@meta.data[[clusters_col]]), colnames(obj)))
    }
    sc <- if (!is.null(contamination_fraction)) {
      SoupX::setContaminationFraction(sc, contamination_fraction)
    } else {
      SoupX::autoEstCont(sc, doPlot = FALSE, verbose = isTRUE(verbose))
    }
    corrected <- SoupX::adjustCounts(sc)
    contamination <- stats::setNames(sc$metaData$rho, rownames(sc$metaData))
    contamination <- contamination[colnames(counts)]
  }

  if (isTRUE(round_to_integer)) corrected <- round(corrected)

  if (isTRUE(verbose)) {
    message(sprintf("  Mean estimated contamination fraction: %.3f (range %.3f-%.3f)",
                    mean(contamination, na.rm = TRUE),
                    min(contamination, na.rm = TRUE),
                    max(contamination, na.rm = TRUE)))
  }

  if (is.null(new_assay)) {
    obj <- Seurat::SetAssayData(obj, assay = a, layer = "counts", new.data = corrected)
  } else {
    obj[[new_assay]] <- Seurat::CreateAssayObject(counts = corrected)
  }
  obj[[paste0(method, "_contamination")]] <- contamination[colnames(obj)]

  obj
}
