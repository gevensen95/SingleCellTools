#' Create Seurat Objects from Parse Biosciences Pipeline Output
#'
#' Reads one or more Parse Biosciences sample directories and constructs a
#' Seurat object for each, applying basic cell- and feature-level filtering.
#' When given a single path, a single \code{Seurat} object is returned; when
#' given multiple paths, a named list of \code{Seurat} objects is returned.
#'
#' @param paths Character vector of paths to Parse Biosciences sample
#'   directories. Each path must contain a \code{DGE_filtered/} subdirectory
#'   with \code{count_matrix.mtx}, \code{all_genes.csv}, and
#'   \code{cell_metadata.csv}.
#' @param sample_names Optional character vector of names to assign to the
#'   resulting Seurat objects. Must be the same length as \code{paths}. If
#'   \code{NULL} (default), \code{basename(paths)} is used.
#' @param treatments Optional vector of treatment labels, one per path, in the
#'   same order as \code{paths}. When supplied, a metadata column (named by
#'   \code{treatment_col}) is added to every cell of the corresponding object
#'   with that treatment label. May be a character, factor, or other atomic
#'   vector. Default \code{NULL} (no treatment column added).
#' @param treatment_col Name of the metadata column used to store treatment
#'   labels. Default \code{"Treatment"}.
#' @param mincellfrac Numeric in \code{[0, 1]}. Minimum fraction of cells in
#'   which a feature must be detected to be retained. Passed to
#'   \code{min.cells} in \code{\link[SeuratObject]{CreateAssayObject}} after
#'   multiplication by the number of cells. Default \code{0.0005}.
#' @param mincellfeat Integer. Minimum number of features a cell must express
#'   to be retained. Passed to \code{min.features} in
#'   \code{\link[SeuratObject]{CreateAssayObject}}. Default \code{50}.
#' @param mt_pattern Regex passed to
#'   \code{\link[Seurat]{PercentageFeatureSet}} to identify mitochondrial
#'   genes. Default \code{"^mt-"} (mouse gene symbol convention, e.g.
#'   \code{"mt-Nd1"}, matching the default used by \code{CreateRNAObjects()}
#'   and friends elsewhere in this package); pass \code{"^MT-"} for human
#'   data, or \code{NULL} to skip the calculation entirely.
#' @param mt_col Name of the metadata column in which to store the
#'   mitochondrial percentage when \code{mt_pattern} is non-\code{NULL}.
#'   Default \code{"percent.mt"}.
#' @param rb_pattern Regex passed to
#'   \code{\link[Seurat]{PercentageFeatureSet}} to identify ribosomal-protein
#'   genes. Default \code{"^(Rp[sl]|RP[SL])"} (matches both mouse
#'   \code{"Rps"}/\code{"Rpl"} and human \code{"RPS"}/\code{"RPL"} gene
#'   symbol conventions, same default as elsewhere in this package); pass
#'   \code{NULL} to skip the calculation entirely.
#' @param rb_col Name of the metadata column in which to store the
#'   ribosomal-protein percentage when \code{rb_pattern} is non-\code{NULL}.
#'   Default \code{"percent.rb"}.
#' @param hb_pattern Regex passed to
#'   \code{\link[Seurat]{PercentageFeatureSet}} to identify hemoglobin
#'   genes. Default \code{"^(Hb[^p]|HB[^P])"} (matches mouse
#'   \code{"Hba"}/\code{"Hbb"} and human \code{"HBA"}/\code{"HBB"}
#'   hemoglobin genes while excluding the unrelated \code{"Hbp1"}/
#'   \code{"HBP1"} gene, a well-known false positive for naive
#'   \code{"^Hb"} patterns; same default as elsewhere in this package);
#'   pass \code{NULL} to skip the calculation entirely.
#' @param hb_col Name of the metadata column in which to store the
#'   hemoglobin percentage when \code{hb_pattern} is non-\code{NULL}.
#'   Default \code{"percent.hb"}.
#' @param run_doublet_finder Logical; if TRUE (default), run \code{calldoublet}
#'   on every object and add a \code{doublet_finder} metadata column.
#' @param doublet_normalization Passed to \code{calldoublet}: one of
#'   \code{"LogNormalize"} (default) or \code{"SCT"}.
#' @param doublet_vars_to_regress Passed to \code{calldoublet} as
#'   \code{vars.to.regress}. Default: regresses out \code{mt_col} (i.e.
#'   \code{"percent.mt"}) whenever \code{mt_pattern} is non-\code{NULL}
#'   (true by default), matching \code{CreateRNAObjects()}'s existing
#'   \code{doublet_vars_to_regress = "percent.mt"} default; falls back to
#'   \code{NULL} (no regression) if \code{mt_pattern} is \code{NULL}, since
#'   there'd be no \code{percent.mt} column to regress. \code{percent.rb}/
#'   \code{percent.hb} are NOT regressed by default even though they're
#'   computed by default -- unlike mitochondrial content, ribosomal/
#'   hemoglobin percentage often carries real biological signal, so
#'   regressing it out isn't standard practice. Pass \code{NULL} explicitly
#'   to disable regression entirely, or your own character vector (e.g.
#'   \code{c(mt_col, rb_col)}) to regress something else or something
#'   additional. Requires the corresponding \code{*_pattern} argument to be
#'   non-\code{NULL} for any column you request.
#' @param doublet_cluster_resolution Passed to \code{calldoublet} as
#'   \code{cluster_resolution}. Default \code{0.1}.
#' @param filter_doublets Logical; if TRUE, subset each object to
#'   \code{doublet_finder == "Singlet"} after doublet calling. Default
#'   \code{FALSE} so the doublet labels are preserved for downstream review.
#' @param workers Number of parallel workers to use (via
#'   \code{future.apply}) for reading/creating each sample's Seurat object
#'   and for the per-sample \code{calldoublet} calls -- the two most
#'   expensive steps, and both fully independent across samples. Default
#'   \code{1} runs sequentially exactly as before (with per-sample progress
#'   messages); \code{workers > 1} spins up that many background R sessions
#'   via \code{future::plan(multisession)}, restored on exit. Note each
#'   worker holds its own copy of that sample's data, so peak memory scales
#'   with \code{workers}.
#'
#' @return If \code{length(paths) == 1}, a single \code{Seurat} object.
#'   Otherwise, a named list of \code{Seurat} objects, one per input path,
#'   with names taken from \code{sample_names} (or \code{basename(paths)}).
#'
#' @details
#' For each input path the function:
#' \enumerate{
#'   \item Reads the sparse count matrix from
#'     \code{DGE_filtered/count_matrix.mtx} and transposes it so genes are
#'     rows and cells are columns.
#'   \item Reads gene and cell metadata from \code{all_genes.csv} and
#'     \code{cell_metadata.csv}.
#'   \item Assigns gene names to the matrix rows and cell barcodes
#'     (\code{bc_wells}) to the columns.
#'   \item Builds an RNA assay with \code{\link[SeuratObject]{CreateAssayObject}},
#'     applying the \code{mincellfeat} and \code{mincellfrac} filters.
#'   \item Wraps the assay in a \code{Seurat} object with the cell metadata
#'     attached, optionally tagging every cell with a treatment label.
#'   \item Optionally runs \code{calldoublet} to add a \code{doublet_finder}
#'     metadata column (and, if requested, drops doublets).
#' }
#'
#' @examples
#' \dontrun{
#' # Single sample
#' obj <- MakeParseObj("/data/parse/S1")
#'
#' # Multiple samples — returns a named list
#' sample_paths <- file.path("/data/parse", paste0("S", 1:16))
#' obj_list <- MakeParseObj(sample_paths)
#'
#' # With explicit sample names
#' obj_list <- MakeParseObj(
#'   sample_paths,
#'   sample_names = paste0("Sample", 1:16)
#' )
#'
#' # With treatments tagged in metadata
#' obj_list <- MakeParseObj(
#'   sample_paths,
#'   treatments = rep(c("Vehicle", "DrugA", "DrugB", "DrugC"), each = 4)
#' )
#'
#' # percent.mt/percent.rb/percent.hb are computed by default (mouse
#' # patterns), and doublet calling regresses out percent.mt by default too:
#' obj_list <- MakeParseObj(sample_paths)
#'
#' # Also regress out percent.rb/percent.hb during doublet calling
#' obj_list <- MakeParseObj(
#'   sample_paths,
#'   doublet_vars_to_regress = c("percent.mt", "percent.rb", "percent.hb")
#' )
#'
#' # Compute percent.mt/rb/hb but skip regressing any of them out
#' obj_list <- MakeParseObj(
#'   sample_paths,
#'   doublet_vars_to_regress = NULL
#' )
#'
#' # Human data -- override the mouse-convention defaults
#' obj_list <- MakeParseObj(
#'   sample_paths,
#'   mt_pattern = "^MT-",
#'   rb_pattern = "^RP[SL]",
#'   hb_pattern = "^HB[^P]"
#' )
#'
#' # Skip these calculations entirely
#' obj_list <- MakeParseObj(
#'   sample_paths,
#'   mt_pattern = NULL, rb_pattern = NULL, hb_pattern = NULL
#' )
#' }
#'
#' @export
MakeParseObj <- function(paths,
                         sample_names               = NULL,
                         treatments                 = NULL,
                         treatment_col              = "Treatment",
                         mincellfrac                = 0.0005,
                         mincellfeat                = 50,
                         mt_pattern                 = '^mt-',
                         mt_col                     = "percent.mt",
                         rb_pattern                 = '^(Rp[sl]|RP[SL])',
                         rb_col                     = "percent.rb",
                         hb_pattern                 = '^(Hb[^p]|HB[^P])',
                         hb_col                     = "percent.hb",
                         run_doublet_finder         = TRUE,
                         doublet_normalization      = c("LogNormalize", "SCT"),
                         doublet_vars_to_regress    = NULL,
                         doublet_cluster_resolution = 0.1,
                         filter_doublets            = FALSE,
                         workers                    = 1) {

  doublet_normalization <- match.arg(doublet_normalization)

  if (workers > 1) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Package 'future.apply' is required for workers > 1. ",
           "install.packages('future.apply')")
    }
    old_plan <- future::plan(future::multisession, workers = workers)
    on.exit(future::plan(old_plan), add = TRUE)
  }

  # `doublet_vars_to_regress` defaults to regressing out `mt_col` whenever
  # `mt_pattern` is non-NULL (i.e. that column actually gets computed),
  # matching CreateRNAObjects()'s existing `doublet_vars_to_regress =
  # "percent.mt"` default. Only kicks in when the caller didn't pass this
  # argument at all, so an explicit `NULL` still means "no regression" and
  # an explicit vector still passes straight through untouched.
  if (missing(doublet_vars_to_regress)) {
    doublet_vars_to_regress <- if (!is.null(mt_pattern)) mt_col else NULL
  }

  # ---- Argument checks -----------------------------------------------------
  if (!is.character(paths) || length(paths) < 1) {
    stop("`paths` must be a non-empty character vector.")
  }
  missing_dirs <- !dir.exists(paths)
  if (any(missing_dirs)) {
    stop("The following path(s) do not exist:\n  ",
         paste(paths[missing_dirs], collapse = "\n  "))
  }

  if (is.null(sample_names)) {
    sample_names <- basename(paths)
  } else if (length(sample_names) != length(paths)) {
    stop("`sample_names` must be the same length as `paths`.")
  }

  if (!is.null(treatments)) {
    if (!is.atomic(treatments)) {
      stop("`treatments` must be an atomic vector (character, factor, etc.).")
    }
    if (length(treatments) != length(paths)) {
      stop("`treatments` must be the same length as `paths` ",
           "(got ", length(treatments), " treatments for ",
           length(paths), " paths).")
    }
    if (!is.character(treatment_col) || length(treatment_col) != 1L ||
        !nzchar(treatment_col)) {
      stop("`treatment_col` must be a single non-empty string.")
    }
  }

  # If the user asked to regress out percent.mt/percent.rb/percent.hb during
  # doublet calling but never asked for the corresponding column to be
  # computed, fail fast rather than have calldoublet() error later on a
  # missing column.
  .check_regress_col <- function(col, pattern, pattern_arg) {
    if (isTRUE(run_doublet_finder) &&
        !is.null(doublet_vars_to_regress) &&
        col %in% doublet_vars_to_regress &&
        is.null(pattern)) {
      stop("`doublet_vars_to_regress` requests '", col, "' but `", pattern_arg,
           "` is NULL, so ", col, " won't exist. Set `", pattern_arg,
           "` or remove '", col, "' from `doublet_vars_to_regress`.")
    }
  }
  .check_regress_col(mt_col, mt_pattern, "mt_pattern")
  .check_regress_col(rb_col, rb_pattern, "rb_pattern")
  .check_regress_col(hb_col, hb_pattern, "hb_pattern")

  # ---- Helper: build a single Seurat object from one path ------------------
  .make_one <- function(path, treatment) {
    dge_dir <- file.path(path, "DGE_filtered")
    if (!dir.exists(dge_dir)) {
      stop("No 'DGE_filtered' directory found at: ", path)
    }

    # 1. Read and orient the sparse count matrix (genes x cells).
    counts <- Matrix::readMM(file.path(dge_dir, "count_matrix.mtx"))
    counts <- Matrix::t(counts)
    counts <- methods::as(counts, "CsparseMatrix")

    # 2. Read gene and cell metadata.
    genes <- utils::read.csv(file.path(dge_dir, "all_genes.csv"),
                             stringsAsFactors = FALSE)
    cells <- utils::read.csv(file.path(dge_dir, "cell_metadata.csv"),
                             stringsAsFactors = FALSE)

    # 3. Assign dimnames (ensure they are unique).
    rownames(counts) <- make.unique(genes$gene_name)
    colnames(counts) <- cells$bc_wells

    # 4. Build the RNA assay and wrap in a Seurat object.
    assay <- SeuratObject::CreateAssayObject(
      counts,
      assay        = "RNA",
      min.features = mincellfeat,
      min.cells    = ceiling(nrow(cells) * mincellfrac)
    )

    meta <- tibble::column_to_rownames(cells, "bc_wells")
    # Tag every cell with the treatment label, if provided. Done on the full
    # cell metadata so that any cells dropped by the assay filters are
    # automatically dropped here too when CreateSeuratObject reconciles them.
    if (!is.null(treatment)) {
      meta[[treatment_col]] <- treatment
    }

    obj <- SeuratObject::CreateSeuratObject(assay, meta.data = meta)

    # Convert the RNA assay to Seurat v5's Assay5 class so downstream code
    # (e.g. multi-layer merges, IntegrateLayers) operates on the v5 API.
    obj[["RNA"]] <- methods::as(obj[["RNA"]], Class = "Assay5")

    # Compute mitochondrial/ribosomal/hemoglobin percentages. Done after
    # object creation so each pattern is matched against the post-filter
    # feature set and the result is automatically aligned to the surviving
    # cells. Each is independently opt-in (NULL pattern = skip) and warns
    # rather than silently adding an all-zero column when nothing matches.
    .add_pct <- function(obj, pattern, col, label) {
      if (is.null(pattern)) return(obj)
      n_match <- sum(grepl(pattern, rownames(obj)))
      if (n_match > 0) {
        obj[[col]] <- Seurat::PercentageFeatureSet(obj, pattern = pattern)
      } else {
        warning("No ", label, " features matched pattern '", pattern,
                "' in '", path, "'; ", col, " not added.")
      }
      obj
    }
    obj <- .add_pct(obj, mt_pattern, mt_col, "mitochondrial")
    obj <- .add_pct(obj, rb_pattern, rb_col, "ribosomal-protein")
    obj <- .add_pct(obj, hb_pattern, hb_col, "hemoglobin")

    obj
  }

  # ---- Build object(s) -----------------------------------------------------
  message(sprintf('--- Reading %d Parse Biosciences sample(s)%s ---',
                  length(paths),
                  if (workers > 1) sprintf(', %d parallel workers', workers) else ''))
  # Each path's read+build is fully independent, so this parallelizes
  # cleanly. Zip over paths/treatments directly (future_mapply/mapply)
  # rather than indexing into captured vectors from inside the worker, so
  # each worker only ever receives the one path/treatment pair it needs.
  treatment_list <- if (is.null(treatments)) {
    vector("list", length(paths))
  } else {
    as.list(treatments)
  }
  objs <- if (workers > 1) {
    future.apply::future_mapply(.make_one, paths, treatment_list,
                                SIMPLIFY = FALSE, future.seed = TRUE)
  } else {
    mapply(.make_one, paths, treatment_list, SIMPLIFY = FALSE)
  }
  names(objs) <- sample_names

  # ---- Doublet detection --------------------------------------------------
  if (isTRUE(run_doublet_finder)) {
    message(sprintf('--- Calling doublets with DoubletFinder (%s)%s ---',
                    doublet_normalization,
                    if (workers > 1) sprintf(', %d parallel workers', workers) else ''))
    sample_labels <- names(objs)
    n_samples     <- length(objs)

    # As above: mapply/future_mapply over `objs` directly rather than
    # indexing into a captured list from inside the worker, so each worker
    # only receives the one sample it's processing instead of the whole list.
    .call_one <- function(obj, i, lab) {
      # Per-sample progress messages only make sense when running
      # sequentially -- with workers > 1 these run in background sessions
      # and wouldn't surface here in order anyway.
      if (workers == 1) {
        message(sprintf('  [%d/%d] %s', i, n_samples, lab))
      }
      out <- calldoublet(obj,
                         samplenameIndex    = i,
                         normalization      = doublet_normalization,
                         vars.to.regress    = doublet_vars_to_regress,
                         cluster_resolution = doublet_cluster_resolution)
      if (isTRUE(filter_doublets)) {
        n_before <- ncol(out)
        out      <- subset(out, doublet_finder == "Singlet")
        if (workers == 1) {
          message(sprintf('    %s: dropped %d doublets (%d singlets remaining)',
                          lab, n_before - ncol(out), ncol(out)))
        }
      }
      out
    }

    results <- if (workers > 1) {
      future.apply::future_mapply(
        .call_one, objs, seq_len(n_samples), sample_labels,
        SIMPLIFY = FALSE, future.seed = TRUE
      )
    } else {
      mapply(.call_one, objs, seq_len(n_samples), sample_labels,
             SIMPLIFY = FALSE)
    }
    objs <- setNames(results, sample_labels)
  }

  if (length(objs) == 1L) {
    return(objs[[1L]])
  }
  objs
}
