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
#' @param mt_pattern Optional regex passed to
#'   \code{\link[Seurat]{PercentageFeatureSet}} to identify mitochondrial
#'   genes. \code{NULL} (default) skips the calculation entirely. Pass a
#'   pattern such as \code{"^MT-"} (human), \code{"^mt-"} (mouse), or
#'   \code{"^[Mm][Tt]-"} (either) to add a \code{percent.mt} column.
#' @param mt_col Name of the metadata column in which to store the
#'   mitochondrial percentage when \code{mt_pattern} is supplied. Default
#'   \code{"percent.mt"}.
#' @param rb_pattern Optional regex passed to
#'   \code{\link[Seurat]{PercentageFeatureSet}} to identify ribosomal-protein
#'   genes. \code{NULL} (default) skips the calculation entirely, matching
#'   \code{mt_pattern}'s behavior. Pass a pattern such as
#'   \code{"^(Rp[sl]|RP[SL])"} to add a \code{percent.rb} column (matches
#'   both mouse and human gene symbol conventions).
#' @param rb_col Name of the metadata column in which to store the
#'   ribosomal-protein percentage when \code{rb_pattern} is supplied.
#'   Default \code{"percent.rb"}.
#' @param hb_pattern Optional regex passed to
#'   \code{\link[Seurat]{PercentageFeatureSet}} to identify hemoglobin
#'   genes. \code{NULL} (default) skips the calculation entirely, matching
#'   \code{mt_pattern}'s behavior. Pass a pattern such as
#'   \code{"^(Hb[^p]|HB[^P])"} to add a \code{percent.hb} column -- that
#'   pattern excludes the unrelated \code{"Hbp1"}/\code{"HBP1"} gene, a
#'   well-known false positive for naive \code{"^Hb"} patterns.
#' @param hb_col Name of the metadata column in which to store the
#'   hemoglobin percentage when \code{hb_pattern} is supplied. Default
#'   \code{"percent.hb"}.
#' @param run_doublet_finder Logical; if TRUE (default), run \code{calldoublet}
#'   on every object and add a \code{doublet_finder} metadata column.
#' @param doublet_normalization Passed to \code{calldoublet}: one of
#'   \code{"LogNormalize"} (default) or \code{"SCT"}.
#' @param doublet_vars_to_regress Passed to \code{calldoublet} as
#'   \code{vars.to.regress}. Default \code{NULL}; set to \code{mt_col} (e.g.
#'   \code{"percent.mt"}), \code{rb_col}, and/or \code{hb_col} if you
#'   supplied the corresponding \code{*_pattern} argument(s) and want to
#'   regress that content out during the doublet workflow's normalization
#'   step.
#' @param doublet_cluster_resolution Passed to \code{calldoublet} as
#'   \code{cluster_resolution}. Default \code{0.1}.
#' @param filter_doublets Logical; if TRUE, subset each object to
#'   \code{doublet_finder == "Singlet"} after doublet calling. Default
#'   \code{FALSE} so the doublet labels are preserved for downstream review.
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
#' # With doublet calling that regresses out percent.mt
#' obj_list <- MakeParseObj(
#'   sample_paths,
#'   mt_pattern              = "^[Mm][Tt]-",
#'   doublet_vars_to_regress = "percent.mt"
#' )
#'
#' # Also compute ribosomal/hemoglobin percentages
#' obj_list <- MakeParseObj(
#'   sample_paths,
#'   mt_pattern = "^[Mm][Tt]-",
#'   rb_pattern = "^(Rp[sl]|RP[SL])",
#'   hb_pattern = "^(Hb[^p]|HB[^P])"
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
                         mt_pattern                 = NULL,
                         mt_col                     = "percent.mt",
                         rb_pattern                 = NULL,
                         rb_col                     = "percent.rb",
                         hb_pattern                 = NULL,
                         hb_col                     = "percent.hb",
                         run_doublet_finder         = TRUE,
                         doublet_normalization      = c("LogNormalize", "SCT"),
                         doublet_vars_to_regress    = NULL,
                         doublet_cluster_resolution = 0.1,
                         filter_doublets            = FALSE) {

  doublet_normalization <- match.arg(doublet_normalization)

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
  # Iterate by index so we can pair each path with its treatment label.
  objs <- lapply(seq_along(paths), function(i) {
    .make_one(paths[[i]],
              treatment = if (is.null(treatments)) NULL else treatments[[i]])
  })
  names(objs) <- sample_names

  # ---- Doublet detection --------------------------------------------------
  if (isTRUE(run_doublet_finder)) {
    message(sprintf('--- Calling doublets with DoubletFinder (%s) ---',
                    doublet_normalization))
    objs <- setNames(lapply(seq_along(objs), function(i) {
      lab <- names(objs)[i]
      message(sprintf('  [%d/%d] %s', i, length(objs), lab))
      out <- calldoublet(objs[[i]],
                         samplenameIndex    = i,
                         normalization      = doublet_normalization,
                         vars.to.regress    = doublet_vars_to_regress,
                         cluster_resolution = doublet_cluster_resolution)
      if (isTRUE(filter_doublets)) {
        n_before <- ncol(out)
        out      <- subset(out, doublet_finder == "Singlet")
        message(sprintf('    %s: dropped %d doublets (%d singlets remaining)',
                        lab, n_before - ncol(out), ncol(out)))
      }
      out
    }), names(objs))
  }

  if (length(objs) == 1L) {
    return(objs[[1L]])
  }
  objs
}
