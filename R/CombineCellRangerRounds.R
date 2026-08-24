#' Combine Resequenced CellRanger Rounds by Cell Barcode
#'
#' The CellRanger counterpart to \code{\link{CombineParseRounds}}: combines
#' two sequencing rounds of the \emph{same} 10X/CellRanger library (i.e.
#' round 2 is additional sequencing depth on the same physical GEM run, not a
#' fresh library prep) into a single synthetic barcodes/features/matrix
#' triplet per sample, suitable for feeding straight into
#' \code{\link{CreateRNAObjects}}. The two functions share the same
#' union/reindex/sum core internally and differ only in how each platform's
#' files are read and written.
#'
#' For each sample, cell barcodes present in \strong{both} rounds have their
#' raw UMI counts summed gene-by-gene (as if the library had simply been
#' sequenced deeper in one run); barcodes present in only one round are
#' carried through unchanged. Genes are unioned across rounds (zero-filled
#' where a gene is present in one round's reference but not the other's).
#'
#' @section Reading, writing, and what's different from Parse:
#' Each round is read with this package's existing CellRanger reader (the
#' same one \code{\link{CreateRNAObjects}} uses), so it accepts the same
#' inputs that already work there: a directory containing a bare
#' barcodes/features/matrix triplet, a \code{filtered_feature_bc_matrix/}
#' subdirectory, or an \code{outs/filtered_feature_bc_matrix/} subdirectory
#' -- gzipped or not, with or without a GEO-style sample-name prefix on the
#' filenames. Note that \code{.h5} input, which \code{CreateRNAObjects()}
#' does accept, is \emph{not} supported here: unpack it to a triplet first.
#' Unlike Parse's
#' \code{cell_metadata.csv}, standard CellRanger output carries no per-cell
#' metadata table beyond the count matrix itself, so there is nothing to
#' reconcile there -- \code{metadata_conflict} is accepted for API symmetry
#' with \code{CombineParseRounds()} but will rarely have anything to act on
#' in practice.
#'
#' The gene key used to match genes between rounds is whatever this
#' package's CellRanger reader already uses as the count matrix's row names
#' -- gene \emph{symbol} by default (CellRanger's \code{features.tsv}
#' column 2), not the Ensembl ID in column 1, since the underlying reader
#' doesn't override \code{Seurat::ReadMtx()}'s default \code{feature.column}.
#' Two distinct genes sharing a display symbol in either round's reference
#' would be incorrectly merged as one, exactly as the analogous
#' \code{gene_name} fallback in \code{CombineParseRounds()} warns about.
#' Combined output is written with matching \code{gene_id} and
#' \code{gene_name} columns (both set to that same symbol) since the
#' original Ensembl ID isn't retained through this path.
#'
#' Combined output is written as an uncompressed
#' \code{matrix.mtx}/\code{features.tsv}/\code{barcodes.tsv} triplet under
#' \code{output_dir/<sample_name>/filtered_feature_bc_matrix/} -- this
#' package's CellRanger reader matches those filenames whether or not
#' they're gzipped, so no compression step is needed for the result to be
#' readable again.
#'
#' @section A caveat on combining already-deduplicated counts:
#' CellRanger's \code{matrix.mtx} is already UMI-deduplicated \emph{within
#' each round independently}. If the very same transcript molecule happened
#' to be sequenced (and therefore counted) in both round 1 and round 2's
#' read sets, summing the two rounds' post-dedup counts will count that
#' molecule twice. The only way to avoid this entirely is to concatenate the
#' raw FASTQs from both rounds and re-run alignment/dedup once on the
#' combined read pool -- which this function exists precisely because that
#' isn't an option here. In practice the double-counted fraction is small
#' unless per-cell sequencing saturation is already high, but it is a real
#' approximation and not a mathematically exact reconstruction of "one
#' deeper sequencing run." Keep this in mind if you see cells with
#' implausibly high UMI counts after combining.
#'
#' @param round1_paths Character vector of paths to round 1 CellRanger
#'   sample directories (anything this package's CellRanger reader already
#'   accepts -- see Details).
#' @param round2_paths Character vector of paths to round 2 CellRanger
#'   sample directories, the same length as \code{round1_paths} and in the
#'   same sample order (i.e. \code{round1_paths[i]} and
#'   \code{round2_paths[i]} must be the same biological sample / same GEM
#'   run, just two rounds of sequencing).
#' @param sample_names Optional character vector of sample names, same
#'   length as \code{round1_paths}. If \code{NULL} (default),
#'   \code{basename(round1_paths)} is used. Used to name the output
#'   subdirectories and as the names of the returned list elements.
#' @param output_dir Directory in which to write one synthetic sample
#'   directory per sample
#'   (\code{output_dir/<sample_name>/filtered_feature_bc_matrix/}). Created
#'   if it doesn't already exist.
#' @param strip_barcode_suffix Logical, default \code{TRUE}. CellRanger
#'   barcodes carry a GEM-well suffix (e.g. \code{"-1"}); two independent
#'   runs of the same physical library will typically both emit \code{"-1"}
#'   regardless of round, and a leftover suffix on a now-merged library
#'   isn't meaningful anyway. When \code{TRUE}, the suffix is stripped
#'   before matching barcodes across rounds and in the combined output (an
#'   error is raised if stripping it creates duplicate barcodes within a
#'   round). Set \code{FALSE} to match/keep barcodes exactly as CellRanger
#'   wrote them.
#' @param metadata_conflict One of \code{"warn"} (default) or \code{"error"}.
#'   Accepted for symmetry with \code{CombineParseRounds()}; see Details --
#'   standard CellRanger output has no per-cell metadata columns for this to
#'   act on.
#' @param overwrite Logical; if \code{FALSE} (default) and
#'   \code{output_dir/<sample_name>/filtered_feature_bc_matrix/} already
#'   exists, the function stops rather than silently overwriting a previous
#'   combine run.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{paths}}{Character vector of the newly-written sample
#'     directories, in the same order as \code{sample_names} -- pass this
#'     straight to \code{CreateRNAObjects(data_dirs = combined$paths, ...)}.}
#'   \item{\code{sample_names}}{The sample names used, same order as
#'     \code{paths}.}
#'   \item{\code{summary}}{A data frame, one row per sample -- same columns
#'     as \code{CombineParseRounds()}'s summary -- use this to sanity-check
#'     that the overlap between rounds looks like what you'd expect for a
#'     resequenced (not re-prepped) library.}
#' }
#'
#' @examples
#' \dontrun{
#' combined <- CombineCellRangerRounds(
#'   round1_paths = file.path("/data/cellranger/round1", c("SampleA", "SampleB")),
#'   round2_paths = file.path("/data/cellranger/round2", c("SampleA", "SampleB")),
#'   sample_names = c("SampleA", "SampleB"),
#'   output_dir   = "/data/cellranger/combined"
#' )
#' print(combined$summary)
#'
#' objs <- CreateRNAObjects(
#'   data_dirs = combined$paths,
#'   workers   = length(combined$paths)
#' )
#' }
#'
#' @seealso \code{\link{CombineParseRounds}}, \code{\link{CreateRNAObjects}}
#' @export
CombineCellRangerRounds <- function(round1_paths,
                                    round2_paths,
                                    sample_names          = NULL,
                                    output_dir,
                                    strip_barcode_suffix  = TRUE,
                                    metadata_conflict     = c("warn", "error"),
                                    overwrite             = FALSE) {

  metadata_conflict <- match.arg(metadata_conflict)

  # ---- Argument checks ------------------------------------------------------
  if (!is.character(round1_paths) || length(round1_paths) < 1) {
    stop("`round1_paths` must be a non-empty character vector.")
  }
  if (!is.character(round2_paths) || length(round2_paths) != length(round1_paths)) {
    stop("`round2_paths` must be a character vector the same length as ",
         "`round1_paths` (got ", length(round2_paths), " vs ",
         length(round1_paths), "). round1_paths[i] and round2_paths[i] must ",
         "be the same sample.")
  }
  missing1 <- !dir.exists(round1_paths)
  missing2 <- !dir.exists(round2_paths)
  if (any(missing1)) {
    stop("The following round 1 path(s) do not exist:\n  ",
         paste(round1_paths[missing1], collapse = "\n  "))
  }
  if (any(missing2)) {
    stop("The following round 2 path(s) do not exist:\n  ",
         paste(round2_paths[missing2], collapse = "\n  "))
  }

  if (is.null(sample_names)) {
    sample_names <- basename(round1_paths)
  } else if (length(sample_names) != length(round1_paths)) {
    stop("`sample_names` must be the same length as `round1_paths`.")
  }
  if (anyDuplicated(sample_names)) {
    stop("`sample_names` must be unique (duplicates would overwrite each ",
         "other's output directory).")
  }

  if (missing(output_dir) || !is.character(output_dir) || length(output_dir) != 1L) {
    stop("`output_dir` must be a single path.")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Check every sample's output directory UP FRONT rather than only inside
  # .write_one() -- see the matching comment in CombineParseRounds().
  if (!overwrite) {
    existing <- file.path(output_dir, sample_names, "filtered_feature_bc_matrix")
    clash    <- dir.exists(existing)
    if (any(clash)) {
      stop("Output director", if (sum(clash) == 1) "y" else "ies",
           " already exist:\n  ",
           paste(existing[clash], collapse = "\n  "),
           "\nPass overwrite = TRUE to replace them.")
    }
  }

  # ---- Helper: read one round's raw CellRanger output ------------------------
  # Reuses .read_10x_triplet() (defined in CreateRNAObjects.R, same package)
  # rather than re-implementing barcode/features/matrix parsing -- same
  # bare-dir / filtered_feature_bc_matrix / outs/filtered_feature_bc_matrix
  # search order CreateRNAObjects() itself uses. Returns the generic
  # list(counts, genes, cells, path) shape .combine_one_core() (in
  # combine-rounds-utils.R) expects.
  .read_round_cellranger <- function(path) {
    mat <- .read_10x_triplet(path)
    if (is.null(mat)) {
      mat <- .read_10x_triplet(file.path(path, "filtered_feature_bc_matrix"))
    }
    if (is.null(mat)) {
      mat <- .read_10x_triplet(file.path(path, "outs", "filtered_feature_bc_matrix"))
    }
    if (is.null(mat)) {
      stop("Could not find a barcodes/features/matrix triplet in '", path,
           "' (checked the directory itself, its 'filtered_feature_bc_matrix' ",
           "subdirectory, and 'outs/filtered_feature_bc_matrix'). ",
           ".h5 input isn't supported by this function -- unpack it to a ",
           "triplet first if that's what you have.")
    }

    barcodes <- colnames(mat)
    if (isTRUE(strip_barcode_suffix)) {
      stripped <- sub("-[0-9]+$", "", barcodes)
      if (anyDuplicated(stripped)) {
        stop("Stripping the GEM-well suffix (e.g. '-1') from barcodes at '",
             path, "' produced duplicate cell IDs -- pass ",
             "strip_barcode_suffix = FALSE and investigate before combining.")
      }
      barcodes <- stripped
      colnames(mat) <- barcodes
    }

    list(
      counts = mat,
      genes  = data.frame(gene_id = rownames(mat), stringsAsFactors = FALSE),
      cells  = data.frame(cell_id = barcodes, stringsAsFactors = FALSE),
      path   = path
    )
  }

  # ---- Helper: combine one sample's two rounds -------------------------------
  .combine_one <- function(p1, p2, sample_name) {
    message(sprintf("--- Combining rounds for sample '%s' ---", sample_name))
    r1 <- .read_round_cellranger(p1)
    r2 <- .read_round_cellranger(p2)

    .combine_one_core(r1, r2, sample_name,
                      gene_key_col      = "gene_id",
                      cell_id_col       = "cell_id",
                      metadata_conflict = metadata_conflict)
  }

  # ---- Helper: write one sample's combined data as a native CellRanger-style
  # triplet -------------------------------------------------------------------
  .write_one <- function(sample_name, combined) {
    sample_dir <- file.path(output_dir, sample_name)
    mtx_dir    <- file.path(sample_dir, "filtered_feature_bc_matrix")
    if (dir.exists(mtx_dir) && !overwrite) {
      stop("Output directory already exists: ", mtx_dir,
           "\nPass overwrite = TRUE to replace it.")
    }
    dir.create(mtx_dir, recursive = TRUE, showWarnings = FALSE)

    # Native 10X orientation is genes x cells, which combined$counts already
    # is (unlike Parse's cells x genes) -- no transpose needed here.
    Matrix::writeMM(combined$counts, file.path(mtx_dir, "matrix.mtx"))

    # gene_id and gene_name are written identically since the original
    # Ensembl ID isn't retained through this reading path -- see Details.
    utils::write.table(
      data.frame(gene_id      = combined$genes$gene_id,
                gene_name    = combined$genes$gene_id,
                feature_type = "Gene Expression"),
      file.path(mtx_dir, "features.tsv"),
      sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
    )
    utils::write.table(
      data.frame(barcode = combined$cells$cell_id),
      file.path(mtx_dir, "barcodes.tsv"),
      sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
    )
    sample_dir
  }

  # ---- Run for every sample --------------------------------------------------
  out_paths <- character(length(sample_names))
  summaries <- vector("list", length(sample_names))
  for (i in seq_along(sample_names)) {
    combined      <- .combine_one(round1_paths[i], round2_paths[i], sample_names[i])
    out_paths[i]  <- .write_one(sample_names[i], combined)
    summaries[[i]] <- combined$summary
  }
  summary_df <- do.call(rbind, summaries)

  message("--- Done. Overlap summary: ---")
  message(paste(utils::capture.output(print(summary_df, row.names = FALSE)),
                collapse = "\n"))

  list(paths = out_paths, sample_names = sample_names, summary = summary_df)
}
