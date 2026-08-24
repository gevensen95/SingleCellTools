#' Combine Resequenced Parse Biosciences Rounds by Cell Barcode
#'
#' Combines two sequencing rounds of the *same* Parse Biosciences sublibrary
#' (i.e. round 2 is additional sequencing depth on the same physical
#' barcoded library, not a fresh library prep) into a single synthetic
#' \code{DGE_filtered/} directory per sample, suitable for feeding straight
#' into \code{\link{MakeParseObj}}. This lets filtering, QC, and doublet
#' calling run exactly once on the combined data, instead of independently
#' per round.
#'
#' For each sample, cell barcodes (\code{bc_wells}) present in \strong{both}
#' rounds have their raw UMI counts summed gene-by-gene (as if the library
#' had simply been sequenced deeper in one run); barcodes present in only one
#' round are carried through unchanged. Genes are unioned across rounds
#' (zero-filled where a gene is present in one round's reference but not the
#' other's).
#'
#' For the CellRanger equivalent of this same "resequenced, not re-prepped"
#' combining problem, see \code{\link{CombineCellRangerRounds}} -- the two
#' functions share the same union/reindex/sum core internally, differing
#' only in how each platform's files are read and written.
#'
#' @section A caveat on combining already-deduplicated counts:
#' \code{DGE_filtered/count_matrix.mtx} is already UMI-deduplicated
#' \emph{within each round independently}. If the very same transcript
#' molecule happened to be sequenced (and therefore counted) in both round 1
#' and round 2's read sets, summing the two rounds' post-dedup counts will
#' count that molecule twice. The only way to avoid this entirely is to
#' concatenate the raw FASTQs from both rounds and re-run alignment/dedup
#' once on the combined read pool -- which this function exists precisely
#' because that isn't an option here. In practice the double-counted
#' fraction is small unless per-cell sequencing saturation is already high
#' (i.e. you were close to having sequenced every molecule in the library
#' at least once even before adding round 2), but it is a real approximation
#' and not a mathematically exact reconstruction of "one deeper sequencing
#' run." Keep this in mind if you see cells with implausibly high UMI counts
#' after combining.
#'
#' @param round1_paths Character vector of paths to round 1 Parse sample
#'   directories (each must contain a \code{DGE_filtered/} subdirectory, same
#'   as \code{\link{MakeParseObj}}'s \code{paths} argument).
#' @param round2_paths Character vector of paths to round 2 Parse sample
#'   directories, the same length as \code{round1_paths} and in the same
#'   sample order (i.e. \code{round1_paths[i]} and \code{round2_paths[i]}
#'   must be the same biological sample / same sublibrary, just two rounds
#'   of sequencing).
#' @param sample_names Optional character vector of sample names, same
#'   length as \code{round1_paths}. If \code{NULL} (default),
#'   \code{basename(round1_paths)} is used. Used to name the output
#'   subdirectories and as the names of the returned list elements.
#' @param output_dir Directory in which to write one synthetic sample
#'   directory per sample (\code{output_dir/<sample_name>/DGE_filtered/}).
#'   Created if it doesn't already exist.
#' @param gene_key_col Column in \code{all_genes.csv} used to match genes
#'   between the two rounds. Default \code{NULL}, which uses
#'   \code{"gene_id"} when that column is present in both rounds' files
#'   (safer -- unique per gene even if two genes share a display name), and
#'   falls back to \code{"gene_name"} with a warning otherwise.
#' @param metadata_conflict One of \code{"warn"} (default) or \code{"error"}.
#'   For a barcode present in both rounds, non-count \code{cell_metadata.csv}
#'   columns (e.g. \code{species}) are taken from round 1. If round 2 has a
#'   different value for that column on that barcode, \code{"warn"} prints a
#'   message and keeps round 1's value; \code{"error"} stops immediately so
#'   you can look into why the two rounds disagree.
#' @param overwrite Logical; if \code{FALSE} (default) and
#'   \code{output_dir/<sample_name>/DGE_filtered/} already exists, the
#'   function stops rather than silently overwriting a previous combine run.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{\code{paths}}{Character vector of the newly-written sample
#'     directories, in the same order as \code{sample_names} -- pass this
#'     straight to \code{MakeParseObj(paths = combined$paths, ...)}.}
#'   \item{\code{sample_names}}{The sample names used, same order as
#'     \code{paths}.}
#'   \item{\code{summary}}{A data frame, one row per sample, with columns
#'     \code{sample_name}, \code{n_cells_round1}, \code{n_cells_round2},
#'     \code{n_cells_shared}, \code{n_cells_round1_only},
#'     \code{n_cells_round2_only}, \code{n_cells_combined},
#'     \code{n_genes_round1}, \code{n_genes_round2}, \code{n_genes_combined}
#'     -- use this to sanity-check that the overlap between rounds looks
#'     like what you'd expect for a resequenced (not re-prepped) library.}
#' }
#'
#' @examples
#' \dontrun{
#' combined <- CombineParseRounds(
#'   round1_paths = file.path("/data/parse/round1", c("SampleA", "SampleB")),
#'   round2_paths = file.path("/data/parse/round2", c("SampleA", "SampleB")),
#'   sample_names = c("SampleA", "SampleB"),
#'   output_dir   = "/data/parse/combined"
#' )
#' print(combined$summary)
#'
#' obj_list <- MakeParseObj(
#'   combined$paths,
#'   sample_names = combined$sample_names
#' )
#' }
#'
#' @seealso \code{\link{CombineCellRangerRounds}}
#' @export
CombineParseRounds <- function(round1_paths,
                               round2_paths,
                               sample_names      = NULL,
                               output_dir,
                               gene_key_col      = NULL,
                               metadata_conflict = c("warn", "error"),
                               overwrite         = FALSE) {

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
  # .write_one(). Combining is expensive, and discovering on sample 5 that its
  # output already exists would mean throwing away the work already done on
  # samples 1-4 (which are, by then, already written to disk).
  if (!overwrite) {
    existing <- file.path(output_dir, sample_names, "DGE_filtered")
    clash    <- dir.exists(existing)
    if (any(clash)) {
      stop("Output director", if (sum(clash) == 1) "y" else "ies",
           " already exist:\n  ",
           paste(existing[clash], collapse = "\n  "),
           "\nPass overwrite = TRUE to replace them.")
    }
  }

  # ---- Helper: read one round's raw DGE_filtered output ---------------------
  # Mirrors MakeParseObj()'s internal reader (same orientation, same
  # dimname assignment) so downstream behavior stays consistent. Returns the
  # generic list(counts, genes, cells, path) shape .combine_one_core()
  # (in combine-rounds-utils.R) expects.
  .read_round <- function(path) {
    dge_dir <- file.path(path, "DGE_filtered")
    if (!dir.exists(dge_dir)) {
      stop("No 'DGE_filtered' directory found at: ", path)
    }
    # Name the missing file explicitly. Without this, a DGE_filtered/ that is
    # present but incomplete fails inside readMM()/read.csv() with a bare
    # "cannot open the connection", which names neither the file nor the
    # sample it belongs to.
    required <- c("count_matrix.mtx", "all_genes.csv", "cell_metadata.csv")
    absent   <- required[!file.exists(file.path(dge_dir, required))]
    if (length(absent)) {
      stop("Incomplete 'DGE_filtered' directory at ", path, " -- missing: ",
           paste(absent, collapse = ", "), ".")
    }
    counts <- Matrix::readMM(file.path(dge_dir, "count_matrix.mtx"))
    counts <- Matrix::t(counts)
    counts <- methods::as(counts, "CsparseMatrix")

    genes <- utils::read.csv(file.path(dge_dir, "all_genes.csv"),
                             stringsAsFactors = FALSE)
    cells <- utils::read.csv(file.path(dge_dir, "cell_metadata.csv"),
                             stringsAsFactors = FALSE)

    if (nrow(genes) != nrow(counts)) {
      stop("'all_genes.csv' at ", path, " has ", nrow(genes),
           " rows but count_matrix.mtx has ", nrow(counts), " genes.")
    }
    if (nrow(cells) != ncol(counts)) {
      stop("'cell_metadata.csv' at ", path, " has ", nrow(cells),
           " rows but count_matrix.mtx has ", ncol(counts), " cells.")
    }
    if (!"bc_wells" %in% colnames(cells)) {
      stop("'cell_metadata.csv' at ", path, " has no 'bc_wells' column.")
    }
    if (anyDuplicated(cells$bc_wells)) {
      stop("'cell_metadata.csv' at ", path, " has duplicate bc_wells values ",
           "-- can't use them as unique cell identifiers.")
    }
    colnames(counts) <- cells$bc_wells

    list(counts = counts, genes = genes, cells = cells, path = path)
  }

  # ---- Helper: combine one sample's two rounds -------------------------------
  # Resolves the Parse-specific gene_id/gene_name fallback (platform-specific,
  # so it lives here rather than in the shared core), then hands off to
  # .combine_one_core() for the actual union/reindex/sum work.
  .combine_one <- function(p1, p2, sample_name) {
    message(sprintf("--- Combining rounds for sample '%s' ---", sample_name))
    r1 <- .read_round(p1)
    r2 <- .read_round(p2)

    key_col <- gene_key_col
    if (is.null(key_col)) {
      if ("gene_id" %in% colnames(r1$genes) && "gene_id" %in% colnames(r2$genes)) {
        key_col <- "gene_id"
      } else {
        key_col <- "gene_name"
        warning(sample_name, ": no 'gene_id' column in all_genes.csv; ",
                "matching genes between rounds by 'gene_name' instead. If ",
                "two different genes share a display name in either round, ",
                "they will be incorrectly merged as one.")
      }
    }
    if (!key_col %in% colnames(r1$genes) || !key_col %in% colnames(r2$genes)) {
      stop(sample_name, ": gene_key_col '", key_col,
           "' not found in both rounds' all_genes.csv.")
    }

    .combine_one_core(r1, r2, sample_name,
                      gene_key_col      = key_col,
                      cell_id_col       = "bc_wells",
                      metadata_conflict = metadata_conflict)
  }

  # ---- Helper: write one sample's combined data as native Parse output -----
  .write_one <- function(sample_name, combined) {
    sample_dir <- file.path(output_dir, sample_name)
    dge_dir    <- file.path(sample_dir, "DGE_filtered")
    if (dir.exists(dge_dir) && !overwrite) {
      stop("Output directory already exists: ", dge_dir,
           "\nPass overwrite = TRUE to replace it.")
    }
    dir.create(dge_dir, recursive = TRUE, showWarnings = FALSE)

    # Native Parse orientation is cells x genes -- transpose back.
    Matrix::writeMM(Matrix::t(combined$counts),
                    file.path(dge_dir, "count_matrix.mtx"))
    utils::write.csv(combined$genes, file.path(dge_dir, "all_genes.csv"),
                     row.names = FALSE)
    utils::write.csv(combined$cells, file.path(dge_dir, "cell_metadata.csv"),
                     row.names = FALSE)
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
