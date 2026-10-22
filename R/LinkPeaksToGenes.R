#' Link ATAC peaks to nearby genes by co-variation with expression
#'
#' Thin wrapper around \code{Signac::RegionStats()} +
#' \code{Signac::LinkPeaks()}: computes each peak's GC content/length
#' (needed for \code{LinkPeaks()}'s background-matched null), then tests
#' each gene's expression against every peak within \code{distance} of it
#' for a significant correlation across cells. A natural companion to
#' \code{\link{RunATACWrapper}} for objects that have both a peak assay and
#' an RNA assay (e.g. multiome, or ATAC + RNA integrated onto the same
#' cells).
#'
#' @param obj A Seurat object with a \code{ChromatinAssay} (peaks) and a
#'   gene expression assay.
#' @param genome A \code{BSgenome} object matching the peaks' genome build
#'   (e.g. \code{BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38}),
#'   passed to \code{Signac::RegionStats()}. For a genome with no packaged
#'   \code{BSgenome} (e.g. a custom macaque T2T assembly built only from a
#'   cellranger-ref FASTA), pass an \code{Rsamtools::FaFile} instead --
#'   \code{RegionStats()}'s GC-content step only calls
#'   \code{Biostrings::getSeq()}/\code{seqlevels()} on \code{genome}, and
#'   both have \code{FaFile} methods, so no \code{BSgenome} package is
#'   required for this to work. A plain character path to a FASTA file is
#'   also accepted directly and wrapped in \code{Rsamtools::FaFile()} for
#'   you (building its \code{.fai} index if one doesn't already exist
#'   alongside it).
#' @param peak_assay Name of the \code{ChromatinAssay}. \code{NULL}
#'   (default) auto-detects it -- errors if \code{obj} has zero or more
#'   than one.
#' @param expression_assay Name of the gene expression assay. \code{NULL}
#'   (default) uses \code{"RNA"} if present.
#' @param genes.use Optional character vector restricting which genes to
#'   test. \code{NULL} (default) tests every gene in
#'   \code{expression_assay}.
#' @param distance Maximum peak-to-gene distance (bp) to test. Default
#'   \code{5e5} (500kb), \code{Signac::LinkPeaks()}'s own default.
#' @param verbose Passed to \code{Signac::LinkPeaks()}. Default \code{TRUE}.
#' @return \code{obj} with peak-gene links stored in
#'   \code{Signac::Links(obj[[peak_assay]])} (a \code{GRanges}), same as
#'   calling \code{Signac::LinkPeaks()} directly.
#' @examples
#' \dontrun{
#' multiome <- LinkPeaksToGenes(multiome, genome = BSgenome.Hsapiens.UCSC.hg38)
#' Signac::Links(multiome[["ATAC"]])
#' }
#' @export
LinkPeaksToGenes <- function(obj,
                             genome,
                             peak_assay       = NULL,
                             expression_assay = NULL,
                             genes.use        = NULL,
                             distance         = 5e5,
                             verbose          = TRUE) {

  .assert_seurat(obj)

  # `genome` can be a BSgenome object (the common case), an already-built
  # Rsamtools::FaFile (for a custom genome with no packaged BSgenome --
  # see the `genome` doc above for why that's enough for RegionStats()),
  # or a plain path to a FASTA, wrapped in FaFile() here for convenience.
  if (is.character(genome)) {
    if (!requireNamespace("Rsamtools", quietly = TRUE)) {
      stop("Passing `genome` as a file path requires the 'Rsamtools' ",
           "package (BiocManager::install('Rsamtools')) to build an ",
           "FaFile from it. Install it, or pass an already-built ",
           "Rsamtools::FaFile or BSgenome object instead.")
    }
    if (length(genome) != 1 || !file.exists(genome)) {
      stop("`genome` ('", paste(genome, collapse = ", "), "') is neither ",
           "a BSgenome object, an Rsamtools::FaFile, nor the path to an ",
           "existing FASTA file.")
    }
    fai_path <- paste0(genome, ".fai")
    if (!file.exists(fai_path)) {
      if (isTRUE(verbose)) {
        message(sprintf("--- Indexing FASTA (no .fai found): %s ---", genome))
      }
      tryCatch(
        Rsamtools::indexFa(genome),
        error = function(e) {
          stop("Failed to build a .fai index for '", genome, "': ",
               conditionMessage(e), ". If the containing directory is ",
               "read-only, index it yourself first (samtools faidx, or ",
               "Rsamtools::indexFa() somewhere writable) and pass the ",
               "resulting Rsamtools::FaFile object instead of a path.")
        }
      )
    }
    if (isTRUE(verbose)) {
      message(sprintf("--- Loading FASTA as Rsamtools::FaFile: %s ---", genome))
    }
    genome <- Rsamtools::FaFile(genome)
  }

  pa <- peak_assay
  if (is.null(pa)) {
    is_chromatin <- vapply(Seurat::Assays(obj),
                           function(a) inherits(obj[[a]], "ChromatinAssay"),
                           logical(1))
    candidates <- Seurat::Assays(obj)[is_chromatin]
    if (length(candidates) != 1) {
      stop("Could not auto-detect a single ChromatinAssay in `obj` ",
           "(found: ", paste(candidates, collapse = ", "),
           "). Pass `peak_assay` explicitly.")
    }
    pa <- candidates
  } else if (!inherits(obj[[pa]], "ChromatinAssay")) {
    stop("`peak_assay` ('", pa, "') is not a ChromatinAssay.")
  }

  ea <- expression_assay
  if (is.null(ea)) {
    if (!"RNA" %in% Seurat::Assays(obj)) {
      stop("Could not find an 'RNA' assay in `obj`. Pass ",
           "`expression_assay` explicitly.")
    }
    ea <- "RNA"
  } else if (!ea %in% Seurat::Assays(obj)) {
    stop("`expression_assay` ('", ea, "') not found in `obj`.")
  }

  if (isTRUE(verbose)) {
    message(sprintf("--- Computing peak GC content/length (assay = '%s') ---", pa))
  }
  obj <- Signac::RegionStats(obj, genome = genome, assay = pa)

  if (isTRUE(verbose)) {
    message(sprintf(
      "--- Linking peaks ('%s') to genes ('%s'), distance <= %s bp ---",
      pa, ea, format(distance, big.mark = ",", scientific = FALSE)))
  }
  obj <- Signac::LinkPeaks(obj, peak.assay = pa, expression.assay = ea,
                           genes.use = genes.use, distance = distance,
                           verbose = isTRUE(verbose))

  links <- tryCatch(Signac::Links(obj[[pa]]), error = function(e) NULL)
  if (isTRUE(verbose) && !is.null(links)) {
    message(sprintf("  %d significant peak-gene link(s) found.", length(links)))
  }

  obj
}
