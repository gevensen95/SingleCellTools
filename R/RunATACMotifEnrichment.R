#' Motif scanning + chromVAR accessibility deviations for ATAC data
#'
#' Wraps the standard Signac motif workflow -- scan peaks for
#' transcription-factor motif matches (\code{Signac::AddMotifs()}, via
#' \code{motifmatchr}/\code{TFBSTools}), then compute per-cell,
#' per-motif accessibility deviation scores (\code{Signac::RunChromVAR()},
#' via \code{chromVAR}) -- into one call. A natural companion to
#' \code{\link{RunATACWrapper}} (normalization + LSI) for going from peaks
#' to "which TF motifs are more/less accessible in which cells."
#'
#' @param obj A Seurat object with a \code{ChromatinAssay}.
#' @param genome A \code{BSgenome} object matching the peaks' genome build,
#'   passed to both \code{Signac::AddMotifs()} and
#'   \code{Signac::RunChromVAR()}.
#' @param peak_assay Name of the \code{ChromatinAssay}. \code{NULL}
#'   (default) auto-detects it -- errors if \code{obj} has zero or more
#'   than one.
#' @param pfm A \code{TFBSTools::PFMatrixList} of position-frequency
#'   matrices to scan for. \code{NULL} (default) fetches JASPAR2020's
#'   motif set for \code{species} (requires the \code{JASPAR2020}
#'   package) -- pass your own to use a different motif database.
#' @param species NCBI taxonomy ID used to filter the default JASPAR2020
#'   set when \code{pfm = NULL}. Default \code{9606} (human); use
#'   \code{10090} for mouse.
#' @param compute_motifs If \code{TRUE} (default), run
#'   \code{Signac::AddMotifs()} first. Set \code{FALSE} to skip (e.g. if
#'   \code{obj} already has a motif object attached) and go straight to
#'   \code{run_chromvar}.
#' @param run_chromvar If \code{TRUE} (default), run
#'   \code{Signac::RunChromVAR()} after motif scanning, adding a
#'   \code{"chromvar"} assay of per-cell motif deviation scores. Set
#'   \code{FALSE} to only scan motifs without the chromVAR step.
#' @param verbose Message progress. Default \code{TRUE}.
#' @return \code{obj} with a motif object attached to \code{peak_assay}
#'   (if \code{compute_motifs = TRUE}) and, if \code{run_chromvar = TRUE},
#'   a new \code{"chromvar"} assay of per-cell/per-motif deviation scores
#'   -- same as calling the underlying Signac functions directly.
#' @examples
#' \dontrun{
#' atac <- RunATACMotifEnrichment(atac, genome = BSgenome.Hsapiens.UCSC.hg38)
#' DefaultAssay(atac) <- "chromvar"
#' FeaturePlot(atac, features = "MA0139.1")  # CTCF, e.g.
#' }
#' @export
RunATACMotifEnrichment <- function(obj,
                                   genome,
                                   peak_assay     = NULL,
                                   pfm            = NULL,
                                   species        = 9606,
                                   compute_motifs = TRUE,
                                   run_chromvar   = TRUE,
                                   verbose        = TRUE) {

  .assert_seurat(obj)

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

  if (isTRUE(compute_motifs)) {
    if (!requireNamespace("motifmatchr", quietly = TRUE) ||
        !requireNamespace("TFBSTools", quietly = TRUE)) {
      stop("'motifmatchr' and 'TFBSTools' are required to compute motifs. ",
           "Install with BiocManager::install(c('motifmatchr', 'TFBSTools')).")
    }
    if (is.null(pfm)) {
      if (!requireNamespace("JASPAR2020", quietly = TRUE)) {
        stop("Pass `pfm` explicitly (a TFBSTools::PFMatrixList), or install ",
             "'JASPAR2020' (BiocManager::install('JASPAR2020')) to fetch a ",
             "default motif set.")
      }
      pfm <- TFBSTools::getMatrixSet(
        x = JASPAR2020::JASPAR2020,
        opts = list(species = species, all_versions = FALSE))
    }
    if (isTRUE(verbose)) {
      message(sprintf(
        "--- Scanning %d motif(s) against peaks (assay = '%s') ---",
        length(pfm), pa))
    }
    obj <- Signac::AddMotifs(obj, genome = genome, pfm = pfm, assay = pa)
  }

  if (isTRUE(run_chromvar)) {
    if (!requireNamespace("chromVAR", quietly = TRUE)) {
      stop("'chromVAR' is required for run_chromvar = TRUE. Install with ",
           "BiocManager::install('chromVAR').")
    }
    if (isTRUE(verbose)) {
      message(sprintf("--- Running chromVAR (assay = '%s') ---", pa))
    }
    obj <- Signac::RunChromVAR(obj, genome = genome, assay = pa)
  }

  obj
}
