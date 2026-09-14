#' TF footprinting -- physical binding evidence at motif-matched sites
#'
#' Wraps \code{Signac::Footprint()} + \code{Signac::PlotFootprint()}: looks
#' at the actual Tn5 insertion pattern \emph{within} each motif-matched
#' region for the characteristic dip left by a bound protein physically
#' blocking transposition, rather than just relying on
#' \code{\link{RunATACMotifEnrichment}}'s sequence-match ("this DNA matches
#' the TF's motif") evidence. Meaningfully stronger evidence of direct
#' binding, at the cost of needing more reads per region and being slower
#' -- hence \code{motifs} is required rather than defaulting to "every motif
#' in the object."
#'
#' @param obj A Seurat object that has already been through
#'   \code{\link{RunATACMotifEnrichment}} (\code{compute_motifs = TRUE}).
#' @param motifs Character vector of TF names or motif IDs to footprint.
#'   Required -- footprinting every motif in the object is expensive and
#'   rarely what you want; pick the TF(s) you actually care about (e.g.
#'   ones that came out of \code{\link{RunSCENIC}} or
#'   \code{\link{RunERegulons}}).
#' @param genome A \code{BSgenome} object matching the peaks' genome build,
#'   same as \code{\link{RunATACMotifEnrichment}}'s \code{genome} argument.
#' @param peak_assay Name of the \code{ChromatinAssay}. \code{NULL}
#'   (default) auto-detects it -- errors if \code{obj} has zero or more
#'   than one.
#' @param in.peaks Passed to \code{Signac::Footprint()}: restrict to
#'   motif matches inside called peaks (\code{TRUE}, default) rather than
#'   anywhere in the genome.
#' @param plot Logical; if \code{TRUE} (default), also build the per-TF
#'   footprint plot via \code{Signac::PlotFootprint()}.
#' @param verbose Message progress. Default \code{TRUE}.
#' @return \code{obj} with footprinting data attached to \code{peak_assay}
#'   (same as calling \code{Signac::Footprint()} directly) and, if
#'   \code{plot = TRUE}, \code{obj@misc$footprint_plot} (the
#'   \code{PlotFootprint()} \code{ggplot}/patchwork object).
#' @examples
#' \dontrun{
#' atac <- RunATACMotifEnrichment(atac, genome = BSgenome.Hsapiens.UCSC.hg38)
#' atac <- RunTFFootprinting(atac, motifs = c("CTCF", "STAT1"),
#'                           genome = BSgenome.Hsapiens.UCSC.hg38)
#' atac@misc$footprint_plot
#' }
#' @export
RunTFFootprinting <- function(obj,
                              motifs,
                              genome,
                              peak_assay = NULL,
                              in.peaks   = TRUE,
                              plot       = TRUE,
                              verbose    = TRUE) {

  .assert_seurat(obj)
  if (missing(motifs) || is.null(motifs) || length(motifs) == 0) {
    stop("`motifs` (a character vector of TF names/motif IDs) is required ",
         "-- footprinting every motif in the object by default would be ",
         "needlessly expensive.")
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

  motif_obj <- tryCatch(Signac::Motifs(obj[[pa]]), error = function(e) NULL)
  if (is.null(motif_obj)) {
    stop("No motif data found on assay '", pa, "'. Run ",
         "RunATACMotifEnrichment(obj, ..., compute_motifs = TRUE) first.")
  }

  if (isTRUE(verbose)) {
    message(sprintf("--- Footprinting %d motif(s) (assay = '%s') ---",
                    length(motifs), pa))
  }
  obj <- Signac::Footprint(obj, assay = pa, motif.name = motifs,
                           genome = genome, in.peaks = in.peaks,
                           verbose = isTRUE(verbose))

  if (isTRUE(plot)) {
    obj@misc$footprint_plot <- Signac::PlotFootprint(obj, features = motifs)
  }
  obj
}
