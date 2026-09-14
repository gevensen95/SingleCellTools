#' Peak-peak co-accessibility network (Cicero) -- a second, RNA-free line of
#' enhancer-promoter evidence
#'
#' Wraps \code{cicero::make_cicero_cds()} + \code{cicero::run_cicero()} to
#' build a co-accessibility network directly from ATAC data: which peaks
#' tend to be open together across cells, independent of any paired RNA
#' measurement. This complements \code{\link{LinkPeaksToGenes}} (which links
#' peaks to genes via expression correlation) rather than replacing it --
#' use this when you don't trust, or don't have, paired RNA to link off of,
#' or want a second source of evidence before calling a peak pair a real
#' enhancer-promoter link.
#'
#' If \code{link_promoters = TRUE} (default) and \code{\link{AnnotatePeaks}}
#' has already been run on \code{obj} (so \code{obj@misc$peak_annotation}
#' exists), highly co-accessible pairs where exactly one peak is
#' \code{"promoter_proximal"} are additionally resolved to a gene (that
#' promoter peak's \code{nearest_gene}) and reported separately as
#' \code{obj@misc$cicero$enhancer_promoter_links} -- the distal member of
#' each such pair is the candidate enhancer.
#'
#' @param obj A Seurat object with a \code{ChromatinAssay}.
#' @param peak_assay Name of the \code{ChromatinAssay}. \code{NULL}
#'   (default) auto-detects it -- errors if \code{obj} has zero or more
#'   than one.
#' @param genome_size A 2-column data.frame/matrix of chromosome name and
#'   length (bp) for the genome build your peaks are called against (e.g.
#'   \code{data.frame(V1 = c("chr1","chr2",...), V2 = c(248956422, ...))}).
#'   Required -- there is no default, since guessing a genome size from
#'   peak coordinates alone risks silently using the wrong build.
#' @param reduction A 2D cell embedding to seed Cicero's UMAP/tSNE-based
#'   cell aggregation. \code{NULL} (default) uses \code{"lsi"} if present
#'   on \code{obj}, else \code{"umap"}; errors if neither is found.
#' @param distance_threshold Max peak-peak distance (bp) Cicero will test.
#'   Default 5e5 (500kb), matching \code{\link{LinkPeaksToGenes}}'s default.
#' @param link_promoters If \code{TRUE} (default), also resolve high-
#'   co-accessibility enhancer-promoter pairs into
#'   \code{obj@misc$cicero$enhancer_promoter_links} -- requires
#'   \code{\link{AnnotatePeaks}} to have been run first; skipped with a
#'   message otherwise.
#' @param coaccess_cutoff Minimum co-accessibility score for a pair to be
#'   included in \code{enhancer_promoter_links}. Default 0.25 (Cicero's own
#'   commonly-used threshold).
#' @param verbose Message progress. Default \code{TRUE}.
#' @return \code{obj} with \code{obj@misc$cicero} populated:
#'   \code{list(cicero_cds, connections, enhancer_promoter_links)} --
#'   \code{connections} is Cicero's full \code{Peak1}/\code{Peak2}/
#'   \code{coaccess} table; \code{enhancer_promoter_links} (if computed) adds
#'   \code{gene}/\code{enhancer_peak}/\code{promoter_peak}/\code{coaccess}.
#' @examples
#' \dontrun{
#' atac <- AnnotatePeaks(atac)
#' hg38_sizes <- data.frame(V1 = names(seqlengths(BSgenome.Hsapiens.UCSC.hg38)),
#'                          V2 = seqlengths(BSgenome.Hsapiens.UCSC.hg38))
#' atac <- RunCicero(atac, genome_size = hg38_sizes)
#' head(atac@misc$cicero$enhancer_promoter_links)
#' }
#' @importFrom Seurat Embeddings DefaultAssay
#' @export
RunCicero <- function(obj,
                      peak_assay         = NULL,
                      genome_size        = NULL,
                      reduction          = NULL,
                      distance_threshold = 5e5,
                      link_promoters     = TRUE,
                      coaccess_cutoff    = 0.25,
                      verbose            = TRUE) {

  .assert_seurat(obj)
  if (!requireNamespace("cicero", quietly = TRUE)) {
    stop("'cicero' is required. Install with BiocManager::install('cicero').")
  }
  if (!requireNamespace("monocle3", quietly = TRUE) ||
      !requireNamespace("SeuratWrappers", quietly = TRUE)) {
    stop("'monocle3' and 'SeuratWrappers' are required. Install with ",
         "remotes::install_github(c('cole-trapnell-lab/monocle3', ",
         "'satijalab/seurat-wrappers')).")
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("'SingleCellExperiment' is required. Install with ",
         "BiocManager::install('SingleCellExperiment').")
  }
  if (is.null(genome_size)) {
    stop("`genome_size` is required: a 2-column data.frame of chromosome ",
         "name/length matching your peaks' genome build. There is no ",
         "default -- guessing a build from peak coordinates risks silently ",
         "using the wrong one.")
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

  rd <- reduction
  if (is.null(rd)) {
    rd <- if ("lsi" %in% names(obj@reductions)) "lsi" else
          if ("umap" %in% names(obj@reductions)) "umap" else NULL
    if (is.null(rd)) {
      stop("No 'lsi' or 'umap' reduction found on `obj`. Compute one, or ",
           "pass `reduction` explicitly.")
    }
  } else if (!(rd %in% names(obj@reductions))) {
    stop("Reduction '", rd, "' not found.")
  }
  emb <- Seurat::Embeddings(obj, reduction = rd)[, 1:2, drop = FALSE]

  if (isTRUE(verbose)) message(sprintf("--- Building cell_data_set (assay = '%s') ---", pa))
  old_default <- Seurat::DefaultAssay(obj)
  Seurat::DefaultAssay(obj) <- pa
  cds <- SeuratWrappers::as.cell_data_set(obj)
  Seurat::DefaultAssay(obj) <- old_default
  SingleCellExperiment::reducedDim(cds, "UMAP") <- emb

  if (isTRUE(verbose)) message("--- Aggregating cells (cicero::make_cicero_cds) ---")
  cicero_cds <- cicero::make_cicero_cds(cds, reduced_coordinates = emb)

  if (isTRUE(verbose)) {
    message(sprintf("--- Running cicero::run_cicero (distance <= %s bp) ---",
                    format(distance_threshold, big.mark = ",", scientific = FALSE)))
  }
  conns <- cicero::run_cicero(cicero_cds, genomic_coords = genome_size,
                              window = distance_threshold)
  conns <- as.data.frame(conns)
  conns <- conns[!is.na(conns$coaccess), , drop = FALSE]
  conns <- conns[order(-conns$coaccess), ]
  rownames(conns) <- NULL

  out <- list(cicero_cds = cicero_cds, connections = conns)

  if (isTRUE(link_promoters)) {
    ann <- obj@misc$peak_annotation
    if (is.null(ann)) {
      if (isTRUE(verbose)) {
        message("  `link_promoters = TRUE` but obj@misc$peak_annotation ",
               "not found -- run AnnotatePeaks(obj) first. Skipping ",
               "enhancer-promoter resolution; raw `connections` still returned.")
      }
    } else {
      sig <- conns[conns$coaccess >= coaccess_cutoff, , drop = FALSE]
      type1 <- ann$peak_type[match(sig$Peak1, ann$peak)]
      type2 <- ann$peak_type[match(sig$Peak2, ann$peak)]
      is_prom1 <- !is.na(type1) & type1 == "promoter_proximal"
      is_prom2 <- !is.na(type2) & type2 == "promoter_proximal"
      keep <- xor(is_prom1, is_prom2)
      sig <- sig[keep, , drop = FALSE]
      is_prom1 <- is_prom1[keep]

      if (nrow(sig) > 0) {
        promoter_peak  <- ifelse(is_prom1, sig$Peak1, sig$Peak2)
        enhancer_peak  <- ifelse(is_prom1, sig$Peak2, sig$Peak1)
        gene           <- ann$nearest_gene[match(promoter_peak, ann$peak)]
        out$enhancer_promoter_links <- data.frame(
          gene           = gene,
          enhancer_peak  = enhancer_peak,
          promoter_peak  = promoter_peak,
          coaccess       = sig$coaccess,
          stringsAsFactors = FALSE
        )
        out$enhancer_promoter_links <- out$enhancer_promoter_links[
          order(-out$enhancer_promoter_links$coaccess), ]
        rownames(out$enhancer_promoter_links) <- NULL
        if (isTRUE(verbose)) {
          message(sprintf("  %d enhancer-promoter link(s) resolved (coaccess >= %.2f).",
                          nrow(out$enhancer_promoter_links), coaccess_cutoff))
        }
      } else if (isTRUE(verbose)) {
        message("  No enhancer-promoter pairs cleared `coaccess_cutoff`.")
      }
    }
  }

  obj@misc$cicero <- out
  obj
}
