#' Classify ATAC peaks by genomic context (promoter / genic / distal-enhancer)
#'
#' The foundation piece for anything that wants to talk about "enhancers"
#' rather than just "peaks": labels every peak in a \code{ChromatinAssay} by
#' its position relative to the nearest gene, so downstream functions
#' (\code{\link{RunERegulons}}, \code{\link{RunCicero}},
#' \code{\link{CallSuperEnhancers}}) can restrict themselves to genuinely
#' distal, non-promoter regions instead of treating every accessible peak
#' the same way.
#'
#' \strong{\code{method = "signac"} (default).} Uses
#' \code{Signac::ClosestFeature()} against whatever gene annotation is
#' already attached to the \code{ChromatinAssay} (\code{Signac::Annotation()}
#' -- set when the object was built, e.g. via \code{\link{CreateATACObjects}}'s
#' \code{genome} argument). No extra package or annotation database needed.
#' Peaks are classified as \code{"promoter_proximal"} (within
#' \code{promoter_window} bp of, but not overlapping, a gene),
#' \code{"exonic"}/\code{"intronic"} (overlapping a gene body), or
#' \code{"distal"} (everything else -- the putative-enhancer set). This is a
#' coarser scheme than \code{"chipseeker"} below (no 5'UTR/3'UTR/downstream
#' distinction), but needs nothing beyond what a Signac ChromatinAssay
#' already carries.
#'
#' \strong{\code{method = "chipseeker"}.} Uses
#' \code{ChIPseeker::annotatePeak()} against a \code{TxDb} you supply, for
#' the fuller genomic-feature breakdown (5'UTR, 3'UTR, downstream, etc.)
#' collapsed down to the same \code{promoter_proximal}/\code{exonic}/
#' \code{intronic}/\code{distal} labels for consistency with
#' \code{"signac"}'s output shape -- inspect \code{obj@misc$peak_annotation}'s
#' \code{full_annotation} column for ChIPseeker's original, finer category
#' if you want it.
#'
#' @param obj A Seurat object with a \code{ChromatinAssay}.
#' @param peak_assay Name of the \code{ChromatinAssay}. \code{NULL}
#'   (default) auto-detects it -- errors if \code{obj} has zero or more
#'   than one.
#' @param method \code{"signac"} (default) or \code{"chipseeker"}. See
#'   Details.
#' @param txdb \code{method = "chipseeker"} only: a \code{TxDb} object
#'   matching your genome build (e.g. \code{TxDb.Hsapiens.UCSC.hg38.knownGene}).
#'   Required for that method -- no default, no guessing the species.
#' @param promoter_window Max distance (bp) from a gene for a non-overlapping
#'   peak to be called \code{"promoter_proximal"} rather than
#'   \code{"distal"}. Default 2000.
#' @param verbose Message progress. Default \code{TRUE}.
#' @return \code{obj} with \code{obj@misc$peak_annotation} (a data frame,
#'   one row per peak: \code{peak}, \code{peak_type}, \code{is_enhancer},
#'   \code{nearest_gene}, \code{distance_to_tss}) and matching per-peak
#'   feature-metadata columns on the \code{ChromatinAssay} itself
#'   (\code{obj[[peak_assay]][["peak_type"]]}, etc.), so both
#'   \code{obj[[peak_assay]][[]]} and \code{obj@misc$peak_annotation} carry
#'   the same information.
#' @examples
#' \dontrun{
#' atac <- AnnotatePeaks(atac)
#' table(atac@misc$peak_annotation$peak_type)
#' enhancer_peaks <- rownames(atac@misc$peak_annotation)[atac@misc$peak_annotation$is_enhancer]
#' }
#' @export
AnnotatePeaks <- function(obj,
                          peak_assay      = NULL,
                          method          = c("signac", "chipseeker"),
                          txdb            = NULL,
                          promoter_window = 2000,
                          verbose         = TRUE) {

  method <- match.arg(method)
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

  peaks_gr <- Signac::granges(obj[[pa]])
  peak_names <- rownames(obj[[pa]])

  if (method == "signac") {
    ann <- tryCatch(Signac::Annotation(obj[[pa]]), error = function(e) NULL)
    if (is.null(ann) || length(ann) == 0) {
      stop("No gene annotation attached to assay '", pa, "' ",
           "(Signac::Annotation(obj[['", pa, "']]) is empty). Set one with ",
           "Signac::Annotation(obj[['", pa, "']]) <- ..., or use ",
           "method = 'chipseeker' with an explicit `txdb` instead.")
    }
    if (isTRUE(verbose)) {
      message(sprintf("--- Annotating %d peak(s) via Signac::ClosestFeature() ---",
                      length(peaks_gr)))
    }
    closest <- Signac::ClosestFeature(obj, regions = peaks_gr)
    closest <- closest[match(seq_along(peaks_gr), seq_len(nrow(closest))), , drop = FALSE]

    peak_type <- ifelse(
      closest$distance == 0 & !is.na(closest$type) & closest$type == "exon", "exonic",
      ifelse(closest$distance == 0, "intronic",
      ifelse(closest$distance > 0 & closest$distance <= promoter_window,
             "promoter_proximal", "distal")))

    df <- data.frame(
      peak            = peak_names,
      peak_type       = peak_type,
      is_enhancer     = peak_type == "distal",
      nearest_gene    = closest$gene_name,
      distance_to_tss = closest$distance,
      stringsAsFactors = FALSE
    )

  } else {
    if (is.null(txdb)) {
      stop("`txdb` is required for method = 'chipseeker' (e.g. ",
           "TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene) ",
           "-- there is no default, since it must match your genome build ",
           "and species exactly.")
    }
    if (!requireNamespace("ChIPseeker", quietly = TRUE)) {
      stop("'ChIPseeker' is required for method = 'chipseeker'. Install with ",
           "BiocManager::install('ChIPseeker').")
    }
    if (isTRUE(verbose)) {
      message(sprintf("--- Annotating %d peak(s) via ChIPseeker::annotatePeak() ---",
                      length(peaks_gr)))
    }
    names(peaks_gr) <- peak_names
    peak_anno <- ChIPseeker::annotatePeak(
      peaks_gr, TxDb = txdb,
      tssRegion = c(-promoter_window, promoter_window),
      verbose = FALSE)
    anno_df <- as.data.frame(peak_anno)
    anno_df <- anno_df[match(peak_names, rownames(anno_df)), , drop = FALSE]

    full_anno <- anno_df$annotation
    peak_type <- ifelse(grepl("Promoter", full_anno), "promoter_proximal",
                 ifelse(grepl("Exon",     full_anno), "exonic",
                 ifelse(grepl("Intron",   full_anno), "intronic", "distal")))

    gene_col <- intersect(c("SYMBOL", "geneId"), colnames(anno_df))[1]
    df <- data.frame(
      peak             = peak_names,
      peak_type        = peak_type,
      is_enhancer      = peak_type == "distal",
      nearest_gene      = if (!is.na(gene_col)) anno_df[[gene_col]] else NA_character_,
      distance_to_tss  = anno_df$distanceToTSS,
      full_annotation  = full_anno,
      stringsAsFactors = FALSE
    )
  }
  rownames(df) <- df$peak

  obj[[pa]][["peak_type"]]       <- df$peak_type
  obj[[pa]][["is_enhancer"]]     <- df$is_enhancer
  obj[[pa]][["nearest_gene"]]    <- df$nearest_gene
  obj[[pa]][["distance_to_tss"]] <- df$distance_to_tss
  obj@misc$peak_annotation <- df

  if (isTRUE(verbose)) {
    tab <- table(df$peak_type)
    message("  ", paste(names(tab), tab, sep = "=", collapse = ", "))
  }
  obj
}
