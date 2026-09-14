#' Demultiplex pooled samples from hashtag (HTO) or CMO tags
#'
#' When multiple samples are pooled into one 10x run, each is stained with a
#' distinct antibody- or lipid-tagged oligo ("hashtag") before pooling.
#' After sequencing, each cell's hashtag UMI counts identify which original
#' sample it came from. This wraps \code{Seurat::HTODemux()} (clustering-
#' based) or \code{Seurat::MULTIseqDemux()} (semi-supervised, threshold-
#' based) to make that call, then standardizes both methods' output into
#' one set of columns so downstream code doesn't need to know which was
#' used.
#'
#' \strong{Which method.} \code{"htodemux"} (default) fits a k-means model
#' per hashtag and calls a cell positive/negative for each; cells positive
#' for exactly one are singlets, more than one are doublets, none are
#' negative. \code{"multiseqdemux"} instead sweeps a quantile threshold
#' range (\code{qrange}) per hashtag and picks the one that minimizes
#' negative/doublet calls -- generally a bit more conservative about
#' calling doublets. Both ship inside Seurat; no extra package needed.
#'
#' @param obj A Seurat object with a hashtag/CMO assay already attached
#'   (e.g. via \code{CreateAssayObject(counts = hto_counts)}) and raw --
#'   \strong{not log-normalized} -- counts in it. This function CLR-
#'   normalizes it internally.
#' @param hto_assay Name of the hashtag assay. Default \code{"HTO"}.
#' @param method \code{"htodemux"} (default) or \code{"multiseqdemux"}.
#' @param positive_quantile \code{method = "htodemux"} only: passed to
#'   \code{Seurat::HTODemux()}. Default 0.99.
#' @param autoThresh,maxiter,qrange \code{method = "multiseqdemux"} only:
#'   passed to \code{Seurat::MULTIseqDemux()}. Defaults match Seurat's own.
#' @param plot Logical; if \code{TRUE} (default), also compute a ridge plot
#'   of each hashtag's signal split by the final call, stored in
#'   \code{obj@misc$hashtag_qc_plot} (a \code{ggplot}/patchwork object).
#' @param remove_negative,remove_doublets Logical; if \code{TRUE}, drop
#'   cells called Negative / Doublet after classification. Both default
#'   \code{FALSE} -- classify first, decide what to filter afterward.
#' @param verbose Message progress and call counts. Default \code{TRUE}.
#' @return \code{obj} with the method's native column(s) plus two
#'   standardized columns: \code{hash_call} (the assigned sample name, or
#'   \code{"Doublet"}/\code{"Negative"}) and \code{hash_global}
#'   (\code{"Singlet"}/\code{"Doublet"}/\code{"Negative"}), and, if
#'   \code{plot = TRUE}, \code{obj@misc$hashtag_qc_plot}.
#' @examples
#' \dontrun{
#' obj <- DemultiplexHashtags(obj, hto_assay = "HTO")
#' table(obj$hash_global)
#' obj@misc$hashtag_qc_plot
#'
#' # Keep only confidently-assigned singlets before downstream analysis
#' obj_clean <- subset(obj, hash_global == "Singlet")
#' }
#' @importFrom Seurat Assays HTODemux MULTIseqDemux NormalizeData RidgePlot
#' @export
DemultiplexHashtags <- function(obj,
                                hto_assay         = "HTO",
                                method            = c("htodemux", "multiseqdemux"),
                                positive_quantile = 0.99,
                                autoThresh        = TRUE,
                                maxiter           = 5,
                                qrange            = seq(0.1, 0.9, by = 0.05),
                                plot              = TRUE,
                                remove_negative   = FALSE,
                                remove_doublets   = FALSE,
                                verbose           = TRUE) {

  method <- match.arg(method)
  .assert_seurat(obj)
  if (!hto_assay %in% Seurat::Assays(obj)) {
    stop("Assay '", hto_assay, "' not found in `obj`. Available: ",
         paste(Seurat::Assays(obj), collapse = ", "))
  }

  obj <- Seurat::NormalizeData(obj, assay = hto_assay,
                               normalization.method = "CLR", verbose = FALSE)

  if (method == "htodemux") {
    if (isTRUE(verbose)) message(sprintf("--- Running HTODemux (assay = '%s') ---", hto_assay))
    obj <- Seurat::HTODemux(obj, assay = hto_assay,
                            positive.quantile = positive_quantile,
                            verbose = isTRUE(verbose))
    call_col   <- paste0(hto_assay, "_classification")
    global_col <- paste0(hto_assay, "_classification.global")
    obj$hash_call   <- obj@meta.data[[call_col]]
    obj$hash_global <- obj@meta.data[[global_col]]

  } else {
    if (isTRUE(verbose)) message(sprintf("--- Running MULTIseqDemux (assay = '%s') ---", hto_assay))
    obj <- Seurat::MULTIseqDemux(obj, assay = hto_assay, autoThresh = autoThresh,
                                 maxiter = maxiter, qrange = qrange,
                                 verbose = isTRUE(verbose))
    call_col <- "MULTI_ID"
    obj$hash_call <- obj@meta.data[[call_col]]
    obj$hash_global <- ifelse(
      obj$hash_call %in% c("Doublet", "Negative"),
      as.character(obj$hash_call), "Singlet")
  }

  if (isTRUE(verbose)) {
    tab <- table(obj$hash_global)
    message("  ", paste(names(tab), tab, sep = "=", collapse = ", "))
  }

  if (isTRUE(plot)) {
    p <- Seurat::RidgePlot(obj, assay = hto_assay,
                           features = rownames(obj[[hto_assay]]),
                           group.by = "hash_call", combine = TRUE)
    obj@misc$hashtag_qc_plot <- p
  }

  if (isTRUE(remove_negative) || isTRUE(remove_doublets)) {
    drop_vals <- character(0)
    if (isTRUE(remove_negative)) drop_vals <- c(drop_vals, "Negative")
    if (isTRUE(remove_doublets)) drop_vals <- c(drop_vals, "Doublet")
    keep <- !obj$hash_global %in% drop_vals
    if (isTRUE(verbose)) {
      message(sprintf("  Removing %d cell(s) called %s.",
                      sum(!keep), paste(drop_vals, collapse = "/")))
    }
    obj <- subset(obj, cells = colnames(obj)[keep])
  }

  obj
}
