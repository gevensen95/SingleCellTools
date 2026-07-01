#' Reference-based cell-type annotation with Azimuth
#'
#' Wraps \code{Azimuth::RunAzimuth} to project a query Seurat object onto
#' one of Azimuth's curated references and transfer cell-type labels.
#' Complements \code{AnnotateClusters} (marker / SingleR-based): Azimuth
#' typically gives finer, more consistent labels for common tissues
#' (PBMC, lung, kidney, adipose, ...), and is much faster than
#' scoring against curated marker lists.
#'
#' \strong{Which reference to use.} Common Azimuth references (as of
#' writing): \code{"pbmcref"}, \code{"lungref"}, \code{"kidneyref"},
#' \code{"adiposeref"}, \code{"fetusref"}, \code{"heartref"},
#' \code{"tonsilref"}, \code{"bonemarrowref"}, \code{"mousecortexref"},
#' \code{"humancortexref"}. Full list at
#' \url{https://azimuth.hubmapconsortium.org}.
#'
#' After the projection, the reference's cell-type levels are written into
#' \code{obj@meta.data} under columns prefixed with \code{prefix} (e.g.
#' \code{predicted.celltype.l1}, \code{predicted.celltype.l2}). The
#' matching prediction score columns are preserved so downstream code can
#' threshold on confidence.
#'
#' @param obj A Seurat object (query).
#' @param reference One of Azimuth's reference names; see Details. Default
#'   \code{"pbmcref"}.
#' @param assay Assay to project. Default DefaultAssay.
#' @param annotation_levels Character vector of prediction levels to keep
#'   (e.g. \code{c("l1", "l2")} for coarse and fine PBMC labels).
#'   \code{NULL} (default) keeps whatever Azimuth returns.
#' @param min_score Optional numeric threshold; cells with prediction
#'   score below this get relabelled to \code{unassigned_label}.
#'   \code{NULL} (default) keeps all Azimuth calls.
#' @param unassigned_label Label applied to low-confidence cells when
#'   \code{min_score} is set. Default \code{"Unknown"}.
#' @param verbose Passed to Azimuth. Default FALSE.
#' @return The Seurat object with new \code{predicted.*} metadata columns
#'   (and a \code{prediction.score.*} column per level), plus the
#'   \code{"ref.umap"} reduction.
#' @examples
#' \dontrun{
#' # PBMC
#' obj <- AnnotateWithReference(obj, reference = "pbmcref",
#'                              annotation_levels = c("l1", "l2"))
#' table(obj$predicted.celltype.l2)
#'
#' # Threshold low-confidence calls
#' obj <- AnnotateWithReference(obj, reference = "lungref",
#'                              min_score = 0.75)
#' }
#' @importFrom Seurat DefaultAssay
#' @export
AnnotateWithReference <- function(obj,
                                  reference          = "pbmcref",
                                  assay              = NULL,
                                  annotation_levels  = NULL,
                                  min_score          = NULL,
                                  unassigned_label   = "Unknown",
                                  verbose            = FALSE) {

  if (!requireNamespace("Azimuth", quietly = TRUE)) {
    stop("'Azimuth' is required. Install with ",
         "remotes::install_github('satijalab/azimuth').")
  }
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")

  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay

  message(sprintf("--- Running Azimuth (reference = '%s') ---", reference))
  obj <- Azimuth::RunAzimuth(
    query     = obj,
    reference = reference,
    assay     = a,
    verbose   = isTRUE(verbose)
  )

  # Which predicted.* columns exist? Filter to requested levels if any.
  pred_cols <- grep("^predicted\\.", colnames(obj@meta.data), value = TRUE)
  score_cols <- grep("^prediction\\.score\\.", colnames(obj@meta.data),
                     value = TRUE)
  if (!is.null(annotation_levels)) {
    pat <- paste0("(", paste(annotation_levels, collapse = "|"), ")$")
    keep_pred <- grep(pat, pred_cols, value = TRUE)
    if (length(keep_pred) == 0) {
      warning("Azimuth returned no predicted columns matching levels: ",
              paste(annotation_levels, collapse = ", "),
              ". Keeping all predicted columns.")
    } else {
      drop_pred <- setdiff(pred_cols, keep_pred)
      obj@meta.data[, drop_pred] <- NULL
      # Also drop the corresponding score columns
      drop_scores <- grep(
        paste0("^prediction\\.score\\.(", paste(annotation_levels,
                                                collapse = "|"), ")$"),
        score_cols, value = TRUE, invert = TRUE)
      obj@meta.data[, drop_scores] <- NULL
      pred_cols <- keep_pred
    }
  }

  # Optional confidence threshold
  if (!is.null(min_score)) {
    for (col in pred_cols) {
      score_col <- sub("^predicted\\.", "prediction.score.", col)
      if (!(score_col %in% colnames(obj@meta.data))) next
      low <- obj@meta.data[[score_col]] < min_score
      low[is.na(low)] <- TRUE
      obj@meta.data[[col]][low] <- unassigned_label
    }
    message(sprintf("  Applied min_score = %.2f to %d prediction column(s).",
                    min_score, length(pred_cols)))
  }

  message(sprintf("  Wrote %d predicted column(s): %s",
                  length(pred_cols),
                  paste(pred_cols, collapse = ", ")))
  obj
}
