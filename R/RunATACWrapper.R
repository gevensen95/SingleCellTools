#' Run Signac's standard ATAC normalization + dimensionality reduction (LSI)
#'
#' Thin wrapper around Signac's standard TF-IDF -> top-features -> SVD
#' recipe (\code{RunTFIDF()} -> \code{FindTopFeatures()} -> \code{RunSVD()}),
#' the ATAC-appropriate substitute for RNA's Normalize+PCA. Peak accessibility
#' counts are not gene-expression counts, so log-normalization/SCT and PCA
#' aren't the right tools here -- TF-IDF corrects for per-cell sequencing
#' depth and per-peak accessibility, and SVD on the TF-IDF matrix ("LSI",
#' latent semantic indexing) is Signac's drop-in replacement for PCA.
#'
#' \strong{Caveat (Signac convention)}: the first LSI component usually
#' correlates strongly with sequencing depth rather than biological variation
#' (see Signac's own \code{DepthCor()}/vignette guidance). Callers typically
#' drop it downstream (e.g. \code{FindNeighbors(dims = 2:30)}) rather than
#' this wrapper dropping it itself, since whether to drop it -- and which
#' component(s) -- is a per-dataset judgment call best left to
#' \code{DepthCor()} inspection. \code{\link{MergeSeurat}} makes this a
#' first-class parameter (\code{atac_lsi_first_dim}, default \code{2}) for
#' exactly this reason.
#'
#' @param obj A Seurat object with a \code{ChromatinAssay} (from Signac).
#' @param assay Assay to run TF-IDF/LSI on. Defaults to
#'   \code{DefaultAssay(obj)}.
#' @param min_cutoff Passed to \code{Signac::FindTopFeatures()}; the minimum
#'   fraction of cells a peak must be detected in to be retained for LSI.
#'   Default \code{"q0"} (Signac's own default -- keep all peaks).
#' @param n_components Number of SVD components to compute. Default
#'   \code{30}. Passed to \code{Signac::RunSVD()} as \code{n}.
#' @param reduction_name Name for the resulting reduction. Default
#'   \code{"lsi"}.
#' @param verbose Passed through to \code{RunTFIDF()}/\code{RunSVD()}.
#' @return The Seurat object with a new \code{"lsi"} (or \code{reduction_name})
#'   reduction, computed on TF-IDF-normalized data.
#' @examples
#' \dontrun{
#' obj <- RunATACWrapper(obj, n_components = 30)
#' Seurat::DepthCor(obj)  # check whether component 1 tracks sequencing depth
#' obj <- Seurat::FindNeighbors(obj, reduction = "lsi", dims = 2:30)
#' }
#' @importFrom Seurat DefaultAssay
#' @export
RunATACWrapper <- function(obj,
                           assay          = NULL,
                           min_cutoff     = 'q0',
                           n_components   = 30,
                           reduction_name = 'lsi',
                           verbose        = TRUE) {

  # Cheap, local checks first so they surface their own errors even on a
  # machine without Signac installed, rather than being masked by the
  # package-availability error below.
  if (!inherits(obj, 'Seurat')) stop('`obj` must be a Seurat object.')
  if (!requireNamespace('Signac', quietly = TRUE)) {
    stop("'Signac' is required. Install with: ",
         "BiocManager::install('Signac')")
  }

  if (is.null(assay)) assay <- Seurat::DefaultAssay(obj)

  if (!inherits(obj[[assay]], 'ChromatinAssay')) {
    stop("Assay '", assay, "' is not a ChromatinAssay. RunATACWrapper() ",
         "expects ATAC (Signac) peak-by-cell data -- pass `assay` explicitly ",
         "if the ChromatinAssay isn't the default assay.")
  }

  counts_mat <- tryCatch(Seurat::GetAssayData(obj, assay = assay, layer = 'counts'),
                         error = function(e) NULL)
  if (is.null(counts_mat) || nrow(counts_mat) == 0 || ncol(counts_mat) == 0) {
    stop("The 'counts' layer of assay '", assay, "' is empty (",
         if (is.null(counts_mat)) '0 x 0' else paste(nrow(counts_mat), 'x', ncol(counts_mat)),
         "). Nothing to run TF-IDF/LSI on.")
  }

  # ---- TF-IDF ---------------------------------------------------------------
  message(sprintf('--- Running TF-IDF (assay = %s) ---', assay))
  obj <- Signac::RunTFIDF(obj, assay = assay, verbose = verbose)

  # ---- Top features ----------------------------------------------------------
  message(sprintf('--- Finding top features (min.cutoff = %s) ---', min_cutoff))
  obj <- Signac::FindTopFeatures(obj, assay = assay, min.cutoff = min_cutoff)

  # ---- SVD (LSI) --------------------------------------------------------------
  message(sprintf('--- Running SVD (n = %d, reduction.name = %s) ---',
                  n_components, reduction_name))
  obj <- Signac::RunSVD(obj, assay = assay, n = n_components,
                        reduction.name = reduction_name, verbose = verbose)

  obj
}
