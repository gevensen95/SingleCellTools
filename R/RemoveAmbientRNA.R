#' Remove ambient RNA contamination from a Seurat object's counts
#'
#' Wraps \code{decontX} (\code{celda} package) or \code{SoupX} to estimate
#' and subtract ambient/background RNA contamination, matching this
#' package's existing "known upstream tool, one call" wrappers
#' (\code{\link{RunLIANA}}, \code{\link{RunCellChat}},
#' \code{\link{RunRCTD}}, \code{\link{RunBanksyWrapper}}). Particularly
#' relevant for tissues where one dominant, fragile cell type contributes
#' outsized ambient signal -- hepatocytes in liver are a textbook case.
#'
#' \strong{Which method.} \code{method = "decontx"} (default) estimates
#' contamination directly from the filtered counts matrix already in
#' \code{obj} -- no extra input needed, at the cost of being a fully
#' data-driven estimate. \code{method = "soupx"} additionally needs the
#' raw, \emph{unfiltered} droplet matrix (every barcode SoupX should treat
#' as potential background, not just the cells Cell Ranger/your pipeline
#' kept) via \code{raw_counts}, and generally benefits more from being
#' given cluster labels (\code{clusters_col}) to estimate contamination
#' per cluster rather than globally.
#'
#' \strong{This overwrites the assay's counts layer by default.} Pass
#' \code{new_assay} to instead write the corrected counts into a new assay
#' and leave the original untouched.
#'
#' @param obj A Seurat object.
#' @param method \code{"decontx"} (default) or \code{"soupx"}. See Details.
#' @param assay Assay to correct. Default \code{DefaultAssay(obj)}.
#' @param clusters_col Metadata column of per-cell cluster/cell-type labels.
#'   Passed to decontX's \code{z} argument, or to SoupX's
#'   \code{setClusters()}. \code{NULL} (default): decontX clusters
#'   internally; SoupX runs unclustered (less accurate per its own docs).
#' @param raw_counts \code{method = "soupx"} only (required there): the
#'   raw, unfiltered droplet x gene counts matrix, or a path to a
#'   directory \code{Seurat::Read10X()} can read.
#' @param contamination_fraction \code{method = "soupx"} only: skip
#'   \code{autoEstCont()} and use this fixed contamination fraction for
#'   every cell instead. \code{NULL} (default) estimates it automatically.
#' @param new_assay If set, write corrected counts to a new assay of this
#'   name instead of overwriting \code{assay}'s counts layer in place.
#'   \code{NULL} (default) overwrites in place.
#' @param round_to_integer If \code{TRUE} (default), round corrected counts
#'   to the nearest integer -- both tools return non-integer "corrected"
#'   values, and most downstream steps (DoubletFinder, DE) expect integer
#'   counts.
#' @param verbose Message progress and the mean estimated contamination
#'   fraction. Default \code{TRUE}.
#' @return \code{obj}, with corrected counts (in place or in
#'   \code{new_assay}) and a new \code{<method>_contamination} metadata
#'   column recording each cell's estimated contamination fraction.
#' @examples
#' \dontrun{
#' # decontX -- no extra input needed
#' obj <- RemoveAmbientRNA(obj, clusters_col = "seurat_clusters")
#'
#' # SoupX -- needs the raw, unfiltered matrix too
#' obj <- RemoveAmbientRNA(obj, method = "soupx",
#'                         raw_counts = "raw_feature_bc_matrix/",
#'                         clusters_col = "seurat_clusters")
#' }
#' @importFrom Seurat DefaultAssay GetAssayData SetAssayData CreateAssayObject Read10X
#' @export
RemoveAmbientRNA <- function(obj,
                             method                  = c("decontx", "soupx"),
                             assay                   = NULL,
                             clusters_col            = NULL,
                             raw_counts              = NULL,
                             contamination_fraction  = NULL,
                             new_assay               = NULL,
                             round_to_integer        = TRUE,
                             verbose                 = TRUE) {

  method <- match.arg(method)
  .assert_seurat(obj)
  if (!is.null(clusters_col) && !clusters_col %in% colnames(obj@meta.data)) {
    stop("`clusters_col` ('", clusters_col, "') not found in obj@meta.data.")
  }
  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay
  counts <- Seurat::GetAssayData(obj, assay = a, layer = "counts")

  if (method == "decontx") {
    if (!requireNamespace("celda", quietly = TRUE)) {
      stop("'celda' is required for method = 'decontx'. Install with ",
           "BiocManager::install('celda').")
    }
    z <- if (!is.null(clusters_col)) as.character(obj@meta.data[[clusters_col]]) else NULL
    if (isTRUE(verbose)) {
      message(sprintf("--- Running decontX (assay = '%s'%s) ---", a,
                      if (!is.null(z)) ", using provided clusters" else ""))
    }
    res <- celda::decontX(as.matrix(counts), z = z)
    corrected <- res$decontXcounts
    contamination <- stats::setNames(res$contamination, colnames(counts))

  } else {
    if (!requireNamespace("SoupX", quietly = TRUE)) {
      stop("'SoupX' is required for method = 'soupx'. Install with ",
           "install.packages('SoupX').")
    }
    if (is.null(raw_counts)) {
      stop("`raw_counts` (the raw, unfiltered droplet matrix, or a path ",
           "Seurat::Read10X() can read) is required for method = 'soupx'.")
    }
    tod <- if (is.character(raw_counts)) {
      Seurat::Read10X(data.dir = raw_counts)
    } else {
      raw_counts
    }
    if (isTRUE(verbose)) message(sprintf("--- Running SoupX (assay = '%s') ---", a))
    sc <- SoupX::SoupChannel(tod, counts)
    if (!is.null(clusters_col)) {
      sc <- SoupX::setClusters(
        sc, stats::setNames(as.character(obj@meta.data[[clusters_col]]), colnames(obj)))
    }
    sc <- if (!is.null(contamination_fraction)) {
      SoupX::setContaminationFraction(sc, contamination_fraction)
    } else {
      SoupX::autoEstCont(sc, doPlot = FALSE, verbose = isTRUE(verbose))
    }
    corrected <- SoupX::adjustCounts(sc)
    contamination <- stats::setNames(sc$metaData$rho, rownames(sc$metaData))
    contamination <- contamination[colnames(counts)]
  }

  if (isTRUE(round_to_integer)) corrected <- round(corrected)

  if (isTRUE(verbose)) {
    message(sprintf("  Mean estimated contamination fraction: %.3f (range %.3f-%.3f)",
                    mean(contamination, na.rm = TRUE),
                    min(contamination, na.rm = TRUE),
                    max(contamination, na.rm = TRUE)))
  }

  if (is.null(new_assay)) {
    obj <- Seurat::SetAssayData(obj, assay = a, layer = "counts", new.data = corrected)
  } else {
    obj[[new_assay]] <- Seurat::CreateAssayObject(counts = corrected)
  }
  obj[[paste0(method, "_contamination")]] <- contamination[colnames(obj)]

  obj
}
