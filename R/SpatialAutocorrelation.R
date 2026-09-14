#' Spatially variable genes via Moran's I
#'
#' Cleaner interface to \code{Seurat::FindSpatiallyVariableFeatures()}: runs
#' it, pulls the resulting per-gene statistics into a tidy ranked table
#' (rather than leaving them buried in \code{obj[[assay]][[]]}), and plots
#' the top hits with \code{\link{SpatialFeaturePlotFixed}}.
#'
#' @param obj A Visium (or other spatial) Seurat object.
#' @param assay Assay to test. Default \code{DefaultAssay(obj)}.
#' @param features Genes to test. \code{NULL} (default) uses
#'   \code{Seurat::VariableFeatures(obj)} if set, otherwise every gene in
#'   \code{assay} restricted to the top \code{top_n} by mean expression
#'   (testing the full transcriptome is slow and rarely necessary).
#' @param image_name Image to plot on for the QC plot. \code{NULL}
#'   (default) uses the first image on \code{obj}.
#' @param selection.method Passed to
#'   \code{Seurat::FindSpatiallyVariableFeatures()}. Default
#'   \code{"moransi"} (fast); \code{"markvariogram"} is the slower
#'   alternative Seurat also supports.
#' @param top_n Cap on genes tested when \code{features} isn't supplied.
#'   Default 2000.
#' @param plot Logical; if \code{TRUE} (default), also return a spatial
#'   feature plot of the top \code{top_n_plot} genes.
#' @param top_n_plot Number of top genes to plot. Default 9.
#' @param verbose Message progress. Default \code{TRUE}.
#' @return A list with \code{results} (data frame: \code{gene} plus
#'   whatever Moran's-I/mark-variogram statistic columns Seurat wrote --
#'   column names vary slightly by Seurat version, so every matching
#'   column is kept rather than one hard-coded name), \code{ranked_features}
#'   (character vector, most-to-least spatially variable), and, if
#'   \code{plot = TRUE}, \code{plot}.
#' @examples
#' \dontrun{
#' svg <- SpatialAutocorrelation(visium)
#' head(svg$results)
#' svg$plot
#' }
#' @importFrom Seurat DefaultAssay FindSpatiallyVariableFeatures SpatiallyVariableFeatures VariableFeatures
#' @export
SpatialAutocorrelation <- function(obj,
                                   assay            = NULL,
                                   features         = NULL,
                                   image_name       = NULL,
                                   selection.method = "moransi",
                                   top_n            = 2000,
                                   plot             = TRUE,
                                   top_n_plot       = 9,
                                   verbose          = TRUE) {

  .assert_seurat(obj)
  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay

  feats <- features
  if (is.null(feats)) {
    feats <- tryCatch(Seurat::VariableFeatures(obj), error = function(e) character(0))
    if (length(feats) == 0) {
      means <- Matrix::rowMeans(Seurat::GetAssayData(obj, assay = a, layer = "data"))
      feats <- names(sort(means, decreasing = TRUE))[seq_len(min(top_n, length(means)))]
    } else {
      feats <- feats[seq_len(min(top_n, length(feats)))]
    }
  }

  if (isTRUE(verbose)) {
    message(sprintf("--- FindSpatiallyVariableFeatures (%s, %d genes) ---",
                    selection.method, length(feats)))
  }
  obj <- Seurat::FindSpatiallyVariableFeatures(
    obj, assay = a, features = feats, selection.method = selection.method
  )
  ranked <- Seurat::SpatiallyVariableFeatures(obj, assay = a, selection.method = selection.method)

  stat_cols <- grep("moran|markvariogram|spatial", colnames(obj[[a]][[]]),
                    ignore.case = TRUE, value = TRUE)
  results <- obj[[a]][[]][ranked, stat_cols, drop = FALSE]
  results$gene <- ranked
  results <- results[, c("gene", stat_cols)]
  rownames(results) <- NULL

  out <- list(results = results, ranked_features = ranked)

  if (isTRUE(plot)) {
    top_feats <- ranked[seq_len(min(top_n_plot, length(ranked)))]
    if (isTRUE(verbose)) {
      message(sprintf("  Top spatially variable gene(s): %s",
                      paste(head(top_feats, 5), collapse = ", ")))
    }
    out$plot <- SpatialFeaturePlotFixed(obj, features = top_feats, image = image_name)
  }

  out
}
