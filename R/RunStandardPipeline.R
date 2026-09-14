#' An opinionated, one-call standard processing pipeline
#'
#' Chains normalization, feature selection, scaling, PCA, optional batch
#' integration, neighbor graph, clustering, and UMAP into a single call
#' with sensible defaults -- the sequence a new user otherwise has to know
#' to run by hand across ~6 separate function calls, in the right order.
#' Intended as a fast path to "a clustered, visualizable object"; for
#' anything beyond the defaults (custom feature selection, SCTransform,
#' non-linear integration), call the underlying Seurat functions directly
#' -- this function is deliberately not configurable beyond what's exposed
#' here.
#'
#' @param obj A Seurat object with raw counts.
#' @param normalize Logical; run \code{Seurat::NormalizeData()}. Default
#'   \code{TRUE} -- set \code{FALSE} if already normalized.
#' @param n_variable_features Genes kept by \code{FindVariableFeatures()}.
#'   Default 2000.
#' @param n_pcs Principal components computed and used downstream. Default 30.
#' @param integration_method \code{"none"} (default), \code{"harmony"}, or
#'   \code{"cca"} (via \code{Seurat::IntegrateLayers}). Requires
#'   \code{batch_col} when not \code{"none"}.
#' @param batch_col Metadata column identifying the batch/sample to
#'   integrate over. Required if \code{integration_method != "none"}.
#' @param cluster_resolution Passed to \code{Seurat::FindClusters()}.
#'   Default 0.8.
#' @param run_umap Logical; run \code{Seurat::RunUMAP()}. Default \code{TRUE}.
#' @param assay Assay to process. Default \code{DefaultAssay(obj)}.
#' @param save_provenance Optional path; if set, calls
#'   \code{\link{SaveWithProvenance}} on the final object with the pipeline
#'   parameters recorded under \code{extra}.
#' @param verbose Message progress per step. Default \code{TRUE}.
#' @return The processed Seurat object (PCA -- and, if requested, harmony/
#'   integrated.cca and UMAP -- reductions; \code{seurat_clusters} in
#'   metadata).
#' @examples
#' \dontrun{
#' obj <- RunStandardPipeline(obj)
#' DimPlot(obj, group.by = "seurat_clusters")
#'
#' # With batch integration + a provenance sidecar
#' obj <- RunStandardPipeline(obj, integration_method = "harmony",
#'                            batch_col = "orig.ident",
#'                            save_provenance = "processed.rds")
#' }
#' @importFrom Seurat DefaultAssay FindClusters FindNeighbors FindVariableFeatures IntegrateLayers NormalizeData RunPCA RunUMAP ScaleData
#' @export
RunStandardPipeline <- function(obj,
                                normalize            = TRUE,
                                n_variable_features  = 2000,
                                n_pcs                = 30,
                                integration_method   = c("none", "harmony", "cca"),
                                batch_col            = NULL,
                                cluster_resolution   = 0.8,
                                run_umap             = TRUE,
                                assay                = NULL,
                                save_provenance      = NULL,
                                verbose              = TRUE) {

  integration_method <- match.arg(integration_method)
  .assert_seurat(obj)
  if (integration_method != "none" &&
      (is.null(batch_col) || !batch_col %in% colnames(obj@meta.data))) {
    stop("`batch_col` must name a column of obj@meta.data when ",
         "integration_method != 'none'.")
  }
  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay

  if (isTRUE(normalize)) {
    if (isTRUE(verbose)) message("--- NormalizeData ---")
    obj <- Seurat::NormalizeData(obj, assay = a, verbose = FALSE)
  }
  if (isTRUE(verbose)) message(sprintf("--- FindVariableFeatures (%d) ---", n_variable_features))
  obj <- Seurat::FindVariableFeatures(obj, assay = a, nfeatures = n_variable_features,
                                      verbose = FALSE)
  if (isTRUE(verbose)) message("--- ScaleData ---")
  obj <- Seurat::ScaleData(obj, assay = a, verbose = FALSE)
  if (isTRUE(verbose)) message(sprintf("--- RunPCA (%d pcs) ---", n_pcs))
  obj <- Seurat::RunPCA(obj, assay = a, npcs = n_pcs, verbose = FALSE)

  reduction <- "pca"
  if (integration_method == "harmony") {
    if (!requireNamespace("harmony", quietly = TRUE)) {
      stop("'harmony' is required for integration_method = 'harmony'. ",
           "Install with install.packages('harmony').")
    }
    if (isTRUE(verbose)) message(sprintf("--- RunHarmony (batch = '%s') ---", batch_col))
    obj <- harmony::RunHarmony(obj, group.by.vars = batch_col, verbose = isTRUE(verbose))
    reduction <- "harmony"
  } else if (integration_method == "cca") {
    if (isTRUE(verbose)) message(sprintf("--- IntegrateLayers (CCA, batch = '%s') ---", batch_col))
    obj[[a]] <- split(obj[[a]], f = obj@meta.data[[batch_col]])
    obj <- Seurat::IntegrateLayers(obj, method = Seurat::CCAIntegration,
                                   orig.reduction = "pca",
                                   new.reduction = "integrated.cca",
                                   verbose = isTRUE(verbose))
    obj[[a]] <- SeuratObject::JoinLayers(obj[[a]])
    reduction <- "integrated.cca"
  }

  if (isTRUE(verbose)) message(sprintf("--- FindNeighbors/FindClusters (reduction = '%s') ---", reduction))
  obj <- Seurat::FindNeighbors(obj, reduction = reduction, dims = 1:n_pcs, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = cluster_resolution, verbose = FALSE)

  if (isTRUE(run_umap)) {
    if (isTRUE(verbose)) message("--- RunUMAP ---")
    obj <- Seurat::RunUMAP(obj, reduction = reduction, dims = 1:n_pcs, verbose = FALSE)
  }

  if (!is.null(save_provenance)) {
    if (isTRUE(verbose)) message(sprintf("--- Saving with provenance: %s ---", save_provenance))
    SaveWithProvenance(obj, save_provenance, extra = list(
      pipeline = "RunStandardPipeline",
      params = list(n_variable_features = n_variable_features, n_pcs = n_pcs,
                    integration_method = integration_method, batch_col = batch_col,
                    cluster_resolution = cluster_resolution)
    ))
  }

  obj
}
