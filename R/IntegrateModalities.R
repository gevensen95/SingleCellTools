#' Weighted-nearest-neighbor integration of two or more modalities
#'
#' Thin wrapper around \code{Seurat::FindMultiModalNeighbors()} (WNN) for
#' objects with more than one modality already processed to a reduction
#' each -- CITE-seq (RNA + ADT) or multiome (RNA + ATAC). WNN learns, per
#' cell, how much to weight each modality's neighbor graph rather than
#' concatenating reductions or picking one modality as primary.
#'
#' Each assay must already have its own reduction computed (e.g.
#' \code{RunPCA} for RNA/ADT, \code{Signac::RunTFIDF} + \code{RunSVD} for
#' ATAC) -- this function does not compute those for you, since the right
#' preprocessing differs per modality.
#'
#' @param obj A Seurat object with >= 2 assays, each already reduced.
#' @param reduction_list Named list mapping assay name -> its reduction
#'   name (e.g. \code{list(RNA = "pca", ADT = "apca")}). \code{NULL}
#'   (default) auto-detects by matching each reduction's
#'   \code{Seurat::Assays()} slot to \code{assays}.
#' @param assays Which assays to integrate. \code{NULL} (default) uses
#'   every assay in \code{obj} (errors if that's fewer than 2).
#' @param dims_list Named list of dimensions to use per assay, matching
#'   the names of \code{reduction_list}. \code{NULL} (default) uses
#'   \code{1:30} for each, except a reduction whose name contains
#'   \code{"lsi"} (the ATAC/Signac convention), which uses \code{2:30} --
#'   the first LSI component is conventionally dropped as a
#'   sequencing-depth proxy.
#' @param k.nn Neighbors per modality. Default 20 (Seurat's own default).
#' @param reduction_name Name for the joint UMAP this also computes.
#'   Default \code{"wnn.umap"}.
#' @param graph_name_prefix Prefix for the new neighbor/graph objects
#'   (\code{<prefix>.nn}, \code{<prefix>}, \code{<prefix>.snn}). Default
#'   \code{"wsnn"}.
#' @param run_umap Logical; if \code{TRUE} (default), also run
#'   \code{Seurat::RunUMAP} on the joint neighbor graph.
#' @param verbose Message progress. Default \code{TRUE}.
#' @return \code{obj} with the joint neighbor/graph objects and, if
#'   \code{run_umap}, a \code{reduction_name} UMAP embedding.
#' @examples
#' \dontrun{
#' cite <- RunPCA(cite, npcs = 30)
#' DefaultAssay(cite) <- "ADT"
#' cite <- RunPCA(cite, npcs = 18, reduction.name = "apca")
#' cite <- IntegrateModalities(cite,
#'                             reduction_list = list(RNA = "pca", ADT = "apca"))
#' DimPlot(cite, reduction = "wnn.umap")
#' }
#' @importFrom Seurat Assays DefaultAssay FindMultiModalNeighbors RunUMAP
#' @export
IntegrateModalities <- function(obj,
                                reduction_list    = NULL,
                                assays            = NULL,
                                dims_list         = NULL,
                                k.nn              = 20,
                                reduction_name    = "wnn.umap",
                                graph_name_prefix = "wsnn",
                                run_umap          = TRUE,
                                verbose           = TRUE) {

  .assert_seurat(obj)
  if (is.null(assays)) assays <- Seurat::Assays(obj)
  if (length(assays) < 2) {
    stop("Need >= 2 assays for multimodal integration (e.g. RNA + ADT, or ",
         "RNA + ATAC); found only: ", paste(assays, collapse = ", "))
  }

  if (is.null(reduction_list)) {
    reduction_list <- list()
    for (a in assays) {
      candidates <- names(obj@reductions)[vapply(obj@reductions, function(r) {
        isTRUE(SeuratObject::DefaultAssay(r) == a)
      }, logical(1))]
      if (length(candidates) == 0) {
        stop("Could not auto-detect a reduction for assay '", a, "'. Pass ",
             "`reduction_list` explicitly, e.g. list(", a, " = \"pca\").")
      }
      reduction_list[[a]] <- candidates[1]
    }
  }

  if (is.null(dims_list)) {
    dims_list <- lapply(reduction_list, function(r) {
      if (grepl("lsi", r, ignore.case = TRUE)) 2:30 else 1:30
    })
  }

  if (isTRUE(verbose)) {
    message(sprintf("--- FindMultiModalNeighbors (%s) ---",
                    paste(names(reduction_list), unlist(reduction_list),
                         sep = "=", collapse = ", ")))
  }
  obj <- Seurat::FindMultiModalNeighbors(
    obj,
    reduction.list      = unname(reduction_list),
    dims.list           = unname(dims_list),
    k.nn                = k.nn,
    weighted.nn.name     = paste0(graph_name_prefix, ".nn"),
    knn.graph.name       = graph_name_prefix,
    snn.graph.name       = paste0(graph_name_prefix, ".snn")
  )

  if (isTRUE(run_umap)) {
    if (isTRUE(verbose)) message(sprintf("--- RunUMAP on '%s' ---", paste0(graph_name_prefix, ".nn")))
    obj <- Seurat::RunUMAP(obj, nn.name = paste0(graph_name_prefix, ".nn"),
                           reduction.name = reduction_name,
                           verbose = isTRUE(verbose))
  }

  obj
}
