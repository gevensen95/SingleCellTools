#' Call doublets on a Seurat object with DoubletFinder
#'
#' Runs the standard DoubletFinder workflow on a single Seurat object:
#' normalize → PCA → choose # PCs → UMAP/clusters → pK sweep → doubletFinder.
#'
#' On return, the normalized \code{data} and \code{scale.data} layers, the
#' variable features, and the \code{pca} / \code{umap} reductions created
#' during the workflow are stripped from the object so the result carries
#' only counts plus the new \code{doublet_finder} (and \code{seurat_clusters})
#' metadata columns.
#'
#' @param obj A Seurat object.
#' @param samplenameIndex Passed through unchanged from the original API
#'   (currently unused inside the function; retained for call-site
#'   compatibility).
#' @param normalization One of \code{"LogNormalize"} (default) or \code{"SCT"}.
#'   Selects the \code{NormalizeData / FindVariableFeatures / ScaleData}
#'   pipeline vs \code{SCTransform} AND flips the \code{sct} flag passed to
#'   DoubletFinder's \code{paramSweep} and \code{doubletFinder}.
#' @param vars.to.regress Optional character vector of variables to regress
#'   out during normalization (passed to either \code{SCTransform} or
#'   \code{ScaleData}). Default \code{NULL}.
#' @param cluster_resolution Resolution passed to \code{FindClusters}.
#'   Default \code{0.1}.
#' @return The input Seurat object with a \code{doublet_finder} metadata
#'   column ("Doublet" / "Singlet").
#' @importFrom Seurat SCTransform NormalizeData FindVariableFeatures ScaleData RunPCA RunUMAP FindNeighbors FindClusters
#' @importFrom SeuratObject DefaultAssay DefaultAssay<- Assays Layers Reductions VariableFeatures VariableFeatures<-
#' @importFrom DoubletFinder paramSweep summarizeSweep find.pK modelHomotypic doubletFinder
#' @export
calldoublet <- function(obj,
                        samplenameIndex,
                        normalization      = c("LogNormalize", "SCT"),
                        vars.to.regress    = NULL,
                        cluster_resolution = 0.1) {
  normalization <- match.arg(normalization)
  use_sct       <- identical(normalization, "SCT")

  # ---- Normalize ----------------------------------------------------------
  if (use_sct) {
    obj <- Seurat::SCTransform(obj, vars.to.regress = vars.to.regress,
                               verbose = FALSE)
  } else {
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, verbose = FALSE)
    obj <- Seurat::ScaleData(obj, vars.to.regress = vars.to.regress,
                             verbose = FALSE)
  }

  obj <- Seurat::RunPCA(obj, verbose = FALSE)

  # ---- Find significant PCs -----------------------------------------------
  stdv         <- obj[["pca"]]@stdev
  percent.stdv <- (stdv / sum(stdv)) * 100
  cumulative   <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:(length(percent.stdv) - 1)] -
                       percent.stdv[2:length(percent.stdv)]) > 0.1),
              decreasing = TRUE)[1] + 1
  min.pc <- min(co1, co2)

  # ---- Finish pre-processing ----------------------------------------------
  # NB: irlba and RSpectra are pulled in via the package Imports so that the
  # PCA/UMAP code that needs them can find them — no library() calls here.
  obj <- Seurat::RunUMAP(obj, dims = 1:min.pc, verbose = FALSE)
  obj <- Seurat::FindNeighbors(obj, dims = 1:min.pc, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = cluster_resolution, verbose = FALSE)

  # ---- pK identification (no ground-truth) --------------------------------
  sweep.list  <- DoubletFinder::paramSweep(obj, PCs = 1:min.pc, sct = use_sct)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.list)
  bcmvn       <- DoubletFinder::find.pK(sweep.stats)
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  bcmvn.max   <- bcmvn[which.max(bcmvn$BCmetric), ]
  optimal.pk  <- bcmvn.max$pK
  optimal.pk  <- as.numeric(levels(optimal.pk))[optimal.pk]

  # ---- Homotypic doublet proportion estimate ------------------------------
  annotations    <- obj@meta.data$seurat_clusters
  homotypic.prop <- DoubletFinder::modelHomotypic(annotations)
  nExp.poi       <- round(optimal.pk * nrow(obj@meta.data))
  nExp.poi.adj   <- round(nExp.poi * (1 - homotypic.prop))

  # ---- Run DoubletFinder --------------------------------------------------
  obj <- DoubletFinder::doubletFinder(seu  = obj,
                                      PCs  = 1:min.pc,
                                      pK   = optimal.pk,
                                      nExp = nExp.poi.adj,
                                      sct  = use_sct)

  metadata <- obj@meta.data
  colnames(metadata)[length(colnames(metadata))] <- "doublet_finder"
  obj@meta.data <- metadata

  cat("Doublet estimation with DoubletFinder() ...\n")
  cat("Normalization: ", normalization, "\n")
  cat("Doublet: ", as.vector(table(metadata[, length(colnames(metadata))]))[1], "\n")
  cat("Singlet: ", as.vector(table(metadata[, length(colnames(metadata))]))[2], "\n")

  # ---- Cleanup: strip everything we created during the workflow -----------
  # For SCT, drop the entire SCT assay — all of its derived layers go with it.
  if (use_sct && "SCT" %in% SeuratObject::Assays(obj)) {
    SeuratObject::DefaultAssay(obj) <- 'RNA'
    obj[['SCT']]                    <- NULL
  }

  # Drop the normalized 'data' and 'scale.data' layers from the RNA assay.
  # No-op for any layer that isn't present.
  for (lyr in c("data", "scale.data")) {
    if (lyr %in% SeuratObject::Layers(obj[["RNA"]])) {
      obj[["RNA"]][[lyr]] <- NULL
    }
  }

  # Clear the variable-features set chosen by FindVariableFeatures /
  # SCTransform. tryCatch keeps this safe across Seurat versions where the
  # setter signature has varied.
  tryCatch({
    SeuratObject::VariableFeatures(obj, assay = "RNA") <- character(0)
  }, error = function(e) invisible(NULL))

  # Drop the PCA and UMAP reductions we computed for the pK sweep.
  for (red in intersect(c("pca", "umap"), SeuratObject::Reductions(obj))) {
    obj[[red]] <- NULL
  }

  return(obj)
}
