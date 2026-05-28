#' Call doublets on a Seurat object with DoubletFinder
#'
#' Runs the standard DoubletFinder workflow on a single Seurat object:
#' normalize → PCA → choose # PCs → UMAP/clusters → pK sweep → doubletFinder.
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
calldoublet <- function(obj,
                        samplenameIndex,
                        normalization      = c("LogNormalize", "SCT"),
                        vars.to.regress    = NULL,
                        cluster_resolution = 0.1) {
  normalization <- match.arg(normalization)
  use_sct       <- identical(normalization, "SCT")

  # ---- Normalize ----------------------------------------------------------
  if (use_sct) {
    obj <- SCTransform(obj, vars.to.regress = vars.to.regress, verbose = FALSE)
  } else {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, verbose = FALSE)
    obj <- ScaleData(obj, vars.to.regress = vars.to.regress, verbose = FALSE)
  }

  obj <- RunPCA(obj, verbose = FALSE)

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
  library("irlba")
  library("RSpectra")
  obj <- RunUMAP(obj, dims = 1:min.pc, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:min.pc, verbose = FALSE)
  obj <- FindClusters(obj, resolution = cluster_resolution, verbose = FALSE)

  # ---- pK identification (no ground-truth) --------------------------------
  sweep.list  <- paramSweep(obj, PCs = 1:min.pc, sct = use_sct)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn       <- find.pK(sweep.stats)
  # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
  bcmvn.max   <- bcmvn[which.max(bcmvn$BCmetric), ]
  optimal.pk  <- bcmvn.max$pK
  optimal.pk  <- as.numeric(levels(optimal.pk))[optimal.pk]

  # ---- Homotypic doublet proportion estimate ------------------------------
  annotations    <- obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp.poi       <- round(optimal.pk * nrow(obj@meta.data))
  nExp.poi.adj   <- round(nExp.poi * (1 - homotypic.prop))

  # ---- Run DoubletFinder --------------------------------------------------
  obj <- doubletFinder(seu  = obj,
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

  # ---- Cleanup: only strip the SCT assay when we created it ---------------
  if (use_sct) {
    DefaultAssay(obj) <- 'RNA'
    obj[['SCT']]      <- NULL
  }

  return(obj)
}
