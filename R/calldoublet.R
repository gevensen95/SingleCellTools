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
  # For SCT, drop the entire SCT assay — all its derived layers go with it.
  if (use_sct && "SCT" %in% SeuratObject::Assays(obj)) {
    SeuratObject::DefaultAssay(obj) <- "RNA"
    obj[["SCT"]]                    <- NULL
  }

  # Drop normalized layers from the RNA assay, including any split variants
  # (data.1, data.2, scale.data.1, ...) that show up after a previous merge.
  rna <- obj[["RNA"]]
  target_patterns <- c("^data(\\.|$)", "^scale\\.data(\\.|$)")
  lyrs_to_drop <- unique(unlist(lapply(
    target_patterns,
    function(p) grep(p, names(rna@layers), value = TRUE)
  )))
  for (lyr in lyrs_to_drop) {
    rna@layers[[lyr]] <- NULL
  }

  # Clear variable features. Assay5 stores them as vf_* columns in @meta.data;
  # the legacy Assay class uses the @var.features slot. Handle both.
  if (inherits(rna, "Assay5")) {
    vf_cols <- grep("^vf_", colnames(rna@meta.data), value = TRUE)
    if (length(vf_cols)) rna@meta.data[vf_cols] <- list(NULL)
  } else {
    rna@var.features <- character(0)
  }
  obj[["RNA"]] <- rna

  # Drop pca / umap reductions, allowing for assay-suffixed names like
  # "pca.RNA" or "umap.SCT".
  red_to_drop <- grep("^(pca|umap)(\\.|$)",
                      names(obj@reductions), value = TRUE)
  for (red in red_to_drop) {
    obj@reductions[[red]] <- NULL
  }

  return(obj)
}
