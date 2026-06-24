#' Subset and re-cluster a Seurat object
#'
#' Subsets cells (by identity, metadata column value, or explicit cell list),
#' cleans up the resulting object, then re-runs the rest of the standard
#' pipeline on the subset: (optional) normalize, PCA, integrate (optional),
#' UMAP, cluster. Useful for cell-type-specific re-analyses where the
#' first-pass clustering is too coarse.
#'
#' \strong{Normalization (\code{normalize}).} \code{RunPCA()} requires
#' \code{data} and \code{scale.data} layers; if the working assay has only
#' a \code{counts} layer, PCA errors with
#' \code{'arg' should be "counts"}. To avoid that, the default
#' \code{normalize = "auto"} inspects the layers in the working assay and:
#' \itemize{
#'   \item if \code{data} and \code{scale.data} layers already exist, does
#'     nothing;
#'   \item if either is missing, runs \code{NormalizeData()} +
#'     \code{FindVariableFeatures()} + \code{ScaleData()} (or
#'     \code{SCTransform()} when \code{normalization_method = "SCT"}).
#' }
#' Set \code{normalize = "none"} to keep the historical behavior (assume
#' the caller already normalized the input) or \code{"LogNormalize"} /
#' \code{"SCT"} to force a re-normalization regardless of what layers
#' already exist.
#'
#' \strong{Layers are not split or re-split here.} If \code{obj}'s assay has
#' per-sample split layers (e.g. \code{data.1}, \code{data.2}, ... from the
#' standard Seurat v5 integration workflow), \code{subset()} preserves that
#' split structure, and Seurat's normalization functions run per-layer
#' automatically. If you need split layers for integration and \code{obj}
#' doesn't already have them, split it (e.g.
#' \code{obj[[assay]] <- split(obj[[assay]], f = obj$orig.ident)}) before
#' calling this function.
#'
#' Exactly one of \code{idents}, \code{metadata_col}+\code{metadata_value},
#' or \code{cells} must be supplied.
#'
#' After subsetting, the object is cleaned up before re-clustering:
#' \enumerate{
#'   \item \strong{Assay selection.} The working assay (\code{a}) is chosen
#'     via \code{assay}: if supplied, it is used directly (and must exist on
#'     \code{obj}); otherwise \code{"RNA"} is used if present, falling back
#'     to \code{"Spatial"}, and finally to \code{DefaultAssay(obj)}. The
#'     chosen assay becomes the object's default assay.
#'   \item \strong{Join \code{counts} layers (for QC only).} If the working
#'     assay has more than one \code{counts} layer (e.g. a Seurat v5 object
#'     split by sample), \code{SeuratObject::JoinLayers()} is called on just
#'     the \code{counts} layers so they can be read as a single matrix for
#'     the QC steps below.
#'   \item \strong{Recalculate QC metadata.} \code{nCount_<assay>} and
#'     \code{nFeature_<assay>} are recomputed from the subsetted counts.
#'   \item \strong{Drop empty cells.} Any cell with \code{nCount_<assay> == 0}
#'     after subsetting is removed.
#'   \item \strong{Drop empty features.} Any gene detected in fewer than
#'     \code{min_cells_per_gene} cells across all retained cells is removed
#'     (prevents downstream SCT \emph{subscript out of bounds} errors).
#' }
#'
#' @param obj A Seurat object.
#' @param idents Identity labels (active idents) to keep.
#' @param metadata_col Metadata column to filter on.
#' @param metadata_value Value(s) of \code{metadata_col} to keep.
#' @param cells Explicit character vector of cell barcodes to keep.
#' @param assay Assay to use for QC recalculation and as the object's default
#'   assay going forward. \code{NULL} (default) auto-selects: \code{"RNA"}
#'   if present, else \code{"Spatial"}, else \code{DefaultAssay(obj)}.
#' @param min_cells_per_gene After subsetting, drop genes detected (count > 0)
#'   in fewer than this many cells. Default 3.
#' @param normalize One of \code{"auto"} (default), \code{"none"},
#'   \code{"LogNormalize"}, \code{"SCT"}. See Details.
#' @param integrate Logical; if TRUE, run \code{IntegrateLayers}. Default TRUE.
#' @param integration_method Integration method (e.g. \code{"HarmonyIntegration"}).
#' @param normalization_method Normalization method that was \emph{already}
#'   applied to \code{obj} (\code{"LogNormalize"} or \code{"SCT"}), passed
#'   through to \code{IntegrateLayers(normalization.method = ...)} when
#'   \code{integrate = TRUE}. Default \code{"LogNormalize"}. Also drives the
#'   choice between LogNormalize and SCT in \code{normalize = "auto"} mode.
#' @param new_reduction Name of the integrated reduction (default
#'   \code{"harmony"}).
#' @param dims Number of PCs to use for neighbors/UMAP. Default 15.
#' @param resolution Cluster resolution. Default 0.3.
#' @return The subsetted, re-clustered Seurat object.
#' @importFrom Seurat Idents RunPCA RunUMAP FindNeighbors FindClusters IntegrateLayers NormalizeData FindVariableFeatures ScaleData SCTransform
#' @importFrom SeuratObject DefaultAssay Layers JoinLayers LayerData
#' @export
SubsetAndRecluster <- function(obj,
                               idents               = NULL,
                               metadata_col         = NULL,
                               metadata_value       = NULL,
                               cells                = NULL,
                               assay                = NULL,
                               min_cells_per_gene   = 3,
                               normalize            = c("auto", "none",
                                                        "LogNormalize", "SCT"),
                               integrate            = TRUE,
                               integration_method   = "HarmonyIntegration",
                               normalization_method = c("LogNormalize", "SCT"),
                               new_reduction        = "harmony",
                               dims                 = 15,
                               resolution           = 0.3) {

  normalize            <- match.arg(normalize)
  normalization_method <- match.arg(normalization_method)

  # ---- Validate that exactly one subset spec is given ---------------------
  modes <- c(!is.null(idents),
             !is.null(metadata_col) || !is.null(metadata_value),
             !is.null(cells))
  if (sum(modes) != 1) {
    stop("Provide exactly one of: `idents`, `metadata_col`+`metadata_value`, ",
         "or `cells`.")
  }

  # ---- Subset --------------------------------------------------------------
  message("--- Subsetting ---")
  if (!is.null(idents)) {
    obj <- subset(obj, idents = idents)
  } else if (!is.null(metadata_col)) {
    if (is.null(metadata_value)) stop("`metadata_value` required.")
    keep <- rownames(obj@meta.data)[obj@meta.data[[metadata_col]] %in% metadata_value]
    if (length(keep) == 0) stop("No cells match the metadata filter.")
    obj <- subset(obj, cells = keep)
  } else {
    obj <- subset(obj, cells = cells)
  }
  message(sprintf("  %d cells retained", ncol(obj)))

  # ---- Choose working assay -------------------------------------------------
  if (is.null(assay)) {
    if ("RNA" %in% names(obj@assays)) {
      a <- "RNA"
    } else if ("Spatial" %in% names(obj@assays)) {
      a <- "Spatial"
    } else {
      a <- SeuratObject::DefaultAssay(obj)
    }
  } else {
    if (!assay %in% names(obj@assays)) {
      stop("Assay '", assay, "' not found in object. Available assays: ",
           paste(names(obj@assays), collapse = ", "))
    }
    a <- assay
  }
  message(sprintf("--- Using assay '%s' ---", a))
  SeuratObject::DefaultAssay(obj) <- a

  # ---- Join 'counts' layers if the working assay was split (Seurat v5) ----
  counts_layers <- grep("^counts", SeuratObject::Layers(obj[[a]]), value = TRUE)
  if (length(counts_layers) > 1) {
    message(sprintf("--- Joining %d 'counts' layer(s) in assay '%s' ---",
                    length(counts_layers), a))
    obj <- tryCatch(
      SeuratObject::JoinLayers(obj, assay = a, layers = counts_layers),
      error = function(e) {
        warning("JoinLayers() failed (", conditionMessage(e),
                "); continuing with separate layers.")
        obj
      }
    )
  }

  # ---- Recalculate nCount_<assay> / nFeature_<assay> -----------------------
  message(sprintf("--- Recalculating nCount_%s / nFeature_%s ---", a, a))
  counts_layer_names <- grep("^counts", SeuratObject::Layers(obj[[a]]), value = TRUE)
  if (length(counts_layer_names) == 0) {
    stop("Assay '", a, "' has no 'counts' layer(s); cannot recalculate ",
         "nCount/nFeature.")
  }
  all_cells <- colnames(obj)
  ncount   <- stats::setNames(numeric(length(all_cells)), all_cells)
  nfeature <- stats::setNames(numeric(length(all_cells)), all_cells)
  for (ly in counts_layer_names) {
    m <- SeuratObject::LayerData(obj, assay = a, layer = ly)
    if (is.null(m) || ncol(m) == 0) next
    cs <- Matrix::colSums(m)
    fs <- Matrix::colSums(m > 0)
    ncount[names(cs)]   <- ncount[names(cs)]   + cs
    nfeature[names(fs)] <- nfeature[names(fs)] + fs
  }
  obj[[paste0("nCount_",   a)]] <- unname(ncount[all_cells])
  obj[[paste0("nFeature_", a)]] <- unname(nfeature[all_cells])

  # ---- Drop cells with 0 total counts --------------------------------------
  zero_cells <- names(ncount)[ncount == 0]
  if (length(zero_cells) > 0) {
    message(sprintf("--- Removing %d cell(s) with 0 total counts in assay '%s' ---",
                    length(zero_cells), a))
    obj <- subset(obj, cells = setdiff(colnames(obj), zero_cells))
  }

  # ---- Drop sparsely-detected features (genes) -----------------------------
  counts_layer_names <- grep("^counts", SeuratObject::Layers(obj[[a]]), value = TRUE)
  all_genes    <- rownames(obj[[a]])
  gene_ncells  <- stats::setNames(numeric(length(all_genes)), all_genes)
  for (ly in counts_layer_names) {
    m <- SeuratObject::LayerData(obj, assay = a, layer = ly)
    if (is.null(m) || ncol(m) == 0) next
    fs <- Matrix::rowSums(m > 0)
    gene_ncells[names(fs)] <- gene_ncells[names(fs)] + fs
  }
  drop_genes <- names(gene_ncells)[gene_ncells < min_cells_per_gene]
  if (length(drop_genes) > 0 && length(drop_genes) < length(all_genes)) {
    message(sprintf(
      "--- Removing %d feature(s) detected in < %d cell(s) in assay '%s' ---",
      length(drop_genes), min_cells_per_gene, a))
    obj <- subset(obj, features = setdiff(all_genes, drop_genes))
  }

  # ---- Normalize (auto / forced) ------------------------------------------
  # RunPCA() requires 'data' and 'scale.data' layers. When only 'counts'
  # exists (the original object was never normalized, or pruning above
  # collapsed it), PCA errors with `'arg' should be "counts"`. This block
  # rebuilds them.
  current_layers <- SeuratObject::Layers(obj[[a]])
  has_data       <- any(grepl("^data",       current_layers))
  has_scaledata  <- any(grepl("^scale.data", current_layers))

  do_norm <- switch(
    normalize,
    auto         = !(has_data && has_scaledata),
    none         = FALSE,
    LogNormalize = TRUE,
    SCT          = TRUE
  )

  if (do_norm) {
    # Explicit normalize value wins; otherwise in auto mode follow
    # normalization_method (so SCT users get SCT re-normalization).
    method <- if (normalize %in% c("LogNormalize", "SCT")) normalize
              else normalization_method

    if (method == "SCT") {
      message("--- Normalizing (SCTransform) ---")
      obj <- Seurat::SCTransform(obj, assay = a, verbose = FALSE)
      # SCTransform creates / sets a new assay called "SCT" and makes it
      # the default; downstream PCA / integration should use it.
      a <- SeuratObject::DefaultAssay(obj)
      message(sprintf("  (working assay is now '%s')", a))
    } else {
      message("--- Normalizing (LogNormalize + FindVariableFeatures + ScaleData) ---")
      obj <- Seurat::NormalizeData(obj, assay = a, verbose = FALSE)
      obj <- Seurat::FindVariableFeatures(obj, assay = a, verbose = FALSE)
      obj <- Seurat::ScaleData(obj, assay = a, verbose = FALSE)
    }
  } else if (normalize == "auto") {
    message("--- Skipping normalization (data + scale.data already present) ---")
  }

  # ---- PCA -----------------------------------------------------------------
  message("--- PCA ---")
  obj <- Seurat::RunPCA(obj, assay = a, verbose = FALSE)

  # ---- Integrate (optional) ------------------------------------------------
  reduction_for_downstream <- "pca"
  if (isTRUE(integrate)) {
    message(sprintf("--- Integrating (%s) ---", integration_method))
    obj <- Seurat::IntegrateLayers(
      object               = obj,
      method               = integration_method,
      orig.reduction       = "pca",
      normalization.method = normalization_method,
      new.reduction        = new_reduction,
      verbose              = FALSE
    )
    reduction_for_downstream <- new_reduction
  }

  # ---- Cluster + UMAP ------------------------------------------------------
  message("--- Neighbors / clusters / UMAP ---")
  obj <- Seurat::FindNeighbors(obj, reduction = reduction_for_downstream,
                               dims = seq_len(dims), verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = resolution, verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, reduction = reduction_for_downstream,
                         dims = seq_len(dims), verbose = FALSE)

  obj
}
