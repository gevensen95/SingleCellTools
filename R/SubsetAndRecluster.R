#' Subset and re-cluster a Seurat object
#'
#' Subsets cells (by identity, metadata column value, or explicit cell list),
#' cleans up the resulting object, then re-runs the rest of the standard
#' pipeline on the subset: PCA, integrate (optional), UMAP, cluster. Useful
#' for cell-type-specific re-analyses where the first-pass clustering is too
#' coarse.
#'
#' \strong{Layers are not split or re-split here.} If \code{obj}'s assay has
#' per-sample split layers (e.g. \code{data.1}, \code{data.2}, ... from the
#' standard Seurat v5 integration workflow), \code{subset()} preserves that
#' split structure on the retained cells/features, and it is left as-is for
#' \code{RunPCA()}/\code{IntegrateLayers()} below. If you need split layers
#' for integration and \code{obj} doesn't already have them, split it (e.g.
#' \code{obj[[assay]] <- split(obj[[assay]], f = obj$orig.ident)}) and
#' re-run normalization on the split object \emph{before} calling this
#' function.
#'
#' \strong{No normalization is performed.} \code{obj} must already have
#' normalized data, variable features, and scaled data (\code{scale.data})
#' for the chosen assay -- e.g. from \code{NormalizeData()} +
#' \code{FindVariableFeatures()} + \code{ScaleData()}, or
#' \code{SCTransform()} -- run on the object \emph{before} subsetting.
#' \code{subset()} carries these layers over for the retained cells/features,
#' and \code{RunPCA()} below requires them. If you need to (re-)normalize the
#' subset from raw counts, do so on the returned object before calling this
#' again, or normalize beforehand.
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
#'     the QC steps below. Other layers (\code{data}, \code{scale.data},
#'     etc.) are left untouched, preserving any existing per-sample split.
#'     If joining fails (this can happen when \code{subset()} leaves a
#'     per-sample layer with zero cells), the layers are left split and QC
#'     metadata is recomputed by summing across them instead.
#'   \item \strong{Recalculate QC metadata.} \code{nCount_<assay>} and
#'     \code{nFeature_<assay>} are recomputed from the subsetted counts,
#'     since Seurat does not update these automatically after \code{subset()}.
#'   \item \strong{Drop empty cells.} Any cell with \code{nCount_<assay> == 0}
#'     after subsetting is removed (can arise when a subset by identity or
#'     metadata leaves a cell with no counts for the chosen assay).
#'   \item \strong{Drop empty features.} Any gene with zero total counts
#'     across all retained cells is removed. This commonly arises after
#'     subsetting to a few clusters and, if left in, can make downstream
#'     normalization (e.g. \code{SCTransform()}) fail with a "subscript out
#'     of bounds" error.
#' }
#'
#' @param obj A Seurat object.
#' @param idents Identity labels (active idents) to keep.
#' @param metadata_col Metadata column to filter on.
#' @param metadata_value Value(s) of \code{metadata_col} to keep.
#' @param cells Explicit character vector of cell barcodes to keep.
#' @param assay Assay to use for QC recalculation and as the object's default
#'   assay going forward. \code{NULL} (default) auto-selects: \code{"RNA"}
#'   if present, else \code{"Spatial"}, else \code{DefaultAssay(obj)}. If
#'   supplied explicitly, it must be one of \code{names(obj@assays)}.
#' @param min_cells_per_gene After subsetting, drop genes detected (count > 0)
#'   in fewer than this many cells. Default 3. Genes that are extremely
#'   sparse in a small subset can otherwise make downstream normalization
#'   (e.g. \code{SCTransform()}) fail with a "subscript out of bounds" error.
#' @param integrate Logical; if TRUE, run \code{IntegrateLayers}. Default TRUE.
#' @param integration_method Integration method (e.g. \code{"HarmonyIntegration"}).
#' @param normalization_method Normalization method that was \emph{already}
#'   applied to \code{obj} (\code{"LogNormalize"} or \code{"SCT"}), passed
#'   through to \code{IntegrateLayers(normalization.method = ...)} when
#'   \code{integrate = TRUE}. Default \code{"LogNormalize"}. Ignored if
#'   \code{integrate = FALSE}.
#' @param new_reduction Name of the integrated reduction (default
#'   \code{"harmony"}).
#' @param dims Number of PCs to use for neighbors/UMAP. Default 15.
#' @param resolution Cluster resolution. Default 0.3.
#' @return The subsetted, re-clustered Seurat object.
#' @importFrom Seurat Idents RunPCA RunUMAP FindNeighbors FindClusters IntegrateLayers
#' @importFrom SeuratObject DefaultAssay Layers JoinLayers LayerData
#' @export
SubsetAndRecluster <- function(obj,
                               idents               = NULL,
                               metadata_col         = NULL,
                               metadata_value       = NULL,
                               cells                = NULL,
                               assay                = NULL,
                               min_cells_per_gene   = 3,
                               integrate            = TRUE,
                               integration_method   = "HarmonyIntegration",
                               normalization_method = c("LogNormalize", "SCT"),
                               new_reduction        = "harmony",
                               dims                 = 15,
                               resolution           = 0.3) {

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
  # Default to RNA, falling back to Spatial, falling back to whatever the
  # object's current default assay is. A user-supplied `assay` is used
  # directly (and validated).
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
  # Best-effort: after subset(), one or more per-sample layers can end up
  # with zero cells, and SeuratObject::JoinLayers() can error with
  # "... error in evaluating the argument 'x' in selecting a method for
  # function 'as.matrix': subscript out of bounds" when asked to join a
  # zero-column layer. If that happens, fall back to leaving the layers
  # split -- the nCount/nFeature recalculation below sums across whatever
  # 'counts*' layers remain, so it works either way.
  #
  # Only the 'counts' layers are joined (via `layers =`), so that any
  # existing per-sample split of 'data'/'scale.data' (e.g. from the
  # standard Seurat v5 integration workflow) is left intact -- this
  # function does not split or re-split layers itself.
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
  # `subset()` does not update these, so they can be stale (reflecting the
  # pre-subset object) unless recomputed here. Summing across all
  # 'counts*' layers (rather than reading a single "counts" layer) works
  # whether or not the join above succeeded.
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
  # After subsetting to a handful of clusters/cells, many genes can end up
  # with zero (or near-zero) detection across every retained cell. Leaving
  # these in causes downstream errors -- in particular SCTransform's vst()
  # can fail with "... error in evaluating the argument 'x' in selecting a
  # method for function 'as.matrix': subscript out of bounds" when its
  # internal gene-binning step (vst.flavor = 'v2') indexes genes that are
  # detected in too few of the retained cells. Recompute, on the
  # (cell-filtered) object, the number of cells each gene is detected in
  # (count > 0) and drop genes below `min_cells_per_gene`.
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

  # ---- PCA -----------------------------------------------------------------
  # No normalization is performed here, and layers are not split/re-split --
  # `obj` is expected to already carry normalized data, variable features,
  # and scale.data (split per-sample if needed for integration) for assay
  # `a` (see the function-level documentation).
  message("--- PCA ---")
  obj <- Seurat::RunPCA(obj, verbose = FALSE)

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
