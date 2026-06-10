#' Strip normalized layers, variable features, and reductions
#'
#' Removes the intermediate artifacts that normalization + PCA + UMAP +
#' integration leave on a Seurat object: the \code{data} and
#' \code{scale.data} layers, the variable-features set, and **every**
#' reduction (\code{pca}, \code{umap}, \code{harmony}, \code{tsne}, ...).
#' By default the result carries only counts plus metadata.
#'
#' Handles single objects or lists, v5 split layers (\code{data.1},
#' \code{scale.data.2}, ...), and assay-suffixed reductions
#' (\code{pca.RNA}, \code{umap.SCT}).
#'
#' @param obj A Seurat object OR a (optionally named) list of Seurat objects.
#' @param assay Assay whose layers to strip. Default \code{"RNA"}.
#' @param keep_reductions Character vector of reduction names to **preserve**.
#'   Everything else is dropped. Default \code{NULL} drops all reductions.
#'   Example: \code{keep_reductions = "harmony"} keeps harmony but drops
#'   \code{pca}, \code{umap}, \code{tsne}, etc.
#' @return The same shape as the input.
#' @importFrom SeuratObject Assays Layers Reductions VariableFeatures<- DefaultAssay<-
#' @export
strip_workflow_artifacts <- function(obj,
                                     assay           = "RNA",
                                     keep_reductions = NULL) {
  single_input <- inherits(obj, "Seurat")
  if (single_input) obj <- list(obj)

  if (!is.list(obj) ||
      !all(vapply(obj, inherits, logical(1), "Seurat"))) {
    stop("`obj` must be a Seurat object or a list of Seurat objects.")
  }

  out <- lapply(obj, function(o) {
    # If an SCT assay exists, drop it whole (its layers go with it).
    if ("SCT" %in% SeuratObject::Assays(o)) {
      SeuratObject::DefaultAssay(o) <- assay
      o[["SCT"]] <- NULL
    }

    # Drop normalized layers, including any split variants like data.1, data.2
    rna <- o[[assay]]
    target_patterns <- c("^data(\\.|$)", "^scale\\.data(\\.|$)")
    lyrs_to_drop <- unique(unlist(lapply(
      target_patterns,
      function(p) grep(p, names(rna@layers), value = TRUE)
    )))
    for (lyr in lyrs_to_drop) {
      rna@layers[[lyr]] <- NULL
    }

    # Clear variable features. Assay5 stores them as vf_* meta.data columns;
    # legacy Assay uses @var.features.
    if (inherits(rna, "Assay5")) {
      vf_cols <- grep("^vf_", colnames(rna@meta.data), value = TRUE)
      if (length(vf_cols)) rna@meta.data[vf_cols] <- list(NULL)
    } else if (.hasSlot(rna, "var.features")) {
      rna@var.features <- character(0)
    }
    o[[assay]] <- rna

    # Drop ALL reductions except any the caller explicitly asked to keep.
    all_reds <- names(o@reductions)
    red_to_drop <- setdiff(all_reds, keep_reductions)
    for (red in red_to_drop) {
      o@reductions[[red]] <- NULL
    }
    o
  })

  if (single_input) out[[1]] else out
}
