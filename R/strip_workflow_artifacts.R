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

    # Clear variable features. On Assay5 (SeuratObject >= 5), a plain
    # VariableFeatures(obj) <- genes call stores the result in two
    # feature-level meta.data columns literally named "var.features" (the
    # gene names) and "var.features.rank" (their rank) -- NOT as
    # "vf_<method>_variable" boolean columns (that naming is specific to
    # FindVariableFeatures()'s per-method diagnostics) and NOT in an
    # @var.features slot (Assay5 has no such slot at all).
    #
    # These columns must be DROPPED entirely, not just NA'd out.
    # SeuratObject:::.GetVariableFeatures() -- what VariableFeatures()
    # actually reads through -- has a real bug for the "everything is NA"
    # case: it computes `nfeatures <- nfeatures %||% sum(label_mask)`,
    # which is 0 when every label is NA/FALSE, and then does
    # `variable_features[1:nfeatures]`. In R, `1:0` is `c(1, 0)`, not an
    # empty sequence, so indexing a character(0) result at position 1
    # returns a single out-of-bounds NA -- VariableFeatures() reports
    # length 1 (one NA), never length 0, as long as the columns exist at
    # all, however they're populated. Removing the columns instead makes
    # .GetVariableFeatures() hit its very first check (neither column
    # present -> return NULL) and skip that broken path entirely.
    vf_cols <- c(
      grep("^vf_", colnames(rna@meta.data), value = TRUE),
      intersect(c("var.features", "var.features.rank"), colnames(rna@meta.data))
    )
    if (length(vf_cols)) rna@meta.data[vf_cols] <- list(NULL)
    if (.hasSlot(rna, "var.features")) rna@var.features <- character(0)
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
