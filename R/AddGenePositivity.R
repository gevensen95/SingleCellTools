#' Add per-gene positivity metadata to one or more Seurat objects
#'
#' For every gene in `genes`, adds a logical metadata column to the Seurat
#' object (or every Seurat object in a list) indicating whether each cell has
#' expression strictly greater than `threshold` in the specified assay/layer.
#'
#' Behavior:
#' * `seurat_objects` may be a single Seurat object OR a list of them — the
#'   return type matches the input shape.
#' * `genes` is pre-filtered to the genes present in **every** supplied object
#'   (the intersection of feature names across the chosen assay). Any gene
#'   missing from even one object is dropped with a warning, so the resulting
#'   metadata schema is consistent across the list.
#'
#' @param seurat_objects A Seurat object, or a (optionally named) list of them.
#' @param genes Character vector of gene symbols to score.
#' @param assay Assay to read from. Defaults to each object's DefaultAssay.
#' @param layer Layer to read ('counts', 'data', ...). Default 'counts'.
#' @param threshold Numeric; cells with expression > threshold are positive.
#'   Use the default (0) for counts; bump it (e.g. 1) for log-normalized data.
#' @param suffix String appended to gene name for the metadata column.
#'   Gene 'CD3D' with default suffix becomes column 'CD3D_pos'.
#' @return The same shape as `seurat_objects` (single object or list) with new
#'   logical metadata columns added.
#' @export
AddGenePositivity <- function(seurat_objects,
                              genes,
                              assay     = NULL,
                              layer     = 'counts',
                              threshold = 0,
                              suffix    = '_pos') {
  stopifnot(is.character(genes), length(genes) >= 1,
            is.numeric(threshold), length(threshold) == 1)

  # ---- Normalize input: accept a single Seurat object or a list ----------
  single_input <- inherits(seurat_objects, 'Seurat')
  if (single_input) {
    seurat_objects <- list(seurat_objects)
  }
  if (!is.list(seurat_objects) ||
      !all(vapply(seurat_objects, inherits, logical(1), 'Seurat'))) {
    stop("`seurat_objects` must be a Seurat object or a list of Seurat objects.")
  }

  resolved_assay <- function(o) if (is.null(assay)) Seurat::DefaultAssay(o) else assay

  # ---- Filter `genes` to those present in EVERY object -------------------
  # Intersection across all objects' feature sets so the resulting metadata
  # schema is consistent across the list.
  feature_sets <- lapply(seurat_objects, function(o) {
    rownames(o[[resolved_assay(o)]])
  })
  shared_features <- Reduce(intersect, feature_sets)
  kept    <- intersect(genes, shared_features)
  dropped <- setdiff(genes, kept)
  if (length(dropped)) {
    warning(sprintf("Dropping %d gene(s) not present in every object: %s",
                    length(dropped), paste(dropped, collapse = ', ')))
  }
  if (!length(kept)) {
    warning("No requested genes are present in every object; ",
            "returning input unchanged.")
    return(if (single_input) seurat_objects[[1]] else seurat_objects)
  }
  genes <- kept

  # ---- Add positivity columns to each object -----------------------------
  out <- lapply(seurat_objects, function(o) {
    a    <- resolved_assay(o)
    # FetchData handles v5 split layers and returns a cells x genes data.frame
    expr <- Seurat::FetchData(o, vars = genes, layer = layer, assay = a)
    pos  <- as.data.frame(expr > threshold)
    colnames(pos) <- paste0(genes, suffix)
    SeuratObject::AddMetaData(o, metadata = pos)
  })

  if (single_input) out[[1]] else out
}
