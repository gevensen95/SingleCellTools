#' Visium spot deconvolution with SPOTlight
#'
#' Wraps \code{SPOTlight::SPOTlight} (NMF-based) to estimate per-spot
#' cell-type proportions from a reference single-cell Seurat object.
#' Faster than \code{\link{RunRCTD}}/\code{\link{RunCARD}} and a useful
#' third opinion for \code{\link{CompareDeconvolution}}; marker genes per
#' cell type are computed internally via \code{Seurat::FindAllMarkers}
#' unless you supply your own.
#'
#' @param obj A Visium Seurat object, or a (optionally named) list of them.
#' @param reference A single-cell Seurat object with a cell-type column.
#' @param celltype_col Reference metadata column with cell-type labels.
#' @param assay_query Assay to read counts from on the Visium object.
#'   Default \code{"Spatial"}.
#' @param assay_ref Assay on the reference. Default \code{"RNA"}.
#' @param markers Optional marker data frame (as from
#'   \code{Seurat::FindAllMarkers}, needing \code{gene}, \code{cluster},
#'   \code{avg_log2FC} columns) to skip computing markers internally.
#'   \code{NULL} (default) runs \code{FindAllMarkers} on \code{reference}.
#' @param n_top_markers Top markers per cell type to use (by
#'   \code{avg_log2FC}) when computing them internally. Default 100.
#' @param write_metadata Logical; if \code{TRUE} (default), writes
#'   \code{spotlight_<celltype>} proportion columns plus
#'   \code{spotlight_dominant} and \code{spotlight_max_weight}.
#' @param verbose Message progress. Default \code{TRUE}.
#' @return The Visium object (or list, matching \code{obj}'s shape) with
#'   \code{obj@misc$spotlight_weights} (spot x cell-type matrix) and, if
#'   \code{write_metadata}, the associated metadata columns.
#' @examples
#' \dontrun{
#' visium <- RunSPOTlight(visium, reference = ref, celltype_col = "cell_type")
#' SpatialFeaturePlot(visium, features = "spotlight_T_cell")
#' }
#' @importFrom Seurat FindAllMarkers GetAssayData Idents Idents<-
#' @export
RunSPOTlight <- function(obj,
                         reference      = NULL,
                         celltype_col   = NULL,
                         assay_query    = "Spatial",
                         assay_ref      = "RNA",
                         markers        = NULL,
                         n_top_markers  = 100,
                         write_metadata = TRUE,
                         verbose        = TRUE) {

  if (!inherits(reference, "Seurat")) stop("`reference` must be a Seurat object.")
  if (is.null(celltype_col) || !celltype_col %in% colnames(reference@meta.data)) {
    stop("`celltype_col` must name a column of reference@meta.data.")
  }
  if (!requireNamespace("SPOTlight", quietly = TRUE)) {
    stop("'SPOTlight' is required. Install with ",
         "remotes::install_github('MarcElosua/SPOTlight').")
  }

  parsed <- .as_seurat_list(obj)
  objs   <- parsed$objs
  orig_names <- names(objs)

  ref_counts <- as.matrix(Seurat::GetAssayData(reference, assay = assay_ref, layer = "counts"))
  ref_groups <- as.character(reference@meta.data[[celltype_col]])

  mgs <- markers
  if (is.null(mgs)) {
    if (isTRUE(verbose)) message("--- Computing marker genes (FindAllMarkers) ---")
    Seurat::Idents(reference) <- factor(ref_groups)
    mgs <- Seurat::FindAllMarkers(reference, assay = assay_ref, only.pos = TRUE,
                                  verbose = isTRUE(verbose))
    mgs <- do.call(rbind, lapply(split(mgs, mgs$cluster), function(d) {
      d <- d[order(-d$avg_log2FC), ]
      d[seq_len(min(n_top_markers, nrow(d))), ]
    }))
  }

  objs <- lapply(seq_along(objs), function(i) {
    o <- objs[[i]]
    tag <- if (length(objs) > 1) sprintf(" ('%s')", orig_names[i]) else ""
    .run_spotlight_one(o, ref_counts = ref_counts, ref_groups = ref_groups,
                       mgs = mgs, assay_query = assay_query,
                       write_metadata = write_metadata, verbose = verbose, tag = tag)
  })
  names(objs) <- orig_names

  if (parsed$was_single) return(objs[[1]])
  objs
}

#' @keywords internal
#' @noRd
.run_spotlight_one <- function(obj, ref_counts, ref_groups, mgs, assay_query,
                               write_metadata, verbose, tag) {

  q_counts <- as.matrix(Seurat::GetAssayData(obj, assay = assay_query, layer = "counts"))
  if (isTRUE(verbose)) {
    message(sprintf("--- Running SPOTlight%s (%d spots) ---", tag, ncol(q_counts)))
  }
  decon <- SPOTlight::SPOTlight(
    x = ref_counts, y = q_counts, groups = ref_groups,
    mgs = mgs, weight_id = "avg_log2FC", group_id = "cluster", gene_id = "gene"
  )
  weights <- as.matrix(decon$mat)
  weights <- sweep(weights, 1, pmax(rowSums(weights), 1e-8), "/")

  obj@misc$spotlight_weights <- weights

  if (isTRUE(write_metadata)) {
    ct_names <- colnames(weights)
    ct_safe  <- make.names(ct_names)
    all_cells <- rownames(obj@meta.data)
    match_idx <- match(all_cells, rownames(weights))

    weight_block <- weights[match_idx, , drop = FALSE]
    rownames(weight_block) <- all_cells
    colnames(weight_block) <- paste0("spotlight_", ct_safe)
    weight_df <- as.data.frame(weight_block)

    dom <- ct_names[apply(weights, 1, which.max)]
    weight_df$spotlight_dominant <- dom[match_idx]
    max_weight <- apply(weights, 1, max)
    weight_df$spotlight_max_weight <- max_weight[match_idx]

    obj@meta.data[, colnames(weight_df)] <- weight_df
    if (isTRUE(verbose)) {
      message(sprintf("  Wrote %d per-cell-type columns + spotlight_dominant/spotlight_max_weight.",
                      length(ct_names)))
    }
  }
  obj
}
