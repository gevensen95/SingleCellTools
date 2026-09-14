#' Visium spot deconvolution with CARD
#'
#' Wraps \code{CARD::createCARDObject} + \code{CARD::CARD_deconvolution} to
#' estimate per-spot cell-type proportions from a reference single-cell
#' Seurat object. Unlike \code{\link{RunRCTD}}, CARD models spatial
#' correlation between neighboring spots (a conditional autoregressive
#' prior), which tends to give smoother, more spatially coherent
#' proportion maps -- a natural second opinion to compare against RCTD via
#' \code{\link{CompareDeconvolution}}.
#'
#' @param obj A Visium Seurat object, or a (optionally named) list of them.
#' @param reference A single-cell Seurat object with a cell-type column.
#' @param celltype_col Reference metadata column with cell-type labels.
#' @param assay_query Assay to read counts from on the Visium object.
#'   Default \code{"Spatial"}.
#' @param assay_ref Assay on the reference. Default \code{"RNA"}.
#' @param minCountGene,minCountSpot Passed to
#'   \code{CARD::createCARDObject}'s own filtering. Defaults 100 / 5.
#' @param write_metadata Logical; if \code{TRUE} (default), writes
#'   \code{card_<celltype>} proportion columns plus \code{card_dominant}
#'   and \code{card_max_weight}, mirroring \code{\link{RunRCTD}}'s columns.
#' @param verbose Message progress. Default \code{TRUE}.
#' @return The Visium object (or list, matching \code{obj}'s shape) with
#'   \code{obj@misc$card_weights} (spot x cell-type matrix) and, if
#'   \code{write_metadata}, the associated metadata columns.
#' @examples
#' \dontrun{
#' visium <- RunCARD(visium, reference = ref, celltype_col = "cell_type")
#' SpatialFeaturePlot(visium, features = "card_T_cell")
#' }
#' @importFrom Seurat GetTissueCoordinates
#' @export
RunCARD <- function(obj,
                    reference,
                    celltype_col   = NULL,
                    assay_query    = "Spatial",
                    assay_ref      = "RNA",
                    minCountGene   = 100,
                    minCountSpot   = 5,
                    write_metadata = TRUE,
                    verbose        = TRUE) {

  if (!inherits(reference, "Seurat")) stop("`reference` must be a Seurat object.")
  if (is.null(celltype_col) || !celltype_col %in% colnames(reference@meta.data)) {
    stop("`celltype_col` must name a column of reference@meta.data.")
  }
  if (!requireNamespace("CARD", quietly = TRUE)) {
    stop("'CARD' is required. Install with ",
         "remotes::install_github('YMa-Lab/CARD').")
  }

  parsed <- .as_seurat_list(obj)
  objs   <- parsed$objs
  orig_names <- names(objs)

  ref_counts <- Seurat::GetAssayData(reference, assay = assay_ref, layer = "counts")
  ref_meta <- data.frame(
    cellID   = colnames(reference),
    cellType = as.character(reference@meta.data[[celltype_col]]),
    sampleInfo = "ref",
    row.names = colnames(reference),
    stringsAsFactors = FALSE
  )

  objs <- lapply(seq_along(objs), function(i) {
    o <- objs[[i]]
    tag <- if (length(objs) > 1) sprintf(" ('%s')", orig_names[i]) else ""
    .run_card_one(o, ref_counts = ref_counts, ref_meta = ref_meta,
                 assay_query = assay_query, minCountGene = minCountGene,
                 minCountSpot = minCountSpot, write_metadata = write_metadata,
                 verbose = verbose, tag = tag)
  })
  names(objs) <- orig_names

  if (parsed$was_single) return(objs[[1]])
  objs
}

#' @keywords internal
#' @noRd
.run_card_one <- function(obj, ref_counts, ref_meta, assay_query, minCountGene,
                          minCountSpot, write_metadata, verbose, tag) {

  q_counts <- Seurat::GetAssayData(obj, assay = assay_query, layer = "counts")
  coords <- as.data.frame(Seurat::GetTissueCoordinates(obj))
  coord_cols <- intersect(c("x", "y"), colnames(coords))
  if (length(coord_cols) < 2) coord_cols <- colnames(coords)[1:2]
  spatial_location <- data.frame(
    x = coords[[coord_cols[1]]], y = coords[[coord_cols[2]]],
    row.names = rownames(coords)
  )
  common <- intersect(colnames(q_counts), rownames(spatial_location))
  q_counts <- q_counts[, common, drop = FALSE]
  spatial_location <- spatial_location[common, , drop = FALSE]

  if (isTRUE(verbose)) {
    message(sprintf("--- Running CARD%s (%d spots) ---", tag, length(common)))
  }
  card_obj <- CARD::createCARDObject(
    sc_count = ref_counts, sc_meta = ref_meta,
    spatial_count = q_counts, spatial_location = spatial_location,
    ct.varname = "cellType", ct.select = unique(ref_meta$cellType),
    sample.varname = "sampleInfo",
    minCountGene = minCountGene, minCountSpot = minCountSpot
  )
  card_obj <- CARD::CARD_deconvolution(CARD_object = card_obj)
  weights <- as.matrix(card_obj@Proportion_CARD)
  weights <- sweep(weights, 1, pmax(rowSums(weights), 1e-8), "/")

  obj@misc$card_weights <- weights

  if (isTRUE(write_metadata)) {
    ct_names <- colnames(weights)
    ct_safe  <- make.names(ct_names)
    all_cells <- rownames(obj@meta.data)
    match_idx <- match(all_cells, rownames(weights))

    weight_block <- weights[match_idx, , drop = FALSE]
    rownames(weight_block) <- all_cells
    colnames(weight_block) <- paste0("card_", ct_safe)
    weight_df <- as.data.frame(weight_block)

    dom <- ct_names[apply(weights, 1, which.max)]
    weight_df$card_dominant <- dom[match_idx]
    max_weight <- apply(weights, 1, max)
    weight_df$card_max_weight <- max_weight[match_idx]

    obj@meta.data[, colnames(weight_df)] <- weight_df
    if (isTRUE(verbose)) {
      message(sprintf("  Wrote %d per-cell-type columns + card_dominant/card_max_weight.",
                      length(ct_names)))
    }
  }
  obj
}
