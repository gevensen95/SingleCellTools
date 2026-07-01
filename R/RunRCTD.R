#' Visium spot deconvolution with RCTD (spacexr)
#'
#' Wraps \code{spacexr::create.RCTD} + \code{spacexr::run.RCTD} to
#' estimate per-spot cell-type proportions from a reference single-cell
#' Seurat object. RCTD is the recommended primary annotation strategy for
#' Visium because each spot is a mixture of 1-10 cells of different types;
#' a winner-takes-all classifier (like \code{AnnotateClusters}) loses
#' minority populations by construction.
#'
#' Three modes are supported (RCTD's \code{doublet_mode}):
#' \describe{
#'   \item{\code{"full"} (default)}{Assumes each spot is a mixture of many
#'     cell types. Returns full proportion matrix. Best for standard
#'     Visium.}
#'   \item{\code{"doublet"}}{Assumes each spot has one or two cell types.
#'     Faster, cleaner, but less appropriate for high-density tissue.}
#'   \item{\code{"multi"}}{Iteratively fits up to \code{max_cores_multi}
#'     types per spot. Middle ground.}
#' }
#'
#' Adds the resulting weights matrix (per-cell-type proportions) to the
#' Visium object under \code{obj@misc$rctd_weights}, and — when
#' \code{write_metadata = TRUE} — writes each cell type's per-spot
#' proportion as a metadata column named \code{rctd_<celltype>}. The
#' dominant cell type per spot is written to
#' \code{obj@meta.data$rctd_dominant}.
#'
#' @param obj A Visium Seurat object.
#' @param reference A single-cell Seurat object with a cell-type column
#'   for the reference.
#' @param celltype_col Reference metadata column with cell-type labels.
#' @param assay_query Assay to read counts from on the Visium object.
#'   Default \code{"Spatial"}.
#' @param assay_ref Assay on the reference. Default \code{"RNA"}.
#' @param mode RCTD \code{doublet_mode}. See Details. Default \code{"full"}.
#' @param max_cells_per_ref_celltype Downsample the reference to this
#'   many cells per cell type before building the RCTD reference. Default
#'   10000. Set to \code{Inf} to disable.
#' @param CELL_MIN_INSTANCE Minimum cells per reference cell type required
#'   to keep it (passed through to \code{create.RCTD}). Default 25.
#' @param write_metadata Logical; if TRUE (default) writes proportion
#'   columns and dominant-type column into obj@meta.data.
#' @param n_cores Number of cores for RCTD. Default 4.
#' @return The Visium object with \code{obj@misc$rctd_weights} (spot x
#'   cell-type matrix) and the associated metadata columns (if
#'   \code{write_metadata}).
#' @examples
#' \dontrun{
#' visium <- RunRCTD(visium,
#'                   reference    = pbmc_ref,
#'                   celltype_col = "cell_type",
#'                   mode         = "full",
#'                   n_cores      = 8)
#' colnames(visium@misc$rctd_weights)
#' # Downstream: SpatialFeaturePlot on rctd_<celltype> columns
#' SpatialFeaturePlot(visium, features = "rctd_T_cell")
#' }
#' @importFrom Seurat DefaultAssay GetTissueCoordinates
#' @importFrom SeuratObject LayerData
#' @export
RunRCTD <- function(obj,
                    reference,
                    celltype_col               = NULL,
                    assay_query                = "Spatial",
                    assay_ref                  = "RNA",
                    mode                       = c("full", "doublet", "multi"),
                    max_cells_per_ref_celltype = 10000,
                    CELL_MIN_INSTANCE          = 25,
                    write_metadata             = TRUE,
                    n_cores                    = 4) {

  mode <- match.arg(mode)
  if (!requireNamespace("spacexr", quietly = TRUE)) {
    stop("'spacexr' is required. Install with ",
         "remotes::install_github('dmcable/spacexr').")
  }
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (!inherits(reference, "Seurat")) {
    stop("`reference` must be a Seurat object.")
  }
  if (is.null(celltype_col) ||
      !celltype_col %in% colnames(reference@meta.data)) {
    stop("`celltype_col` must name a column of reference@meta.data.")
  }

  # ---- Build RCTD reference ----------------------------------------------
  ref_counts <- SeuratObject::LayerData(reference,
                                        assay = assay_ref, layer = "counts")
  ref_types  <- factor(as.character(reference@meta.data[[celltype_col]]))
  names(ref_types) <- colnames(reference)

  # Optional downsample per cell type
  if (is.finite(max_cells_per_ref_celltype)) {
    keep <- unlist(lapply(levels(ref_types), function(ct) {
      cells <- names(ref_types)[ref_types == ct]
      if (length(cells) > max_cells_per_ref_celltype) {
        cells <- sample(cells, max_cells_per_ref_celltype)
      }
      cells
    }))
    ref_counts <- ref_counts[, keep, drop = FALSE]
    ref_types  <- ref_types[keep]
  }

  message(sprintf("--- Building RCTD reference (%d cells, %d types) ---",
                  length(ref_types), nlevels(ref_types)))
  rctd_ref <- spacexr::Reference(
    counts     = ref_counts,
    cell_types = ref_types,
    nUMI       = as.numeric(Matrix::colSums(ref_counts))
  )

  # ---- Build spatial "puck" ----------------------------------------------
  # RCTD wants a Puck object: coords + counts + nUMI. Pull from the first
  # image (Visium objects usually have exactly one).
  img_name <- names(obj@images)[1]
  coords <- Seurat::GetTissueCoordinates(obj[[img_name]])
  # Normalize the coords columns depending on Seurat version
  if ("cell" %in% colnames(coords)) {
    rownames(coords) <- coords$cell
    coords <- coords[, c("x", "y")]
  } else if (all(c("imagerow", "imagecol") %in% colnames(coords))) {
    coords <- coords[, c("imagerow", "imagecol")]
    colnames(coords) <- c("x", "y")
  } else {
    coords <- coords[, 1:2]
    colnames(coords) <- c("x", "y")
  }
  cells_in_both <- intersect(rownames(coords), colnames(obj))
  coords <- coords[cells_in_both, , drop = FALSE]

  q_counts <- SeuratObject::LayerData(obj, assay = assay_query,
                                      layer = "counts")[, cells_in_both,
                                                        drop = FALSE]
  message(sprintf("--- Building query puck (%d spots) ---", ncol(q_counts)))
  puck <- spacexr::SpatialRNA(
    coords = coords,
    counts = q_counts,
    nUMI   = as.numeric(Matrix::colSums(q_counts))
  )

  # ---- Run RCTD -----------------------------------------------------------
  message(sprintf("--- Running RCTD (mode = %s, cores = %d) ---",
                  mode, n_cores))
  rctd <- spacexr::create.RCTD(
    spatialRNA        = puck,
    reference         = rctd_ref,
    max_cores         = n_cores,
    CELL_MIN_INSTANCE = CELL_MIN_INSTANCE
  )
  rctd <- spacexr::run.RCTD(rctd, doublet_mode = mode)

  # ---- Extract weights and normalize --------------------------------------
  if (mode == "full") {
    weights <- as.matrix(rctd@results$weights)
    # Row-normalize to proportions (they can be counts otherwise)
    weights <- sweep(weights, 1, pmax(rowSums(weights), 1e-8), "/")
  } else {
    # For doublet / multi mode, results is a data frame; construct a matrix
    df <- rctd@results$results_df
    ct <- levels(factor(c(as.character(df$first_type),
                          as.character(df$second_type))))
    weights <- matrix(0, nrow = nrow(df), ncol = length(ct),
                      dimnames = list(rownames(df), ct))
    for (i in seq_len(nrow(df))) {
      w1 <- df$weight_first[i]
      w2 <- df$weight_second[i]
      weights[i, as.character(df$first_type[i])]  <- w1
      if (!is.na(df$second_type[i])) {
        weights[i, as.character(df$second_type[i])] <- w2
      }
    }
    weights <- sweep(weights, 1, pmax(rowSums(weights), 1e-8), "/")
  }

  obj@misc$rctd_weights <- weights

  # ---- Optional metadata columns -----------------------------------------
  if (isTRUE(write_metadata)) {
    # One column per cell type
    ct_names <- colnames(weights)
    ct_safe  <- make.names(ct_names)
    for (i in seq_along(ct_names)) {
      col <- paste0("rctd_", ct_safe[i])
      obj@meta.data[[col]] <- NA_real_
      obj@meta.data[rownames(weights), col] <- weights[, i]
    }
    # Dominant type
    dom <- ct_names[apply(weights, 1, which.max)]
    obj@meta.data$rctd_dominant <- NA_character_
    obj@meta.data[rownames(weights), "rctd_dominant"] <- dom
    message(sprintf("  Wrote %d per-cell-type columns + rctd_dominant.",
                    length(ct_names)))
  }

  obj
}
