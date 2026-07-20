#' Run BANKSY spatial clustering directly on a Seurat object
#'
#' Thin wrapper around \code{SeuratWrappers::RunBanksy()} (which itself
#' wraps \code{Banksy::computeBanksy()}) so BANKSY can be run on a Seurat
#' object without converting to a \code{SpatialExperiment} first.
#'
#' For imaging-based spatial objects (Xenium, CosMx, etc. -- one or more
#' FOV/centroids images under \code{obj@images}), spatial x/y coordinates
#' are pulled automatically via \code{get_all_coords()} and stashed in
#' metadata columns before calling \code{RunBanksy()}. For objects using
#' Seurat's native spatial framework (Visium) with no FOV images and no
#' explicit \code{dimx}/\code{dimy}, coordinates are left to
#' \code{RunBanksy()} to resolve via \code{GetTissueCoordinates()}
#' directly, matching the "within Seurat's spatial framework" usage shown
#' in the BANKSY-Seurat vignette. You can also bypass coordinate lookup
#' entirely by passing existing metadata column names via \code{dimx}/
#' \code{dimy}.
#'
#' \strong{Caveat (inherited from \code{RunBanksy()})}: do not call
#' \code{ScaleData()} on the resulting BANKSY assay -- \code{RunBanksy()}
#' already fills \code{scale.data} with the lambda-weighted BANKSY matrix,
#' and re-scaling would negate the effect of \code{lambda}.
#'
#' @param obj A Seurat object with spatial information (FOV images and/or
#'   coordinate metadata columns).
#' @param lambda Numeric in `[0,1]`; spatial weight. Low values (~0.2) favor
#'   cell-typing, high values (~0.8) favor spatial domain segmentation.
#' @param assay Assay to compute BANKSY on. Defaults to
#'   \code{DefaultAssay(obj)}.
#' @param slot Slot/layer to pull expression data from. Default
#'   \code{"data"}.
#' @param features \code{"all"}, \code{"variable"}, or a character vector
#'   of feature names.
#' @param k_geom Number of spatial neighbors (numeric, or length-2 for
#'   \code{use_agf = TRUE}).
#' @param image_name Optional single image/FOV name to restrict coordinate
#'   lookup to (passed to \code{get_all_coords()}'s output). If NULL and
#'   the object has multiple FOVs, coordinates from every FOV are used.
#'   Ignored when \code{dimx}/\code{dimy} are supplied directly.
#' @param dimx,dimy Optional names of existing metadata columns already
#'   holding spatial x/y coordinates. If supplied, automatic coordinate
#'   lookup via \code{get_all_coords()} is skipped entirely.
#' @param group Optional metadata column name for multi-sample analysis;
#'   staggers coordinates so cells from different samples don't spatially
#'   overlap (see \code{RunBanksy()} docs).
#' @param split.scale Logical; scale each group's BANKSY matrix
#'   separately. Only relevant when \code{group} is supplied.
#' @param run_pca Logical; if TRUE (default), also runs \code{RunPCA()} on
#'   the resulting BANKSY assay. Skipped automatically if \code{lazy =
#'   TRUE} is passed via \code{...}, since that path already returns a PCA
#'   reduction directly.
#' @param npcs Number of PCs to compute when \code{run_pca = TRUE}.
#' @param assay_name Name for the new BANKSY assay (or, if
#'   \code{lazy = TRUE}, the resulting reduction). Default \code{"BANKSY"}.
#' @param verbose Passed through to \code{RunBanksy()} / \code{RunPCA()}.
#' @param ... Additional arguments passed through to
#'   \code{SeuratWrappers::RunBanksy()} (e.g. \code{use_agf},
#'   \code{spatial_mode}, \code{lazy}).
#' @return The Seurat object with a new BANKSY assay (default assay set to
#'   it) and, if \code{run_pca = TRUE}, a \code{"pca_banksy"} reduction.
#' @examples
#' \dontrun{
#' obj <- RunBanksyWrapper(obj, lambda = 0.2, k_geom = 15)
#' obj <- FindNeighbors(obj, reduction = "pca_banksy", dims = 1:30)
#' obj <- FindClusters(obj, resolution = 0.5)
#' }
#' @importFrom Seurat DefaultAssay RunPCA
#' @export
RunBanksyWrapper <- function(obj,
                             lambda,
                             assay       = NULL,
                             slot        = 'data',
                             features    = 'variable',
                             k_geom      = 15,
                             image_name  = NULL,
                             dimx        = NULL,
                             dimy        = NULL,
                             group       = NULL,
                             split.scale = TRUE,
                             run_pca     = TRUE,
                             npcs        = 30,
                             assay_name  = 'BANKSY',
                             verbose     = TRUE,
                             ...) {

  # Cheap, local checks first so they surface their own errors even on a
  # machine without Banksy/SeuratWrappers installed, rather than being
  # masked by the package-availability errors below.
  if (!inherits(obj, 'Seurat')) stop('`obj` must be a Seurat object.')
  if (is.null(lambda) || lambda < 0 || lambda > 1) {
    stop('`lambda` must be a single numeric value between 0 and 1.')
  }
  if (!requireNamespace('Banksy', quietly = TRUE)) {
    stop("'Banksy' is required. Install with: ",
         "BiocManager::install('Banksy')")
  }
  if (!requireNamespace('SeuratWrappers', quietly = TRUE)) {
    stop("'SeuratWrappers' is required. Install with: ",
         "remotes::install_github('satijalab/seurat-wrappers')")
  }

  if (is.null(assay)) assay <- Seurat::DefaultAssay(obj)

  # A fresh Assay5 only has a "counts" layer until NormalizeData()/
  # SCTransform() populate "data" -- GetAssayData(layer = "data") on that
  # empty layer silently returns a 0x0 matrix rather than erroring, which
  # then surfaces many calls later, deep inside SeuratWrappers::RunBanksy(),
  # as a cryptic Matrix "subscript out of bounds" once it tries to index
  # that empty matrix by feature name. Catch it here with an actionable
  # message instead.
  layer_data <- tryCatch(Seurat::GetAssayData(obj, assay = assay, layer = slot),
                         error = function(e) NULL)
  if (is.null(layer_data) || nrow(layer_data) == 0 || ncol(layer_data) == 0) {
    stop("The '", slot, "' layer of assay '", assay, "' is empty (",
         if (is.null(layer_data)) '0 x 0' else paste(nrow(layer_data), 'x', ncol(layer_data)),
         "). Run NormalizeData() (or SCTransform()) on `obj` first, or pass ",
         "slot = 'counts' to run BANKSY directly on raw counts.")
  }

  dots <- list(...)
  lazy_requested <- isTRUE(dots$lazy)

  # ---- Resolve spatial x/y columns ----------------------------------------
  if (is.null(dimx) || is.null(dimy)) {
    has_fovs <- length(Seurat::Images(obj)) > 0
    if (has_fovs) {
      message('--- Pulling spatial coordinates via get_all_coords() ---')
      coords <- get_all_coords(obj)

      # Normalize column names across spatial technologies -- same
      # imagecol/imagerow fallback AnnotateRegions() uses.
      if (all(c('imagecol', 'imagerow') %in% colnames(coords))) {
        coords$x <- coords$imagecol
        coords$y <- coords$imagerow
      }
      if (!all(c('x', 'y') %in% colnames(coords))) {
        stop("Couldn't find x/y (or imagecol/imagerow) coordinates via get_all_coords().")
      }
      if (!'cell' %in% colnames(coords)) coords$cell <- rownames(coords)

      if (!is.null(image_name)) {
        if (!'fov' %in% colnames(coords) || !image_name %in% coords$fov) {
          stop('`image_name` "', image_name, '" not found among obj images: ',
               paste(Seurat::Images(obj), collapse = ', '))
        }
        coords <- coords[coords$fov == image_name, , drop = FALSE]
      } else if ('fov' %in% colnames(coords) && length(unique(coords$fov)) > 1) {
        message('  Multiple FOVs found and no `image_name` given -- using ',
                'coordinates from all ', length(unique(coords$fov)),
                ' FOVs. Pass `group` if these are separate samples that ',
                'should not be spatially staggered together.')
      }

      dup <- duplicated(coords$cell)
      if (any(dup)) {
        warning(sum(dup), ' cell(s) appeared in more than one FOV; keeping the first occurrence.')
        coords <- coords[!dup, , drop = FALSE]
      }
      rownames(coords) <- coords$cell

      common <- intersect(colnames(obj), rownames(coords))
      if (length(common) == 0) {
        stop('No cells in `obj` matched coordinates pulled from get_all_coords().')
      }
      if (length(common) < ncol(obj)) {
        warning(ncol(obj) - length(common), ' / ', ncol(obj),
                ' cells had no spatial coordinates and will be dropped by RunBanksy().')
      }

      obj <- Seurat::AddMetaData(obj, metadata = coords[common, 'x', drop = FALSE],
                                 col.name = 'banksy_x')
      obj <- Seurat::AddMetaData(obj, metadata = coords[common, 'y', drop = FALSE],
                                 col.name = 'banksy_y')
      dimx <- 'banksy_x'
      dimy <- 'banksy_y'
    } else {
      message('--- No FOV images found on `obj`; leaving coordinates to ',
              'RunBanksy() to resolve via GetTissueCoordinates() directly ',
              '(Seurat native spatial framework, e.g. Visium). Pass `dimx`/',
              '`dimy` explicitly if that lookup does not apply here. ---')
    }
  }

  # ---- Run BANKSY -----------------------------------------------------------
  message(sprintf('--- Running BANKSY (lambda = %s, k_geom = %s, assay = %s) ---',
                  lambda, paste(k_geom, collapse = ','), assay))
  obj <- SeuratWrappers::RunBanksy(obj, lambda = lambda, assay = assay, slot = slot,
                                   dimx = dimx, dimy = dimy, features = features,
                                   k_geom = k_geom, group = group,
                                   split.scale = split.scale,
                                   assay_name = assay_name, verbose = verbose, ...)

  # ---- Optional PCA on the BANKSY matrix -----------------------------------
  if (isTRUE(run_pca) && !lazy_requested) {
    message('--- Running PCA on BANKSY assay ---')
    obj <- Seurat::RunPCA(obj, assay = assay_name, features = rownames(obj),
                          npcs = npcs, reduction.name = 'pca_banksy',
                          reduction.key = 'PCABANKSY_', verbose = verbose)
  } else if (isTRUE(run_pca) && lazy_requested) {
    message('--- `lazy = TRUE` already produced a PCA reduction ("', assay_name,
            '"); skipping the extra RunPCA() call. ---')
  }

  obj
}
