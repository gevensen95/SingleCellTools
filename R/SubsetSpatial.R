#' Subset a spatial Seurat object without corrupting \code{@images}
#'
#' \code{subset()} on a Seurat object is supposed to subset every attached
#' image (\code{obj@images}) in lockstep with the rest of the object, but a
#' known, still-unresolved upstream Seurat bug
#' (\url{https://github.com/satijalab/seurat/issues/8848}) can leave one or
#' more images holding cells that no longer exist in the subsetted object.
#' \code{subset()} itself won't complain -- the inconsistency only surfaces
#' later, whenever the next function happens to call
#' \code{methods::validObject()} (\code{RunPCA()}, \code{PrepSCTFindMarkers()},
#' \code{ScaleData()}, and others all do), with \code{invalid class "Seurat"
#' object: All cells in images must be present in the Seurat object}.
#'
#' \code{SubsetSpatial()} sidesteps this by letting \code{base::subset()}
#' handle everything about the subset \emph{except} \code{@images}, then
#' explicitly re-subsetting each image against the object's own,
#' already-subsetted cell names -- the exact invariant the validity check
#' enforces -- rather than relying on \code{subset()}'s own per-image
#' handling. Unlike \code{\link{DropSpatialImage}}, this keeps images
#' attached (correctly re-subset, not removed), so a subsetted object is
#' still usable for spatial plotting (\code{\link{SpatialDimPlotFixed}},
#' \code{\link{SpatialFeaturePlotFixed}}, \code{\link{SpatialConcordance}},
#' etc.) afterward.
#'
#' Works on any image type attached to \code{obj@images} -- pixel-backed
#' Visium images (\code{VisiumV1}/\code{VisiumV2}) and coordinate-only FOV
#' images (Xenium/CosMx/MERFISH-style) alike -- since the fix (re-subset
#' each image directly, by cell name) doesn't depend on image class. Unlike
#' \code{\link{subset_opt}}, which is written specifically for CosMx/Xenium
#' \code{FOV} objects (it reaches into a \code{centroids} sub-slot and a
#' \code{molecules} slot that pixel-backed Visium images don't have),
#' \code{SubsetSpatial()} makes no assumptions about image internals beyond
#' "has cells, can be subset by cell name" -- true of every
#' \code{SpatialImage}-derived class.
#'
#' @param object A Seurat object.
#' @param subset Logical expression indicating cells to keep, evaluated
#'   against \code{object}'s metadata -- same non-standard-evaluation
#'   convention as \code{base::subset()}'s own \code{subset} argument (not
#'   a string). Optional.
#' @param cells Optional character vector of cell names (or integer indices
#'   into \code{Cells(object)}) to keep. Combined with \code{subset}/
#'   \code{idents} if more than one is given, same as \code{base::subset()}.
#' @param idents Optional vector of identity classes to keep.
#' @param features Optional character vector of features (genes) to keep.
#'   Only affects assay data, not \code{@images}.
#' @param ... Passed through to \code{base::subset()}.
#' @return \code{object} subsetted, with every attached image explicitly
#'   re-subset to match the object's final cell set (an image left with
#'   zero cells is dropped, with a message).
#' @examples
#' \dontrun{
#' # Keeps images, correctly re-subset -- safe to plot afterward
#' hep <- SubsetSpatial(visium, subset = Annotation == "Hepatocytes")
#' SpatialDimPlotFixed(hep, group.by = "Zone_final")
#' }
#' @export
SubsetSpatial <- function(object, subset = NULL, cells = NULL, idents = NULL,
                          features = NULL, ...) {
  .assert_seurat(object)
  subset_q <- rlang::enquo(subset)

  resolved_cells <- if (!rlang::quo_is_null(subset_q) || !is.null(idents)) {
    # WhichCells() treats a NULL `expression` (whether passed literally or
    # as a quosure wrapping NULL) differently from `expression` simply not
    # being supplied at all -- passing NULL explicitly makes it try to
    # FetchData() an empty variable set and error with "None of the
    # requested variables were found", even though the *only* filter
    # actually wanted here is `idents`. So only forward `expression` when
    # the caller actually gave a `subset=` to preserve (confirmed via a
    # standalone repro against SeuratObject::WhichCells() directly).
    if (rlang::quo_is_null(subset_q)) {
      SeuratObject::WhichCells(object, cells = cells, idents = idents,
                               return.null = TRUE)
    } else {
      SeuratObject::WhichCells(object, cells = cells, idents = idents,
                               expression = subset_q, return.null = TRUE)
    }
  } else if (is.null(cells)) {
    SeuratObject::Cells(object)
  } else if (is.numeric(cells)) {
    SeuratObject::Cells(object)[cells]
  } else {
    cells
  }

  obj_subset <- base::subset(object, cells = resolved_cells, features = features, ...)

  if (length(object@images) == 0L) {
    return(obj_subset)
  }

  kept_cells <- colnames(obj_subset)
  new_images <- list()
  for (img in names(object@images)) {
    orig_image <- object@images[[img]]
    img_cells  <- intersect(SeuratObject::Cells(orig_image), kept_cells)
    if (length(img_cells) == 0L) {
      message("SubsetSpatial(): image '", img, "' has no cells left after ",
              "subsetting -- dropping it.")
      next
    }
    new_images[[img]] <- base::subset(orig_image, cells = img_cells)
  }
  obj_subset@images <- new_images

  obj_subset
}
