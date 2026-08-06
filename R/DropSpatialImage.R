#' Free Image/Molecule Memory on Already-Built Spatial Objects
#'
#' Retroactively frees the memory \code{@images} holds on Seurat object(s) --
#' useful right before a step that doesn't need images/coordinates at all
#' (merge, integration, pseudobulk DE, saving to disk). Handles both kinds
#' of spatial image this package builds: pixel-backed \code{VisiumV1}
#' images (from \code{\link{CreateVisiumObjects}}) and coordinate-only
#' \code{FOV} images (Xenium/CosMx/MERFISH-style, e.g. from
#' \code{\link{LoadXenium2}} or \code{Seurat::LoadNanostring()}), detected
#' by capability rather than class name (same approach as
#' \code{\link{SpatialObjectInfo}}).
#'
#' @param obj A Seurat object, or a (optionally named) list of them.
#' @param mode One of \code{"remove"} (default) or \code{"downgrade"}.
#'
#'   \code{"remove"} drops \code{@images} entirely, regardless of image
#'   type -- the most memory you can free, but afterward nothing that needs
#'   \code{Seurat::GetTissueCoordinates()} (e.g. re-running
#'   \code{\link{EdgeDetectionVisium}}, \code{\link{GenerateQCReport}}'s
#'   spatial section, or any FOV-based function in this package like
#'   \code{\link{NeighborhoodEnrichment}}) will work on this object anymore;
#'   only appropriate once you're done with anything spatial-coordinate-based.
#'
#'   \code{"downgrade"} only applies to pixel-backed images -- it rebuilds
#'   them at lowres resolution (the same result \code{image_backend =
#'   "deferred"} would have produced at construction time), requiring
#'   \code{obj@misc$visium_image_dir} to be stashed (every object built by
#'   the current \code{\link{CreateVisiumObjects}} has it regardless of
#'   \code{image_backend}). Coordinate-only (FOV) images have no lowres
#'   equivalent to downgrade to, so they're left untouched with an
#'   explanatory message -- their attached \code{molecules} are already the
#'   QV-filtered subset from construction time, not the full raw table (see
#'   \code{\link{QueryXeniumMolecules}} for re-querying that on disk without
#'   holding it in memory); use \code{mode = "remove"} on those objects
#'   instead if you want the FOV gone entirely. A no-op (with a message) if
#'   a pixel-backed image is already deferred, or if no images are attached.
#' @return The updated Seurat object, or list of them -- matches the shape
#'   of \code{obj} (a single object in, a single object out; a list in, a
#'   list out).
#' @export
DropSpatialImage <- function(obj, mode = c("remove", "downgrade")) {
  mode <- match.arg(mode)
  parsed <- .as_seurat_list(obj)
  objs <- parsed$objs
  was_single <- parsed$was_single

  .is_pixel <- function(img) {
    arr <- tryCatch(Seurat::GetImage(img, mode = "raw"), error = function(e) NULL)
    !is.null(arr) && length(dim(arr)) >= 2L
  }

  objs <- lapply(objs, function(o) {
    if (mode == "remove") {
      o@images <- list()
      return(o)
    }

    # mode == "downgrade"
    img_names <- names(o@images)
    is_pixel_vec <- vapply(img_names, function(in_) .is_pixel(o@images[[in_]]), logical(1))
    fov_names   <- img_names[!is_pixel_vec]
    pixel_names <- img_names[is_pixel_vec]

    if (length(fov_names) > 0L) {
      message(
        "Skipping ", length(fov_names), " coordinate-only image(s) (",
        paste(fov_names, collapse = ", "), ") -- 'downgrade' only applies ",
        "to pixel-backed images; there's no lowres equivalent for FOV-style ",
        "objects. Their attached molecules are already the QV-filtered ",
        "subset from construction time (see LoadXenium2(microns_lazy = ...) ",
        "/ QueryXeniumMolecules() to re-query the full table on disk instead ",
        "of holding it in memory). Use `mode = 'remove'` if you want these ",
        "dropped entirely."
      )
    }

    if (length(img_names) > 0L && length(pixel_names) == 0L) {
      # Only coordinate-only images were present (already messaged above) --
      # nothing pixel-backed to downgrade, and no Visium prerequisites
      # (visium_image_dir) are relevant for an object like this.
      return(o)
    }

    if (!is.null(o@misc$hires_image_path)) {
      message("Already deferred (lowres attached, hires stashed) -- nothing to do.")
      return(o)
    }
    image_dir <- o@misc$visium_image_dir
    if (is.null(image_dir)) {
      stop("`obj@misc$visium_image_dir` not found, so there's no known ",
           "location to re-read a lowres image from. This is stashed ",
           "automatically by the current CreateVisiumObjects() -- rebuild ",
           "with it, or use `mode = 'remove'` instead if you don't need ",
           "images at all.")
    }
    if (length(o@images) == 0L) {
      message("No images attached -- nothing to downgrade.")
      return(o)
    }

    for (in_ in pixel_names) {
      old <- o@images[[in_]]
      vis.image <- Seurat::Read10X_Image(
        image.dir  = image_dir,
        image.name = "tissue_lowres_image.png",
        filter.matrix = TRUE
      )
      vis.image@assay <- old@assay
      vis.image@key   <- old@key
      o@images[[in_]] <- vis.image
    }
    o@misc$hires_image_path <- file.path(image_dir, "tissue_hires_image.png")
    o
  })

  if (was_single) objs[[1]] else objs
}
