#' Retrieve or Attach the Full-Resolution Visium Image for a Deferred Sample
#'
#' Complement to \code{CreateVisiumObjects(..., image_backend = "deferred")}:
#' that mode attaches the small \code{tissue_lowres_image.png} to every
#' sample (so a list of N samples doesn't mean N times ~100MB of decoded
#' hires PNGs before you've done anything with them) and stashes the hires
#' image's path. This function decodes that hires image on demand, for just
#' the sample(s) you actually need full detail on.
#'
#' @param obj A Seurat object built with
#'   \code{CreateVisiumObjects(..., image_backend = "deferred")}.
#' @param attach Logical; if \code{FALSE} (default), \code{obj} is left
#'   untouched and the decoded image array is returned directly (as
#'   \code{png::readPNG()} returns it) -- useful for a one-off custom plot
#'   without permanently swapping what's attached to the object. If
#'   \code{TRUE}, rebuilds \code{obj@images$image} as a full-resolution
#'   \code{VisiumV1} image (via \code{Seurat::Read10X_Image()} on the
#'   stashed hires path, exactly as \code{image_backend = "eager"} would
#'   have built it originally) and returns the updated Seurat object.
#' @return If \code{attach = FALSE} (default), a numeric array (the decoded
#'   image). If \code{attach = TRUE}, the updated Seurat object.
#' @export
GetHiresVisiumImage <- function(obj, attach = FALSE) {
  .assert_seurat(obj)

  hires_path <- obj@misc$hires_image_path
  if (is.null(hires_path)) {
    stop("`obj@misc$hires_image_path` not found. Rebuild with ",
         "CreateVisiumObjects(..., image_backend = 'deferred') to enable ",
         "on-demand hires image loading.")
  }
  if (!file.exists(hires_path)) {
    stop("Stashed hires image path no longer exists on disk: ", hires_path)
  }

  if (!isTRUE(attach)) {
    if (!requireNamespace("png", quietly = TRUE)) {
      stop("Package 'png' is required to decode '", hires_path, "'.")
    }
    message("--- Decoding hires image: ", hires_path, " ---")
    return(png::readPNG(hires_path))
  }

  # Rebuild in place, preserving the currently-attached pixel-backed image's
  # own name/@assay/@key -- these may no longer be "image"/"Spatial"/
  # "slice1" (CreateVisiumObjects()'s original defaults) if
  # RenameSpatialImages() has run since construction. Same reasoning
  # DropSpatialImage(mode = "downgrade") already applies per-image; hard-
  # coding the old defaults here would silently attach a stray new
  # obj@images$image entry instead of updating the real, renamed one.
  img_names <- names(obj@images)
  is_pixel <- vapply(img_names, function(nm) {
    arr <- tryCatch(Seurat::GetImage(obj@images[[nm]], mode = "raw"),
                    error = function(e) NULL)
    !is.null(arr) && length(dim(arr)) >= 2L
  }, logical(1))
  pixel_names <- img_names[is_pixel]
  if (length(pixel_names) != 1L) {
    stop("Expected exactly one pixel-backed image in obj@images to rebuild ",
         "at full resolution (to match obj@misc$hires_image_path), found ",
         length(pixel_names), ". ",
         if (length(pixel_names) == 0L) {
           "No pixel-backed image is currently attached."
         } else {
           paste0("Ambiguous which one to rebuild: ",
                  paste(pixel_names, collapse = ", "), ".")
         })
  }
  img_name <- pixel_names
  old <- obj@images[[img_name]]

  message("--- Rebuilding obj@images$", img_name, " at full resolution from ", hires_path, " ---")
  vis.image <- Seurat::Read10X_Image(
    image.dir  = dirname(hires_path),
    image.name = basename(hires_path),
    filter.matrix = TRUE
  )
  vis.image@assay <- old@assay
  vis.image@key   <- old@key
  obj@images[[img_name]] <- vis.image
  obj
}
