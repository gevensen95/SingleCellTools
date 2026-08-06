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

  message("--- Rebuilding obj@images$image at full resolution from ", hires_path, " ---")
  vis.image <- Seurat::Read10X_Image(
    image.dir  = dirname(hires_path),
    image.name = basename(hires_path),
    filter.matrix = TRUE
  )
  vis.image@assay <- 'Spatial'
  vis.image@key   <- 'slice1'
  obj@images$image <- vis.image
  obj
}
