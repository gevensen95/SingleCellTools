#' Register spatial coordinates across serial sections
#'
#' Aligns each sample's tissue coordinates onto a common reference frame
#' using user-identified corresponding landmark points (e.g. tissue
#' corners, fiducial marks, or the same anatomical structure annotated in
#' each section) -- automatic registration without any known correspondence
#' isn't attempted here, since it isn't reliable in general. Writes the
#' transformed coordinates as new metadata columns so aligned samples can
#' be overlaid or compared spot-for-spot.
#'
#' @param obj_list A (optionally named) list of >= 2 spatial Seurat objects.
#' @param reference_index Index (or name) into \code{obj_list} to treat as
#'   the fixed reference frame. Default 1.
#' @param landmarks A list, same length/order as \code{obj_list}, of
#'   n x 2 matrices/data frames (x, y) of corresponding landmark points --
#'   point \code{i} in every element must be the same physical location.
#'   Required.
#' @param method \code{"procrustes"} (default; requires the \code{vegan}
#'   package -- similarity transform: rotation + uniform scale +
#'   translation, minimizing landmark distance) or \code{"affine"} (no
#'   extra dependency -- a general 2D affine fit via least squares; use
#'   this if \code{vegan} isn't available, or if samples genuinely need
#'   independent x/y scaling or shear).
#' @param coord_prefix Prefix for the new metadata columns. Default
#'   \code{"aligned_"} (columns \code{aligned_x}, \code{aligned_y}).
#' @param image_name Image to read each object's tissue coordinates from.
#'   \code{NULL} (default) uses the first image on each object.
#' @param verbose Message progress. Default \code{TRUE}.
#' @return \code{obj_list}, each object with \code{<coord_prefix>x} /
#'   \code{<coord_prefix>y} metadata columns added (the reference object's
#'   coordinates copied through unchanged).
#' @examples
#' \dontrun{
#' lm <- list(
#'   matrix(c(100, 120, 400, 410, 90, 500, 95, 505), ncol = 2),  # sample 1 (reference)
#'   matrix(c(105, 130, 395, 420, 85, 495, 100, 510), ncol = 2)  # sample 2
#' )
#' aligned <- AlignSpatialSamples(list(s1 = visium1, s2 = visium2), landmarks = lm)
#' ggplot() +
#'   geom_point(data = aligned$s1[[]], aes(aligned_x, aligned_y), color = "blue") +
#'   geom_point(data = aligned$s2[[]], aes(aligned_x, aligned_y), color = "red")
#' }
#' @importFrom Seurat GetTissueCoordinates
#' @export
AlignSpatialSamples <- function(obj_list,
                                reference_index = 1,
                                landmarks       = NULL,
                                method          = c("procrustes", "affine"),
                                coord_prefix    = "aligned_",
                                image_name      = NULL,
                                verbose         = TRUE) {

  method <- match.arg(method)
  if (!is.list(obj_list) || length(obj_list) < 2 ||
      !all(vapply(obj_list, inherits, logical(1), "Seurat"))) {
    stop("`obj_list` must be a list of >= 2 Seurat objects.")
  }
  if (is.null(landmarks) || length(landmarks) != length(obj_list)) {
    stop("`landmarks` is required: a list of the same length as ",
         "`obj_list`, one n x 2 (x, y) matrix/data frame per object, ",
         "with matching points in the same row order across all of them.")
  }
  if (method == "procrustes" && !requireNamespace("vegan", quietly = TRUE)) {
    stop("'vegan' is required for method = 'procrustes'. Install with ",
         "install.packages('vegan'), or use method = 'affine' instead.")
  }

  ref_i <- reference_index
  ref_lm <- as.matrix(landmarks[[ref_i]])

  .get_coords <- function(obj) {
    im <- if (is.null(image_name)) names(obj@images)[1] else image_name
    coords <- as.data.frame(Seurat::GetTissueCoordinates(obj[[im]]))
    coord_cols <- intersect(c("x", "y"), colnames(coords))
    if (length(coord_cols) < 2) coord_cols <- colnames(coords)[1:2]
    if ("cell" %in% colnames(coords)) rownames(coords) <- coords$cell
    as.matrix(coords[, coord_cols])
  }

  .fit_affine <- function(src, dst) {
    # dst ~ [src, 1] %*% M  (M is 3x2: 2 linear params + translation, per axis)
    design <- cbind(src, 1)
    M <- solve(t(design) %*% design, t(design) %*% dst)
    function(pts) cbind(pts, 1) %*% M
  }

  out <- vector("list", length(obj_list))
  names(out) <- names(obj_list)

  for (i in seq_along(obj_list)) {
    obj <- obj_list[[i]]
    all_coords <- .get_coords(obj)

    if (i == ref_i) {
      transformed <- all_coords
      if (isTRUE(verbose)) message(sprintf("--- Sample %d: reference, no transform ---", i))
    } else {
      lm_i <- as.matrix(landmarks[[i]])
      if (nrow(lm_i) != nrow(ref_lm)) {
        stop("landmarks[[", i, "]] has ", nrow(lm_i), " point(s); reference ",
             "has ", nrow(ref_lm), ". Landmark counts must match.")
      }
      if (isTRUE(verbose)) {
        message(sprintf("--- Sample %d: aligning to reference (%s, %d landmarks) ---",
                        i, method, nrow(lm_i)))
      }
      if (method == "procrustes") {
        fit <- vegan::procrustes(X = ref_lm, Y = lm_i, scale = TRUE)
        transformed <- as.matrix(predict(fit, newdata = all_coords))
      } else {
        xform <- .fit_affine(lm_i, ref_lm)
        transformed <- xform(all_coords)
      }
    }

    common <- intersect(rownames(transformed), colnames(obj))
    ax <- setNames(rep(NA_real_, ncol(obj)), colnames(obj))
    ay <- ax
    ax[common] <- transformed[common, 1]
    ay[common] <- transformed[common, 2]
    obj@meta.data[[paste0(coord_prefix, "x")]] <- ax
    obj@meta.data[[paste0(coord_prefix, "y")]] <- ay

    out[[i]] <- obj
  }

  out
}
