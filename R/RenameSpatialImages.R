#' Rename the images attached to a spatial Seurat object
#'
#' \code{names(obj@images) <- new_names} is the standard way to rename
#' spatial images, but a plain positional rename is fragile: nothing
#' guarantees \code{new_names} (e.g. \code{unique(obj$orig.ident)}) is in
#' the same order \code{obj@images} was built in, so a naive rename can
#' silently swap which name ends up attached to which image's actual data.
#' \code{RenameSpatialImages()} avoids that either by deriving each image's
#' new name directly from its own cells (\code{group_col}) or by validating
#' an explicit, caller-supplied mapping (\code{new_names}) before applying it.
#'
#' @details
#' Exactly one of \code{group_col} or \code{new_names} must be supplied.
#'
#' \strong{\code{group_col} (auto-derive).} For each image, looks up which
#' value of \code{obj@meta.data[[group_col]]} its own cells actually carry
#' (via \code{\link[SeuratObject]{Cells}}, not position) and uses that as
#' the new name. Errors clearly, per image, if: the image has no cells
#' present in \code{obj@meta.data}; every cell's \code{group_col} value is
#' \code{NA}; or the image's cells span more than one distinct value (can't
#' derive a single name automatically -- use \code{new_names} instead).
#'
#' \strong{\code{new_names} (manual).} Either a plain character vector the
#' same length as \code{names(obj@images)}, matched by \strong{position}
#' (same order as \code{names(obj@images)}); or a named character vector
#' (\code{c(old_name = "new_name", ...)}) for a partial rename -- images not
#' mentioned keep their current name. Use this when \code{group_col} isn't
#' the right axis, or there's no metadata column to derive names from at all.
#'
#' Either way, the resolved names are checked for collisions (two images
#' ending up with the same name) before anything is applied.
#'
#' @param obj A Seurat object with at least one image in \code{obj@images}.
#' @param group_col Metadata column to auto-derive each image's new name
#'   from (e.g. \code{"orig.ident"}). \code{NULL} by default.
#' @param new_names Explicit new name(s) -- plain vector (positional) or
#'   named vector (old -> new, partial rename allowed). \code{NULL} by
#'   default.
#' @return \code{obj} with \code{names(obj@images)} updated.
#' @examples
#' \dontrun{
#' # Auto-derive from a metadata column
#' visium <- RenameSpatialImages(visium, group_col = "orig.ident")
#'
#' # Manual, positional (must supply one name per image, in
#' # names(visium@images) order)
#' visium <- RenameSpatialImages(visium, new_names = c("anterior", "posterior"))
#'
#' # Manual, partial rename by name
#' visium <- RenameSpatialImages(visium, new_names = c(slice1 = "anterior"))
#' }
#' @export
RenameSpatialImages <- function(obj, group_col = NULL, new_names = NULL) {
  .assert_seurat(obj)
  if (length(obj@images) == 0) {
    stop("`obj` has no images in obj@images -- nothing to rename.")
  }
  if (is.null(group_col) && is.null(new_names)) {
    stop("Exactly one of `group_col` or `new_names` must be supplied -- both are NULL.")
  }
  if (!is.null(group_col) && !is.null(new_names)) {
    stop("Exactly one of `group_col` or `new_names` must be supplied -- both were given.")
  }

  current_names <- names(obj@images)

  if (!is.null(new_names)) {
    if (!is.character(new_names)) {
      stop("`new_names` must be a character vector.")
    }
    if (!is.null(names(new_names))) {
      missing_images <- setdiff(names(new_names), current_names)
      if (length(missing_images)) {
        stop("`new_names` has name(s) not found in obj@images: ",
             paste(missing_images, collapse = ", "), ". Available: ",
             paste(current_names, collapse = ", "))
      }
      resolved <- current_names
      resolved[match(names(new_names), current_names)] <- unname(new_names)
    } else {
      if (length(new_names) != length(current_names)) {
        stop("`new_names` has length ", length(new_names), " but obj@images has ",
             length(current_names), " image(s) -- supply one name per image ",
             "(in the same order as names(obj@images)), or a named vector for ",
             "a partial rename.")
      }
      resolved <- new_names
    }
  } else {
    if (!group_col %in% colnames(obj@meta.data)) {
      stop("`group_col` column '", group_col, "' not found in obj@meta.data.")
    }
    group_vals <- obj@meta.data[[group_col]]
    names(group_vals) <- rownames(obj@meta.data)

    resolved <- vapply(current_names, function(img) {
      cells_in_img <- intersect(SeuratObject::Cells(obj[[img]]), names(group_vals))
      if (length(cells_in_img) == 0) {
        stop("Image '", img, "' has no cells present in obj@meta.data -- can't ",
             "derive a name from `group_col`.")
      }
      vals <- unique(as.character(group_vals[cells_in_img]))
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0) {
        stop("Image '", img, "': every cell has NA for `group_col` '", group_col,
             "' -- can't derive a name.")
      }
      if (length(vals) != 1) {
        stop("Image '", img, "' spans ", length(vals), " distinct value(s) of '",
             group_col, "' (", paste(vals, collapse = ", "), ") -- can't derive ",
             "a single name automatically. Use `new_names` instead.")
      }
      vals
    }, character(1), USE.NAMES = FALSE)
  }

  dup <- unique(resolved[duplicated(resolved)])
  if (length(dup)) {
    stop("Resolved name(s) would collide -- more than one image would be named: ",
         paste(dup, collapse = ", "))
  }

  names(obj@images) <- resolved
  obj
}
