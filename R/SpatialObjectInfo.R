#' Report Per-Sample Image/FOV Size and Memory Footprint for Spatial Objects
#'
#' Audits what's actually attached at \code{obj@images} across one or more
#' Seurat objects, for either kind of spatial image this package builds:
#' pixel-backed images (\code{VisiumV1}, from \code{\link{CreateVisiumObjects}})
#' and coordinate-only \code{FOV} images (Xenium/CosMx/MERFISH-style, e.g.
#' from \code{\link{LoadXenium2}}, \code{Seurat::LoadNanostring()}, or any
#' hand-built \code{Seurat::CreateFOV()} object). Which columns are
#' meaningful depends on which kind of image a given row is -- see Value.
#' Useful for seeing where spatial memory is going across a multi-sample
#' list before deciding whether to \code{\link{DropSpatialImage}}.
#'
#' Pixel-backed vs. coordinate-only is detected by capability
#' (\code{Seurat::GetImage(img, mode = "raw")} succeeding and returning an
#' array), not by class name -- so this also covers any other
#' \code{SpatialImage} subclass that behaves like one or the other, not
#' just \code{VisiumV1}/\code{FOV} specifically.
#'
#' @param obj A Seurat object, or a (optionally named) list of them -- e.g.
#'   the list returned by \code{\link{CreateVisiumObjects}} or
#'   \code{\link{LoadXenium2}}.
#' @return A data frame with one row per image (sample x image name):
#'   \describe{
#'     \item{sample, image_name, class}{Identify the row.}
#'     \item{n_cells}{Number of cells/spots the image indexes
#'       (\code{Seurat::Cells(img)}, falling back to the \code{centroids}
#'       boundary's cell count for \code{FOV} images where \code{Cells()}
#'       doesn't resolve). \code{NA} if neither works.}
#'     \item{width, height, size_mb}{Pixel dimensions and in-memory size of
#'       the decoded array. \code{NA} for coordinate-only images (nothing
#'       decoded).}
#'     \item{deferred, hires_path}{Visium-specific: \code{TRUE}/path if
#'       \code{obj@misc$hires_image_path} is stashed (lowres attached,
#'       hires available on demand). \code{NA} for coordinate-only images
#'       -- the eager/deferred distinction doesn't apply to them.}
#'     \item{boundary_sets}{FOV-specific: comma-joined names of attached
#'       boundary sets (e.g. \code{"centroids, segmentation"}, from
#'       \code{Seurat::Boundaries()}). \code{NA} for pixel-backed images.}
#'     \item{has_molecules, n_molecule_features, n_molecules}{FOV-specific,
#'       best-effort: whether a \code{molecules} boundary is attached, how
#'       many features (genes) have molecule data, and the total transcript
#'       count across them. \code{NA}/\code{FALSE} for pixel-backed images
#'       or if introspection fails.}
#'     \item{molecules_lazy}{Object-level (repeated per row for that
#'       sample): \code{TRUE} if \code{obj@misc$molecules_lazy} is stashed,
#'       i.e. an \code{arrow} Dataset handle exists for re-querying the full
#'       (pre-QV-filter) transcript table -- see
#'       \code{\link{QueryXeniumMolecules}}. This is independent of
#'       \code{has_molecules}: the molecules already attached to the FOV are
#'       the QV-filtered subset either way, in-memory or not.}
#'   }
#'   Samples with no images attached get a single row of \code{NA}s past
#'   \code{sample}.
#' @export
SpatialObjectInfo <- function(obj) {
  objs <- .as_seurat_list(obj)$objs
  if (is.null(names(objs)) || any(!nzchar(names(objs)))) {
    names(objs) <- paste0("sample", seq_along(objs))
  }

  .empty_row <- function(nm) {
    data.frame(sample = nm, image_name = NA_character_, class = NA_character_,
              n_cells = NA_integer_, width = NA_integer_, height = NA_integer_,
              size_mb = NA_real_, deferred = NA, hires_path = NA_character_,
              boundary_sets = NA_character_, has_molecules = NA,
              n_molecule_features = NA_integer_, n_molecules = NA_integer_,
              molecules_lazy = NA, stringsAsFactors = FALSE)
  }

  .n_cells <- function(img) {
    n <- tryCatch(length(SeuratObject::Cells(img)), error = function(e) NA_integer_)
    if (is.na(n)) {
      n <- tryCatch(length(img$centroids@cells), error = function(e) NA_integer_)
    }
    as.integer(n)
  }

  .boundary_sets <- function(img) {
    tryCatch({
      # Boundaries() is exported by SeuratObject, not Seurat -- confirmed
      # directly: Seurat::Boundaries(fov) errors with "'Boundaries' is not
      # an exported object from 'namespace:Seurat'", which is exactly why
      # this was silently swallowed to NA_character_ by the tryCatch below.
      b <- SeuratObject::Boundaries(img)
      if (length(b) == 0L) NA_character_ else paste(b, collapse = ", ")
    }, error = function(e) NA_character_)
  }

  .molecule_stats <- function(img) {
    mols <- tryCatch(img$molecules, error = function(e) NULL)
    if (is.null(mols) || length(mols) == 0L) {
      return(list(has = FALSE, n_features = NA_integer_, n_molecules = NA_integer_))
    }
    n_mol <- tryCatch(
      sum(vapply(mols, function(m) nrow(m@coords), integer(1))),
      error = function(e) NA_integer_
    )
    list(has = TRUE, n_features = length(mols), n_molecules = n_mol)
  }

  rows <- lapply(names(objs), function(nm) {
    o <- objs[[nm]]
    img_names <- names(o@images)
    if (length(img_names) == 0L) {
      return(.empty_row(nm))
    }
    lazy <- !is.null(o@misc$molecules_lazy)
    do.call(rbind, lapply(img_names, function(in_) {
      img <- o@images[[in_]]
      arr <- tryCatch(Seurat::GetImage(img, mode = "raw"), error = function(e) NULL)
      is_pixel <- !is.null(arr) && length(dim(arr)) >= 2L

      if (is_pixel) {
        dims <- dim(arr)[1:2]
        size_mb <- as.numeric(utils::object.size(arr)) / 1e6
        hires_path <- o@misc$hires_image_path
        data.frame(
          sample = nm, image_name = in_, class = class(img)[1],
          n_cells = .n_cells(img), width = dims[2], height = dims[1],
          size_mb = round(size_mb, 1), deferred = !is.null(hires_path),
          hires_path = if (is.null(hires_path)) NA_character_ else hires_path,
          boundary_sets = NA_character_, has_molecules = FALSE,
          n_molecule_features = NA_integer_, n_molecules = NA_integer_,
          molecules_lazy = lazy, stringsAsFactors = FALSE
        )
      } else {
        mol <- .molecule_stats(img)
        data.frame(
          sample = nm, image_name = in_, class = class(img)[1],
          n_cells = .n_cells(img), width = NA_integer_, height = NA_integer_,
          size_mb = NA_real_, deferred = NA, hires_path = NA_character_,
          boundary_sets = .boundary_sets(img), has_molecules = mol$has,
          n_molecule_features = mol$n_features, n_molecules = mol$n_molecules,
          molecules_lazy = lazy, stringsAsFactors = FALSE
        )
      }
    }))
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
