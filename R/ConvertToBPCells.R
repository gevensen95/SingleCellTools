#' Move a Seurat Assay Layer to an On-Disk BPCells Matrix
#'
#' Writes one or more layers of a Seurat assay (typically \code{"counts"})
#' to disk via \href{https://bnprks.github.io/BPCells/}{BPCells} and swaps
#' the in-memory matrix for the on-disk one. Because Seurat v5's
#' \code{Assay5}/\code{Layers} system accepts any matrix-like backend that
#' implements the necessary generics -- which is exactly what BPCells'
#' \code{IterableMatrix} does -- this doesn't require a custom object type;
#' native Seurat functions (\code{NormalizeData}, \code{FindVariableFeatures},
#' \code{ScaleData}, \code{RunPCA}, plotting, ...) keep working unmodified on
#' the returned object, just against a matrix that streams from disk instead
#' of living fully in RAM.
#'
#' Two important caveats:
#' \itemize{
#'   \item BPCells matrices are effectively immutable/lazy -- operations
#'     build a new lazy transform pipeline rather than mutating in place.
#'     This is fine for the normal analysis pipeline, but any code that
#'     forces \code{as.matrix()} on the \emph{whole} assay (rather than a
#'     deliberately-subsetted piece of it) will still fully materialize it
#'     into RAM, silently defeating the point of converting it in the first
#'     place.
#'   \item BPCells' compiled backend has historically had weaker/unofficial
#'     Windows support. It's Suggests-only here specifically so that's a
#'     choice, not a hard requirement, for everyone using this package.
#'   \item Signac \code{ChromatinAssay} objects (from
#'     \code{\link{CreateATACObjects}}/\code{\link{CreateATACObjectsFilter}})
#'     aren't supported -- they carry fragments/ranges/annotation slots a
#'     generic \code{Assay5} coercion can't preserve, so this errors clearly
#'     rather than attempting it. This is also why those two functions don't
#'     have an \code{on_disk} argument.
#' }
#'
#' @param obj A Seurat object, or a (optionally named) list of them.
#' @param assay Assay to convert. Default \code{NULL} uses
#'   \code{DefaultAssay(obj)} (resolved separately per object if \code{obj}
#'   is a list).
#' @param layers Character vector of layer name(s) to convert. Default
#'   \code{"counts"}. A layer not present on a given object's assay is
#'   skipped with a warning rather than an error, so one call can cover a
#'   mixed list.
#' @param path Directory to write the on-disk matrix/matrices to. A single
#'   object gets one subdirectory per layer (\code{path/<layer>}); a list of
#'   objects also gets one subdirectory per object first
#'   (\code{path/<name>/<layer>}), named from \code{names(obj)} if present,
#'   otherwise \code{Project(obj)}. Created if it doesn't exist. Default
#'   \code{file.path(getwd(), "bpcells")}.
#' @param overwrite Logical; if a target layer's directory already exists,
#'   error unless \code{TRUE} (in which case it's deleted and rewritten).
#'   Default \code{FALSE} -- protects an existing on-disk matrix from being
#'   silently discarded by a repeat call.
#' @return The updated Seurat object, or list of them -- matches the shape
#'   of \code{obj} (a single object in, a single object out; a list in, a
#'   list out with the original names preserved).
#' @export
ConvertToBPCells <- function(obj, assay = NULL, layers = "counts",
                             path = file.path(getwd(), "bpcells"),
                             overwrite = FALSE) {
  if (!requireNamespace("BPCells", quietly = TRUE)) {
    stop("Package 'BPCells' is required for ConvertToBPCells(). Install with: ",
         "remotes::install_github('bnprks/BPCells/r')")
  }

  parsed      <- .as_seurat_list(obj)
  objs        <- parsed$objs
  was_single  <- parsed$was_single
  orig_names  <- names(objs)  # NULL if `obj` was a single object or an unnamed list

  # Directory names used for a list of objects -- prefer the list's own
  # names, falling back to each object's Project() when unnamed (matches
  # how these objects are typically created, e.g. CreateRNAObjects() names
  # its list from `object_names`/basename(data_dirs)).
  dir_names <- orig_names
  if (is.null(dir_names)) {
    dir_names <- vapply(objs, function(o) SeuratObject::Project(o), character(1))
  }
  blank <- !nzchar(dir_names)
  if (any(blank)) dir_names[blank] <- paste0("sample", which(blank))

  message(sprintf('--- Converting %d object(s) to on-disk BPCells matrices (path = %s) ---',
                  length(objs), path))

  # ---- Pass 1: validate everything before writing anything to disk --------
  # Assay existence, the ChromatinAssay check, and "target directory already
  # exists without overwrite = TRUE" are all checked across every object x
  # layer combination here, before any BPCells::write_matrix_dir() call.
  # Without this, a failure partway through a many-object list (say, object
  # 15 of 20) would abort the whole call *after* objects 1-14 had already
  # had their (potentially many-GB) matrices written to disk -- work the
  # caller gets no object reference to and has to notice/clean up manually.
  plan <- lapply(seq_along(objs), function(i) {
    o   <- objs[[i]]
    a   <- if (is.null(assay)) SeuratObject::DefaultAssay(o) else assay
    tag <- if (length(objs) > 1) paste0(" ('", dir_names[i], "')") else ""
    if (!a %in% names(o@assays)) {
      stop("Assay '", a, "' not found in object", tag, ".")
    }
    if (inherits(o[[a]], "ChromatinAssay")) {
      stop("Assay '", a, "'", tag, " is a Signac `ChromatinAssay`, not a ",
           "standard Assay -- it carries fragments/ranges/annotation slots ",
           "that a generic Assay5 coercion doesn't know how to preserve, so ",
           "ConvertToBPCells() doesn't support it. This is why `on_disk` ",
           "isn't offered on CreateATACObjects()/CreateATACObjectsFilter().")
    }
    avail_layers <- SeuratObject::Layers(o[[a]])
    obj_path     <- if (length(objs) > 1) file.path(path, dir_names[i]) else path
    ly_paths     <- stats::setNames(file.path(obj_path, layers), layers)

    for (ly in layers) {
      if (!ly %in% avail_layers) next  # warned about (not fatal) in pass 2
      ly_path <- ly_paths[[ly]]
      if (dir.exists(ly_path) && !isTRUE(overwrite)) {
        stop("'", ly_path, "' already exists. Set overwrite = TRUE to ",
             "replace it, or choose a different `path`.")
      }
    }
    list(a = a, tag = tag, avail_layers = avail_layers, ly_paths = ly_paths)
  })

  # ---- Pass 2: do the actual conversion/writing ----------------------------
  objs <- lapply(seq_along(objs), function(i) {
    o    <- objs[[i]]
    info <- plan[[i]]
    a    <- info$a
    tag  <- info$tag

    if (!inherits(o[[a]], "Assay5")) {
      o[[a]] <- methods::as(o[[a]], Class = "Assay5")
    }

    for (ly in layers) {
      if (!ly %in% info$avail_layers) {
        warning("Layer '", ly, "' not found in assay '", a, "'", tag,
                "; skipping.", call. = FALSE)
        next
      }
      ly_path <- info$ly_paths[[ly]]
      if (dir.exists(ly_path)) {
        unlink(ly_path, recursive = TRUE)
      }
      # Only the parent is created here -- BPCells::write_matrix_dir() below
      # creates `ly_path` itself, and (per its own convention) expects that
      # leaf directory to not already exist.
      dir.create(dirname(ly_path), recursive = TRUE, showWarnings = FALSE)

      # BPCells stores genes x cells sparse matrices on disk -- the same
      # orientation Seurat already uses, so no transpose is needed either
      # direction.
      mat <- SeuratObject::LayerData(o, assay = a, layer = ly)
      BPCells::write_matrix_dir(mat, ly_path)
      on_disk_mat <- BPCells::open_matrix_dir(ly_path)
      SeuratObject::LayerData(o, assay = a, layer = ly) <- on_disk_mat

      message(sprintf("  %s%s layer '%s' -> %s", if (nzchar(tag)) paste0(dir_names[i], ": ") else "",
                      a, ly, ly_path))
    }
    o
  })

  if (was_single) return(objs[[1]])
  setNames(objs, orig_names)
}
