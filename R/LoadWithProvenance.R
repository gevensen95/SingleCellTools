#' Load a Seurat object saved by SaveWithProvenance, safely
#'
#' Companion to \code{\link{SaveWithProvenance}}: reads the \code{.rds}, and
#' if a \code{<name>_provenance.json} sidecar sits next to it, uses it to
#' report what was in the file (\code{saved_at}, package versions at save
#' time, git state) without needing to load the object first.
#'
#' \strong{Why this exists.} SeuratObject occasionally adds slots to its S4
#' classes between versions (recent examples: \code{misc} added to
#' \code{SpatialImage} in 5.2.0, \code{coords_x_orientation} added to
#' \code{FOV} in 5.3.0). An object serialized under an older SeuratObject
#' doesn't have those slots, and \code{methods::validObject()} -- which
#' Seurat/SeuratObject call constantly, deep inside otherwise-unrelated
#' functions -- then fails with a cryptic \code{invalid class ... slots in
#' class definition but not in object} error, often long after the object
#' was actually loaded (whatever function happens to call
#' \code{validObject()} first "owns" the confusing error, even though it
#' didn't cause it). This function catches that at load time instead: it
#' validates the object right after \code{readRDS()}, and if that fails,
#' runs \code{UpdateSeuratObject()} automatically.
#'
#' \code{UpdateSeuratObject()}'s own \code{coords_x_orientation} handling is
#' documented as Visium-specific, so it isn't guaranteed to reach every
#' \code{FOV} object on non-Visium spatial platforms (Xenium/CosMx/MERFISH).
#' If the object is still invalid after \code{UpdateSeuratObject()}, each
#' attached image is re-initialized against its own, currently-loaded class
#' definition (\code{methods::new(class(image), image)}, which copies
#' existing slot values across and fills any newly-added slots with that
#' class's prototype defaults) as a fallback. If it's \emph{still} invalid
#' after both, a warning says so rather than silently handing back a broken
#' object.
#'
#' @param file Path to a \code{.rds} file, typically one written by
#'   \code{\link{SaveWithProvenance}}. Works on a plain \code{saveRDS()}
#'   file too (with no sidecar, or one that isn't a Seurat object) -- the
#'   provenance summary and update logic are skipped, respectively, rather
#'   than erroring.
#' @param update If \code{TRUE} (default) and the loaded object fails its
#'   validity check, attempt to fix it automatically (see Details). If
#'   \code{FALSE}, only warn.
#' @param verbose Message the provenance summary and any update steps
#'   taken. Default \code{TRUE}.
#' @return The loaded object. For a Seurat object, the parsed provenance
#'   sidecar (if found) is attached as \code{attr(obj, "provenance")} --
#'   best-effort only; ordinary Seurat operations (subset, merge, ...) are
#'   not guaranteed to preserve arbitrary attributes.
#' @examples
#' \dontrun{
#' obj <- LoadWithProvenance("results/obj_annotated.rds")
#' attr(obj, "provenance")$seurat_state$n_cells
#' }
#' @export
LoadWithProvenance <- function(file, update = TRUE, verbose = TRUE) {

  if (!file.exists(file)) stop("File not found: ", file)

  sidecar <- sub("\\.rds$", "_provenance.json", file, ignore.case = TRUE)
  prov <- NULL
  if (file.exists(sidecar)) {
    if (requireNamespace("jsonlite", quietly = TRUE)) {
      prov <- tryCatch(
        jsonlite::fromJSON(sidecar, simplifyVector = TRUE),
        error = function(e) {
          warning("Could not parse provenance sidecar '", sidecar, "': ",
                  conditionMessage(e))
          NULL
        }
      )
    } else if (isTRUE(verbose)) {
      message("'jsonlite' not installed -- loading without reading the ",
              "provenance sidecar.")
    }
  } else if (isTRUE(verbose)) {
    message("No provenance sidecar found at ", sidecar, ".")
  }

  obj <- readRDS(file)

  if (inherits(obj, "Seurat")) {
    saved_ver <- tryCatch(prov$package_versions$SeuratObject,
                          error = function(e) NULL)
    obj <- .update_if_invalid(obj, update = update, verbose = verbose,
                              saved_ver = saved_ver)
  } else if (isTRUE(verbose)) {
    message("Loaded object is not a Seurat object -- skipping the ",
            "validity/UpdateSeuratObject check.")
  }

  if (!is.null(prov) && isTRUE(verbose)) {
    message(sprintf(
      "Loaded %s (saved %s%s)", file,
      if (!is.null(prov$saved_at)) prov$saved_at else "unknown time",
      if (isTRUE(prov$git$dirty)) ", from a dirty git state" else ""
    ))
  }

  if (!is.null(prov)) attr(obj, "provenance") <- prov
  obj
}

#' @keywords internal
#' @noRd
.update_if_invalid <- function(obj, update, verbose, saved_ver) {
  valid <- tryCatch({ methods::validObject(obj); TRUE },
                    error = function(e) FALSE)
  if (valid) return(obj)

  current_ver <- tryCatch(as.character(utils::packageVersion("SeuratObject")),
                          error = function(e) NA_character_)
  ver_msg <- sprintf(
    "saved under SeuratObject %s, currently running %s",
    if (is.null(saved_ver) || is.na(saved_ver)) "an unknown version" else saved_ver,
    if (is.na(current_ver)) "an unknown version" else current_ver)

  if (!isTRUE(update)) {
    warning("Loaded object fails Seurat/SeuratObject validity checks (", ver_msg,
            "). Pass update = TRUE (default) to attempt an automatic fix, or ",
            "run UpdateSeuratObject() yourself.")
    return(obj)
  }

  if (isTRUE(verbose)) {
    message("Object fails validity checks (", ver_msg,
            ") -- running UpdateSeuratObject()...")
  }
  obj <- tryCatch(
    Seurat::UpdateSeuratObject(obj),
    error = function(e) {
      warning("UpdateSeuratObject() failed: ", conditionMessage(e))
      obj
    }
  )

  valid <- tryCatch({ methods::validObject(obj); TRUE },
                    error = function(e) FALSE)
  if (!valid) {
    # UpdateSeuratObject()'s own new-slot handling (e.g. coords_x_orientation,
    # added in SeuratObject 5.3.0) is documented as Visium-specific, so it
    # doesn't necessarily reach every FOV object on non-Visium spatial
    # platforms (Xenium/CosMx/MERFISH). Fall back to re-initializing each
    # image against its own, currently-loaded class definition --
    # methods::new(Class, sourceObject) copies over slots that already exist
    # on the source and fills any newly-added ones with that class's
    # prototype defaults, which is what plain methods::as() does NOT do when
    # the source is already (nominally) of the target class.
    imgs <- tryCatch(SeuratObject::Images(obj), error = function(e) character(0))
    for (img in imgs) {
      obj[[img]] <- tryCatch(
        methods::new(class(obj[[img]])[1], obj[[img]]),
        error = function(e) obj[[img]]
      )
    }
    valid <- tryCatch({ methods::validObject(obj); TRUE },
                      error = function(e) FALSE)
  }

  if (!valid) {
    warning("Object still fails validity checks after UpdateSeuratObject() ",
            "and image re-initialization -- inspect manually (see ",
            "?UpdateSeuratObject).")
  }
  obj
}
