# Fail fast, with a clear message, when an optional (Suggests) package a code
# path needs isn't installed. Several packages Seurat functions rely on --
# fitdistrplus (HTODemux), spatstat.explore/spatstat.geom (markvariogram),
# lmtest (test.use = "LR"), ROCR, ... -- moved from Seurat's Imports to its
# Suggests in Seurat 6, so a fresh install no longer brings them along.
# `why` names what needs them, e.g. "DemultiplexHashtags(method = 'htodemux')".
#' @keywords internal
#' @noRd
.require_suggested <- function(pkgs, why) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    quoted <- paste0("'", missing, "'")
    stop(why, " needs the ", if (length(missing) > 1) "packages " else "package ",
         paste(quoted, collapse = ", "),
         " (not installed automatically with Seurat >= 6). Install with: ",
         "install.packages(c(", paste(quoted, collapse = ", "), "))",
         call. = FALSE)
  }
  invisible(TRUE)
}
