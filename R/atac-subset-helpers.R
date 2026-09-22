# ============================================================================
# Shared internal helper for safely cell-subsetting a Seurat object whose
# active assay is a Signac ChromatinAssay.
#
# Background: Signac's own `subset.ChromatinAssay()` (the S3 method that
# `subset.Seurat()` -> `base::subset(x[[assay]], ...)` dispatches to for any
# ChromatinAssay-backed assay) calls `GetAssayData(object, slot = "bias")`
# and `GetAssayData(object, slot = "positionEnrichment")` internally. As of
# SeuratObject 5.0.0 the `slot` argument of GetAssayData() is deprecated in
# favor of `layer`, and in current SeuratObject releases (>= 5.4.0,
# confirmed against a real HPC install) it is fully defunct -- that call now
# throws a hard error instead of a warning. SeuratObject's own
# `subset.Seurat()` wraps that call in a tryCatch whose error handler checks
# `e$message %in% c("Cannot find features provided", "None of the features
# provided found in this assay")` -- but `e$message` for this particular
# lifecycle/cli-formatted "defunct argument" error is a multi-element
# character vector (not a single string), so that `%in%` comparison itself
# returns a vector, and the `if()` around it crashes with "the condition has
# length > 1" -- masking the real error entirely.
#
# None of this is fixable from SingleCellTools by changing Signac or
# SeuratObject (confirmed on the real HPC install: both `install.packages`
# and a GitHub `develop`-branch reinstall of Signac fail to even compile
# under the available toolchain, independent of Signac version -- so a
# Signac upgrade isn't a viable near-term fix either). `bias` and
# `positionEnrichment` aren't Assay5 "layers" at all (they're not
# counts/data/scale.data-style data) -- they're ChromatinAssay's own custom
# S4 slots, and the old `GetAssayData(slot = ...)` was really just a
# generic S4-slot accessor for them. There's no `layer =` equivalent to
# swap in, but there doesn't need to be: reading those slots directly
# (`methods::slot(x, "bias")` / `methods::slot(x, "positionEnrichment")`)
# works under any SeuratObject version and sidesteps the whole
# slot-vs-layer deprecation.
#
# .subset_chromatin_assay() below is a faithful copy of
# Signac:::subset.ChromatinAssay (confirmed via
# `getS3method("subset", "ChromatinAssay")` against the installed Signac
# 1.9.0.9000), with only those two calls changed, plus its internal-only
# `SetIfNull()` helper inlined (Signac doesn't export it). Everything else
# -- the standard-Assay subset, FindTopFeatures(), ranges/motifs/fragments
# handling, seqinfo/annotation carry-through -- is unchanged from Signac's
# own logic. .subset_atac_seurat() uses it to rebuild a whole single-assay
# ATAC Seurat object, bypassing subset.Seurat()'s broken dispatch chain
# entirely. .safe_subset_cells() is the dispatcher used by callers: it
# routes ChromatinAssay-backed objects through this workaround and leaves
# every other object type on ordinary subset() (unaffected by this bug).
#
# NOTE: the very first line of .subset_chromatin_assay() (`standardassay <-
# subset(x = standardassay, ...)`, subsetting the *plain* Assay produced by
# `as(x, "Assay")`) is left calling ordinary subset() -- that's
# SeuratObject's core classic-Assay subsetting, not Signac's
# ChromatinAssay-specific code, and isn't known to have this bug. If a
# *different* error ever surfaces from that specific line, that would be a
# new, separate issue to chase, not this one.
# ============================================================================

#' @keywords internal
#' @noRd
.find_top_features_assay <- function(object, assay = NULL, min.cutoff = "q5",
                                      verbose = TRUE, ...) {
  # Faithful copy of Signac:::FindTopFeatures.Assay (confirmed via
  # getS3method("FindTopFeatures", "Assay") against the installed Signac),
  # with its one GetAssayData(object, slot = "counts") call -- now defunct
  # under current SeuratObject, same root cause as everywhere else in this
  # file, see the file header -- replaced by direct methods::slot() access.
  # `counts` is a standard slot on the classic (non-Assay5) `Assay` class
  # produced by `methods::as(x, "Assay")`, so this is an exact, safe
  # equivalent; everything below this line is unchanged from Signac's own
  # logic and operates on `data.use`, a plain matrix, not an Assay object,
  # so it never hits the broken GetAssayData() call.
  data.use <- methods::slot(object, "counts")
  if (Signac:::IsMatrixEmpty(x = data.use)) {
    if (verbose) {
      message("Count slot empty")
    }
    return(object)
  }
  hvf.info <- Signac::FindTopFeatures(object = data.use, assay = assay,
                                       min.cutoff = min.cutoff, verbose = verbose, ...)
  if (utils::packageVersion(pkg = "SeuratObject") < "4.9.9") {
    object[[names(x = hvf.info)]] <- hvf.info
  } else {
    object[names(x = hvf.info)] <- hvf.info
  }
  if (is.null(x = min.cutoff)) {
    SeuratObject::VariableFeatures(object = object) <- rownames(x = hvf.info)
  } else if (is.numeric(x = min.cutoff)) {
    SeuratObject::VariableFeatures(object = object) <- rownames(
      x = hvf.info[hvf.info$count > min.cutoff, ])
  } else if (is.na(x = min.cutoff)) {
    return(object)
  } else {
    percentile.use <- as.numeric(x = sub(pattern = "q", replacement = "",
                                          x = as.character(x = min.cutoff))) / 100
    SeuratObject::VariableFeatures(object = object) <- rownames(
      x = hvf.info[hvf.info$percentile > percentile.use, ])
  }
  return(object)
}

#' @keywords internal
#' @noRd
.subset_chromatin_assay <- function(x, features = NULL, cells = NULL) {
  standardassay <- methods::as(object = x, Class = "Assay")
  standardassay <- subset(x = standardassay, features = features, cells = cells)
  standardassay <- .find_top_features_assay(object = standardassay,
                                             min.cutoff = NA, verbose = FALSE)
  ranges.keep <- GenomicRanges::granges(x = x)
  if (!is.null(x = features)) {
    idx.keep <- rownames(x = x) %in% features
    ranges.keep <- ranges.keep[idx.keep]
  }
  motifs <- Signac::Motifs(object = x)
  if (!is.null(x = motifs)) {
    motifs <- subset(x = motifs, features = features)
  }
  cells <- if (is.null(x = cells)) colnames(x = x) else cells

  # Direct slot access instead of GetAssayData(x, slot = "positionEnrichment")
  # -- see file header.
  posmat <- methods::slot(x, "positionEnrichment")
  for (i in seq_along(along.with = posmat)) {
    posmat[[i]] <- posmat[[i]][cells, ]
  }

  frags <- Signac::Fragments(object = x)
  for (i in seq_along(along.with = frags)) {
    frags[[i]] <- subset(x = frags[[i]], cells = cells)
  }

  Signac::as.ChromatinAssay(
    x = standardassay,
    ranges = ranges.keep,
    seqinfo = GenomeInfoDb::seqinfo(x = x),
    annotation = Signac::Annotation(object = x),
    motifs = motifs,
    fragments = frags,
    # Direct slot access instead of GetAssayData(x, slot = "bias") -- see
    # file header.
    bias = methods::slot(x, "bias"),
    positionEnrichment = posmat
  )
}

#' @keywords internal
#' @noRd
.subset_atac_seurat <- function(so, cells) {
  assays_present <- SeuratObject::Assays(so)
  if (length(assays_present) > 1) {
    stop(
      "Safe ATAC-aware cell subsetting (the workaround for Signac's ",
      "subset.ChromatinAssay() calling the now-defunct ",
      "GetAssayData(slot = ...) -- see .subset_chromatin_assay()) only ",
      "supports a single-assay Seurat object. Found ",
      length(assays_present), " assays: ",
      paste(assays_present, collapse = ", "), ".")
  }
  a <- SeuratObject::DefaultAssay(so)
  new_assay <- .subset_chromatin_assay(so[[a]], cells = cells)
  meta <- so@meta.data[cells, , drop = FALSE]
  Seurat::CreateSeuratObject(counts = new_assay, assay = a,
                             meta.data = meta, project = so@project.name)
}

#' @keywords internal
#' @noRd
.safe_subset_cells <- function(so, cells) {
  a <- SeuratObject::DefaultAssay(so)
  if (methods::is(so[[a]], "ChromatinAssay")) {
    .subset_atac_seurat(so, cells)
  } else {
    subset(so, cells = cells)
  }
}
