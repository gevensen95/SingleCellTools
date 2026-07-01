#' Save a Seurat object with a provenance sidecar
#'
#' Writes a Seurat object as \code{.rds} and, next to it, a
#' \code{<name>_provenance.json} sidecar recording the analysis state:
#' package versions, DefaultAssay, layers per assay, reductions,
#' metadata column names, cell / gene counts, and (optionally) the git
#' SHA / dirty flag of the calling script directory. Future-you (or a
#' reviewer) can look at the sidecar to understand what was in the RDS
#' without loading it into R.
#'
#' @param obj A Seurat object.
#' @param file Output \code{.rds} path. The sidecar is written to
#'   \code{sub(".rds$", "_provenance.json", file)}.
#' @param git_dir Optional path to a git repository whose HEAD to record
#'   in the sidecar. \code{NULL} (default) skips git.
#' @param extra Optional named list of additional key-value pairs to
#'   include in the sidecar under \code{extra}.
#' @return Invisible \code{file} path.
#' @examples
#' \dontrun{
#' SaveWithProvenance(obj, file = "results/obj_annotated.rds",
#'                    git_dir = getwd(),
#'                    extra = list(analyst = "K. Evensen",
#'                                 project = "study42"))
#' }
#' @importFrom Seurat DefaultAssay Reductions Assays Idents
#' @export
SaveWithProvenance <- function(obj,
                               file,
                               git_dir = NULL,
                               extra   = NULL) {

  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("'jsonlite' is required for the provenance sidecar.")
  }

  # Ensure output dir exists
  out_dir <- dirname(file)
  if (nzchar(out_dir) && !dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Save the object
  saveRDS(obj, file = file)

  # Assemble provenance
  layers_by_assay <- lapply(
    Seurat::Assays(obj),
    function(a) tryCatch(SeuratObject::Layers(obj[[a]]),
                         error = function(e) NULL)
  )
  names(layers_by_assay) <- Seurat::Assays(obj)

  # Package versions of key deps that are actually loaded
  pkgs <- c("Seurat", "SeuratObject", "Signac", "Matrix", "DESeq2",
            "harmony", "UCell", "SingleR", "SummarizedExperiment",
            "SingleCellExperiment")
  versions <- lapply(pkgs, function(p) {
    tryCatch(as.character(utils::packageVersion(p)),
             error = function(e) NA_character_)
  })
  names(versions) <- pkgs

  # Optional git info
  git_info <- NULL
  if (!is.null(git_dir) && dir.exists(git_dir)) {
    git_info <- tryCatch({
      sha <- suppressWarnings(system2(
        "git", c("-C", shQuote(git_dir), "rev-parse", "HEAD"),
        stdout = TRUE, stderr = FALSE))
      dirty <- suppressWarnings(system2(
        "git", c("-C", shQuote(git_dir), "status", "--porcelain"),
        stdout = TRUE, stderr = FALSE))
      list(sha   = if (length(sha) == 1L) sha else NA_character_,
           dirty = length(dirty) > 0L,
           dir   = normalizePath(git_dir, mustWork = FALSE))
    }, error = function(e) NULL)
  }

  prov <- list(
    saved_at        = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    file            = normalizePath(file, mustWork = FALSE),
    hostname        = Sys.info()[["nodename"]],
    R_version       = as.character(getRversion()),
    package_versions = versions,
    seurat_state    = list(
      n_cells       = ncol(obj),
      n_genes_by_assay = vapply(Seurat::Assays(obj),
                                function(a) tryCatch(nrow(obj[[a]]),
                                                     error = function(e) NA_integer_),
                                numeric(1)),
      default_assay = Seurat::DefaultAssay(obj),
      assays        = Seurat::Assays(obj),
      layers        = layers_by_assay,
      reductions    = Seurat::Reductions(obj),
      images        = names(obj@images),
      metadata_cols = colnames(obj@meta.data),
      idents_levels = levels(as.factor(Seurat::Idents(obj)))
    ),
    git             = git_info,
    extra           = extra
  )

  sidecar <- sub("\\.rds$", "_provenance.json", file, ignore.case = TRUE)
  jsonlite::write_json(prov, sidecar, pretty = TRUE, auto_unbox = TRUE,
                       null = "null")
  message(sprintf("Wrote %s and %s", file, sidecar))
  invisible(file)
}
