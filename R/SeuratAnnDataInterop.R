#' Convert a Seurat object to an AnnData (h5ad) file
#'
#' Round-trips data with the Python / scanpy world via the
#' \code{zellkonverter} Bioconductor package. The Seurat object is first
#' converted to a \code{SingleCellExperiment} (using SeuratObject's
#' \code{as.SingleCellExperiment}) and then written as \code{.h5ad}.
#'
#' @param obj A Seurat object.
#' @param file Output path. Should end in \code{.h5ad}.
#' @param assay Assay to export. Default DefaultAssay.
#' @return Invisibly, the path written to.
#' @importFrom Seurat as.SingleCellExperiment DefaultAssay
#' @export
ToAnnData <- function(obj, file, assay = NULL) {
  if (!requireNamespace("zellkonverter", quietly = TRUE)) {
    stop("Package 'zellkonverter' is required. Install with: ",
         "BiocManager::install('zellkonverter')")
  }
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (!grepl("\\.h5ad$", file)) {
    warning("Output file does not end in '.h5ad'; writing anyway.")
  }
  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay

  message("--- Converting Seurat -> SingleCellExperiment ---")
  sce <- Seurat::as.SingleCellExperiment(obj, assay = a)

  message(sprintf("--- Writing AnnData to %s ---", file))
  zellkonverter::writeH5AD(sce, file = file)
  invisible(file)
}


#' Read an AnnData (h5ad) file into a Seurat object
#'
#' Reads the file as a \code{SingleCellExperiment} via \code{zellkonverter}
#' and converts to a Seurat object with \code{Seurat::as.Seurat}. Note that
#' some scanpy-specific bits (e.g. \code{obsm} entries beyond UMAP/PCA,
#' \code{varm}, \code{uns}) may not round-trip cleanly.
#'
#' After reading, each assay is coerced to a column-compressed sparse
#' matrix (\code{dgCMatrix}) with a double-precision \code{@x} slot. This
#' avoids two common errors with h5ads written by scanpy / scvi-tools:
#' \itemize{
#'   \item \code{invalid class "dgRMatrix" object: 'x' slot is not of type "double"}
#'     — counts stored as integers in the source file.
#'   \item Seurat refusing to construct an assay from a row-compressed
#'     sparse matrix.
#' }
#'
#' @param file Path to an \code{.h5ad} file.
#' @param counts Name of the layer in the AnnData object to treat as raw
#'   counts. Default \code{"counts"}; falls back to \code{"X"} if not found.
#' @param data Name of the layer to treat as normalized data. Default
#'   \code{"logcounts"}.
#' @param reader Which zellkonverter backend to use. \code{"python"}
#'   (default) goes through the \code{anndata} Python module via basilisk
#'   — required for files where counts are stored as integer sparse
#'   matrices (e.g. Tabula Sapiens, scvi-tools exports), because the
#'   native R reader builds an invalid \code{dgRMatrix} with integer
#'   \code{@x} slot. \code{"R"} is the native reader; faster startup but
#'   may fail on integer-typed sparse data.
#' @return A Seurat object.
#' @importFrom Seurat as.Seurat
#' @importFrom methods as
#' @export
FromAnnData <- function(file,
                        counts = "counts",
                        data   = "logcounts",
                        reader = c("python", "R")) {
  reader <- match.arg(reader)
  if (!requireNamespace("zellkonverter", quietly = TRUE)) {
    stop("Package 'zellkonverter' is required. Install with: ",
         "BiocManager::install('zellkonverter')")
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required.")
  }
  if (!file.exists(file)) stop("File not found: ", file)

  message(sprintf("--- Reading AnnData from %s (%s reader) ---", file, reader))
  sce <- tryCatch(
    zellkonverter::readH5AD(file, reader = reader),
    error = function(e) {
      if (reader == "python") {
        message("Python reader failed (", conditionMessage(e), ").")
        message("Falling back to native R reader.")
        zellkonverter::readH5AD(file, reader = "R")
      } else {
        stop(e)
      }
    }
  )

  available_assays <- SummarizedExperiment::assayNames(sce)
  if (!counts %in% available_assays) {
    if ("X" %in% available_assays) {
      counts <- "X"
      message("  No 'counts' assay; using 'X' as counts.")
    } else {
      stop("Neither '", counts, "' nor 'X' found in the AnnData assays. ",
           "Available: ", paste(available_assays, collapse = ", "))
    }
  }
  if (!data %in% available_assays) data <- counts

  # ---- Sanitize each assay -----------------------------------------------
  # 1) Convert dgRMatrix (row-compressed) to dgCMatrix (col-compressed),
  #    which is the format Seurat expects.
  # 2) Coerce the @x slot to double — newer Matrix versions strictly
  #    require this; integer-stored counts (as in Tabula Sapiens) trip
  #    validObject otherwise.
  message("--- Sanitizing assays (dgCMatrix + double) ---")
  for (a in available_assays) {
    m <- SummarizedExperiment::assay(sce, a)
    if (inherits(m, "sparseMatrix")) {
      m <- methods::as(m, "CsparseMatrix")
      if (!is.double(m@x)) m@x <- as.numeric(m@x)
    } else if (is.matrix(m) && !is.double(m)) {
      storage.mode(m) <- "double"
    }
    SummarizedExperiment::assay(sce, a) <- m
  }

  message("--- Converting SingleCellExperiment -> Seurat ---")
  Seurat::as.Seurat(sce, counts = counts, data = data)
}
