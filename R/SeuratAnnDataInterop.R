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
#' and converts to a Seurat object. After reading, each assay is coerced
#' to a column-compressed sparse matrix (\code{dgCMatrix}) with a
#' double-precision \code{@x} slot. This avoids two common errors with
#' h5ads written by scanpy / scvi-tools:
#' \itemize{
#'   \item \code{invalid class "dgRMatrix" object: 'x' slot is not of type "double"}
#'     — counts stored as integers in the source file.
#'   \item Seurat refusing to construct an assay from a row-compressed
#'     sparse matrix.
#' }
#'
#' For \code{as.Seurat} failures in Seurat v5 when the source AnnData
#' has only a single expression matrix (so \code{counts == data}), this
#' function transparently falls back to building the Seurat object
#' directly from the SCE assays, colData, and reducedDims.
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
  # Try Seurat::as.Seurat first; it handles reductions, alt assays, etc.
  # In Seurat v5 it fails when `counts == data` because internally it
  # validates `data` against the Seurat layer namespace ("counts",
  # "data", "scale.data"), not against SCE assay names. Fall back to
  # constructing the object manually in that case.
  tryCatch(
    Seurat::as.Seurat(sce, counts = counts, data = data),
    error = function(e) {
      message("  as.Seurat failed (", conditionMessage(e), ")")
      message("  Falling back to manual Seurat construction.")

      cmat <- SummarizedExperiment::assay(sce, counts)
      cmat <- methods::as(cmat, "CsparseMatrix")
      if (!is.double(cmat@x)) cmat@x <- as.numeric(cmat@x)
      obj <- SeuratObject::CreateSeuratObject(counts = cmat)

      # If `data` is a distinct assay, write it into the data layer
      if (!identical(data, counts) &&
          data %in% SummarizedExperiment::assayNames(sce)) {
        dmat <- SummarizedExperiment::assay(sce, data)
        dmat <- methods::as(dmat, "CsparseMatrix")
        if (!is.double(dmat@x)) dmat@x <- as.numeric(dmat@x)
        SeuratObject::LayerData(obj, assay = "RNA", layer = "data") <- dmat
      }

      # Carry over colData -> @meta.data (skip cols Seurat already added)
      cd <- as.data.frame(SummarizedExperiment::colData(sce))
      cd <- cd[, !colnames(cd) %in% colnames(obj@meta.data), drop = FALSE]
      if (ncol(cd)) obj@meta.data <- cbind(obj@meta.data, cd)

      # Carry over reducedDims as DimReducs (best-effort)
      if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        for (rn in SingleCellExperiment::reducedDimNames(sce)) {
          embed <- SingleCellExperiment::reducedDim(sce, rn)
          if (is.null(colnames(embed))) {
            colnames(embed) <- paste0(toupper(rn), "_", seq_len(ncol(embed)))
          }
          obj[[tolower(rn)]] <- SeuratObject::CreateDimReducObject(
            embeddings = embed, key = paste0(toupper(rn), "_"),
            assay = SeuratObject::DefaultAssay(obj)
          )
        }
      }
      obj
    }
  )
}
