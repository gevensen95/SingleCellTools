#' Create Seurat Objects from Parse Biosciences Pipeline Output
#'
#' Reads one or more Parse Biosciences sample directories and constructs a
#' Seurat object for each, applying basic cell- and feature-level filtering.
#' When given a single path, a single \code{Seurat} object is returned; when
#' given multiple paths, a named list of \code{Seurat} objects is returned.
#'
#' @param paths Character vector of paths to Parse Biosciences sample
#'   directories. Each path must contain a \code{DGE_filtered/} subdirectory
#'   with \code{count_matrix.mtx}, \code{all_genes.csv}, and
#'   \code{cell_metadata.csv}.
#' @param sample_names Optional character vector of names to assign to the
#'   resulting Seurat objects. Must be the same length as \code{paths}. If
#'   \code{NULL} (default), \code{basename(paths)} is used.
#' @param mincellfrac Numeric in \code{[0, 1]}. Minimum fraction of cells in
#'   which a feature must be detected to be retained. Passed to
#'   \code{min.cells} in \code{\link[SeuratObject]{CreateAssayObject}} after
#'   multiplication by the number of cells. Default \code{0.0005}.
#' @param mincellfeat Integer. Minimum number of features a cell must express
#'   to be retained. Passed to \code{min.features} in
#'   \code{\link[SeuratObject]{CreateAssayObject}}. Default \code{50}.
#'
#' @return If \code{length(paths) == 1}, a single \code{Seurat} object.
#'   Otherwise, a named list of \code{Seurat} objects, one per input path,
#'   with names taken from \code{sample_names} (or \code{basename(paths)}).
#'
#' @details
#' For each input path the function:
#' \enumerate{
#'   \item Reads the sparse count matrix from
#'     \code{DGE_filtered/count_matrix.mtx} and transposes it so genes are
#'     rows and cells are columns.
#'   \item Reads gene and cell metadata from \code{all_genes.csv} and
#'     \code{cell_metadata.csv}.
#'   \item Assigns gene names to the matrix rows and cell barcodes
#'     (\code{bc_wells}) to the columns.
#'   \item Builds an RNA assay with \code{\link[SeuratObject]{CreateAssayObject}},
#'     applying the \code{mincellfeat} and \code{mincellfrac} filters.
#'   \item Wraps the assay in a \code{Seurat} object with the cell metadata
#'     attached.
#' }
#'
#' @examples
#' \dontrun{
#' # Single sample
#' obj <- MakeParseObj("/data/parse/S1")
#'
#' # Multiple samples — returns a named list
#' sample_paths <- file.path("/data/parse", paste0("S", 1:16))
#' obj_list <- MakeParseObj(sample_paths)
#'
#' # With explicit sample names
#' obj_list <- MakeParseObj(
#'   sample_paths,
#'   sample_names = paste0("Sample", 1:16)
#' )
#' }
#'
#' @importFrom Matrix readMM t
#' @importFrom methods as
#' @importFrom utils read.csv
#' @importFrom tibble column_to_rownames
#' @importFrom SeuratObject CreateAssayObject CreateSeuratObject
#' @export
MakeParseObj <- function(paths,
                         sample_names = NULL,
                         mincellfrac = 0.0005,
                         mincellfeat = 50) {

  # ---- Argument checks -----------------------------------------------------
  if (!is.character(paths) || length(paths) < 1) {
    stop("`paths` must be a non-empty character vector.")
  }
  missing_dirs <- !dir.exists(paths)
  if (any(missing_dirs)) {
    stop("The following path(s) do not exist:\n  ",
         paste(paths[missing_dirs], collapse = "\n  "))
  }

  if (is.null(sample_names)) {
    sample_names <- basename(paths)
  } else if (length(sample_names) != length(paths)) {
    stop("`sample_names` must be the same length as `paths`.")
  }

  # ---- Helper: build a single Seurat object from one path ------------------
  .make_one <- function(path) {
    dge_dir <- file.path(path, "DGE_filtered")
    if (!dir.exists(dge_dir)) {
      stop("No 'DGE_filtered' directory found at: ", path)
    }

    # 1. Read and orient the sparse count matrix (genes x cells).
    counts <- Matrix::readMM(file.path(dge_dir, "count_matrix.mtx"))
    counts <- Matrix::t(counts)
    counts <- methods::as(counts, "CsparseMatrix")

    # 2. Read gene and cell metadata.
    genes <- utils::read.csv(file.path(dge_dir, "all_genes.csv"),
                             stringsAsFactors = FALSE)
    cells <- utils::read.csv(file.path(dge_dir, "cell_metadata.csv"),
                             stringsAsFactors = FALSE)

    # 3. Assign dimnames (ensure they are unique).
    rownames(counts) <- make.unique(genes$gene_name)
    colnames(counts) <- cells$bc_wells

    # 4. Build the RNA assay and wrap in a Seurat object.
    assay <- SeuratObject::CreateAssayObject(
      counts,
      assay        = "RNA",
      min.features = mincellfeat,
      min.cells    = ceiling(nrow(cells) * mincellfrac)
    )
    SeuratObject::CreateSeuratObject(
      assay,
      meta.data = tibble::column_to_rownames(cells, "bc_wells")
    )
  }

  # ---- Build object(s) -----------------------------------------------------
  objs <- lapply(paths, .make_one)
  names(objs) <- sample_names

  if (length(objs) == 1L) {
    return(objs[[1L]])
  }
  objs
}
