#' Copy-number inference from expression with copyKAT
#'
#' Wraps \code{copykat::copykat} to infer per-cell copy-number profiles and
#' classify cells as aneuploid (likely malignant) vs. diploid (likely
#' normal) directly from expression counts -- no matched normal/tumor pair
#' or genotyping needed, just a set of cells you're confident are normal to
#' anchor the baseline.
#'
#' \code{copykat} writes several files (plots, the segmented CNA matrix, a
#' prediction table) into the current working directory with no way to
#' redirect that via an argument -- this function creates \code{output_dir},
#' \code{setwd()}s into it for the duration of the call, and restores the
#' original working directory afterward (even on error), so a package
#' function never leaves stray files in the caller's working directory.
#'
#' @param obj A Seurat object.
#' @param assay Assay to read counts from. Default \code{DefaultAssay(obj)}.
#' @param normal_cells Character vector of cell barcodes known/believed to
#'   be normal (non-malignant), used to anchor copyKAT's baseline. \code{NULL}
#'   (default) lets copyKAT estimate the baseline itself from the most
#'   stable cells -- less reliable than supplying real normal cells if you
#'   have them (e.g. immune cells in a tumor sample).
#' @param genome \code{"hg20"} (default) or \code{"mm10"}.
#' @param output_dir Directory copyKAT writes its output files to. Default
#'   a fresh \code{tempfile()} directory.
#' @param id_type Gene ID type passed to copyKAT (\code{"S"} = gene symbol,
#'   default; \code{"E"} = Ensembl).
#' @param ngene_chr Minimum genes per chromosome to keep it in the analysis.
#'   copyKAT's own default, 5.
#' @param sam_name Sample name prefix for copyKAT's output files. Default
#'   \code{"sample"}.
#' @param n_cores Cores for copyKAT's internal parallelization. Default 1.
#' @param write_metadata Logical; if \code{TRUE} (default), writes a
#'   \code{copykat_call} metadata column (\code{"aneuploid"} / \code{"diploid"}
#'   / \code{"not.defined"}).
#' @param verbose Message progress. Default \code{TRUE}.
#' @return \code{obj} with \code{obj@misc$copykat} (list with
#'   \code{prediction} and \code{cna_matrix}, copyKAT's own output) and,
#'   if \code{write_metadata}, a \code{copykat_call} metadata column.
#' @examples
#' \dontrun{
#' obj <- InferCNVWrapper(obj, normal_cells = colnames(obj)[obj$cell_type == "T cell"])
#' table(obj$copykat_call)
#' DimPlot(obj, group.by = "copykat_call")
#' }
#' @importFrom Seurat DefaultAssay GetAssayData
#' @export
InferCNVWrapper <- function(obj,
                            assay          = NULL,
                            normal_cells   = NULL,
                            genome         = c("hg20", "mm10"),
                            output_dir     = tempfile("copykat_"),
                            id_type        = "S",
                            ngene_chr      = 5,
                            sam_name       = "sample",
                            n_cores        = 1,
                            write_metadata = TRUE,
                            verbose        = TRUE) {

  genome <- match.arg(genome)
  .assert_seurat(obj)
  if (!requireNamespace("copykat", quietly = TRUE)) {
    stop("'copykat' is required. Install with ",
         "remotes::install_github('navinlabcode/copykat').")
  }

  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay
  counts <- as.matrix(Seurat::GetAssayData(obj, assay = a, layer = "counts"))

  norm_names <- if (!is.null(normal_cells)) {
    keep <- intersect(normal_cells, colnames(counts))
    if (length(keep) == 0) {
      warning("None of `normal_cells` matched obj's cells; letting copyKAT ",
              "estimate its own baseline.")
    }
    keep
  } else {
    ""
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(output_dir)

  if (isTRUE(verbose)) {
    message(sprintf(
      "--- Running copyKAT (%d cells, %d genes, genome = %s%s) ---",
      ncol(counts), nrow(counts), genome,
      if (length(norm_names) > 0) sprintf(", %d normal cell(s) supplied", length(norm_names)) else ""))
    message(sprintf("  Output files: %s", output_dir))
  }

  ck <- copykat::copykat(
    rawmat          = counts,
    id.type         = id_type,
    ngene.chr       = ngene_chr,
    sam.name        = sam_name,
    norm.cell.names = norm_names,
    genome          = genome,
    n.cores         = n_cores
  )

  obj@misc$copykat <- list(
    prediction = ck$prediction,
    cna_matrix = ck$CNAmat
  )

  if (isTRUE(write_metadata)) {
    pred <- ck$prediction
    call_by_cell <- setNames(as.character(pred$copykat.pred), pred$cell.names)
    obj$copykat_call <- call_by_cell[colnames(obj)]
    if (isTRUE(verbose)) {
      message("  Calls: ", paste(names(table(obj$copykat_call)),
                                 table(obj$copykat_call), sep = "=", collapse = ", "))
    }
  }

  obj
}
