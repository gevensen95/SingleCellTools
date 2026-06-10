#' Combine per-cell metadata from a list of Seurat objects
#'
#' Stacks the \code{@meta.data} of every Seurat object in a list into a single
#' long-format data frame, with the original list name carried in a
#' \code{sample} column and the original cell barcode in a \code{cell_id}
#' column. Tolerates objects with different metadata columns
#' (missing ones become \code{NA}) via \code{dplyr::bind_rows}.
#'
#' @param seurat_list A list of Seurat objects.
#' @param sample_col Name of the column in the result that records which list
#'   element each row came from. Default \code{"sample"}.
#' @param cell_col Name of the column that holds the original cell barcode.
#'   Default \code{"cell_id"}.
#' @return A tibble with all per-cell metadata rows from every input object.
#' @importFrom dplyr bind_rows
#' @importFrom tibble rownames_to_column
#' @export
combine_metadata <- function(seurat_list,
                             sample_col = "sample",
                             cell_col   = "cell_id") {
  if (inherits(seurat_list, "Seurat")) seurat_list <- list(seurat_list)
  if (is.null(names(seurat_list))) {
    names(seurat_list) <- paste0("obj_", seq_along(seurat_list))
  }

  per_obj <- lapply(seurat_list, function(o) {
    tibble::rownames_to_column(o@meta.data, var = cell_col)
  })

  dplyr::bind_rows(per_obj, .id = sample_col)
}
