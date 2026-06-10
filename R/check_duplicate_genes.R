#' Check for duplicate gene names across a list of Seurat objects
#'
#' Reports, per object, how many features in the chosen assay have duplicated
#' rownames and which they are. Catches the most common cause of the
#' \code{"duplicate 'row.names' are not allowed"} error during \code{merge()}.
#'
#' @param seurat_list A Seurat object or a (optionally named) list of them.
#' @param assay Assay whose rownames to inspect. Defaults to each object's
#'   DefaultAssay.
#' @param max_show How many duplicated genes to print per object. Default 20.
#' @return Invisibly, a list with one character vector per input object
#'   containing the duplicated gene names (empty if none).
#' @importFrom Seurat DefaultAssay
#' @importFrom utils head
#' @export
check_duplicate_genes <- function(seurat_list, assay = NULL, max_show = 20) {
  # Accept single object too
  if (inherits(seurat_list, "Seurat")) seurat_list <- list(seurat_list)

  nms <- names(seurat_list)
  if (is.null(nms)) nms <- paste0("obj_", seq_along(seurat_list))

  cat(sprintf("%-20s %10s %10s %10s\n",
              "object", "n_genes", "n_unique", "n_dup"))
  cat(strrep("-", 53), "\n", sep = "")

  out <- vector("list", length(seurat_list))
  names(out) <- nms

  for (i in seq_along(seurat_list)) {
    o  <- seurat_list[[i]]
    a  <- if (is.null(assay)) Seurat::DefaultAssay(o) else assay
    rn <- rownames(o[[a]])
    dup_mask  <- duplicated(rn) | duplicated(rn, fromLast = TRUE)
    dup_genes <- sort(unique(rn[dup_mask]))

    cat(sprintf("%-20s %10d %10d %10d\n",
                nms[i], length(rn), length(unique(rn)), length(dup_genes)))

    if (length(dup_genes)) {
      shown <- utils::head(dup_genes, max_show)
      cat("  ", paste(shown, collapse = ", "), sep = "")
      if (length(dup_genes) > max_show) {
        cat(sprintf("  ...and %d more", length(dup_genes) - max_show))
      }
      cat("\n")
    }
    out[[i]] <- dup_genes
  }

  invisible(out)
}
