#' Cell-type composition per sample / condition
#'
#' Returns counts and proportions of cells in each \code{group_col} (typically
#' cluster or cell type) per \code{sample_col}, optionally with a chi-square
#' or Fisher's exact test comparing distributions across a
#' \code{condition_col} grouping of the samples.
#'
#' @param obj A Seurat object.
#' @param group_col Metadata column with the grouping to analyze (cluster
#'   or cell-type labels).
#' @param sample_col Metadata column identifying samples / replicates.
#' @param condition_col Optional metadata column grouping samples into
#'   conditions for testing.
#' @param test One of \code{"none"} (default), \code{"chisq"}, or
#'   \code{"fisher"}. Tests whether group proportions differ across
#'   \code{condition_col}.
#' @return A list with elements:
#'   \describe{
#'     \item{\code{counts}}{Long-format tibble of cell counts.}
#'     \item{\code{proportions}}{Long-format tibble of within-sample fractions.}
#'     \item{\code{test}}{Test result (if requested) or \code{NULL}.}
#'   }
#' @importFrom dplyr group_by mutate ungroup count
#' @importFrom tibble as_tibble
#' @importFrom stats chisq.test fisher.test
#' @export
CompositionAnalysis <- function(obj,
                                group_col,
                                sample_col,
                                condition_col = NULL,
                                test          = c("none", "chisq", "fisher")) {
  test <- match.arg(test)
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  md <- obj@meta.data
  for (col in c(group_col, sample_col)) {
    if (!col %in% colnames(md)) {
      stop("Column '", col, "' not found in obj@meta.data.")
    }
  }

  group  <- factor(md[[group_col]])
  sample <- factor(md[[sample_col]])

  counts <- as.data.frame(table(sample = sample, group = group),
                          stringsAsFactors = FALSE)
  colnames(counts)[3] <- "n_cells"

  # Within-sample proportions
  totals <- tapply(counts$n_cells, counts$sample, sum)
  counts$prop <- counts$n_cells / totals[as.character(counts$sample)]

  # Attach condition mapping if provided
  if (!is.null(condition_col)) {
    if (!condition_col %in% colnames(md)) {
      stop("Column '", condition_col, "' not found in obj@meta.data.")
    }
    sample_to_cond <- unique(md[, c(sample_col, condition_col)])
    colnames(sample_to_cond) <- c("sample", "condition")
    counts <- merge(counts, sample_to_cond, by = "sample",
                    all.x = TRUE, sort = FALSE)
  }

  # Tidy
  counts_tib <- tibble::as_tibble(counts)
  props_tib  <- counts_tib

  # Optional test: build contingency of group x condition (collapse samples)
  test_result <- NULL
  if (test != "none") {
    if (is.null(condition_col)) {
      stop("`condition_col` must be supplied to run a test.")
    }
    cond_lookup <- unique(md[, c(sample_col, condition_col)])
    cond <- cond_lookup[match(as.character(sample), cond_lookup[[sample_col]]),
                        condition_col]
    tab <- table(group = group, condition = cond)
    test_result <- switch(
      test,
      chisq  = stats::chisq.test(tab),
      fisher = stats::fisher.test(tab, simulate.p.value = TRUE, B = 10000)
    )
  }

  list(
    counts      = counts_tib,
    proportions = props_tib,
    test        = test_result
  )
}
