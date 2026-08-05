#' Cell-type composition per sample / condition
#'
#' Returns counts and proportions of cells in each \code{group_col} (typically
#' cluster or cell type) per \code{sample_col}, optionally with a chi-square
#' or Fisher's exact test comparing distributions across a
#' \code{condition_col} grouping of the samples.
#'
#' @details
#' The optional \code{test} is computed on cells pooled across all samples
#' within each condition (a single \code{group_col} x \code{condition_col}
#' contingency table). That treats every \emph{cell} as an independent
#' replicate rather than every \emph{sample} -- with few samples and many
#' cells per sample (the usual case in scRNA-seq), this inflates the
#' effective sample size the test sees and can report significance that
#' isn't actually supported at the sample level (the "pseudoreplication"
#' problem well documented for single-cell composition testing). Treat it as
#' a rough, descriptive check of the pooled data, not a substitute for a
#' proper sample-level test. \code{\link{CompositionalTest}} (propeller /
#' betareg / Wilcoxon, all operating on per-sample proportions) is the
#' function to reach for when you need a test whose p-value actually
#' reflects the number of biological replicates.
#'
#' @param obj A Seurat object.
#' @param group_col Metadata column with the grouping to analyze (cluster
#'   or cell-type labels).
#' @param sample_col Metadata column identifying samples / replicates.
#' @param condition_col Optional metadata column grouping samples into
#'   conditions for testing.
#' @param test One of \code{"none"} (default), \code{"chisq"}, or
#'   \code{"fisher"}. Tests whether group proportions differ across
#'   \code{condition_col}. See Details for an important caveat about what
#'   this test's p-value does and doesn't account for -- emits a
#'   \code{warning()} to that effect whenever it's requested.
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
    warning(
      "This ", test, " test is computed on cells pooled across all samples ",
      "within each condition -- it treats every cell as an independent ",
      "replicate, not just every sample. With few samples and many cells ",
      "per sample (the common case), this inflates the effective sample ",
      "size and can report significance that isn't actually supported at ",
      "the sample level (pseudoreplication). Use CompositionalTest() ",
      "instead for a test that operates on per-sample proportions.",
      call. = FALSE
    )
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
