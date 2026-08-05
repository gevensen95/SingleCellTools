#' Gene-positivity rates per sample / condition
#'
#' Companion to \code{\link{CompositionAnalysis}}, but for per-gene positivity
#' (from \code{\link{AddGenePositivity}}) rather than cluster/cell-type
#' composition. Returns, for each requested gene, the fraction of cells
#' positive for that gene per \code{sample_col} (optionally further split by
#' \code{group_col}, e.g. cell type / cluster / niche), and optionally a
#' chi-square or Fisher's exact test comparing positivity rates across a
#' \code{condition_col} grouping of the samples.
#'
#' @details
#' This function does not compute positivity itself -- it reads the
#' \code{<gene><suffix>} logical metadata columns that
#' \code{\link{AddGenePositivity}} already writes, the same way
#' \code{\link{PlotGenePositivity}} does. Run \code{AddGenePositivity()} first;
#' this errors out (naming the missing columns) if it wasn't.
#'
#' When \code{test} is requested, one test is run per gene (and, if
#' \code{group_col} is supplied, per gene x group combination) on a 2x2
#' positive/negative x condition table built directly from cell-level data --
#' the same level of aggregation \code{\link{CompositionAnalysis}}'s own test
#' uses, and the same caveat applies: cells, not samples, are the test's
#' unit of replication, so with few samples and many cells per sample the
#' p-value can report significance that isn't actually supported at the
#' sample level (pseudoreplication) -- emits a \code{warning()} to that
#' effect whenever \code{test != "none"}. Unlike \code{CompositionAnalysis}
#' there's no ready-made sample-level significance test for gene positivity
#' in this package; \code{\link{GenePositivityEstimationPlot}}'s bootstrap CI
#' on the per-\emph{sample} rates is the sample-respecting view to lean on,
#' or run your own paired test (e.g. \code{wilcox.test()}) on
#' \code{gpa$proportions$prop_pos} per gene if you need a p-value that
#' reflects the number of biological replicates.
#'
#' @param obj A Seurat object with \code{AddGenePositivity()} already run.
#' @param genes Character vector of gene symbols to analyze. Each must have a
#'   corresponding \code{paste0(gene, suffix)} metadata column.
#' @param sample_col Metadata column identifying samples / replicates.
#' @param condition_col Optional metadata column grouping samples into
#'   conditions for testing / downstream estimation plots.
#' @param group_col Optional metadata column to additionally stratify by
#'   (e.g. cell type, cluster, or niche). \code{NULL} (default) computes one
#'   positivity rate per gene per sample only.
#' @param suffix Suffix used by \code{AddGenePositivity} when creating the
#'   positivity columns. Default \code{"_pos"}.
#' @param test One of \code{"none"} (default), \code{"chisq"}, or
#'   \code{"fisher"}. Tests whether each gene's positivity rate differs
#'   across \code{condition_col}. See Details for an important caveat about
#'   what this test's p-value does and doesn't account for -- emits a
#'   \code{warning()} to that effect whenever it's requested.
#' @return A list with elements:
#'   \describe{
#'     \item{\code{counts}}{Long-format tibble: \code{sample}, \code{gene},
#'       optionally \code{group}, \code{n_pos}, \code{n_total},
#'       \code{prop_pos}, and \code{condition} if \code{condition_col} was
#'       supplied.}
#'     \item{\code{proportions}}{Identical to \code{counts} -- kept as a
#'       separate name for symmetry with \code{\link{CompositionAnalysis}}
#'       and so \code{\link{GenePositivityEstimationPlot}} has a stable
#'       \code{$proportions} entry to read.}
#'     \item{\code{test}}{Named list of \code{htest} results (if requested),
#'       one per gene (name = gene) or per gene x group combination (name =
#'       \code{"<gene> | <group>"} if \code{group_col} was supplied), or
#'       \code{NULL} if \code{test = "none"}.}
#'   }
#' @examples
#' \dontrun{
#' obj <- AddGenePositivity(obj, genes = c("CD3D", "CD4"))
#' res <- GenePositivityAnalysis(obj, genes = c("CD3D", "CD4"),
#'                               sample_col = "orig.ident",
#'                               condition_col = "treatment",
#'                               group_col = "cell_type",
#'                               test = "chisq")
#' res$proportions
#' res$test[["CD3D | T cell"]]
#' }
#' @importFrom tibble as_tibble
#' @importFrom stats chisq.test fisher.test
#' @export
GenePositivityAnalysis <- function(obj,
                                   genes,
                                   sample_col,
                                   condition_col = NULL,
                                   group_col     = NULL,
                                   suffix        = "_pos",
                                   test          = c("none", "chisq", "fisher")) {
  test <- match.arg(test)
  stopifnot(is.character(genes), length(genes) >= 1)
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  md <- obj@meta.data

  for (col in c(sample_col, group_col)) {
    if (!col %in% colnames(md)) {
      stop("Column '", col, "' not found in obj@meta.data.")
    }
  }

  pos_cols <- paste0(genes, suffix)
  missing_cols <- setdiff(pos_cols, colnames(md))
  if (length(missing_cols)) {
    stop("Missing positivity column(s): ", paste(missing_cols, collapse = ", "),
         ". Run AddGenePositivity(obj, genes = ...) first",
         if (suffix != "_pos") paste0(" (with suffix = '", suffix, "')"),
         ".")
  }

  sample <- factor(md[[sample_col]])
  group  <- if (!is.null(group_col)) factor(md[[group_col]]) else factor(rep("all", nrow(md)))

  counts_list <- lapply(genes, function(gene) {
    pos <- as.logical(md[[paste0(gene, suffix)]])

    total_tab <- as.data.frame(table(sample = sample, group = group),
                               stringsAsFactors = FALSE)
    colnames(total_tab)[3] <- "n_total"

    pos_tab <- as.data.frame(table(sample = sample[pos], group = group[pos]),
                             stringsAsFactors = FALSE)
    colnames(pos_tab)[3] <- "n_pos"

    merged <- merge(total_tab, pos_tab, by = c("sample", "group"),
                    all.x = TRUE, sort = FALSE)
    merged$n_pos[is.na(merged$n_pos)] <- 0
    merged$gene <- gene
    merged
  })
  counts <- do.call(rbind, counts_list)
  counts$prop_pos <- counts$n_pos / pmax(1, counts$n_total)

  if (is.null(group_col)) counts$group <- NULL

  if (!is.null(condition_col)) {
    if (!condition_col %in% colnames(md)) {
      stop("Column '", condition_col, "' not found in obj@meta.data.")
    }
    sample_to_cond <- unique(md[, c(sample_col, condition_col)])
    colnames(sample_to_cond) <- c("sample", "condition")
    counts <- merge(counts, sample_to_cond, by = "sample", all.x = TRUE, sort = FALSE)
  }

  counts_tib <- tibble::as_tibble(counts)
  props_tib  <- counts_tib

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
      "the sample level (pseudoreplication). There's no ready-made ",
      "sample-level significance test for gene positivity in this package ",
      "-- use GenePositivityEstimationPlot() for a sample-respecting ",
      "effect-size view, or run your own test (e.g. wilcox.test()) on the ",
      "per-sample `prop_pos` values if you need a p-value that reflects ",
      "the number of biological replicates.",
      call. = FALSE
    )
    cond_lookup <- unique(md[, c(sample_col, condition_col)])
    cond <- cond_lookup[match(as.character(sample), cond_lookup[[sample_col]]),
                        condition_col]

    group_levels <- if (is.null(group_col)) "all" else levels(group)
    test_result <- list()
    for (gene in genes) {
      pos <- as.logical(md[[paste0(gene, suffix)]])
      for (g in group_levels) {
        idx_cells <- if (is.null(group_col)) rep(TRUE, length(pos)) else as.character(group) == g
        tab <- table(positive = pos[idx_cells], condition = cond[idx_cells])
        key <- if (is.null(group_col)) gene else paste(gene, g, sep = " | ")
        test_result[[key]] <- switch(
          test,
          chisq  = stats::chisq.test(tab),
          fisher = stats::fisher.test(tab, simulate.p.value = TRUE, B = 10000)
        )
      }
    }
  }

  list(
    counts      = counts_tib,
    proportions = props_tib,
    test        = test_result
  )
}
