#' Bootstrap effect-size ("estimation") plots for gene-positivity shifts
#'
#' Takes the list returned by \code{\link{GenePositivityAnalysis}} and, for
#' one or more genes (and, if computed with \code{group_col}, gene x group
#' combinations), builds a \code{dabestr} estimation plot comparing
#' per-sample positivity rates between two conditions: a bootstrap 95\%
#' confidence interval on the effect size, shown alongside the raw per-sample
#' values, rather than (or in addition to) the chi-square/Fisher p-value in
#' \code{gpa$test}.
#'
#' @details
#' \code{GenePositivityAnalysis()}'s chi-square/Fisher test answers "is there
#' a difference," operating on pooled cell-level counts. This function
#' instead asks "by how much, with what uncertainty," using per-sample
#' positivity rates as the unit of observation -- the two are complementary,
#' not redundant, in exactly the same way \code{\link{CompositionAnalysis}}
#' and \code{\link{CompositionEstimationPlot}} are.
#'
#' \code{"cohens_h"} is the effect size specifically designed for comparing
#' two proportions and may be a more principled choice than the default
#' \code{"mean_diff"} for this proportion-valued data -- see the dabestr
#' documentation for details on each.
#'
#' @param gpa The list returned by \code{\link{GenePositivityAnalysis}},
#'   called with \code{condition_col} set. Only \code{gpa$proportions} is
#'   used; \code{gpa$test} (if present) is ignored here.
#' @param genes Optional character vector restricting to these genes.
#'   \code{NULL} (default) uses every gene present in \code{gpa$proportions}.
#' @param group_levels Optional character vector restricting to these
#'   \code{group_col} levels. Only valid if \code{gpa} was computed with
#'   \code{group_col} set -- errors if supplied otherwise.
#'   \code{NULL} (default) uses every level present.
#' @param idx Length-2 character vector giving the two conditions to
#'   compare, reference first. \code{NULL} (default) auto-detects the 2
#'   levels present, sorted alphabetically, with a message -- errors if more
#'   than 2 conditions are present and \code{idx} wasn't supplied explicitly.
#' @param effect One of \code{"mean_diff"} (default), \code{"median_diff"},
#'   \code{"cohens_d"}, \code{"hedges_g"}, \code{"cliffs_delta"}, or
#'   \code{"cohens_h"}. See Details.
#' @return A single \code{dabestr} estimation plot if the gene (and group,
#'   if applicable) filter resolves to exactly one combination, otherwise a
#'   named list of plots (names = the gene, or \code{"<gene> | <group>"} if
#'   \code{group_col} was used). A combination present in only one of the
#'   two \code{idx} conditions is skipped with a \code{warning()} naming it
#'   and the missing condition, rather than erroring the whole call -- this
#'   is expected with imbalanced data and errors only if \emph{no} requested
#'   combination has both conditions represented.
#' @examples
#' \dontrun{
#' obj <- AddGenePositivity(obj, genes = c("CD3D", "CD4"))
#' res <- GenePositivityAnalysis(obj, genes = c("CD3D", "CD4"),
#'                               sample_col = "orig.ident",
#'                               condition_col = "treatment")
#'
#' # One gene
#' GenePositivityEstimationPlot(res, genes = "CD3D", idx = c("Vehicle", "DrugA"))
#'
#' # Every gene present, as a named list
#' plots <- GenePositivityEstimationPlot(res, idx = c("Vehicle", "DrugA"))
#' plots[["CD3D"]]
#'
#' # Cohen's h -- the proportion-specific effect size
#' GenePositivityEstimationPlot(res, idx = c("Vehicle", "DrugA"),
#'                              effect = "cohens_h")
#' }
#' @export
GenePositivityEstimationPlot <- function(gpa,
                                         genes        = NULL,
                                         group_levels = NULL,
                                         idx          = NULL,
                                         effect       = c("mean_diff", "median_diff",
                                                          "cohens_d", "hedges_g",
                                                          "cliffs_delta", "cohens_h")) {
  effect <- match.arg(effect)

  if (!is.list(gpa) || !"proportions" %in% names(gpa)) {
    stop("`gpa` must be the list returned by GenePositivityAnalysis().")
  }
  df <- gpa$proportions
  if (!"condition" %in% colnames(df)) {
    stop("`gpa$proportions` has no `condition` column -- re-run ",
         "GenePositivityAnalysis() with `condition_col` set.")
  }
  if (!"gene" %in% colnames(df)) {
    stop("`gpa$proportions` has no `gene` column -- was this really produced ",
         "by GenePositivityAnalysis()?")
  }

  has_group <- "group" %in% colnames(df)
  if (!has_group && !is.null(group_levels)) {
    stop("`group_levels` was supplied but `gpa` wasn't computed with ",
         "`group_col` set -- rerun GenePositivityAnalysis() with `group_col` ",
         "if you want to filter/stratify by it.")
  }

  if (!is.null(genes)) df <- df[df$gene %in% genes, , drop = FALSE]
  if (has_group && !is.null(group_levels)) {
    df <- df[df$group %in% group_levels, , drop = FALSE]
  }
  if (nrow(df) == 0) {
    stop("No rows left after filtering by `genes`",
         if (has_group) " / `group_levels`", ". Check that the requested ",
         "gene(s) (and group level(s)) are present in `gpa$proportions`.")
  }

  df$feature <- if (has_group) paste(df$gene, df$group, sep = " | ") else df$gene

  dabest_objs <- .dabest_from_long(df,
                                   group_col     = "feature",
                                   y_col         = "prop_pos",
                                   condition_col = "condition",
                                   idx           = idx,
                                   effect        = effect)

  plots <- lapply(dabest_objs, dabestr::dabest_plot)
  if (length(plots) == 1) plots[[1]] else plots
}
