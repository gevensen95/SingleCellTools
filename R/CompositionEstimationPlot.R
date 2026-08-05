#' Bootstrap effect-size ("estimation") plots for composition shifts
#'
#' Takes the \code{proportions} tibble from \code{\link{CompositionAnalysis}}
#' (must have been run with \code{condition_col} set) and, for one or more
#' \code{group_col} levels (cell types / clusters), builds a
#' \code{dabestr} estimation plot comparing per-sample proportions between
#' two conditions: a bootstrap 95\% confidence interval on the effect size,
#' shown alongside the raw per-sample values, rather than (or in addition
#' to) a bare p-value from \code{CompositionAnalysis(test = ...)}.
#'
#' @details
#' \code{CompositionAnalysis()}'s chi-square/Fisher tests answer "do these
#' distributions differ." This function instead asks "by how much, with
#' what uncertainty" for one cell type/cluster at a time -- the two are
#' complementary, not redundant. See
#' \url{https://acclab.github.io/dabestr/} for the estimation-statistics
#' framework this wraps.
#'
#' @param comp The list returned by \code{\link{CompositionAnalysis}},
#'   called with \code{condition_col} set. Only \code{comp$proportions} is
#'   used; \code{comp$test} (if present) is ignored here.
#' @param group_levels Character vector of \code{group_col} levels (e.g.
#'   cell types) to build a plot for. \code{NULL} (default) uses every
#'   level present in \code{comp$proportions}.
#' @param idx Length-2 character vector giving the two \code{condition_col}
#'   levels to compare, reference first (e.g.
#'   \code{c("Vehicle", "DrugA")}). \code{NULL} (default) auto-detects the
#'   2 levels present, sorted alphabetically, with a message -- errors if
#'   more than 2 conditions are present and \code{idx} wasn't supplied
#'   explicitly, since silently picking 2 out of 3+ would be a
#'   consequential, easy-to-miss decision.
#' @param effect One of \code{"mean_diff"} (default), \code{"median_diff"},
#'   \code{"cohens_d"}, \code{"hedges_g"}, \code{"cliffs_delta"}, or
#'   \code{"cohens_h"}. \code{"cohens_h"} is the effect size specifically
#'   designed for comparing two proportions and may be a more principled
#'   choice than the default \code{"mean_diff"} for this proportion-valued
#'   data -- see the dabestr documentation for details on each.
#' @return A single \code{dabestr} estimation plot if \code{group_levels}
#'   resolves to exactly one level, otherwise a named list of plots (one
#'   per \code{group_levels} entry, names = the group level).
#' @examples
#' \dontrun{
#' comp <- CompositionAnalysis(obj, group_col = "cell_type",
#'                             sample_col = "orig.ident",
#'                             condition_col = "treatment")
#'
#' # One cell type
#' CompositionEstimationPlot(comp, group_levels = "T cell",
#'                           idx = c("Vehicle", "DrugA"))
#'
#' # Every cell type present, as a named list
#' plots <- CompositionEstimationPlot(comp, idx = c("Vehicle", "DrugA"))
#' plots[["T cell"]]
#'
#' # Cohen's h -- the proportion-specific effect size
#' CompositionEstimationPlot(comp, idx = c("Vehicle", "DrugA"),
#'                           effect = "cohens_h")
#' }
#' @export
CompositionEstimationPlot <- function(comp,
                                      group_levels = NULL,
                                      idx          = NULL,
                                      effect       = c("mean_diff", "median_diff",
                                                       "cohens_d", "hedges_g",
                                                       "cliffs_delta", "cohens_h")) {
  effect <- match.arg(effect)

  if (!is.list(comp) || !"proportions" %in% names(comp)) {
    stop("`comp` must be the list returned by CompositionAnalysis().")
  }
  df <- comp$proportions
  if (!"condition" %in% colnames(df)) {
    stop("`comp$proportions` has no `condition` column -- re-run ",
         "CompositionAnalysis() with `condition_col` set.")
  }

  dabest_objs <- .dabest_from_long(df,
                                   group_col     = "group",
                                   y_col         = "prop",
                                   condition_col = "condition",
                                   idx           = idx,
                                   effect        = effect,
                                   group_levels  = group_levels)

  plots <- lapply(dabest_objs, dabestr::dabest_plot)
  if (length(plots) == 1) plots[[1]] else plots
}
