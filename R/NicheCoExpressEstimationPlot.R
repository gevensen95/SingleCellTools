#' Bootstrap effect-size ("estimation") plots for niche co-expression shifts
#'
#' Takes the list returned by \code{\link{NicheCoExpress}} and, for one or
#' more (niche, gene-pair) combinations, builds a \code{dabestr} estimation
#' plot comparing per-sample co-expression scores between the two
#' conditions: a bootstrap 95\% confidence interval on the effect size,
#' shown alongside the raw per-sample values, rather than (or in addition
#' to) the Wilcoxon/t-test p-value in \code{res$stats}.
#'
#' @details
#' \code{res$stats$delta} is always \code{mean(test) - mean(reference)},
#' using whichever order \code{conditions} resolved to inside
#' \code{NicheCoExpress()} (stored as \code{attr(res$stats, "conditions")}).
#' By default this function reuses that exact order for its own effect-size
#' calculation, rather than re-deriving its own (e.g. alphabetical) default
#' -- if it didn't, and \code{NicheCoExpress()} had been called with a
#' non-alphabetical \code{conditions} argument, the two functions could
#' silently disagree on which condition is the reference, flipping the
#' sign of the effect size relative to \code{stats$delta}. Pass \code{idx}
#' explicitly to override.
#'
#' @param res The list returned by \code{\link{NicheCoExpress}}.
#' @param niches Optional character vector restricting to these niches.
#'   \code{NULL} (default) uses every niche present in \code{res$per_sample}.
#' @param pairs Optional character vector of pair IDs (\code{"geneA_geneB"},
#'   matching \code{res$per_sample$pair}) restricting to these gene pairs.
#'   \code{NULL} (default) uses every pair present.
#' @param idx Length-2 character vector giving the two conditions to
#'   compare, reference first. \code{NULL} (default) uses
#'   \code{attr(res$stats, "conditions")} -- see Details.
#' @param effect One of \code{"mean_diff"} (default), \code{"median_diff"},
#'   \code{"cohens_d"}, \code{"hedges_g"}, \code{"cliffs_delta"}, or
#'   \code{"cohens_h"}. See \code{\link{CompositionEstimationPlot}} (which
#'   shares this same set of options) for notes on each.
#' @return A single \code{dabestr} estimation plot if the (niche, pair)
#'   filter resolves to exactly one combination, otherwise a named list of
#'   plots (names formatted as \code{"<niche> | <pair>"}).
#' @examples
#' \dontrun{
#' res <- NicheCoExpress(obj, genes = c("Vegfa", "Kdr"), niche_col = "niche",
#'                       sample_col = "sample", condition_col = "condition")
#'
#' # One niche x pair combination
#' NicheCoExpressEstimationPlot(res, niches = "N1", pairs = "Vegfa_Kdr")
#'
#' # Every niche x pair combination present, as a named list
#' plots <- NicheCoExpressEstimationPlot(res)
#' plots[["N1 | Vegfa_Kdr"]]
#' }
#' @export
NicheCoExpressEstimationPlot <- function(res,
                                         niches = NULL,
                                         pairs  = NULL,
                                         idx    = NULL,
                                         effect = c("mean_diff", "median_diff",
                                                   "cohens_d", "hedges_g",
                                                   "cliffs_delta", "cohens_h")) {
  effect <- match.arg(effect)

  if (!is.list(res) || !all(c("per_sample", "stats") %in% names(res))) {
    stop("`res` must be the list returned by NicheCoExpress().")
  }
  df <- res$per_sample
  if (!"condition" %in% colnames(df)) {
    stop("`res$per_sample` has no `condition` column -- was this really ",
         "produced by NicheCoExpress()?")
  }

  if (is.null(idx)) {
    idx <- attr(res$stats, "conditions")  # see Details -- keep in sync with
                                          # stats$delta's reference/test order
  }

  if (!is.null(niches)) df <- df[df$niche %in% niches, , drop = FALSE]
  if (!is.null(pairs))  df <- df[df$pair  %in% pairs,  , drop = FALSE]
  df <- df[!is.na(df$coexpr), , drop = FALSE]
  if (nrow(df) == 0) {
    stop("No per-sample scores left after filtering by `niches`/`pairs`. ",
         "Check that the requested niche(s)/pair(s) are present in ",
         "`res$per_sample`.")
  }

  df$niche_pair <- paste(df$niche, df$pair, sep = " | ")

  dabest_objs <- .dabest_from_long(df,
                                   group_col     = "niche_pair",
                                   y_col         = "coexpr",
                                   condition_col = "condition",
                                   idx           = idx,
                                   effect        = effect)

  plots <- lapply(dabest_objs, dabestr::dabest_plot)
  if (length(plots) == 1) plots[[1]] else plots
}
