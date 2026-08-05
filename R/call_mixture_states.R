#' Call data-driven states from a Gaussian mixture model
#'
#' Fits a 1D or multivariate Gaussian mixture model (via \code{mclust}) on
#' one or more numeric metadata columns and returns a ranked, relabeled
#' state call per row -- state 1 is always the lowest on the score axis,
#' state G is always the highest. Useful anywhere you'd otherwise hand-pick
#' a quantile cutoff on a module score or composite score: BIC selects the
#' number of states for you, and the returned posterior probabilities give
#' a principled confidence/severity score instead of a hard threshold.
#'
#' Works on a Seurat object (using \code{meta.data}) or a plain data.frame/
#' tibble, and an optional \code{group_col}/\code{group_value} lets you
#' restrict the fit to a subset (e.g. one cell type) -- or omit both to fit
#' on every row, for data that isn't stratified by cell type at all.
#'
#' @param data A Seurat object or a data.frame/tibble.
#' @param score_col One or more numeric metadata column names to fit the
#'   mixture on. Length 1 fits a standard 1D mixture; length > 1 fits a
#'   joint multivariate mixture across those columns.
#' @param group_col Optional metadata column to filter/stratify by (e.g.
#'   \code{'annotation_first_pass'}). \code{NULL} (default) fits on every
#'   row of \code{data}.
#' @param group_value Value(s) of \code{group_col} to keep. Required if
#'   \code{group_col} is set; ignored otherwise.
#' @param G Range of component counts to try; BIC picks the winner. Passed
#'   directly to \code{mclust::Mclust}.
#' @param modelNames Optional \code{mclust} covariance model spec (e.g.
#'   \code{"E"}/\code{"V"} for 1D, \code{"EII"}/\code{"VVV"} for
#'   multivariate). \code{NULL} lets Mclust search all applicable models.
#' @param decreasing If \code{FALSE} (default), state 1 = lowest mean on the
#'   score axis, state G = highest. Set \code{TRUE} for scores where LOWER
#'   means more severe (e.g. a vitality score), so state G still means
#'   "most severe".
#' @param min_n Minimum number of eligible rows required to attempt a fit.
#'   Below this, returns \code{NULL} with a message instead of erroring out
#'   on a rare group. Default 10.
#' @param label Prefix for the output column names, so the function reads
#'   cleanly for whatever it's being used to call (e.g.
#'   \code{label = 'senescence'} produces \code{senescence_state}, etc.).
#'   Default \code{'state'}.
#' @param verbose Print a one-line summary (n, BIC-selected G, model) per
#'   call. Default \code{TRUE}.
#' @param ... Passed through to \code{mclust::Mclust} (e.g. \code{prior = }).
#'
#' @return \code{NULL} if skipped (too few rows), otherwise a data.frame
#'   with columns \code{id}, \code{<label>_state}, \code{<label>_confidence}
#'   (posterior of the ASSIGNED component -- a confidence score, not a
#'   severity score), and \code{<label>_prob_high} (posterior probability of
#'   NOT being in the lowest-ranked state -- the continuous score to use
#'   for ROC/AUC-style evaluation). \code{G}, \code{bic}, and
#'   \code{modelName} are attached as attributes.
#'
#' @examples
#' \dontrun{
#' # Single cell type, one score column (the original use case)
#' call_mixture_states(mouse, score_col = 'stress_composite',
#'                     group_col = 'annotation_first_pass',
#'                     group_value = 'Hepatocytes')
#'
#' # No stratification, multivariate mixture across several module scores
#' call_mixture_states(mouse, score_col = c('sasp1', 'cdki1', 'sdtmc1'),
#'                     label = 'senescence')
#' }
#'
#' @importFrom stats complete.cases
#' @importFrom mclust mclustBIC
#' @export
call_mixture_states <- function(data,
                                score_col,
                                group_col   = NULL,
                                group_value = NULL,
                                G           = 2:5,
                                modelNames  = NULL,
                                decreasing  = FALSE,
                                min_n       = 10,
                                label       = 'state',
                                verbose     = TRUE,
                                ...) {

  # ---- Accept a Seurat object or a plain data.frame -----------------------
  if (inherits(data, 'Seurat')) {
    meta <- data@meta.data
    ids  <- rownames(meta)
  } else if (is.data.frame(data)) {
    meta <- data
    ids  <- if (!is.null(rownames(data))) rownames(data) else as.character(seq_len(nrow(data)))
  } else {
    stop('`data` must be a Seurat object or a data.frame.')
  }

  if (!all(score_col %in% colnames(meta))) {
    stop('`score_col` not found in metadata: ',
        paste(setdiff(score_col, colnames(meta)), collapse = ', '))
  }

  # ---- Optional grouping/filtering ----------------------------------------
  keep <- rep(TRUE, nrow(meta))
  group_label <- 'all rows'
  if (!is.null(group_col)) {
    if (!group_col %in% colnames(meta)) {
      stop(sprintf("`group_col` '%s' not found in metadata.", group_col))
    }
    if (is.null(group_value)) {
      stop('`group_value` must be supplied when `group_col` is set.')
    }
    keep <- meta[[group_col]] %in% group_value
    group_label <- paste(group_value, collapse = ', ')
  }

  score_vals <- meta[, score_col, drop = FALSE]
  keep <- keep & stats::complete.cases(score_vals)
  idx  <- which(keep)

  if (length(idx) < min_n) {
    if (verbose) message(sprintf('  Skipping %s (n=%d < min_n=%d)',
                                 group_label, length(idx), min_n))
    return(NULL)
  }

  # 1D: plain numeric vector for Mclust. Multivariate: keep as a matrix so
  # Mclust fits a joint mixture across all requested score columns.
  # `[[1]][idx]` (rather than `[idx, 1]`) so this works whether `meta` is a
  # base data.frame or a tibble -- tibble's `[` never drops to a vector.
  x <- if (length(score_col) == 1) {
    score_vals[[1]][idx]
  } else {
    as.matrix(score_vals[idx, , drop = FALSE])
  }
  # mclust::Mclust() builds a call to mclustBIC() and evaluates it in ITS
  # caller's frame (i.e. this function's environment) rather than mclust's
  # own namespace -- so `mclustBIC` must be resolvable from here, which is
  # exactly what the @importFrom mclust mclustBIC above guarantees. Without
  # it (and without mclust attached via library()), this errors with
  # "could not find function \"mclustBIC\"" even though `mclust::Mclust`
  # itself resolves fine.
  fit <- mclust::Mclust(x, G = G, modelNames = modelNames, ...)

  if (verbose) {
    message(sprintf('  %s: n=%d | BIC selected G=%d (model=%s)',
                    group_label, length(idx), fit$G, fit$modelName))
  }

  # Rank components along the score axis. For multivariate fits,
  # fit$parameters$mean is a (variables x components) matrix -- average
  # across variables to get one ranking value per component. This assumes
  # the score columns are on comparable scales; scale them first if not.
  comp_means <- fit$parameters$mean
  if (is.matrix(comp_means)) comp_means <- colMeans(comp_means)
  ord     <- order(comp_means, decreasing = decreasing)
  relabel <- setNames(seq_along(ord), ord)
  states  <- relabel[as.character(fit$classification)]

  lowest_comp <- ord[1]
  prob_high   <- 1 - fit$z[, lowest_comp]

  out <- data.frame(
    id         = ids[idx],
    state      = factor(states, levels = sort(unique(states))),
    confidence = apply(fit$z, 1, max),
    prob_high  = prob_high,
    stringsAsFactors = FALSE
  )
  colnames(out) <- c('id', paste0(label, c('_state', '_confidence', '_prob_high')))

  structure(out, G = fit$G, bic = fit$bic, modelName = fit$modelName)
}

#' Call stress states for one or more cell types (stress-specific wrapper)
#'
#' Thin wrapper around \code{\link{call_mixture_states}} that reproduces the
#' original stress-calling signature/output shape used across the
#' aging/senescence scripts, for cases where you don't need the more general
#' interface directly.
#'
#' @param obj A Seurat object.
#' @param cell_type Value(s) of \code{annotation_first_pass} to restrict to.
#' @param score_col Metadata column(s) to fit the mixture on. Default
#'   \code{'stress_composite'}.
#' @param G Range of component counts to try. Default \code{2:5}.
#' @param ... Passed through to \code{\link{call_mixture_states}}.
#'
#' @return \code{NULL} if skipped, otherwise \code{list(calls, G, bic)} where
#'   \code{calls} is a data.frame with columns \code{cell}, \code{stress_state},
#'   \code{stress_prob} (assignment confidence), and \code{prob_stressed}
#'   (posterior probability of not being in the most homeostatic state).
#'
#' @seealso \code{\link{call_mixture_states}}
#' @export
call_stress_states <- function(obj, cell_type, score_col = 'stress_composite',
                               G = 2:5, ...) {
  res <- call_mixture_states(obj, score_col = score_col,
                             group_col = 'annotation_first_pass',
                             group_value = cell_type, G = G,
                             label = 'stress', ...)
  if (is.null(res)) return(NULL)

  list(
    calls = data.frame(cell          = res$id,
                       stress_state  = res$stress_state,
                       stress_prob   = res$stress_confidence,
                       prob_stressed = res$stress_prob_high,
                       stringsAsFactors = FALSE),
    G   = attr(res, 'G'),
    bic = attr(res, 'bic')
  )
}
