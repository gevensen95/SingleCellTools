#' Recursive function to get all child terms for a given GO term
#'
#' This function will get all child GO terms
#'
#' @param go_term Parent GO term
#' @param verbose If TRUE, print a single high-level message before recursing.
#'   Defaults to FALSE because the function calls itself; emitting a message
#'   in every recursive call would flood the console. Recursive calls always
#'   pass \code{verbose = FALSE} internally regardless of the outer call's
#'   value, so passing \code{verbose = TRUE} at the top level prints exactly
#'   one message, not one per descendant term.
#' @return Vector containing child GO terms
#' @details Results are memoized in a session-scoped cache keyed by GO term
#'   id. The GO DAG lets a term have more than one parent, so a naive
#'   recursive walk can revisit -- and fully re-expand -- the same subtree
#'   many times both within one call (a term reachable from two different
#'   children of \code{go_term}) and across separate calls (two unrelated
#'   \code{go_term}s that share descendants). The cache persists for the
#'   life of the R session and assumes \code{GO.db}'s \code{GOBPCHILDREN}
#'   mapping doesn't change mid-session, which holds unless the annotation
#'   package itself is reinstalled/reloaded.
#' @export
#'
get_all_children <- function(go_term, verbose = FALSE) {
  if (verbose) message(sprintf('--- Collecting all child GO terms for %s ---', go_term))

  result <- .get_all_children_memo(go_term)

  if (verbose) message(sprintf('  Total descendant terms collected: %d',
                               length(result)))

  if (length(result) == 0) NULL else result
}

# Session-scoped memoization cache for .get_all_children_memo() below --
# keyed by GO term id, value is that term's full (already-deduped)
# descendant vector (character(0) for a leaf with no children). See the
# @details note on get_all_children() above for why this is safe and worth
# doing.
.go_children_cache <- new.env(parent = emptyenv())

#' @keywords internal
#' @noRd
.get_all_children_memo <- function(go_term) {
  # exists(..., inherits = FALSE) rather than a plain double-bracket lookup
  # against NULL -- character(0) (a real, cacheable "no children" result)
  # and "not cached yet" both look like NULL/empty under a double-bracket
  # lookup, so only an explicit existence check tells them apart.
  if (exists(go_term, envir = .go_children_cache, inherits = FALSE)) {
    return(get(go_term, envir = .go_children_cache, inherits = FALSE))
  }

  # Get direct children of the GO term. GOBPCHILDREN is a Bioconductor
  # AnnDbBimap object, not a plain environment -- base::mget() requires an
  # actual environment for its second argument and errors on this class
  # ("second argument must be an environment") on current GO.db releases.
  # AnnotationDbi provides an mget() S4 method that dispatches correctly
  # on Bimap objects like this one.
  children <- unlist(AnnotationDbi::mget(go_term, GO.db::GOBPCHILDREN, ifnotfound = NA))  # 'ifnotfound = NA' prevents errors if no children

  result <- if (all(is.na(children))) {
    character(0)
  } else {
    # Recursively get children of each child term
    all_children <- children[!is.na(children)]  # Remove NAs
    for (child in all_children) {
      child_descendants <- .get_all_children_memo(child)  # Recursive call (memoized)
      if (length(child_descendants)) {
        all_children <- c(all_children, child_descendants)
      }
    }
    unique(all_children)
  }

  assign(go_term, result, envir = .go_children_cache)
  result
}
