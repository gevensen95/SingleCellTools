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
#' @export
#'
get_all_children <- function(go_term, verbose = FALSE) {
  if (verbose) message(sprintf('--- Collecting all child GO terms for %s ---', go_term))

  # Get direct children of the GO term. GOBPCHILDREN is a Bioconductor
  # AnnDbBimap object, not a plain environment -- base::mget() requires an
  # actual environment for its second argument and errors on this class
  # ("second argument must be an environment") on current GO.db releases.
  # AnnotationDbi provides an mget() S4 method that dispatches correctly
  # on Bimap objects like this one.
  children <- unlist(AnnotationDbi::mget(go_term, GO.db::GOBPCHILDREN, ifnotfound = NA))  # 'ifnotfound = NA' prevents errors if no children

  # If there are no children, return an empty vector
  if (all(is.na(children))) {
    return(NULL)
  }

  # Otherwise, recursively get children of each child term
  all_children <- children[!is.na(children)]  # Remove NAs
  for (child in all_children) {
    child_descendants <- get_all_children(child, verbose = FALSE)  # Recursive call (always silent)
    if (!is.null(child_descendants)) {
      all_children <- c(all_children, child_descendants)
    }
  }

  if (verbose) message(sprintf('  Total descendant terms collected: %d',
                               length(unique(all_children))))

  return(unique(all_children))
}
