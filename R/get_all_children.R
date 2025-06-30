#' Recursive function to get all child terms for a given GO term
#'
#' This function will get all child GO terms
#'
#' @param go_term Parent GO term
#' @return Vector containing child GO terms
#' @export
#'
get_all_children <- function(go_term) {
  # Get direct children of the GO term
  children <- unlist(mget(go_term, GOBPCHILDREN, ifnotfound = NA))  # 'ifnotfound = NA' prevents errors if no children

  # If there are no children, return an empty vector
  if (all(is.na(children))) {
    return(NULL)
  }

  # Otherwise, recursively get children of each child term
  all_children <- children[!is.na(children)]  # Remove NAs
  for (child in all_children) {
    child_descendants <- get_all_children(child)  # Recursive call for each child
    if (!is.null(child_descendants)) {
      all_children <- c(all_children, child_descendants)
    }
  }

  return(unique(all_children))
}
