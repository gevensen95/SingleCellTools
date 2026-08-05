#' Attach SingleCellTools' core dependencies on `library(SingleCellTools)`
#'
#' By default, R only \emph{loads} (but does not \emph{attach}) the packages
#' listed under \code{Imports} in DESCRIPTION -- their exported functions
#' are available to SingleCellTools' own code (which calls them via
#' \code{::}), but not directly to the user's session. This package is meant
#' to set up a core single-cell analysis environment, so on
#' \code{library(SingleCellTools)} this hook also attaches (i.e. puts on
#' \code{search()}) the small set of packages listed in \code{.scta_deps}
#' below -- the ones actually worth having unqualified in an interactive
#' session (\code{Seurat}/\code{Signac}, the tidyverse core, \code{ggplot2}
#' and friends) -- so things like \code{ggplot()}, \code{%>%},
#' \code{filter()}, etc. are available without qualifying them.
#'
#' \strong{This is a deliberately curated subset of \code{Imports}, not all
#' of it.} The many single-function dependencies (\code{DESeq2},
#' \code{slingshot}, \code{scmap}, \code{GO.db}, \code{UCell},
#' \code{DoubletFinder}, \code{MASS}, etc.) are intentionally left
#' unattached -- they're only ever used inside one or two specific
#' SingleCellTools functions (via \code{::}), so forcing them onto every
#' user's search path would slow down \code{library(SingleCellTools)} and
#' risk name-masking (e.g. \code{MASS::select}, \code{edgeR}/\code{limma}
#' exports) for no benefit to sessions that never call those functions.
#'
#' \strong{Notes / caveats:}
#' \itemize{
#'   \item Packages already on the search path (e.g. attached earlier by the
#'     user, or already attached by another package) are left as-is.
#'   \item If a listed dependency isn't installed, a single
#'     \code{packageStartupMessage()} names it and it is skipped --
#'     \code{library(SingleCellTools)} itself still succeeds.
#'   \item R prints the usual conflict messages for any masking these
#'     packages introduce (e.g. \code{dplyr::filter} vs \code{stats::filter});
#'     these are left visible (not suppressed).
#'   \item To attach an additional package yourself, just add it to
#'     \code{.scta_deps} below (it will still be available via \code{::}
#'     either way, as long as it remains in DESCRIPTION's \code{Imports}).
#' }
#'
#' @param libname,pkgname Standard \code{.onAttach} arguments; supplied by R.
#' @keywords internal
#' @noRd
.scta_deps <- c(
  # Core Seurat / Signac
  "Seurat",
  "SeuratObject",
  "Signac",
  # Base R namespaces we lean on
  "Matrix",
  "methods",
  "grDevices",
  "grid",
  # tidyverse tools
  "dplyr",
  "tibble",
  "tidyr",
  "magrittr",
  "readr",
  "stringr",
  "rlang",
  "purrr",
  # Numerical / plotting used pervasively across plotting functions
  "RANN",
  "ggplot2",
  "patchwork",         # QCComparePlots grid layout
  "RColorBrewer"
)

.onAttach <- function(libname, pkgname) {
  attached <- character(0)
  missing  <- character(0)

  # Verbose per-package progress -- lets a hang or crash during this loop
  # be pinpointed to a single dependency instead of just "library() never
  # returns" / "R died somewhere in here". Opt out with
  # options(SingleCellTools.quiet_attach = TRUE).
  verbose_attach <- !isTRUE(getOption("SingleCellTools.quiet_attach", FALSE))

  for (pkg in .scta_deps) {
    if (paste0("package:", pkg) %in% search()) next

    if (verbose_attach) packageStartupMessage("  ...loading ", pkg)

    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
      next
    }

    # attachNamespace() errors if a package is already attached; the
    # search()-path check above already guards against that, but stay
    # defensive in case of races with other startup hooks.
    attached_ok <- tryCatch({
      attachNamespace(pkg)
      TRUE
    }, error = function(e) FALSE)

    if (attached_ok) attached <- c(attached, pkg)
  }

  if (length(missing)) {
    packageStartupMessage(
      pkgname, ": ", length(missing),
      " optional dependenc", if (length(missing) == 1) "y is" else "ies are",
      " not installed and will not be attached: ",
      paste(missing, collapse = ", "))
  }

  if (length(attached)) {
    packageStartupMessage(
      pkgname, ": attached ", length(attached),
      " dependenc", if (length(attached) == 1) "y" else "ies",
      ":\n  ", paste(strwrap(paste(attached, collapse = ", "),
                             width = 70, exdent = 2), collapse = "\n"))
  }
}
