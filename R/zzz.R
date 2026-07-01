#' Attach SingleCellTools' dependencies on `library(SingleCellTools)`
#'
#' By default, R only \emph{loads} (but does not \emph{attach}) the packages
#' listed under \code{Imports} in DESCRIPTION -- their exported functions
#' are available to SingleCellTools' own code (which calls them via
#' \code{::}), but not directly to the user's session. This package is meant
#' to set up a full single-cell analysis environment, so on
#' \code{library(SingleCellTools)} this hook also attaches (i.e. puts on
#' \code{search()}) every package listed under \code{Imports} in
#' DESCRIPTION, so things like \code{ggplot()}, \code{%>%}, \code{filter()},
#' etc. are available without qualifying them.
#'
#' \strong{Notes / caveats:}
#' \itemize{
#'   \item Packages already on the search path (e.g. attached earlier by the
#'     user, or already attached by another package) are left as-is.
#'   \item If a listed dependency isn't installed, a single
#'     \code{packageStartupMessage()} names it and it is skipped --
#'     \code{library(SingleCellTools)} itself still succeeds.
#'   \item Attaching this many packages at once can introduce masking
#'     (e.g. \code{dplyr::filter} vs \code{stats::filter}); R prints the
#'     usual conflict messages, which are left visible (not suppressed)
#'     for that reason.
#'   \item \code{EnsDb.Mmusculus.v79} and \code{GO.db} are annotation/DB
#'     packages and can take a moment to load. To exclude any package from
#'     this auto-attach behavior, just remove it from \code{.scta_deps}
#'     below (it will still be available via \code{::} as long as it
#'     remains in DESCRIPTION's \code{Imports}).
#' }
#'
#' @param libname,pkgname Standard \code{.onAttach} arguments; supplied by R.
#' @keywords internal
#' @noRd
.scta_deps <- c(
  # Core Seurat / SeuratObject / Signac
  "Seurat",
  "SeuratObject",
  "Signac",
  "EnsDb.Mmusculus.v79",
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
  # Numerical / plotting
  "RANN",
  "ggplot2",
  "patchwork",         # QCComparePlots grid layout
  "RColorBrewer",
  "ClusterR",
  "irlba",
  "RSpectra",
  # Single-cell specific analysis
  "glmGamPoi",
  "GO.db",
  "UCell",
  "DoubletFinder",
  # Differential expression (PseudobulkDE)
  "DESeq2",
  "SummarizedExperiment",
  "S4Vectors",
  "edgeR",             # optional pseudobulk backend
  # Composition testing (CellComposition / CompositionalTest)
  "speckle",           # propeller
  "limma",             # speckle dependency; useful directly
  # Ligand-receptor (RunLIANA)
  "liana",
  "OmnipathR",         # liana pulls resources through OmnipathR
  # Reference-based annotation (AnnotateWithReference)
  "Azimuth",
  # Spatial deconvolution (RunRCTD)
  "spacexr",
  # Trajectory (PseudotimeWrapper)
  "slingshot",
  # Density plots (PlotFeatureDensity)
  "ks",
  "MASS",              # fallback KDE backend
  # Silhouette / integration QC (BatchEffectQC)
  "cluster",
  # Provenance sidecars (SaveWithProvenance)
  "jsonlite"
)

.onAttach <- function(libname, pkgname) {
  attached <- character(0)
  missing  <- character(0)

  for (pkg in .scta_deps) {
    if (paste0("package:", pkg) %in% search()) next

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
