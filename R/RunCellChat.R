#' Cell-cell communication inference with CellChat
#'
#' Wraps the standard \code{CellChat} pipeline (\code{createCellChat} ->
#' \code{subsetDB} -> \code{subsetData} ->
#' \code{identifyOverExpressedGenes}/\code{identifyOverExpressedInteractions}
#' -> \code{computeCommunProb} -> \code{filterCommunication} ->
#' \code{computeCommunProbPathway} -> \code{aggregateNet}) into a single
#' call, for one Seurat subset (e.g. one condition, sample, or region) at a
#' time.
#'
#' @param seurat_sub A Seurat object (typically a subset already restricted
#'   to the condition/sample/region you want to analyze).
#' @param label A short label used only in the summary message printed at
#'   the end (e.g. the sample or condition name), so repeated calls in a
#'   loop are easy to tell apart in the console.
#' @param group_col Metadata column naming the groups CellChat should treat
#'   as signaling sources/targets. Default \code{"signaling_group"}.
#' @param assay Assay to read expression from. Default
#'   \code{DefaultAssay(seurat_sub)}.
#' @param species \code{"mouse"} (default) or \code{"human"}; selects
#'   \code{CellChatDB.mouse} / \code{CellChatDB.human} as the
#'   ligand-receptor database.
#' @param signaling_type Value(s) passed to \code{CellChat::subsetDB(...,
#'   search = signaling_type)} to restrict the database to one or more
#'   signaling categories (e.g. \code{"Secreted Signaling"},
#'   \code{"ECM-Receptor"}, \code{"Cell-Cell Contact"}). Pass \code{NULL}
#'   to use the entire database unfiltered. Default \code{"Secreted
#'   Signaling"}.
#' @param comm_prob_type Method passed to \code{computeCommunProb(type =
#'   ...)}. Default \code{"triMean"}.
#' @param nboot Number of bootstrap iterations used by
#'   \code{computeCommunProb()} for significance testing. Default 25.
#' @param min_cells Minimum cells per group required to keep an
#'   interaction, passed to \code{filterCommunication(min.cells = ...)}.
#'   Default 10.
#' @param verbose Logical; print the per-call summary message (groups /
#'   significant L-R pairs). Default \code{TRUE}.
#' @return The fitted \code{CellChat} object (post \code{aggregateNet()}),
#'   with \code{cc@net}, \code{cc@netP}, and \code{cc@LR$LRsig} populated.
#' @examples
#' \dontrun{
#' cc_healthy <- RunCellChat(subset(obj, condition == "healthy"),
#'                           label = "healthy")
#' cc_tumor   <- RunCellChat(subset(obj, condition == "tumor"),
#'                           label = "tumor", nboot = 100)
#' }
#' @importFrom Seurat DefaultAssay
#' @importFrom utils data
#' @export
RunCellChat <- function(seurat_sub,
                        label,
                        group_col      = "signaling_group",
                        assay          = NULL,
                        species        = c("mouse", "human"),
                        signaling_type = "Secreted Signaling",
                        comm_prob_type = "triMean",
                        nboot          = 25,
                        min_cells      = 10,
                        verbose        = TRUE) {

  species <- match.arg(species)

  # Cheap, local checks first so they surface their own errors even on a
  # machine without CellChat installed.
  if (!inherits(seurat_sub, "Seurat")) {
    stop("`seurat_sub` must be a Seurat object.")
  }
  if (!group_col %in% colnames(seurat_sub@meta.data)) {
    stop("`group_col` '", group_col, "' not found in seurat_sub@meta.data.")
  }
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("'CellChat' is required. Install with ",
         "devtools::install_github('sqjin/CellChat').")
  }

  a <- if (is.null(assay)) Seurat::DefaultAssay(seurat_sub) else assay

  # ---- Load the requested ligand-receptor database -------------------------
  # utils::data(..., envir = environment()) loads the package dataset into
  # this call's local environment without attaching CellChat to the search
  # path, so CellChatDB.mouse/.human don't need to already be exported.
  db_name <- paste0("CellChatDB.", species)
  utils::data(list = db_name, package = "CellChat", envir = environment())
  db <- get(db_name, envir = environment())

  message(sprintf(
    "--- Running CellChat (%s, group.by = '%s', assay = '%s', species = %s) ---",
    label, group_col, a, species))

  cc <- CellChat::createCellChat(object = seurat_sub, group.by = group_col,
                                 assay = a)
  cc@DB <- if (is.null(signaling_type)) {
    db
  } else {
    CellChat::subsetDB(db, search = signaling_type)
  }
  cc <- CellChat::subsetData(cc)
  cc <- CellChat::identifyOverExpressedGenes(cc)
  cc <- CellChat::identifyOverExpressedInteractions(cc)
  cc <- CellChat::computeCommunProb(cc, type = comm_prob_type, nboot = nboot)
  cc <- CellChat::filterCommunication(cc, min.cells = min_cells)
  cc <- CellChat::computeCommunProbPathway(cc)
  cc <- CellChat::aggregateNet(cc)

  if (isTRUE(verbose)) {
    message(sprintf("%s: %d signaling groups, %d significant L-R pairs.",
                    label, length(levels(cc@idents)), nrow(cc@LR$LRsig)))
  }

  cc
}
