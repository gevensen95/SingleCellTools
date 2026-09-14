#' Compare LIANA and CellChat ligand-receptor results
#'
#' \code{\link{RunLIANA}} and \code{\link{RunCellChat}} both infer
#' cell-cell communication on the same object, independently, with no way
#' to check where they agree -- the same gap \code{\link{CompareMarkers}}
#' already fills for two DE/marker result sets, applied here to
#' ligand-receptor inference instead.
#'
#' \strong{Matching caveat.} LIANA reports \code{ligand.complex}/
#' \code{receptor.complex} as possibly multi-subunit names joined with
#' \code{"_"} (e.g. \code{"TGFBR1_TGFBR2"}); CellChat's
#' \code{ligand}/\code{receptor} columns (from
#' \code{CellChat::subsetCommunication()}) are typically the primary
#' subunit only, sometimes with a \code{"+"}-joined complex name elsewhere.
#' This function matches on \code{source}, \code{target}, and each side's
#' \emph{first} gene symbol (everything before the first \code{"_"} or
#' \code{"+"}) -- a good-enough heuristic for "same interaction," but not a
#' guaranteed exact reconciliation when the two tools disagree about which
#' subunit is primary. Inspect \code{merged} directly if that distinction
#' matters for your analysis.
#'
#' @param liana_result The data frame returned by \code{\link{RunLIANA}}
#'   (i.e. \code{return_raw = FALSE}, the default there). Must have
#'   \code{source}, \code{target}, \code{ligand.complex},
#'   \code{receptor.complex}, \code{aggregate_rank}.
#' @param cellchat_result Either the \code{CellChat} object returned by
#'   \code{\link{RunCellChat}} (in which case
#'   \code{CellChat::subsetCommunication()} is called on it internally), or
#'   already the data frame that returns.
#' @param pval_cutoff Keep only CellChat interactions with \code{pval <=}
#'   this, if a \code{pval} column is present. \code{NULL} skips filtering.
#'   Default \code{0.05}.
#' @param top_n Keep only the top \code{top_n} LIANA interactions by
#'   \code{aggregate_rank} (lower = better) before comparing. \code{NULL}
#'   (default) keeps all of them.
#' @param plot Logical; if \code{TRUE} (default), also return a bar plot of
#'   how many interactions were found by both tools vs. only one.
#' @return A list with \code{merged} (one row per unique
#'   source/target/ligand/receptor combination seen in either tool, with
#'   \code{liana_rank}, \code{cellchat_prob}, \code{cellchat_pval}, and a
#'   \code{category} of \code{"both"} / \code{"only LIANA"} /
#'   \code{"only CellChat"}), \code{overlap} (named counts), and
#'   optionally \code{plot}.
#' @examples
#' \dontrun{
#' lr_liana    <- RunLIANA(obj, idents_col = "cell_type")
#' cc          <- RunCellChat(obj, label = "all")
#' cmp <- CompareLRResults(lr_liana, cc)
#' cmp$overlap
#' cmp$plot
#' }
#' @importFrom ggplot2 aes element_text geom_bar ggplot labs scale_fill_manual theme theme_bw
#' @export
CompareLRResults <- function(liana_result,
                             cellchat_result,
                             pval_cutoff = 0.05,
                             top_n       = NULL,
                             plot        = TRUE) {

  if (!is.data.frame(liana_result)) {
    stop("`liana_result` must be the data frame RunLIANA() returns ",
         "(return_raw = FALSE, the default).")
  }
  needed_liana <- c("source", "target", "ligand.complex", "receptor.complex",
                    "aggregate_rank")
  missing_liana <- setdiff(needed_liana, colnames(liana_result))
  if (length(missing_liana) > 0) {
    stop("`liana_result` is missing expected column(s): ",
         paste(missing_liana, collapse = ", "),
         ". Did this come from RunLIANA()?")
  }

  if (inherits(cellchat_result, "CellChat")) {
    if (!requireNamespace("CellChat", quietly = TRUE)) {
      stop("'CellChat' is required to read a CellChat object. Install with ",
           "devtools::install_github('sqjin/CellChat').")
    }
    cellchat_result <- CellChat::subsetCommunication(cellchat_result)
  }
  if (!is.data.frame(cellchat_result)) {
    stop("`cellchat_result` must be a CellChat object (as returned by ",
         "RunCellChat()) or the data frame from ",
         "CellChat::subsetCommunication().")
  }
  needed_cc <- c("source", "target", "ligand", "receptor")
  missing_cc <- setdiff(needed_cc, colnames(cellchat_result))
  if (length(missing_cc) > 0) {
    stop("`cellchat_result` is missing expected column(s): ",
         paste(missing_cc, collapse = ", "))
  }

  # Primary gene symbol from a possibly multi-subunit complex name -- see
  # the "Matching caveat" above.
  .primary_gene <- function(x) sub("[_+].*$", "", as.character(x))

  liana_norm <- data.frame(
    source     = as.character(liana_result$source),
    target     = as.character(liana_result$target),
    ligand     = .primary_gene(liana_result$ligand.complex),
    receptor   = .primary_gene(liana_result$receptor.complex),
    liana_rank = liana_result$aggregate_rank,
    stringsAsFactors = FALSE
  )
  if (!is.null(top_n)) {
    liana_norm <- liana_norm[order(liana_norm$liana_rank), , drop = FALSE]
    liana_norm <- liana_norm[seq_len(min(top_n, nrow(liana_norm))), , drop = FALSE]
  }

  has_prob <- "prob" %in% colnames(cellchat_result)
  has_pval <- "pval" %in% colnames(cellchat_result)
  if (!has_prob && !has_pval) {
    warning("Neither `prob` nor `pval` found in `cellchat_result` -- your ",
           "installed CellChat version may name these columns differently; ",
           "cellchat_prob/cellchat_pval will be NA.")
  }
  cellchat_norm <- data.frame(
    source        = as.character(cellchat_result$source),
    target        = as.character(cellchat_result$target),
    ligand        = .primary_gene(cellchat_result$ligand),
    receptor      = .primary_gene(cellchat_result$receptor),
    cellchat_prob = if (has_prob) cellchat_result$prob else NA_real_,
    cellchat_pval = if (has_pval) cellchat_result$pval else NA_real_,
    stringsAsFactors = FALSE
  )
  if (!is.null(pval_cutoff) && has_pval) {
    cellchat_norm <- cellchat_norm[
      is.na(cellchat_norm$cellchat_pval) | cellchat_norm$cellchat_pval <= pval_cutoff,
      , drop = FALSE]
  }

  key_cols <- c("source", "target", "ligand", "receptor")
  merged <- merge(liana_norm, cellchat_norm, by = key_cols, all = TRUE)
  merged$in_liana    <- !is.na(merged$liana_rank)
  merged$in_cellchat <- !is.na(merged$cellchat_prob) | !is.na(merged$cellchat_pval)
  merged$category <- factor(
    ifelse(merged$in_liana & merged$in_cellchat, "both",
    ifelse(merged$in_liana, "only LIANA", "only CellChat")),
    levels = c("both", "only LIANA", "only CellChat")
  )
  merged$in_liana <- merged$in_cellchat <- NULL
  merged <- merged[order(merged$category, merged$liana_rank), ]
  rownames(merged) <- NULL

  overlap <- c(
    both          = sum(merged$category == "both"),
    only_liana    = sum(merged$category == "only LIANA"),
    only_cellchat = sum(merged$category == "only CellChat"),
    n_total       = nrow(merged)
  )

  message(sprintf(
    "Both tools: %d. LIANA only: %d. CellChat only: %d (of %d total interactions).",
    overlap["both"], overlap["only_liana"], overlap["only_cellchat"],
    overlap["n_total"]))

  out <- list(merged = merged, overlap = overlap)

  if (isTRUE(plot)) {
    category <- NULL  # NSE
    p <- ggplot2::ggplot(merged, ggplot2::aes(x = category, fill = category)) +
      ggplot2::geom_bar() +
      ggplot2::scale_fill_manual(values = c(
        "both"          = "#4393C3",
        "only LIANA"    = "#8AB0D6",
        "only CellChat" = "#F4A261"
      )) +
      ggplot2::labs(x = NULL, y = "L-R interactions", fill = NULL,
                   title = "LIANA vs CellChat agreement") +
      Ol_Reliable()
    out$plot <- p
  }

  out
}
