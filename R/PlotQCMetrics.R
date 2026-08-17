#' Plot QC metrics across samples
#'
#' Builds a multi-panel QC figure from a Seurat object's metadata, one panel
#' per requested metric, grouped along the x-axis by \code{group.by}
#' (default \code{"orig.ident"}). Numeric columns (e.g. \code{nFeature_RNA},
#' \code{percent.mt}) are drawn as boxplots; categorical columns (e.g. a
#' doublet-calling column) are drawn as stacked proportion bars labeled with
#' the per-group count in each segment, matching the doublet panel style
#' used in \code{CreateRNAObjects}.
#'
#' By default (\code{qc_cols = NULL}) the columns to plot are auto-detected:
#' the first metadata column matching \code{^nFeature}, the first matching
#' \code{^nCount}, \code{percent.mt}, and a doublet-calling column
#' (\code{doublet_finder} if present, otherwise the first column matching
#' \code{"doublet"} case-insensitively). Any of these that aren't present
#' are silently skipped rather than erroring, since not every object will
#' have gone through doublet calling, have a percent.mt column, etc.
#'
#' @param obj A Seurat object, or a list of Seurat objects. For a list,
#'   only each object's \code{@meta.data} is row-bound together (via
#'   \code{dplyr::bind_rows}) rather than merging the objects themselves,
#'   since the count matrices aren't needed here -- this keeps things fast
#'   even for many/large samples.
#' @param group.by Metadata column to group samples along the x-axis.
#'   Default \code{"orig.ident"}.
#' @param qc_cols Character vector of metadata column names to plot instead
#'   of the auto-detected defaults. If supplied, only these columns are
#'   used (auto-detection is skipped entirely); any requested column not
#'   found in the object's metadata is skipped with a message.
#' @param ncol Optional number of columns for the panel layout. Default
#'   \code{NULL} lets patchwork lay panels out in a single row.
#' @return A patchwork/ggplot object. Not printed automatically -- assign
#'   it or let it auto-print at the console.
#' @examples
#' \dontrun{
#' PlotQCMetrics(obj)
#' PlotQCMetrics(obj, group.by = "Treatment")
#' PlotQCMetrics(obj, qc_cols = c("nFeature_RNA", "percent.mt", "Phase"))
#' }
#' @importFrom ggplot2 aes element_text geom_boxplot geom_col geom_text ggplot labs position_fill position_stack theme
#' @export
PlotQCMetrics <- function(obj, group.by = "orig.ident", qc_cols = NULL, ncol = NULL) {

  if (is.list(obj) && !inherits(obj, "Seurat")) {
    if (length(obj) == 0) {
      stop("`obj` is an empty list.")
    }
    if (!all(vapply(obj, inherits, logical(1), "Seurat"))) {
      stop("`obj` must be a Seurat object or a list of Seurat objects.")
    }
    # Row-bind just the metadata rather than Seurat::merge()-ing the objects
    # themselves -- merge() would combine the underlying count matrices too,
    # which this function never touches, so it's wasted work (and memory)
    # for large/many-sample lists. bind_rows() also tolerates samples that
    # are missing a given QC column (e.g. no doublet_finder yet), filling
    # NA rather than erroring.
    meta <- dplyr::bind_rows(lapply(obj, function(x) x@meta.data))
  } else {
    if (!inherits(obj, "Seurat")) {
      stop("`obj` must be a Seurat object or a list of Seurat objects.")
    }
    meta <- obj@meta.data
  }

  if (!group.by %in% colnames(meta)) {
    stop("`group.by` column '", group.by, "' was not found in the object's metadata.")
  }

  if (is.null(qc_cols)) {
    feature_col <- grep("^nFeature", colnames(meta), value = TRUE)[1]
    count_col   <- grep("^nCount",   colnames(meta), value = TRUE)[1]
    mito_col    <- if ("percent.mt" %in% colnames(meta)) "percent.mt" else NA_character_
    doublet_col <- if ("doublet_finder" %in% colnames(meta)) {
      "doublet_finder"
    } else {
      hits <- grep("doublet", colnames(meta), ignore.case = TRUE, value = TRUE)
      if (length(hits)) hits[1] else NA_character_
    }
    qc_cols <- c(feature_col, count_col, mito_col, doublet_col)
    qc_cols <- qc_cols[!is.na(qc_cols)]
    if (length(qc_cols) == 0) {
      stop("Could not auto-detect any QC columns (nFeature*/nCount*/percent.mt/doublet*) ",
           "in the object's metadata. Supply `qc_cols` explicitly.")
    }
  } else {
    missing_cols <- setdiff(qc_cols, colnames(meta))
    if (length(missing_cols)) {
      message("Skipping QC column(s) not found in metadata: ",
              paste(missing_cols, collapse = ", "))
    }
    qc_cols <- intersect(qc_cols, colnames(meta))
    if (length(qc_cols) == 0) {
      stop("None of the requested `qc_cols` were found in the object's metadata.")
    }
  }

  message(sprintf("--- Plotting QC metrics (%s) grouped by '%s' ---",
                  paste(qc_cols, collapse = ", "), group.by))

  .numeric_panel <- function(col) {
    ggplot2::ggplot(meta, ggplot2::aes(x = .data[[group.by]], y = .data[[col]])) +
      ggplot2::geom_boxplot() +
      ggplot2::labs(title = col, x = NULL, y = col) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      Ol_Reliable()
  }

  .categorical_panel <- function(col) {
    pct <- Freq <- NULL  # silence R CMD check NSE notes
    # useNA = "ifany" -- meta[[col]] can legitimately be NA for some rows
    # (e.g. a mixed list of objects where only some went through doublet
    # calling; dplyr::bind_rows() fills the missing column with NA rather
    # than erroring, per this function's own docs). table()'s default
    # useNA = "no" would silently drop those rows from both the counts and
    # each group's sum(Freq) denominator below, understating group size
    # with no indication anything was excluded.
    counts_df <- as.data.frame(table(
      group_ = meta[[group.by]],
      level_ = meta[[col]],
      useNA  = "ifany"
    ))
    names(counts_df) <- c(group.by, col, "Freq")
    counts_df <- counts_df |>
      dplyr::group_by(.data[[group.by]]) |>
      dplyr::mutate(pct = Freq / sum(Freq)) |>
      dplyr::ungroup()

    ggplot2::ggplot(counts_df, ggplot2::aes(x = .data[[group.by]], y = pct,
                                            fill = .data[[col]])) +
      ggplot2::geom_col(position = ggplot2::position_stack()) +
      ggplot2::geom_text(
        ggplot2::aes(label = ifelse(Freq > 0, Freq, "")),
        position = ggplot2::position_fill(vjust = 0.5), size = 3
      ) +
      ggplot2::labs(title = col, x = NULL, y = "Proportion", fill = NULL) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      Ol_Reliable()
  }

  panels <- lapply(qc_cols, function(col) {
    if (is.numeric(meta[[col]])) .numeric_panel(col) else .categorical_panel(col)
  })

  combined <- Reduce(`+`, panels)
  if (!is.null(ncol)) {
    combined <- combined + patchwork::plot_layout(ncol = ncol)
  }
  combined
}
