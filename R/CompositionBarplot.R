#' Stacked / grouped bar plot of cell-type composition
#'
#' Companion to \code{\link{CompositionAnalysis}}. Takes either a Seurat
#' object (in which case it calls \code{CompositionAnalysis} internally) or
#' the \code{counts} or \code{proportions} tibble from that function.
#'
#' @param x A Seurat object, OR a tibble of the form returned by
#'   \code{CompositionAnalysis()$proportions}.
#' @param group_col,sample_col,condition_col Used only when \code{x} is a
#'   Seurat object &mdash; passed to \code{CompositionAnalysis()}.
#' @param style One of \code{"stacked"} (default) or \code{"grouped"}.
#' @param y One of \code{"prop"} (default) or \code{"n_cells"}.
#' @param facet_by_condition Logical; if TRUE and a condition column exists,
#'   facet by condition.
#' @param colors Optional named character vector of colors per group level.
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_col facet_wrap theme_minimal theme element_text labs scale_fill_manual position_dodge position_stack
#' @export
CompositionBarplot <- function(x,
                               group_col          = NULL,
                               sample_col         = NULL,
                               condition_col      = NULL,
                               style              = c("stacked", "grouped"),
                               y                  = c("prop", "n_cells"),
                               facet_by_condition = TRUE,
                               colors             = NULL) {
  style <- match.arg(style)
  y     <- match.arg(y)

  # Resolve input to a data frame
  if (inherits(x, "Seurat")) {
    if (is.null(group_col) || is.null(sample_col)) {
      stop("`group_col` and `sample_col` are required when passing a Seurat object.")
    }
    res <- CompositionAnalysis(x, group_col = group_col,
                               sample_col = sample_col,
                               condition_col = condition_col,
                               test = "none")
    df <- res$proportions
  } else {
    df <- as.data.frame(x)
    required <- c("sample", "group", "n_cells", "prop")
    missing  <- setdiff(required, colnames(df))
    if (length(missing)) {
      stop("Input is missing columns: ", paste(missing, collapse = ", "))
    }
  }

  position <- if (style == "stacked") ggplot2::position_stack() else
                                      ggplot2::position_dodge(width = 0.9)

  sample <- group <- n_cells <- prop <- NULL  # NSE
  p <- ggplot2::ggplot(df, ggplot2::aes(x = sample,
                                        y = if (y == "prop") prop else n_cells,
                                        fill = group)) +
    ggplot2::geom_col(position = position) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = NULL,
                  y = if (y == "prop") "Proportion" else "Number of cells",
                  fill = NULL)

  if (!is.null(colors)) p <- p + ggplot2::scale_fill_manual(values = colors)

  if (isTRUE(facet_by_condition) && "condition" %in% colnames(df)) {
    condition <- NULL
    p <- p + ggplot2::facet_wrap(~ condition, scales = "free_x")
  }

  p
}
