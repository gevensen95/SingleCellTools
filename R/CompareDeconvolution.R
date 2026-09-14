#' Compare two spatial deconvolution results
#'
#' \code{\link{RunRCTD}}, \code{\link{RunCARD}}, and \code{\link{RunSPOTlight}}
#' each write a spot x cell-type proportion matrix to \code{obj@misc}
#' (\code{rctd_weights} / \code{card_weights} / \code{spotlight_weights}).
#' This compares any two of them the same way \code{\link{CompareMarkers}}/
#' \code{\link{CompareLRResults}} compare two result sets: a merged table,
#' correlation statistics, and a scatter plot.
#'
#' @param obj A Seurat object that has already been run through at least
#'   two of \code{\link{RunRCTD}}/\code{\link{RunCARD}}/\code{\link{RunSPOTlight}}.
#' @param method_a,method_b Which two methods to compare:
#'   \code{"rctd"}, \code{"card"}, \code{"spotlight"}. Must differ.
#' @param plot Logical; if \code{TRUE} (default), also return a scatter of
#'   matched proportions colored by cell type.
#' @return A list with \code{merged} (long-format data frame: \code{spot},
#'   \code{cell_type}, \code{prop_<method_a>}, \code{prop_<method_b>}),
#'   \code{correlation} (named numeric vector: \code{overall} Pearson r
#'   across every spot x cell-type pair, plus one entry per shared cell
#'   type), and optionally \code{plot}.
#' @examples
#' \dontrun{
#' visium <- RunRCTD(visium, reference = ref, celltype_col = "cell_type")
#' visium <- RunCARD(visium, reference = ref, celltype_col = "cell_type")
#' cmp <- CompareDeconvolution(visium, "rctd", "card")
#' cmp$correlation
#' cmp$plot
#' }
#' @importFrom ggplot2 aes facet_wrap geom_abline geom_point ggplot labs
#' @export
CompareDeconvolution <- function(obj,
                                 method_a = c("rctd", "card", "spotlight"),
                                 method_b = c("rctd", "card", "spotlight"),
                                 plot     = TRUE) {

  method_a <- match.arg(method_a)
  method_b <- match.arg(method_b)
  if (identical(method_a, method_b)) {
    stop("`method_a` and `method_b` must differ.")
  }
  .assert_seurat(obj)

  slot_name <- function(m) paste0(m, "_weights")
  get_weights <- function(m) {
    w <- obj@misc[[slot_name(m)]]
    if (is.null(w)) {
      stop("obj@misc$", slot_name(m), " not found -- run Run",
           toupper(substr(m, 1, 1)), substr(m, 2, nchar(m)), "() first.")
    }
    w
  }
  w_a <- get_weights(method_a)
  w_b <- get_weights(method_b)

  common_spots <- intersect(rownames(w_a), rownames(w_b))
  common_types <- intersect(colnames(w_a), colnames(w_b))
  if (length(common_spots) == 0) stop("No spots in common between the two methods.")
  if (length(common_types) == 0) stop("No cell types in common between the two methods.")

  a_sub <- w_a[common_spots, common_types, drop = FALSE]
  b_sub <- w_b[common_spots, common_types, drop = FALSE]

  merged <- data.frame(
    spot        = rep(common_spots, times = length(common_types)),
    cell_type   = rep(common_types, each = length(common_spots)),
    stringsAsFactors = FALSE
  )
  merged[[paste0("prop_", method_a)]] <- as.vector(a_sub)
  merged[[paste0("prop_", method_b)]] <- as.vector(b_sub)

  overall_r <- stats::cor(merged[[paste0("prop_", method_a)]],
                          merged[[paste0("prop_", method_b)]])
  per_type_r <- vapply(common_types, function(ct) {
    stats::cor(a_sub[, ct], b_sub[, ct])
  }, numeric(1))
  correlation <- c(overall = overall_r, per_type_r)

  message(sprintf("Overall correlation (%s vs %s): r = %.3f across %d spot x cell-type pairs.",
                  method_a, method_b, overall_r, nrow(merged)))

  out <- list(merged = merged, correlation = correlation)

  if (isTRUE(plot)) {
    x_col <- paste0("prop_", method_a)
    y_col <- paste0("prop_", method_b)
    cell_type <- NULL  # NSE
    merged$.x <- merged[[x_col]]
    merged$.y <- merged[[y_col]]
    p <- ggplot2::ggplot(merged, ggplot2::aes(x = .x, y = .y, color = cell_type)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
      ggplot2::geom_point(alpha = 0.4, size = 0.8) +
      ggplot2::labs(x = paste(method_a, "proportion"), y = paste(method_b, "proportion"),
                   color = NULL,
                   title = sprintf("r = %.3f", overall_r)) +
      Ol_Reliable()
    out$plot <- p
  }

  out
}
