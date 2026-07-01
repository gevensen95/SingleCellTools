#' Compare two DE / marker result sets
#'
#' Given two DE result data frames (e.g. from \code{FindMarkers} run in
#' two conditions, or from \code{PseudobulkDE} between two contrasts),
#' returns a merged per-gene table, overlap statistics, and a scatter of
#' log2FC values. Useful for asking questions like "which markers are
#' shared, which are unique, do the shared ones agree in direction."
#'
#' @param a,b Two data frames of DE results. Must contain the columns
#'   named by \code{fc_col}, \code{p_col}, \code{gene_col} (or the
#'   auto-detected equivalents; see \code{PlotVolcano}).
#' @param labels Length-2 character labels for \code{a} / \code{b}. Used
#'   as column suffixes and axis labels.
#' @param fc_col Log2 fold change column. \code{NULL} auto-detects.
#' @param p_col Adjusted p-value column. \code{NULL} auto-detects.
#' @param gene_col Gene column. \code{NULL} auto-detects.
#' @param fc_threshold |log2FC| threshold for calling a gene significant.
#'   Default 1.
#' @param p_threshold Adjusted p-value threshold. Default 0.05.
#' @param plot Logical; if TRUE (default) also returns a scatter of
#'   log2FC values from a vs b, colored by significance category.
#' @return A list with \code{merged} (per-gene merged data frame with
#'   \code{log2FC_<a>}, \code{log2FC_<b>}, \code{padj_<a>}, \code{padj_<b>},
#'   \code{category}), \code{overlap} (integer vector of counts:
#'   \code{shared_up}, \code{shared_down}, \code{opposite_sign},
#'   \code{only_a}, \code{only_b}), \code{fisher} (Fisher's exact test
#'   result on the overlap contingency), and optionally \code{plot}.
#' @examples
#' \dontrun{
#' res_drug    <- PseudobulkDE(obj, sample_col = "orig.ident",
#'                             condition_col = "treatment",
#'                             ident_1 = "drug", ident_2 = "vehicle")
#' res_disease <- PseudobulkDE(obj, sample_col = "orig.ident",
#'                             condition_col = "status",
#'                             ident_1 = "disease", ident_2 = "healthy")
#' cmp <- CompareMarkers(res_drug$results, res_disease$results,
#'                       labels = c("drug_vs_vehicle", "disease_vs_healthy"))
#' cmp$overlap
#' cmp$plot
#' }
#' @importFrom ggplot2 aes element_text geom_hline geom_point geom_smooth geom_vline ggplot labs scale_color_manual theme theme_bw
#' @export
CompareMarkers <- function(a, b,
                           labels       = c("a", "b"),
                           fc_col       = NULL,
                           p_col        = NULL,
                           gene_col     = NULL,
                           fc_threshold = 1,
                           p_threshold  = 0.05,
                           plot         = TRUE) {

  if (!is.data.frame(a) || !is.data.frame(b)) {
    stop("`a` and `b` must both be data frames.")
  }
  if (length(labels) != 2) {
    stop("`labels` must be a length-2 character vector.")
  }

  fc_a <- .detect_col(fc_col, a,
                      c("log2FC", "avg_log2FC", "logFC", "log2FoldChange"),
                      "log2 fold change")
  fc_b <- .detect_col(fc_col, b,
                      c("log2FC", "avg_log2FC", "logFC", "log2FoldChange"),
                      "log2 fold change")
  p_a  <- .detect_col(p_col, a,
                      c("padj", "p_val_adj", "FDR", "adj.P.Val"),
                      "adjusted p-value")
  p_b  <- .detect_col(p_col, b,
                      c("padj", "p_val_adj", "FDR", "adj.P.Val"),
                      "adjusted p-value")
  g_a <- if (is.null(gene_col)) {
    intersect(c("gene", "Gene", "gene_symbol"), colnames(a))[1]
  } else gene_col
  g_b <- if (is.null(gene_col)) {
    intersect(c("gene", "Gene", "gene_symbol"), colnames(b))[1]
  } else gene_col
  if (is.na(g_a)) { a$.gene <- rownames(a); g_a <- ".gene" }
  if (is.na(g_b)) { b$.gene <- rownames(b); g_b <- ".gene" }

  df_a <- data.frame(gene = a[[g_a]], log2FC = a[[fc_a]], padj = a[[p_a]],
                     stringsAsFactors = FALSE)
  df_b <- data.frame(gene = b[[g_b]], log2FC = b[[fc_b]], padj = b[[p_b]],
                     stringsAsFactors = FALSE)
  colnames(df_a) <- c("gene",
                      paste0("log2FC_", labels[1]),
                      paste0("padj_",   labels[1]))
  colnames(df_b) <- c("gene",
                      paste0("log2FC_", labels[2]),
                      paste0("padj_",   labels[2]))

  merged <- merge(df_a, df_b, by = "gene", all = FALSE)

  # Significance calls
  sig_a <- merged[[paste0("padj_", labels[1])]] < p_threshold &
           abs(merged[[paste0("log2FC_", labels[1])]]) > fc_threshold
  sig_b <- merged[[paste0("padj_", labels[2])]] < p_threshold &
           abs(merged[[paste0("log2FC_", labels[2])]]) > fc_threshold
  sig_a[is.na(sig_a)] <- FALSE
  sig_b[is.na(sig_b)] <- FALSE

  same_dir <- sign(merged[[paste0("log2FC_", labels[1])]]) ==
              sign(merged[[paste0("log2FC_", labels[2])]])

  category <- rep("n.s. in both", nrow(merged))
  category[sig_a & !sig_b] <- paste0("only ", labels[1])
  category[!sig_a & sig_b] <- paste0("only ", labels[2])
  category[sig_a & sig_b & same_dir &
           merged[[paste0("log2FC_", labels[1])]] > 0]  <- "shared up"
  category[sig_a & sig_b & same_dir &
           merged[[paste0("log2FC_", labels[1])]] < 0]  <- "shared down"
  category[sig_a & sig_b & !same_dir]                   <- "opposite sign"
  merged$category <- factor(category, levels = c(
    "n.s. in both", paste0("only ", labels[1]), paste0("only ", labels[2]),
    "shared up", "shared down", "opposite sign"
  ))

  overlap <- c(
    shared_up     = sum(category == "shared up"),
    shared_down   = sum(category == "shared down"),
    opposite_sign = sum(category == "opposite sign"),
    only_a        = sum(category == paste0("only ", labels[1])),
    only_b        = sum(category == paste0("only ", labels[2])),
    n_tested_both = nrow(merged)
  )

  # Fisher on 2x2 shared vs unique
  mat <- matrix(c(sum(sig_a & sig_b),  sum(sig_a & !sig_b),
                  sum(!sig_a & sig_b), sum(!sig_a & !sig_b)),
                nrow = 2)
  fisher <- stats::fisher.test(mat)

  out <- list(merged = merged, overlap = overlap, fisher = fisher)

  if (isTRUE(plot)) {
    log2FC_a_name <- paste0("log2FC_", labels[1])
    log2FC_b_name <- paste0("log2FC_", labels[2])
    log2FC_a <- log2FC_b <- category <- NULL  # NSE
    merged$log2FC_a <- merged[[log2FC_a_name]]
    merged$log2FC_b <- merged[[log2FC_b_name]]
    p <- ggplot2::ggplot(merged,
                         ggplot2::aes(x = log2FC_a, y = log2FC_b,
                                      color = category)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "grey60") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                          color = "grey60") +
      ggplot2::geom_point(alpha = 0.6, size = 0.8) +
      ggplot2::geom_smooth(inherit.aes = FALSE,
                           mapping = ggplot2::aes(x = log2FC_a,
                                                  y = log2FC_b),
                           method = "lm", se = FALSE,
                           color = "black", linewidth = 0.4,
                           linetype = "dotted") +
      ggplot2::scale_color_manual(values = c(
        "n.s. in both"                      = "grey80",
        setNames("#8AB0D6", paste0("only ", labels[1])),
        setNames("#F4A261", paste0("only ", labels[2])),
        "shared up"                         = "#D6604D",
        "shared down"                       = "#4393C3",
        "opposite sign"                     = "purple"
      )) +
      ggplot2::labs(x = log2FC_a_name, y = log2FC_b_name, color = NULL) +
      ggplot2::theme_bw()
    out$plot <- p
  }

  out
}
