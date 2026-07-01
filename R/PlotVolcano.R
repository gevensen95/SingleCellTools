#' Volcano plot for DE results
#'
#' Standard volcano plot with configurable significance / effect
#' thresholds, top-gene labeling, and automatic column-name detection so
#' it accepts either \code{FindMarkers}-style output (\code{avg_log2FC},
#' \code{p_val_adj}) or \code{PseudobulkDE}-style output (\code{log2FC},
#' \code{padj}).
#'
#' @param de A data frame of DE results with, at minimum, a log2-fold-change
#'   column and an adjusted p-value column.
#' @param fc_col Name of the log2 fold change column. \code{NULL} (default)
#'   auto-detects among \code{"log2FC"}, \code{"avg_log2FC"},
#'   \code{"logFC"}, \code{"log2FoldChange"}.
#' @param p_col Name of the adjusted p-value column. \code{NULL} (default)
#'   auto-detects among \code{"padj"}, \code{"p_val_adj"}, \code{"FDR"},
#'   \code{"adj.P.Val"}.
#' @param gene_col Name of the gene-symbol column. \code{NULL} (default)
#'   auto-detects among \code{"gene"}, \code{"Gene"}, \code{"gene_symbol"}
#'   and falls back to rownames.
#' @param fc_threshold Absolute log2FC threshold for the "significant"
#'   coloring / labeling. Default 1 (2-fold change).
#' @param p_threshold Adjusted p-value threshold. Default 0.05.
#' @param top_n Number of top genes (by combined significance + effect
#'   size) to label. Default 20. Set to 0 to disable labeling.
#' @param label_genes Character vector of specific genes to label
#'   regardless of significance. Combined with the top-N selection.
#' @param colors Length-3 vector of colors: c(non-significant, up, down).
#'   Default \code{c("grey70", "#D6604D", "#4393C3")}.
#' @param point_size Base point size. Default 1.
#' @param label_size Text size for the labeled gene names. Default 3.
#' @param title Plot title. \code{NULL} auto-generates one from the DE
#'   metadata if present (\code{contrast_num}, \code{contrast_den},
#'   \code{cluster} columns).
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' # From PseudobulkDE
#' res <- PseudobulkDE(obj, sample_col = "orig.ident",
#'                     condition_col = "treatment",
#'                     ident_1 = "drug", ident_2 = "vehicle")
#' PlotVolcano(res$results)
#'
#' # From FindMarkers
#' fm <- FindMarkers(obj, ident.1 = "cluster1", ident.2 = "cluster2")
#' PlotVolcano(fm, gene_col = NULL)   # uses rownames
#'
#' # Custom labeling
#' PlotVolcano(res$results,
#'             label_genes  = c("MYC", "TP53"),
#'             fc_threshold = 1.5,
#'             top_n        = 10)
#' }
#' @importFrom ggplot2 aes element_text geom_hline geom_point geom_text geom_vline ggplot labs scale_color_manual theme theme_bw
#' @export
PlotVolcano <- function(de,
                        fc_col       = NULL,
                        p_col        = NULL,
                        gene_col     = NULL,
                        fc_threshold = 1,
                        p_threshold  = 0.05,
                        top_n        = 20,
                        label_genes  = NULL,
                        colors       = c("grey70", "#D6604D", "#4393C3"),
                        point_size   = 1,
                        label_size   = 3,
                        title        = NULL) {

  if (!is.data.frame(de)) stop("`de` must be a data frame.")
  if (length(colors) != 3) {
    stop("`colors` must be a length-3 vector: c(non-sig, up, down).")
  }

  # ---- Auto-detect columns ------------------------------------------------
  fc_col <- .detect_col(fc_col, de,
                        c("log2FC", "avg_log2FC", "logFC", "log2FoldChange"),
                        "log2 fold change")
  p_col  <- .detect_col(p_col, de,
                        c("padj", "p_val_adj", "FDR", "adj.P.Val"),
                        "adjusted p-value")
  gene_col_used <- gene_col
  if (is.null(gene_col_used)) {
    for (candidate in c("gene", "Gene", "gene_symbol")) {
      if (candidate %in% colnames(de)) {
        gene_col_used <- candidate
        break
      }
    }
  }
  if (is.null(gene_col_used)) {
    de$.gene <- rownames(de)
    gene_col_used <- ".gene"
  }

  # ---- Prepare plot data --------------------------------------------------
  df <- data.frame(
    gene   = as.character(de[[gene_col_used]]),
    log2FC = as.numeric(de[[fc_col]]),
    padj   = as.numeric(de[[p_col]]),
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$log2FC) & is.finite(df$padj) & df$padj > 0, ,
           drop = FALSE]
  df$negLog10P <- -log10(df$padj)

  df$category <- "n.s."
  df$category[df$padj < p_threshold & df$log2FC >   fc_threshold] <- "up"
  df$category[df$padj < p_threshold & df$log2FC < -fc_threshold] <- "down"
  df$category <- factor(df$category, levels = c("n.s.", "up", "down"))

  # ---- Which genes to label ----------------------------------------------
  labels_df <- df[integer(0), , drop = FALSE]
  if (top_n > 0) {
    sig <- df[df$category %in% c("up", "down"), ]
    if (nrow(sig)) {
      # Rank by combined criterion: -log10(padj) * |log2FC|
      sig$rank_score <- sig$negLog10P * abs(sig$log2FC)
      sig <- sig[order(-sig$rank_score), ]
      labels_df <- sig[seq_len(min(top_n, nrow(sig))), ]
    }
  }
  if (length(label_genes)) {
    extra <- df[df$gene %in% label_genes, ]
    labels_df <- unique(rbind(labels_df, extra))
  }

  # ---- Auto title ---------------------------------------------------------
  if (is.null(title)) {
    parts <- character(0)
    if (all(c("contrast_num", "contrast_den") %in% colnames(de))) {
      parts <- c(parts, sprintf("%s vs %s",
                                de$contrast_num[1], de$contrast_den[1]))
    }
    if ("cluster" %in% colnames(de) && length(unique(de$cluster)) == 1L) {
      parts <- c(parts, sprintf("cluster %s", de$cluster[1]))
    }
    if (length(parts)) title <- paste(parts, collapse = " — ")
  }

  # ---- Build plot ---------------------------------------------------------
  category <- log2FC <- negLog10P <- gene <- NULL  # NSE
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = log2FC, y = negLog10P,
                                    color = category)) +
    ggplot2::geom_point(size = point_size, alpha = 0.7) +
    ggplot2::geom_vline(xintercept = c(-fc_threshold, fc_threshold),
                        linetype = "dashed", color = "grey50") +
    ggplot2::geom_hline(yintercept = -log10(p_threshold),
                        linetype = "dashed", color = "grey50") +
    ggplot2::scale_color_manual(values = c("n.s." = colors[1],
                                           "up"   = colors[2],
                                           "down" = colors[3]),
                                name = NULL) +
    ggplot2::labs(x = expression(log[2]~"fold change"),
                  y = expression(-log[10]~"padj"),
                  title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  if (nrow(labels_df)) {
    p <- p + ggplot2::geom_text(
      data    = labels_df,
      mapping = ggplot2::aes(label = gene),
      size    = label_size,
      color   = "black",
      vjust   = -0.5,
      inherit.aes = TRUE
    )
  }

  p
}


#' @keywords internal
#' @noRd
.detect_col <- function(col, df, candidates, what) {
  if (!is.null(col)) {
    if (!col %in% colnames(df)) {
      stop("Column '", col, "' not found in data.")
    }
    return(col)
  }
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    stop("Could not auto-detect ", what,
         " column; tried: ", paste(candidates, collapse = ", "),
         ". Available columns: ", paste(colnames(df), collapse = ", "))
  }
  hit[1]
}
