#' Dense stacked violin plot
#'
#' One row per gene, one violin per cluster, narrow strip layout. The kind of
#' figure scanpy produces with \code{sc.pl.stacked_violin()} and Seurat
#' doesn't have a clean equivalent for out of the box.
#'
#' @param obj A Seurat object.
#' @param features Character vector of gene symbols to plot.
#' @param group.by Metadata column (or active idents if NULL) to group cells by.
#' @param assay Assay to read from. Default DefaultAssay.
#' @param layer Layer to read from. Default \code{"data"} (normalized).
#' @param colors Optional named character vector of colors per group level.
#' @param scale_per_gene Logical; if TRUE (default), each gene is scaled to
#'   its own max so low-expression genes are still visible.
#' @return A \code{ggplot} object.
#' @importFrom Seurat FetchData Idents DefaultAssay
#' @importFrom ggplot2 ggplot aes geom_violin facet_grid theme_classic theme element_blank element_text scale_fill_manual labs
#' @export
StackedViolinPlot <- function(obj,
                              features,
                              group.by       = NULL,
                              assay          = NULL,
                              layer          = "data",
                              colors         = NULL,
                              scale_per_gene = TRUE) {

  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay

  # Pull expression cells x genes
  features <- intersect(features, rownames(obj[[a]]))
  if (!length(features)) stop("None of the requested features are in assay '", a, "'.")
  expr <- Seurat::FetchData(obj, vars = features, layer = layer, assay = a)

  # Group vector
  if (is.null(group.by)) {
    grp <- as.character(Seurat::Idents(obj))
  } else {
    if (!group.by %in% colnames(obj@meta.data)) {
      stop("group.by '", group.by, "' not found in metadata.")
    }
    grp <- as.character(obj@meta.data[[group.by]])
  }
  expr$.group <- grp

  # Long format
  long <- stats::reshape(
    expr,
    direction = "long",
    varying   = list(features),
    v.names   = "expression",
    timevar   = "gene",
    times     = features,
    idvar     = ".rowid"
  )
  long$.rowid <- NULL
  long$gene <- factor(long$gene, levels = features)
  long$.group <- factor(long$.group)

  # Per-gene scaling so low-expression genes still show structure
  if (isTRUE(scale_per_gene)) {
    max_per_gene <- tapply(long$expression, long$gene, max, na.rm = TRUE)
    long$expression <- long$expression / max_per_gene[as.character(long$gene)]
  }

  expression <- .group <- gene <- NULL  # NSE silencing
  p <- ggplot2::ggplot(long, ggplot2::aes(x = .group, y = expression, fill = .group)) +
    ggplot2::geom_violin(scale = "width", trim = TRUE, color = NA) +
    ggplot2::facet_grid(rows = ggplot2::vars(gene), switch = "y") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1),
      strip.background  = ggplot2::element_blank(),
      axis.text.x       = ggplot2::element_text(angle = 45, hjust = 1),
      axis.title.y      = ggplot2::element_blank(),
      axis.text.y       = ggplot2::element_blank(),
      axis.ticks.y      = ggplot2::element_blank(),
      panel.spacing.y   = grid::unit(0, "lines"),
      legend.position   = "none"
    ) +
    ggplot2::labs(x = NULL)

  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }
  p
}
