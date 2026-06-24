#' Annotated Marker Gene Dot Plot
#'
#' Builds a Seurat \code{DotPlot} in which the features are organized by an
#' annotation (cell type, pathway, etc.) and the resulting plot is labeled
#' with those annotation groups along the right edge. Optionally clusters the
#' identities by correlation across the requested features before plotting.
#'
#' Before plotting, the gene table is filtered down to genes that:
#' \enumerate{
#'   \item are present in the chosen \code{assay},
#'   \item have non-zero expression in at least one cell of the object, and
#'   \item show variation in average expression across the identities being
#'     plotted (genes that scale to \code{NaN} because they're uniformly
#'     expressed would render as a row of grey dots and are dropped).
#' }
#' Each filtering step emits a message listing the dropped genes so all-empty
#' or all-grey rows never sneak into the plot.
#'
#' \strong{Dependency.} Uses the internal helpers \code{.pull_dotplot_data}
#' and \code{.cluster_mat} defined in \code{dotplot-helpers.R}. If sourcing
#' files directly (rather than loading as a package), source that file
#' before this one.
#'
#' @param obj A Seurat object.
#' @param genes A two-column data frame. The first column holds gene names,
#'   the second holds annotation labels (e.g. cell type, pathway). Column
#'   names are ignored — order is what matters.
#' @param margin_factor Scaling factor for plot margins. Increase if the
#'   annotation labels along the right edge are getting clipped.
#' @param maxsize Maximum dot size, passed to \code{scale_size_continuous}.
#' @param label.fontsize Font size for the right-edge annotation labels.
#' @param assay Assay name to read expression from. Default \code{"RNA"}.
#' @param show.annotations Logical; if TRUE (default), draw the annotation
#'   group labels along the right edge of the plot.
#' @param cluster Logical; if TRUE (default), cluster identities by
#'   correlation across the requested features before plotting.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' cellID <- data.frame(
#'   Gene    = c("Bank1", "Cd79a", "Ebf1", "Epcam", "Kcnma1", "Sox4"),
#'   Details = c("B cell", "B cell", "B cell",
#'               "Cholangiocyte", "Cholangiocyte", "Cholangiocyte")
#' )
#' MarkerPlot(harmony, cellID)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr pivot_wider
#' @importFrom ggplot2 aes coord_flip element_text geom_text geom_vline ggplot labs margin scale_color_gradientn scale_size_continuous theme theme_linedraw xlab ylab
#' @importFrom RColorBrewer brewer.pal
#' @importFrom Seurat DotPlot DefaultAssay DefaultAssay<- FetchData Idents Idents<-
#' @export
MarkerPlot <- function(obj,
                       genes,
                       margin_factor    = 0.5,
                       maxsize          = 4,
                       label.fontsize   = 3,
                       assay            = "RNA",
                       show.annotations = TRUE,
                       cluster          = TRUE) {

  # ---- Validate inputs ----------------------------------------------------
  if (!inherits(obj, "Seurat")) {
    stop("`obj` must be a Seurat object.")
  }
  if (!is.data.frame(genes) || ncol(genes) < 2) {
    stop("`genes` must be a data frame with at least two columns ",
         "(gene name, annotation).")
  }
  if (!assay %in% Seurat::Assays(obj)) {
    stop("Assay '", assay, "' not found in object. Available: ",
         paste(Seurat::Assays(obj), collapse = ", "))
  }

  # ---- Normalize the gene table -------------------------------------------
  genes <- stats::setNames(genes[, 1:2], c("Gene", "Details"))
  genes$Gene <- as.character(genes$Gene)
  genes$Details <- as.character(genes$Details)

  # Drop genes not present in the requested assay
  in_assay <- intersect(genes$Gene, rownames(obj[[assay]]))
  dropped_missing <- setdiff(genes$Gene, in_assay)
  if (length(dropped_missing)) {
    message(sprintf("  Dropping %d gene(s) not in assay '%s': %s",
                    length(dropped_missing), assay,
                    paste(dropped_missing, collapse = ", ")))
  }
  genes <- genes[genes$Gene %in% in_assay, , drop = FALSE]

  # Drop genes with zero expression across the entire object
  if (nrow(genes) > 0) {
    expr <- Seurat::FetchData(obj, vars = genes$Gene,
                              layer = "counts", assay = assay)
    has_expr <- colnames(expr)[colSums(expr, na.rm = TRUE) > 0]
    dropped_zero <- setdiff(genes$Gene, has_expr)
    if (length(dropped_zero)) {
      message(sprintf("  Dropping %d gene(s) with zero expression: %s",
                      length(dropped_zero),
                      paste(dropped_zero, collapse = ", ")))
    }
    genes <- genes[genes$Gene %in% has_expr, , drop = FALSE]
  }

  if (nrow(genes) == 0) {
    stop("No genes remain after filtering for assay presence and ",
         "non-zero expression.")
  }

  # ---- Pull dot plot summary data ----------------------------------------
  # Needed for the variance filter below and (optionally) for clustering.
  # This computes the same z-scored avg expression that Seurat::DotPlot will
  # later render, so we can preempt any features that scale to NaN.
  Seurat::DefaultAssay(obj) <- assay
  dp_data <- .pull_dotplot_data(obj, genes$Gene)

  # Drop genes whose avg.exp.scaled is non-finite across identities. These
  # are genes whose log1p(avg) is identical for every identity (commonly:
  # zero across every cluster that survives CellsByIdentities), so scale()
  # divides by zero and returns NaN. ggplot would render them as a row of
  # grey dots — drop them so the plot is clean.
  uniform_genes <- unique(as.character(
    dp_data$features.plot[!is.finite(dp_data$avg.exp.scaled)]
  ))
  if (length(uniform_genes)) {
    message(sprintf(
      "  Dropping %d gene(s) with no variation across identities (would render grey): %s",
      length(uniform_genes), paste(uniform_genes, collapse = ", ")
    ))
    genes   <- genes[!genes$Gene %in% uniform_genes, , drop = FALSE]
    dp_data <- dp_data[!as.character(dp_data$features.plot) %in% uniform_genes, ,
                       drop = FALSE]
  }

  if (nrow(genes) == 0) {
    stop("No genes remain after filtering for variation across identities.")
  }

  # ---- Annotation-label geometry ------------------------------------------
  tmp <- genes
  tmp$count <- 1
  Details <- count <- NULL  # silence R CMD check NSE notes
  tmp <- tmp %>%
    dplyr::group_by(Details) %>%
    dplyr::summarise(index = sum(count), .groups = "drop")
  tmp <- tmp[order(tmp$Details), ]
  tmp$cumsum  <- cumsum(tmp$index)
  tmp$annot_y <- tmp$cumsum - (0.5 * tmp$index)

  intersects <- cumsum(as.numeric(table(genes$Details))) + 0.5
  intersects <- intersects[seq_len(length(intersects) - 1)]

  # ---- Optionally cluster identities by correlation -----------------------
  if (isTRUE(cluster)) {
    message("Grouping clusters by correlation...")
    features.plot <- avg.exp.scaled <- id <- NULL  # NSE
    testmat <- tidyr::pivot_wider(
      dp_data[, c("id", "features.plot", "avg.exp.scaled")],
      names_from  = features.plot,
      values_from = avg.exp.scaled
    ) %>%
      tibble::column_to_rownames("id")

    # Any remaining non-finite cells (shouldn't be many after the uniform
    # filter, but be defensive) become 0 so hclust doesn't choke.
    testmat <- as.matrix(testmat)
    testmat[!is.finite(testmat)] <- 0

    tree_row  <- .cluster_mat(testmat)
    row_order <- rownames(testmat)[tree_row$order]
    Seurat::Idents(obj) <- factor(Seurat::Idents(obj), levels = row_order)
  }

  # ---- Build the plot -----------------------------------------------------
  ordered_genes <- genes$Gene[order(genes$Details)]
  tmp$xpos <- length(unique(Seurat::Idents(obj))) + 1

  d <- Seurat::DotPlot(
    obj,
    features = factor(ordered_genes, levels = ordered_genes)
  ) +
    ggplot2::scale_color_gradientn(
      colours = rev(RColorBrewer::brewer.pal(9, "RdBu"))
    ) +
    ggplot2::scale_size_continuous(range = c(0, maxsize)) +
    ggplot2::theme_linedraw() +
    ggplot2::geom_vline(xintercept = intersects) +
    ggplot2::xlab("") +
    ggplot2::ylab(" ") +
    ggplot2::coord_flip(clip  = "off",
                        ylim = c(1, length(unique(Seurat::Idents(obj))))) +
    ggplot2::theme(
      legend.position = "left",
      legend.title    = ggplot2::element_text(size = 6),
      axis.text.x     = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.margin     = ggplot2::margin(
        t = 1 * margin_factor,
        r = 5 * margin_factor,
        b = 1 * margin_factor,
        l = 1 * margin_factor,
        unit = "cm"
      )
    )

  if (isTRUE(show.annotations)) {
    annot_y <- xpos <- NULL  # NSE
    d <- d + ggplot2::geom_text(
      data    = tmp,
      mapping = ggplot2::aes(y = xpos, x = annot_y, label = Details),
      size    = label.fontsize,
      hjust   = 0,
      vjust   = 0
    )
  }

  d
}
