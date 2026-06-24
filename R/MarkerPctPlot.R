#' Annotated Percent-Positive Heatmap
#'
#' A sibling to \code{\link{MarkerPlot}} that plots, per gene and per identity,
#' the \strong{percent of cells expressing the gene} (count > 0). The features
#' are organized by an annotation (cell type, pathway, etc.) and the resulting
#' heatmap is labeled with those annotation groups along the right edge.
#' Optionally clusters the identities by percent-positive similarity across
#' the requested features before plotting.
#'
#' This complements \code{MarkerPlot}: where \code{MarkerPlot} encodes
#' \emph{average expression} via color and \emph{percent positive} via dot
#' size on a single dot plot, \code{MarkerPctPlot} isolates percent positive
#' as the only signal so it can be read off directly without being entangled
#' with avg-expression scaling.
#'
#' Before plotting, the gene table is filtered down to genes that:
#' \enumerate{
#'   \item are present in the chosen \code{assay},
#'   \item have non-zero expression in at least one cell of the object, and
#'   \item are detected (count > 0) in at least one cell of at least one
#'     identity being plotted (genes with 0\% positivity across every
#'     identity would render as a row of empty / minimum-color tiles and are
#'     dropped).
#' }
#' Each filtering step emits a message listing the dropped genes.
#'
#' @param obj A Seurat object.
#' @param genes A two-column data frame. The first column holds gene names,
#'   the second holds annotation labels. Column names are ignored — order is
#'   what matters.
#' @param margin_factor Scaling factor for plot margins. Increase if the
#'   annotation labels along the right edge are getting clipped.
#' @param label.fontsize Font size for the right-edge annotation labels.
#' @param assay Assay name to read expression from. Default \code{"RNA"}.
#' @param show.annotations Logical; if TRUE (default), draw the annotation
#'   group labels along the right edge of the plot.
#' @param cluster Logical; if TRUE (default), cluster identities by
#'   percent-positive similarity across the requested features before
#'   plotting.
#' @param min_pct Minimum percent-positive to display (smaller values are
#'   shown as the minimum color). Useful to dampen noise from tiny clusters.
#'   Default 0.
#' @param max_pct Maximum percent-positive for the color scale (values above
#'   are clipped). Default 100. Lower it (e.g. 50) to spread the color
#'   range over the more informative low/mid-percent region.
#' @param style \code{"tile"} (default) draws a heatmap with one rectangle
#'   per gene/identity, fill = percent positive. \code{"dot"} draws a dot
#'   plot whose size \emph{and} color both scale with percent positive.
#' @param colors A length-2 character vector of low/high colors for the
#'   gradient. Default \code{c("white", "firebrick")}.
#' @param maxsize (style = "dot" only) Maximum dot size passed to
#'   \code{scale_size_continuous}. Default 5.
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
#' MarkerPctPlot(harmony, cellID)
#' MarkerPctPlot(harmony, cellID, style = "dot", max_pct = 80)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr pivot_wider
#' @importFrom ggplot2 aes coord_flip element_text geom_text geom_tile geom_point geom_vline ggplot labs margin scale_fill_gradient scale_color_gradient scale_size_continuous theme theme_linedraw xlab ylab
#' @importFrom Seurat DefaultAssay DefaultAssay<- FetchData Idents Idents<- Assays
#' @export
MarkerPctPlot <- function(obj,
                          genes,
                          margin_factor    = 0.5,
                          label.fontsize   = 3,
                          assay            = "RNA",
                          show.annotations = TRUE,
                          cluster          = TRUE,
                          min_pct          = 0,
                          max_pct          = 100,
                          style            = c("tile", "dot"),
                          colors           = c("white", "firebrick"),
                          maxsize          = 5) {

  style <- match.arg(style)

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
  if (length(colors) != 2) {
    stop("`colors` must be a length-2 character vector (low, high).")
  }
  if (min_pct < 0 || max_pct > 100 || min_pct >= max_pct) {
    stop("Need 0 <= min_pct < max_pct <= 100.")
  }

  # ---- Normalize the gene table -------------------------------------------
  genes <- stats::setNames(genes[, 1:2], c("Gene", "Details"))
  genes$Gene    <- as.character(genes$Gene)
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

  # ---- Pull per-(ident, gene) percent-positive ----------------------------
  # Reuse the DotPlot data-prep helper so the percent-positive numbers are
  # computed exactly the way Seurat::DotPlot would compute them.
  Seurat::DefaultAssay(obj) <- assay
  dp_data <- .pull_dotplot_data(obj, genes$Gene)

  # pct.exp comes back on a 0-100 scale and is NA for any (gene, ident)
  # combination dropped by Seurat's own dot.min filter. Treat NA as 0
  # positivity for our purposes.
  dp_data$pct.exp[is.na(dp_data$pct.exp)] <- 0

  # Drop genes that are 0% positive in every identity (would render as a
  # row of empty tiles).
  zero_genes <- unique(as.character(
    dp_data$features.plot[ave(dp_data$pct.exp,
                              dp_data$features.plot,
                              FUN = max) == 0]
  ))
  if (length(zero_genes)) {
    message(sprintf(
      "  Dropping %d gene(s) with 0%% positivity across all identities: %s",
      length(zero_genes), paste(zero_genes, collapse = ", ")
    ))
    genes   <- genes[!genes$Gene %in% zero_genes, , drop = FALSE]
    dp_data <- dp_data[!as.character(dp_data$features.plot) %in% zero_genes, ,
                       drop = FALSE]
  }

  if (nrow(genes) == 0) {
    stop("No genes remain after filtering for non-zero positivity.")
  }

  # Clip to [min_pct, max_pct] for display
  dp_data$pct.exp[dp_data$pct.exp < min_pct] <- min_pct
  dp_data$pct.exp[dp_data$pct.exp > max_pct] <- max_pct

  # ---- Annotation-label geometry -----------------------------------------
  tmp <- genes
  tmp$count <- 1
  Details <- count <- NULL  # NSE
  tmp <- tmp %>%
    dplyr::group_by(Details) %>%
    dplyr::summarise(index = sum(count), .groups = "drop")
  tmp <- tmp[order(tmp$Details), ]
  tmp$cumsum  <- cumsum(tmp$index)
  tmp$annot_y <- tmp$cumsum - (0.5 * tmp$index)

  intersects <- cumsum(as.numeric(table(genes$Details))) + 0.5
  intersects <- intersects[seq_len(length(intersects) - 1)]

  # ---- Optionally cluster identities by pct similarity --------------------
  if (isTRUE(cluster)) {
    message("Grouping clusters by percent-positive correlation...")
    features.plot <- pct.exp <- id <- NULL  # NSE
    testmat <- tidyr::pivot_wider(
      dp_data[, c("id", "features.plot", "pct.exp")],
      names_from  = features.plot,
      values_from = pct.exp
    ) %>%
      tibble::column_to_rownames("id")

    testmat <- as.matrix(testmat)
    testmat[!is.finite(testmat)] <- 0

    tree_row  <- .cluster_mat(testmat)
    row_order <- rownames(testmat)[tree_row$order]
    dp_data$id <- factor(as.character(dp_data$id), levels = row_order)
    Seurat::Idents(obj) <- factor(Seurat::Idents(obj), levels = row_order)
  } else {
    # Preserve current Idents() ordering on the y axis.
    dp_data$id <- factor(as.character(dp_data$id),
                         levels = levels(Seurat::Idents(obj)))
  }

  # ---- Order features by annotation group --------------------------------
  ordered_genes <- genes$Gene[order(genes$Details)]
  dp_data$features.plot <- factor(as.character(dp_data$features.plot),
                                  levels = ordered_genes)

  n_idents <- length(unique(dp_data$id))
  tmp$xpos <- n_idents + 1

  # ---- Build the plot -----------------------------------------------------
  features.plot <- pct.exp <- id <- annot_y <- xpos <- NULL  # NSE

  if (style == "tile") {
    p <- ggplot2::ggplot(
      dp_data,
      ggplot2::aes(x = id, y = features.plot, fill = pct.exp)
    ) +
      ggplot2::geom_tile(color = "grey70") +
      ggplot2::scale_fill_gradient(
        low    = colors[1],
        high   = colors[2],
        limits = c(min_pct, max_pct),
        name   = "% positive"
      )
  } else {
    p <- ggplot2::ggplot(
      dp_data,
      ggplot2::aes(x = id, y = features.plot,
                   size = pct.exp, color = pct.exp)
    ) +
      ggplot2::geom_point() +
      ggplot2::scale_size_continuous(
        range  = c(0, maxsize),
        limits = c(min_pct, max_pct),
        name   = "% positive"
      ) +
      ggplot2::scale_color_gradient(
        low    = colors[1],
        high   = colors[2],
        limits = c(min_pct, max_pct),
        name   = "% positive",
        guide  = "none"
      )
  }

  p <- p +
    ggplot2::theme_linedraw() +
    ggplot2::geom_vline(xintercept = intersects) +
    ggplot2::xlab("") +
    ggplot2::ylab(" ") +
    ggplot2::coord_flip(clip = "off", ylim = c(1, n_idents)) +
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

  # Note: coord_flip() swaps x/y, so the per-identity axis ends up horizontal
  # and the per-gene axis ends up vertical — matching MarkerPlot's layout.
  # The annotation labels go on the (post-flip) right edge, which means
  # mapping to (y = xpos, x = annot_y) in pre-flip coordinates.
  if (isTRUE(show.annotations)) {
    p <- p + ggplot2::geom_text(
      data    = tmp,
      mapping = ggplot2::aes(x = xpos, y = annot_y, label = Details),
      size    = label.fontsize,
      hjust   = 0,
      vjust   = 0,
      inherit.aes = FALSE
    )
  }

  p
}
