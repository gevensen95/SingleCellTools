#' Annotated Marker Gene Dot Plot
#'
#' Builds a Seurat \code{DotPlot} in which the features are organized by an
#' annotation (cell type, pathway, etc.) and the resulting plot is labeled
#' with those annotation groups along the right edge. Optionally clusters the
#' identities by correlation across the requested features before plotting.
#'
#' Before plotting, the gene table is filtered down to genes that:
#' \enumerate{
#'   \item are present in the chosen \code{assay}, and
#'   \item have non-zero expression in at least one cell of the object.
#' }
#' Both filtering steps emit a message listing the dropped genes so all-empty
#' rows never sneak into the plot.
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
  Seurat::DefaultAssay(obj) <- assay
  if (isTRUE(cluster)) {
    message("Grouping clusters by correlation...")
    dp_data <- .pull_dotplot_data(obj, genes$Gene)
    features.plot <- avg.exp.scaled <- id <- NULL  # NSE
    testmat <- tidyr::pivot_wider(
      dp_data[, c("id", "features.plot", "avg.exp.scaled")],
      names_from  = features.plot,
      values_from = avg.exp.scaled
    ) %>%
      tibble::column_to_rownames("id")

    # Features with zero variance across identities (e.g. uniformly zero
    # expression for the surviving cells) come out of scale() as NaN and
    # would break hclust. Treat them as "no contribution" by zeroing them out.
    testmat <- as.matrix(testmat)
    testmat[!is.finite(testmat)] <- 0

    tree_row  <- .cluster_mat(testmat)
    row_order <- rownames(testmat)[tree_row$order]
    Seurat::Idents(obj) <- factor(Seurat::Idents(obj), levels = row_order)
  }

  # ---- Build the plot -----------------------------------------------------
  ordered_genes <- genes$Gene[order(genes$Details)]
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


# ============================================================================
# Internal helpers
# ============================================================================

#' Extract DotPlot summary data without rendering
#'
#' Internal helper that reproduces the data-prep half of \code{Seurat::DotPlot}
#' so the average-expression matrix can be used for clustering identities
#' before plotting. Not exported.
#'
#' @keywords internal
#' @noRd
.pull_dotplot_data <- function(object,
                               features,
                               assay           = NULL,
                               cols            = c("lightgrey", "blue"),
                               col.min         = -2.5,
                               col.max         = 2.5,
                               dot.min         = 0,
                               dot.scale       = 6,
                               idents          = NULL,
                               group.by        = NULL,
                               split.by        = NULL,
                               cluster.idents  = FALSE,
                               scale           = TRUE,
                               scale.by        = "radius",
                               scale.min       = NA,
                               scale.max       = NA) {
  if (is.null(assay)) assay <- Seurat::DefaultAssay(object)
  Seurat::DefaultAssay(object) <- assay

  split.colors <- !is.null(split.by) &&
    !any(cols %in% rownames(RColorBrewer::brewer.pal.info))
  scale.func <- switch(
    scale.by,
    size   = ggplot2::scale_size,
    radius = ggplot2::scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )

  feature.groups <- NULL
  if (is.list(features) || any(!is.na(names(features)))) {
    feature.groups <- unlist(sapply(seq_along(features), function(x) {
      rep(names(features)[x], each = length(features[[x]]))
    }))
    if (any(is.na(feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, immediate. = TRUE)
    }
    features <- unlist(features)
    names(feature.groups) <- features
  }

  cells <- unlist(Seurat::CellsByIdentities(
    object = object,
    cells  = colnames(object[[assay]]),
    idents = idents
  ))
  data.features <- Seurat::FetchData(object, vars = features, cells = cells)
  data.features$id <- if (is.null(group.by)) {
    Seurat::Idents(object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(data.features$id)) {
    data.features$id <- factor(data.features$id)
  }
  id.levels <- levels(data.features$id)
  data.features$id <- as.vector(data.features$id)

  if (!is.null(split.by)) {
    splits <- Seurat::FetchData(object, vars = split.by)[cells, split.by]
    if (split.colors) {
      if (length(unique(splits)) > length(cols)) {
        stop("Need to specify at least ", length(unique(splits)),
             " colors using the cols parameter")
      }
      cols <- cols[seq_along(unique(splits))]
      names(cols) <- unique(splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(splits)
    id.levels <- paste0(
      rep(id.levels, each = length(unique.splits)), "_",
      rep(unique(splits), times = length(id.levels))
    )
  }

  data.plot <- lapply(unique(data.features$id), function(ident) {
    data.use <- data.features[data.features$id == ident,
                              seq_len(ncol(data.features) - 1),
                              drop = FALSE]
    avg.exp <- apply(data.use, 2, function(x) mean(expm1(x)))
    pct.exp <- apply(data.use, 2, Seurat::PercentAbove, threshold = 0)
    list(avg.exp = avg.exp, pct.exp = pct.exp)
  })
  names(data.plot) <- unique(data.features$id)

  if (isTRUE(cluster.idents)) {
    mat <- do.call(rbind, lapply(data.plot, unlist))
    mat <- scale(mat)
    id.levels <- id.levels[hclust(dist(mat))$order]
  }

  data.plot <- lapply(names(data.plot), function(x) {
    data.use <- as.data.frame(data.plot[[x]])
    data.use$features.plot <- rownames(data.use)
    data.use$id <- x
    data.use
  })
  data.plot <- do.call(rbind, data.plot)
  if (!is.null(id.levels)) {
    data.plot$id <- factor(data.plot$id, levels = id.levels)
  }

  ngroup <- length(levels(data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will not be scaled",
            call. = FALSE, immediate. = TRUE)
  } else if (ngroup < 5 && scale) {
    warning("Scaling data with a low number of groups may produce misleading results",
            call. = FALSE, immediate. = TRUE)
  }

  avg.exp.scaled <- sapply(unique(data.plot$features.plot), function(x) {
    data.use <- data.plot[data.plot$features.plot == x, "avg.exp"]
    if (scale) {
      data.use <- scale(log1p(data.use))
      data.use <- Seurat::MinMax(data.use, min = col.min, max = col.max)
    } else {
      data.use <- log1p(data.use)
    }
    data.use
  })
  avg.exp.scaled <- as.vector(t(avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(cut(avg.exp.scaled, breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot  <- factor(data.plot$features.plot, levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100

  if (split.colors) {
    splits.use <- unlist(lapply(data.plot$id, function(x) {
      sub(paste0(".*_(",
                 paste(sort(unique(splits), decreasing = TRUE), collapse = "|"),
                 ")$"), "\\1", x)
    }))
    data.plot$colors <- mapply(function(color, value) {
      grDevices::colorRampPalette(colors = c("grey", color))(20)[value]
    }, color = cols[splits.use], value = avg.exp.scaled)
  }

  if (!is.na(scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(feature.groups)) {
    data.plot$feature.groups <- factor(
      feature.groups[data.plot$features.plot],
      levels = unique(feature.groups)
    )
  }

  data.plot
}


#' Hierarchically cluster a numeric matrix
#'
#' Internal helper supporting \code{MarkerPlot}. Defaults to correlation
#' distance with Ward's linkage. Not exported.
#'
#' @keywords internal
#' @noRd
.cluster_mat <- function(mat, distance = "correlation", method = "ward.D2") {
  valid_methods <- c("ward.D", "ward.D2", "ward", "single", "complete",
                     "average", "mcquitty", "median", "centroid")
  if (!(method %in% valid_methods)) {
    stop("clustering method has to be one from the list: ",
         paste(shQuote(valid_methods), collapse = ", "))
  }

  valid_distances <- c("correlation", "euclidean", "maximum", "manhattan",
                       "canberra", "binary", "minkowski")
  if (!(distance[1] %in% valid_distances) && !inherits(distance, "dist")) {
    stop("`distance` must be a 'dist' object or one of: ",
         paste(shQuote(valid_distances), collapse = ", "))
  }

  if (identical(distance[1], "correlation")) {
    cor_mat <- stats::cor(t(mat), use = "pairwise.complete.obs")
    # Any remaining NA (e.g. a row with zero variance vs. everything) becomes
    # 0 correlation so hclust doesn't choke on NA/NaN/Inf.
    cor_mat[!is.finite(cor_mat)] <- 0
    d <- stats::as.dist(1 - cor_mat)
  } else if (inherits(distance, "dist")) {
    d <- distance
  } else {
    d <- stats::dist(mat, method = distance)
  }

  stats::hclust(d, method = method)
}
