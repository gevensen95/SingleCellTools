# ============================================================================
# Shared internal helpers for the MarkerPlot / MarkerPctPlot family of
# functions. Both MarkerPlot.R and MarkerPctPlot.R depend on the two
# functions defined here; keeping them in a single file avoids duplication
# and the need for `if (!exists(...))` guards.
#
# Source order (when not using as a package):
#   source("dotplot-helpers.R")
#   source("MarkerPlot.R")
#   source("MarkerPctPlot.R")
#
# Inside an R package the file order doesn't matter — every file in R/ is
# loaded before any function is called.
# ============================================================================


#' Extract DotPlot summary data without rendering
#'
#' Internal helper that reproduces the data-prep half of
#' \code{Seurat::DotPlot} so the per-(ident, gene) summary table is
#' available for downstream plotting. Returns a data frame with one row per
#' (identity, gene) combination and columns \code{avg.exp},
#' \code{avg.exp.scaled}, \code{pct.exp}, \code{features.plot}, \code{id}.
#'
#' Not exported. Used by \code{\link{MarkerPlot}} and
#' \code{\link{MarkerPctPlot}}.
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
      warning("Some feature groups are unnamed.", call. = FALSE,
              immediate. = TRUE)
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
  data.plot$features.plot  <- factor(data.plot$features.plot,
                                     levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100

  if (split.colors) {
    splits.use <- unlist(lapply(data.plot$id, function(x) {
      sub(paste0(".*_(",
                 paste(sort(unique(splits), decreasing = TRUE),
                       collapse = "|"),
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
#' Internal helper supporting \code{\link{MarkerPlot}} and
#' \code{\link{MarkerPctPlot}}. Defaults to correlation distance with
#' Ward's linkage. Not exported.
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


#' Normalize and validate a (Gene, Details) marker table against an assay
#'
#' Shared prep step for \code{\link{MarkerPlot}} and
#' \code{\link{MarkerPctPlot}} -- previously ~150 lines duplicated verbatim
#' across both files (gene-table validation, duplicate-gene dedup,
#' assay-presence filtering, zero-expression filtering, with identical
#' comments), extracted here so a future fix to this logic can't drift
#' between the two copies the way it previously could.
#'
#' Validates that \code{obj} is a Seurat object, \code{genes} is a data
#' frame with at least two columns, and \code{assay} exists in \code{obj}.
#' Keeps only \code{genes}'s first two columns (renamed \code{Gene}/
#' \code{Details}), drops genes listed under more than one Details group
#' (keeping the first occurrence -- see inline comment below for why a
#' duplicate gene would otherwise crash deep inside
#' \code{.pull_dotplot_data()}), drops genes not present in \code{assay},
#' and drops genes with zero expression across the whole object (via the
#' \code{"counts"} layer). Errors if nothing survives. Each dropped
#' category is reported via \code{message()}.
#'
#' Not exported.
#'
#' @keywords internal
#' @noRd
.prepare_marker_gene_table <- function(obj, genes, assay) {
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

  genes <- stats::setNames(genes[, 1:2], c("Gene", "Details"))
  genes$Gene    <- as.character(genes$Gene)
  genes$Details <- as.character(genes$Details)

  # `genes` maps each gene to exactly one Details group. A gene listed under
  # two different Details groups (e.g. a marker genuinely shared between two
  # subtypes/cell types) reaches .pull_dotplot_data()'s own
  # `factor(features.plot, levels = features)` call with a duplicated level
  # and crashes with a cryptic, hard-to-trace "factor level [n] is
  # duplicated" error -- fail predictably here instead: keep the first
  # occurrence (its first-listed Details group) and warn about the rest,
  # rather than silently guessing or crashing deep inside plot internals.
  dup_genes <- unique(genes$Gene[duplicated(genes$Gene)])
  if (length(dup_genes) > 0) {
    message(sprintf(
      "  %d gene(s) listed under more than one Details group -- keeping only the first occurrence: %s",
      length(dup_genes), paste(dup_genes, collapse = ", ")))
    genes <- genes[!duplicated(genes$Gene), , drop = FALSE]
  }

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

  genes
}
