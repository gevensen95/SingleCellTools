#' Plot gene-positivity results from AddGenePositivity
#'
#' Visualizes the per-cell logical positivity columns added by
#' \code{\link{AddGenePositivity}}. Three plot styles:
#'
#' \describe{
#'   \item{\code{"bar"} (default)}{Dodged bar chart of percent-positive per
#'     group, one bar per gene. Best when comparing a few genes across
#'     several groups.}
#'   \item{\code{"heatmap"}}{Group x gene tile chart, fill = percent
#'     positive. Best for many genes.}
#'   \item{\code{"combo"}}{Stacked bar chart of co-expression combination
#'     categories per group (e.g. \code{CD3D+/CD4-}, \code{CD3D+/CD4+},
#'     ...). Best for 2-3 genes.}
#' }
#'
#' Accepts a single Seurat object or a named list of them (same shape as
#' \code{AddGenePositivity}). For a list, per-sample percent-positive is
#' computed and rendered as facets (\code{"bar"} / \code{"heatmap"}) or as
#' faceted stacks (\code{"combo"}). Pass \code{sample_col} to instead
#' aggregate by a metadata column in a single combined object.
#'
#' @param obj A Seurat object or a named list of Seurat objects that
#'   already carry the per-gene \code{<gene><suffix>} logical metadata
#'   columns (i.e. run \code{AddGenePositivity} first).
#' @param genes Character vector of gene symbols. The function looks for
#'   \code{paste0(gene, suffix)} in each object's metadata.
#' @param group.by Metadata column used as the grouping axis (usually
#'   \code{"seurat_clusters"} or a cell-type label). Default
#'   \code{"seurat_clusters"}.
#' @param sample_col When \code{obj} is a single Seurat object containing
#'   multiple samples in one metadata column (e.g. \code{"orig.ident"}),
#'   the name of that column. Each sample gets its own facet. Ignored
#'   when \code{obj} is a list. Default \code{NULL}.
#' @param style One of \code{"bar"}, \code{"heatmap"}, \code{"combo"}.
#' @param suffix Suffix used by \code{AddGenePositivity} when creating the
#'   columns. Default \code{"_pos"}.
#' @param colors For \code{"bar"} and \code{"combo"}, a color palette for
#'   the discrete fill (gene name or combination category). \code{NULL}
#'   uses ggplot2's default. For \code{"heatmap"}, a length-2 vector of
#'   low/high gradient colors. Default \code{c("white", "firebrick")}.
#' @param drop_absent Logical; drop gene columns that are missing from
#'   \emph{any} object with a warning. Default TRUE. Set FALSE to instead
#'   fail with an error when a column is missing.
#' @param facet_ncol Number of facet columns for the sample facets. Default
#'   \code{NULL} (ggplot chooses).
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' # Single object
#' obj <- AddGenePositivity(obj, genes = c("CD3D", "CD4", "CD8A"))
#' PlotGenePositivity(obj, genes = c("CD3D", "CD4", "CD8A"))
#'
#' # Heatmap style, many genes
#' PlotGenePositivity(obj,
#'                    genes  = big_gene_vec,
#'                    style  = "heatmap",
#'                    colors = c("white", "steelblue4"))
#'
#' # Co-expression combinations per cluster
#' PlotGenePositivity(obj, genes = c("CD3D", "CD4"), style = "combo")
#'
#' # List of samples, faceted
#' list_of_objs <- AddGenePositivity(list_of_objs, c("CD3D", "CD4"))
#' PlotGenePositivity(list_of_objs, c("CD3D", "CD4"))
#'
#' # Single combined object, sample column
#' PlotGenePositivity(obj, c("CD3D", "CD4"),
#'                    sample_col = "orig.ident",
#'                    facet_ncol = 3)
#' }
#'
#' @importFrom ggplot2 aes coord_flip element_text facet_wrap geom_col geom_tile ggplot labs scale_fill_gradient scale_fill_manual theme theme_bw theme_linedraw
#' @export
PlotGenePositivity <- function(obj,
                               genes,
                               group.by    = "seurat_clusters",
                               sample_col  = NULL,
                               style       = c("bar", "heatmap", "combo"),
                               suffix      = "_pos",
                               colors      = NULL,
                               drop_absent = TRUE,
                               facet_ncol  = NULL) {

  style <- match.arg(style)
  stopifnot(is.character(genes), length(genes) >= 1)
  pos_cols <- paste0(genes, suffix)

  # ---- Build a long tidy data frame: sample, group, gene, positive -------
  df <- .collect_positivity(obj, pos_cols, genes,
                            group.by = group.by,
                            sample_col = sample_col,
                            drop_absent = drop_absent)

  if (nrow(df) == 0) {
    stop("No positivity data available; check that AddGenePositivity has ",
         "been run and that `genes` / `suffix` match the metadata columns.")
  }

  # ---- Dispatch to style --------------------------------------------------
  switch(
    style,
    bar     = .plot_bar(df,     colors, facet_ncol),
    heatmap = .plot_heatmap(df, colors, facet_ncol),
    combo   = .plot_combo(obj,  pos_cols, genes, group.by, sample_col,
                          drop_absent, colors, facet_ncol, suffix)
  )
}


# ============================================================================
# Internal: collect a long per-sample per-group per-gene percent-positive
# data frame from either a single object or a list.
# Columns: sample, group, gene, n_cells, n_pos, pct_pos
# ============================================================================

#' @keywords internal
#' @noRd
.collect_positivity <- function(obj, pos_cols, genes,
                                group.by, sample_col, drop_absent) {

  .one <- function(o, sample_name) {
    md <- o@meta.data
    if (!group.by %in% colnames(md)) {
      stop("`group.by` column '", group.by, "' not found in sample '",
           sample_name, "'.")
    }
    missing_cols <- setdiff(pos_cols, colnames(md))
    if (length(missing_cols)) {
      msg <- sprintf(
        "Sample '%s' is missing positivity column(s): %s. Did you run AddGenePositivity()?",
        sample_name, paste(missing_cols, collapse = ", "))
      if (!isTRUE(drop_absent)) stop(msg)
      warning(msg, call. = FALSE)
    }
    use_cols <- intersect(pos_cols, colnames(md))
    if (!length(use_cols)) return(NULL)

    grp <- as.character(md[[group.by]])
    out <- do.call(rbind, lapply(use_cols, function(pc) {
      gene <- sub(paste0(regmatches_escape(suffix_from(pc, genes)), "$"), "", pc)
      v <- md[[pc]]
      # aggregate per group
      by_grp <- split(as.logical(v), grp)
      data.frame(
        sample  = sample_name,
        group   = names(by_grp),
        gene    = gene,
        n_cells = vapply(by_grp, length, numeric(1)),
        n_pos   = vapply(by_grp, function(x) sum(x, na.rm = TRUE), numeric(1)),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }))
    out$pct_pos <- 100 * out$n_pos / pmax(1, out$n_cells)
    out
  }

  if (inherits(obj, "Seurat")) {
    if (!is.null(sample_col)) {
      if (!sample_col %in% colnames(obj@meta.data)) {
        stop("`sample_col` '", sample_col, "' not found in obj@meta.data.")
      }
      samples <- unique(as.character(obj@meta.data[[sample_col]]))
      dfs <- lapply(samples, function(samp) {
        cells <- rownames(obj@meta.data)[
          as.character(obj@meta.data[[sample_col]]) == samp
        ]
        .one(subset(obj, cells = cells), samp)
      })
      return(do.call(rbind, dfs))
    } else {
      out <- .one(obj, "all")
      out$sample <- NA_character_   # marker meaning "no sample facet"
      return(out)
    }
  }

  if (is.list(obj)) {
    if (is.null(names(obj))) {
      stop("List of Seurat objects must be named.")
    }
    dfs <- lapply(names(obj), function(nm) .one(obj[[nm]], nm))
    return(do.call(rbind, dfs))
  }

  stop("`obj` must be a Seurat object or a named list of Seurat objects.")
}


# Small utility to convert a suffixed column name back to a gene name by
# stripping whichever suffix produced it — safer than a fixed regex when
# the caller has an unusual suffix.
#' @keywords internal
#' @noRd
suffix_from <- function(colname, genes) {
  hit <- vapply(genes, function(g) startsWith(colname, g), logical(1))
  if (!any(hit)) return("")
  g <- genes[which(hit)[1]]
  sub(paste0("^", regmatches_escape(g)), "", colname)
}

# Escape regex meta-characters in a literal string.
#' @keywords internal
#' @noRd
regmatches_escape <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}


# ============================================================================
# Style: bar chart of percent-positive per group, dodged by gene.
# ============================================================================
#' @keywords internal
#' @noRd
.plot_bar <- function(df, colors, facet_ncol) {
  sample <- group <- gene <- pct_pos <- NULL  # NSE
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = group, y = pct_pos, fill = gene)) +
    ggplot2::geom_col(position = "dodge", color = "black", linewidth = 0.2) +
    ggplot2::labs(x = NULL, y = "% positive", fill = "Gene") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }
  .maybe_facet(p, df, facet_ncol)
}


# ============================================================================
# Style: heatmap tile, group x gene, fill = pct positive.
# ============================================================================
#' @keywords internal
#' @noRd
.plot_heatmap <- function(df, colors, facet_ncol) {
  sample <- group <- gene <- pct_pos <- NULL  # NSE
  if (is.null(colors)) colors <- c("white", "firebrick")
  if (length(colors) < 2) {
    stop("`colors` must be a length-2 low/high vector for style = 'heatmap'.")
  }
  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = group, y = gene, fill = pct_pos)) +
    ggplot2::geom_tile(color = "grey80") +
    ggplot2::scale_fill_gradient(low = colors[1], high = colors[2],
                                 limits = c(0, 100), name = "% positive") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_linedraw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  .maybe_facet(p, df, facet_ncol)
}


# ============================================================================
# Style: stacked bar of co-expression category (A+/B+, A+/B-, ...) per
# group. Requires re-computing the per-cell combination, so it takes the
# original obj rather than the aggregated df.
# ============================================================================
#' @keywords internal
#' @noRd
.plot_combo <- function(obj, pos_cols, genes, group.by, sample_col,
                        drop_absent, colors, facet_ncol, suffix) {

  .one_combo <- function(o, sample_name) {
    md <- o@meta.data
    missing_cols <- setdiff(pos_cols, colnames(md))
    if (length(missing_cols)) {
      msg <- sprintf(
        "Sample '%s' is missing positivity column(s): %s.",
        sample_name, paste(missing_cols, collapse = ", "))
      if (!isTRUE(drop_absent)) stop(msg)
      warning(msg, call. = FALSE)
    }
    use_cols <- intersect(pos_cols, colnames(md))
    if (!length(use_cols)) return(NULL)

    # Build a "combination" label per cell: e.g. "CD3D+/CD4-"
    combo_strs <- do.call(paste, c(lapply(use_cols, function(pc) {
      gene <- sub(paste0(regmatches_escape(suffix), "$"), "", pc)
      ifelse(as.logical(md[[pc]]), paste0(gene, "+"), paste0(gene, "-"))
    }), sep = " / "))

    data.frame(
      sample = sample_name,
      group  = as.character(md[[group.by]]),
      combo  = combo_strs,
      stringsAsFactors = FALSE
    )
  }

  if (inherits(obj, "Seurat")) {
    if (!is.null(sample_col)) {
      samples <- unique(as.character(obj@meta.data[[sample_col]]))
      combo_df <- do.call(rbind, lapply(samples, function(samp) {
        cells <- rownames(obj@meta.data)[
          as.character(obj@meta.data[[sample_col]]) == samp
        ]
        .one_combo(subset(obj, cells = cells), samp)
      }))
    } else {
      combo_df <- .one_combo(obj, "all")
      combo_df$sample <- NA_character_
    }
  } else if (is.list(obj)) {
    combo_df <- do.call(rbind, lapply(names(obj),
                                      function(nm) .one_combo(obj[[nm]], nm)))
  } else {
    stop("`obj` must be a Seurat object or a named list of Seurat objects.")
  }

  # Per (sample, group) proportions of each combo
  agg <- do.call(rbind, lapply(
    split(combo_df, list(combo_df$sample, combo_df$group), drop = TRUE),
    function(d) {
      tab <- table(d$combo)
      data.frame(
        sample = d$sample[1],
        group  = d$group[1],
        combo  = names(tab),
        prop   = as.numeric(tab) / sum(tab),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }))

  sample <- group <- combo <- prop <- NULL  # NSE
  p <- ggplot2::ggplot(agg,
                       ggplot2::aes(x = group, y = prop, fill = combo)) +
    ggplot2::geom_col(color = "black", linewidth = 0.2) +
    ggplot2::labs(x = NULL, y = "Proportion", fill = "Combination") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_fill_manual(values = colors)
  }
  .maybe_facet(p, agg, facet_ncol)
}


# ============================================================================
# Add facet_wrap(~ sample) when we have multiple samples; no-op otherwise.
# ============================================================================
#' @keywords internal
#' @noRd
.maybe_facet <- function(p, df, facet_ncol) {
  samples_present <- unique(df$sample[!is.na(df$sample)])
  if (length(samples_present) >= 2L) {
    p <- p + ggplot2::facet_wrap(~ sample, ncol = facet_ncol)
  }
  p
}
