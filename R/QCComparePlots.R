#' Compare QC metrics before and after filtering
#'
#' Draws violin overlays of QC metrics from paired pre- and post-filter
#' Seurat objects (or lists of them) so the effect of \code{ApplyQCFilters}
#' — or any other cell-level filtering step — can be inspected at a
#' glance. One violin per (sample, filter-state) shows the distribution
#' shift; a text annotation reports \code{n_before -> n_after} per sample.
#'
#' Pairs naturally with \code{GenerateQCReport} +
#' \code{ApplyQCFilters}: run the report to pick cutoffs, apply them, then
#' call this on the (pre, post) pair to confirm the filter did what you
#' expected.
#'
#' @param pre A Seurat object or (named) list of Seurat objects
#'   \emph{before} filtering.
#' @param post A Seurat object or list of Seurat objects \emph{after}
#'   filtering. Must have the same shape and names as \code{pre}.
#' @param metrics Character vector of metadata column names to plot. If
#'   \code{NULL} (default), auto-detects any of \code{"nCount_RNA"},
#'   \code{"nFeature_RNA"}, \code{"percent.mt"}, \code{"percent.rb"},
#'   \code{"percent.hb"}, \code{"nCount_Spatial"},
#'   \code{"nFeature_Spatial"} that exist in the objects.
#' @param sample_col When \code{pre}/\code{post} are single Seurat objects
#'   with multiple samples in one column, that column's name. Default
#'   \code{NULL} (treat each single-object input as one sample).
#' @param log_y Metrics whose distributions to display on \code{log10}
#'   scale (skewed by nature). Default \code{c("nCount_RNA",
#'   "nFeature_RNA", "nCount_Spatial", "nFeature_Spatial")}.
#' @param colors A length-2 vector of colors for \code{pre} / \code{post}.
#'   Default \code{c("#8AB0D6", "#F4A261")}.
#' @param ncol Number of columns in the patchwork grid over metrics.
#'   Default \code{NULL} (patchwork chooses).
#' @param show_counts Logical; if TRUE (default) annotate each sample with
#'   \code{n_before -> n_after (pct_kept\%)}.
#' @return A \code{ggplot} object (single metric) or a patchwork of them
#'   (multiple metrics).
#' @examples
#' \dontrun{
#' # Report + filter workflow
#' GenerateQCReport(sample_list, output_file = "qc/qc_report.html")
#' sample_list_filt <- ApplyQCFilters(sample_list,
#'                                    cutoffs = "qc/qc_report_cutoffs.csv")
#' QCComparePlots(sample_list, sample_list_filt)
#'
#' # Single combined object
#' QCComparePlots(pre_obj, post_obj, sample_col = "orig.ident")
#' }
#' @importFrom ggplot2 aes annotate element_text geom_violin ggplot labs scale_fill_manual scale_y_log10 theme theme_bw
#' @export
QCComparePlots <- function(pre,
                           post,
                           metrics     = NULL,
                           sample_col  = NULL,
                           log_y       = c("nCount_RNA", "nFeature_RNA",
                                           "nCount_Spatial", "nFeature_Spatial"),
                           colors      = c("#8AB0D6", "#F4A261"),
                           ncol        = NULL,
                           show_counts = TRUE) {

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("`patchwork` not installed; multi-metric grid will be returned ",
            "as a list of ggplots.", call. = FALSE)
  }

  df <- .qc_gather(pre, post, sample_col)
  if (nrow(df) == 0) stop("No metadata rows found in `pre` / `post`.")

  # Determine metrics
  standard <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb",
                "percent.hb", "nCount_Spatial", "nFeature_Spatial")
  if (is.null(metrics)) {
    metrics <- intersect(standard, colnames(df))
    if (!length(metrics)) {
      stop("No standard QC metrics found; pass `metrics` explicitly.")
    }
  } else {
    metrics <- intersect(metrics, colnames(df))
    if (!length(metrics)) {
      stop("None of the requested metrics were found in `pre`/`post`.")
    }
  }

  # Per-sample retention annotation
  counts <- .qc_counts(df)

  plots <- lapply(metrics, function(m) {
    .qc_violin(df, m, log_y, colors, counts, show_counts)
  })
  names(plots) <- metrics

  if (length(plots) == 1L) return(plots[[1]])
  if (requireNamespace("patchwork", quietly = TRUE)) {
    Reduce(`+`, plots) + patchwork::plot_layout(ncol = ncol,
                                                guides = "collect")
  } else {
    plots
  }
}


# ============================================================================
# Internal: gather metadata from pre + post into one long data frame with
# a `state` column ("pre" / "post") and a `sample` column.
# ============================================================================
#' @keywords internal
#' @noRd
.qc_gather <- function(pre, post, sample_col) {

  .one <- function(o, sample_name) {
    if (!inherits(o, "Seurat")) {
      stop("Object for '", sample_name, "' is not a Seurat object.")
    }
    md <- o@meta.data
    md$sample <- sample_name
    md$.cell  <- rownames(md)
    md
  }

  .collect <- function(x, state) {
    if (inherits(x, "Seurat")) {
      if (!is.null(sample_col)) {
        if (!sample_col %in% colnames(x@meta.data)) {
          stop("`sample_col` '", sample_col, "' not in metadata.")
        }
        samples <- unique(as.character(x@meta.data[[sample_col]]))
        parts <- lapply(samples, function(samp) {
          cells <- rownames(x@meta.data)[
            as.character(x@meta.data[[sample_col]]) == samp
          ]
          .one(subset(x, cells = cells), samp)
        })
        d <- do.call(rbind, parts)
      } else {
        d <- .one(x, "all")
      }
    } else if (is.list(x)) {
      if (is.null(names(x))) stop("List inputs must be named.")
      d <- do.call(rbind,
                   lapply(names(x), function(nm) .one(x[[nm]], nm)))
    } else {
      stop("`pre`/`post` must be a Seurat object or named list.")
    }
    d$state <- state
    d
  }

  d1 <- .collect(pre,  "pre")
  d2 <- .collect(post, "post")
  # rbind with column intersection so mismatched extra metadata cols don't
  # prevent stacking.
  common <- intersect(colnames(d1), colnames(d2))
  d1 <- d1[, common, drop = FALSE]
  d2 <- d2[, common, drop = FALSE]
  d <- rbind(d1, d2)
  d$state <- factor(d$state, levels = c("pre", "post"))
  d
}


# ============================================================================
# Internal: n_before / n_after / pct_kept per sample
# ============================================================================
#' @keywords internal
#' @noRd
.qc_counts <- function(df) {
  tab <- as.data.frame(table(sample = df$sample, state = df$state),
                       stringsAsFactors = FALSE)
  wide <- reshape(tab, idvar = "sample", timevar = "state",
                  direction = "wide")
  colnames(wide) <- sub("^Freq\\.", "", colnames(wide))
  wide$pct_kept <- round(100 * wide$post / pmax(1, wide$pre), 1)
  wide$label <- sprintf("%d -> %d\n(%.1f%%)",
                        wide$pre, wide$post, wide$pct_kept)
  wide
}


# ============================================================================
# Internal: violin plot for one metric, pre/post side by side, per sample
# ============================================================================
#' @keywords internal
#' @noRd
.qc_violin <- function(df, metric, log_y, colors, counts, show_counts) {
  sample <- state <- NULL  # NSE
  df$.y <- df[[metric]]
  df <- df[is.finite(df$.y), , drop = FALSE]

  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = sample, y = .y, fill = state)) +
    ggplot2::geom_violin(scale = "width", position = ggplot2::position_dodge(0.8),
                         width = 0.7, color = "black", linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = c(pre = colors[1], post = colors[2]),
                               name = NULL) +
    ggplot2::labs(x = NULL, y = metric) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  if (metric %in% log_y) {
    p <- p + ggplot2::scale_y_log10()
  }

  if (isTRUE(show_counts) && nrow(counts)) {
    y_pos <- if (metric %in% log_y) {
      max(df$.y, na.rm = TRUE) ^ 1.02
    } else {
      max(df$.y, na.rm = TRUE) * 1.05
    }
    p <- p + ggplot2::annotate(
      "text",
      x     = counts$sample,
      y     = y_pos,
      label = counts$label,
      size  = 2.5, vjust = 0
    )
  }
  p
}
