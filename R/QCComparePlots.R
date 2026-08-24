#' Compare QC metrics before and after filtering
#'
#' Draws violin overlays of QC metrics from paired pre- and post-filter
#' Seurat objects (or lists of them) so the effect of \code{ApplyQCFilters}
#' -- or any other cell-level filtering step -- can be inspected at a
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
#'   \code{n_before}, \code{-> n_after} and \code{(pct_kept\%)} stacked on
#'   three lines above that sample's violins.
#' @param counts_size Text size for the \code{show_counts} annotation, in
#'   \code{ggplot2} size units. Default \code{2.5}. The y axis is expanded to
#'   fit the annotation automatically, but with many samples in one panel the
#'   labels can still crowd each other horizontally -- lower this (or pass
#'   \code{show_counts = FALSE}) if they do.
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
#' @importFrom ggplot2 aes annotate coord_cartesian element_text expansion geom_violin ggplot labs margin scale_fill_manual scale_y_continuous scale_y_log10 theme theme_bw
#' @export
QCComparePlots <- function(pre,
                           post,
                           metrics     = NULL,
                           sample_col  = NULL,
                           log_y       = c("nCount_RNA", "nFeature_RNA",
                                           "nCount_Spatial", "nFeature_Spatial"),
                           colors      = c("#8AB0D6", "#F4A261"),
                           ncol        = NULL,
                           show_counts = TRUE,
                           counts_size = 2.5) {

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
    .qc_violin(df, m, log_y, colors, counts, show_counts, counts_size)
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

  # rbind a list of per-object metadata frames on their COMMON columns only.
  # Seurat metadata routinely carries object-specific column NAMES -- most
  # commonly DoubletFinder's "pANN_<pN>_<pK>_<nExp>", whose pK and nExp are
  # fit per sample, so every object in a list gets a differently-named
  # column even though the column counts all match. Plain rbind() rejects
  # that with match.names()'s "names do not match previous names", which
  # points at neither the offending column nor the sample it came from.
  # Dropped columns are accumulated here and reported in ONE warning at the
  # end, rather than warned about inside each call: `pre` and `post` are
  # collected separately but almost always carry the same object-specific
  # columns, so warning per call means the user reads the same list twice.
  dropped_cols <- character(0)

  .rbind_common <- function(parts, what) {
    if (!length(parts)) return(NULL)
    all_cols <- lapply(parts, colnames)
    common   <- Reduce(intersect, all_cols)
    if (!length(common)) {
      stop("The ", what, " share no metadata columns in common.")
    }
    dropped_cols <<- union(dropped_cols,
                           setdiff(Reduce(union, all_cols), common))
    do.call(rbind, lapply(parts, function(p) p[, common, drop = FALSE]))
  }

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
        d <- .rbind_common(parts, paste0("'", sample_col, "' groups"))
      } else {
        d <- .one(x, "all")
      }
    } else if (is.list(x)) {
      if (is.null(names(x))) stop("List inputs must be named.")
      d <- .rbind_common(lapply(names(x), function(nm) .one(x[[nm]], nm)),
                         "object in the list")
    } else {
      stop("`pre`/`post` must be a Seurat object or named list.")
    }
    d$state <- state
    d
  }

  d1 <- .collect(pre,  "pre")
  d2 <- .collect(post, "post")
  if (length(dropped_cols)) {
    warning("Dropping ", length(dropped_cols), " metadata column(s) not ",
            "present in every object: ",
            paste(utils::head(dropped_cols, 6L), collapse = ", "),
            if (length(dropped_cols) > 6L) ", ..." else "",
            call. = FALSE)
  }
  # `pre`/`post` are documented as needing the same sample set -- if they
  # don't (e.g. a sample dropped entirely, or one added between calls), the
  # per-sample pct_kept annotation in .qc_counts() below silently divides by
  # a 0-vs-0 (pre = 0) or produces a nonsensical >>100% figure rather than
  # erroring, so catch the mismatch here instead.
  mismatched <- union(setdiff(unique(d1$sample), unique(d2$sample)),
                      setdiff(unique(d2$sample), unique(d1$sample)))
  if (length(mismatched)) {
    stop("`pre` and `post` have mismatched sample names -- present in only ",
         "one of the two: ", paste(mismatched, collapse = ", "),
         ". They must cover the same samples.")
  }
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
  # Stacked rather than "n_before -> n_after (pct%)" on one line: at 8+
  # samples the single-line form is roughly twice the width of its x slot,
  # so adjacent labels run into each other. Breaking after each number makes
  # a label about as wide as its longest number. Kept ASCII deliberately --
  # a Unicode arrow here would be a non-ASCII character in package source
  # (an R CMD check warning) and is device-dependent to render.
  wide$label <- sprintf("%d\n-> %d\n(%.1f%%)",
                        wide$pre, wide$post, wide$pct_kept)
  wide
}


# ============================================================================
# Internal: violin plot for one metric, pre/post side by side, per sample
# ============================================================================
#' @keywords internal
#' @noRd
.qc_violin <- function(df, metric, log_y, colors, counts, show_counts,
                       counts_size = 2.5) {
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
    Ol_Reliable() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  annotating <- isTRUE(show_counts) && nrow(counts) > 0

  # The count annotation is text drawn in absolute (pt) units, so it occupies
  # NO height in data space -- the y scale expands to fit the violins and
  # nothing else, and the label's top line gets clipped at the panel edge.
  # (The old code tried to buy room by pushing the label's y position to
  # max * 1.05/1.1, but moving the anchor up inside a panel that grows to
  # match only pushes the text further past the top.) Two things fix it
  # properly: extra expansion at the top of the scale, in scale units so it
  # works the same on log10 and linear, and clip = "off" so anything still
  # protruding is drawn into the margin rather than cut.
  # Scale the headroom to how tall the label actually is, so this keeps
  # working if the label format gains or loses a line.
  n_lines <- if (annotating) {
    max(lengths(strsplit(as.character(counts$label), "\n", fixed = TRUE)), 1L)
  } else {
    0L
  }
  y_expand <- ggplot2::expansion(mult = c(0.05, if (annotating) {
    0.06 + 0.06 * n_lines
  } else {
    0.05
  }))
  p <- p + if (metric %in% log_y) {
    ggplot2::scale_y_log10(expand = y_expand)
  } else {
    ggplot2::scale_y_continuous(expand = y_expand)
  }

  if (annotating) {
    y_max <- suppressWarnings(max(df$.y, na.rm = TRUE))
    if (is.finite(y_max)) {
      p <- p +
        ggplot2::annotate(
          "text",
          x          = counts$sample,
          y          = y_max,
          label      = counts$label,
          size       = counts_size,
          vjust      = 0,
          lineheight = 0.9
        ) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 5.5,
                                                     b = 5.5, l = 5.5))
    }
  }
  p
}
