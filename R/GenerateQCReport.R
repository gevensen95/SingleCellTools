#' Generate an HTML QC report for a Seurat object or list
#'
#' Builds an \code{rmarkdown} report with a comprehensive set of per-sample
#' single-cell QC plots and diagnostics. Sections that don't apply to the
#' input (cell-cycle phase, spatial maps, edge-spot summary, sample
#' correlation) are silently skipped.
#'
#' Sample identification:
#' \itemize{
#'   \item If \code{obj} is a single Seurat object, the per-cell
#'     \code{sample_col} (default \code{"orig.ident"}) is used as the
#'     sample identifier (one row in the summary per unique value). If
#'     \code{sample_col} is missing from the object, every cell is
#'     treated as the single \code{"sample"}.
#'   \item If \code{obj} is a list of Seurat objects, each list element is
#'     treated as one sample, using the list names (or \code{obj_1},
#'     \code{obj_2}, ... if unnamed). \code{sample_col} is ignored.
#' }
#'
#' Sections included:
#' \enumerate{
#'   \item Parameters appendix (top of report) — records the arguments.
#'   \item Overview table — cells, genes, median QC metrics per sample.
#'   \item Cell-count bar plot — quickest possible visual outlier check.
#'   \item Violin distributions of \code{metadata_cols} per sample.
#'   \item Density overlay — same metrics, all samples overlaid in one panel.
#'   \item QC scatters: \code{nCount_RNA} vs \code{nFeature_RNA}, and
#'     \code{percent.mt} vs \code{nFeature_RNA}.
#'   \item Top \code{top_n_genes} expressed genes per sample.
#'   \item Suggested filtering cutoffs (median ± \code{mad_multiplier} * MAD)
#'     and the cells that would survive each cutoff.
#'   \item Cell-cycle phase breakdown (if a \code{Phase} column exists).
#'   \item Doublet calls (if \code{doublet_col} exists).
#'   \item Edge-spot summary (if an \code{is_edge} column exists).
#'   \item Sample-sample pseudobulk correlation heatmap (>= 2 samples).
#'   \item Spatial QC maps per FOV (if \code{obj@images} is populated).
#'   \item Session info appendix.
#' }
#'
#' @param obj A Seurat object or a (optionally named) list of them.
#' @param output_file Path to write the HTML report to. May be relative
#'   (resolved against the current working directory) or absolute.
#' @param title Title for the report.
#' @param metadata_cols Character vector of metadata columns to plot as
#'   violin distributions and use for the cutoff table.
#' @param doublet_col Metadata column holding doublet calls. \code{NULL}
#'   skips the doublet section.
#' @param top_n_genes Number of top-expressed genes per sample to display.
#'   Default 20.
#' @param mad_multiplier Multiplier for MAD when suggesting cutoffs.
#'   Default 3 (the typical "outlier" threshold).
#' @param assay Assay to read counts from for top-expressed genes and the
#'   correlation heatmap. Default DefaultAssay of each input.
#' @param sample_col Metadata column used as the sample identifier when
#'   \code{obj} is a single Seurat object. Default \code{"orig.ident"}.
#'   Ignored when \code{obj} is a list (the list names are used instead).
#' @param log_skewed Logical; if TRUE (default), automatically log10-transform
#'   metrics whose distribution has a heavy right tail (e.g.,
#'   \code{nCount_RNA}). A metric is log-transformed when its
#'   \code{max / median > log_threshold}.
#' @param log_threshold Ratio threshold for triggering log transformation.
#'   Default 10 (i.e., the metric's max is more than 10x its median).
#' @param spatial_max_cols Maximum number of spatial-QC panels per row.
#'   Default 3 — caps facet density to keep panels wide enough that the
#'   FOV/sample names aren't truncated.
#' @return Invisibly, the absolute path of the rendered report.
#' @export
GenerateQCReport <- function(obj,
                             output_file      = "qc_report.html",
                             title            = "Single-cell QC report",
                             metadata_cols    = c("nFeature_RNA",
                                                  "nCount_RNA",
                                                  "percent.mt"),
                             sample_col       = "orig.ident",
                             doublet_col      = "doublet_finder",
                             top_n_genes      = 20,
                             mad_multiplier   = 3,
                             assay            = NULL,
                             log_skewed       = TRUE,
                             log_threshold    = 10,
                             spatial_max_cols = 3) {

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required. install.packages('rmarkdown')")
  }
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("Package 'knitr' is required. install.packages('knitr')")
  }

  call_text <- paste(deparse(match.call(), width.cutoff = 80L),
                     collapse = "\n")

  # ---- Normalize input ---------------------------------------------------
  # samples: named list of lists, each entry { md, counts (sparse), images }
  # n_genes_per_sample: integer vector
  if (inherits(obj, "Seurat")) {
    src_obj <- obj
    a       <- if (is.null(assay)) Seurat::DefaultAssay(src_obj) else assay
    md_full <- src_obj@meta.data
    if (sample_col %in% colnames(md_full)) {
      split_key <- as.character(md_full[[sample_col]])
    } else {
      if (sample_col != "orig.ident") {
        warning("`sample_col` ('", sample_col, "') not found in metadata; ",
                "treating every cell as one sample.")
      }
      split_key <- rep("sample", nrow(md_full))
    }
    sample_to_cells <- split(rownames(md_full), split_key)
    full_counts <- tryCatch(
      SeuratObject::LayerData(src_obj, assay = a, layer = "counts"),
      error = function(e) NULL
    )
    samples <- lapply(names(sample_to_cells), function(nm) {
      cells <- sample_to_cells[[nm]]
      list(
        name   = nm,
        md     = md_full[cells, , drop = FALSE],
        counts = if (!is.null(full_counts)) full_counts[, cells, drop = FALSE] else NULL,
        images = src_obj@images,
        cells  = cells
      )
    })
    names(samples) <- names(sample_to_cells)
    n_genes_per_sample <- setNames(rep(nrow(src_obj), length(samples)),
                                    names(samples))
  } else if (is.list(obj) &&
             all(vapply(obj, inherits, logical(1), "Seurat"))) {
    if (is.null(names(obj))) names(obj) <- paste0("obj_", seq_along(obj))
    samples <- lapply(names(obj), function(nm) {
      o <- obj[[nm]]
      a <- if (is.null(assay)) Seurat::DefaultAssay(o) else assay
      list(
        name   = nm,
        md     = o@meta.data,
        counts = tryCatch(
          SeuratObject::LayerData(o, assay = a, layer = "counts"),
          error = function(e) NULL
        ),
        images = o@images,
        cells  = colnames(o)
      )
    })
    names(samples) <- names(obj)
    # Seurat v5 note: nrow() on some assay classes returns a `double`, not
    # `integer`, so a strict integer(1) template fails with
    #   "values must be type 'integer', but FUN(X[[1]]) result is type 'double'"
    # Use numeric(1) which accepts both.
    n_genes_per_sample <- vapply(obj, nrow, numeric(1))
    names(n_genes_per_sample) <- names(obj)
  } else {
    stop("`obj` must be a Seurat object or a list of Seurat objects.")
  }

  # ---- Resolve metadata_cols to actual column names ---------------------
  # Substitute hardcoded nFeature_RNA / nCount_RNA with the first column
  # that matches the pattern in the data, so the report works for Spatial,
  # Xenium, ATAC, etc. without requiring the user to pass metadata_cols.
  all_cols <- unique(unlist(lapply(samples, function(s) colnames(s$md))))
  resolve_one <- function(req) {
    if (req %in% all_cols) return(req)
    if (startsWith(req, "nFeature")) {
      m <- grep("^nFeature", all_cols, value = TRUE)
      return(if (length(m)) m[1] else NA_character_)
    }
    if (startsWith(req, "nCount")) {
      m <- grep("^nCount", all_cols, value = TRUE)
      return(if (length(m)) m[1] else NA_character_)
    }
    NA_character_
  }
  resolved <- vapply(metadata_cols, resolve_one, character(1))
  resolved <- unique(na.omit(resolved))
  if (!identical(resolved, metadata_cols)) {
    message("Resolved metadata_cols to: ",
            paste(resolved, collapse = ", "))
  }
  metadata_cols <- resolved

  # ---- Long-format QC frame for violin/density --------------------------
  message("--- Collecting metadata across ", length(samples), " sample(s) ---")
  long_qc <- do.call(rbind, lapply(samples, function(s) {
    cols <- intersect(metadata_cols, colnames(s$md))
    if (!length(cols)) return(NULL)
    do.call(rbind, lapply(cols, function(c) {
      data.frame(sample = s$name, metric = c, value = s$md[[c]],
                 stringsAsFactors = FALSE, row.names = NULL)
    }))
  }))
  if (!is.null(long_qc)) rownames(long_qc) <- NULL

  # ---- Optionally log10-transform heavy-tailed metrics ------------------
  # Count distributions (nCount_*) routinely span 2-3 orders of magnitude;
  # plotting them on a linear axis squeezes 90 percent of cells into a thin
  # strip at the bottom. Auto-transform any metric whose max/median exceeds
  # `log_threshold`. The metric name is suffixed with " (log10)" so the
  # axis labels reflect what's actually shown.
  if (isTRUE(log_skewed) && !is.null(long_qc)) {
    metric_skew <- vapply(unique(long_qc$metric), function(m) {
      v <- long_qc$value[long_qc$metric == m]
      v <- v[is.finite(v) & v > 0]
      if (length(v) < 10) return(FALSE)
      med <- stats::median(v)
      if (med <= 0) return(FALSE)
      (max(v) / med) > log_threshold
    }, logical(1))
    metrics_to_log <- names(metric_skew)[metric_skew]
    if (length(metrics_to_log)) {
      message("Log10-transforming heavy-tailed metric(s): ",
              paste(metrics_to_log, collapse = ", "))
      for (m in metrics_to_log) {
        idx <- long_qc$metric == m & long_qc$value > 0
        long_qc$value[idx] <- log10(long_qc$value[idx])
        long_qc$metric[long_qc$metric == m] <- paste0(m, " (log10)")
      }
    }
  }

  # ---- Summary table -----------------------------------------------------
  summary_df <- data.frame(
    sample  = names(samples),
    n_cells = vapply(samples, function(s) nrow(s$md), numeric(1)),
    n_genes = unname(n_genes_per_sample[names(samples)]),
    stringsAsFactors = FALSE
  )
  for (mc in metadata_cols) {
    summary_df[[paste0("median_", mc)]] <- vapply(samples, function(s) {
      if (mc %in% colnames(s$md)) stats::median(s$md[[mc]], na.rm = TRUE)
      else NA_real_
    }, numeric(1))
  }
  rownames(summary_df) <- NULL

  # ---- Scatter frame: nCount vs nFeature, percent.mt vs nFeature --------
  # Pick the first nFeature* and nCount* columns we can find — works for
  # RNA, Spatial, Xenium, ATAC, etc.
  scatter_df <- do.call(rbind, lapply(samples, function(s) {
    md <- s$md
    feat_col  <- grep("^nFeature", colnames(md), value = TRUE)[1]
    count_col <- grep("^nCount",   colnames(md), value = TRUE)[1]
    if (is.na(feat_col) || is.na(count_col)) return(NULL)
    data.frame(
      sample     = s$name,
      nFeature   = md[[feat_col]],
      nCount     = md[[count_col]],
      percent.mt = if ("percent.mt" %in% colnames(md)) md$percent.mt else NA_real_,
      feat_col   = feat_col,
      count_col  = count_col,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))

  # ---- Top expressed genes per sample -----------------------------------
  top_genes_df <- do.call(rbind, lapply(samples, function(s) {
    if (is.null(s$counts) || nrow(s$counts) == 0) return(NULL)
    means <- Matrix::rowMeans(s$counts)
    top_idx <- order(means, decreasing = TRUE)[seq_len(min(top_n_genes, length(means)))]
    data.frame(
      sample    = s$name,
      gene      = rownames(s$counts)[top_idx],
      mean_expr = as.numeric(means[top_idx]),
      rank      = seq_along(top_idx),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))

  # ---- Suggested filter cutoffs (median +/- mad_multiplier * MAD) -------
  cutoffs_df <- do.call(rbind, lapply(samples, function(s) {
    do.call(rbind, lapply(intersect(metadata_cols, colnames(s$md)), function(mc) {
      v   <- s$md[[mc]]
      med <- stats::median(v, na.rm = TRUE)
      m   <- stats::mad(v, na.rm = TRUE)
      lo  <- max(0, med - mad_multiplier * m)
      hi  <- med + mad_multiplier * m
      n_total <- sum(!is.na(v))
      n_pass  <- sum(v >= lo & v <= hi, na.rm = TRUE)
      data.frame(
        sample      = s$name,
        metric      = mc,
        median      = round(med, 2),
        mad         = round(m,   2),
        suggest_lo  = round(lo,  2),
        suggest_hi  = round(hi,  2),
        n_total     = n_total,
        n_pass      = n_pass,
        pct_pass    = round(100 * n_pass / max(1, n_total), 1),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }))
  }))

  # ---- Cell cycle (if Phase column exists in any sample) -----------------
  cc_df <- do.call(rbind, lapply(samples, function(s) {
    if (!"Phase" %in% colnames(s$md)) return(NULL)
    tab <- table(s$md$Phase, useNA = "ifany")
    data.frame(
      sample = s$name,
      phase  = names(tab),
      n      = as.integer(tab),
      prop   = as.numeric(tab) / max(1, sum(tab)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
  if (!is.null(cc_df) && nrow(cc_df) == 0) cc_df <- NULL

  # ---- Doublet summary ---------------------------------------------------
  doublet_summary <- NULL
  if (!is.null(doublet_col)) {
    has_d <- vapply(samples, function(s) doublet_col %in% colnames(s$md),
                    logical(1))
    if (any(has_d)) {
      doublet_summary <- do.call(rbind, lapply(samples[has_d], function(s) {
        tab <- table(s$md[[doublet_col]])
        data.frame(
          sample      = s$name,
          singlet     = as.integer(tab["Singlet"]),
          doublet     = as.integer(tab["Doublet"]),
          pct_doublet = round(100 * as.integer(tab["Doublet"]) / sum(tab), 2),
          stringsAsFactors = FALSE
        )
      }))
      rownames(doublet_summary) <- NULL
    }
  }

  # ---- Edge-spot summary -------------------------------------------------
  edge_summary <- do.call(rbind, lapply(samples, function(s) {
    if (!"is_edge" %in% colnames(s$md)) return(NULL)
    n_total <- nrow(s$md)
    n_edge  <- sum(isTRUE(s$md$is_edge) | s$md$is_edge == TRUE, na.rm = TRUE)
    data.frame(
      sample   = s$name,
      n_total  = n_total,
      n_edge   = n_edge,
      pct_edge = round(100 * n_edge / max(1, n_total), 2),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
  if (!is.null(edge_summary) && nrow(edge_summary) == 0) edge_summary <- NULL

  # ---- Pseudobulk sample-sample correlation -----------------------------
  cor_mat <- NULL
  if (length(samples) >= 2) {
    have_counts <- vapply(samples, function(s) !is.null(s$counts), logical(1))
    if (all(have_counts)) {
      gene_sets    <- lapply(samples, function(s) rownames(s$counts))
      common_genes <- Reduce(intersect, gene_sets)
      if (length(common_genes) >= 100) {
        pb <- vapply(samples, function(s) {
          Matrix::rowSums(s$counts[common_genes, , drop = FALSE])
        }, FUN.VALUE = numeric(length(common_genes)))
        colnames(pb) <- names(samples)
        # CPM-like normalization, then log1p
        size_factors <- colSums(pb)
        size_factors[size_factors == 0] <- 1
        pb_norm <- log1p(t(t(pb) / (size_factors / 1e4)))
        cor_mat <- stats::cor(pb_norm, method = "pearson")
      }
    }
  }

  # ---- Spatial QC --------------------------------------------------------
  # Detect whether *any* sample has any spatial images attached. Avoids the
  # `make.unique(NULL)` crash that fires when no sample is spatial.
  spatial_df <- NULL
  has_images <- any(vapply(samples,
                           function(s) length(s$images) > 0,
                           logical(1)))
  if (has_images) {
    spatial_df <- do.call(rbind, lapply(samples, function(s) {
      do.call(rbind, lapply(names(s$images), function(img_name) {
        coords <- tryCatch(
          Seurat::GetTissueCoordinates(s$images[[img_name]], which = "centroids"),
          error = function(e) NULL
        )
        if (is.null(coords)) return(NULL)
        coords <- as.data.frame(coords)
        if (all(c("imagecol", "imagerow") %in% colnames(coords))) {
          coords$x <- coords$imagecol; coords$y <- coords$imagerow
        }
        if (!all(c("x", "y") %in% colnames(coords))) return(NULL)
        if (!"cell" %in% colnames(coords)) coords$cell <- rownames(coords)
        # Only keep cells that belong to this sample
        keep <- coords$cell %in% s$cells
        coords <- coords[keep, c("cell", "x", "y"), drop = FALSE]
        if (!nrow(coords)) return(NULL)
        # Attach QC metrics
        for (mc in intersect(metadata_cols, colnames(s$md))) {
          coords[[mc]] <- s$md[coords$cell, mc]
        }
        coords$sample <- s$name
        coords$image  <- img_name
        coords
      }))
    }))
    if (!is.null(spatial_df) && nrow(spatial_df) == 0) spatial_df <- NULL
  }

  # ---- Compute dynamic figure heights for facetted panels ---------------
  # Each top-genes panel needs enough vertical room to render `top_n_genes`
  # gene labels along the y-axis after coord_flip. Cap panels per row at 3
  # so each panel stays wide enough to read.
  if (!is.null(top_genes_df)) {
    n_smp        <- length(unique(top_genes_df$sample))
    top_ncol     <- min(3, n_smp)
    top_nrow     <- ceiling(n_smp / top_ncol)
    # ~0.18 in per gene label + 0.8 in for axis/title/padding
    top_genes_fig_h <- max(4, top_nrow * (top_n_genes * 0.18 + 0.8))
  } else {
    top_ncol        <- 1
    top_genes_fig_h <- 4
  }

  # Spatial layout: cap columns per row and scale strip text so FOV names
  # don't get truncated when there are many panels. With more panels per
  # row, each strip is narrower, so we shrink the font.
  #
  # Figure height needs to match what `theme(aspect.ratio = 1)` will
  # actually consume — square panels of width fig.width/ncol stacked
  # `nrow` deep, plus a small buffer for the title / legend / strips.
  # Over-allocating fig.height shows up as huge whitespace before/after
  # the figure in the rendered HTML.
  spatial_ncol        <- 1L
  spatial_fig_h       <- 6
  spatial_strip_size  <- 11
  fig_w               <- 9   # matches the chunk-level fig.width default
  if (!is.null(spatial_df)) {
    facet_var <- if (length(unique(spatial_df$image)) > 1) "image" else "sample"
    n_facets  <- length(unique(spatial_df[[facet_var]]))
    spatial_ncol       <- min(spatial_max_cols, max(1L, n_facets))
    spatial_nrow       <- ceiling(n_facets / spatial_ncol)
    # Effective per-panel side after axis/legend overhead (~0.5 in).
    panel_in <- max(2, (fig_w - 0.5) / spatial_ncol)
    # +1.5 in for title + subtitle + legend + strip strip.
    spatial_fig_h      <- panel_in * spatial_nrow + 1.5
    # Font size: 1 col -> 14pt; each extra col drops 2pt; floor at 7pt.
    spatial_strip_size <- max(7, 16 - 2 * spatial_ncol)
  }

  # ---- Stash inputs for the Rmd ------------------------------------------
  data_path <- tempfile(fileext = ".rds")
  saveRDS(list(
    title            = title,
    call_text        = call_text,
    metadata_cols    = metadata_cols,
    long_qc          = long_qc,
    summary_df       = summary_df,
    scatter_df       = scatter_df,
    top_genes_df       = top_genes_df,
    top_ncol           = top_ncol,
    cutoffs_df         = cutoffs_df,
    cc_df              = cc_df,
    doublet_summary    = doublet_summary,
    edge_summary       = edge_summary,
    cor_mat            = cor_mat,
    spatial_df         = spatial_df,
    spatial_ncol       = spatial_ncol,
    spatial_strip_size = spatial_strip_size
  ), data_path)

  # ---- Build Rmd ---------------------------------------------------------
  rmd <- c(
    "---",
    sprintf("title: \"%s\"", title),
    "output: html_document",
    "params:",
    sprintf("  data_path: \"%s\"", data_path),
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE,",
    "                      fig.width = 9, fig.height = 5)",
    "b <- readRDS(params$data_path)",
    "```",
    "",
    "## Run details",
    "",
    "```{r call}",
    "cat('```\\n', b$call_text, '\\n```', sep = '')",
    "```",
    "",
    "## Overview",
    "",
    "```{r overview}",
    "knitr::kable(b$summary_df, row.names = FALSE)",
    "```",
    "",
    "## Cell counts",
    "",
    "```{r counts_bar, fig.height=3.5}",
    "ggplot2::ggplot(b$summary_df, ggplot2::aes(x = sample, y = n_cells, fill = sample)) +",
    "  ggplot2::geom_col() +",
    "  ggplot2::geom_text(ggplot2::aes(label = n_cells), vjust = -0.3, size = 3) +",
    "  ggplot2::theme_bw() +",
    "  ggplot2::theme(legend.position = 'none',",
    "                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +",
    "  ggplot2::labs(x = NULL, y = 'Number of cells')",
    "```",
    "",
    "## QC distributions",
    "",
    "```{r vln, fig.height=5}",
    "if (is.null(b$long_qc)) {",
    "  cat('No requested metadata columns were found in any object.')",
    "} else {",
    "  ggplot2::ggplot(b$long_qc, ggplot2::aes(x = sample, y = value, fill = sample)) +",
    "    ggplot2::geom_violin(scale = 'width', trim = TRUE, color = NA) +",
    "    ggplot2::geom_boxplot(width = 0.1, outlier.size = 0.2, fill = 'white') +",
    "    ggplot2::facet_wrap(~ metric, scales = 'free_y', ncol = 1) +",
    "    ggplot2::theme_bw() +",
    "    ggplot2::theme(legend.position = 'none',",
    "                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +",
    "    ggplot2::labs(x = NULL, y = NULL)",
    "}",
    "```",
    "",
    "## QC density overlay",
    "",
    "```{r density, fig.height=4}",
    "if (!is.null(b$long_qc)) {",
    "  ggplot2::ggplot(b$long_qc, ggplot2::aes(x = value, color = sample)) +",
    "    ggplot2::geom_density(linewidth = 0.6) +",
    "    ggplot2::facet_wrap(~ metric, scales = 'free', ncol = 1) +",
    "    ggplot2::theme_bw() +",
    "    ggplot2::labs(x = NULL, y = 'Density')",
    "}",
    "```",
    "",
    "## QC scatter (nCount vs nFeature)",
    "",
    "```{r scatter1, fig.height=5}",
    "if (!is.null(b$scatter_df)) {",
    "  feat_lab  <- paste0(unique(b$scatter_df$feat_col),  collapse = ' / ')",
    "  count_lab <- paste0(unique(b$scatter_df$count_col), collapse = ' / ')",
    "  ggplot2::ggplot(b$scatter_df, ggplot2::aes(x = nCount, y = nFeature)) +",
    "    ggplot2::geom_point(size = 0.2, alpha = 0.3) +",
    "    ggplot2::facet_wrap(~ sample) +",
    "    ggplot2::scale_x_log10() +",
    "    ggplot2::scale_y_log10() +",
    "    ggplot2::theme_bw() +",
    "    ggplot2::labs(x = paste0(count_lab, ' (log10)'),",
    "                  y = paste0(feat_lab,  ' (log10)'))",
    "} else cat('No nCount* / nFeature* columns found.')",
    "```",
    "",
    "## QC scatter (percent.mt vs nFeature)",
    "",
    "```{r scatter2, fig.height=5}",
    "if (!is.null(b$scatter_df) && any(!is.na(b$scatter_df$percent.mt))) {",
    "  feat_lab <- paste0(unique(b$scatter_df$feat_col), collapse = ' / ')",
    "  ggplot2::ggplot(b$scatter_df, ggplot2::aes(x = nFeature, y = percent.mt)) +",
    "    ggplot2::geom_point(size = 0.2, alpha = 0.3) +",
    "    ggplot2::facet_wrap(~ sample) +",
    "    ggplot2::scale_x_log10() +",
    "    ggplot2::theme_bw() +",
    "    ggplot2::labs(x = paste0(feat_lab, ' (log10)'), y = 'percent.mt')",
    "} else cat('percent.mt not present.')",
    "```",
    "",
    "## Top expressed genes per sample",
    "",
    sprintf("```{r top_genes, fig.height=%.1f}", top_genes_fig_h),
    "if (!is.null(b$top_genes_df)) {",
    "  ggplot2::ggplot(b$top_genes_df,",
    "                  ggplot2::aes(x = stats::reorder(gene, mean_expr),",
    "                               y = mean_expr, fill = sample)) +",
    "    ggplot2::geom_col() +",
    "    ggplot2::coord_flip() +",
    "    ggplot2::facet_wrap(~ sample, scales = 'free', ncol = b$top_ncol) +",
    "    ggplot2::theme_bw() +",
    "    ggplot2::theme(legend.position = 'none',",
    "                   axis.text.y = ggplot2::element_text(size = 7)) +",
    "    ggplot2::labs(x = NULL, y = 'Mean expression')",
    "} else cat('No counts available; skipping top-gene plot.')",
    "```",
    "",
    "## Suggested filtering cutoffs",
    "",
    "```{r cutoffs}",
    "if (!is.null(b$cutoffs_df)) {",
    "  knitr::kable(b$cutoffs_df, row.names = FALSE)",
    "}",
    "```",
    "",
    "## Cell cycle",
    "",
    "```{r ccphase, fig.height=3.5}",
    "if (is.null(b$cc_df)) {",
    "  cat('No \"Phase\" column found in any sample.')",
    "} else {",
    "  ggplot2::ggplot(b$cc_df, ggplot2::aes(x = sample, y = prop, fill = phase)) +",
    "    ggplot2::geom_col() +",
    "    ggplot2::theme_bw() +",
    "    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +",
    "    ggplot2::labs(x = NULL, y = 'Proportion', fill = 'Phase')",
    "}",
    "```",
    "",
    "## Doublet calls",
    "",
    "```{r doublets}",
    "if (is.null(b$doublet_summary)) {",
    "  cat('No doublet_finder column present in any sample.')",
    "} else {",
    "  knitr::kable(b$doublet_summary, row.names = FALSE)",
    "}",
    "```",
    "",
    "## Edge spots",
    "",
    "```{r edges}",
    "if (is.null(b$edge_summary)) {",
    "  cat('No is_edge column present in any sample (skip if non-Visium).')",
    "} else {",
    "  knitr::kable(b$edge_summary, row.names = FALSE)",
    "}",
    "```",
    "",
    "## Sample-sample correlation (pseudobulk)",
    "",
    "```{r cor, fig.height=5}",
    "if (is.null(b$cor_mat)) {",
    "  cat('Not enough samples (or shared genes) for a correlation heatmap.')",
    "} else {",
    "  cm <- b$cor_mat",
    "  long <- data.frame(",
    "    sample_a = rep(rownames(cm), times = ncol(cm)),",
    "    sample_b = rep(colnames(cm), each  = nrow(cm)),",
    "    value    = as.vector(cm),",
    "    stringsAsFactors = FALSE",
    "  )",
    "  ggplot2::ggplot(long, ggplot2::aes(x = sample_a, y = sample_b, fill = value)) +",
    "    ggplot2::geom_tile() +",
    "    ggplot2::geom_text(ggplot2::aes(label = sprintf('%.2f', value)), size = 3) +",
    "    ggplot2::scale_fill_gradient2(midpoint = mean(long$value, na.rm = TRUE),",
    "                                  low = 'blue', mid = 'white', high = 'red') +",
    "    ggplot2::theme_bw() +",
    "    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +",
    "    ggplot2::labs(x = NULL, y = NULL, fill = 'Pearson r')",
    "}",
    "```",
    "",
    "## Spatial QC",
    "",
    sprintf("```{r spatial, fig.height=%.1f}", spatial_fig_h),
    "if (is.null(b$spatial_df)) {",
    "  cat('No spatial images found (skip if non-spatial data).')",
    "} else {",
    "  metric_for_color <- intersect(b$metadata_cols, colnames(b$spatial_df))[1]",
    "  if (is.na(metric_for_color)) {",
    "    cat('No metric available to color spatial plot.')",
    "  } else {",
    "    facet_var <- if (length(unique(b$spatial_df$image)) > 1) 'image' else 'sample'",
    "    # NB: dropped coord_equal() because newer ggplot2 refuses to combine",
    "    # fixed aspect ratio with free scales. Within-panel x/y proportions",
    "    # are no longer guaranteed 1:1, but each FOV fills its panel.",
    "    ggplot2::ggplot(b$spatial_df,",
    "                    ggplot2::aes(x = x, y = y, color = .data[[metric_for_color]])) +",
    "      ggplot2::geom_point(size = 0.3) +",
    "      ggplot2::facet_wrap(stats::as.formula(paste('~', facet_var)),",
    "                          scales = 'free', ncol = b$spatial_ncol) +",
    "      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.02)) +",
    "      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.02)) +",
    "      ggplot2::scale_color_viridis_c() +",
    "      ggplot2::theme_void() +",
    "      ggplot2::theme(aspect.ratio = 1,",
    "                     strip.text = ggplot2::element_text(size = b$spatial_strip_size)) +",
    "      ggplot2::labs(color = metric_for_color)",
    "  }",
    "}",
    "```",
    "",
    "## Session info",
    "",
    "```{r session}",
    "utils::sessionInfo()",
    "```"
  )

  rmd_path <- tempfile(fileext = ".Rmd")
  writeLines(rmd, rmd_path)

  # ---- Render -------------------------------------------------------------
  output_file <- normalizePath(output_file, mustWork = FALSE)
  out_dir  <- dirname(output_file)
  out_name <- basename(output_file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  message(sprintf("--- Rendering report to %s ---", output_file))
  rmarkdown::render(
    input       = rmd_path,
    output_file = out_name,
    output_dir  = out_dir,
    quiet       = TRUE,
    envir       = new.env()
  )

  if (!file.exists(output_file)) {
    warning("rmarkdown finished but '", output_file, "' is missing.")
  } else {
    message(sprintf("--- Wrote %s (%.1f KB) ---", output_file,
                    file.info(output_file)$size / 1024))
  }

  # ---- Sidecar: machine-readable cutoffs for ApplyQCFilters() -----------
  # Write the suggested-cutoffs table next to the HTML as CSV so it can be
  # loaded programmatically (e.g. by ApplyQCFilters) without parsing the
  # rendered report.
  if (!is.null(cutoffs_df) && nrow(cutoffs_df) > 0) {
    cutoffs_csv <- sub("\\.html?$", "", output_file, ignore.case = TRUE)
    cutoffs_csv <- paste0(cutoffs_csv, "_cutoffs.csv")
    utils::write.csv(cutoffs_df, cutoffs_csv, row.names = FALSE)
    message(sprintf("--- Wrote sidecar cutoffs to %s ---", cutoffs_csv))
  }

  invisible(output_file)
}
