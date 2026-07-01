#' Apply QC filter cutoffs from a GenerateQCReport sidecar
#'
#' Reads the recommended per-sample, per-metric filter cutoffs that
#' \code{\link{GenerateQCReport}} writes as a sidecar CSV
#' (\code{<report>_cutoffs.csv}) and applies them to a Seurat object or a
#' named list of Seurat objects. For each sample, cells are kept only if
#' \emph{every} selected metric's value falls within
#' \code{[suggest_lo, suggest_hi]}.
#'
#' The sidecar CSV is expected to have columns:
#' \code{sample}, \code{metric}, \code{suggest_lo}, \code{suggest_hi}.
#' Any other columns (\code{median}, \code{mad}, \code{pct_pass}, ...)
#' are ignored during filtering but pass through in \code{$cutoffs} when
#' \code{return_report = TRUE}.
#'
#' @param obj A Seurat object, or a named list of Seurat objects. When a
#'   list, the names are matched against the \code{sample} column of the
#'   cutoffs table.
#' @param cutoffs Either a path to a cutoffs CSV (typically
#'   \code{<report>_cutoffs.csv} written by \code{GenerateQCReport}), or a
#'   \code{data.frame} with the columns described above, or a path to the
#'   HTML report itself (the function will look for a sibling
#'   \code{_cutoffs.csv} file).
#' @param metrics Character vector of metric names to apply. \code{NULL}
#'   (default) uses every metric present in the cutoffs table. Metric names
#'   must match metadata columns on each object. The " (log10)" suffix that
#'   GenerateQCReport appends to skewed metrics is stripped automatically
#'   for matching.
#' @param sample_col When \code{obj} is a single Seurat object that
#'   contains multiple samples in one metadata column, the name of that
#'   column (e.g. \code{"orig.ident"}). Each unique value is matched
#'   against a row of the cutoffs table and filtered independently.
#'   \code{NULL} (default) treats the whole object as one sample and
#'   requires the cutoffs table to have exactly one \code{sample} value
#'   (or the value passed via \code{single_sample_name}).
#' @param single_sample_name Only used when \code{obj} is a single Seurat
#'   object and \code{sample_col} is \code{NULL}. Which \code{sample} row
#'   of the cutoffs table to use. Defaults to the sole sample in the
#'   cutoffs table (errors if there are multiple and this is not set).
#' @param override Named list for tweaking cutoffs before applying, of the
#'   form \code{list(sample_name = list(metric_name = c(lo, hi)))}. Any
#'   entry replaces the corresponding cutoff pair; useful for tightening
#'   or loosening a single metric on a single sample without editing the
#'   CSV.
#' @param verbose Logical; print per-sample retention counts. Default TRUE.
#' @param return_report If TRUE, returns a list with \code{obj} (filtered
#'   input) and \code{report} (data frame of per-sample per-metric
#'   retention stats). If FALSE (default), returns just \code{obj}.
#' @return The filtered Seurat object or list of Seurat objects, matching
#'   the shape of the input. If \code{return_report = TRUE}, a
#'   \code{list(obj, report, cutoffs)} instead.
#'
#' @examples
#' \dontrun{
#' # Typical two-step workflow
#' GenerateQCReport(sample_list, output_file = "qc/qc_report.html")
#' sample_list_filtered <- ApplyQCFilters(
#'   sample_list,
#'   cutoffs = "qc/qc_report_cutoffs.csv"
#' )
#'
#' # Or point at the HTML and let it find the sidecar:
#' ApplyQCFilters(sample_list, cutoffs = "qc/qc_report.html")
#'
#' # Restrict to specific metrics
#' ApplyQCFilters(sample_list,
#'                cutoffs = "qc/qc_report_cutoffs.csv",
#'                metrics = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
#'
#' # Override one sample's mitochondrial cutoff
#' ApplyQCFilters(sample_list,
#'                cutoffs = "qc/qc_report_cutoffs.csv",
#'                override = list(sample1 = list(percent.mt = c(0, 15))))
#' }
#'
#' @importFrom utils read.csv
#' @export
ApplyQCFilters <- function(obj,
                           cutoffs,
                           metrics             = NULL,
                           sample_col          = NULL,
                           single_sample_name  = NULL,
                           override            = NULL,
                           verbose             = TRUE,
                           return_report       = FALSE) {

  # ---- Resolve cutoffs input ----------------------------------------------
  if (is.data.frame(cutoffs)) {
    cutoffs_df <- cutoffs
  } else if (is.character(cutoffs) && length(cutoffs) == 1L) {
    csv_path <- cutoffs
    # If the user passed the HTML report, look for the sidecar CSV next to it
    if (grepl("\\.html?$", csv_path, ignore.case = TRUE)) {
      csv_path <- sub("\\.html?$", "", csv_path, ignore.case = TRUE)
      csv_path <- paste0(csv_path, "_cutoffs.csv")
    }
    if (!file.exists(csv_path)) {
      stop("Cutoffs file not found: ", csv_path,
           ". Did GenerateQCReport() finish successfully?")
    }
    cutoffs_df <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  } else {
    stop("`cutoffs` must be a data.frame or a path to a CSV / HTML file.")
  }

  required <- c("sample", "metric", "suggest_lo", "suggest_hi")
  missing_cols <- setdiff(required, colnames(cutoffs_df))
  if (length(missing_cols) > 0) {
    stop("Cutoffs table is missing required column(s): ",
         paste(missing_cols, collapse = ", "))
  }

  # Strip the " (log10)" suffix that GenerateQCReport adds to log-transformed
  # metrics, so the metric name matches the underlying metadata column.
  cutoffs_df$metric_key <- sub(" \\(log10\\)$", "", cutoffs_df$metric)

  # ---- Optional metric restriction ----------------------------------------
  if (!is.null(metrics)) {
    keep <- cutoffs_df$metric_key %in% metrics ||
            cutoffs_df$metric     %in% metrics  # accept either form
    cutoffs_df <- cutoffs_df[cutoffs_df$metric_key %in% metrics |
                             cutoffs_df$metric     %in% metrics, ,
                             drop = FALSE]
    if (nrow(cutoffs_df) == 0) {
      stop("No cutoffs remain after `metrics` filter. Available metrics: ",
           paste(unique(cutoffs_df$metric_key), collapse = ", "))
    }
  }

  # ---- Apply user overrides -----------------------------------------------
  if (!is.null(override)) {
    for (samp in names(override)) {
      for (met in names(override[[samp]])) {
        pair <- override[[samp]][[met]]
        if (!(is.numeric(pair) && length(pair) == 2L)) {
          stop("override[[", samp, "]][[", met,
               "]] must be a length-2 numeric c(lo, hi).")
        }
        idx <- which(cutoffs_df$sample == samp &
                     (cutoffs_df$metric_key == met |
                      cutoffs_df$metric     == met))
        if (length(idx) == 0) {
          # Not present — append a fresh row
          cutoffs_df <- rbind(
            cutoffs_df,
            data.frame(sample = samp, metric = met,
                       suggest_lo = pair[1], suggest_hi = pair[2],
                       metric_key = met, stringsAsFactors = FALSE)
          )
        } else {
          cutoffs_df$suggest_lo[idx] <- pair[1]
          cutoffs_df$suggest_hi[idx] <- pair[2]
        }
      }
    }
  }

  # ---- Dispatch: list of Seurat / single Seurat with sample_col / single -
  if (is.list(obj) && !inherits(obj, "Seurat")) {
    if (is.null(names(obj))) {
      stop("List of Seurat objects must be named; names are matched ",
           "against the `sample` column of the cutoffs table.")
    }
    out_list <- lapply(names(obj), function(nm) {
      .filter_one(obj[[nm]], nm, cutoffs_df, verbose)
    })
    names(out_list) <- names(obj)
    report_df <- do.call(rbind, lapply(out_list, `[[`, "report"))
    result    <- lapply(out_list, `[[`, "obj")
  } else if (inherits(obj, "Seurat")) {
    if (!is.null(sample_col)) {
      if (!sample_col %in% colnames(obj@meta.data)) {
        stop("`sample_col` '", sample_col, "' not found in obj@meta.data.")
      }
      samples_in_obj <- unique(as.character(obj@meta.data[[sample_col]]))
      # Filter each sample's cells independently, then combine.
      keep_all <- character(0)
      report_list <- list()
      for (samp in samples_in_obj) {
        cells_samp <- rownames(obj@meta.data)[
          as.character(obj@meta.data[[sample_col]]) == samp
        ]
        sub <- subset(obj, cells = cells_samp)
        one <- .filter_one(sub, samp, cutoffs_df, verbose)
        keep_all    <- c(keep_all, colnames(one$obj))
        report_list <- c(report_list, list(one$report))
      }
      result    <- subset(obj, cells = keep_all)
      report_df <- do.call(rbind, report_list)
    } else {
      # Single sample — figure out which cutoff row to use
      if (is.null(single_sample_name)) {
        uniq_samp <- unique(cutoffs_df$sample)
        if (length(uniq_samp) == 1L) {
          single_sample_name <- uniq_samp
        } else {
          stop("Cutoffs table has ", length(uniq_samp),
               " sample(s); pass `single_sample_name` or use `sample_col`.")
        }
      }
      one <- .filter_one(obj, single_sample_name, cutoffs_df, verbose)
      result    <- one$obj
      report_df <- one$report
    }
  } else {
    stop("`obj` must be a Seurat object or a named list of Seurat objects.")
  }

  if (isTRUE(return_report)) {
    return(list(obj = result, report = report_df, cutoffs = cutoffs_df))
  }
  result
}


# ============================================================================
# Internal: apply cutoffs for one named sample to one Seurat object.
# Returns list(obj = filtered, report = per-metric retention data frame).
# ============================================================================

#' @keywords internal
#' @noRd
.filter_one <- function(so, sample_name, cutoffs_df, verbose) {
  if (!inherits(so, "Seurat")) {
    stop("Object for sample '", sample_name, "' is not a Seurat object.")
  }
  rows <- cutoffs_df[cutoffs_df$sample == sample_name, , drop = FALSE]
  n_start <- ncol(so)

  if (nrow(rows) == 0) {
    if (verbose) {
      message(sprintf(
        "  [%s] no matching sample in cutoffs table; keeping all %d cells.",
        sample_name, n_start))
    }
    return(list(
      obj    = so,
      report = data.frame(sample = sample_name, metric = NA_character_,
                          n_before = n_start, n_after = n_start,
                          pct_kept = 100,
                          stringsAsFactors = FALSE)
    ))
  }

  # Combined mask: keep a cell only if every metric passes.
  cells    <- colnames(so)
  keep     <- rep(TRUE, length(cells))
  metric_report <- list()

  for (i in seq_len(nrow(rows))) {
    mkey <- rows$metric_key[i]
    lo   <- rows$suggest_lo[i]
    hi   <- rows$suggest_hi[i]
    if (!(mkey %in% colnames(so@meta.data))) {
      if (verbose) {
        message(sprintf(
          "  [%s] metric '%s' not in metadata; skipping.",
          sample_name, mkey))
      }
      next
    }
    v <- so@meta.data[[mkey]]
    passes <- !is.na(v) & v >= lo & v <= hi
    keep <- keep & passes
    metric_report[[length(metric_report) + 1]] <- data.frame(
      sample   = sample_name,
      metric   = mkey,
      lo       = lo,
      hi       = hi,
      n_before = length(v),
      n_after  = sum(passes),
      pct_kept = round(100 * sum(passes) / length(v), 1),
      stringsAsFactors = FALSE
    )
  }

  n_end <- sum(keep)
  if (verbose) {
    message(sprintf("  [%s] %d -> %d cells (%.1f%% retained across %d metric(s))",
                    sample_name, n_start, n_end,
                    100 * n_end / max(1, n_start),
                    length(metric_report)))
  }

  filtered <- subset(so, cells = cells[keep])
  list(
    obj    = filtered,
    report = do.call(rbind, metric_report)
  )
}
