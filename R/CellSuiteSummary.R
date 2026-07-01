#' One-command project summary of a Seurat object
#'
#' Returns a compact list of everything you'd want at the top of a README
#' or project handoff for a Seurat object: cell / gene counts, per-cluster
#' counts, QC-metric summaries, top markers per cluster, reductions
#' present, and (if a sample column is supplied) per-sample cell counts.
#' Printed nicely by \code{print.CellSuiteSummary}.
#'
#' @param obj A Seurat object.
#' @param cluster_col Metadata column holding cluster / cell-type ids.
#'   Default \code{"seurat_clusters"}.
#' @param sample_col Optional metadata column for the per-sample
#'   breakdown. Default \code{NULL}.
#' @param qc_metrics Character vector of QC metric metadata columns to
#'   summarize (median + IQR). \code{NULL} (default) auto-detects the
#'   standard set (\code{nCount_RNA}, \code{nFeature_RNA}, \code{percent.mt},
#'   ...).
#' @param top_markers Number of top marker genes per cluster to include.
#'   \code{0} skips marker computation. Default 0 (skip; markers are
#'   expensive).
#' @param verbose Print progress. Default TRUE.
#' @return An object of class \code{"CellSuiteSummary"} (a list) with
#'   pretty print method.
#' @examples
#' \dontrun{
#' s <- CellSuiteSummary(obj, sample_col = "orig.ident", top_markers = 5)
#' print(s)
#' # Access individual pieces
#' s$cluster_counts
#' s$qc_summary
#' }
#' @importFrom Seurat DefaultAssay Assays Reductions Idents FindAllMarkers
#' @export
CellSuiteSummary <- function(obj,
                             cluster_col  = "seurat_clusters",
                             sample_col   = NULL,
                             qc_metrics   = NULL,
                             top_markers  = 0,
                             verbose      = TRUE) {

  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")

  # Auto-detect QC metrics if not specified
  standard <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb",
                "percent.hb", "nCount_Spatial", "nFeature_Spatial")
  if (is.null(qc_metrics)) {
    qc_metrics <- intersect(standard, colnames(obj@meta.data))
  }

  # Cluster counts
  cluster_counts <- NULL
  if (cluster_col %in% colnames(obj@meta.data)) {
    tab <- table(as.character(obj@meta.data[[cluster_col]]))
    cluster_counts <- data.frame(
      cluster = names(tab),
      n_cells = as.integer(tab),
      pct     = round(100 * as.numeric(tab) / sum(tab), 1),
      stringsAsFactors = FALSE
    )
    cluster_counts <- cluster_counts[order(-cluster_counts$n_cells), ]
    rownames(cluster_counts) <- NULL
  }

  # Sample counts
  sample_counts <- NULL
  if (!is.null(sample_col) && sample_col %in% colnames(obj@meta.data)) {
    tab <- table(as.character(obj@meta.data[[sample_col]]))
    sample_counts <- data.frame(
      sample  = names(tab),
      n_cells = as.integer(tab),
      pct     = round(100 * as.numeric(tab) / sum(tab), 1),
      stringsAsFactors = FALSE
    )
    sample_counts <- sample_counts[order(-sample_counts$n_cells), ]
    rownames(sample_counts) <- NULL
  }

  # QC summary
  qc_summary <- NULL
  if (length(qc_metrics)) {
    qc_summary <- do.call(rbind, lapply(qc_metrics, function(m) {
      v <- obj@meta.data[[m]]
      data.frame(
        metric = m,
        median = round(stats::median(v, na.rm = TRUE), 2),
        q25    = round(stats::quantile(v, 0.25, na.rm = TRUE), 2),
        q75    = round(stats::quantile(v, 0.75, na.rm = TRUE), 2),
        min    = round(min(v, na.rm = TRUE), 2),
        max    = round(max(v, na.rm = TRUE), 2),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }))
  }

  # Top markers per cluster (optional; expensive)
  markers <- NULL
  if (top_markers > 0 && cluster_col %in% colnames(obj@meta.data)) {
    if (isTRUE(verbose)) {
      message("--- Running FindAllMarkers (may take a while) ---")
    }
    Seurat::Idents(obj) <- as.factor(as.character(obj@meta.data[[cluster_col]]))
    all_m <- Seurat::FindAllMarkers(obj, only.pos = TRUE,
                                    min.pct = 0.25,
                                    logfc.threshold = 0.25,
                                    verbose = FALSE)
    all_m <- all_m[all_m$p_val_adj < 0.05, ]
    all_m <- do.call(rbind, lapply(split(all_m, all_m$cluster), function(d) {
      d[order(-d$avg_log2FC), ][seq_len(min(top_markers, nrow(d))), ]
    }))
    markers <- all_m[, c("cluster", "gene", "avg_log2FC", "p_val_adj")]
    rownames(markers) <- NULL
  }

  out <- list(
    n_cells        = ncol(obj),
    n_genes        = tryCatch(nrow(obj[[Seurat::DefaultAssay(obj)]]),
                              error = function(e) NA_integer_),
    default_assay  = Seurat::DefaultAssay(obj),
    assays         = Seurat::Assays(obj),
    reductions     = Seurat::Reductions(obj),
    images         = names(obj@images),
    cluster_col    = cluster_col,
    n_clusters     = if (!is.null(cluster_counts)) nrow(cluster_counts) else NA,
    cluster_counts = cluster_counts,
    sample_col     = sample_col,
    sample_counts  = sample_counts,
    qc_summary     = qc_summary,
    top_markers    = markers
  )
  class(out) <- c("CellSuiteSummary", "list")
  out
}


#' @export
print.CellSuiteSummary <- function(x, ...) {
  cat("Seurat object summary\n")
  cat(rep("-", 45), "\n", sep = "")
  cat(sprintf("  cells:         %d\n", x$n_cells))
  cat(sprintf("  genes:         %d\n", x$n_genes))
  cat(sprintf("  default assay: %s\n", x$default_assay))
  if (length(x$assays) > 1L) {
    cat(sprintf("  all assays:    %s\n", paste(x$assays, collapse = ", ")))
  }
  if (length(x$reductions)) {
    cat(sprintf("  reductions:    %s\n",
                paste(x$reductions, collapse = ", ")))
  }
  if (length(x$images)) {
    cat(sprintf("  images:        %s\n",
                paste(x$images, collapse = ", ")))
  }
  cat(sprintf("  clusters:      %s (col: %s)\n",
              x$n_clusters, x$cluster_col))

  if (!is.null(x$cluster_counts)) {
    cat("\nCluster sizes\n")
    print(utils::head(x$cluster_counts, 12), row.names = FALSE)
    if (nrow(x$cluster_counts) > 12L) {
      cat(sprintf("... (%d more clusters)\n",
                  nrow(x$cluster_counts) - 12L))
    }
  }
  if (!is.null(x$sample_counts)) {
    cat("\nSamples\n")
    print(x$sample_counts, row.names = FALSE)
  }
  if (!is.null(x$qc_summary)) {
    cat("\nQC metric summaries\n")
    print(x$qc_summary, row.names = FALSE)
  }
  if (!is.null(x$top_markers)) {
    cat("\nTop markers per cluster\n")
    print(x$top_markers, row.names = FALSE)
  }
  invisible(x)
}
