#' Pseudobulk differential expression with DESeq2
#'
#' Aggregates single-cell counts per (sample, condition) pair with
#' \code{Seurat::AggregateExpression} and runs DESeq2 to test the contrast
#' between conditions. This is the statistically correct way to compare
#' conditions across donors with single-cell data — per-cell Wilcoxon DE
#' treats each cell as independent and inflates type-I error when there
#' are multiple cells per donor.
#'
#' Two usage modes, controlled by \code{cluster_col}:
#' \describe{
#'   \item{Single-group mode (\code{cluster_col = NULL}, default)}{
#'     Runs DE on the whole object (or on the subset specified by
#'     \code{group_by} / \code{group_value}, e.g. one cell type). Returns
#'     a two-element list: \code{results} + \code{normalized_counts}.}
#'   \item{Multi-cluster mode (\code{cluster_col} supplied)}{
#'     Loops the analysis over each level of \code{cluster_col} (or the
#'     subset in \code{clusters}), returns a named list where each
#'     element is a two-element list per the single-group mode.}
#' }
#'
#' Custom \code{design} and \code{contrast} are supported so batch
#' adjustment, three-level factors, and interaction models all work. The
#' default is \code{design = ~ condition} and a two-level contrast built
#' from \code{ident_1} / \code{ident_2}.
#'
#' Gene filtering uses \code{edgeR::filterByExpr} when available, falling
#' back to a simple row-sum threshold (\code{min_counts_per_gene}).
#'
#' @param obj A Seurat object.
#' @param sample_col Metadata column identifying the biological replicate
#'   (e.g. donor, mouse). Required.
#' @param condition_col Metadata column identifying the condition. Required.
#' @param ident_1,ident_2 Two levels of \code{condition_col} to contrast.
#'   \code{ident_1} is the numerator (positive log2FC = higher in
#'   \code{ident_1}). Required if \code{contrast} is not set.
#' @param cluster_col Optional metadata column holding cluster / cell-type
#'   ids. When supplied, DE runs per level and returns a named list.
#'   \code{NULL} (default) runs a single analysis on the whole object.
#' @param clusters Character vector of \code{cluster_col} levels to run on;
#'   defaults to every level. Ignored when \code{cluster_col} is \code{NULL}.
#' @param group_by,group_value Optional pre-filter to a single group before
#'   pseudobulking (e.g. one cluster). Ignored when \code{cluster_col} is
#'   supplied.
#' @param assay Assay to aggregate. Default \code{DefaultAssay(obj)}.
#' @param design Model design as formula or string. Default
#'   \code{~ condition}.
#' @param contrast DESeq2 contrast spec. \code{NULL} builds
#'   \code{c("condition", ident_1, ident_2)}.
#' @param min_cells_per_sample Drop pseudobulk samples below this cell count.
#'   Default 10.
#' @param min_samples_per_condition Refuse the cluster / analysis if any
#'   condition has fewer than this many pseudobulk samples. Default 2.
#' @param min_counts_per_gene Fallback low-count-gene threshold when
#'   \code{edgeR} isn't available. Default 10.
#' @param verbose Print progress per cluster. Default TRUE.
#' @return In single-group mode, a list \code{list(results, normalized_counts)}.
#'   In multi-cluster mode, a named list of those lists, one per cluster.
#' @examples
#' \dontrun{
#' # Single cluster
#' res <- PseudobulkDE(obj,
#'                     sample_col    = "orig.ident",
#'                     condition_col = "treatment",
#'                     ident_1       = "drug",
#'                     ident_2       = "vehicle",
#'                     group_by      = "cell_type",
#'                     group_value   = "T cell")
#' head(res$results)
#'
#' # Every cluster at once
#' res <- PseudobulkDE(obj,
#'                     sample_col    = "orig.ident",
#'                     condition_col = "treatment",
#'                     ident_1       = "drug",
#'                     ident_2       = "vehicle",
#'                     cluster_col   = "seurat_clusters")
#' names(res)
#' head(res$`0`$results)
#'
#' # Batch-adjusted design
#' res <- PseudobulkDE(obj,
#'                     sample_col    = "orig.ident",
#'                     condition_col = "treatment",
#'                     ident_1       = "drug",
#'                     ident_2       = "vehicle",
#'                     cluster_col   = "seurat_clusters",
#'                     design        = ~ batch + condition,
#'                     contrast      = c("condition", "drug", "vehicle"))
#' }
#' @importFrom Seurat AggregateExpression DefaultAssay
#' @importFrom SeuratObject LayerData
#' @export
PseudobulkDE <- function(obj,
                         sample_col,
                         condition_col,
                         ident_1                   = NULL,
                         ident_2                   = NULL,
                         cluster_col               = NULL,
                         clusters                  = NULL,
                         group_by                  = NULL,
                         group_value               = NULL,
                         assay                     = NULL,
                         design                    = NULL,
                         contrast                  = NULL,
                         min_cells_per_sample      = 10,
                         min_samples_per_condition = 2,
                         min_counts_per_gene       = 10,
                         verbose                   = TRUE) {

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required. Install with: ",
         "BiocManager::install('DESeq2').")
  }
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  for (col in c(sample_col, condition_col)) {
    if (!col %in% colnames(obj@meta.data)) {
      stop("Column '", col, "' not found in obj@meta.data.")
    }
  }
  if (is.null(ident_1) || is.null(ident_2)) {
    if (is.null(contrast)) {
      stop("Provide `ident_1` and `ident_2`, or supply a full `contrast`.")
    }
    ident_1 <- contrast[2]
    ident_2 <- contrast[3]
  }
  if (is.null(design)) {
    design <- stats::as.formula("~ condition")
  } else if (is.character(design)) {
    design <- stats::as.formula(design)
  }
  if (is.null(contrast)) {
    contrast <- c("condition", ident_1, ident_2)
  }

  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay

  # ---- Multi-cluster dispatch ---------------------------------------------
  if (!is.null(cluster_col)) {
    if (!cluster_col %in% colnames(obj@meta.data)) {
      stop("`cluster_col` '", cluster_col, "' not found in metadata.")
    }
    all_clu <- sort(unique(as.character(obj@meta.data[[cluster_col]])))
    if (is.null(clusters)) clusters <- all_clu
    clusters <- intersect(clusters, all_clu)

    out <- list()
    for (cl in clusters) {
      if (isTRUE(verbose)) {
        message(sprintf("--- Cluster '%s' ---", cl))
      }
      one <- tryCatch(
        .pseudobulk_one(
          obj, sample_col, condition_col, ident_1, ident_2,
          group_by = cluster_col, group_value = cl,
          assay = a, design = design, contrast = contrast,
          min_cells_per_sample = min_cells_per_sample,
          min_samples_per_condition = min_samples_per_condition,
          min_counts_per_gene = min_counts_per_gene,
          verbose = verbose
        ),
        error = function(e) {
          message(sprintf("  Skipping cluster '%s': %s", cl,
                          conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(one)) out[[as.character(cl)]] <- one
    }
    return(out)
  }

  # ---- Single-group mode --------------------------------------------------
  .pseudobulk_one(
    obj, sample_col, condition_col, ident_1, ident_2,
    group_by = group_by, group_value = group_value,
    assay = a, design = design, contrast = contrast,
    min_cells_per_sample = min_cells_per_sample,
    min_samples_per_condition = min_samples_per_condition,
    min_counts_per_gene = min_counts_per_gene,
    verbose = verbose
  )
}


# ============================================================================
# Internal: run DESeq2 on one (sub)set of the Seurat object.
# Returns list(results, normalized_counts).
# ============================================================================
#' @keywords internal
#' @noRd
.pseudobulk_one <- function(obj, sample_col, condition_col, ident_1, ident_2,
                            group_by, group_value, assay, design, contrast,
                            min_cells_per_sample, min_samples_per_condition,
                            min_counts_per_gene, verbose) {

  md <- obj@meta.data

  # Optional single-group subset
  if (!is.null(group_by)) {
    if (is.null(group_value)) {
      stop("`group_value` required when `group_by` is supplied.")
    }
    if (!group_by %in% colnames(md)) {
      stop("Column '", group_by, "' not found in obj@meta.data.")
    }
    keep <- rownames(md)[as.character(md[[group_by]]) == group_value]
    if (length(keep) == 0) {
      stop("No cells with ", group_by, " == '", group_value, "'.")
    }
    obj <- subset(obj, cells = keep)
    md  <- obj@meta.data
  }

  # Restrict to the two conditions being contrasted
  keep_cond <- md[[condition_col]] %in% c(ident_1, ident_2)
  if (!any(keep_cond)) {
    stop("No cells in either '", ident_1, "' or '", ident_2, "'.")
  }
  obj <- subset(obj, cells = rownames(md)[keep_cond])
  md  <- obj@meta.data

  # Aggregate to pseudobulk. AggregateExpression can return a plain matrix,
  # a dgCMatrix, or an Assay / Assay5 depending on Seurat version — coerce
  # to a base integer matrix so downstream operations work uniformly.
  agg_raw <- Seurat::AggregateExpression(
    obj,
    assays        = assay,
    group.by      = c(sample_col, condition_col),
    return.seurat = FALSE,
    slot          = "counts"
  )[[assay]]
  if (inherits(agg_raw, c("Assay", "Assay5"))) {
    agg_raw <- SeuratObject::LayerData(agg_raw, layer = "counts")
  }
  agg <- as.matrix(agg_raw)
  storage.mode(agg) <- "integer"

  # Reconstruct (sample, condition) from column names. AggregateExpression
  # joins group.by values with "_"; if sample IDs contain "_" themselves,
  # anchor on the two known condition labels at the end.
  condition_pat <- paste0("_(", ident_1, "|", ident_2, ")$")
  if (all(grepl(condition_pat, colnames(agg)))) {
    sample_ids    <- sub(condition_pat, "", colnames(agg))
    condition_ids <- sub("^.*_", "", colnames(agg))
  } else {
    parts <- strsplit(colnames(agg), "_", fixed = TRUE)
    sample_ids    <- vapply(parts, `[`, character(1), 1)
    condition_ids <- vapply(parts, `[`, character(1), 2)
  }
  coldata <- data.frame(
    sample    = sample_ids,
    condition = factor(condition_ids, levels = c(ident_2, ident_1)),
    row.names = colnames(agg),
    stringsAsFactors = FALSE
  )
  # Expose the original condition column name too (in case `design` uses it)
  coldata[[condition_col]] <- coldata$condition

  # Preserve additional per-sample metadata (batch, sex, ...) for custom
  # designs like `~ batch + condition`.
  extra_cols <- setdiff(colnames(md),
                        c(sample_col, condition_col, "nCount_RNA",
                          "nFeature_RNA", "orig.ident"))
  if (length(extra_cols)) {
    per_sample <- unique(md[, c(sample_col, extra_cols), drop = FALSE])
    per_sample <- per_sample[!duplicated(per_sample[[sample_col]]), ,
                             drop = FALSE]
    ix <- match(sample_ids, per_sample[[sample_col]])
    for (ec in extra_cols) {
      if (!(ec %in% colnames(coldata))) coldata[[ec]] <- per_sample[[ec]][ix]
    }
  }

  # Drop pseudobulk samples below the cell-count threshold
  cells_per_pb <- table(paste(md[[sample_col]], md[[condition_col]],
                              sep = "_"))
  cells_per_col <- cells_per_pb[colnames(agg)]
  keep_cols <- !is.na(cells_per_col) & cells_per_col >= min_cells_per_sample
  if (isTRUE(verbose) && sum(keep_cols) < ncol(agg)) {
    message(sprintf("  Dropping %d pseudobulk sample(s) with < %d cells",
                    ncol(agg) - sum(keep_cols), min_cells_per_sample))
  }
  agg     <- agg[, keep_cols, drop = FALSE]
  coldata <- coldata[keep_cols, , drop = FALSE]

  n_per_cond <- table(coldata$condition)
  if (any(n_per_cond < min_samples_per_condition)) {
    stop("Condition counts after filtering: ",
         paste(names(n_per_cond), "=", n_per_cond, collapse = ", "),
         "; need >= ", min_samples_per_condition, " per condition.")
  }

  # Low-count gene filter
  if (requireNamespace("edgeR", quietly = TRUE)) {
    keep_g <- edgeR::filterByExpr(agg, group = factor(coldata$condition))
  } else {
    keep_g <- Matrix::rowSums(agg) >= min_counts_per_gene
  }
  agg <- agg[keep_g, , drop = FALSE]

  if (isTRUE(verbose)) {
    message(sprintf("  Running DESeq2 on %d genes x %d pseudobulk samples",
                    nrow(agg), ncol(agg)))
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = agg,
    colData   = coldata,
    design    = design
  )
  dds <- DESeq2::DESeq(dds, quiet = !isTRUE(verbose))
  res <- DESeq2::results(dds, contrast = contrast)

  results_df <- as.data.frame(res)
  results_df$gene <- rownames(results_df)
  names(results_df)[names(results_df) == "log2FoldChange"] <- "log2FC"
  results_df <- results_df[, c("gene", "baseMean", "log2FC", "lfcSE",
                               "stat", "pvalue", "padj")]
  results_df <- results_df[order(results_df$padj, results_df$pvalue,
                                 na.last = TRUE), ]
  rownames(results_df) <- NULL

  # DESeq2 size-factor-normalized pseudobulk counts, for downstream plotting
  # (per-gene boxplots etc.) without re-running the aggregation.
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)
  norm_counts_df <- as.data.frame(norm_counts)
  norm_counts_df$gene <- rownames(norm_counts_df)
  norm_counts_df <- norm_counts_df[, c("gene",
                                       setdiff(colnames(norm_counts_df), "gene"))]
  rownames(norm_counts_df) <- NULL

  if (isTRUE(verbose)) {
    n_sig <- sum(results_df$padj < 0.05, na.rm = TRUE)
    message(sprintf("  %d genes; %d significant at padj < 0.05.",
                    nrow(results_df), n_sig))
  }

  list(
    results           = results_df,
    normalized_counts = norm_counts_df
  )
}
