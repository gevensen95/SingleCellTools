#' Pseudobulk differential expression between two conditions
#'
#' Aggregates single-cell counts per (sample, group) pair and runs DESeq2 to
#' test the contrast between two levels of a condition column. This is the
#' statistically correct way to compare conditions across donors with
#' single-cell data &mdash; per-cell Wilcoxon DE treats each cell as
#' independent and dramatically inflates type-I error when there are multiple
#' cells per donor.
#'
#' @param obj A Seurat object.
#' @param sample_col Metadata column identifying the biological replicate
#'   (e.g. donor, mouse, well). Required.
#' @param condition_col Metadata column identifying the condition to test
#'   (e.g. treatment vs vehicle).
#' @param ident_1,ident_2 The two levels of \code{condition_col} to contrast.
#'   \code{ident_1} is the "numerator" in log-fold-changes
#'   (positive logFC = higher in \code{ident_1}).
#' @param group_by Optional metadata column to subset to a single group
#'   before pseudobulking (e.g. a cluster or cell type). When supplied,
#'   only cells with \code{group_by == group_value} are used.
#' @param group_value The value within \code{group_by} to keep. Required if
#'   \code{group_by} is supplied.
#' @param assay Assay to aggregate. Defaults to DefaultAssay.
#' @param min_cells_per_sample Drop pseudobulk samples that aggregate fewer
#'   than this many cells (the pseudobulk count for that sample is unstable).
#'   Default 10.
#' @param min_samples_per_condition Refuse to run if either condition has
#'   fewer than this many usable pseudobulk samples. Default 2.
#' @return A list with two elements:
#' \describe{
#'   \item{\code{results}}{A data frame of DESeq2 results sorted by adjusted
#'     p-value, with gene names in the \code{gene} column.}
#'   \item{\code{normalized_counts}}{A data frame of DESeq2
#'     size-factor-normalized pseudobulk counts (\code{counts(dds,
#'     normalized = TRUE)}), one row per gene (\code{gene} column) and one
#'     column per pseudobulk sample retained after filtering.}
#' }
#' @importFrom Seurat AggregateExpression DefaultAssay
#' @importFrom SeuratObject LayerData
#' @export
PseudobulkDE <- function(obj,
                         sample_col,
                         condition_col,
                         ident_1,
                         ident_2,
                         group_by                   = NULL,
                         group_value                = NULL,
                         assay                      = NULL,
                         min_cells_per_sample       = 10,
                         min_samples_per_condition  = 2) {

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required for PseudobulkDE. Install with: ",
         "BiocManager::install('DESeq2')")
  }
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")

  md <- obj@meta.data
  for (col in c(sample_col, condition_col)) {
    if (!col %in% colnames(md)) {
      stop("Column '", col, "' not found in obj@meta.data")
    }
  }

  # Optional subset to a single group (cluster / celltype)
  if (!is.null(group_by)) {
    if (is.null(group_value)) {
      stop("`group_value` required when `group_by` is supplied.")
    }
    if (!group_by %in% colnames(md)) {
      stop("Column '", group_by, "' not found in obj@meta.data")
    }
    keep <- rownames(md)[md[[group_by]] == group_value]
    if (length(keep) == 0) {
      stop("No cells with ", group_by, " == '", group_value, "'.")
    }
    obj <- subset(obj, cells = keep)
    md  <- obj@meta.data
  }

  # Keep only the two conditions being contrasted
  keep_cond <- md[[condition_col]] %in% c(ident_1, ident_2)
  if (!any(keep_cond)) {
    stop("No cells in either '", ident_1, "' or '", ident_2,
         "' for condition column '", condition_col, "'.")
  }
  obj <- subset(obj, cells = rownames(md)[keep_cond])
  md  <- obj@meta.data

  # Aggregate to pseudobulk. AggregateExpression returns a list keyed by
  # assay; the value can be a dense matrix, a dgCMatrix, or even an Assay
  # depending on Seurat version. Coerce defensively to a base matrix so
  # downstream subsetting works regardless of the underlying class.
  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay
  agg_raw <- Seurat::AggregateExpression(
    obj,
    assays        = a,
    group.by      = c(sample_col, condition_col),
    return.seurat = FALSE,
    slot          = "counts"
  )[[a]]

  if (inherits(agg_raw, c("Assay", "Assay5"))) {
    agg_raw <- SeuratObject::LayerData(agg_raw, layer = "counts")
  }
  agg <- as.matrix(agg_raw)
  storage.mode(agg) <- "integer"

  # Reconstruct sample / condition from the column names.
  # AggregateExpression joins group.by values with "_"; if sample IDs already
  # contain "_", we anchor on the known condition values at the end.
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

  # Drop pseudobulk samples below the cell-count threshold
  cells_per_pb <- table(paste(md[[sample_col]], md[[condition_col]], sep = "_"))
  cells_per_col <- cells_per_pb[colnames(agg)]
  keep_cols <- !is.na(cells_per_col) & cells_per_col >= min_cells_per_sample
  if (sum(keep_cols) < ncol(agg)) {
    message(sprintf("Dropping %d pseudobulk sample(s) with < %d cells",
                    ncol(agg) - sum(keep_cols), min_cells_per_sample))
  }
  agg     <- agg[, keep_cols, drop = FALSE]
  coldata <- coldata[keep_cols, , drop = FALSE]

  n_per_cond <- table(coldata$condition)
  if (any(n_per_cond < min_samples_per_condition)) {
    stop("After filtering, condition counts are ",
         paste(names(n_per_cond), "=", n_per_cond, collapse = ", "),
         "; need >= ", min_samples_per_condition, " per condition.")
  }

  # Run DESeq2
  message(sprintf("Running DESeq2 on %d genes x %d pseudobulk samples",
                  nrow(agg), ncol(agg)))
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = agg,
    colData   = coldata,
    design    = ~ condition
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- DESeq2::results(dds, contrast = c("condition", ident_1, ident_2))

  out <- as.data.frame(res)
  out$gene <- rownames(out)
  out <- out[order(out$padj, out$pvalue), ]
  rownames(out) <- NULL
  out <- out[, c("gene", setdiff(colnames(out), "gene"))]

  # ---- Normalized counts ---------------------------------------------------
  # DESeq2 size-factor-normalized pseudobulk counts, for downstream plotting
  # (e.g. per-gene boxplots) without re-running the aggregation/DESeq2 steps.
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)
  norm_counts_df <- as.data.frame(norm_counts)
  norm_counts_df$gene <- rownames(norm_counts_df)
  norm_counts_df <- norm_counts_df[, c("gene", setdiff(colnames(norm_counts_df), "gene"))]
  rownames(norm_counts_df) <- NULL

  list(
    results           = out,
    normalized_counts = norm_counts_df
  )
}
