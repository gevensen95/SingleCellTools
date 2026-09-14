#' Single-cell-level gene-set enrichment (no pseudobulk aggregation)
#'
#' Scores every cell for every gene set in \code{gene_sets} (via
#' \code{UCell} or \code{AUCell}), then tests whether each gene set's
#' per-cell score distribution differs between two conditions
#' (\code{ident_1}/\code{ident_2}) or between one cluster and the rest
#' (\code{cluster_col}) -- entirely at cell resolution, with no
#' sample-level aggregation step.
#'
#' \strong{When to use this instead of \code{\link{PseudobulkDE}} +
#' \code{\link{RunPathwayEnrichment}}.} \code{PseudobulkDE} requires at
#' least \code{min_samples_per_condition} (default 2) independent samples
#' per condition -- the statistically correct way to compare conditions
#' across donors. This function tests per-\emph{cell} score distributions
#' directly, so it works even with a single animal/sample per condition,
#' where pseudobulk literally cannot be run at all. The tradeoff is real:
#' cells from the same sample are not independent observations, so a
#' Wilcoxon test across thousands of cells is more liberal than the same
#' test across a handful of pseudobulk samples -- p-values here should be
#' read as "detectable difference in this dataset," not treated as
#' generalizing to new donors the way \code{PseudobulkDE}'s do. Prefer
#' \code{PseudobulkDE} whenever you have >= 2 samples per condition.
#'
#' \strong{Scoring method.} \code{"ucell"} (default) uses
#' \code{UCell::AddModuleScore_UCell()} -- already a dependency of this
#' package (see \code{\link{ClassifyByReferenceCutoff}}), rank-based, and
#' robust to dropout. \code{"aucell"} uses
#' \code{AUCell::AUCell_buildRankings()} + \code{AUCell_calcAUC()}, the
#' classic SCENIC-ecosystem method.
#'
#' \strong{Comparison mode} is chosen by which arguments you supply, not a
#' separate flag: pass \code{cluster_col} for per-cluster one-vs-rest
#' testing, or \code{condition_col}/\code{ident_1}/\code{ident_2} (and
#' optionally \code{group_by}/\code{group_value} to restrict to a subset
#' first, e.g. one cell type) for a two-group contrast. Passing both, or
#' neither, is an error -- exactly one mode per call.
#'
#' @param obj A Seurat object.
#' @param gene_sets A named list of character vectors, gene-set name ->
#'   gene symbols (e.g. from \code{msigdbr}). Genes not present in
#'   \code{assay} are dropped from each set; sets left with fewer than
#'   \code{min_size} genes, or more than \code{max_size}, are excluded
#'   (matching \code{\link{RunPathwayEnrichment}}'s conventions). Supply
#'   this OR \code{collection}, not both.
#' @param collection One or more named collections to fetch via
#'   \code{\link{FetchGeneSets}} instead of supplying \code{gene_sets}
#'   yourself -- e.g. \code{"go_bp"}, \code{"go_mf"}, \code{"kegg"},
#'   \code{"reactome"}, \code{"hallmark"}. See \code{?FetchGeneSets} for the
#'   full list. \code{NULL} (default) requires \code{gene_sets} instead.
#' @param species \code{collection} only: passed to
#'   \code{\link{FetchGeneSets}} -- \strong{required} when \code{collection}
#'   is supplied (e.g. \code{"Homo sapiens"}, \code{"Mus musculus"}). There
#'   is no default and no auto-detection from \code{obj}'s gene symbols:
#'   casing alone can't reliably distinguish some species (mouse vs. rat,
#'   for instance), so this is never guessed for you.
#' @param gene_set_source \code{collection} only: \code{"msigdbr"} (default)
#'   or \code{"godb"}, passed to \code{\link{FetchGeneSets}} as \code{source}.
#' @param method \code{"ucell"} (default) or \code{"aucell"}. See Details.
#' @param assay Assay to score. Default \code{DefaultAssay(obj)}.
#' @param condition_col,ident_1,ident_2 Two-group mode: metadata column and
#'   the two levels to contrast (\code{ident_1} is the positive-\code{mean_diff}
#'   side, matching \code{\link{PseudobulkDE}}'s \code{ident_1}/\code{ident_2}).
#' @param group_by,group_value Two-group mode only: optional pre-filter to
#'   one subset (e.g. one cell type) before contrasting, exactly like
#'   \code{\link{PseudobulkDE}}'s same-named arguments.
#' @param cluster_col Per-cluster mode: metadata column of cluster/cell-type
#'   labels. Triggers one-vs-rest testing instead of a two-group contrast.
#' @param clusters Per-cluster mode only: which \code{cluster_col} levels to
#'   test; \code{NULL} (default) tests every level.
#' @param min_size,max_size Gene-set size bounds after filtering to genes
#'   present in \code{assay}. Defaults 5 / 500.
#' @param store_scores Logical; if \code{TRUE} (default), also write each
#'   gene set's per-cell score as a \code{<gene_set>_<method>} metadata
#'   column (in addition to the full matrix always stored -- see Return).
#' @param top_n_plot Number of top gene sets shown in the plot. Default 20.
#' @param plot Logical; if \code{TRUE} (default), also build a plot -- a
#'   bar chart of top gene sets by significance (two-group mode) or a
#'   cluster x gene-set heatmap of the union of each cluster's top hits
#'   (per-cluster mode), styled like \code{\link{RCTDCompositionHeatmap}}.
#' @param verbose Message progress. Default \code{TRUE}.
#' @return \code{obj}, with \code{obj@misc$sc_gsea} populated:
#'   \code{list(mode, scores, results, plot)} where \code{scores} is the
#'   full cell x gene-set score matrix; \code{results} is a data frame
#'   (two-group mode: \code{gene_set}, \code{mean_<ident_1>},
#'   \code{mean_<ident_2>}, \code{mean_diff}, \code{pvalue}, \code{padj}) or
#'   a named list of such data frames, one per cluster (per-cluster mode);
#'   and \code{plot} is the bar plot or heatmap described above. If
#'   \code{store_scores = TRUE}, per-cell \code{<gene_set>_<method>}
#'   metadata columns are also added.
#' @examples
#' \dontrun{
#' # Fetch the gene sets directly -- no separate FetchGeneSets() call needed
#' # (species is always required when using `collection`, never guessed)
#' obj <- RunSingleCellGSEA(obj, collection = "go_bp", species = "Homo sapiens",
#'                          condition_col = "treatment",
#'                          ident_1 = "drug", ident_2 = "vehicle")
#' obj@misc$sc_gsea$results
#' obj@misc$sc_gsea$plot
#'
#' # Combine collections, or supply your own gene_sets instead
#' obj <- RunSingleCellGSEA(obj, collection = c("kegg", "reactome"),
#'                          species = "Homo sapiens",
#'                          condition_col = "treatment",
#'                          ident_1 = "drug", ident_2 = "vehicle")
#'
#' # Per-cluster one-vs-rest instead
#' obj <- RunSingleCellGSEA(obj, collection = "hallmark", species = "Mus musculus",
#'                          cluster_col = "seurat_clusters")
#' obj@misc$sc_gsea$results$`0`
#' obj@misc$sc_gsea$plot
#' }
#' @importFrom Seurat DefaultAssay GetAssayData
#' @importFrom ggplot2 aes element_blank element_text geom_col geom_tile ggplot labs scale_fill_gradientn theme
#' @export
RunSingleCellGSEA <- function(obj,
                              gene_sets       = NULL,
                              collection      = NULL,
                              species         = NULL,
                              gene_set_source = c("msigdbr", "godb"),
                              method          = c("ucell", "aucell"),
                              assay           = NULL,
                              condition_col   = NULL,
                              ident_1         = NULL,
                              ident_2         = NULL,
                              group_by        = NULL,
                              group_value     = NULL,
                              cluster_col     = NULL,
                              clusters        = NULL,
                              min_size        = 5,
                              max_size        = 500,
                              store_scores    = TRUE,
                              top_n_plot      = 20,
                              plot            = TRUE,
                              verbose         = TRUE) {

  method <- match.arg(method)
  gene_set_source <- match.arg(gene_set_source)
  .assert_seurat(obj)

  if (is.null(gene_sets) && is.null(collection)) {
    stop("Supply either `gene_sets` (a named list) or `collection` (e.g. ",
         "'go_bp', 'kegg', 'reactome') -- see ?FetchGeneSets for valid names.")
  }
  if (!is.null(gene_sets) && !is.null(collection)) {
    stop("Supply only one of `gene_sets` or `collection`, not both.")
  }

  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay
  present <- rownames(obj[[a]])

  if (!is.null(collection)) {
    if (is.null(species)) {
      stop("`species` is required when `collection` is supplied (e.g. ",
           "'Homo sapiens', 'Mus musculus') -- there is no default and no ",
           "auto-detection, since gene symbol casing alone can't reliably ",
           "distinguish some species (mouse vs. rat, for instance).")
    }
    gene_sets <- FetchGeneSets(collection, species = species,
                               source = gene_set_source, verbose = verbose)
  }
  if (!is.list(gene_sets) || is.null(names(gene_sets)) || any(names(gene_sets) == "")) {
    stop("`gene_sets` must be a named list of character vectors (gene-set ",
         "name -> gene symbols).")
  }

  two_group_mode <- !is.null(condition_col) && !is.null(ident_1) && !is.null(ident_2)
  cluster_mode    <- !is.null(cluster_col)
  if (two_group_mode == cluster_mode) {
    stop("Pass exactly one of: cluster_col (per-cluster one-vs-rest mode), ",
         "or condition_col + ident_1 + ident_2 (two-group contrast mode). ",
         "Got ", if (two_group_mode) "both" else "neither", ".")
  }
  if (cluster_mode && !cluster_col %in% colnames(obj@meta.data)) {
    stop("`cluster_col` '", cluster_col, "' not found in obj@meta.data.")
  }
  if (two_group_mode && !condition_col %in% colnames(obj@meta.data)) {
    stop("`condition_col` '", condition_col, "' not found in obj@meta.data.")
  }

  sets <- .sc_gsea_filter_sets(gene_sets, present, min_size, max_size, verbose)

  scores <- .sc_gsea_score_matrix(obj, sets, method, a, verbose)

  if (isTRUE(store_scores)) {
    cols <- paste0(colnames(scores), "_", method)
    score_df <- as.data.frame(scores[colnames(obj), , drop = FALSE])
    colnames(score_df) <- cols
    obj@meta.data[, cols] <- score_df
  }

  if (two_group_mode) {
    res <- .sc_gsea_two_group(obj, scores, condition_col, ident_1, ident_2,
                              group_by, group_value, verbose)
    out <- list(mode = "two_group", scores = res$scores, results = res$results)
    if (isTRUE(plot) && nrow(res$results) > 0) {
      out$plot <- .sc_gsea_two_group_plot(res$results, ident_1, ident_2, top_n_plot)
    }
  } else {
    res_list <- .sc_gsea_per_cluster(obj, scores, cluster_col, clusters, verbose)
    out <- list(mode = "per_cluster", scores = scores, results = res_list)
    if (isTRUE(plot) && length(res_list) > 0) {
      out$plot <- .sc_gsea_per_cluster_heatmap(scores, obj@meta.data[[cluster_col]],
                                               res_list, top_n_plot)
    }
  }

  obj@misc$sc_gsea <- out
  obj
}


# ============================================================================
# Internal helpers for RunSingleCellGSEA()
# ============================================================================

#' @keywords internal
#' @noRd
.sc_gsea_filter_sets <- function(gene_sets, present_genes, min_size, max_size, verbose) {
  filtered <- lapply(gene_sets, function(g) intersect(g, present_genes))
  sizes <- vapply(filtered, length, integer(1))
  keep <- sizes >= min_size & sizes <= max_size
  if (isTRUE(verbose)) {
    message(sprintf("--- %d/%d gene set(s) kept after filtering to genes in assay (size %d-%d) ---",
                    sum(keep), length(gene_sets), min_size, max_size))
  }
  if (sum(keep) == 0) {
    stop("No gene sets survived filtering to genes present in `obj`'s assay ",
         "and the [min_size, max_size] bounds.")
  }
  filtered[keep]
}

#' @keywords internal
#' @noRd
.sc_gsea_score_matrix <- function(obj, gene_sets, method, assay, verbose) {
  if (method == "ucell") {
    if (!requireNamespace("UCell", quietly = TRUE)) {
      stop("'UCell' is required for method = 'ucell'. Install with ",
           "BiocManager::install('UCell').")
    }
    if (isTRUE(verbose)) {
      message(sprintf("--- Scoring %d gene set(s) with UCell (assay = '%s') ---",
                      length(gene_sets), assay))
    }
    scored <- UCell::AddModuleScore_UCell(obj, features = gene_sets, assay = assay,
                                          name = "_scgsea_tmp")
    cols <- paste0(names(gene_sets), "_scgsea_tmp")
    scores <- as.matrix(scored@meta.data[, cols, drop = FALSE])
    colnames(scores) <- names(gene_sets)
    scores

  } else {
    if (!requireNamespace("AUCell", quietly = TRUE)) {
      stop("'AUCell' is required for method = 'aucell'. Install with ",
           "BiocManager::install('AUCell').")
    }
    if (isTRUE(verbose)) {
      message(sprintf("--- Scoring %d gene set(s) with AUCell (assay = '%s') ---",
                      length(gene_sets), assay))
    }
    mat <- as.matrix(Seurat::GetAssayData(obj, assay = assay, layer = "counts"))
    rankings <- AUCell::AUCell_buildRankings(mat, plotStats = FALSE,
                                             verbose = isTRUE(verbose))
    auc <- AUCell::AUCell_calcAUC(gene_sets, rankings, verbose = isTRUE(verbose))
    t(AUCell::getAUC(auc))
  }
}

#' @keywords internal
#' @noRd
.sc_gsea_wilcox_one <- function(x, y) {
  p <- tryCatch(stats::wilcox.test(x, y)$p.value, error = function(e) NA_real_)
  c(mean_x = mean(x, na.rm = TRUE), mean_y = mean(y, na.rm = TRUE), pvalue = p)
}

#' @keywords internal
#' @noRd
.sc_gsea_two_group <- function(obj, scores, condition_col, ident_1, ident_2,
                               group_by, group_value, verbose) {
  md <- obj@meta.data

  if (!is.null(group_by)) {
    if (is.null(group_value)) stop("`group_value` required when `group_by` is supplied.")
    if (!group_by %in% colnames(md)) stop("Column '", group_by, "' not found in obj@meta.data.")
    keep <- rownames(md)[as.character(md[[group_by]]) == group_value]
    if (length(keep) == 0) stop("No cells with ", group_by, " == '", group_value, "'.")
    md <- md[keep, , drop = FALSE]
  }

  keep_cond <- md[[condition_col]] %in% c(ident_1, ident_2)
  if (!any(keep_cond)) stop("No cells in either '", ident_1, "' or '", ident_2, "'.")
  md <- md[keep_cond, , drop = FALSE]

  scores_sub <- scores[rownames(md), , drop = FALSE]
  grp <- as.character(md[[condition_col]])

  if (isTRUE(verbose)) {
    message(sprintf("--- Testing %d gene set(s), %s (n=%d) vs %s (n=%d) ---",
                    ncol(scores_sub), ident_1, sum(grp == ident_1),
                    ident_2, sum(grp == ident_2)))
  }

  rows <- lapply(colnames(scores_sub), function(gs) {
    x <- scores_sub[grp == ident_1, gs]
    y <- scores_sub[grp == ident_2, gs]
    .sc_gsea_wilcox_one(x, y)
  })
  results <- as.data.frame(do.call(rbind, rows))
  results$gene_set <- colnames(scores_sub)
  colnames(results)[colnames(results) == "mean_x"] <- paste0("mean_", ident_1)
  colnames(results)[colnames(results) == "mean_y"] <- paste0("mean_", ident_2)
  results$mean_diff <- results[[paste0("mean_", ident_1)]] - results[[paste0("mean_", ident_2)]]
  results$padj <- stats::p.adjust(results$pvalue, method = "BH")
  results <- results[, c("gene_set", paste0("mean_", ident_1), paste0("mean_", ident_2),
                        "mean_diff", "pvalue", "padj")]
  results <- results[order(results$padj, results$pvalue), ]
  rownames(results) <- NULL

  list(scores = scores_sub, results = results)
}

#' @keywords internal
#' @noRd
.sc_gsea_per_cluster <- function(obj, scores, cluster_col, clusters, verbose) {
  all_clusters <- as.character(obj@meta.data[[cluster_col]])
  levels_use <- if (is.null(clusters)) sort(unique(all_clusters)) else intersect(clusters, unique(all_clusters))

  out <- list()
  for (cl in levels_use) {
    if (isTRUE(verbose)) message(sprintf("--- Cluster '%s' (one-vs-rest) ---", cl))
    in_cl <- all_clusters == cl
    rows <- lapply(colnames(scores), function(gs) {
      .sc_gsea_wilcox_one(scores[in_cl, gs], scores[!in_cl, gs])
    })
    res <- as.data.frame(do.call(rbind, rows))
    res$gene_set <- colnames(scores)
    colnames(res)[colnames(res) == "mean_x"] <- "mean_in"
    colnames(res)[colnames(res) == "mean_y"] <- "mean_out"
    res$mean_diff <- res$mean_in - res$mean_out
    res$padj <- stats::p.adjust(res$pvalue, method = "BH")
    res <- res[, c("gene_set", "mean_in", "mean_out", "mean_diff", "pvalue", "padj")]
    res <- res[order(res$padj, res$pvalue), ]
    rownames(res) <- NULL
    out[[cl]] <- res
  }
  out
}

#' @keywords internal
#' @noRd
.sc_gsea_two_group_plot <- function(results, ident_1, ident_2, top_n_plot) {
  top <- results[seq_len(min(top_n_plot, nrow(results))), , drop = FALSE]
  top$gene_set <- factor(top$gene_set, levels = rev(top$gene_set))
  top$neglog10padj <- -log10(pmax(top$padj, 1e-300))
  top$direction <- ifelse(top$mean_diff > 0, paste0("higher in ", ident_1),
                          paste0("higher in ", ident_2))
  direction <- neglog10padj <- gene_set <- NULL  # NSE
  ggplot2::ggplot(top, ggplot2::aes(x = neglog10padj, y = gene_set, fill = direction)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = "-log10(padj)", y = NULL, fill = NULL,
                 title = sprintf("Top %d gene sets: %s vs %s", nrow(top), ident_1, ident_2)) +
    Ol_Reliable()
}

#' @keywords internal
#' @noRd
.sc_gsea_per_cluster_heatmap <- function(scores, cluster_vec, res_list, top_n_plot) {
  cluster_vec <- as.character(cluster_vec)
  top_sets <- unique(unlist(lapply(res_list, function(d) {
    d$gene_set[seq_len(min(top_n_plot, nrow(d)))]
  })))

  clusters <- names(res_list)
  mat <- matrix(NA_real_, nrow = length(clusters), ncol = length(top_sets),
               dimnames = list(clusters, top_sets))
  for (cl in clusters) {
    in_cl <- cluster_vec == cl
    mat[cl, ] <- colMeans(scores[in_cl, top_sets, drop = FALSE], na.rm = TRUE)
  }
  mat_z <- scale(mat)
  mat_z[!is.finite(mat_z)] <- 0

  long <- data.frame(
    cluster   = rep(rownames(mat_z), times = ncol(mat_z)),
    gene_set  = rep(colnames(mat_z), each = nrow(mat_z)),
    value     = as.vector(mat_z),
    stringsAsFactors = FALSE
  )
  long$cluster  <- factor(long$cluster, levels = rownames(mat_z))
  long$gene_set <- factor(long$gene_set, levels = colnames(mat_z))

  cluster <- gene_set <- value <- NULL  # NSE
  ggplot2::ggplot(long, ggplot2::aes(x = cluster, y = gene_set, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11, "RdBu"))) +
    Ol_Reliable() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                  panel.grid  = ggplot2::element_blank()) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Z-score",
                 title = "Top gene sets per cluster (one-vs-rest)")
}
