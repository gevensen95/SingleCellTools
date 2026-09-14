#' Gene-set / pathway enrichment on DE results
#'
#' Takes a DE results table -- typically from \code{\link{PseudobulkDE}} or
#' \code{Seurat::FindMarkers} -- and runs gene-set enrichment via
#' \code{fgsea} (rank-based GSEA over every tested gene) or \code{enrichR}
#' (hypergeometric test against a curated significant-gene list, via the
#' Enrichr web API).
#'
#' \strong{Choosing a method.} \code{"fgsea"} (default) ranks every gene in
#' \code{de_result} by \code{sign(log2FC) * -log10(pvalue)} and tests
#' whether genes in each set from \code{gene_sets} are skewed toward either
#' end of that ranking -- it uses the whole result table, so genes with a
#' modest but consistent effect still contribute signal, and it needs no
#' network access. \code{"enrichr"} instead takes only the genes passing
#' \code{p_threshold}/\code{fc_threshold} and tests them against Enrichr's
#' hosted gene-set libraries (\code{databases}) over the network -- simpler
#' to call (no \code{gene_sets} object to assemble, e.g. via
#' \code{msigdbr}) but throws away everything below the significance cutoff
#' and requires internet access.
#'
#' \code{de_result} may also be the named list \code{\link{PseudobulkDE}}
#' returns in multi-cluster mode (\code{cluster_col} set there); each
#' cluster's \code{results} table is enriched separately and the return
#' value is a named list of the same shape, matching
#' \code{\link{PseudobulkDE}}'s own convention for that mode.
#'
#' @param de_result A DE results data frame, a two-element list with a
#'   \code{results} element (as \code{\link{PseudobulkDE}} returns in
#'   single-group mode), or a named list of either (multi-cluster mode).
#' @param gene_sets \code{method = "fgsea"} only: a named list of character
#'   vectors, pathway name -> gene symbols. Supply this OR \code{collection}
#'   (one of them is required for \code{method = "fgsea"}); ignored for
#'   \code{method = "enrichr"}, which uses \code{databases} instead.
#' @param collection \code{method = "fgsea"} only: one or more named
#'   collections to fetch via \code{\link{FetchGeneSets}} instead of
#'   supplying \code{gene_sets} yourself -- e.g. \code{"go_bp"},
#'   \code{"kegg"}, \code{"reactome"}, \code{"hallmark"}. See
#'   \code{?FetchGeneSets} for the full list.
#' @param species \code{collection} only: passed to
#'   \code{\link{FetchGeneSets}} -- \strong{required} when \code{collection}
#'   is supplied (e.g. \code{"Homo sapiens"}, \code{"Mus musculus"}). There
#'   is no default and no auto-detection from \code{de_result}'s gene
#'   symbols: casing alone can't reliably distinguish some species (mouse
#'   vs. rat, for instance), so this is never guessed for you.
#' @param gene_set_source \code{collection} only: \code{"msigdbr"} (default)
#'   or \code{"godb"}, passed to \code{\link{FetchGeneSets}} as \code{source}.
#' @param method \code{"fgsea"} (default) or \code{"enrichr"}. See Details.
#' @param fc_col,p_col,gene_col Column names for log2FC / adjusted p-value /
#'   gene symbol. \code{NULL} (default) auto-detects the same way
#'   \code{\link{CompareMarkers}} does.
#' @param p_threshold,fc_threshold \code{method = "enrichr"} only: cutoffs
#'   defining the "significant gene" list sent to Enrichr. Defaults 0.05
#'   and 0 (any direction).
#' @param databases \code{method = "enrichr"} only: Enrichr library name(s).
#'   Default \code{c("GO_Biological_Process_2023", "KEGG_2021_Human")}.
#' @param min_size,max_size \code{method = "fgsea"} only: minimum/maximum
#'   gene-set size to test. Defaults 10 / 500.
#' @param plot Logical; if \code{TRUE} (default), also return a bar plot of
#'   the top \code{top_n_plot} gene sets.
#' @param top_n_plot Number of top gene sets to show in the plot. Default 20.
#' @param verbose Message progress. Default \code{TRUE}.
#' @return A list with \code{results} (a tidy data frame of gene sets and
#'   their enrichment statistics) and, if \code{plot = TRUE}, \code{plot}.
#'   In multi-cluster mode, a named list of those lists, one per cluster.
#' @examples
#' \dontrun{
#' de <- PseudobulkDE(obj, sample_col = "orig.ident",
#'                    condition_col = "treatment",
#'                    ident_1 = "drug", ident_2 = "vehicle")
#'
#' # fgsea, offline -- fetch the pathway list directly
#' # (species is always required when using `collection`, never guessed)
#' pe <- RunPathwayEnrichment(de, collection = "hallmark", species = "Homo sapiens")
#' pe$plot
#'
#' # Combine collections, or supply your own gene_sets instead
#' pe2 <- RunPathwayEnrichment(de, collection = c("go_bp", "kegg"),
#'                             species = "Homo sapiens")
#'
#' # enrichR, needs internet, no gene_sets object required
#' pe3 <- RunPathwayEnrichment(de, method = "enrichr")
#' }
#' @importFrom ggplot2 aes coord_flip element_text geom_col ggplot labs reorder theme
#' @export
RunPathwayEnrichment <- function(de_result,
                                 gene_sets       = NULL,
                                 collection      = NULL,
                                 species         = NULL,
                                 gene_set_source = c("msigdbr", "godb"),
                                 method          = c("fgsea", "enrichr"),
                                 fc_col          = NULL,
                                 p_col           = NULL,
                                 gene_col        = NULL,
                                 p_threshold     = 0.05,
                                 fc_threshold    = 0,
                                 databases       = c("GO_Biological_Process_2023",
                                                     "KEGG_2021_Human"),
                                 min_size        = 10,
                                 max_size        = 500,
                                 plot            = TRUE,
                                 top_n_plot      = 20,
                                 verbose         = TRUE) {

  method <- match.arg(method)
  gene_set_source <- match.arg(gene_set_source)
  if (method == "fgsea") {
    if (!is.null(gene_sets) && !is.null(collection)) {
      stop("Supply only one of `gene_sets` or `collection`, not both.")
    }
    if (is.null(gene_sets) && is.null(collection)) {
      stop("`gene_sets` (a named list of pathway -> gene symbols) or ",
           "`collection` (e.g. 'go_bp', 'kegg') is required for ",
           "method = 'fgsea'. See ?FetchGeneSets for valid collection names.")
    }
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
  }

  # ---- Multi-cluster dispatch (mirrors PseudobulkDE's own list-of-lists) --
  is_pbde_result <- is.list(de_result) && !is.data.frame(de_result) &&
                    "results" %in% names(de_result)
  is_multi <- is.list(de_result) && !is.data.frame(de_result) && !is_pbde_result

  if (is_multi) {
    out <- list()
    for (nm in names(de_result)) {
      if (isTRUE(verbose)) message(sprintf("--- Cluster '%s' ---", nm))
      one <- tryCatch(
        RunPathwayEnrichment(de_result[[nm]], gene_sets = gene_sets,
                             method = method, fc_col = fc_col, p_col = p_col,
                             gene_col = gene_col, p_threshold = p_threshold,
                             fc_threshold = fc_threshold, databases = databases,
                             min_size = min_size, max_size = max_size,
                             plot = plot, top_n_plot = top_n_plot,
                             verbose = verbose),
        error = function(e) {
          message(sprintf("  Skipping '%s': %s", nm, conditionMessage(e)))
          NULL
        }
      )
      if (!is.null(one)) out[[nm]] <- one
    }
    return(out)
  }

  df <- if (is_pbde_result) de_result$results else de_result
  if (!is.data.frame(df)) {
    stop("`de_result` must be a DE results data frame, the list ",
         "PseudobulkDE() returns, or a named list of either.")
  }

  fc_c   <- .detect_col(fc_col, df,
                        c("log2FC", "avg_log2FC", "logFC", "log2FoldChange"),
                        "log2 fold change")
  p_c    <- .detect_col(p_col, df,
                        c("padj", "p_val_adj", "FDR", "adj.P.Val"),
                        "adjusted p-value")
  gene_c <- if (is.null(gene_col)) {
    cand <- intersect(c("gene", "Gene", "gene_symbol"), colnames(df))[1]
    if (is.na(cand)) { df$.gene <- rownames(df); ".gene" } else cand
  } else gene_col

  if (method == "fgsea") {
    if (!requireNamespace("fgsea", quietly = TRUE)) {
      stop("'fgsea' is required for method = 'fgsea'. Install with ",
           "BiocManager::install('fgsea').")
    }
    rank_stat <- sign(df[[fc_c]]) * -log10(pmax(df[[p_c]], 1e-300))
    names(rank_stat) <- df[[gene_c]]
    rank_stat <- sort(rank_stat[!is.na(rank_stat)], decreasing = TRUE)

    if (isTRUE(verbose)) {
      message(sprintf("--- Running fgsea (%d gene sets, %d ranked genes) ---",
                      length(gene_sets), length(rank_stat)))
    }
    res <- fgsea::fgsea(pathways = gene_sets, stats = rank_stat,
                        minSize = min_size, maxSize = max_size)
    res <- as.data.frame(res)
    res$leadingEdge <- vapply(res$leadingEdge, paste, character(1), collapse = ",")
    res <- res[order(res$padj, res$pval), ]
    rownames(res) <- NULL
    score_col <- "NES"

  } else {
    if (!requireNamespace("enrichR", quietly = TRUE)) {
      stop("'enrichR' is required for method = 'enrichr'. Install with ",
           "install.packages('enrichR').")
    }
    sig <- df[[p_c]] < p_threshold & abs(df[[fc_c]]) > fc_threshold
    sig[is.na(sig)] <- FALSE
    genes <- unique(df[[gene_c]][sig])
    if (length(genes) == 0) {
      stop("No genes pass p_threshold/fc_threshold to send to Enrichr.")
    }
    if (isTRUE(verbose)) {
      message(sprintf("--- Running enrichR (%d genes, %d database(s)) ---",
                      length(genes), length(databases)))
    }
    res_list <- enrichR::enrichr(genes, databases)
    res <- do.call(rbind, lapply(names(res_list), function(dbn) {
      d <- res_list[[dbn]]
      if (nrow(d) == 0) return(NULL)
      d$database <- dbn
      d
    }))
    res <- res[order(res$Adjusted.P.value), ]
    rownames(res) <- NULL
    score_col <- "Combined.Score"
  }

  out <- list(results = res)

  if (isTRUE(plot) && nrow(res) > 0) {
    top <- res[seq_len(min(top_n_plot, nrow(res))), , drop = FALSE]
    name_col <- if (method == "fgsea") "pathway" else "Term"
    top[[name_col]] <- factor(top[[name_col]], levels = rev(top[[name_col]]))
    score <- name <- NULL  # NSE
    top$.score <- top[[score_col]]
    top$.name  <- top[[name_col]]
    p <- ggplot2::ggplot(top, ggplot2::aes(x = .score, y = .name)) +
      ggplot2::geom_col(fill = "#4393C3") +
      ggplot2::labs(x = score_col, y = NULL,
                   title = sprintf("Top %d gene sets (%s)", nrow(top), method)) +
      Ol_Reliable()
    out$plot <- p
  }

  out
}
