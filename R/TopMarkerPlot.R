#' Top-N Marker Genes per Cluster, as an Annotated Dot Plot
#'
#' One call from a clustered Seurat object to a labeled marker dot plot: runs
#' \code{\link[Seurat]{FindAllMarkers}}, keeps the top \code{n} genes per
#' identity, and hands the result to \code{\link{MarkerPlot}} so the styling,
#' filtering, auto-sizing, and right-edge annotation labels all match every
#' other marker panel in this package.
#'
#' \code{FindAllMarkers} is the slow part. If you already have its output --
#' or want to tweak its arguments beyond what's exposed here -- pass it as
#' \code{markers} and this skips straight to the ranking and plotting.
#'
#' @section How genes are chosen:
#' Marker rows are filtered to \code{p_val_adj < max_padj}, then the top
#' \code{n} per cluster are taken by \code{rank_by}. A gene that lands in
#' more than one cluster's top-N is assigned to whichever cluster it scored
#' highest in and dropped from the others -- duplicated genes would otherwise
#' plot as repeated, identical rows. That means a cluster can end up with
#' fewer than \code{n} genes, which is usually a signal in itself (that
#' cluster has no markers of its own).
#'
#' @section A note on cluster ordering:
#' \code{MarkerPlot()} orders its annotation groups \emph{alphabetically}. With
#' plain numeric cluster IDs that puts cluster 10 between 1 and 2. Cluster
#' labels are therefore zero-padded to equal width ("00", "01", ... "11") when
#' the identities are numeric and there are 10 or more, so alphabetical order
#' matches numeric order. Pass your own labels via \code{group_by} (pointing
#' at a metadata column of cell-type names, say) and they are used verbatim.
#'
#' @param obj A Seurat object with identities set (or see \code{group_by}).
#' @param n Number of top marker genes to keep per cluster. Default \code{5}.
#' @param markers Optional pre-computed \code{FindAllMarkers()} data frame
#'   (needs \code{gene}, \code{cluster}, \code{avg_log2FC}, \code{p_val_adj}).
#'   \code{NULL} (default) runs \code{FindAllMarkers()} for you.
#' @param group_by Optional metadata column to set as identity before
#'   finding markers. \code{NULL} (default) uses the object's current
#'   \code{Idents()}. The object's identities are restored on exit.
#' @param rank_by One of \code{"avg_log2FC"} (default) or \code{"p_val_adj"}.
#'   Effect size or significance. Effect size is the usual choice for a
#'   visual panel; p-values in scRNA-seq are inflated by treating cells as
#'   independent replicates, so ranking by them mostly ranks by cluster size.
#' @param max_padj Adjusted p-value cutoff applied before ranking. Default
#'   \code{0.05}. Set to \code{1} to disable.
#' @param min_pct,logfc_threshold,only_pos Passed to
#'   \code{\link[Seurat]{FindAllMarkers}} when \code{markers} is \code{NULL}.
#'   Defaults \code{0.1}, \code{0.25}, \code{TRUE}.
#' @param cluster Passed to \code{\link{MarkerPlot}}. Default \code{FALSE}
#'   here (unlike \code{MarkerPlot()}'s own \code{TRUE}): re-ordering
#'   identities by correlation destroys the block-diagonal structure that
#'   makes a top-N-per-cluster panel readable in the first place.
#' @param ... Further arguments passed to \code{\link{MarkerPlot}} --
#'   \code{save_path}, \code{assay}, \code{maxsize}, \code{width}/
#'   \code{height}, \code{label.fontsize}, and so on.
#'
#' @return The \code{ggplot} returned by \code{\link{MarkerPlot}}, with two
#'   extra attributes: \code{attr(p, "markers")} (the full
#'   \code{FindAllMarkers} table) and \code{attr(p, "top_markers")} (the
#'   gene/cluster table actually plotted). Retrieve them with
#'   \code{attr(p, "top_markers")} rather than re-running the search.
#'
#' @examples
#' \dontrun{
#' # Simplest case
#' p <- TopMarkerPlot(obj)
#'
#' # 10 per cluster, saved at the auto-computed figure size
#' p <- TopMarkerPlot(obj, n = 10, save_path = "top10_markers.pdf")
#'
#' # Reuse an existing FindAllMarkers() run
#' mk <- Seurat::FindAllMarkers(obj, only.pos = TRUE)
#' p  <- TopMarkerPlot(obj, n = 5, markers = mk)
#'
#' # Group by annotated cell type instead of numeric clusters
#' p <- TopMarkerPlot(obj, group_by = "celltype")
#'
#' # Pull the table back out
#' utils::write.csv(attr(p, "top_markers"), "top_markers.csv", row.names = FALSE)
#' }
#'
#' @seealso \code{\link{MarkerPlot}}, \code{\link{MarkerPctPlot}}
#' @importFrom Seurat FindAllMarkers Idents Idents<-
#' @export
TopMarkerPlot <- function(obj,
                          n               = 5,
                          markers         = NULL,
                          group_by        = NULL,
                          rank_by         = c("avg_log2FC", "p_val_adj"),
                          max_padj        = 0.05,
                          min_pct         = 0.1,
                          logfc_threshold = 0.25,
                          only_pos        = TRUE,
                          cluster         = FALSE,
                          ...) {

  rank_by <- match.arg(rank_by)

  # ---- Validate -----------------------------------------------------------
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1) {
    stop("`n` must be a single positive number.")
  }
  n <- as.integer(n)

  if (!is.null(group_by)) {
    if (!group_by %in% colnames(obj@meta.data)) {
      stop("`group_by` column '", group_by, "' not found in metadata.")
    }
    old_idents <- Seurat::Idents(obj)
    on.exit(Seurat::Idents(obj) <- old_idents, add = TRUE)
    Seurat::Idents(obj) <- obj@meta.data[[group_by]]
  }
  if (length(unique(Seurat::Idents(obj))) < 2) {
    stop("Need at least 2 identities to find markers between. Set ",
         "`group_by`, or run clustering first.")
  }

  # ---- Markers ------------------------------------------------------------
  if (is.null(markers)) {
    message(sprintf("--- FindAllMarkers across %d identities ---",
                    length(unique(Seurat::Idents(obj)))))
    markers <- Seurat::FindAllMarkers(obj,
                                      only.pos        = only_pos,
                                      min.pct         = min_pct,
                                      logfc.threshold = logfc_threshold)
  }
  need <- c("gene", "cluster", "avg_log2FC", "p_val_adj")
  miss <- setdiff(need, colnames(markers))
  if (length(miss)) {
    stop("`markers` is missing column(s): ", paste(miss, collapse = ", "),
         ". Expected a FindAllMarkers() data frame.")
  }
  if (!nrow(markers)) stop("No markers found. Loosen logfc_threshold/min_pct.")

  # ---- Filter and rank ----------------------------------------------------
  m <- markers[!is.na(markers$p_val_adj) & markers$p_val_adj < max_padj, ,
               drop = FALSE]
  if (!nrow(m)) {
    stop("No markers passed p_val_adj < ", max_padj,
         ". Raise `max_padj`, or check the marker table.")
  }
  m$gene    <- as.character(m$gene)
  m$cluster <- as.character(m$cluster)

  # Rank within cluster: descending effect size, or ascending p-value.
  score <- if (rank_by == "avg_log2FC") -m$avg_log2FC else m$p_val_adj
  m     <- m[order(m$cluster, score), , drop = FALSE]
  top   <- do.call(rbind, lapply(split(m, m$cluster), utils::head, n))

  # A gene topping more than one cluster would plot as duplicate identical
  # rows, so keep only its best-scoring cluster. Sorting by the same score
  # first makes "best" well defined regardless of rank_by.
  top   <- top[order(if (rank_by == "avg_log2FC") -top$avg_log2FC
                     else top$p_val_adj), , drop = FALSE]
  n_dup <- sum(duplicated(top$gene))
  top   <- top[!duplicated(top$gene), , drop = FALSE]
  if (n_dup) {
    message(sprintf("  %d gene(s) topped more than one cluster; kept the ",
                    n_dup), "best-scoring cluster for each.")
  }

  # ---- Cluster labels -----------------------------------------------------
  # MarkerPlot() orders annotation groups alphabetically, so numeric IDs need
  # zero-padding or cluster 10 lands between 1 and 2.
  labs <- unique(top$cluster)
  if (all(grepl("^[0-9]+$", labs)) && length(labs) >= 10) {
    w <- max(nchar(labs))
    top$cluster <- formatC(as.integer(top$cluster), width = w, flag = "0")
  }

  # Rows ordered by cluster, then by rank within cluster. MarkerPlot() does
  # its own stable order(Details), so within-group order is preserved.
  top <- top[order(top$cluster,
                   if (rank_by == "avg_log2FC") -top$avg_log2FC
                   else top$p_val_adj), , drop = FALSE]

  n_per <- table(top$cluster)
  if (any(n_per < n)) {
    short <- names(n_per)[n_per < n]
    message(sprintf("  %d cluster(s) yielded fewer than %d genes: %s",
                    length(short), n, paste(short, collapse = ", ")))
  }
  message(sprintf("--- Plotting %d genes across %d cluster(s) ---",
                  nrow(top), length(n_per)))

  # ---- Plot ---------------------------------------------------------------
  genes_df <- data.frame(Gene    = top$gene,
                         Details = top$cluster,
                         stringsAsFactors = FALSE)

  p <- MarkerPlot(obj, genes_df, cluster = cluster, ...)
  attr(p, "markers")     <- markers
  attr(p, "top_markers") <- top
  p
}
