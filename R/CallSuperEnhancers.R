#' Super-enhancer calling from pseudobulk ATAC signal (ROSE-style)
#'
#' Classic bulk-epigenomics super-enhancer calling (ROSE: Whyte/Loven et al.
#' 2013), adapted to per-cell-type pseudobulk ATAC signal instead of a bulk
#' H3K27ac track: peaks are pseudobulked per level of \code{cluster_col},
#' nearby peaks are stitched into candidate regions, and regions are ranked
#' by total signal to find the "hockey stick" inflection point separating
#' ordinary enhancers from a small set of unusually large, unusually active
#' stitched regions (super-enhancers) -- master regulatory loci often
#' associated with cell-identity genes.
#'
#' \strong{The cutoff here is a simplified stand-in for ROSE's own
#' tangent-line method}, not a re-implementation of it: it takes the
#' rank-normalized point that maximizes distance above the diagonal
#' (\code{signal_norm - rank_norm}), which behaves similarly for a typical
#' right-skewed signal distribution but is not guaranteed identical to
#' ROSE's slope-1-tangent construction. Treat the cutoff as a reasonable
#' default to inspect on the hockey-stick plot (\code{plot = TRUE}), not as
#' a validated statistical threshold.
#'
#' If \code{exclude_promoters = TRUE} (default) and
#' \code{\link{AnnotatePeaks}} has been run (\code{obj@misc$peak_annotation}
#' present), promoter-proximal peaks are dropped before stitching -- ROSE's
#' own convention, since promoter (not enhancer) signal at highly-expressed
#' genes would otherwise dominate the ranking.
#'
#' @param obj A Seurat object with a \code{ChromatinAssay}.
#' @param cluster_col Metadata column of cluster/cell-type labels; signal is
#'   pseudobulked (summed) separately per level. Required.
#' @param peak_assay Name of the \code{ChromatinAssay}. \code{NULL}
#'   (default) auto-detects it -- errors if \code{obj} has zero or more
#'   than one.
#' @param clusters Which \code{cluster_col} levels to run on. \code{NULL}
#'   (default) runs every level.
#' @param stitch_distance Peaks within this many bp of each other are
#'   merged into one candidate region before ranking. Default 12500 (ROSE's
#'   own default).
#' @param exclude_promoters Drop promoter-proximal peaks before stitching.
#'   Default \code{TRUE}; requires \code{\link{AnnotatePeaks}} to have been
#'   run, otherwise a warning is issued and it's treated as \code{FALSE}.
#' @param min_cells Minimum cells in a \code{cluster_col} level to attempt
#'   pseudobulking. Default 20; smaller clusters are skipped with a message.
#' @param plot Logical; if \code{TRUE} (default), also build a hockey-stick
#'   plot (rank vs. signal, super-enhancers highlighted) per cluster.
#' @param verbose Message progress. Default \code{TRUE}.
#' @return \code{obj} with \code{obj@misc$super_enhancers} populated: a
#'   named list (one element per cluster) of
#'   \code{list(regions, table, cutoff_rank, plot)} -- \code{regions} is a
#'   \code{GRanges} of stitched candidate regions with a \code{signal}
#'   column; \code{table} is the same as a ranked data frame with
#'   \code{is_super_enhancer}.
#' @examples
#' \dontrun{
#' atac <- AnnotatePeaks(atac)
#' atac <- CallSuperEnhancers(atac, cluster_col = "cell_type")
#' atac@misc$super_enhancers[["Hepatocyte"]]$plot
#' se <- atac@misc$super_enhancers[["Hepatocyte"]]$table
#' subset(se, is_super_enhancer)
#' }
#' @importFrom Seurat DefaultAssay GetAssayData
#' @importFrom ggplot2 aes geom_point geom_vline ggplot labs scale_color_manual
#' @export
CallSuperEnhancers <- function(obj,
                               cluster_col,
                               peak_assay        = NULL,
                               clusters          = NULL,
                               stitch_distance   = 12500,
                               exclude_promoters = TRUE,
                               min_cells         = 20,
                               plot              = TRUE,
                               verbose           = TRUE) {

  .assert_seurat(obj)
  if (!cluster_col %in% colnames(obj@meta.data)) {
    stop("`cluster_col` '", cluster_col, "' not found in obj@meta.data.")
  }

  pa <- peak_assay
  if (is.null(pa)) {
    is_chromatin <- vapply(Seurat::Assays(obj),
                           function(a) inherits(obj[[a]], "ChromatinAssay"),
                           logical(1))
    candidates <- Seurat::Assays(obj)[is_chromatin]
    if (length(candidates) != 1) {
      stop("Could not auto-detect a single ChromatinAssay in `obj` ",
           "(found: ", paste(candidates, collapse = ", "),
           "). Pass `peak_assay` explicitly.")
    }
    pa <- candidates
  } else if (!inherits(obj[[pa]], "ChromatinAssay")) {
    stop("`peak_assay` ('", pa, "') is not a ChromatinAssay.")
  }

  peaks_gr <- Signac::granges(obj[[pa]])
  peak_names <- rownames(obj[[pa]])
  names(peaks_gr) <- peak_names

  keep_peaks <- peak_names
  if (isTRUE(exclude_promoters)) {
    ann <- obj@misc$peak_annotation
    if (is.null(ann)) {
      warning("`exclude_promoters = TRUE` but obj@misc$peak_annotation not ",
             "found -- run AnnotatePeaks(obj) first. Proceeding without ",
             "promoter exclusion.")
    } else {
      keep_peaks <- ann$peak[is.na(ann$peak_type) | ann$peak_type != "promoter_proximal"]
    }
  }
  peaks_gr <- peaks_gr[keep_peaks]

  counts <- Seurat::GetAssayData(obj, assay = pa, layer = "counts")
  cl_levels <- if (is.null(clusters)) {
    sort(unique(as.character(obj@meta.data[[cluster_col]])))
  } else as.character(clusters)

  out <- list()
  for (cl in cl_levels) {
    cells <- colnames(obj)[as.character(obj@meta.data[[cluster_col]]) == cl]
    if (length(cells) < min_cells) {
      if (isTRUE(verbose)) {
        message(sprintf("  Skipping '%s': %d cell(s) < min_cells (%d)",
                        cl, length(cells), min_cells))
      }
      next
    }
    if (isTRUE(verbose)) {
      message(sprintf("--- Cluster '%s': %d cells, %d peaks (stitch <= %d bp) ---",
                      cl, length(cells), length(peaks_gr), stitch_distance))
    }
    signal <- Matrix::rowSums(counts[keep_peaks, cells, drop = FALSE])
    gr <- peaks_gr
    gr$signal <- signal[names(gr)]
    gr <- sort(gr)

    stitched <- GenomicRanges::reduce(gr, min.gapwidth = stitch_distance)
    ov <- GenomicRanges::findOverlaps(gr, stitched)
    sig_by_region <- tapply(gr$signal[S4Vectors::queryHits(ov)],
                            S4Vectors::subjectHits(ov), sum)
    stitched$signal <- as.numeric(sig_by_region[as.character(seq_along(stitched))])
    stitched$signal[is.na(stitched$signal)] <- 0

    df <- data.frame(region = as.character(stitched), signal = stitched$signal,
                     stringsAsFactors = FALSE)
    df <- df[order(df$signal), ]
    df$rank <- seq_len(nrow(df))

    x <- df$rank / max(df$rank)
    y <- if (max(df$signal) > 0) df$signal / max(df$signal) else rep(0, nrow(df))
    idx_cut <- which.max(y - x)
    cutoff_rank <- df$rank[idx_cut]
    df$is_super_enhancer <- df$rank > cutoff_rank
    rownames(df) <- NULL
    stitched$is_super_enhancer <- stitched$signal %in% df$signal[df$is_super_enhancer]

    p <- NULL
    if (isTRUE(plot)) {
      is_super_enhancer <- rank <- signal <- NULL  # NSE
      p <- ggplot2::ggplot(df, ggplot2::aes(x = rank, y = signal,
                                            color = is_super_enhancer)) +
        ggplot2::geom_vline(xintercept = cutoff_rank, linetype = 2, color = "grey50") +
        ggplot2::geom_point(size = 1.2) +
        ggplot2::scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "#B2182B"),
                                    name = "super-enhancer") +
        ggplot2::labs(x = "stitched region rank", y = "pseudobulk signal",
                     title = sprintf("Super-enhancers: %s (%d/%d regions)",
                                    cl, sum(df$is_super_enhancer), nrow(df))) +
        Ol_Reliable()
    }

    if (isTRUE(verbose)) {
      message(sprintf("  %d/%d stitched region(s) called super-enhancers.",
                      sum(df$is_super_enhancer), nrow(df)))
    }
    out[[cl]] <- list(regions = stitched, table = df,
                      cutoff_rank = cutoff_rank, plot = p)
  }

  if (length(out) == 0) {
    stop("No cluster had >= min_cells (", min_cells, ") cells to pseudobulk.")
  }
  obj@misc$super_enhancers <- out
  obj
}
