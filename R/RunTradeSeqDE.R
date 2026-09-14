#' Differential expression along pseudotime with tradeSeq
#'
#' Fits a negative-binomial GAM per gene along pseudotime
#' (\code{tradeSeq::fitGAM}) and tests it, closing the loop between fitting
#' a trajectory (\code{\link{PseudotimeWrapper}}) and asking which genes
#' actually change along it or differ between branches. Consumes the
#' standardized \code{list(pseudotime, weights)} structure
#' \code{\link{PseudotimeWrapper}} writes to \code{obj@misc$pseudotime}, so
#' it works the same way regardless of which trajectory method produced it.
#'
#' @param obj A Seurat object that has already been run through
#'   \code{\link{PseudotimeWrapper}} (or has an equivalent structure passed
#'   via \code{pseudotime}).
#' @param assay Assay to read counts from. Default \code{DefaultAssay(obj)}.
#' @param pseudotime A list with \code{pseudotime} and \code{weights}
#'   matrices (cells x lineages), as written by
#'   \code{\link{PseudotimeWrapper}}. \code{NULL} (default) reads
#'   \code{obj@misc$pseudotime}.
#' @param genes.use Genes to fit. \code{NULL} (default) uses
#'   \code{Seurat::VariableFeatures(obj)} if set, otherwise every gene in
#'   \code{assay} (can be slow for a full transcriptome -- a warning is
#'   issued in that case).
#' @param n_knots Number of knots for the GAM. Default 6, tradeSeq's own
#'   recommended starting point (see \code{tradeSeq::evaluateK} to tune).
#' @param test Which test(s) to run after fitting.
#'   \code{"associationTest"} (default) tests whether each gene changes at
#'   all along pseudotime; \code{"startVsEndTest"} compares expression at
#'   the start vs. end of each lineage; \code{"diffEndTest"} compares end
#'   points across lineages (only meaningful with >1 lineage). Pass a
#'   vector to run more than one.
#' @param parallel Logical; use \code{BiocParallel} to fit genes in
#'   parallel. Default \code{FALSE}.
#' @param workers Workers for \code{parallel = TRUE}. Default 1.
#' @param verbose Passed to \code{tradeSeq::fitGAM}. Default \code{TRUE}.
#' @return A list with \code{sce} (the fitted \code{SingleCellExperiment}
#'   from \code{fitGAM}, for \code{tradeSeq::plotSmoothers} etc.) and
#'   \code{results}, a named list of per-gene result data frames (with a
#'   \code{gene} column and \code{padj} added via FDR correction), one per
#'   entry in \code{test}.
#' @examples
#' \dontrun{
#' obj <- PseudotimeWrapper(obj, method = "slingshot", start_cluster = "3")
#' de <- RunTradeSeqDE(obj, genes.use = VariableFeatures(obj)[1:500])
#' head(de$results$associationTest)
#' tradeSeq::plotSmoothers(de$sce, counts(de$sce), gene = "Sox9")
#' }
#' @importFrom Seurat DefaultAssay GetAssayData VariableFeatures
#' @export
RunTradeSeqDE <- function(obj,
                          assay      = NULL,
                          pseudotime = NULL,
                          genes.use  = NULL,
                          n_knots    = 6,
                          test       = c("associationTest", "startVsEndTest", "diffEndTest"),
                          parallel   = FALSE,
                          workers    = 1,
                          verbose    = TRUE) {

  .assert_seurat(obj)
  if (!requireNamespace("tradeSeq", quietly = TRUE)) {
    stop("'tradeSeq' is required. Install with ",
         "BiocManager::install('tradeSeq').")
  }

  pt_info <- if (is.null(pseudotime)) obj@misc$pseudotime else pseudotime
  if (is.null(pt_info) || is.null(pt_info$pseudotime) || is.null(pt_info$weights)) {
    stop("No pseudotime structure found. Run PseudotimeWrapper(obj, ...) ",
         "first, or pass `pseudotime = list(pseudotime = ..., weights = ...)` ",
         "directly.")
  }

  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay
  genes <- genes.use
  if (is.null(genes)) {
    genes <- tryCatch(Seurat::VariableFeatures(obj), error = function(e) character(0))
    if (length(genes) == 0) {
      genes <- rownames(obj[[a]])
      warning("`genes.use` not set and no VariableFeatures() found -- ",
              "fitting all ", length(genes), " genes. This can be slow; ",
              "consider passing a smaller `genes.use`.")
    }
  }

  counts <- as.matrix(Seurat::GetAssayData(obj, assay = a, layer = "counts")[genes, , drop = FALSE])
  pt <- pt_info$pseudotime[colnames(counts), , drop = FALSE]
  wt <- pt_info$weights[colnames(counts), , drop = FALSE]
  pt[is.na(pt)] <- 0  # fitGAM expects numeric, not NA, where weight is 0

  BPPARAM <- if (isTRUE(parallel)) {
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      stop("'BiocParallel' is required for parallel = TRUE. Install with ",
           "BiocManager::install('BiocParallel').")
    }
    BiocParallel::MulticoreParam(workers)
  } else {
    if (requireNamespace("BiocParallel", quietly = TRUE)) {
      BiocParallel::SerialParam()
    } else {
      NULL
    }
  }

  if (isTRUE(verbose)) {
    message(sprintf("--- Fitting GAMs (%d genes, %d lineage(s), %d knots) ---",
                    nrow(counts), ncol(pt), n_knots))
  }
  fit_args <- list(counts = counts, pseudotime = pt, cellWeights = wt,
                   nknots = n_knots, verbose = isTRUE(verbose))
  if (!is.null(BPPARAM)) fit_args$parallel <- isTRUE(parallel)
  if (!is.null(BPPARAM)) fit_args$BPPARAM  <- BPPARAM
  sce <- do.call(tradeSeq::fitGAM, fit_args)

  test_fun <- list(
    associationTest = tradeSeq::associationTest,
    startVsEndTest  = tradeSeq::startVsEndTest,
    diffEndTest     = tradeSeq::diffEndTest
  )
  results <- list()
  for (t in test) {
    if (isTRUE(verbose)) message(sprintf("--- Running %s ---", t))
    res <- as.data.frame(test_fun[[t]](sce))
    res$gene <- rownames(res)
    p_col <- intersect(c("pvalue", "pvalue_lineage1"), colnames(res))[1]
    if (!is.na(p_col)) res$padj <- stats::p.adjust(res[[p_col]], method = "fdr")
    res <- res[, c("gene", setdiff(colnames(res), "gene"))]
    results[[t]] <- res
  }

  list(sce = sce, results = results)
}
