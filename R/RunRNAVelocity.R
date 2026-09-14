#' RNA velocity from spliced/unspliced counts
#'
#' Estimates per-cell transcriptional direction ("RNA velocity") from
#' spliced/unspliced read counts, via \code{scVelo} (Python, through
#' \code{reticulate}) or \code{velocyto.R} (pure R). Input is typically a
#' \code{.loom} file from \code{velocyto}, \code{kb-python}, or STARsolo's
#' \code{--soloFeatures Velocyto} output; a matched pair of spliced/
#' unspliced matrices works too.
#'
#' \strong{Choosing a method.} \code{"scvelo"} (default) is the actively
#' maintained, more complete option -- it supports \code{mode =
#' "dynamical"} (EM-based, gives per-cell latent time, the mode scVelo's
#' authors recommend, but slower) as well as the cheaper legacy
#' \code{"stochastic"}/\code{"deterministic"} modes -- but requires a
#' working Python environment with \code{scvelo} installed
#' (\code{reticulate::py_install("scvelo")} or \code{pip install scvelo}).
#' \code{"velocyto.R"} needs no Python at all, at the cost of being older
#' and less actively maintained, and its streamline plot
#' (\code{velocyto.R::show.velocity.on.embedding.cor}) is a base-R plot,
#' not a \code{ggplot} -- call it yourself on \code{obj@misc$velocity} for
#' that visualization; this function does not attempt to wrap it as one.
#'
#' scVelo's own plotting (streamline/grid velocity plots) is Python-side
#' and not reproduced here -- \code{obj@misc$velocity$h5ad} points at a
#' written-out \code{.h5ad} file with the full result (including the
#' velocity graph) for use with \code{scv.pl.velocity_embedding_stream()}
#' in Python, or via \code{reticulate} directly.
#'
#' @param obj A Seurat object whose cells match the spliced/unspliced input.
#' @param loom_file Path to a \code{.loom} file containing spliced/
#'   unspliced (and typically ambiguous) layers. Either this or
#'   \code{spliced}/\code{unspliced} is required.
#' @param spliced,unspliced Gene x cell spliced / unspliced count matrices,
#'   as an alternative to \code{loom_file}. Column names must match
#'   \code{colnames(obj)}.
#' @param method \code{"scvelo"} (default) or \code{"velocyto.R"}.
#' @param mode \code{method = "scvelo"} only: \code{"dynamical"} (default),
#'   \code{"stochastic"}, or \code{"deterministic"}. See Details.
#' @param reduction Embedding to attach the velocity graph to / compute
#'   cell-cell distances from. Default \code{"umap"}.
#' @param n_top_genes \code{method = "scvelo"} only: genes kept after
#'   \code{scv.pp.filter_and_normalize}. Default 2000.
#' @param h5ad_out \code{method = "scvelo"} only: path to write the full
#'   AnnData result to. Default a fresh \code{tempfile(fileext = ".h5ad")}.
#' @param k_cells \code{method = "velocyto.R"} only: neighborhood size for
#'   \code{gene.relative.velocity.estimates}. Default 25.
#' @param verbose Message progress. Default \code{TRUE}.
#' @return \code{obj} with \code{obj@misc$velocity} populated:
#'   for \code{"scvelo"}, \code{list(method, mode, h5ad, latent_time)}
#'   (\code{latent_time} only in dynamical mode, also written as a
#'   \code{velocity_latent_time} metadata column); for \code{"velocyto.R"},
#'   the raw \code{gene.relative.velocity.estimates()} result.
#' @examples
#' \dontrun{
#' obj <- RunRNAVelocity(obj, loom_file = "sample.loom", mode = "dynamical")
#' FeaturePlot(obj, features = "velocity_latent_time")
#' # obj@misc$velocity$h5ad -- open in Python for scv.pl.velocity_embedding_stream()
#'
#' # No Python available
#' obj <- RunRNAVelocity(obj, loom_file = "sample.loom", method = "velocyto.R")
#' }
#' @importFrom Seurat Embeddings
#' @export
RunRNAVelocity <- function(obj,
                           loom_file   = NULL,
                           spliced     = NULL,
                           unspliced   = NULL,
                           method      = c("scvelo", "velocyto.R"),
                           mode        = c("dynamical", "stochastic", "deterministic"),
                           reduction   = "umap",
                           n_top_genes = 2000,
                           h5ad_out    = tempfile(fileext = ".h5ad"),
                           k_cells     = 25,
                           verbose     = TRUE) {

  method <- match.arg(method)
  mode   <- match.arg(mode)
  .assert_seurat(obj)
  if (is.null(loom_file) && (is.null(spliced) || is.null(unspliced))) {
    stop("Provide either `loom_file`, or both `spliced` and `unspliced`.")
  }
  if (!(reduction %in% names(obj@reductions))) {
    stop("Reduction '", reduction, "' not found.")
  }

  if (method == "scvelo") {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      stop("'reticulate' is required for method = 'scvelo'.")
    }
    scv <- tryCatch(reticulate::import("scvelo"), error = function(e) NULL)
    if (is.null(scv)) {
      stop("Python package 'scvelo' not found in the active reticulate ",
           "environment. Install with reticulate::py_install('scvelo') ",
           "or pip install scvelo, or use method = 'velocyto.R' instead.")
    }
    ad <- reticulate::import("anndata")

    if (!is.null(loom_file)) {
      if (isTRUE(verbose)) message(sprintf("--- Reading %s ---", loom_file))
      adata <- scv$read_loom(loom_file)
      common <- intersect(adata$obs_names$values, colnames(obj))
      adata <- adata[common, ]
    } else {
      common <- intersect(intersect(colnames(spliced), colnames(unspliced)), colnames(obj))
      if (length(common) == 0) stop("No common cells between spliced/unspliced and `obj`.")
      s <- t(as.matrix(spliced[, common, drop = FALSE]))
      u <- t(as.matrix(unspliced[, common, drop = FALSE]))
      adata <- ad$AnnData(X = s, layers = list(spliced = s, unspliced = u))
      adata$obs_names <- common
      adata$var_names <- colnames(s)
    }
    emb <- Seurat::Embeddings(obj, reduction = reduction)[reticulate::py_to_r(adata$obs_names$values), , drop = FALSE]
    adata$obsm[["X_umap"]] <- emb

    if (isTRUE(verbose)) message("--- scv.pp.filter_and_normalize / moments ---")
    scv$pp$filter_and_normalize(adata, min_shared_counts = 20L,
                                n_top_genes = as.integer(n_top_genes))
    scv$pp$moments(adata, n_pcs = 30L, n_neighbors = 30L)

    if (isTRUE(verbose)) message(sprintf("--- scv.tl.velocity (mode = %s) ---", mode))
    if (mode == "dynamical") scv$tl$recover_dynamics(adata, verbose = isTRUE(verbose))
    scv$tl$velocity(adata, mode = mode)
    scv$tl$velocity_graph(adata)

    latent_time <- NULL
    if (mode == "dynamical") {
      scv$tl$velocity_pseudotime(adata)
      scv$tl$recover_latent_time(adata)
      latent_time <- setNames(reticulate::py_to_r(adata$obs[["latent_time"]]$values),
                              reticulate::py_to_r(adata$obs_names$values))
      obj$velocity_latent_time <- latent_time[colnames(obj)]
    }

    adata$write(h5ad_out)
    obj@misc$velocity <- list(method = "scvelo", mode = mode,
                              h5ad = h5ad_out, latent_time = latent_time)
    if (isTRUE(verbose)) message(sprintf("  Full result written to %s", h5ad_out))

  } else {
    if (!requireNamespace("velocyto.R", quietly = TRUE)) {
      stop("'velocyto.R' is required for method = 'velocyto.R'. Install ",
           "with remotes::install_github('velocyto-team/velocyto.R').")
    }
    if (!is.null(loom_file)) {
      if (isTRUE(verbose)) message(sprintf("--- Reading %s ---", loom_file))
      dat <- velocyto.R::read.loom.matrices(loom_file)
      spliced   <- dat$spliced
      unspliced <- dat$unspliced
    }
    common <- intersect(intersect(colnames(spliced), colnames(unspliced)), colnames(obj))
    if (length(common) == 0) stop("No common cells between spliced/unspliced and `obj`.")
    emb <- Seurat::Embeddings(obj, reduction = reduction)[common, , drop = FALSE]
    cell_dist <- as.dist(1 - stats::cor(t(emb)))

    if (isTRUE(verbose)) {
      message(sprintf("--- Running velocyto.R (%d cells, k = %d) ---", length(common), k_cells))
    }
    rvel <- velocyto.R::gene.relative.velocity.estimates(
      spliced[, common, drop = FALSE], unspliced[, common, drop = FALSE],
      kCells = k_cells, cell.dist = cell_dist, fit.quantile = 0.02
    )
    obj@misc$velocity <- rvel
  }

  obj
}
