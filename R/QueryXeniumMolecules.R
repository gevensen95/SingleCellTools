#' Query Xenium Molecule Positions From a Lazy Arrow Dataset
#'
#' Windowed / gene-subset accessor for Seurat objects built with
#' \code{\link{LoadXenium2}(..., microns_lazy = TRUE)}. Reads only the rows
#' that match the requested filters straight from \code{transcripts.parquet}
#' via \code{arrow}'s query engine, instead of re-reading (or having kept in
#' memory) the full transcript table -- for a whole-slide Xenium run this is
#' the difference between touching 10^8+ rows and touching only the ones you
#' actually asked for.
#'
#' @param obj A Seurat object built with
#'   \code{LoadXenium2(..., microns_lazy = TRUE)}.
#' @param genes Optional character vector of gene/feature names to restrict
#'   to. \code{NULL} (default) doesn't restrict by gene.
#' @param x_range,y_range Optional length-2 numeric vectors \code{c(min,
#'   max)} giving a bounding-box window, in the same pixel coordinate space
#'   as \code{x_location}/\code{y_location} in \code{transcripts.parquet}.
#'   \code{NULL} (default) doesn't restrict that axis.
#' @param qv_threshold Minimum QV to keep. Default \code{20} (10x's own
#'   recommended cutoff) -- pass the same threshold \code{obj} was
#'   originally loaded with if you want consistency with the molecules
#'   already attached to \code{obj[["fov"]]}.
#' @return A data frame with columns \code{x}, \code{y}, \code{gene} for the
#'   matching transcripts.
#' @export
QueryXeniumMolecules <- function(obj, genes = NULL, x_range = NULL,
                                 y_range = NULL, qv_threshold = 20) {
  qv <- feature_name <- x_location <- y_location <- NULL  # silence R CMD check NSE notes

  if (!inherits(obj, "Seurat")) {
    stop("`obj` must be a Seurat object.")
  }
  if (!is.null(genes) && !is.character(genes)) {
    stop("`genes` must be a character vector, or NULL.")
  }
  if (!is.null(x_range) && (!is.numeric(x_range) || length(x_range) != 2L)) {
    stop("`x_range` must be a numeric vector of length 2 (min, max), or NULL.")
  }
  if (!is.null(y_range) && (!is.numeric(y_range) || length(y_range) != 2L)) {
    stop("`y_range` must be a numeric vector of length 2 (min, max), or NULL.")
  }
  if (!is.numeric(qv_threshold) || length(qv_threshold) != 1L) {
    stop("`qv_threshold` must be a single number.")
  }

  ds <- obj@misc$molecules_lazy
  if (is.null(ds)) {
    stop("`obj@misc$molecules_lazy` not found. Rebuild with ",
         "LoadXenium2(..., microns_lazy = TRUE) to enable windowed/",
         "gene-subset queries without re-reading transcripts.parquet.")
  }

  message(sprintf('--- Querying transcripts.parquet (qv >= %g%s%s%s) ---',
                  qv_threshold,
                  if (!is.null(genes)) sprintf(', %d gene(s)', length(genes)) else '',
                  if (!is.null(x_range)) ', x-windowed' else '',
                  if (!is.null(y_range)) ', y-windowed' else ''))

  q <- dplyr::filter(ds, qv >= qv_threshold)
  if (!is.null(genes))   q <- dplyr::filter(q, feature_name %in% genes)
  if (!is.null(x_range)) q <- dplyr::filter(q, x_location >= x_range[1], x_location <= x_range[2])
  if (!is.null(y_range)) q <- dplyr::filter(q, y_location >= y_range[1], y_location <= y_range[2])
  q <- dplyr::select(q, x = x_location, y = y_location, gene = feature_name)

  as.data.frame(dplyr::collect(q))
}
