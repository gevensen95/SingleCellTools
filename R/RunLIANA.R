#' Ligand-receptor analysis with LIANA
#'
#' Wraps \code{liana::liana_wrap} + \code{liana::liana_aggregate} to run
#' one or more ligand-receptor inference methods (CellPhoneDB, NATMI,
#' Connectome, SingleCellSignalR, logFC, etc.) and produce a single
#' aggregated / consensus-scored results table with a standardized
#' schema.
#'
#' By default runs LIANA's consensus method set and aggregates their
#' rankings, which is more robust than any single method. Pass
#' \code{method = "cellphonedb"} (or similar) to run a single method.
#'
#' \strong{Species and resources.} LIANA pulls its ligand-receptor
#' resource from OmnipathR. Pass \code{resource = "Consensus"} (default)
#' for human OR use \code{liana::show_resources()} to explore. For mouse,
#' set \code{use_ortho = TRUE} to have LIANA translate the human resource
#' to mouse gene symbols. Genes are matched case-sensitively — check that
#' \code{rownames(obj[[assay]])} match the resource's naming convention.
#'
#' @param obj A Seurat object with clusters / cell-type labels ready to
#'   set as the \code{Idents()} (via \code{idents_col}). Uses the
#'   \code{"data"} layer of \code{assay} for expression.
#' @param idents_col Metadata column to use as the source/target cell-type
#'   identity. Default \code{"seurat_clusters"}.
#' @param assay Assay to read from. Default DefaultAssay(obj).
#' @param method LIANA method(s) to run. \code{"consensus"} (default) runs
#'   LIANA's default consensus set (\code{sca}, \code{natmi},
#'   \code{connectome}, \code{logfc}, \code{cellphonedb}) and aggregates.
#'   Pass a character vector to run a specific subset.
#' @param resource Resource name for LR pairs, passed to
#'   \code{liana_wrap(resource = ...)}. Default \code{"Consensus"}.
#' @param use_ortho If TRUE, translate the human resource to mouse
#'   orthologs via OmnipathR. Default FALSE.
#' @param source_cells Character vector of \code{idents_col} values to use
#'   as sources (ligand-expressing). \code{NULL} (default) uses all.
#' @param target_cells Same as \code{source_cells} but for targets
#'   (receptor-expressing).
#' @param min_cells Minimum cells per (cluster, method) required to score
#'   a cell type. LIANA's default is 5. Default 10.
#' @param verbose Passed to liana. Default FALSE.
#' @param return_raw If TRUE, returns the full aggregated LIANA object
#'   instead of the flattened data frame. Default FALSE.
#' @return A data frame with one row per (source, target, ligand, receptor)
#'   interaction, sorted by aggregated rank. Columns include \code{source},
#'   \code{target}, \code{ligand.complex}, \code{receptor.complex},
#'   \code{aggregate_rank}, and the per-method scores that were run.
#' @examples
#' \dontrun{
#' Idents(obj) <- obj$cell_type
#' lr <- RunLIANA(obj, idents_col = "cell_type")
#' head(lr)
#'
#' # Just interactions FROM T cells TO everything
#' lr_t <- RunLIANA(obj, idents_col = "cell_type",
#'                  source_cells = "T cell")
#'
#' # Single method, mouse
#' lr_cpdb <- RunLIANA(obj, idents_col = "cell_type",
#'                     method    = "cellphonedb",
#'                     use_ortho = TRUE)
#' }
#' @importFrom Seurat DefaultAssay Idents Idents<-
#' @export
RunLIANA <- function(obj,
                     idents_col   = "seurat_clusters",
                     assay        = NULL,
                     method       = "consensus",
                     resource     = "Consensus",
                     use_ortho    = FALSE,
                     source_cells = NULL,
                     target_cells = NULL,
                     min_cells    = 10,
                     verbose      = FALSE,
                     return_raw   = FALSE) {

  if (!requireNamespace("liana", quietly = TRUE)) {
    stop("'liana' is required. Install with ",
         "remotes::install_github('saezlab/liana').")
  }
  if (!inherits(obj, "Seurat")) stop("`obj` must be a Seurat object.")
  if (!idents_col %in% colnames(obj@meta.data)) {
    stop("Idents column '", idents_col, "' not found in obj@meta.data.")
  }

  # Set Idents from the chosen column
  Seurat::Idents(obj) <- as.factor(as.character(obj@meta.data[[idents_col]]))
  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay

  # LIANA v0.1.x uses method="consensus" as a special string; newer
  # versions expect a vector via `method =`. Handle both.
  method_arg <- if (identical(method, "consensus")) {
    c("connectome", "logfc", "natmi", "sca", "cellphonedb")
  } else {
    method
  }

  message(sprintf(
    "--- Running LIANA (%d method(s), resource '%s'%s) ---",
    length(method_arg), resource,
    if (use_ortho) ", ortholog-translated" else ""))

  liana_res <- liana::liana_wrap(
    sce         = obj,
    method      = method_arg,
    resource    = resource,
    assay       = a,
    idents_col  = "seurat_clusters",  # LIANA reads Idents() directly
    min_cells   = min_cells,
    verbose     = isTRUE(verbose),
    # ortholog handling — only applies when the resource is human-origin
    resource_orthologs = if (use_ortho) TRUE else FALSE
  )

  agg <- liana::liana_aggregate(liana_res, verbose = isTRUE(verbose))
  if (isTRUE(return_raw)) return(agg)

  # ---- Flatten to a tidy data frame ---------------------------------------
  df <- as.data.frame(agg)

  # Restrict to source / target if requested
  if (!is.null(source_cells)) {
    df <- df[as.character(df$source) %in% source_cells, , drop = FALSE]
  }
  if (!is.null(target_cells)) {
    df <- df[as.character(df$target) %in% target_cells, , drop = FALSE]
  }

  # Sort by aggregate_rank if present
  if ("aggregate_rank" %in% colnames(df)) {
    df <- df[order(df$aggregate_rank), ]
  }
  rownames(df) <- NULL
  df
}
