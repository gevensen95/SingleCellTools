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
#' set \code{use_ortho = TRUE}: with the default \code{resource =
#' "Consensus"} this switches to LIANA's built-in, already-mouse-symbol
#' \code{"MouseConsensus"} resource; for any other named resource it
#' translates that human resource to mouse gene symbols via
#' \code{liana::generate_homologs()}. Genes are matched case-sensitively
#' — check that \code{rownames(obj[[assay]])} match the resource's naming
#' convention (mouse symbols are Title Case, e.g. \code{"Tnf"}, not
#' \code{"TNF"}).
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
#'   \code{aggregate_rank} occasionally can't be computed (see the
#'   \code{get_agrank} fallback note in the source) -- in that case the
#'   data frame is instead sorted by \code{mean_rank}, and
#'   \code{aggregate_rank} is simply absent.
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
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("'SingleCellExperiment' is required. Install with ",
         "BiocManager::install('SingleCellExperiment').")
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

  # ---- Convert Seurat -> SingleCellExperiment ourselves ------------------
  # liana::liana_wrap() dispatches on `sce`'s class via liana::liana_prep().
  # As of this writing, liana::liana_prep.Seurat() still calls
  # SeuratObject::GetAssayData(object, assay, slot = "counts"/"data"), and
  # the `slot=` argument was made fully defunct in SeuratObject 5.0.0 — so
  # handing liana_wrap() a Seurat object directly hard-errors on any
  # SeuratObject >= 5.0.0 install. This is a real, currently open upstream
  # bug (https://github.com/saezlab/liana/issues/194), not something on our
  # side to "wait out". liana::liana_prep.SingleCellExperiment() doesn't
  # touch GetAssayData() at all, so we sidestep the broken code path by
  # doing the Seurat -> SCE conversion here ourselves with the modern
  # `layer=` argument, and pass liana_wrap() an SCE instead of the raw
  # Seurat object. Remove this workaround once liana ships a real fix.
  sce <- SingleCellExperiment::SingleCellExperiment(
    list(
      counts    = Seurat::GetAssayData(obj, assay = a, layer = "counts"),
      logcounts = Seurat::GetAssayData(obj, assay = a, layer = "data")
    ),
    colData = obj@meta.data
  )
  SingleCellExperiment::colLabels(sce) <- Seurat::Idents(obj)
  # Seurat's "data" layer is natural-log-normalized by default (unlike SCE's
  # log2 default) -- match what liana_prep.Seurat() would have set so the
  # logFC method's base-conversion stays correct.
  sce@int_metadata$base <- exp(1)

  # ---- use_ortho: translate the resource to mouse gene symbols -----------
  # `resource_orthologs` (passed here previously) is not a real parameter
  # of liana::liana_wrap() -- it doesn't appear in liana_wrap()'s formals
  # or in any of the liana_defaults() sub-lists its `...` forwards to, so
  # passing it was silently absorbed and did *nothing*: use_ortho = TRUE
  # was matching the human resource against mouse gene symbols exactly as
  # typed, case-sensitively, the whole time. Confirmed empirically: even
  # with use_ortho = TRUE, every method returned 0 rows against a mouse
  # object using the human "Consensus" resource (the "no non-missing
  # arguments to min/max" warnings liana_aggregate() throws afterward are
  # that all-zero-rows condition propagating through its own code -- not a
  # separate bug in the aggregate_rank fallback above).
  #
  # liana ships a curated "MouseConsensus" resource already in mouse gene
  # symbols (see liana::show_resources()) -- used directly whenever the
  # (default) "Consensus" resource is requested, no translation needed.
  # For any other named resource, fall back to the actual (exported)
  # liana::generate_homologs(), which translates a human OmniPath resource
  # to another organism's symbols via OmnipathR's HomoloGene mapping.
  # .missing_fun = stringr::str_to_title is liana's own documented
  # suggestion for murine data: genes with no exact HomoloGene match get
  # Title-Cased from their human symbol instead of being dropped, trading
  # a little mismatch risk for substantially better coverage.
  resource_arg <- resource
  external_arg <- NULL
  if (isTRUE(use_ortho)) {
    if (identical(resource, "Consensus") &&
        "MouseConsensus" %in% liana::show_resources()) {
      resource_arg <- "MouseConsensus"
    } else {
      if (!requireNamespace("stringr", quietly = TRUE)) {
        stop("'stringr' is required for use_ortho = TRUE with a resource ",
             "other than 'Consensus'.")
      }
      message("--- Translating '", resource, "' resource to mouse gene ",
              "symbols via liana::generate_homologs() ---")
      human_resource <- liana::select_resource(resource)[[1]]
      external_arg   <- liana::generate_homologs(
        human_resource,
        target_organism = 10090L,  # NCBI taxid for Mus musculus
        .missing_fun    = stringr::str_to_title,
        verbose         = isTRUE(verbose)
      )
      resource_arg <- "custom"
    }
  }

  liana_res <- liana::liana_wrap(
    sce               = sce,
    method            = method_arg,
    resource          = resource_arg,
    external_resource = external_arg,
    idents_col        = NULL,  # NULL -> liana matches colLabels(sce), set above
    min_cells         = min_cells,
    verbose           = isTRUE(verbose)
  )

  # liana::liana_aggregate()'s default get_agrank = TRUE runs a hand-rolled
  # RobustRankAggreg reimplementation (liana's own comments say this is a
  # stopgap since RobustRankAggreg was removed from CRAN) that unites
  # source/target/ligand/receptor into a single "interaction" string,
  # pivots it to rownames, computes ranks, then splits "interaction" back
  # apart via separate(). If that intermediate matrix ends up with no
  # rownames to split -- e.g. very few surviving interactions, which is
  # exactly what heavy per-cell-type downsampling upstream produces -- it
  # fails with "Column `interaction` doesn't exist", aborting the whole
  # call even though every method's actual results computed fine. Try the
  # real default first (it's a legitimate, more informative ranking when
  # it works), and fall back to get_agrank = FALSE -- which just skips the
  # `aggregate_rank` column and keeps the plain mean-of-ranks `mean_rank`
  # column from get_ranks = TRUE (still on by default, unaffected by this
  # bug) -- rather than losing the whole result over one fragile column.
  agg <- tryCatch(
    liana::liana_aggregate(liana_res, verbose = isTRUE(verbose)),
    error = function(e) {
      message("liana::liana_aggregate()'s RobustRankAggreg-based ",
              "`aggregate_rank` step failed (", conditionMessage(e),
              "); falling back to get_agrank = FALSE. Results below are ",
              "sorted by `mean_rank` instead of `aggregate_rank`.")
      liana::liana_aggregate(liana_res, verbose = isTRUE(verbose),
                             get_agrank = FALSE)
    }
  )
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

  # Sort by aggregate_rank if present (falls back to mean_rank when the
  # RobustRankAggreg step above was skipped -- see the tryCatch note).
  if ("aggregate_rank" %in% colnames(df)) {
    df <- df[order(df$aggregate_rank), ]
  } else if ("mean_rank" %in% colnames(df)) {
    df <- df[order(df$mean_rank), ]
  }
  rownames(df) <- NULL
  df
}
