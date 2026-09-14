#' RNA-only regulon inference (SCENIC: GENIE3/GRNBoost2 + RcisTarget + AUCell)
#'
#' Infers transcription-factor regulons from expression data alone, the
#' classic SCENIC workflow: (1) co-expression modules between candidate
#' regulators (\code{tf_list}) and every gene, via \code{GENIE3} (R) or
#' \code{GRNBoost2} (Python, through \code{reticulate}); (2) motif-based
#' pruning of each module down to a "regulon" -- only targets whose
#' upstream region is independently enriched for the TF's own motif, via
#' \code{RcisTarget} -- so a regulon means "co-expressed AND has direct
#' motif support," not just correlated; (3) per-cell regulon activity
#' scoring via \code{AUCell} (reusing the same scorer as
#' \code{\link{RunSingleCellGSEA}}).
#'
#' \strong{Why \code{tf_list} is required, with no bundled default.} Unlike
#' \code{\link{FetchGeneSets}}'s curated collections, this package does not
#' ship a species-specific transcription-factor list -- silently defaulting
#' to a human or mouse TF list would misclassify runs on other species
#' exactly the way un-checked gene-symbol casing did (see
#' \code{\link{FetchGeneSets}}'s species handling), except worse, since a
#' wrong TF list doesn't just mismatch gene sets, it changes which genes are
#' even considered as regulators. Get one from a source appropriate to your
#' species (e.g. AnimalTFDB, a curated GO:0003700 "DNA-binding transcription
#' factor activity" gene list, or \code{RcisTarget}'s own bundled TF lists
#' for human/mouse) and pass it explicitly.
#'
#' \strong{Motif pruning is optional but strongly recommended.} If
#' \code{motif_rankings}/\code{motif_annotations} are omitted, this function
#' skips the \code{RcisTarget} step entirely and returns the raw GENIE3/
#' GRNBoost2 co-expression modules as "regulons" -- these are \emph{not}
#' SCENIC regulons in the usual sense (no motif support required), just
#' TF-anchored co-expression modules, and are labeled as such in the
#' returned object and console messages. Real \code{RcisTarget} motif
#' ranking databases are large external \code{.feather} files (per-species,
#' from \url{https://resources.aertslab.org/cistarget/}) that this package
#' does not bundle or download automatically.
#'
#' @param obj A Seurat object.
#' @param tf_list Character vector of candidate regulator gene symbols.
#'   Required -- see Details for why there is no default.
#' @param assay Assay to use. Default \code{DefaultAssay(obj)}.
#' @param method \code{"genie3"} (default; pure R, via the \code{GENIE3}
#'   Bioconductor package) or \code{"grnboost2"} (Python, via
#'   \code{reticulate} + the \code{arboreto} package -- faster on large
#'   datasets, needs \code{reticulate::py_install("arboreto")} or
#'   \code{pip install arboreto}).
#' @param motif_rankings Optional path to (or already-loaded
#'   \code{RcisTarget}-compatible object for) a cisTarget motif ranking
#'   database matching your species/genome build. \code{NULL} (default)
#'   skips motif pruning -- see Details.
#' @param motif_annotations Optional motif-to-TF annotation table (e.g. from
#'   \code{RcisTarget::importAnnotations()}). Required alongside
#'   \code{motif_rankings} if either is supplied.
#' @param top_n_targets Maximum targets kept per regulator from the raw
#'   co-expression ranking, before motif pruning. Default 50.
#' @param min_targets Minimum targets a regulon/module must retain (after
#'   motif pruning, if run) to be kept. Default 10.
#' @param n_cores \code{method = "genie3"} only: passed to
#'   \code{GENIE3::GENIE3(nCores = )}. Default 1.
#' @param verbose Message progress. Default \code{TRUE}.
#' @return \code{obj} with \code{obj@misc$scenic} populated:
#'   \code{list(method, motif_pruned, modules, regulons, activity)} --
#'   \code{modules} is the named list of raw TF -> co-expression targets;
#'   \code{regulons} is the named list of TF -> final target genes (equal to
#'   \code{modules} if \code{motif_pruned = FALSE}); \code{activity} is the
#'   cell x regulon AUCell score matrix. Per-regulon
#'   \code{<TF>_regulon_aucell} metadata columns are also added.
#' @examples
#' \dontrun{
#' # Co-expression modules only (no motif database available)
#' obj <- RunSCENIC(obj, tf_list = human_tfs)
#' obj@misc$scenic$regulons[["STAT1"]]
#' FeaturePlot(obj, features = "STAT1_regulon_aucell")
#'
#' # Full SCENIC with motif pruning
#' rankings <- RcisTarget::importRankings("hg19-500bp-upstream-7species.mc9nr.feather")
#' annot    <- RcisTarget::importAnnotations("motifs-v9-nr.hgnc-m0.001-o0.0.tbl")
#' obj <- RunSCENIC(obj, tf_list = human_tfs,
#'                  motif_rankings = rankings, motif_annotations = annot)
#' }
#' @importFrom Seurat DefaultAssay GetAssayData
#' @export
RunSCENIC <- function(obj,
                      tf_list           = NULL,
                      assay             = NULL,
                      method            = c("genie3", "grnboost2"),
                      motif_rankings    = NULL,
                      motif_annotations = NULL,
                      top_n_targets     = 50,
                      min_targets       = 10,
                      n_cores           = 1,
                      verbose           = TRUE) {

  method <- match.arg(method)
  .assert_seurat(obj)
  if (is.null(tf_list) || length(tf_list) == 0) {
    stop("`tf_list` (a character vector of candidate regulator gene ",
         "symbols) is required -- there is no bundled default. See ?RunSCENIC.")
  }
  if (xor(is.null(motif_rankings), is.null(motif_annotations))) {
    stop("Supply both `motif_rankings` and `motif_annotations`, or neither.")
  }

  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay
  expr <- as.matrix(Seurat::GetAssayData(obj, assay = a, layer = "data"))
  tfs_present <- intersect(tf_list, rownames(expr))
  if (length(tfs_present) == 0) {
    stop("None of `tf_list` are present in assay '", a, "'.")
  }
  if (isTRUE(verbose) && length(tfs_present) < length(tf_list)) {
    message(sprintf("  %d/%d `tf_list` genes present in assay '%s'.",
                    length(tfs_present), length(tf_list), a))
  }

  # ---- Step 1: co-expression (GENIE3 or GRNBoost2) ------------------------
  if (method == "genie3") {
    if (!requireNamespace("GENIE3", quietly = TRUE)) {
      stop("'GENIE3' is required for method = 'genie3'. Install with ",
           "BiocManager::install('GENIE3').")
    }
    if (isTRUE(verbose)) {
      message(sprintf("--- Running GENIE3 (%d regulators x %d genes, %d cells) ---",
                      length(tfs_present), nrow(expr), ncol(expr)))
    }
    weight_mat <- GENIE3::GENIE3(expr, regulators = tfs_present, nCores = n_cores,
                                 verbose = isTRUE(verbose))
    link_list <- GENIE3::getLinkList(weight_mat)
    colnames(link_list) <- c("regulatoryGene", "targetGene", "weight")
  } else {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      stop("'reticulate' is required for method = 'grnboost2'.")
    }
    arboreto <- tryCatch(reticulate::import("arboreto.algo"), error = function(e) NULL)
    if (is.null(arboreto)) {
      stop("Python package 'arboreto' not found in the active reticulate ",
           "environment. Install with reticulate::py_install('arboreto') ",
           "or pip install arboreto, or use method = 'genie3' instead.")
    }
    if (isTRUE(verbose)) {
      message(sprintf("--- Running GRNBoost2 (%d regulators x %d genes, %d cells) ---",
                      length(tfs_present), nrow(expr), ncol(expr)))
    }
    expr_df <- as.data.frame(t(expr))
    link_list <- arboreto$grnboost2(expression_data = expr_df,
                                    tf_names = tfs_present)
    link_list <- reticulate::py_to_r(link_list)
    colnames(link_list) <- c("regulatoryGene", "targetGene", "weight")
  }

  link_list <- link_list[order(link_list$regulatoryGene, -link_list$weight), ]
  modules <- lapply(split(link_list, link_list$regulatoryGene), function(d) {
    utils::head(d$targetGene, top_n_targets)
  })
  modules <- modules[vapply(modules, length, integer(1)) >= min_targets]
  if (length(modules) == 0) {
    stop("No regulator retained >= min_targets (", min_targets,
         ") co-expressed targets. Lower `min_targets` or supply more data.")
  }
  if (isTRUE(verbose)) {
    message(sprintf("  %d co-expression module(s) built (>= %d targets each).",
                    length(modules), min_targets))
  }

  # ---- Step 2: motif-based pruning (RcisTarget), optional -----------------
  motif_pruned <- FALSE
  regulons <- modules
  if (!is.null(motif_rankings)) {
    if (!requireNamespace("RcisTarget", quietly = TRUE)) {
      stop("'RcisTarget' is required when `motif_rankings`/`motif_annotations` ",
           "are supplied. Install with BiocManager::install('RcisTarget').")
    }
    if (isTRUE(verbose)) {
      message(sprintf("--- Pruning %d module(s) with RcisTarget motif enrichment ---",
                      length(modules)))
    }
    motif_enrich <- RcisTarget::cisTarget(
      geneSets        = modules,
      motifRankings   = motif_rankings,
      motifAnnot      = motif_annotations,
      verbose         = isTRUE(verbose))
    motif_enrich <- as.data.frame(motif_enrich)

    regulons <- lapply(names(modules), function(tf) {
      rows <- motif_enrich[motif_enrich$geneSet == tf, , drop = FALSE]
      if (nrow(rows) == 0) return(character(0))
      targets <- unique(unlist(strsplit(rows$enrichedGenes, ";")))
      intersect(targets, modules[[tf]])
    })
    names(regulons) <- names(modules)
    regulons <- regulons[vapply(regulons, length, integer(1)) >= min_targets]
    motif_pruned <- TRUE
    if (length(regulons) == 0) {
      stop("No regulon retained >= min_targets (", min_targets,
           ") motif-supported targets after RcisTarget pruning. Lower ",
           "`min_targets`, raise `top_n_targets`, or check that ",
           "`motif_rankings`/`motif_annotations` match your species/genome.")
    }
    if (isTRUE(verbose)) {
      message(sprintf("  %d/%d module(s) retained a motif-supported regulon.",
                      length(regulons), length(modules)))
    }
  } else if (isTRUE(verbose)) {
    message("  No `motif_rankings`/`motif_annotations` supplied -- returning ",
           "raw co-expression modules as 'regulons' (no motif support checked).")
  }

  # ---- Step 3: per-cell activity scoring (AUCell, reused from RunSingleCellGSEA) --
  if (isTRUE(verbose)) message("--- Scoring regulon activity (AUCell) ---")
  activity <- .sc_gsea_score_matrix(obj, regulons, "aucell", a, verbose)
  cols <- paste0(colnames(activity), "_regulon_aucell")
  score_df <- as.data.frame(activity[colnames(obj), , drop = FALSE])
  colnames(score_df) <- cols
  obj@meta.data[, cols] <- score_df

  obj@misc$scenic <- list(method = method, motif_pruned = motif_pruned,
                          modules = modules, regulons = regulons,
                          activity = activity)
  obj
}
