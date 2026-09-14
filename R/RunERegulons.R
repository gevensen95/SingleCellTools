#' Joint RNA+ATAC regulon inference (SCENIC+-style eRegulons)
#'
#' Extends \code{\link{RunSCENIC}}'s TF -> target regulons with the
#' region-level evidence SCENIC+ adds: a transcription factor's eRegulon
#' requires not just co-expression/motif support at the gene level, but a
#' motif match in an accessible peak (\code{\link{RunATACMotifEnrichment}})
#' that is itself statistically linked to the target gene
#' (\code{\link{LinkPeaksToGenes}}) -- i.e. TF -> peak -> gene, not just
#' TF -> gene.
#'
#' \strong{\code{method = "lite"} (default).} Assembles eRegulons from two
#' functions you've likely already run on \code{obj}: \code{\link{LinkPeaksToGenes}}
#' (peak-to-gene links, read from \code{Signac::Links()}) and
#' \code{\link{RunATACMotifEnrichment}} (per-peak motif matches, read from
#' \code{Signac::GetMotifData()}). For each TF motif, this takes the peaks
#' it matches, intersects those with linked peaks, and keeps the linked
#' genes as that TF's eRegulon targets -- pure R, no new heavy dependency,
#' but it is a simplified approximation of real SCENIC+ (no topic modeling,
#' no region-based co-accessibility beyond what \code{LinkPeaks()} already
#' computed).
#'
#' \strong{\code{method = "scenicplus"}.} The real \pkg{scenicplus} Python
#' package is not a single callable function -- its actual workflow is a
#' multi-stage pipeline (\code{pycisTopic} topic modeling on the ATAC data,
#' \code{pycistarget} motif enrichment on those topics, then SCENIC+'s own
#' eRegulon-linking step), normally run across several notebooks with a
#' curated cisTarget database download. This method therefore does the R-side
#' preparation only -- exporting \code{obj}'s RNA to \code{.h5ad}
#' (\code{\link{ToAnnData}}) and locating the ATAC fragments file -- and, if
#' you supply \code{python_driver} (a function or path to a Python script
#' that runs the real pipeline end-to-end and returns/writes an eRegulon
#' table), attempts to call it and parse the result into the same
#' standardized shape \code{"lite"} produces. If \code{python_driver} is
#' omitted (the common case), it does \strong{not} attempt to run
#' scenicplus itself -- it returns the prepared file paths on \code{obj}
#' and points you at \url{https://scenicplus.readthedocs.io} to run the
#' pipeline yourself, rather than silently doing nothing while looking like
#' it succeeded.
#'
#' @param obj A Seurat object with a \code{ChromatinAssay} that has already
#'   been through \code{\link{RunATACMotifEnrichment}} (\code{compute_motifs
#'   = TRUE}) and \code{\link{LinkPeaksToGenes}} -- \code{method = "lite"}
#'   reads both of those results directly off \code{obj}.
#' @param method \code{"lite"} (default) or \code{"scenicplus"}. See Details.
#' @param peak_assay Name of the \code{ChromatinAssay}. \code{NULL}
#'   (default) auto-detects it -- errors if \code{obj} has zero or more
#'   than one.
#' @param min_targets Minimum linked target genes an eRegulon must retain to
#'   be kept. Default 10.
#' @param score_method \code{"lite"} only, and the export for
#'   \code{"scenicplus"}'s activity scoring once a driver returns a result:
#'   \code{"aucell"} (default, the SCENIC-ecosystem convention) or
#'   \code{"ucell"}, passed to the same internal scorer
#'   \code{\link{RunSingleCellGSEA}} uses.
#' @param python_driver \code{method = "scenicplus"} only: an R function
#'   taking \code{(rna_h5ad, fragments_path, out_dir)} and returning a data
#'   frame with columns \code{TF}, \code{peak}, \code{gene} (one row per
#'   TF-peak-gene triplet in the final eRegulon set), or a path to a Python
#'   script that writes such a table to \code{file.path(out_dir,
#'   "eregulons.csv")}. \code{NULL} (default) skips running anything and
#'   just prepares/returns the input files. See Details.
#' @param fragments_file \code{method = "scenicplus"} only: path to the ATAC
#'   fragments file. \code{NULL} (default) reads it off
#'   \code{Signac::Fragments(obj[[peak_assay]])}.
#' @param out_dir \code{method = "scenicplus"} only: directory to write
#'   prepared inputs (and, if \code{python_driver} is a script path, its
#'   outputs) to. Default a fresh \code{tempfile()} directory.
#' @param verbose Message progress. Default \code{TRUE}.
#' @return \code{obj} with \code{obj@misc$eregulons} populated. On success
#'   (\code{"lite"}, or \code{"scenicplus"} with a working
#'   \code{python_driver}): \code{list(method, eregulons, triplets,
#'   activity)} -- \code{eregulons} is a named list, TF -> target genes;
#'   \code{triplets} is the full TF/peak/gene evidence table; \code{activity}
#'   is the cell x eRegulon score matrix (also added as
#'   \code{<TF>_eregulon_<score_method>} metadata columns). For
#'   \code{"scenicplus"} with no (or a failed) \code{python_driver}, instead
#'   \code{list(method = "scenicplus", prepared = list(rna_h5ad,
#'   fragments_file, out_dir))} -- no activity scores, since no eRegulons
#'   were actually computed.
#' @examples
#' \dontrun{
#' multiome <- RunATACMotifEnrichment(multiome, genome = BSgenome.Hsapiens.UCSC.hg38)
#' multiome <- LinkPeaksToGenes(multiome, genome = BSgenome.Hsapiens.UCSC.hg38)
#' multiome <- RunERegulons(multiome)
#' multiome@misc$eregulons$eregulons[["CTCF"]]
#' FeaturePlot(multiome, features = "CTCF_eregulon_aucell")
#'
#' # Real scenicplus, no driver -- just prepares inputs and tells you what's next
#' multiome <- RunERegulons(multiome, method = "scenicplus")
#' multiome@misc$eregulons$prepared
#' }
#' @importFrom Seurat DefaultAssay Assays
#' @export
RunERegulons <- function(obj,
                         method        = c("lite", "scenicplus"),
                         peak_assay    = NULL,
                         min_targets   = 10,
                         score_method  = c("aucell", "ucell"),
                         python_driver = NULL,
                         fragments_file = NULL,
                         out_dir       = tempfile("eregulons_"),
                         verbose       = TRUE) {

  method <- match.arg(method)
  score_method <- match.arg(score_method)
  .assert_seurat(obj)

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

  if (method == "lite") {
    links <- tryCatch(Signac::Links(obj[[pa]]), error = function(e) NULL)
    if (is.null(links) || length(links) == 0) {
      stop("No peak-gene links found on assay '", pa, "'. Run ",
           "LinkPeaksToGenes(obj, ...) first.")
    }
    links_df <- as.data.frame(links)
    if (!all(c("peak", "gene") %in% colnames(links_df))) {
      stop("Links on assay '", pa, "' are missing 'peak'/'gene' columns -- ",
           "expected Signac::LinkPeaks() output.")
    }

    motif_mat <- tryCatch(
      Signac::GetMotifData(obj, assay = pa, slot = "data"),
      error = function(e) NULL)
    motif_obj <- tryCatch(Signac::Motifs(obj[[pa]]), error = function(e) NULL)
    if (is.null(motif_mat) || is.null(motif_obj)) {
      stop("No motif data found on assay '", pa, "'. Run ",
           "RunATACMotifEnrichment(obj, ..., compute_motifs = TRUE) first.")
    }
    motif_names <- tryCatch(Signac::GetMotifData(obj, assay = pa, slot = "motif.names"),
                            error = function(e) NULL)
    tf_ids <- colnames(motif_mat)
    tf_labels <- if (!is.null(motif_names)) {
      vapply(tf_ids, function(id) {
        nm <- motif_names[[id]]
        if (is.null(nm) || length(nm) == 0) id else nm[1]
      }, character(1))
    } else {
      stats::setNames(tf_ids, tf_ids)
    }

    if (isTRUE(verbose)) {
      message(sprintf(
        "--- Assembling eRegulons: %d TF motif(s), %d peak-gene link(s) ---",
        length(tf_ids), nrow(links_df)))
    }

    triplet_list <- vector("list", length(tf_ids))
    for (i in seq_along(tf_ids)) {
      id <- tf_ids[i]
      tf_peaks <- rownames(motif_mat)[motif_mat[, id] > 0]
      hits <- links_df[links_df$peak %in% tf_peaks, , drop = FALSE]
      if (nrow(hits) == 0) { triplet_list[[i]] <- NULL; next }
      triplet_list[[i]] <- data.frame(TF = tf_labels[i], peak = hits$peak,
                                      gene = hits$gene, score = hits$score,
                                      stringsAsFactors = FALSE)
    }
    triplets <- do.call(rbind, triplet_list)
    if (is.null(triplets) || nrow(triplets) == 0) {
      stop("No TF motif matched any linked peak -- check that `obj`'s motif ",
           "scan and peak-gene links were run on the same peak set.")
    }

    eregulons <- lapply(split(triplets$gene, triplets$TF), unique)
    eregulons <- eregulons[vapply(eregulons, length, integer(1)) >= min_targets]
    triplets <- triplets[triplets$TF %in% names(eregulons), , drop = FALSE]
    rownames(triplets) <- NULL
    if (length(eregulons) == 0) {
      stop("No TF retained >= min_targets (", min_targets, ") linked target ",
           "genes. Lower `min_targets`, or check upstream peak-gene link/",
           "motif-scan quality.")
    }
    if (isTRUE(verbose)) {
      message(sprintf("  %d eRegulon(s) built (>= %d targets each).",
                      length(eregulons), min_targets))
    }

    if (isTRUE(verbose)) message(sprintf("--- Scoring eRegulon activity (%s) ---", score_method))
    a_expr <- Seurat::DefaultAssay(obj)
    ea <- if ("RNA" %in% Seurat::Assays(obj)) "RNA" else a_expr
    activity <- .sc_gsea_score_matrix(obj, eregulons, score_method, ea, verbose)
    cols <- paste0(colnames(activity), "_eregulon_", score_method)
    score_df <- as.data.frame(activity[colnames(obj), , drop = FALSE])
    colnames(score_df) <- cols
    obj@meta.data[, cols] <- score_df

    obj@misc$eregulons <- list(method = "lite", eregulons = eregulons,
                               triplets = triplets, activity = activity)
    return(obj)
  }

  # ---- method == "scenicplus" ----------------------------------------------
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("'reticulate' is required for method = 'scenicplus'.")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ff <- fragments_file
  if (is.null(ff)) {
    frags <- tryCatch(Signac::Fragments(obj[[pa]]), error = function(e) list())
    if (length(frags) == 0) {
      stop("No fragments file found on assay '", pa, "'. Pass ",
           "`fragments_file` explicitly.")
    }
    ff <- frags[[1]]@path
  }
  if (isTRUE(verbose)) message(sprintf("--- Fragments file: %s ---", ff))

  rna_h5ad <- file.path(out_dir, "rna.h5ad")
  if (isTRUE(verbose)) message(sprintf("--- Exporting RNA to %s ---", rna_h5ad))
  ToAnnData(obj, file = rna_h5ad)

  prepared <- list(rna_h5ad = rna_h5ad, fragments_file = ff, out_dir = out_dir)

  if (is.null(python_driver)) {
    message(
      "No `python_driver` supplied -- real scenicplus is a multi-stage ",
      "Python pipeline (pycisTopic topic modeling -> pycistarget motif ",
      "enrichment -> SCENIC+ eRegulon linking), not one callable function, ",
      "so it was not run. Inputs are prepared at:\n",
      "  RNA h5ad:   ", rna_h5ad, "\n",
      "  Fragments:  ", ff, "\n",
      "  Out dir:    ", out_dir, "\n",
      "Run the pipeline yourself (see https://scenicplus.readthedocs.io) ",
      "using these files, or pass `python_driver` to attempt it from here.")
    obj@misc$eregulons <- list(method = "scenicplus", prepared = prepared)
    return(obj)
  }

  result <- tryCatch({
    if (is.function(python_driver)) {
      python_driver(rna_h5ad, ff, out_dir)
    } else {
      reticulate::source_python(python_driver)
      out_csv <- file.path(out_dir, "eregulons.csv")
      if (!file.exists(out_csv)) {
        stop("`python_driver` script ran but did not write '", out_csv, "'.")
      }
      utils::read.csv(out_csv, stringsAsFactors = FALSE)
    }
  }, error = function(e) {
    message("`python_driver` failed: ", conditionMessage(e),
           "\nFalling back to returning prepared input paths only.")
    NULL
  })

  if (is.null(result) || !all(c("TF", "peak", "gene") %in% colnames(result))) {
    if (!is.null(result)) {
      message("`python_driver` result is missing required columns ",
             "(TF, peak, gene) -- falling back to prepared input paths only.")
    }
    obj@misc$eregulons <- list(method = "scenicplus", prepared = prepared)
    return(obj)
  }

  eregulons <- lapply(split(result$gene, result$TF), unique)
  eregulons <- eregulons[vapply(eregulons, length, integer(1)) >= min_targets]
  triplets <- result[result$TF %in% names(eregulons), , drop = FALSE]
  rownames(triplets) <- NULL

  if (isTRUE(verbose)) message(sprintf("--- Scoring eRegulon activity (%s) ---", score_method))
  ea <- if ("RNA" %in% Seurat::Assays(obj)) "RNA" else Seurat::DefaultAssay(obj)
  activity <- .sc_gsea_score_matrix(obj, eregulons, score_method, ea, verbose)
  cols <- paste0(colnames(activity), "_eregulon_", score_method)
  score_df <- as.data.frame(activity[colnames(obj), , drop = FALSE])
  colnames(score_df) <- cols
  obj@meta.data[, cols] <- score_df

  obj@misc$eregulons <- list(method = "scenicplus", eregulons = eregulons,
                            triplets = triplets, activity = activity,
                            prepared = prepared)
  obj
}
