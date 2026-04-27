#' Detect Gene Identifier Type in a Seurat Object
#'
#' Inspects the rownames of a Seurat object's assay and guesses whether the
#' gene identifiers are HGNC/MGI symbols, Ensembl gene IDs, RefSeq IDs, or
#' Entrez IDs. Useful as a sanity check before merging or integrating
#' objects that may have been processed against different references.
#'
#' Pattern matching is performed on a sample of up to 2000 rownames to keep
#' the call cheap on large objects. The function reports the dominant guess
#' along with the fraction of sampled rownames matching each pattern, and a
#' handful of example IDs so the result can be eyeballed.
#'
#' @param seurat_obj A Seurat object.
#' @param assay Assay whose rownames should be inspected. Defaults to
#'   \code{DefaultAssay(seurat_obj)}.
#' @return A list with three elements:
#'   \describe{
#'     \item{\code{guess}}{Single string: one of \code{"symbol"},
#'       \code{"ensembl"}, \code{"refseq"}, or \code{"entrez"}.}
#'     \item{\code{fractions}}{Named numeric vector giving the fraction of
#'       sampled rownames that match each pattern.}
#'     \item{\code{examples}}{First five rownames, for visual inspection.}
#'   }
#' @export

detect_gene_id_type <- function(seurat_obj, assay = NULL) {
  if (is.null(assay)) assay <- Seurat::DefaultAssay(seurat_obj)
  ids <- rownames(seurat_obj[[assay]])

  # Sample a chunk to keep this cheap on large objects
  sample_ids <- utils::head(ids, 2000)

  # Pattern checks
  ensembl_human   <- grepl("^ENSG\\d{11}(\\.\\d+)?$", sample_ids)
  ensembl_mouse   <- grepl("^ENSMUSG\\d{11}(\\.\\d+)?$", sample_ids)
  ensembl_generic <- grepl("^ENS[A-Z]*G\\d{11}(\\.\\d+)?$", sample_ids)
  refseq          <- grepl("^(NM_|NR_|XM_|XR_)\\d+", sample_ids)
  entrez          <- grepl("^\\d+$", sample_ids)
  symbol_like     <- grepl("^[A-Za-z][A-Za-z0-9\\.\\-]*$", sample_ids) &
    !ensembl_generic & !refseq & !entrez

  fracs <- c(
    ensembl = mean(ensembl_human | ensembl_mouse | ensembl_generic),
    refseq  = mean(refseq),
    entrez  = mean(entrez),
    symbol  = mean(symbol_like)
  )

  list(
    guess     = names(which.max(fracs)),
    fractions = fracs,
    examples  = utils::head(ids, 5)
  )
}


#' Check Gene Identifier Types Across a List of Seurat Objects
#'
#' Applies \code{\link{detect_gene_id_type}} to each Seurat object in a list
#' and returns both the per-object detection results and a tidy summary
#' data frame. Helpful for spotting mismatched gene identifier conventions
#' (e.g. some objects using Ensembl IDs, others using gene symbols) before
#' merging or integration, where mismatches would silently cause an
#' empty intersection of features.
#'
#' If the input list is unnamed, objects are labeled \code{obj_1},
#' \code{obj_2}, etc. in the summary.
#'
#' @param seurat_list A list of Seurat objects.
#' @param assay Assay whose rownames should be inspected in each object.
#'   Defaults to \code{DefaultAssay()} of each object.
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{per_object}}{List of per-object results from
#'       \code{detect_gene_id_type()}.}
#'     \item{\code{summary}}{Data frame with one row per object, including
#'       the dominant guess, the per-pattern fractions, and example IDs.}
#'   }
#' @seealso \code{\link{detect_gene_id_type}}
#' @export

check_gene_ids_across_objects <- function(seurat_list, assay = NULL) {
  results <- lapply(seurat_list, detect_gene_id_type, assay = assay)
  nms <- names(seurat_list)
  if (is.null(nms)) nms <- paste0("obj_", seq_along(seurat_list))

  summary_df <- do.call(rbind, lapply(seq_along(results), function(i) {
    r <- results[[i]]
    data.frame(
      object    = nms[i],
      guess     = r$guess,
      ensembl   = round(r$fractions["ensembl"], 3),
      refseq    = round(r$fractions["refseq"],  3),
      entrez    = round(r$fractions["entrez"],  3),
      symbol    = round(r$fractions["symbol"],  3),
      examples  = paste(r$examples, collapse = ", "),
      row.names = NULL
    )
  }))

  list(per_object = results, summary = summary_df)
}
