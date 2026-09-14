#' Allele-specific copy-number + clonal phylogeny (numbat)
#'
#' Wraps \code{numbat::run_numbat()}, which extends
#' \code{\link{InferCNVWrapper}}'s expression-only copy-number calling
#' (copyKAT) with allele-specific expression -- letting it detect
#' copy-neutral LOH that expression alone cannot see -- and additionally
#' reconstructs a clonal phylogeny from the inferred CNA profiles.
#'
#' \strong{The allele-counting/phasing step cannot be run from inside this
#' package.} Unlike copyKAT (raw counts in, everything else internal),
#' numbat requires a per-cell allele count table built by its own
#' \code{pileup_and_phase.R} script, which shells out to \code{bcftools} and
#' \code{eagle2} against a 1000 Genomes reference panel -- external
#' binaries and reference data this package does not bundle or invoke. Run
#' that preprocessing once per sample yourself, per numbat's own
#' instructions (\url{https://kharchenkolab.github.io/numbat/}), and pass
#' its output data frame in as \code{allele_df}; this function starts from
#' there. This mirrors \code{\link{RunERegulons}}'s \code{"scenicplus"}
#' method being upfront about the same kind of non-R preprocessing gap
#' rather than silently skipping it.
#'
#' @param obj A Seurat object.
#' @param allele_df The per-cell allele count data frame produced by
#'   numbat's \code{pileup_and_phase.R} preprocessing script. Required --
#'   see Details for why there is no way to compute this internally.
#' @param assay Assay to read counts from. Default \code{DefaultAssay(obj)}.
#' @param ref_expression Optional named numeric vector (gene -> reference
#'   expression level) or matrix (genes x reference cell types) to use as
#'   numbat's \code{lambdas_ref}. \code{NULL} (default): if \code{normal_cells}
#'   is supplied, builds one from those cells' mean expression; otherwise
#'   falls back to numbat's own bundled HCA reference (\code{numbat::ref_hca}),
#'   with a message, since that's a general (not sample-matched) reference.
#' @param normal_cells Optional character vector of cell barcodes believed
#'   normal, used to build \code{ref_expression} when it isn't supplied
#'   directly.
#' @param genome \code{"hg38"} (default) or \code{"hg19"}.
#' @param out_dir Directory numbat writes its (many) output files to.
#'   Default a fresh \code{tempfile()} directory.
#' @param t Transition probability between copy-number states in numbat's
#'   HMM. Default \code{1e-5}, numbat's own default.
#' @param n_cores Cores for numbat's internal parallelization. Default 1.
#' @param verbose Message progress. Default \code{TRUE}.
#' @return \code{obj} with \code{obj@misc$numbat} populated:
#'   \code{list(out_dir, nb, clone_post, segs_consensus)} -- \code{nb} is
#'   the full \code{numbat::Numbat} R6 result object (for
#'   \code{nb$plot_phylo_heatmap()} etc.); a \code{numbat_clone} metadata
#'   column is also added from \code{clone_post}.
#' @examples
#' \dontrun{
#' # allele_df built externally via numbat's pileup_and_phase.R
#' obj <- RunNumbat(obj, allele_df = allele_df,
#'                  normal_cells = colnames(obj)[obj$cell_type == "T cell"])
#' obj@misc$numbat$nb$plot_phylo_heatmap()
#' table(obj$numbat_clone)
#' }
#' @importFrom Seurat DefaultAssay GetAssayData
#' @export
RunNumbat <- function(obj,
                      allele_df      = NULL,
                      assay          = NULL,
                      ref_expression = NULL,
                      normal_cells   = NULL,
                      genome         = c("hg38", "hg19"),
                      out_dir        = tempfile("numbat_"),
                      t              = 1e-5,
                      n_cores        = 1,
                      verbose        = TRUE) {

  genome <- match.arg(genome)
  .assert_seurat(obj)
  if (!requireNamespace("numbat", quietly = TRUE)) {
    stop("'numbat' is required. Install per ",
         "https://kharchenkolab.github.io/numbat/ (available via ",
         "remotes::install_github('kharchenkolab/numbat')).")
  }
  if (is.null(allele_df)) {
    stop("`allele_df` is required: numbat's allele-specific pileup/phasing ",
         "step (pileup_and_phase.R -- bcftools + eagle2 + a 1000 Genomes ",
         "reference panel) is an external, non-R preprocessing pipeline ",
         "that can't be run from inside this package. Run it once per ",
         "sample per numbat's own instructions, then pass its output data ",
         "frame here. See ?RunNumbat Details.")
  }

  a <- if (is.null(assay)) Seurat::DefaultAssay(obj) else assay
  counts <- as.matrix(Seurat::GetAssayData(obj, assay = a, layer = "counts"))

  lambdas_ref <- ref_expression
  if (is.null(lambdas_ref)) {
    if (!is.null(normal_cells)) {
      keep <- intersect(normal_cells, colnames(counts))
      if (length(keep) == 0) {
        stop("None of `normal_cells` matched obj's cells.")
      }
      if (isTRUE(verbose)) {
        message(sprintf("--- Building reference expression from %d normal cell(s) ---",
                        length(keep)))
      }
      lambdas_ref <- Matrix::rowMeans(counts[, keep, drop = FALSE])
    } else if (!is.null(numbat::ref_hca)) {
      if (isTRUE(verbose)) {
        message("  No `ref_expression`/`normal_cells` supplied -- using ",
               "numbat's bundled HCA reference (numbat::ref_hca). This is a ",
               "general reference, not matched to your sample; supply ",
               "`normal_cells` for a sample-matched baseline if you have one.")
      }
      lambdas_ref <- numbat::ref_hca
    } else {
      stop("Supply `ref_expression` or `normal_cells` -- no bundled ",
           "reference (numbat::ref_hca) was found either.")
    }
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (isTRUE(verbose)) {
    message(sprintf("--- Running numbat (%d cells, %d genes, genome = %s) ---",
                    ncol(counts), nrow(counts), genome))
    message(sprintf("  Output files: %s", out_dir))
  }

  numbat::run_numbat(
    count_mat   = counts,
    lambdas_ref = lambdas_ref,
    df_allele   = allele_df,
    genome      = genome,
    t           = t,
    ncores      = n_cores,
    out_dir     = out_dir,
    plot        = FALSE
  )

  nb <- numbat::Numbat$new(out_dir = out_dir)
  obj@misc$numbat <- list(out_dir = out_dir, nb = nb,
                          clone_post = nb$clone_post,
                          segs_consensus = nb$segs_consensus)

  if (!is.null(nb$clone_post) && "cell" %in% colnames(nb$clone_post)) {
    clone_by_cell <- setNames(as.character(nb$clone_post$clone_opt),
                              nb$clone_post$cell)
    obj$numbat_clone <- clone_by_cell[colnames(obj)]
    if (isTRUE(verbose)) {
      message("  Clones: ", paste(names(table(obj$numbat_clone)),
                                  table(obj$numbat_clone), sep = "=", collapse = ", "))
    }
  }

  obj
}
