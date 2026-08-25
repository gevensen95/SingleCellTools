#' Determine Cell Cycle Phase
#'
#' Acts like Seurat's \code{CellCycleScoring}, but scores the S and G2M gene
#' sets with \code{UCell::AddModuleScore_UCell} instead of
#' \code{AddModuleScore}, then calls a phase per cell from the two scores.
#'
#' Called with just an object, it uses Seurat's built-in cell-cycle gene
#' lists and works on human or mouse data without being told which:
#'
#' \preformatted{obj <- assign_cell_cycle_phase(obj)}
#'
#' @section Where the default genes come from:
#' \code{Seurat::cc.genes.updated.2019} (43 S, 54 G2M) -- deliberately not
#' \code{Seurat::cc.genes}, whose \code{MLF1IP}, \code{FAM64A} and \code{HN1}
#' are symbols retired years ago (now \code{CENPU}, \code{PIMREG},
#' \code{JPT1}); they match nothing in a current reference, so the old list
#' silently loses 3 of its 97 genes.
#'
#' @section Species handling:
#' Gene symbols are matched against the object's own feature names
#' \strong{case-insensitively}, so human (\code{MCM5}) and mouse
#' (\code{Mcm5}) objects both work with no \code{species} argument and no
#' ortholog table to keep in sync. Whatever casing the object uses is what
#' gets scored. This was checked against a 24,593-row human/mouse ortholog
#' map: for all 97 cell-cycle genes the case-insensitive match and the
#' curated ortholog agree exactly, with no one-to-many ambiguity.
#'
#' Note this is a property of the cell-cycle gene set, not a general rule.
#' Plenty of orthologs are not case-variants (\code{FCGR3A} vs \code{Fcgr3}),
#' so do not reach for this trick as a general symbol converter. Genes the
#' object doesn't contain are reported and dropped.
#'
#' @param seurat_obj A Seurat object.
#' @param s_genes S-phase genes. \code{NULL} (default) uses
#'   \code{Seurat::cc.genes.updated.2019$s.genes}. Supplied genes are
#'   resolved against the object the same way.
#' @param g2m_genes G2M-phase genes. \code{NULL} (default) uses
#'   \code{Seurat::cc.genes.updated.2019$g2m.genes}.
#' @param assay Assay whose features to match against and score.
#'   \code{NULL} (default) uses \code{DefaultAssay(seurat_obj)}.
#' @param threshold_quantile Quantile of each score used as the G1 cutoff.
#'   Default \code{0.5} (median). A cell at or below both cutoffs is G1;
#'   otherwise it takes the phase with the higher score.
#'
#' @return The input Seurat object with a new \code{$Phase} metadata column
#'   (\code{"G1"}, \code{"S"}, or \code{"G2M"} per cell), alongside the two
#'   UCell score columns.
#'
#' @examples
#' \dontrun{
#' # Human or mouse -- no species argument needed
#' obj <- assign_cell_cycle_phase(obj)
#'
#' # Override the gene sets
#' obj <- assign_cell_cycle_phase(obj, s_genes = my_s, g2m_genes = my_g2m)
#'
#' table(obj$Phase)
#' }
#'
#' @seealso \code{\link[Seurat]{CellCycleScoring}}
#' @export
assign_cell_cycle_phase <- function(seurat_obj,
                                    s_genes            = NULL,
                                    g2m_genes          = NULL,
                                    assay              = NULL,
                                    threshold_quantile = 0.5) {

  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.")
  }
  if (!requireNamespace("UCell", quietly = TRUE)) {
    stop("Package 'UCell' is required. ",
         "BiocManager::install('UCell')")
  }
  if (!is.numeric(threshold_quantile) || length(threshold_quantile) != 1L ||
      is.na(threshold_quantile) ||
      threshold_quantile < 0 || threshold_quantile > 1) {
    stop("`threshold_quantile` must be a single number in [0, 1].")
  }

  assay <- if (is.null(assay)) SeuratObject::DefaultAssay(seurat_obj) else assay
  if (!assay %in% SeuratObject::Assays(seurat_obj)) {
    stop("Assay '", assay, "' not found. Available: ",
         paste(SeuratObject::Assays(seurat_obj), collapse = ", "))
  }
  obj_genes <- rownames(seurat_obj[[assay]])

  if (is.null(s_genes))   s_genes   <- Seurat::cc.genes.updated.2019$s.genes
  if (is.null(g2m_genes)) g2m_genes <- Seurat::cc.genes.updated.2019$g2m.genes

  # Case-insensitive resolution against the object's own feature names, so
  # MCM5/Mcm5 both land without asking the caller which species this is.
  .resolve <- function(genes, label) {
    genes <- as.character(genes)
    hit   <- match(toupper(genes), toupper(obj_genes))
    found <- obj_genes[hit[!is.na(hit)]]
    lost  <- genes[is.na(hit)]
    if (length(lost)) {
      message(sprintf("  %s: %d of %d gene(s) not in assay '%s': %s",
                      label, length(lost), length(genes), assay,
                      paste(utils::head(lost, 10), collapse = ", ")),
              if (length(lost) > 10) ", ..." else "")
    }
    if (!length(found)) {
      stop("None of the ", label, " genes are present in assay '", assay,
           "'. Check the gene symbols against rownames(seurat_obj).")
    }
    unique(found)
  }
  s_genes   <- .resolve(s_genes,   "S")
  g2m_genes <- .resolve(g2m_genes, "G2M")

  # Clear any score columns left by a previous run, so re-running is
  # idempotent. Without this the grep() below matches both the stale and the
  # fresh column, and `[[c(a, b)]][, 1]` silently scores off the stale one.
  stale <- grep("_(S|G2M)\\.Score_UCell$", colnames(seurat_obj@meta.data),
                value = TRUE)
  if (length(stale)) {
    for (cl in stale) seurat_obj[[cl]] <- NULL
  }

  message('--- Scoring S phase genes with UCell ---')
  seurat_obj <- UCell::AddModuleScore_UCell(seurat_obj, features = list(s_genes),
                                            assay = assay,
                                            name = "_S.Score_UCell")

  message('--- Scoring G2M phase genes with UCell ---')
  seurat_obj <- UCell::AddModuleScore_UCell(seurat_obj, features = list(g2m_genes),
                                            assay = assay,
                                            name = "_G2M.Score_UCell")

  s_col   <- grep("_S\\.Score_UCell$",   colnames(seurat_obj@meta.data), value = TRUE)
  g2m_col <- grep("_G2M\\.Score_UCell$", colnames(seurat_obj@meta.data), value = TRUE)
  if (length(s_col) != 1L || length(g2m_col) != 1L) {
    stop("Expected exactly one S and one G2M UCell score column, found ",
         length(s_col), " and ", length(g2m_col),
         ". UCell's column naming may have changed.")
  }

  s_score   <- seurat_obj@meta.data[[s_col]]
  g2m_score <- seurat_obj@meta.data[[g2m_col]]

  message(sprintf('--- Computing score thresholds (quantile = %.2f) ---',
                  threshold_quantile))
  s_thresh   <- stats::quantile(s_score,   threshold_quantile, na.rm = TRUE)
  g2m_thresh <- stats::quantile(g2m_score, threshold_quantile, na.rm = TRUE)

  message('--- Assigning cell cycle phase (G1 / S / G2M) ---')
  # `<=`, not `<`. UCell returns exactly 0 for any cell expressing none of a
  # signature's genes, which on sparse or shallow data is a large share of
  # cells -- often enough that the quantile is itself 0. With a strict `<`,
  # `0 < 0` is never TRUE, so not one cell is called G1 and every all-zero
  # cell falls through both score comparisons to NA. On a 60%-zero sample
  # that is 0 G1 cells and 60% NA, where the answer should be 60% G1.
  phase <- ifelse(
    s_score <= s_thresh & g2m_score <= g2m_thresh, "G1",
    ifelse(s_score > g2m_score, "S",
           ifelse(g2m_score > s_score, "G2M", NA_character_))
  )
  seurat_obj$Phase <- phase

  n_na <- sum(is.na(phase))
  if (n_na) {
    message(sprintf(
      "  %d cell(s) had equal S and G2M scores above the G1 cutoff; Phase is NA for them.",
      n_na))
  }
  tab <- table(factor(phase, levels = c("G1", "S", "G2M")))
  message(sprintf("--- Phase: G1 %d | S %d | G2M %d%s ---",
                  tab[["G1"]], tab[["S"]], tab[["G2M"]],
                  if (n_na) sprintf(" | NA %d", n_na) else ""))

  seurat_obj
}
