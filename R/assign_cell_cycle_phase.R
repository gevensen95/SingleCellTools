#' Determine Cell Cycle Phase
#'
#' This function acts the same as Seurat's CellCycleScoring, but it uses AddModuleScore_UCell
#'
#' @param seurat_obj Seurat Object
#' @param s_genes S Phase Genes (must be species specific)
#' @param g2m_genes G2M Phase Genes (must be species specific)
#' @param threshold_quantile Threshold for G1 Phase, 0.5 is the recommended threshold
#' @return Vector containing child GO terms
#' @export
#'
assign_cell_cycle_phase <- function(seurat_obj, s_genes, g2m_genes,
                                    threshold_quantile = 0.5) {
  # Calculate UCell scores
  seurat_obj <- AddModuleScore_UCell(seurat_obj, features = list(s_genes), name = "_S.Score_UCell")
  seurat_obj <- AddModuleScore_UCell(seurat_obj, features = list(g2m_genes), name = "_G2M.Score_UCell")

  # Extract column names created by UCell
  s_col <- grep("_S.Score_UCell$", colnames(seurat_obj@meta.data), value = TRUE)
  g2m_col <- grep("_G2M.Score_UCell$", colnames(seurat_obj@meta.data), value = TRUE)

  # Calculate thresholds
  s_thresh <- quantile(seurat_obj[[s_col]][,1], threshold_quantile)
  g2m_thresh <- quantile(seurat_obj[[g2m_col]][,1], threshold_quantile)

  # Assign cell cycle phase
  seurat_obj$Phase <- with(seurat_obj@meta.data, ifelse(
    get(s_col) < s_thresh & get(g2m_col) < g2m_thresh, "G1",
    ifelse(get(s_col) > get(g2m_col), "S",
           ifelse(get(g2m_col) > get(s_col), "G2M", NA)
    )
  ))

  return(seurat_obj)
}
