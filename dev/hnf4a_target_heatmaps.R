# Pseudobulk Hnf4a target expression by Zone x Condition (split by Sex) and
# render as heatmaps -- one pass using the full CUT&RUN-derived target lists,
# a second pass restricted to only the genes actually used in section 5.8
# ("Hnf4a Target Expression by Zone") of Epigenomics_Story_Summary.Rmd.
#
# Rewritten to lean on MarkerHeatmap(genes = ..., pseudobulk = TRUE): it's
# the exact package fit for "pseudobulk a gene list per group, z-score,
# render as a diverging heatmap" -- so the custom pseudobulk-matrix builder
# and pheatmap wrapper from the original version of this script are gone.
# Two real differences from the old pheatmap-based version, both accepted
# per discussion: row clustering uses euclidean distance on the already
# z-scored matrix instead of pheatmap's 'correlation' (nearly equivalent on
# standardized rows -- euclidean distance and correlation are a monotonic
# transform of each other once rows have mean 0 / sd 1), and output is a
# ggplot saved via ggsave() rather than pheatmap's direct-to-file rendering.
library(SingleCellTools)
library(ggplot2)
library(patchwork)

## ---- 0. Zone / Condition annotation colors ----
# Zone: liver zonation gradient, pericentral -> midlobular -> periportal.
# Deliberately NOT drawn from the RdBu family the heatmap body already uses
# (row z-scores), so the annotation strip reads as a distinct signal instead
# of blending into the heatmap's own diverging color scale. RColorBrewer
# "Dark2" triplet -- qualitative, well-separated.
zone_colors <- c(
  pericentral = "#1B9E77",
  midlobular  = "#D95F02",
  periportal  = "#7570B3"
)
# Condition: Condition already encodes Sex ("Female 4mo WT", "Male 4mo WT",
# ...) -- stripped of the Sex prefix below so Female/Male heatmaps share one
# color key. WT ages 4/18/24/30mo get a light-to-dark blue ramp (chronological
# age); Ercc1-/∆ (the accelerated-aging/DNA-repair-deficient model) gets
# "salmon", reusing Epigenomics_Story_Summary.Rmd's own convention of salmon
# for the "old"/damaged condition (skyblue = young, salmon = old throughout
# that document) so this figure reads consistently with the rest of the story.
condition_colors <- c(
  "4mo WT"          = "#C6DBEF",
  "18mo WT"         = "#6BAED6",
  "24mo WT"         = "#3182BD",
  "30mo WT"         = "#08519C",
  "4mo Ercc1-/∆" = "salmon"
)

# Thin geom_tile color strip for one annotation variable, aligned to
# `col_order` -- the exact post-clustering column order MarkerHeatmap()
# landed on (cluster_cols = TRUE reorders columns internally, so this can't
# be known ahead of time; it's read back out of the returned plot's data).
build_annotation_strip <- function(col_order, value_fn, colors, legend_title) {
  values <- vapply(col_order, value_fn, character(1))
  df <- data.frame(
    col   = factor(col_order, levels = col_order),
    value = factor(values, levels = names(colors))
  )
  ggplot(df, aes(x = col, y = 1, fill = value)) +
    geom_tile(color = "black", linewidth = 0.2) +
    scale_fill_manual(values = colors, name = legend_title, drop = FALSE) +
    theme_void() +
    theme(legend.position = "right", plot.margin = margin(0, 0, 0, 0))
}

## ---- 1. Hnf4a target gene lists (same derivation as before) ----
yngTargAll <- read.delim("/gpfs/analyses/kat/PO1/Hnf4a_Epigenetics/YvO_HNF4a_KO_WL/CUTNRUN_Targets_Yng.tsv")
yng_act_targ <- unique(yngTargAll$Gene[yngTargAll$RNA_KO_FDR <= 0.05 & yngTargAll$RNA_KO_log2FC < 0])
yng_inh_targ <- unique(yngTargAll$Gene[yngTargAll$RNA_KO_FDR <= 0.05 & yngTargAll$RNA_KO_log2FC > 0])

## ---- 2. Restrict to hepatocyte-zoned spots, split by Sex, and set Zone x
##         Condition identities -- MarkerHeatmap(pseudobulk = TRUE) groups by
##         Idents(obj) (it has no group.by argument), so that's the one thing
##         the original's explicit group.by = 'Zone_Condition' needs to become. ----
# SubsetSpatial() instead of base subset() -- avoids the stale-image-slot
# Visium validity error from earlier.
hep_zoned <- SubsetSpatial(visium, subset = Zone_final %in% c("pericentral", "periportal", "midlobular"))

set_zone_condition_idents <- function(obj) {
  obj$Zone_Condition <- paste(obj$Zone_final, obj$Condition, sep = " | ")
  SeuratObject::Idents(obj) <- factor(obj$Zone_Condition)
  obj
}
female_obj <- set_zone_condition_idents(SubsetSpatial(hep_zoned, subset = Sex == "Female"))
male_obj   <- set_zone_condition_idents(SubsetSpatial(hep_zoned, subset = Sex == "Male"))

## ---- 3. Reconstruct the exact gene lists used in Rmd section 5.8 ----
# Section 5.8 subsets yng_act_targ/yng_inh_targ down to whatever genes are
# present in the MERFISH panel (rownames of MF_Pseudobulk / corrected_norm).
# Reading the same file here just to get its gene list -- not using its
# expression values, since we're plotting against the visium pseudobulk.
MF_Pseudobulk <- read.delim("/gpfs/analyses/kat/PO1/Hnf4a_Epigenetics/PseudoByZone_MerfishWTHeps.tsv")
merfish_panel_genes <- rownames(MF_Pseudobulk)
sec5_8_act_targ <- intersect(yng_act_targ, merfish_panel_genes)
sec5_8_inh_targ <- intersect(yng_inh_targ, merfish_panel_genes)
message(sprintf(
  "Section 5.8 gene counts -- Activated: %d of %d total, Inhibited: %d of %d total",
  length(sec5_8_act_targ), length(yng_act_targ),
  length(sec5_8_inh_targ), length(yng_inh_targ)
))

## ---- 4. Heatmap driver ----
# Missing-gene handling (drop + message) is now MarkerHeatmap()'s own job;
# the only thing left to handle here is the "nothing to plot at all" case,
# which MarkerHeatmap() errors on rather than silently skipping.
#
# Column labels are "<Zone> | <Condition>" (see set_zone_condition_idents());
# parsed back apart here to drive the Zone/Condition annotation strips.
.zone_of_col <- function(x) strsplit(x, " \\| ")[[1]][1]
.condition_of_col <- function(x) {
  cond <- strsplit(x, " \\| ")[[1]][2]
  sub("^(Female|Male) ", "", cond)
}

save_target_heatmap <- function(obj, genes, title, filename) {
  if (!length(genes)) {
    message(sprintf('  Skipping "%s" -- no genes in this list.', title))
    return(invisible(NULL))
  }
  p <- tryCatch(
    MarkerHeatmap(obj, genes = genes, assay = "SCT",
                 pseudobulk = TRUE, cluster_cols = TRUE),
    error = function(e) {
      message(sprintf('  Skipping "%s" -- %s', title, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(p)) return(invisible(NULL))

  col_order <- levels(p$data$cluster)
  p_zone <- build_annotation_strip(col_order, .zone_of_col, zone_colors, "Zone")
  p_cond <- build_annotation_strip(col_order, .condition_of_col, condition_colors, "Condition")

  # Fixed *absolute* strip height (not a proportion of total figure height --
  # that would grow the strips right along with tall, many-gene heatmaps).
  # patchwork's plot_layout() accepts mixed grid units: "in" pins the two
  # annotation strips to 0.15in each regardless of figure size, "null" lets
  # the heatmap panel absorb all remaining height.
  combined <- (p_cond / p_zone / p) +
    patchwork::plot_layout(
      heights = grid::unit(c(0.15, 0.15, 1), c("in", "in", "null")),
      guides  = "collect"
    ) +
    patchwork::plot_annotation(title = title)

  ggsave(filename, combined, width = 7.5, height = max(4, 0.18 * length(genes) + 2))
  invisible(combined)
}

## ---- 5. Run all 8 combinations: {full target list, section 5.8 genes} x
##         {Activated, Inhibited} x {Female, Male} ----
gene_lists <- list(
  full       = list(Activated = yng_act_targ,     Inhibited = yng_inh_targ),
  section5.8 = list(Activated = sec5_8_act_targ,  Inhibited = sec5_8_inh_targ)
)
sex_objs <- list(female = female_obj, male = male_obj)

for (scope in names(gene_lists)) {
  scope_label <- if (scope == "full") "" else "Section 5.8 Genes: "
  for (direction in names(gene_lists[[scope]])) {
    for (sex in names(sex_objs)) {
      save_target_heatmap(
        obj      = sex_objs[[sex]],
        genes    = gene_lists[[scope]][[direction]],
        title    = sprintf("%sHnf4a-%s Targets by Zone x Condition (%s)",
                           scope_label, direction, tools::toTitleCase(sex)),
        filename = sprintf("hnf4a_%s_%s_pseudobulk_heatmap_%s.pdf",
                           tolower(direction), scope, sex)
      )
    }
  }
}
message("Saved 8 heatmaps: {activated,inhibited} x {female,male} x {full target list, section 5.8 genes}")
