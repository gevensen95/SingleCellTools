# =============================================================================
# Liver zonation: classification + validation figures
#
# Adapts dev/zonation_classification.R + dev/zonation_validation.R into one
# script that classifies hepatocyte zonation (pericentral / periportal /
# midlobular) from marker UCell scores against a fixed sex-specific baseline
# cutoff, then produces a full set of validation figures. Wherever a piece of
# those two scripts was doing something a SingleCellTools function already
# does, it's substituted in here:
#
#   - CompositionAnalysis() / CompositionBarplot() / CompositionalTest() /
#     CompositionEstimationPlot()  replace the manual dplyr group_by/tally +
#     geom_col zone-proportion code (validation sections 4 and 7).
#   - PlotVolcano()  replaces the bare DE table print (validation section 3).
#   - MarkerPlot()   adds a marker-panel dot plot across zone calls.
#   - PlotFeatureDensity()  adds zone-score density on the UMAP embedding
#     (validation section 8).
#   - Ol_Reliable()  replaces theme_classic() everywhere for consistent
#     package styling.
#   - call_mixture_states()  adds a NEW, independent data-driven cross-check
#     (mclust states on the same UCell scores) alongside your original
#     fixed-baseline classification -- Zone_final stays the primary call.
#   - SpatialConcordance()  replaces the hand-rolled get_spatial_concordance()
#     permutation test (validation section 1) -- same statistic (same-label
#     k-NN concordance vs. a shuffled-label null, per sample), now backed by
#     RANN instead of FNN, plus empirical p/padj columns the custom version
#     didn't compute.
#
# Synced against adams_2026_zonation_classification.R / _validation.R:
#   - pericentral panel expanded to 8 genes (added Slc1a2, Gulo, Cyp2a5, Axin2).
#   - cutoff_pct lowered from 0.75 to 0.6.
#   - holdout_pericentral / known_pericentral trimmed to drop Slc1a2/Gulo/Cyp2a5
#     now that they're in the classification panel itself -- otherwise those
#     "independent" checks would be scoring the same genes used to classify.
#     holdout_pericentral is down to a single gene (Aldh1a1) as a result; add
#     more non-overlapping pericentral markers there if you want a sturdier
#     held-out panel.
#
# Run top to bottom. Fill in CONFIG first. All figures are written, in
# order, to a single multi-page PDF at `output_pdf`.
# =============================================================================

suppressPackageStartupMessages({
  library(SingleCellTools)
  library(Seurat)
  library(UCell)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mclust)  # install.packages("mclust") if missing
})

## ---- CONFIG ----------------------------------------------------------------

# Path to an RDS of your merged, multi-sample liver Visium Seurat object.
# Must already carry metadata columns: Annotation (with a "Hepatocytes"
# level), Condition, Sex, orig.ident, plus counts under assay_name below.
visium_path <- "../VisiumLiverInsights/adams_all_genes_merged_seurat_objects.rds"

# Where to write the multi-page PDF of figures.
output_pdf <- "dev/zonation_figures.pdf"

# Assay zonation UCell scores / DE / PCA-UMAP should be computed from.
assay_name <- "SCT"

# Marker panels used for zonation classification (must exist in visium).
pericentral <- c("Glul", "Cyp2e1", "Oat", "Cyp1a2", "Slc1a2", "Gulo", "Cyp2a5", "Axin2")
periportal  <- c("Ass1", "Pck1", "Sds", "Hal", "Cyp2f2")
midlobular  <- c("Hamp2", "Igfbp2", "Ccnd1")
markers     <- list(pericentral = pericentral, periportal = periportal,
                    midlobular = midlobular)

# Held-out marker panels used ONLY for independent validation -- must not
# overlap the classification panels above. Slc1a2/Gulo/Cyp2a5 dropped here
# (and from known_pericentral below) now that they're part of the
# classification `pericentral` panel above -- holdout_pericentral is down to
# one gene as a result; add more non-overlapping pericentral markers if you
# want a sturdier held-out panel.
holdout_pericentral <- c("Aldh1a1")
holdout_periportal  <- c("Arg1", "Ttr", "Aldh1b1", "Ctsc")

# Known zonation genes checked for expected DE direction (section 5).
known_pericentral <- c("Aldh1a1", "Rnase4")
known_periportal  <- c("Arg1", "Ttr", "Aldh1b1", "Ctsc", "Serpina1b")

# Baseline condition(s) that fix the classification cutoff -- the cutoff is
# anchored here, then applied to every condition, including these.
baseline_conditions <- c("Female 4mo WT", "Male 4mo WT")
cutoff_pct <- 0.6

# Sequencing-depth metadata columns (technical confound check, section 7).
depth_col_count   <- "nCount_Spatial"
depth_col_feature <- "nFeature_Spatial"

set.seed(1)

## ---- Load -------------------------------------------------------------------

visium <- readRDS(visium_path)
stopifnot(inherits(visium, "Seurat"))
for (col in c("Annotation", "Condition", "Sex", "orig.ident")) {
  if (!col %in% colnames(visium@meta.data)) {
    stop("Expected metadata column '", col, "' not found in visium@meta.data.")
  }
}

pdf(output_pdf, width = 8, height = 6)
on.exit(dev.off(), add = TRUE)

## =============================================================================
## 1. Classification: UCell scores + sex-specific baseline cutoff -> Zone_final
## =============================================================================

markers_filtered <- lapply(markers, function(g) intersect(g, rownames(visium)))
dropped <- mapply(function(g, kept) setdiff(g, kept), markers, markers_filtered)
dropped <- dropped[lengths(dropped) > 0]
if (length(dropped) > 0) {
  message("Marker genes dropped (not found in visium):")
  print(dropped)
}

visium <- AddModuleScore_UCell(visium, features = markers_filtered)
zone_cols <- paste0(names(markers_filtered), "_UCell")

ref <- visium@meta.data %>%
  filter(Annotation == "Hepatocytes", Condition %in% baseline_conditions)
if (nrow(ref) == 0) {
  stop("No hepatocyte spots found for Condition %in% baseline_conditions -- ",
       "check that column/values exist in visium@meta.data.")
}
cutoffs <- sapply(split(ref[, zone_cols], ref$Sex),
                  function(df) apply(df, 2, quantile, probs = cutoff_pct))
message("Baseline cutoffs (", cutoff_pct, " quantile, by Sex):")
print(cutoffs)

meta       <- visium@meta.data
hep_idx    <- which(meta$Annotation == "Hepatocytes")
zone_final <- as.character(meta$Annotation)
scores     <- as.matrix(meta[hep_idx, zone_cols])
sex_hep    <- meta$Sex[hep_idx]
if (!all(sex_hep %in% colnames(cutoffs))) {
  stop("Hepatocyte spots have Sex values not present in the reference ",
       "cutoffs -- check meta$Sex levels.")
}
rel <- scores - t(cutoffs[, sex_hep])
zone_call <- apply(rel, 1, function(x) {
  if (all(x <= 0)) return("Unclassified")
  sub("_UCell$", "", names(which.max(x)))
})
zone_final[hep_idx] <- zone_call
visium$Zone_final <- zone_final

message("Zone_final counts by Condition:")
print(table(visium$Zone_final, visium$Condition))

# Two hepatocyte subsets, kept distinct on purpose:
#   hep_all        -- every hepatocyte spot, INCLUDING "Unclassified". Used
#                      wherever the Unclassified fraction itself is part of
#                      what's being checked (composition, plausibility,
#                      unsupervised cross-check).
#   hep_classified  -- Unclassified spots dropped. Used wherever a clean
#                      2-3 group comparison is needed (DE, marker plot,
#                      mixture-model cross-check) and an "Unclassified"
#                      group would just add noise.
# SubsetSpatial() (not base subset()) -- base subset() leaves obj@images
# inconsistent in a way validObject() catches (satijalab/seurat#8848 --
# unresolved upstream), surfacing later on whatever function happens to call
# validObject() first (PrepSCTFindMarkers, RunPCA, ScaleData, etc.), not at
# subset() time itself. SubsetSpatial() re-subsets every image explicitly so
# that can't happen, and -- unlike dropping images outright -- keeps them
# attached and correctly subset, so hep_all/hep_classified stay usable with
# SpatialDimPlotFixed()/SpatialConcordance() etc. if a later addition to this
# script wants that.
hep_all        <- SubsetSpatial(visium, subset = Annotation == "Hepatocytes")
hep_classified <- SubsetSpatial(visium, subset = Annotation == "Hepatocytes" & Zone_final != "Unclassified")

## =============================================================================
## 2. Composition: zone proportions per sample / condition
## =============================================================================
# Replaces the manual dplyr group_by/tally + geom_col code (originally
# validation sections 4 and 7) with CompositionAnalysis()/CompositionBarplot()/
# CompositionalTest()/CompositionEstimationPlot(). Uses hep_all so the
# Unclassified fraction is visible, matching the original plausibility check.

comp <- CompositionAnalysis(hep_all,
                            group_col     = "Zone_final",
                            sample_col    = "orig.ident",
                            condition_col = "Condition",
                            test          = "chisq")
message("CompositionAnalysis chi-square test (pooled-cell caveat applies -- see ?CompositionAnalysis):")
print(comp$test)

print(
  CompositionBarplot(comp$proportions, style = "stacked") +
    labs(title = "Zone composition per sample")
)
print(
  CompositionBarplot(comp$proportions, style = "grouped") +
    labs(title = "Zone composition per sample (grouped)")
)

# Sample-level statistical test on the Condition effect (the test to trust,
# per CompositionAnalysis()'s own pseudoreplication caveat above).
comp_test <- CompositionalTest(hep_all,
                               cluster_col   = "Zone_final",
                               sample_col    = "orig.ident",
                               condition_col = "Condition")
message("CompositionalTest() -- per-sample proportions, condition effect:")
print(comp_test)

# Map each Condition level to its Sex, and each baseline_conditions entry to
# the sex it's the baseline for -- shared by the estimation-plot loop below
# AND the sex-matched Unclassified significance check further down. Kept
# independent of whether `dabestr` is installed since the significance
# check doesn't need it.
cond_sex_map <- unique(hep_all@meta.data[, c("Condition", "Sex")])
dup_cond <- cond_sex_map$Condition[duplicated(cond_sex_map$Condition)]
if (length(dup_cond)) {
  stop("Condition value(s) map to more than one Sex -- can't auto-split ",
       "sex-matched comparisons: ", paste(unique(dup_cond), collapse = ", "))
}
cond_sex <- setNames(as.character(cond_sex_map$Sex), as.character(cond_sex_map$Condition))

baseline_sex <- cond_sex[baseline_conditions]
if (anyNA(baseline_sex)) {
  stop("baseline_conditions value(s) not found in Condition: ",
       paste(baseline_conditions[is.na(baseline_sex)], collapse = ", "))
}
if (any(duplicated(baseline_sex))) {
  stop("More than one baseline_conditions entry maps to the same Sex (",
       paste(unique(baseline_sex[duplicated(baseline_sex)]), collapse = ", "),
       ") -- need exactly one baseline per sex to auto-loop comparisons.")
}
baseline_by_sex <- setNames(baseline_conditions, baseline_sex)

# Bootstrap effect-size plot(s), one per (sex, non-baseline condition, zone)
# combination. `Condition` already bakes sex into the label (see
# baseline_conditions above -- "Female 4mo WT" vs "Male 4mo WT" are
# distinct strings), so each sex is compared only against its OWN baseline,
# never against the other sex's: liver zonation differs by sex and pooling
# would be misleading. Requires the optional `dabestr` package.
if (!requireNamespace("dabestr", quietly = TRUE)) {
  message("Package 'dabestr' not installed -- skipping CompositionEstimationPlot(). ",
          "Install with install.packages('dabestr') to include it.")
} else {
  est_plots <- list()
  for (sx in names(baseline_by_sex)) {
    ref_cond <- baseline_by_sex[[sx]]
    test_conds <- setdiff(names(cond_sex)[cond_sex == sx], ref_cond)
    for (tc in test_conds) {
      message("CompositionEstimationPlot: ", ref_cond, " vs ", tc)
      plots_sx <- tryCatch(
        CompositionEstimationPlot(comp, idx = c(ref_cond, tc)),
        error = function(e) {
          message("  Skipped (", conditionMessage(e), ")")
          NULL
        }
      )
      if (!is.null(plots_sx)) {
        if (!is.list(plots_sx) || is.null(names(plots_sx))) plots_sx <- list(zone = plots_sx)
        for (nm in names(plots_sx)) {
          est_plots[[paste(sx, tc, nm, sep = " | ")]] <- plots_sx[[nm]]
        }
      }
    }
  }
  for (nm in names(est_plots)) print(est_plots[[nm]])
}

# Plausibility check: flag samples where the classification looks degenerate.
flagged_samples <- comp$proportions %>%
  filter((group == "Unclassified" & prop > 0.5) |
         (group != "Unclassified" & prop > 0.9))
message("Samples with a degenerate zone call (>50% Unclassified or >90% one zone):")
print(flagged_samples)

# Formal test, not just eyeballing: for each sex, does the Unclassified
# fraction differ from that sex's own baseline in each non-baseline
# condition? Re-runs CompositionalTest() per (sex, condition) pair on a
# 2-level Condition subset (rather than filtering the 9-level comp_test
# above to cluster == "Unclassified") for two reasons: (1) the 9-level
# comp_test's condition effect for hep_all was pooled across BOTH sexes, and
# (2) with >2 condition levels CompositionalTest()/propeller runs an
# omnibus ANOVA (Fstatistic, no PropRatio/effect direction) -- restricting
# to exactly 2 levels here gets the t-test path back (Tstatistic +
# PropRatio), so each row has an actual direction, not just "some
# difference exists somewhere among 9 groups."
unclass_tests <- list()
for (sx in names(baseline_by_sex)) {
  ref_cond <- baseline_by_sex[[sx]]
  test_conds <- setdiff(names(cond_sex)[cond_sex == sx], ref_cond)
  for (tc in test_conds) {
    sub_obj <- subset(hep_all, subset = Condition %in% c(ref_cond, tc))
    res <- tryCatch(
      CompositionalTest(sub_obj,
                        cluster_col   = "Zone_final",
                        sample_col    = "orig.ident",
                        condition_col = "Condition"),
      error = function(e) {
        message("Unclassified test skipped for ", sx, " / ", tc, ": ",
                conditionMessage(e))
        NULL
      }
    )
    if (!is.null(res)) {
      row <- res[res$cluster == "Unclassified", , drop = FALSE]
      if (nrow(row)) {
        unclass_tests[[paste(sx, tc)]] <- cbind(
          sex = sx, reference = ref_cond, test = tc, row,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
if (length(unclass_tests)) {
  unclass_tests <- dplyr::bind_rows(unclass_tests)
  unclass_tests$padj <- stats::p.adjust(unclass_tests$pvalue, method = "BH")
  message("Sex-matched test: does Unclassified proportion differ from baseline?")
  print(unclass_tests)
} else {
  message("No sex-matched Unclassified comparisons could be run.")
}

## =============================================================================
## 3. Spatial neighbor concordance vs. permuted null
## =============================================================================
# Real zonation is laminar -- neighboring spots should share a zone call far
# more than chance. SpatialConcordance() replaces the hand-rolled
# get.knn()-based permutation loop: same statistic (same-label k-NN
# concordance per sample vs. a shuffled-label null), reusing the package's
# own coordinate-normalization helper instead of the get_spatial_concordance()
# `coord_cols <- colnames(coords)[1:2]` positional guess, and RANN instead of
# FNN. `exclude = c(NA, "Unclassified")` reproduces the original custom
# function's `!is.na(zones) & zones != "Unclassified"` filtering.

concordance_results <- SpatialConcordance(visium, group.by = "Zone_final",
                                          k = 6, n_perm = 100,
                                          exclude = c(NA, "Unclassified"))
message("Spatial neighbor concordance vs. permuted null:")
print(concordance_results)

print(
  ggplot(concordance_results, aes(sample, z)) +
    geom_col() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    Ol_Reliable() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Neighbor concordance z-score vs. permuted null", x = "",
         title = "Spatial coherence of zone calls")
)

## =============================================================================
## 4. Held-out marker validation
## =============================================================================
# Genes NOT used to build the classification -- if these track the zone
# calls too, that's independent confirmation rather than circular logic.

visium <- AddModuleScore_UCell(visium, features = list(
  holdout_pericentral = intersect(holdout_pericentral, rownames(visium)),
  holdout_periportal  = intersect(holdout_periportal, rownames(visium))
))

hep <- visium@meta.data %>% filter(Annotation == "Hepatocytes", Zone_final != "Unclassified")

print(
  ggplot(hep, aes(Zone_final, holdout_pericentral_UCell)) +
    geom_boxplot() + Ol_Reliable() +
    labs(title = "Held-out pericentral genes by zone call", x = "", y = "UCell score")
)
print(
  ggplot(hep, aes(Zone_final, holdout_periportal_UCell)) +
    geom_boxplot() + Ol_Reliable() +
    labs(title = "Held-out periportal genes by zone call", x = "", y = "UCell score")
)

message("Kruskal-Wallis, held-out pericentral score by Zone_final:")
print(kruskal.test(holdout_pericentral_UCell ~ Zone_final, data = hep))
message("Kruskal-Wallis, held-out periportal score by Zone_final:")
print(kruskal.test(holdout_periportal_UCell ~ Zone_final, data = hep))

# Dot-plot summary of both marker panels (classification + holdout) across
# zone calls -- a sanity-check that the defining genes actually differ as
# expected between groups.
gene_annot <- rbind(
  data.frame(gene = markers_filtered$pericentral, panel = "pericentral (classification)"),
  data.frame(gene = markers_filtered$periportal,  panel = "periportal (classification)"),
  data.frame(gene = markers_filtered$midlobular,  panel = "midlobular (classification)"),
  data.frame(gene = intersect(holdout_pericentral, rownames(visium)), panel = "pericentral (holdout)"),
  data.frame(gene = intersect(holdout_periportal, rownames(visium)),  panel = "periportal (holdout)")
)
Idents(hep_classified) <- hep_classified$Zone_final
print(
  MarkerPlot(hep_classified, gene_annot, assay = assay_name) +
    labs(title = "Zonation marker panels by zone call")
)

## =============================================================================
## 5. DE recovery of independent zonation genes
## =============================================================================
Idents(visium) <- visium$Zone_final

de_peri_vs_port <- FindMarkers(hep_classified, ident.1 = "pericentral", ident.2 = "periportal",
                               assay = assay_name)
de_peri_vs_port$gene <- rownames(de_peri_vs_port)

de_check <- de_peri_vs_port %>%
  filter(gene %in% c(known_pericentral, known_periportal)) %>%
  select(gene, avg_log2FC, p_val_adj)
message("Known zonation genes in pericentral-vs-periportal DE (expect known_pericentral ",
        "positive logFC, known_periportal negative):")
print(de_check)

print(
  PlotVolcano(de_peri_vs_port,
             label_genes = c(known_pericentral, known_periportal),
             title = "Pericentral vs. periportal DE")
)

## =============================================================================
## 6. Threshold sensitivity
## =============================================================================
# Rerun classification at a few cutoff percentiles -- if a biological
# conclusion (e.g. % periportal shifting with condition) only shows up at
# one specific cutoff, it's a threshold artifact, not biology.

classify_zones_at <- function(visium, cutoff_pct) {
  ref <- visium@meta.data %>%
    filter(Annotation == "Hepatocytes", Condition %in% baseline_conditions)
  cuts <- sapply(split(ref[, zone_cols], ref$Sex),
                function(df) apply(df, 2, quantile, probs = cutoff_pct))

  meta <- visium@meta.data
  hep_idx <- which(meta$Annotation == "Hepatocytes")
  zone_final <- as.character(meta$Annotation)

  scores  <- as.matrix(meta[hep_idx, zone_cols])
  sex_hep <- meta$Sex[hep_idx]
  rel     <- scores - t(cuts[, sex_hep])

  zone_call <- apply(rel, 1, function(x) {
    if (all(x <= 0)) return("Unclassified")
    sub("_UCell$", "", names(which.max(x)))
  })
  zone_final[hep_idx] <- zone_call
  zone_final
}

sensitivity <- lapply(c(0.5, 0.75, 0.9), function(p) {
  df <- visium@meta.data
  df$Zone_test <- classify_zones_at(visium, p)
  df %>% filter(Annotation == "Hepatocytes") %>%
    group_by(Condition, Zone_test) %>% tally() %>%
    group_by(Condition) %>% mutate(pct = 100 * n / sum(n), cutoff_pct = p)
})
sensitivity <- bind_rows(sensitivity)

print(
  ggplot(sensitivity %>% filter(Zone_test == "periportal"),
         aes(factor(cutoff_pct), pct, fill = Condition)) +
    geom_col(position = "dodge") +
    Ol_Reliable() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Cutoff percentile", y = "% Periportal spots",
         title = "Does the Condition effect survive different thresholds?")
)

## =============================================================================
## 7. Technical confound check: sequencing depth
## =============================================================================
hep <- visium@meta.data %>% filter(Annotation == "Hepatocytes")

print(
  ggplot(hep, aes(Zone_final, .data[[depth_col_count]])) +
    geom_boxplot() + Ol_Reliable() +
    labs(title = "Sequencing depth by zone call", x = "", y = depth_col_count)
)

message("Kruskal-Wallis, sequencing depth by Zone_final:")
print(kruskal.test(hep[[depth_col_count]] ~ hep$Zone_final))
print(kruskal.test(hep[[depth_col_feature]] ~ hep$Zone_final))

message("Spearman correlation, depth vs. zone scores:")
print(cor.test(hep[[depth_col_count]], hep$pericentral_UCell, method = "spearman"))
print(cor.test(hep[[depth_col_count]], hep$periportal_UCell, method = "spearman"))

## =============================================================================
## 8. Baseline-shift check: raw score distributions vs. fixed cutoff
## =============================================================================
cutoff_df <- as.data.frame(t(cutoffs)) %>%
  tibble::rownames_to_column("Sex") %>%
  pivot_longer(-Sex, names_to = "score", values_to = "cutoff")

score_long <- hep %>%
  select(Condition, Sex, all_of(zone_cols)) %>%
  pivot_longer(all_of(zone_cols), names_to = "score", values_to = "value")

print(
  ggplot(score_long, aes(Condition, value, fill = Condition)) +
    geom_violin() +
    geom_hline(data = cutoff_df, aes(yintercept = cutoff), linetype = "dashed") +
    facet_grid(score ~ Sex, scales = "free_x") +
    Ol_Reliable() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
    labs(title = "UCell score distributions vs. fixed baseline cutoff",
         x = "", y = "UCell score")
)
# Look for whole-sample shifts relative to the dashed line that would
# systematically over/under-classify a condition, independent of real biology.

## =============================================================================
## 9. Unsupervised cross-check
## =============================================================================
# If marker-based zones don't track any structure in an unbiased clustering
# of the same spots, the panel/resolution may not be capturing the dominant
# biological axis.

# Re-derive from hep_all (built in section 1) rather than reusing it directly,
# since PCA/UMAP/clustering below overwrite it with derived reductions/idents.
# SubsetSpatial() again, same as section 1 -- keeps images intact/consistent
# rather than dropping them, so this hep_all is still plottable if wanted.
hep_all <- SubsetSpatial(visium, subset = Annotation == "Hepatocytes")
DefaultAssay(hep_all) <- assay_name

hep_all <- RunPCA(hep_all, npcs = 30, verbose = FALSE)
hep_all <- RunUMAP(hep_all, dims = 1:30, verbose = FALSE)
hep_all <- FindNeighbors(hep_all, dims = 1:30, verbose = FALSE)
hep_all <- FindClusters(hep_all, resolution = 0.3, verbose = FALSE)

print(
  DimPlot(hep_all, group.by = "Zone_final") +
    labs(title = "Unsupervised hepatocyte embedding colored by zone call")
)

# Where each zone score concentrates on the same embedding.
print(
  PlotFeatureDensity(hep_all, features = zone_cols, reduction = "umap")
)

message("Unsupervised clusters vs. Zone_final:")
print(table(hep_all$seurat_clusters, hep_all$Zone_final))
message("Adjusted Rand Index (closer to 1 = strong agreement, 0 = no agreement):")
print(mclust::adjustedRandIndex(hep_all$seurat_clusters, hep_all$Zone_final))

## =============================================================================
## 10. Data-driven cross-check: mixture-model states on the same UCell scores
## =============================================================================
# Independent second opinion on the classification: call_mixture_states()
# fits a Gaussian mixture (BIC-selected number of components) jointly across
# the three zone score columns, with no baseline-condition anchoring at all.
# This does NOT replace Zone_final (which is deliberately anchored to the
# baseline condition so downstream comparisons are relative to it) -- it's a
# check on whether an unanchored, data-driven view of the same scores agrees.

mix <- call_mixture_states(hep_classified, score_col = zone_cols, label = "zonation")

if (!is.null(mix)) {
  hep_classified$zonation_state <- NA_character_
  hep_classified$zonation_state[match(mix$id, colnames(hep_classified))] <- as.character(mix$zonation_state)

  message("Zone_final vs. data-driven mixture state:")
  print(table(hep_classified$Zone_final, hep_classified$zonation_state, useNA = "ifany"))
  message("Adjusted Rand Index, Zone_final vs. mixture state:")
  print(mclust::adjustedRandIndex(hep_classified$Zone_final, hep_classified$zonation_state))

  print(
    ggplot(data.frame(Zone_final = hep_classified$Zone_final,
                      mixture_state = factor(hep_classified$zonation_state)),
           aes(Zone_final, fill = mixture_state)) +
      geom_bar(position = "fill") +
      Ol_Reliable() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Zone_final composition of each data-driven mixture state",
           x = "", y = "Proportion", fill = "Mixture state")
  )
} else {
  message("call_mixture_states() skipped (too few eligible rows) -- see message above.")
}

message("Done. Figures written to: ", output_pdf)
