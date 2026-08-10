# Critical evaluation of the Zone_final calls from zonation_classification.R
# Run zonation_classification.R first -- this script assumes `visium` already
# has: Annotation, Condition, Sex, pericentral_UCell/periportal_UCell/midlobular_UCell,
# and Zone_final.
#
# Adjustable bits flagged inline: image/orig.ident matching (section 1),
# assay name for depth columns (section 5), coordinate column names (section 1).

library(Seurat)
library(UCell)
library(dplyr)
library(tidyr)
library(ggplot2)
library(FNN)     # install.packages('FNN') if missing
library(mclust)  # install.packages('mclust') if missing -- adjustedRandIndex()

zone_cols <- c('pericentral_UCell', 'periportal_UCell', 'midlobular_UCell')

# recompute the fixed reference cutoffs so this script is self-contained
cutoff_pct <- 0.75
ref <- visium@meta.data %>%
  filter(Annotation == 'Hepatocytes',
         Condition %in% c('Female 4mo WT', 'Male 4mo WT'))
cutoffs <- sapply(split(ref[, zone_cols], ref$Sex),
                   function(df) apply(df, 2, quantile, probs = cutoff_pct))

## =====================================================================
## 1. Spatial neighbor concordance vs. permuted null
## =====================================================================
# Real zonation is laminar -- neighboring spots should share a zone call
# far more than chance. Compares each spot's zone against its k nearest
# spatial neighbors, then compares to a label-shuffled null per sample.

get_spatial_concordance <- function(visium, zone_col = 'Zone_final',
                                     k = 6, n_perm = 100, seed = 1) {
  set.seed(seed)
  out <- list()
  for (img in Images(visium)) {
    coords <- GetTissueCoordinates(visium, image = img)
    if ('cell' %in% colnames(coords)) {
      rownames(coords) <- coords$cell
      coords$cell <- NULL
    }
    coord_cols <- colnames(coords)[1:2]  # adjust if your Seurat version names these differently
    coords <- coords[, coord_cols]

    spots <- rownames(coords)
    spots <- spots[spots %in% rownames(visium@meta.data)]
    zones <- visium@meta.data[spots, zone_col]
    keep <- !is.na(zones) & zones != 'Unclassified'
    coords <- coords[keep, , drop = FALSE]
    zones <- zones[keep]

    if (nrow(coords) < (k + 1)) next

    nn <- get.knn(as.matrix(coords), k = k)$nn.index

    concordance <- function(z) {
      mean(vapply(seq_along(z), function(i) mean(z[nn[i, ]] == z[i]), numeric(1)))
    }

    observed <- concordance(zones)
    null_dist <- replicate(n_perm, concordance(sample(zones)))

    out[[img]] <- data.frame(
      sample = img,
      n_spots = length(zones),
      observed = observed,
      null_mean = mean(null_dist),
      null_sd = sd(null_dist),
      z_score = (observed - mean(null_dist)) / sd(null_dist)
    )
  }
  bind_rows(out)
}

concordance_results <- get_spatial_concordance(visium)
print(concordance_results)

ggplot(concordance_results, aes(sample, z_score)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = 'Neighbor concordance z-score vs. permuted null', x = '',
       title = 'Spatial coherence of zone calls')

## =====================================================================
## 2. Orthogonal (held-out) marker validation
## =====================================================================
# Genes NOT used to build the classification -- if these track the zone
# calls too, that's independent confirmation rather than circular logic.

holdout_pericentral <- intersect(c('Slc1a2', 'Cyp2a5', 'Aldh1a1', 'Gulo'), rownames(visium))
holdout_periportal  <- intersect(c('Arg1', 'Ttr', 'Aldh1b1', 'Ctsc'), rownames(visium))

visium <- AddModuleScore_UCell(visium, features = list(
  holdout_pericentral = holdout_pericentral,
  holdout_periportal  = holdout_periportal
))

hep <- visium@meta.data %>% filter(Annotation == 'Hepatocytes', Zone_final != 'Unclassified')

ggplot(hep, aes(Zone_final, holdout_pericentral_UCell)) +
  geom_boxplot() + theme_classic() +
  labs(title = 'Held-out pericentral genes by zone call', x = '', y = 'UCell score')

ggplot(hep, aes(Zone_final, holdout_periportal_UCell)) +
  geom_boxplot() + theme_classic() +
  labs(title = 'Held-out periportal genes by zone call', x = '', y = 'UCell score')

kruskal.test(holdout_pericentral_UCell ~ Zone_final, data = hep)
kruskal.test(holdout_periportal_UCell ~ Zone_final, data = hep)

## =====================================================================
## 3. DE recovery of independent zonation genes
## =====================================================================
Idents(visium) <- visium$Zone_final
hep_obj <- subset(visium, subset = Annotation == 'Hepatocytes' & Zone_final != 'Unclassified')

de_peri_vs_port <- FindMarkers(hep_obj, ident.1 = 'pericentral', ident.2 = 'periportal',
                                assay = 'SCT')  # adjust assay if needed
de_peri_vs_port$gene <- rownames(de_peri_vs_port)

known_pericentral <- c('Slc1a2', 'Cyp2a5', 'Aldh1a1', 'Gulo', 'Rnase4')
known_periportal  <- c('Arg1', 'Ttr', 'Aldh1b1', 'Ctsc', 'Serpina1b')

de_check <- de_peri_vs_port %>%
  filter(gene %in% c(known_pericentral, known_periportal)) %>%
  select(gene, avg_log2FC, p_val_adj)
print(de_check)
# expect known_pericentral genes positive logFC (pericentral = ident.1),
# known_periportal genes negative logFC

## =====================================================================
## 4. Threshold sensitivity
## =====================================================================
# Rerun classification at a few cutoff percentiles -- if a biological
# conclusion (e.g. % periportal shifting with age) only shows up at one
# specific cutoff, it's a threshold artifact, not biology.

classify_zones_at <- function(visium, cutoff_pct) {
  ref <- visium@meta.data %>%
    filter(Annotation == 'Hepatocytes',
           Condition %in% c('Female 4mo WT', 'Male 4mo WT'))
  cuts <- sapply(split(ref[, zone_cols], ref$Sex),
                 function(df) apply(df, 2, quantile, probs = cutoff_pct))

  meta <- visium@meta.data
  hep_idx <- which(meta$Annotation == 'Hepatocytes')
  zone_final <- as.character(meta$Annotation)

  scores <- as.matrix(meta[hep_idx, zone_cols])
  sex_hep <- meta$Sex[hep_idx]
  rel <- scores - t(cuts[, sex_hep])

  zone_call <- apply(rel, 1, function(x) {
    if (all(x <= 0)) return('Unclassified')
    sub('_UCell$', '', names(which.max(x)))
  })
  zone_final[hep_idx] <- zone_call
  zone_final
}

sensitivity <- lapply(c(0.5, 0.75, 0.9), function(p) {
  df <- visium@meta.data
  df$Zone_test <- classify_zones_at(visium, p)
  df %>% filter(Annotation == 'Hepatocytes') %>%
    group_by(Condition, Zone_test) %>% tally() %>%
    group_by(Condition) %>% mutate(pct = 100 * n / sum(n), cutoff_pct = p)
})
sensitivity <- bind_rows(sensitivity)

ggplot(sensitivity %>% filter(Zone_test == 'periportal'),
       aes(factor(cutoff_pct), pct, fill = Condition)) +
  geom_col(position = 'dodge') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = 'Cutoff percentile', y = '% Periportal spots',
       title = 'Does the Condition effect survive different thresholds?')

## =====================================================================
## 5. Technical confound check: sequencing depth
## =====================================================================
hep <- visium@meta.data %>% filter(Annotation == 'Hepatocytes')

depth_col_count   <- 'nCount_Spatial'    # adjust to match your object (nCount_Spatial / nCount_SCT etc.)
depth_col_feature <- 'nFeature_Spatial'  # same

ggplot(hep, aes(Zone_final, .data[[depth_col_count]])) +
  geom_boxplot() + theme_classic() +
  labs(title = 'Sequencing depth by zone call', x = '', y = depth_col_count)

kruskal.test(hep[[depth_col_count]] ~ hep$Zone_final)
kruskal.test(hep[[depth_col_feature]] ~ hep$Zone_final)

cor.test(hep[[depth_col_count]], hep$pericentral_UCell, method = 'spearman')
cor.test(hep[[depth_col_count]], hep$periportal_UCell, method = 'spearman')

## =====================================================================
## 6. Baseline-shift check: raw score distributions vs. fixed cutoff
## =====================================================================
cutoff_df <- as.data.frame(t(cutoffs)) %>%
  tibble::rownames_to_column('Sex') %>%
  pivot_longer(-Sex, names_to = 'score', values_to = 'cutoff')

score_long <- hep %>%
  select(Condition, Sex, all_of(zone_cols)) %>%
  pivot_longer(all_of(zone_cols), names_to = 'score', values_to = 'value')

ggplot(score_long, aes(Condition, value, fill = Condition)) +
  geom_violin() +
  geom_hline(data = cutoff_df, aes(yintercept = cutoff), linetype = 'dashed') +
  facet_grid(score ~ Sex, scales = 'free_x') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') +
  labs(title = 'UCell score distributions vs. fixed Young WT cutoff',
       x = '', y = 'UCell score')
# look for whole-sample shifts relative to the dashed line that would
# systematically over/under-classify a condition, independent of real biology

## =====================================================================
## 7. Plausibility check: zone proportions per sample
## =====================================================================
zone_props <- visium@meta.data %>%
  filter(Annotation == 'Hepatocytes') %>%
  group_by(orig.ident, Condition, Zone_final) %>%
  tally() %>%
  group_by(orig.ident) %>%
  mutate(pct = 100 * n / sum(n))

ggplot(zone_props, aes(orig.ident, pct, fill = Zone_final)) +
  geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(y = '% of hepatocyte spots', x = '', title = 'Zone composition per sample')

flagged_samples <- zone_props %>%
  filter((Zone_final == 'Unclassified' & pct > 50) |
         (Zone_final != 'Unclassified' & pct > 90))
print(flagged_samples)  # samples where the classification looks degenerate

## =====================================================================
## 8. Unsupervised cross-check
## =====================================================================
# If marker-based zones don't track any structure in an unbiased
# clustering of the same spots, the panel/resolution may not be
# capturing the dominant biological axis.

hep_obj <- subset(visium, subset = Annotation == 'Hepatocytes')
DefaultAssay(hep_obj) <- 'SCT'  # adjust to your default assay

hep_obj <- RunPCA(hep_obj, npcs = 30, verbose = FALSE)
hep_obj <- RunUMAP(hep_obj, dims = 1:30, verbose = FALSE)
hep_obj <- FindNeighbors(hep_obj, dims = 1:30, verbose = FALSE)
hep_obj <- FindClusters(hep_obj, resolution = 0.3, verbose = FALSE)

DimPlot(hep_obj, group.by = 'Zone_final') +
  labs(title = 'Unsupervised hepatocyte embedding colored by zone call')

table(hep_obj$seurat_clusters, hep_obj$Zone_final)
adjustedRandIndex(hep_obj$seurat_clusters, hep_obj$Zone_final)
# closer to 0 = no agreement between unsupervised structure and zone calls
# closer to 1 = strong agreement
