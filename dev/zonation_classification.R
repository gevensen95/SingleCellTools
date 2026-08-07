library(Seurat)
library(UCell)
library(dplyr)

## ---- 1. Marker panels, pre-filtered to genes present in the object ----
pericentral <- c('Glul', 'Cyp2e1', 'Oat', 'Cyp1a2')
periportal  <- c('Ass1', 'Pck1', 'Sds', 'Hal', 'Cyp2f2')
midlobular  <- c('Hamp2', 'Igfbp2', 'Ccnd1')

markers <- list(pericentral = pericentral,
                periportal  = periportal,
                midlobular  = midlobular)

markers_filtered <- lapply(markers, function(g) intersect(g, rownames(visium)))

dropped <- mapply(function(g, kept) setdiff(g, kept), markers, markers_filtered)
dropped <- dropped[lengths(dropped) > 0]
if (length(dropped) > 0) {
  message('Genes dropped (not found in visium):')
  print(dropped)
}

## ---- 2. UCell scores (only on hepatocyte-dominant spots downstream) ----
visium <- AddModuleScore_UCell(visium, features = markers_filtered)
zone_cols <- paste0(names(markers_filtered), '_UCell')
# -> pericentral_UCell, periportal_UCell, midlobular_UCell

## ---- 3. Reference cutoffs: 75th percentile of Young WT hepatocyte spots, per sex ----
# Fixed baseline -- these numbers get reused for every sample/condition below,
# they never get re-derived per sample.
cutoff_pct <- 0.75

ref <- visium@meta.data %>%
  filter(Annotation == 'Hepatocytes',
         Condition %in% c('Female 4mo WT', 'Male 4mo WT'))

if (nrow(ref) == 0) {
  stop('No hepatocyte spots found for Condition %in% c("Female 4mo WT", "Male 4mo WT") -- check that column/values exist in visium@meta.data.')
}

cutoffs <- sapply(split(ref[, zone_cols], ref$Sex),
                   function(df) apply(df, 2, quantile, probs = cutoff_pct))
# cutoffs is a zone_cols x Sex matrix -- print to sanity check before trusting it
print(cutoffs)

## ---- 4. Classify every hepatocyte spot: argmax vs. its sex-specific cutoff ----
meta <- visium@meta.data
hep_idx <- which(meta$Annotation == 'Hepatocytes')

zone_final <- as.character(meta$Annotation)  # non-hepatocyte spots keep their cell-type label

scores <- as.matrix(meta[hep_idx, zone_cols])
sex_hep <- meta$Sex[hep_idx]

if (!all(sex_hep %in% colnames(cutoffs))) {
  stop('Hepatocyte spots have Sex values not present in the reference cutoffs -- check meta$Sex levels.')
}

rel <- scores - t(cutoffs[, sex_hep])  # each spot's score relative to its own sex's cutoff

zone_call <- apply(rel, 1, function(x) {
  if (all(x <= 0)) return('Unclassified')
  sub('_UCell$', '', names(which.max(x)))
})

zone_final[hep_idx] <- zone_call
visium$Zone_final <- zone_final

table(visium$Zone_final, visium$Condition)

## ---- 5. QC: spot-check spatially before trusting the calls ----
# SpatialDimPlot(visium, group.by = 'Zone_final', images = '<sample_id>')
# SpatialFeaturePlot(visium, features = c('Glul', 'Cyp2f2'), images = '<sample_id>')
