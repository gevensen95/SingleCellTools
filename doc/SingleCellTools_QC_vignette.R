## ----eval=FALSE---------------------------------------------------------------
# library(SingleCellTools)
# library(Seurat)
# 
# # objs is a named list of Seurat objects, one per sample
# objs <- CreateRNAObjects(data_dirs, run_doublet_finder = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# PlotQCMetrics(obj, group.by = "orig.ident")
# 
# # Restrict to specific metrics, control the grid
# PlotQCMetrics(obj,
#               group.by = "sample",
#               qc_cols  = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
#               ncol     = 3)

## ----eval=FALSE---------------------------------------------------------------
# GenerateQCReport(objs,
#                  output_file = "qc/qc_report.html",
#                  title       = "Experiment 1 QC",
#                  sample_col  = "orig.ident")

## ----eval=FALSE---------------------------------------------------------------
# GenerateQCReport(
#   objs,
#   output_file      = "qc/qc_report.html",
#   metadata_cols    = c("nFeature_RNA", "nCount_RNA",
#                        "percent.mt", "percent.rb", "percent.hb"),
#   mad_multiplier   = 3,      # width of the suggested interval, in MADs
#   log_skewed       = TRUE,   # log-scale metrics whose skew exceeds...
#   log_threshold    = 10,     # ...this
#   complexity_score = TRUE,   # genes-per-UMI, catches low-complexity barcodes
#   doublet_col      = "doublet_finder"
# )

## ----eval=FALSE---------------------------------------------------------------
# objs_filt <- ApplyQCFilters(objs, cutoffs = "qc/qc_report_cutoffs.csv")

## ----eval=FALSE---------------------------------------------------------------
# res <- ApplyQCFilters(
#   objs,
#   cutoffs         = "qc/qc_report.html",   # sidecar resolved automatically
#   metrics         = c("nFeature_RNA", "percent.mt"),  # subset of the table
#   override        = list(percent.mt = c(0, 5)),       # beats the table
#   filter_doublets = TRUE,
#   return_report   = TRUE
# )
# 
# res$obj        # filtered objects
# res$report     # per-sample, per-metric retention
# res$cutoffs    # what was actually applied

## ----eval=FALSE---------------------------------------------------------------
# ApplyQCFilters(merged, cutoffs = cut_df, sample_col = "sample")

## ----eval=FALSE---------------------------------------------------------------
# p <- QCComparePlots(objs, objs_filt)
# p <- QCComparePlots(objs, objs_filt, metrics = c("nCount_RNA", "percent.mt"))

## ----eval=FALSE---------------------------------------------------------------
# obj <- calldoublet(obj, samplenameIndex = 1)
# 
# # usually run through the loader instead
# objs <- CreateRNAObjects(data_dirs,
#                          run_doublet_finder = TRUE,
#                          doublet_rate       = 0.075)

## ----eval=FALSE---------------------------------------------------------------
# # ~4,000 cells recovered -> ~3%
# objs <- CreateRNAObjects(data_dirs, doublet_rate = 0.03)
# 
# # per-sample rates, if recovery varies a lot
# objs <- Map(function(o, r) calldoublet(o, doublet_rate = r),
#             objs, c(Alt1 = 0.035, Alt2 = 0.05, PBS1 = 0.04))

## ----eval=FALSE---------------------------------------------------------------
# calldoublet(obj,
#             pk_sweep_max_cells = 4000,  # estimate pK on a subsample
#             sweep_cores        = 1)     # paramSweep's own parallelism

## ----eval=FALSE---------------------------------------------------------------
# pre  <- BatchEffectQC(merged, reduction = "pca",     batch_col = "sample")
# post <- BatchEffectQC(merged, reduction = "harmony", batch_col = "sample")
# 
# pre$summary
# post$summary

## ----eval=FALSE---------------------------------------------------------------
# BatchEffectQC(merged, reduction = "harmony",
#               batch_col = "sample", celltype_col = "celltype")

## ----eval=FALSE---------------------------------------------------------------
# obj <- assign_cell_cycle_phase(obj)
# table(obj$Phase)

## ----eval=FALSE---------------------------------------------------------------
# mk  <- intersect(c("Mki67", "Top2a", "Ccnb1", "Birc5"), rownames(obj))
# det <- Matrix::colSums(LayerData(obj, layer = "counts")[mk, , drop = FALSE] > 0)
# 
# for (q in c(0.5, 0.75, 0.9, 0.95)) {
#   p <- assign_cell_cycle_phase(obj, threshold_quantile = q)$Phase
#   cat(sprintf("q=%.2f  cycling=%5.1f%%  marker+ among cycling=%5.1f%%\n",
#               q, 100 * mean(p != "G1"), 100 * mean(det[p != "G1"] > 0)))
# }

## ----eval=FALSE---------------------------------------------------------------
# obj <- SCTransform(obj, vars.to.regress =
#   c("signature_1_S.Score_UCell", "signature_1_G2M.Score_UCell"))

## ----eval=FALSE---------------------------------------------------------------
# # Duplicate gene symbols across a list of objects
# check_duplicate_genes(objs)
# 
# # Gene ID type (symbol vs Ensembl vs Entrez) -- catches a mismatched reference
# check_gene_ids_across_objects(objs)
# 
# # Per-cluster summary: sizes, QC metrics, sample composition
# CellSuiteSummary(merged, cluster_col = "seurat_clusters", sample_col = "sample")

## ----eval=FALSE---------------------------------------------------------------
# vis <- CreateVisiumObjects(data_dirs, file_type = "h5")
# SpatialObjectInfo(vis)

## ----eval=FALSE---------------------------------------------------------------
# coords <- get_all_coords(obj, meta_cols = c("seurat_clusters", "nCount_Spatial"))
# head(coords)

## ----eval=FALSE---------------------------------------------------------------
# GenerateQCReport(vis,
#                  output_file      = "qc/visium_qc.html",
#                  metadata_cols    = c("nFeature_Spatial", "nCount_Spatial",
#                                       "percent.mt"),
#                  spatial_max_cols = 3)
# 
# vis_filt <- ApplyQCFilters(vis, cutoffs = "qc/visium_qc_cutoffs.csv")
# QCComparePlots(vis, vis_filt,
#                metrics = c("nCount_Spatial", "nFeature_Spatial"))

## ----eval=FALSE---------------------------------------------------------------
# edges <- EdgeDetectionVisium(seurat.obj = vis[[1]],
#                              search     = "radius",
#                              radius_factor = 1.5)
# head(edges)   # barcode, coordinates, Filter, Filter2, Filter3, Filter4

## ----eval=FALSE---------------------------------------------------------------
# vis[[1]]$edge <- edges$Filter2[match(colnames(vis[[1]]), edges$barcode)]
# vis_core <- SubsetSpatial(vis[[1]], subset = edge == "Keep")

## ----eval=FALSE---------------------------------------------------------------
# obj <- detect_fov_edges(obj, method = "bbox", n_iterations = 2)
# table(obj$edge_layer)
# 
# # drop the outermost layer only
# obj_core <- SubsetSpatial(obj, subset = edge_layer > 1)

## ----eval=FALSE---------------------------------------------------------------
# obj <- detect_tissue_holes(obj, min_hole_size = 4, sensitivity = 0.75)
# table(obj$hole_layer)

## ----eval=FALSE---------------------------------------------------------------
# detect_tissue_holes(obj, exclude_gene = "Pecam1",
#                     exclude_gene_threshold = 1)

## ----eval=FALSE---------------------------------------------------------------
# SpatialConcordance(obj, group.by = "seurat_clusters", k = 6, n_perm = 100)

## ----eval=FALSE---------------------------------------------------------------
# obj$shuffled <- sample(obj$seurat_clusters)
# SpatialConcordance(obj, group.by = "shuffled")   # z should be ~0

## ----eval=FALSE---------------------------------------------------------------
# obj_sub <- SubsetSpatial(obj, subset = celltype == "Hepatocyte")
# obj_sub <- SubsetSpatial(obj, cells = keep_barcodes)
# obj_sub <- SubsetSpatial(obj, idents = c("0", "3"))
# obj_sub <- SubsetSpatial(obj, features = keep_genes)   # images untouched

## ----eval=FALSE---------------------------------------------------------------
# # Images named image, image.1, image.2... -> meaningful names
# obj <- RenameSpatialImages(obj, group_col = "sample")
# obj <- RenameSpatialImages(obj, new_names = c(image = "slice_A"))
# 
# # Shrink an object for sharing: drop images, or keep coordinates only
# obj <- DropSpatialImage(obj, mode = "downgrade")   # coords stay, pixels go
# obj <- DropSpatialImage(obj, mode = "remove")      # everything goes

## ----eval=FALSE---------------------------------------------------------------
# # 1. What am I holding?
# vis <- CreateVisiumObjects(data_dirs, file_type = "h5")
# SpatialObjectInfo(vis)
# check_gene_ids_across_objects(vis)
# 
# # 2. Standard metric QC
# GenerateQCReport(vis, output_file = "qc/visium_qc.html",
#                  metadata_cols = c("nFeature_Spatial", "nCount_Spatial",
#                                    "percent.mt"))
# #    -> edit qc/visium_qc_cutoffs.csv here <-
# vis_f <- ApplyQCFilters(vis, cutoffs = "qc/visium_qc_cutoffs.csv",
#                         return_report = TRUE)
# QCComparePlots(vis, vis_f$obj)
# 
# # 3. Spatial artifacts -- label, look, then decide
# edges <- EdgeDetectionVisium(seurat.obj = vis_f$obj[[1]])
# vis_f$obj[[1]]$edge <- edges$Filter2[match(colnames(vis_f$obj[[1]]),
#                                            edges$barcode)]
# SpatialDimPlotFixed(vis_f$obj[[1]], group.by = "edge")
# 
# # 4. Cluster, then ask whether the clusters are spatially real
# merged <- MergeSeurat(vis_f$obj, spatial = "Visium")
# SpatialConcordance(merged, group.by = "seurat_clusters")
# 
# # 5. Shrink before saving
# saveRDS(DropSpatialImage(merged, mode = "downgrade"), "visium_qc.rds")

