## ----eval=FALSE---------------------------------------------------------------
# # install.packages("remotes")
# remotes::install_github("gevensen95/SingleCellTools")

## ----eval=FALSE---------------------------------------------------------------
# remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")

## ----eval=FALSE---------------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c(
#   # Core scRNA-seq / spatial
#   "EnsDb.Mmusculus.v79",  # mouse gene annotations
#   "glmGamPoi",            # fast SCTransform fitting
#   "GO.db",                # GO term lookups
#   "UCell",                # module scoring
#   "Signac",               # scATAC-seq
#   # Differential expression and composition
#   "DESeq2",
#   "SummarizedExperiment",
#   "S4Vectors",
#   "edgeR",
#   "speckle",              # propeller composition test
#   "limma",
#   # Trajectory
#   "slingshot",
#   # Reference-based annotation
#   "SingleR",              # optional; used by AnnotateClusters(method="singler")
#   # scanpy interop
#   "zellkonverter"
# ))

## ----eval=FALSE---------------------------------------------------------------
# # Ligand-receptor inference
# remotes::install_github("saezlab/liana")
# 
# # scmap for the R-native branch of AnnotateWithReference
# BiocManager::install("scmap")
# 
# # CellTypist (Python default backend of AnnotateWithReference)
# reticulate::py_install("celltypist")
# # Or: pip install celltypist
# 
# # scANVI (optional Python backend of AnnotateWithReference)
# reticulate::py_install("scvi-tools")
# # Or: pip install scvi-tools
# 
# # Visium deconvolution
# remotes::install_github("dmcable/spacexr")

## ----eval=FALSE---------------------------------------------------------------
# install.packages(c(
#   "patchwork", "ks", "MASS", "cluster",
#   "jsonlite", "betareg"
# ))

## ----eval=FALSE---------------------------------------------------------------
# install.packages("sf")     # required only for get_cells_in_polygon()
# install.packages("anndata"); reticulate::install_python()  # for FromAnnData Python reader

## ----eval=FALSE---------------------------------------------------------------
# library(SingleCellTools)

## ----eval=FALSE---------------------------------------------------------------
# library(SingleCellTools)
# 
# # Paths to per-sample CellRanger output folders
# samples <- c(
#   "data/vehicle_rep1",
#   "data/vehicle_rep2",
#   "data/drugA_rep1",
#   "data/drugA_rep2"
# )
# 
# seurat_list <- CreateRNAObjects(
#   data_dirs    = samples,
# 
#   # Basic QC thresholds applied at object creation
#   cells        = 3,          # a gene must appear in ≥3 cells to be kept
#   features     = 200,        # a cell must express ≥200 genes to be kept
# 
#   # Mitochondrial read percentage
#   mt_pattern   = "^mt-",     # mouse mitochondrial genes; use "^MT-" for human
# 
#   # Optional: add a Treatment column to metadata
#   treatment    = c("Vehicle", "Vehicle", "DrugA", "DrugA"),
# 
#   # Optional: give the list elements readable names
#   object_names = c("Veh1", "Veh2", "Drug1", "Drug2"),
# 
#   # Doublet detection (runs calldoublet() on every sample)
#   run_doublet_finder        = TRUE,
#   doublet_normalization     = "LogNormalize",   # or "SCT"
#   doublet_vars_to_regress   = "percent.mt",
#   doublet_cluster_resolution = 0.1,
# 
#   # Set to TRUE to immediately drop doublets from each object;
#   # FALSE (default) keeps them so you can inspect the labels first
#   filter_doublets = FALSE
# )

## ----eval=FALSE---------------------------------------------------------------
# # Violin plots per sample — decide your nFeature and percent.mt cutoffs here
# library(Seurat)
# library(patchwork)
# 
# VlnPlot(seurat_list[["Veh1"]], features = c("nFeature_RNA", "percent.mt"))

## ----eval=FALSE---------------------------------------------------------------
# seurat_list <- lapply(seurat_list, function(obj) {
#   subset(obj,
#     nFeature_RNA > 500  &
#     nFeature_RNA < 6000 &
#     percent.mt   < 20   &
#     doublet_finder == "Singlet"
#   )
# })

## ----eval=FALSE---------------------------------------------------------------
# seurat_list <- CreateRNAObjects(
#   data_dirs = samples,
#   treatment = c("Vehicle", "Vehicle", "DrugA", "DrugA"),
#   on_disk   = TRUE,
#   bpcells_dir = "bpcells_cache"   # default: "bpcells" under the working directory
# )

## ----eval=FALSE---------------------------------------------------------------
# seurat_list <- CreateRNAObjectsFilter(
#   data_dirs  = samples,
#   mt_pattern = "^mt-",
#   treatment  = c("Vehicle", "Vehicle", "DrugA", "DrugA")
# )

## ----eval=FALSE---------------------------------------------------------------
# # Run on a single object
# obj_with_doublets <- calldoublet(
#   obj                = seurat_list[["Veh1"]],
#   samplenameIndex    = 1,              # unused internally; kept for compatibility
#   normalization      = "LogNormalize", # or "SCT"
#   vars.to.regress    = "percent.mt",
#   cluster_resolution = 0.1,
#   pk_sweep_max_cells = 4000,           # estimate pK from <=4000 cells instead of every cell
#   sweep_cores        = 1               # DoubletFinder's own num.cores, for its 6-value pN sweep
# )
# 
# # Check the result
# table(obj_with_doublets$doublet_finder)
# #> Doublet Singlet
# #>     312    5841

## ----eval=FALSE---------------------------------------------------------------
# integrated <- MergeSeurat(
#   seurat_objects = seurat_list,
# 
#   # Variables to regress during normalization
#   to_regress = "percent.mt",
# 
#   # Normalization method
#   use_SCT = TRUE,           # FALSE uses LogNormalize / FindVariableFeatures / ScaleData
# 
#   # Dimensionality
#   max_dims       = 20,      # number of PCs to use (ignored if use_elbow_plot = TRUE)
#   use_elbow_plot = FALSE,   # set TRUE to pick dims interactively at the console
# 
#   # Integration
#   integration     = "HarmonyIntegration",  # see table below for alternatives
#   new_reduction   = "harmony",
# 
#   # Clustering
#   cluster_resolution = 0.3,
# 
#   # Save outputs
#   save_rds_file = TRUE,
#   file_name     = "lung_experiment",
# 
#   # Automatically run FindAllMarkers after clustering
#   markers = TRUE,
# 
#   # Restrict to genes shared across all samples before merging
#   common_genes_only = FALSE
# )

## ----eval=FALSE---------------------------------------------------------------
# # Example with RPCA
# integrated <- MergeSeurat(
#   seurat_objects  = seurat_list,
#   integration     = "RPCAIntegration",
#   new_reduction   = "integrated.rpca",
#   k_anchor        = 20,
#   k_weight        = 100
# )

## ----eval=FALSE---------------------------------------------------------------
# integrated_visium <- MergeSeurat(
#   seurat_objects = visium_list,
#   spatial        = "Visium",
#   use_SCT        = TRUE,
#   integration    = "HarmonyIntegration",
#   new_reduction  = "harmony"
# )

## ----eval=FALSE---------------------------------------------------------------
# tcell <- SubsetAndRecluster(
#   integrated,
#   metadata_col   = "cell_type",
#   metadata_value = c("CD4 T", "CD8 T", "Treg"),
#   assay          = "RNA",
#   normalize      = "auto",            # detect missing data/scale.data layers
#   integrate      = TRUE,
#   integration_method = "HarmonyIntegration",
#   dims           = 20,
#   resolution     = 0.3
# )

## ----eval=FALSE---------------------------------------------------------------
# # By cluster id
# SubsetAndRecluster(integrated, idents = c("0", "3", "7"))
# 
# # By explicit cells
# SubsetAndRecluster(integrated, cells = my_barcodes)

## ----eval=FALSE---------------------------------------------------------------
# # Define a marker panel as a data frame
# markers <- data.frame(
#   Gene = c(
#     "Sftpc", "Sftpb", "Abca3",        # AT2
#     "Ager",  "Hopx",  "Pdpn",         # AT1
#     "Trp63", "Krt5",  "Krt14",        # Basal
#     "Scgb1a1", "Scgb3a2",             # Club
#     "Pecam1", "Cdh5", "Kdr",          # Endothelial
#     "Col1a1", "Postn", "Acta2"        # Fibroblast
#   ),
#   CellType = c(
#     rep("AT2", 3),
#     rep("AT1", 3),
#     rep("Basal", 3),
#     rep("Club", 2),
#     rep("Endothelial", 3),
#     rep("Fibroblast", 3)
#   )
# )
# 
# # Set identities to the annotated column before plotting
# Idents(integrated) <- "seurat_clusters"
# 
# p <- MarkerPlot(
#   obj              = integrated,
#   genes            = markers,
#   assay            = "RNA",
#   cluster          = TRUE,      # cluster identities by expression correlation
#   show.annotations = TRUE,      # draw cell-type labels along right edge
#   maxsize          = 5,         # maximum dot size
#   label.fontsize   = 3,         # annotation label size
#   margin_factor    = 0.5        # increase if labels are clipped
# )
# print(p)
# 
# # Save
# ggsave("markerplot_lung.pdf", p, width = 12, height = 8)

## ----eval=FALSE---------------------------------------------------------------
# MarkerPlot(integrated, markers, save_path = "markerplot_lung.pdf")

## ----eval=FALSE---------------------------------------------------------------
# Idents(integrated) <- factor(
#   integrated$cell_type,
#   levels = c("AT2", "AT1", "Club", "Basal", "Endothelial", "Fibroblast")
# )
# MarkerPlot(integrated, markers, cluster = FALSE)  # disable correlation clustering

## ----eval=FALSE---------------------------------------------------------------
# MarkerPctPlot2(
#   integrated, cellID,
#   style        = "tile",         # or "dot"
#   colors       = c("white", "firebrick"),
#   max_pct      = 80,             # spread color range over informative 0-80%
#   label_offset = 1.2             # right-margin cell-type labels
# )

## ----eval=FALSE---------------------------------------------------------------
# markers_list <- list(
#   T_cell = c("CD3D", "CD3E", "CD8A", "CD4"),
#   B_cell = c("MS4A1", "CD79A", "CD19"),
#   NK     = c("NKG7", "GNLY", "KLRD1"),
#   Mono   = c("CD14", "LYZ", "FCGR3A")
# )
# 
# integrated <- AnnotateClusters(
#   integrated,
#   method              = "marker",
#   markers             = markers_list,
#   cluster_col         = "seurat_clusters",
#   new_col             = "predicted_cell_type",
#   filter_nonspecific  = TRUE,      # drop broadly-expressed "markers"
#   min_score           = 0.1,
#   min_margin          = 0.05
# )
# table(integrated$predicted_cell_type)

## ----eval=FALSE---------------------------------------------------------------
# obj <- AnnotateClusters(obj, markers = markers_list, tumor_mode = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# res <- AnnotateClusters(
#   visium, markers = markers_list,
#   filter_nonspecific  = FALSE,
#   min_detection_frac  = 0.05,
#   return_scores       = "cluster"
# )
# res$obj               # annotated
# res$scores            # cluster x cell-type matrix; heatmap this

## ----eval=FALSE---------------------------------------------------------------
# integrated <- AnnotateWithReference(
#   integrated,
#   method          = "celltypist",
#   model           = "Immune_All_Low.pkl",   # or Immune_All_High, Lung, etc.
#   majority_voting = TRUE,                    # smoothing via over-clustering
#   min_score       = 0.5
# )
# table(integrated$predicted_cell_type)

## ----eval=FALSE---------------------------------------------------------------
# integrated <- AnnotateWithReference(
#   integrated,
#   method        = "scanvi",
#   reference     = pbmc_ref_seurat,
#   ref_label_col = "cell_type",
#   batch_col     = "orig.ident",
#   n_latent      = 30
# )

## ----eval=FALSE---------------------------------------------------------------
# integrated <- AnnotateWithReference(
#   integrated,
#   method        = "scmap",
#   reference     = pbmc_ref_seurat,
#   ref_label_col = "cell_type",
#   scmap_method  = "cluster",     # or "cell"
#   threshold     = 0.5,
#   n_features    = 500
# )

## ----eval=FALSE---------------------------------------------------------------
# markers_by_sample <- FindMarkersList(
#   seurat_list,
#   idents_col       = "cell_ontology_class",
#   padj_threshold   = 0.05,
#   log2fc_threshold = 1
# )
# 
# # Combine across samples
# all_markers <- do.call(rbind, markers_by_sample)

## ----eval=FALSE---------------------------------------------------------------
# # Single object
# integrated <- AddGenePositivity(
#   seurat_objects = integrated,
#   genes          = c("Sftpc", "Ager", "Trp63"),
#   layer          = "counts",   # use "data" for log-normalized layer
#   threshold      = 0,          # cells with counts > 0 are positive
#   suffix         = "_pos"
# )
# 
# # New metadata columns: Sftpc_pos, Ager_pos, Trp63_pos (logical)
# head(integrated@meta.data[, c("Sftpc_pos", "Ager_pos", "Trp63_pos")])

## ----eval=FALSE---------------------------------------------------------------
# seurat_list <- AddGenePositivity(
#   seurat_objects = seurat_list,
#   genes          = c("Sftpc", "Ager", "Trp63"),
#   layer          = "counts"
# )

## ----eval=FALSE---------------------------------------------------------------
# # Subset to AT2 cells (Sftpc-positive)
# at2 <- subset(integrated, Sftpc_pos == TRUE)
# 
# # Fraction positive per cluster
# tapply(integrated$Sftpc_pos, Idents(integrated), mean)

## ----eval=FALSE---------------------------------------------------------------
# # Percent-positive bar chart per cluster, one bar per gene
# PlotGenePositivity(integrated, c("Sftpc", "Ager", "Trp63"))
# 
# # Many genes → tile heatmap
# PlotGenePositivity(integrated, big_gene_vec, style = "heatmap")
# 
# # Co-expression combinations (CD3D+/CD4+, CD3D+/CD4-, ...)
# PlotGenePositivity(integrated, c("CD3D", "CD4"), style = "combo")
# 
# # For a list of samples, one facet per sample
# PlotGenePositivity(seurat_list, c("Sftpc", "Ager"))

## ----eval=FALSE---------------------------------------------------------------
# gpa <- GenePositivityAnalysis(
#   integrated,
#   genes         = c("Sftpc", "Ager"),
#   sample_col    = "orig.ident",
#   condition_col = "treatment",
#   group_col     = "cell_type",
#   test          = "chisq"
# )
# gpa$proportions       # long tibble: sample, gene, group, n_pos, n_total, prop_pos, condition
# gpa$test[["Sftpc | AT2"]]   # chi-square result for that gene x cell-type combination (pooled-cell caveat above)

## ----eval=FALSE---------------------------------------------------------------
# # One gene x cell-type combination
# GenePositivityEstimationPlot(gpa, genes = "Sftpc", group_levels = "AT2",
#                              idx = c("Vehicle", "DrugA"))
# 
# # Every combination present, as a named list of plots
# plots <- GenePositivityEstimationPlot(gpa, idx = c("Vehicle", "DrugA"))
# plots[["Sftpc | AT2"]]
# 
# # Cohen's h -- the effect size designed specifically for comparing two proportions
# GenePositivityEstimationPlot(gpa, idx = c("Vehicle", "DrugA"), effect = "cohens_h")

## ----eval=FALSE---------------------------------------------------------------
# # A composite ER-stress/UPR module score, as an example of a continuous score
# # worth binning into discrete states
# integrated <- AddModuleScore(
#   integrated,
#   features = list(c("Hspa5", "Ddit3", "Atf4")),
#   name     = "StressScore"
# )
# integrated$stress_composite <- integrated$StressScore1
# 
# states <- call_mixture_states(
#   integrated,
#   score_col   = "stress_composite",
#   group_col   = "cell_type",
#   group_value = "AT2",     # fit the mixture within AT2 cells only
#   label       = "stress"
# )
# head(states)
# attr(states, "G")     # number of states BIC selected
# attr(states, "bic")

## ----eval=FALSE---------------------------------------------------------------
# integrated$stress_state <- states$stress_state[match(colnames(integrated), states$id)]
# VlnPlot(integrated, features = "stress_composite", group.by = "stress_state")

## ----eval=FALSE---------------------------------------------------------------
# # No stratification -- fit on every cell regardless of cell_type
# call_mixture_states(integrated, score_col = "stress_composite")
# 
# # Multivariate: joint mixture across several module scores at once
# call_mixture_states(integrated, score_col = c("stress_composite", "SenescenceScore1"),
#                     label = "combined_state")

## ----eval=FALSE---------------------------------------------------------------
# integrated$annotation_first_pass <- integrated$cell_type   # match call_stress_states()'s fixed column name
# res <- call_stress_states(integrated, cell_type = "AT2")
# res$calls           # data.frame: cell, stress_state, stress_prob, prob_stressed

## ----eval=FALSE---------------------------------------------------------------
# integrated <- CreateAndIntegrateRNA(
#   data_dirs          = samples,
#   mt_pattern         = "^mt-",
#   treatment          = c("Vehicle", "Vehicle", "DrugA", "DrugA"),
#   nFeature_min       = 500,
#   nFeature_max       = 6000,
#   percent_mt_max     = 20,
#   run_doublet_finder = TRUE,
#   filter_doublets    = TRUE,
#   integration        = "HarmonyIntegration",
#   cluster_resolution = 0.3,
#   markers            = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# GenerateQCReport(
#   sample_list,
#   output_file   = "qc/qc_report.html",
#   metadata_cols = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
#   mad_multiplier = 3,          # 3 MAD = "outlier" threshold
#   doublet_col   = "doublet_finder"
# )

## ----eval=FALSE---------------------------------------------------------------
# sample_list_filtered <- ApplyQCFilters(
#   sample_list,
#   cutoffs         = "qc/qc_report_cutoffs.csv",
#   filter_doublets = TRUE,
#   doublet_col     = "doublet_finder",
#   doublet_value   = "Doublet",
#   return_report   = TRUE          # for a per-sample retention table
# )
# sample_list_filtered$report        # sample x metric x pct_kept

## ----eval=FALSE---------------------------------------------------------------
# # Only some metrics — leave doublet score / percent.hb alone
# ApplyQCFilters(sample_list, cutoffs = "qc/qc_report_cutoffs.csv",
#                metrics = c("nCount_RNA", "nFeature_RNA", "percent.mt"))
# 
# # Override a single sample's mito cutoff without editing the CSV
# ApplyQCFilters(sample_list, cutoffs = "qc/qc_report_cutoffs.csv",
#                override = list(sample1 = list(percent.mt = c(0, 15))))

## ----eval=FALSE---------------------------------------------------------------
# QCComparePlots(
#   pre     = sample_list,
#   post    = sample_list_filtered$obj,
#   metrics = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
#   log_y   = c("nCount_RNA", "nFeature_RNA")
# )

## ----eval=FALSE---------------------------------------------------------------
# # Single object, auto-detected columns
# PlotQCMetrics(sample_list$sample1)
# 
# # A list of objects -- metadata is row-bound across samples, no need to merge first
# PlotQCMetrics(sample_list, group.by = "orig.ident")
# 
# # Explicit columns, custom panel layout
# PlotQCMetrics(sample_list,
#              qc_cols = c("nFeature_RNA", "nCount_RNA", "percent.mt", "doublet_finder"),
#              ncol    = 2)

## ----eval=FALSE---------------------------------------------------------------
# # Single cluster / cell type
# res <- PseudobulkDE(
#   integrated,
#   sample_col    = "orig.ident",
#   condition_col = "treatment",
#   ident_1       = "DrugA",
#   ident_2       = "Vehicle",
#   group_by      = "cell_type",
#   group_value   = "T cell"
# )
# head(res$results)                # DESeq2 results table
# head(res$normalized_counts)      # size-factor-normalized pseudobulk counts

## ----eval=FALSE---------------------------------------------------------------
# res_all <- PseudobulkDE(
#   integrated,
#   sample_col    = "orig.ident",
#   condition_col = "treatment",
#   ident_1       = "DrugA",
#   ident_2       = "Vehicle",
#   cluster_col   = "cell_type"
# )
# res_all$T_cell$results            # per-cluster access
# de_long <- do.call(rbind, lapply(res_all, `[[`, "results"))

## ----eval=FALSE---------------------------------------------------------------
# res <- PseudobulkDE(
#   integrated,
#   sample_col    = "orig.ident",
#   condition_col = "treatment",
#   cluster_col   = "cell_type",
#   design        = ~ batch + treatment,
#   contrast      = c("treatment", "DrugA", "Vehicle")
# )

## ----eval=FALSE---------------------------------------------------------------
# PlotVolcano(res$results,
#             fc_threshold = 1,     # |log2FC| threshold
#             p_threshold  = 0.05,
#             top_n        = 20,    # label top-N by combined criterion
#             label_genes  = c("MYC", "TP53"))

## ----eval=FALSE---------------------------------------------------------------
# res_drug    <- PseudobulkDE(obj, sample_col = "orig.ident",
#                             condition_col = "treatment",
#                             ident_1 = "drug", ident_2 = "vehicle")
# res_disease <- PseudobulkDE(obj, sample_col = "orig.ident",
#                             condition_col = "status",
#                             ident_1 = "disease", ident_2 = "healthy")
# 
# cmp <- CompareMarkers(
#   res_drug$results, res_disease$results,
#   labels = c("drug_vs_veh", "disease_vs_healthy")
# )
# cmp$overlap                    # counts per category
# cmp$fisher                     # Fisher's exact test
# cmp$plot                       # log2FC-vs-log2FC scatter
# head(cmp$merged)

## ----eval=FALSE---------------------------------------------------------------
# comp <- CellComposition(
#   integrated,
#   cluster_col = "cell_type",
#   sample_col  = "orig.ident",
#   group_col   = "treatment",
#   style       = "box"
# )
# comp$plot
# head(comp$df)                  # per-sample per-cluster proportions

## ----eval=FALSE---------------------------------------------------------------
# comp_test <- CompositionalTest(
#   integrated,
#   cluster_col   = "cell_type",
#   sample_col    = "orig.ident",
#   condition_col = "treatment",
#   transform     = "asin",
#   method        = "auto"
# )
# subset(comp_test, padj < 0.05)

## ----eval=FALSE---------------------------------------------------------------
# comp <- CompositionAnalysis(
#   integrated,
#   group_col     = "cell_type",
#   sample_col    = "orig.ident",
#   condition_col = "treatment",
#   test          = "chisq"
# )
# comp$counts         # long tibble: sample, group, n_cells
# comp$proportions     # long tibble: sample, group, n_cells, prop, condition
# comp$test            # chi-square result, if test != "none" (pooled-cell caveat above)

## ----eval=FALSE---------------------------------------------------------------
# # One cell type
# CompositionEstimationPlot(comp, group_levels = "T cell",
#                           idx = c("Vehicle", "DrugA"))
# 
# # Every cell type present, as a named list of plots
# plots <- CompositionEstimationPlot(comp, idx = c("Vehicle", "DrugA"))
# plots[["T cell"]]
# 
# # Cohen's h -- the effect size designed specifically for comparing two proportions
# CompositionEstimationPlot(comp, idx = c("Vehicle", "DrugA"), effect = "cohens_h")

## ----eval=FALSE---------------------------------------------------------------
# Idents(integrated) <- integrated$cell_type
# lr <- RunLIANA(
#   integrated,
#   idents_col = "cell_type",
#   method     = "consensus",
#   resource   = "Consensus",
#   min_cells  = 10
# )
# head(lr)
# 
# # Restrict to interactions FROM T cells
# lr_t <- RunLIANA(integrated, idents_col = "cell_type",
#                  source_cells = "T cell")
# 
# # Mouse via ortholog translation
# lr_mouse <- RunLIANA(integrated, idents_col = "cell_type",
#                      use_ortho = TRUE)

## ----eval=FALSE---------------------------------------------------------------
# integrated <- PseudotimeWrapper(
#   integrated,
#   reduction     = "umap",
#   dims          = 2,
#   cluster_col   = "seurat_clusters",
#   start_cluster = "3",             # optional root
#   prefix        = "slingshot"
# )
# # New columns: slingshot_Lineage1, slingshot_Lineage2, ...
# FeaturePlot(integrated, features = "slingshot_Lineage1")
# 
# # The SlingshotDataSet itself, for downstream analysis
# sds <- integrated@misc$slingshot

## ----eval=FALSE---------------------------------------------------------------
# pre  <- BatchEffectQC(integrated, reduction = "pca",
#                       batch_col    = "orig.ident",
#                       celltype_col = "cell_type")
# post <- BatchEffectQC(integrated, reduction = "harmony",
#                       batch_col    = "orig.ident",
#                       celltype_col = "cell_type")
# 
# rbind(pre$summary, post$summary)
# #      batch_asw celltype_asw knn_mixing knn_purity ...
# # pre       0.32          0.28       0.61       0.72
# # post      0.04          0.35       0.94       0.78

## ----eval=FALSE---------------------------------------------------------------
# PlotFeatureDensity(
#   integrated,
#   features  = c("CD3D", "CD4", "CD8A"),
#   reduction = "umap"
# )
# 
# # Joint co-expression panel
# PlotFeatureDensity(integrated, features = c("CD3D", "CD4"), joint = TRUE)
# 
# # From a metadata column (e.g. a module score)
# PlotFeatureDensity(integrated, features = "module_score_Bcell")

## ----eval=FALSE---------------------------------------------------------------
# visium_dirs <- c(
#   "data/visium/sample1/outs",
#   "data/visium/sample2/outs"
# )
# 
# visium_list <- CreateVisiumObjects(
#   data_dirs    = visium_dirs,
#   object_names = c("S1", "S2")
# )

## ----eval=FALSE---------------------------------------------------------------
# visium_hd_list <- CreateVisiumObjects(
#   data_dirs    = visium_hd_dirs,
#   hd_bin_size  = "002um",
#   on_disk      = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# # Run on the first sample
# edge_df <- EdgeDetectionVisium(
#   coord_path = "data/visium/sample1/outs/spatial",  # folder containing tissue_positions_list.csv
#   seurat.obj = visium_list[["S1"]],   # ensures barcode order matches
#   search     = "radius",
#   neighbors  = 7
# )
# 
# # edge_df has columns: barcode, Filter, Filter2, Filter3, Filter4
# # Each Filter column represents one iteration; Filter4 is the most conservative.
# head(edge_df[, c("barcode", "Filter", "Filter4")])

## ----eval=FALSE---------------------------------------------------------------
# library(Seurat)
# 
# # Add the filter column to metadata
# visium_list[["S1"]] <- AddMetaData(
#   visium_list[["S1"]],
#   metadata = setNames(edge_df$Filter4, edge_df$barcode),
#   col.name = "edge_filter"
# )
# 
# # Visualize
# SpatialDimPlot(visium_list[["S1"]], group.by = "edge_filter")
# 
# # Remove edge spots
# visium_list[["S1"]] <- subset(visium_list[["S1"]], edge_filter == "Keep")

## ----eval=FALSE---------------------------------------------------------------
# integrated_visium <- MergeSeurat(
#   seurat_objects     = visium_list,
#   spatial            = "Visium",
#   use_SCT            = TRUE,
#   to_regress         = "percent.mt",
#   integration        = "HarmonyIntegration",
#   new_reduction      = "harmony",
#   cluster_resolution = 0.4,
#   markers            = FALSE
# )

## ----eval=FALSE---------------------------------------------------------------
# xenium_list <- lapply(
#   c("data/xenium/sample1", "data/xenium/sample2"),
#   LoadXenium2
# )
# names(xenium_list) <- c("X1", "X2")
# 
# # Merge and integrate
# integrated_xenium <- MergeSeurat(
#   seurat_objects = xenium_list,
#   spatial        = "Xenium",
#   integration    = "HarmonyIntegration",
#   new_reduction  = "harmony",
#   use_SCT        = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# xenium_obj <- LoadXenium2(
#   data_dir    = "data/xenium/sample1",
#   sample_name = "X1",
#   on_disk     = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# parse_list <- MakeParseObj(
#   data_dirs    = c("data/parse/run1/DGE_filtered",
#                    "data/parse/run2/DGE_filtered"),
#   mt_pattern   = "^mt-",
#   object_names = c("Run1", "Run2")
# )

## ----eval=FALSE---------------------------------------------------------------
# # Parse Biosciences: two DGE_filtered/ rounds per sample
# combined <- CombineParseRounds(
#   round1_paths = c("data/parse/run1/DGE_filtered/sample1",
#                    "data/parse/run1/DGE_filtered/sample2"),
#   round2_paths = c("data/parse/run2/DGE_filtered/sample1",
#                    "data/parse/run2/DGE_filtered/sample2"),
#   sample_names = c("sample1", "sample2"),
#   output_dir   = "data/parse/combined"
# )
# 
# # Then read the combined output as usual
# parse_list <- MakeParseObj(
#   data_dirs = file.path("data/parse/combined", c("sample1", "sample2"), "DGE_filtered")
# )

## ----eval=FALSE---------------------------------------------------------------
# # CellRanger: two filtered_feature_bc_matrix/ rounds per sample
# combined <- CombineCellRangerRounds(
#   round1_paths = c("data/cellranger/run1/sample1/outs/filtered_feature_bc_matrix",
#                    "data/cellranger/run1/sample2/outs/filtered_feature_bc_matrix"),
#   round2_paths = c("data/cellranger/run2/sample1/outs/filtered_feature_bc_matrix",
#                    "data/cellranger/run2/sample2/outs/filtered_feature_bc_matrix"),
#   sample_names = c("sample1", "sample2"),
#   output_dir   = "data/cellranger/combined"
# )
# 
# # Then read the combined output with CreateRNAObjects() as usual
# seurat_list <- CreateRNAObjects(
#   data_dirs = file.path("data/cellranger/combined", c("sample1", "sample2"),
#                         "filtered_feature_bc_matrix")
# )

## ----eval=FALSE---------------------------------------------------------------
# # Bounding-box method (default): fast, deterministic, good for roughly convex FOVs
# obj <- detect_fov_edges(obj,
#                         method       = "bbox",
#                         bbox_factor  = 2,        # ~2 cell-widths of the box edge
#                         n_iterations = 2,
#                         label_col    = "edge_layer")
# 
# # Angular-gap method: catches concave edges / tears the bbox misses
# obj <- detect_fov_edges(obj,
#                         method         = "angular",
#                         k              = 10,
#                         gap_threshold  = 2 * pi / 3,
#                         density_factor = 1)

## ----eval=FALSE---------------------------------------------------------------
# obj <- detect_tissue_holes2(obj,
#                             bin_size     = NULL,      # auto = 2.5 * median NN dist
#                             min_hole_size = 4,
#                             n_iterations  = 2,
#                             label_col     = "hole_layer")

## ----eval=FALSE---------------------------------------------------------------
# obj <- detect_tissue_holes2(obj,
#                             exclude_gene       = "Glul",
#                             sensitivity        = 0.75,   # data-adaptive quantile
#                             exclude_gene_layer = "data")

## ----eval=FALSE---------------------------------------------------------------
# obj <- obj[, obj$edge_layer == 0 & obj$hole_layer == 0]

## ----eval=FALSE---------------------------------------------------------------
# enrich <- NeighborhoodEnrichment(
#   obj,
#   group.by      = "cell_type",
#   k             = 10,
#   n_perm        = 200,
#   assign_niches = TRUE,     # k-means clustering of neighborhood vectors
#   n_niches      = 6
# )
# enrich$z                      # source x target enrichment z-score matrix
# enrich$results                 # same info as a long-form data frame (focal, neighbor, observed, expected, z, p, padj)
# obj <- enrich$obj              # a copy of `obj` with niche labels already written to meta.data$niche

## ----eval=FALSE---------------------------------------------------------------
# co <- NicheCoExpress(
#   obj,
#   genes         = c("Vegfa", "Kdr"),   # or a 2-column data.frame of specific pairs
#   niche_col     = "niche",
#   sample_col    = "sample",            # meta.data column identifying biological samples
#   condition_col = "condition",         # meta.data column with exactly 2 levels to compare
#   layer         = "data"
# )
# co$stats                      # per niche x pair: delta, statistic, p, p_adj (+ comp_diff/comp_flag if celltype_col set)
# co$per_sample                 # long-form per-sample x niche x pair co-expression scores
# 
# plotNicheCoExpress(co, type = "heatmap")   # niche x pair heatmap of delta, significance stars from p_adj
# plotNicheCoExpress(co, type = "scores")    # per-sample score plots for top/selected pairs

## ----eval=FALSE---------------------------------------------------------------
# # One niche x pair combination
# NicheCoExpressEstimationPlot(co, niches = "N1", pairs = "Vegfa_Kdr")
# 
# # Every niche x pair combination present, as a named list of plots
# plots <- NicheCoExpressEstimationPlot(co)
# plots[["N1 | Vegfa_Kdr"]]

## ----eval=FALSE---------------------------------------------------------------
# visium <- RunRCTD(
#   visium,
#   reference    = pbmc_ref,
#   celltype_col = "cell_type",
#   mode         = "full",         # or "doublet" / "multi"
#   max_cells_per_ref_celltype = 10000,
#   n_cores      = 8
# )
# 
# # Per-cell-type proportion columns
# colnames(visium@misc$rctd_weights)
# SpatialFeaturePlot(visium, features = c("rctd_T_cell", "rctd_B_cell"))
# 
# # Winner-takes-all label per spot (for quick visualization)
# SpatialDimPlot(visium, group.by = "rctd_dominant")

## ----eval=FALSE---------------------------------------------------------------
# res <- SummarizeRCTDByCluster(visium)
# res$labels        # named vector: cluster -> most likely cell type
# res$composition    # cluster x cell-type mean-proportion matrix

## ----eval=FALSE---------------------------------------------------------------
# atac_list <- CreateATACObjects(
#   data_dirs     = c("data/atac/sample1", "data/atac/sample2"),
#   add_treatment = TRUE,
#   treatment     = c("Vehicle", "DrugA"),
#   genome        = "mm10"   # or "hg38" for human data
# )
# names(atac_list)   # "sample1", "sample2" -- from basename(data_dirs)
# 
# # Custom names instead of basename(data_dirs)
# atac_list <- CreateATACObjects(
#   data_dirs    = c("data/atac/sample1", "data/atac/sample2"),
#   object_names = c("Vehicle_rep1", "DrugA_rep1")
# )

## ----eval=FALSE---------------------------------------------------------------
# atac_list <- CreateATACObjectsFilter(
#   data_dirs   = c("data/atac/sample1", "data/atac/sample2"),
#   interactive = TRUE,
#   genome      = "mm10"
# )

## ----eval=FALSE---------------------------------------------------------------
# # Export to AnnData (writes .h5ad on disk)
# ToAnnData(integrated, file = "results/pbmc.h5ad", assay = "RNA", layer = "data")
# 
# # Read AnnData back to a Seurat object
# obj <- FromAnnData(file = "data/tabula_sapiens/Lung.h5ad",
#                    reader = "python",   # "python" (via reticulate) or "R" (via zellkonverter)
#                    assay = "RNA")

## ----eval=FALSE---------------------------------------------------------------
# files <- list.files("data/tabula_sapiens", pattern = "\\.h5ad$",
#                     full.names = TRUE)
# objs  <- lapply(files, FromAnnData)
# names(objs) <- sub("\\.h5ad$", "", basename(files))

## ----eval=FALSE---------------------------------------------------------------
# SaveWithProvenance(
#   integrated,
#   file    = "results/pbmc_annotated.rds",
#   git_dir = getwd(),
#   extra   = list(analyst = "K. Evensen", project = "study42")
# )
# # Writes results/pbmc_annotated.rds AND
# #        results/pbmc_annotated_provenance.json

## ----eval=FALSE---------------------------------------------------------------
# s <- CellSuiteSummary(integrated,
#                       cluster_col = "cell_type",
#                       sample_col  = "orig.ident",
#                       top_markers = 5)
# print(s)

## ----eval=FALSE---------------------------------------------------------------
# # Check a single object
# detect_gene_id_type(seurat_list[["Veh1"]])
# #> [1] "MGI symbol"   # or "HGNC symbol", "Ensembl", "RefSeq", "Entrez"
# 
# # Check across an entire list — catches mismatches before they cost you
# check_gene_ids_across_objects(seurat_list)

## ----eval=FALSE---------------------------------------------------------------
# DetectGenes(integrated, genes = c("Sftpc", "Ager", "ENSMUG00000001234"))

## ----eval=FALSE---------------------------------------------------------------
# # Requires mouse S and G2M gene lists
# # The function uses built-in gene sets from Seurat by default
# integrated <- assign_cell_cycle_phase(integrated)
# 
# # New metadata columns: S.Score, G2M.Score, Phase
# table(integrated$Phase)
# #>  G1  G2M    S
# #> 4021  812 1023
# 
# # Regress out cell cycle in a downstream normalization
# integrated <- SCTransform(integrated, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))

## ----eval=FALSE---------------------------------------------------------------
# visium_list <- BuildMultipleNicheAssays(
#   seurat_list     = visium_list,
#   neighbors       = 6,     # spots to include in each neighborhood
#   k               = 8,     # number of niche clusters
#   assay           = "SCT"
# )
# 
# # A new "niche" assay and "niche_cluster" metadata column are added to each object
# SpatialDimPlot(visium_list[["S1"]], group.by = "niche_cluster")

## ----eval=FALSE---------------------------------------------------------------
# # Pull tissue coordinates across all images in an object
# coords <- get_all_coords(integrated_xenium)
# 
# # Parse polygon definitions (e.g., from a JSON or CSV exported from Xenium Explorer)
# polys <- parse_polygons("my_rois.json")
# 
# # Identify which cells fall inside a polygon
# cells_in_roi <- get_cells_in_polygon(
#   coords  = coords,
#   polygon = polys[["ROI_1"]]
# )
# 
# # Subset the object to those cells
# roi_obj <- subset_opt(integrated_xenium, cells = cells_in_roi)

## ----eval=FALSE---------------------------------------------------------------
# # Always prefer subset_opt() over subset() for spatial objects
# clean <- subset_opt(integrated_xenium, seurat_clusters %in% c("0", "1", "3"))

## ----eval=FALSE---------------------------------------------------------------
# # Per-vertex polygon coordinates for one FOV, joined with a metadata column
# poly <- get_polygon_coords(integrated_xenium, image = "fov1", meta_cols = "cell_type")
# 
# # Continuous or discrete feature, auto-detected
# PlotPolygons(poly, feature = "cell_type")
# 
# # It's a plain ggplot, fully chainable
# PlotPolygons(poly, feature = "cell_type") +
#   ggplot2::ggtitle("Cell types, fov1")

## ----eval=FALSE---------------------------------------------------------------
# bg    <- PlotPolygons(poly, background = "grey90")
# tcell <- PlotPolygons(subset(poly, cell_type == "T cell"),
#                       feature = "Cd3e", colors = c("white", "red"))
# bcell <- PlotPolygons(subset(poly, cell_type == "B cell"),
#                       feature = "Cd19", colors = c("white", "blue"))
# 
# overlay <- stack_polygons(bg, poly, first = TRUE) +
#   stack_polygons(tcell, poly) +
#   stack_polygons(bcell, poly)
# 
# legends <- patchwork::wrap_plots(collect_legend(tcell), collect_legend(bcell), ncol = 1)
# overlay + legends + patchwork::plot_layout(widths = c(1, 0.25))

## ----eval=FALSE---------------------------------------------------------------
# # Get all children of "response to oxidative stress" (GO:0006979)
# children <- get_all_children("GO:0006979")
# length(children)
# #> [1] 47
# 
# # Use with a GO annotation database to collect relevant genes
# library(org.Mm.eg.db)
# genes_in_terms <- AnnotationDbi::select(
#   org.Mm.eg.db,
#   keys    = children,
#   columns = c("SYMBOL"),
#   keytype = "GO"
# )

## ----eval=FALSE---------------------------------------------------------------
# sessionInfo()

