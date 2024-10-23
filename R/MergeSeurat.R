#' Merge Seurat RNA Objects
#'
#' This function merges and integrates a list of Seurat objects. By default, it
#' assumes the objects are SCT normalized and will be integrated using Harmony.
#' It can integrate ATAC objects, as well. This function also save multiple
#' Seurat objects along the way and runs FindAllMarkers on the clusters.
#'
#' @param seurat_objects List of Seurat objects
#' #' @param object_names Names for the Seurat objects
#' @param cell_IDs Adds cell identities to the beginning of each barcode
#' @param to_regress Variable to regress out in normalization
#' @param cluster_resolution Resolution for cell clustering
#' @param max_dims Maximum dimensions for FindNeighbors and UMAP
#' @param use_SCT Use SCTransform method or LogNormalize (FALSE)
#' @param save_rds_file Save final Seurat object
#' @param file_name Name of files
#' @param use_elbow_plot Use the ElbowPlot to determine number of componenets
#' for FindNeighbors and UMAP
#' @param spatial If your data is a spatial (i.e., Visium)
#' @param integration Method for integrating data, see IntegrateLayers
#' @param integration_normalization Normalization method used
#' @param integration_assay Assay that has normalized data
#' @param integration_reduction Reduction to use for integration
#' @param new_reduction Name of new integrated reduction
#' @param k_anchor How many neighbors (k) to use when picking anchors
#' @param k_weight Number of neighbors to consider when weighting anchors
#' @param markers Find all markers
#' @return An integrated Seurat object.
#' @export

MergeSeurat <-
  function(seurat_objects, cell_IDs = names(seurat_objects), to_regress = 'percent.mt',
           cluster_resolution = 0.3, max_dims = 15, use_SCT = TRUE,
           save_rds_file = TRUE, file_name = NULL,
           use_elbow_plot = FALSE, spatial = FALSE,
           integration = 'HarmonyIntegration',
           integration_normalization = 'SCT', integration_assay = 'SCT',
           integration_reduction = 'pca', new_reduction = 'harmony',
           k_anchor = NULL, k_weight = NULL,
           markers = TRUE) {
    if(integration != 'HarmonyIntegration' & new_reduction == 'harmony') {
      stop('\n\n  Error: Integration method is not the default (HarmonyIntegration).\nChange new_reduction to match integration method')
    }

    if(integration == 'RPCAIntegration' |
       integration == 'CCAIntegration' |
       integration == 'JointPCAIntegration'  & is.null(k_anchor) == TRUE) {
      stop('\n\n  Error: Integration method is RPCAIntegration.\nSpecifiy k_anchor to an integar value (recommend 20)')
    }

    if(integration == 'RPCAIntegration' |
       integration == 'CCAIntegration' |
       integration == 'JointPCAIntegration'  & is.null(k_weight) == TRUE) {
      stop('\n\n  Error: Integration method is RPCAIntegration.\nSpecifiy k_weight to an integar value (recommend 100)')
    }

    if (spatial){
      obj <- suppressWarnings(merge(seurat_objects[[1]], seurat_objects[-1],
                                    add.cell.ids = cell_IDs))
    } else {
      obj <- merge(seurat_objects[[1]], seurat_objects[-1],
                   add.cell.ids = cell_IDs)
    }

    if (use_SCT){
      obj <- Seurat::SCTransform(obj, vars.to.regress = to_regress)
    }
    else{
      obj <- Seurat::NormalizeData(obj)
      obj <- Seurat::FindVariableFeatures(obj)
      obj <- Seurat::ScaleData(obj, vars.to.regress = to_regress)
    }
    obj <- Seurat::RunPCA(obj)
    if (use_elbow_plot) {

      elbow_plot <- Seurat::ElbowPlot(obj)
      print(elbow_plot)

      pcs <- as.numeric(readline(prompt = 'Enter # of PCs: '))

      obj <- Seurat::IntegrateLayers(object = obj, method = integration,
                             orig.reduction = integration_reduction,
                             assay = integration_assay,
                             normalization.method = integration_normalization,
                             new.reduction = new_reduction,
                             k.anchor = k_anchor, k.weight = k_weight)
      obj <- Seurat::FindNeighbors(obj, reduction = new_reduction, dims = 1:pcs)
      obj <- Seurat::FindClusters(obj, resolution = cluster_resolution)
      obj <- Seurat::RunUMAP(obj, reduction = new_reduction, dims = 1:pcs)

      if (use_SCT == TRUE & markers == TRUE){
        obj <- Seurat::PrepSCTFindMarkers(obj)
      }

      if (save_rds_file == TRUE & is.null(file_name) == TRUE) {
        saveRDS(obj, paste(new_reduction, 'merged_seurat_objects.rds',
                           sep = '_'))
      } else if (save_rds_file == TRUE & is.null(file_name) == FALSE) {
        saveRDS(obj, paste(file_name, 'merged_seurat_objects.rds',
                           sep = '_'))
      }

    } else {
      obj <- Seurat::IntegrateLayers(object = obj, method = integration,
                             orig.reduction = integration_reduction,
                             assay = integration_assay,
                             normalization.method = integration_normalization,
                             new.reduction = new_reduction,
                             k.anchor = k_anchor, k.weight = k_weight)
      obj <- Seurat::FindNeighbors(obj, reduction = new_reduction, dims = 1:max_dims)
      obj <- Seurat::FindClusters(obj, resolution = cluster_resolution)
      obj <- Seurat::RunUMAP(obj, reduction = new_reduction, dims = 1:max_dims)
      if (use_SCT == TRUE & markers == TRUE){
        obj <- Seurat::PrepSCTFindMarkers(obj)
        if (save_rds_file == TRUE & is.null(file_name) == TRUE) {
          saveRDS(obj, paste(new_reduction, 'merged_seurat_objects.rds',
                             sep = '_'))
        } else if (save_rds_file == TRUE & is.null(file_name) == FALSE) {
          saveRDS(obj, paste(file_name, 'merged_seurat_objects.rds',
                             sep = '_'))
        }
        else {
          return(obj)
        }
      } else {
        if (save_rds_file == TRUE & is.null(file_name) == TRUE) {
          saveRDS(obj, paste(new_reduction, 'merged_seurat_objects.rds',
                             sep = '_'))
        } else if (save_rds_file == TRUE & is.null(file_name) == FALSE) {
          saveRDS(obj, paste(file_name, 'merged_seurat_objects.rds',
                             sep = '_'))
        }
      }
    }

    Seurat::DimPlot(obj, label = T)
    ggsave('dimplot_seurat_clusters.pdf', height = 5, width = 7)

    if (markers == TRUE){
        print('Running FindAllMarkers')
    markers <- Seurat::FindAllMarkers(obj, logfc.threshold = 1, only.pos = TRUE,
                              min.pct = 0.25)
    write.csv(markers, 'markers_all.csv')
    return(obj)
      } else {return(obj)}
  }
