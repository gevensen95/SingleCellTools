#' Merge Seurat RNA Objects
#'
#' This function merges and integrates a list of Seurat objects. By default, it
#' assumes the objects are SCT normalized and will be integrated using Harmony.
#' It can integrate ATAC objects, as well. This function also save multiple
#' Seurat objects along the way and runs FindAllMarkers on the clusters.
#'
#' @param seurat_objects List of Seurat objects
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
           k_anchor = NULL, k_weight = NULL) {
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

      if (use_SCT){
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
      if (use_SCT){
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
      }
    }

    print('Running FindAllMarkers')
    markers <- Seurat::FindAllMarkers(obj, logfc.threshold = 1, only.pos = TRUE,
                              min.pct = 0.25)
    write.csv(markers, 'markers_all.csv')
    return(obj)
  }
