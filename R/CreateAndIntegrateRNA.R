#' Create and Integrate Seurat RNA Objects
#'
#' This function creates multiple Seurat objects and then integrates them. It
#' takes a list of directories as input. In each directory, there should be
#' the contents of the filtered feature matrix folder, a filtered feature
#' matrix folder, or a .h5 file. You can directly give the output folder from
#' cellranger into this function. It will preferentially choose the
#' filtered_feature_matrix folder. The function then automatically filters each
#' objects and integrates using Harmony. The user can change these parameters.
#' This function also save multiple Seurat objects along the way and
#' runs FindAllMarkers on the clusters.
#'
#' @param data_dirs Path to directories containing matrix.mtx, features.tsv, and
#'  barcodes.tsv or .h5 files.
#' @param cells Features must be expressed in at least this many cells
#' @param features Cells must have at least this many features
#' @param treatment Treatment metadata column (e.g., Age, chemical, etc.)
#' @param use_quantile Use quantile filtering method for nFeature_RNA and
#' percent.mt
#' @param quantile_value_min Minimum percentile for filtering (max will be 1-min)
#' @param feature_min Minimum cutoff for features
#' @param feature_max Maximum cutoff for features
#' @param percent_mt_max Maximum cutoff for percent.mt
#' @param interactive Run filtering interactively to pick cutoffs. A figure will
#' be shown of the the features and the percent.mt, or you can open the
#' automatically saved file
#' @param object_names Names for the Seurat objects
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
#' @return An integrated Seurat object
#' @export
CreateAndIntegrateRNA <-
  function(data_dirs, cells = 3, features = 200, treatment = NULL,
           use_quantile = TRUE, quantile_value_min = 0.15,
           feature_min = NA, feature_max = NA,
           percent_mt_max = NA, interactive = FALSE,
           object_names = NULL,
           cell_IDs = names(seurat_objects), to_regress = 'percent.mt',
           cluster_resolution = 0.3, max_dims = 15, use_SCT = TRUE,
           save_rds_file = TRUE, file_name = NULL,
           use_elbow_plot = FALSE, spatial = FALSE,
           integration = 'HarmonyIntegration',
           integration_normalization = 'SCT', integration_assay = 'SCT',
           integration_reduction = 'pca', new_reduction = 'harmony',
           k_anchor = NULL, k_weight = NULL) {
    # Ensure thresholds are specified if not using quantiles
    if (!use_quantile) {
      if (!is.numeric(feature_min)) stop("Error: Did not specify threshold for feature_min")
      if (!is.numeric(feature_max)) stop("Error: Did not specify threshold for feature_max")
      if (!is.numeric(percent_mt_max)) stop("Error: Did not specify threshold for percent_mt_max")
    }

    if (use_quantile) {
      if (is.numeric(feature_min)) stop("Error: Set quantile=TRUE and specificed hard cut off for min nFeature_RNA Pick only one.")
      if (is.numeric(feature_max)) stop("Error: Set quantile=TRUE and specificed hard cut off for max nFeature_RNA Pick only one.")
      if (is.numeric(percent_mt_max)) stop("Error: Set quantile=TRUE and specificed hard cut off for percent.mt. Pick only one.")
    }

    # Use lapply to read the data and create Seurat objects

    seurat_objects <- lapply(data_dirs, function(dir) {
      if (rlang::is_empty(
        list.files(dir, 'barcodes.tsv.gz|features.tsv.gz|matrix.mtx.gz')) == FALSE) {
        # Read 10X data
        seurat_data <- Seurat::Read10X(data.dir = dir)
        Seurat::CreateSeuratObject(counts = seurat_data,
                           min.cells = cells,
                           min.features = features,
                           project = basename(dir))
      } else if (rlang::is_empty(
        list.files(paste(dir,'filtered_feature_bc_matrix', sep = '/'),
                   'barcodes.tsv.gz|features.tsv.gz|matrix.mtx.gz')) == FALSE) {

        seurat_data <- Seurat::Read10X(data.dir =
                                         paste(dir,'filtered_feature_bc_matrix',
                                               sep = '/'))

        # Create Seurat object
        Seurat::CreateSeuratObject(counts = seurat_data,
                           min.cells = cells,
                           min.features = features,
                           project = basename(dir))
      } else if (rlang::is_empty(
        list.files(paste(paste(dir, 'outs', sep = '/'),
                         'filtered_feature_bc_matrix', sep = '/'),
                   'barcodes.tsv.gz|features.tsv.gz|matrix.mtx.gz')) == FALSE) {

        seurat_data <- Seurat::Read10X(data.dir =
                                         paste(paste(dir, 'outs', sep = '/'),
                                               'filtered_feature_bc_matrix', sep = '/'))

        # Create Seurat object
        Seurat::CreateSeuratObject(counts = seurat_data,
                                   min.cells = cells,
                                   min.features = features,
                                   project = basename(dir))
      } else if (sum(str_detect(list.files(dir), '.h5'))>0) {
        seurat_data <- Seurat::Read10X_h5(
          paste(dir,
                list.files(dir)[sapply(list.files(dir),
                                       function(x) all(c(grepl("filtered", x),
                                                         grepl(".h5", x))))], sep = '/'))

        # Create Seurat object
        Seurat::CreateSeuratObject(counts = seurat_data,
                           min.cells = cells,
                           min.features = features,
                           project = dirname(dir))
      }

    })
    # Name the list elements with the base names of the directories
    if (is.null(object_names) == TRUE) {
      names(seurat_objects) <- basename(data_dirs)
    } else {names(seurat_objects) <- object_names}

    # Add percent mitochondrial DNA to each Seurat object
    seurat_objects <- lapply(seurat_objects, function(obj) {
      obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^mt-")
      return(obj)
    })

    # Add a to_regress to metadata to specify treatment
    if (is.null(treatment) == FALSE){
      seurat_objects <- setNames(lapply(seq_along(seurat_objects), function(i) {
        seurat_obj <- seurat_objects[[i]]
        seurat_obj[["Treatment"]] <- treatment[i]
        return(seurat_obj)
      }), names(seurat_objects))
    }

    print('Seurat Objects Made')

    saveRDS(seurat_objects, file = 'seurat_objects_unfiltered.rds')

    print('Unfiltered Objects Saved')

    # Merge Seurat objects
    obj <- merge(seurat_objects[[1]], seurat_objects[-1])

    # Create plots
    gene.plot <- ggplot2::ggplot(obj@meta.data, aes(orig.ident, nFeature_RNA)) +
      ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered')
    mt.plot <- ggplot2::ggplot(obj@meta.data, aes(orig.ident, percent.mt)) +
      ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered')
    print(gene.plot + mt.plot)
    ggplot2::ggsave('unfiltered_features_percentMT.pdf', height = 5, width = 7)

    # Interactive mode
    if (interactive) {
      # Ask user for thresholds interactively
      for (param in c("min nFeature_RNA", "max nFeature_RNA", "max percent.mt")) {
        use_quantile <- readline(prompt = paste("Do you want to use quantile for subsetting", param, "? (yes/no): "))
        if (tolower(use_quantile) == "yes") {
          threshold <- as.numeric(readline(prompt = paste("Enter quantile threshold for", param, " (0 to 1): ")))
          seurat_objects <- lapply(seurat_objects, function(obj) {
            if (param == "min nFeature_RNA" ) {
              subset(obj, subset = nFeature_RNA > quantile(obj$nFeature_RNA, threshold))
            } else if (param == "max nFeature_RNA") {
              subset(obj, subset = nFeature_RNA < quantile(obj$nFeature_RNA, threshold))
            } else if (param == "max percent.mt") {
              subset(obj, subset = percent.mt < quantile(obj$percent.mt, threshold))
            }
          })
        } else if (tolower(use_quantile) == "no") {
          threshold <- as.numeric(readline(prompt = paste("Enter threshold for", param, ": ")))
          seurat_objects <- lapply(seurat_objects, function(obj) {
            if (param == "min nFeature_RNA") {
              subset(obj, subset = nFeature_RNA > threshold)
            } else if (param == "max nFeature_RNA") {
              subset(obj, subset = nFeature_RNA < threshold)
            } else if (param == "max percent.mt") {
              subset(obj, subset = percent.mt < threshold)
            }
          })
        }
      }

      obj <- merge(seurat_objects[[1]], seurat_objects[-1])
      gene.plot <- ggplot2::ggplot(obj@meta.data, aes(orig.ident, nFeature_RNA)) +
        ggplot2::geom_boxplot() + ggplot2::labs(title = 'Filtered')
      mt.plot <- ggplot2::ggplot(obj@meta.data, aes(orig.ident, percent.mt)) +
        ggplot2::geom_boxplot() + ggplot2::labs(title = 'Filtered')
      print(gene.plot + mt.plot)

    } else {
      # Non-interactive mode
      seurat_objects <- lapply(seurat_objects, function(obj) {
        if (use_quantile) {
          feature_threshold_min <- quantile(obj$nFeature_RNA, quantile_value_min)
          feature_threshold_max <- quantile(obj$nFeature_RNA, 1 - quantile_value_min)
          mt_threshold_max <- quantile(obj$percent.mt, 1 - quantile_value_min)

          subset(obj, subset = nFeature_RNA > feature_threshold_min &
                   nFeature_RNA < feature_threshold_max &
                   percent.mt < mt_threshold_max)
        } else {
          subset(obj, subset = nFeature_RNA > feature_min &
                   nFeature_RNA < feature_max &
                   percent.mt < percent_mt_max)
        }
      })

      # Name the subsetted list with original names for clarity
      names(seurat_objects) <- names(seurat_objects)

      obj <- merge(seurat_objects[[1]], seurat_objects[-1])
      gene.plot <- ggplot2::ggplot(obj@meta.data, aes(orig.ident, nFeature_RNA)) +
        ggplot2::geom_boxplot() + ggplot2::labs(title = 'Filtered')
      mt.plot <- ggplot2::ggplot(obj@meta.data, aes(orig.ident, percent.mt)) +
        ggplot2::geom_boxplot() + ggplot2::labs(title = 'Filtered')
      print(gene.plot + mt.plot)
    }

    print('Seurat Objects Filtered')

    saveRDS(seurat_objects, file = 'seurat_objects_filtered.rds')

    print('Filtered Objects Saved')

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
      print('Normalizing Data')
      obj <- Seurat::SCTransform(obj, vars.to.regress = to_regress)
    }
    else{
      print('Normalizing Data')
      obj <- Seurat::NormalizeData(obj)
      obj <- Seurat::FindVariableFeatures(obj)
      obj <- Seurat::ScaleData(obj, vars.to.regress = to_regress)
    }
    obj <- Seurat::RunPCA(obj)
    if (use_elbow_plot) {

      elbow_plot <- Seurat::ElbowPlot(obj)
      print(elbow_plot)

      pcs <- as.numeric(readline(prompt = 'Enter # of PCs: '))

      print('Integrating Data')
      obj <- Seurat::IntegrateLayers(object = obj, method = integration,
                             orig.reduction = integration_reduction,
                             assay = integration_assay,
                             normalization.method = integration_normalization,
                             new.reduction = new_reduction,
                             k.anchor = k_anchor, k.weight = k_weight)
      print('Running FindNeighbors')
      obj <- Seurat::FindNeighbors(obj, reduction = new_reduction, dims = 1:pcs)
      print('Running FindClusters')
      obj <- Seurat::FindClusters(obj, resolution = cluster_resolution)
      print('Running FindClusters')
      obj <- Seurat::FindClusters(obj, reduction = new_reduction, dims = 1:pcs)

      if (use_SCT){
        print('Running PrepSCTFindMarkers')
        obj <- Seurat::PrepSCTFindMarkers(obj)
      }

      if (save_rds_file == TRUE & is.null(file_name) == TRUE) {
        print('Saving Integrated Object')
        saveRDS(obj, paste(new_reduction, 'merged_seurat_objects.rds',
                           sep = '_'))
      } else if (save_rds_file == TRUE & is.null(file_name) == FALSE) {
        print('Saving Integrated Object')
        saveRDS(obj, paste(file_name, 'merged_seurat_objects.rds',
                           sep = '_'))
      }

    } else {
      print('Integrating Data')
      obj <- Seurat::IntegrateLayers(object = obj, method = integration,
                             orig.reduction = integration_reduction,
                             assay = integration_assay,
                             normalization.method = integration_normalization,
                             new.reduction = new_reduction,
                             k.anchor = k_anchor, k.weight = k_weight)
      print('Running FindNeighbors')
      obj <- Seurat::FindNeighbors(obj, reduction = new_reduction, dims = 1:max_dims)
      print('Running FindClusters')
      obj <- Seurat::FindClusters(obj, resolution = cluster_resolution)
      print('Running FindClusters')
      obj <- Seurat::FindClusters(obj, reduction = new_reduction, dims = 1:max_dims)
      if (use_SCT){
        print('Running PrepSCTFindMarkers')
        obj <- Seurat::PrepSCTFindMarkers(obj)
        if (save_rds_file == TRUE & is.null(file_name) == TRUE) {
          print('Saving Integrated Object')
          saveRDS(obj, paste(new_reduction, 'merged_seurat_objects.rds',
                             sep = '_'))
        } else if (save_rds_file == TRUE & is.null(file_name) == FALSE) {
          print('Saving Integrated Object')
          saveRDS(obj, paste(file_name, 'merged_seurat_objects.rds',
                             sep = '_'))
        }
        else {
          return(obj)
        }
      }
    }

    Seurat::DimPlot(obj, label = T)
    ggsave('plots/dimplot_seurat_clusters.pdf', height = 5, width = 7)

    print('Running FindAllMarkers')
    markers <- Seurat::FindAllMarkers(obj, logfc.threshold = 1, only.pos = TRUE,
                              min.pct = 0.25)
    write.csv(markers, 'markers_all.csv')
    return(obj)
  }
