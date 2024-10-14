#' Create and Filter Seurat RNA Objects
#'
#' This function creates multiple Seurat objects and filters them. It takes a
#' list of directories as input. In each directory, there should be the contents
#' of the filtered feature matrix folder, a filtered feature matrix folder, or
#' a .h5 file. You can directly give the output folder from cellranger into this
#' function. It will preferentially choose the filtered_feature_matrix folder.
#' The objects are then automatically filtered based on the 15th
#' and 85th percentile for nFeature_RNA and percent.mt. This can function can
#' be run interactively to choose cutoffs without having to re-run it.
#'
#' @param data_dirs Path to the directories
#' @return A list of filtered Seurat objects
#' @export

CreateRNAObjectsFilter <-
  function(data_dirs, cells = 3, features = 200, treatment = NULL,
           use_quantile = TRUE, quantile_value_min = 0.15,
           feature_min = NA, feature_max = NA,
           percent_mt_max = NA, interactive = FALSE) {
    # Ensure thresholds are specified if not using quantiles
    if (!use_quantile) {
      if (!is.numeric(feature_min)) stop("Error: Did not specify threshold for feature_min")
      if (!is.numeric(feature_max)) stop("Error: Did not specify threshold for feature_max")
      if (!is.numeric(percent_mt_max)) stop("Error: Did not specify threshold for percent_mt_max")
    }

    # Use lapply to read the data and create Seurat objects
    seurat_objects <- lapply(data_dirs, function(dir) {
      # Read 10X data
      seurat_data <- Read10X(data.dir = dir)

      # Create Seurat object
      CreateSeuratObject(counts = seurat_data,
                         min.cells = cells,
                         min.features = features,
                         project = basename(dir))
    })


    # Name the list elements with the base names of the directories
    names(seurat_objects) <- basename(data_dirs)

    # Add percent mitochondrial DNA to each Seurat object
    seurat_objects <- lapply(seurat_objects, function(obj) {
      obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
      return(obj)
    })

    print('Seurat Objects Made')

    saveRDS(seurat_objects, 'seurat_objects_unfiltered.rds')

    print('Unfiltered Seurat Objects Saved')

    # Merge Seurat objects
    obj <- merge(seurat_objects[[1]], seurat_objects[-1])

    # Create plots
    gene.plot <- ggplot(obj@meta.data, aes(orig.ident, nFeature_RNA)) +
      geom_boxplot() + labs(title = 'Unfiltered')
    mt.plot <- ggplot(obj@meta.data, aes(orig.ident, percent.mt)) +
      geom_boxplot() + labs(title = 'Unfiltered')
    print(gene.plot + mt.plot)

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
      gene.plot <- ggplot(obj@meta.data, aes(orig.ident, nFeature_RNA)) +
        geom_boxplot() + labs(title = 'Filtered')
      mt.plot <- ggplot(obj@meta.data, aes(orig.ident, percent.mt)) +
        geom_boxplot() + labs(title = 'Filtered')
      print(gene.plot + mt.plot)

      print('Saving Filtered Objects')

      saveRDS(seurat_objects, 'seurat_objects_filtered.rds')

      return(seurat_objects)

    } else {
      # Non-interactive mode
      subsetted_objs <- lapply(seurat_objects, function(obj) {
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
      names(subsetted_objs) <- names(seurat_objects)

      obj <- merge(subsetted_objs[[1]], subsetted_objs[-1])
      gene.plot <- ggplot(obj@meta.data, aes(orig.ident, nFeature_RNA)) +
        geom_boxplot() + labs(title = 'Filtered')
      mt.plot <- ggplot(obj@meta.data, aes(orig.ident, percent.mt)) +
        geom_boxplot() + labs(title = 'Filtered')
      print(gene.plot + mt.plot)

      print('Saving Filtered Objects')

      saveRDS(subsetted_objs, 'seurat_objects_filtered.rds')

      return(subsetted_objs)
    }
  }
