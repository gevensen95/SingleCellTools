#' Create Seurat RNA Objects
#'
#' This function creates multiple Seurat objects. It takes a list of directories
#' as input. In each directory, there should be the contents of the filtered
#' feature matrix folder, a filtered feature matrix folder, or a .h5 file. You
#' can directly give the output folder from cellranger into this function.
#' It will preferentially choose the filtered_feature_matrix folder.
#'
#'
#' @param data_dirs Path to directories containing matrix.mtx, features.tsv, and
#'  barcodes.tsv or .h5 files.
#' @param cells Features must be expressed in at least this many cells
#' @param features Cells must have at least this many features
#' @param treatment Treatment metadata column (e.g., Age, chemical, etc.)
#' @param object_names Names for the Seurat objects
#' @return A list of Seurat objects
#' @export

CreateRNAObjects <- function(data_dirs, cells = 3, features = 200,
                             treatment = NULL,
                             object_names = NULL) {
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

  seurat_objects <- lapply(seurat_objects, function(obj) {
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^mt-")
    return(obj)
  })

  # Add a to_regress to metadata to specify treatment
  if (is.null(treatment)==FALSE){
    seurat_objects <- setNames(lapply(seq_along(seurat_objects), function(i) {
      seurat_obj <- seurat_objects[[i]]
      seurat_obj[["Treatment"]] <- treatment[i]
      return(seurat_obj)
    }), names(seurat_objects))
  }

  obj <- merge(seurat_objects[[1]], seurat_objects[-1])
  gene.plot <- ggplot2::ggplot(obj@meta.data, aes(orig.ident, nFeature_RNA)) +
    ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered')
  mt.plot <- ggplot2::ggplot(obj@meta.data, aes(orig.ident, percent.mt)) +
    ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered')
  print(gene.plot + mt.plot)

  return(seurat_objects)
}
