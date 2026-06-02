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
#' @param mt_pattern Pattern for calculating percent mtDNA
#' @param treatment Treatment metadata column (e.g., Age, chemical, etc.)
#' @param object_names Names for the Seurat objects
#' @param run_doublet_finder Logical; if TRUE (default), run \code{calldoublet}
#'   on every object and add a \code{doublet_finder} metadata column.
#' @param doublet_normalization Passed to \code{calldoublet}: one of
#'   \code{"LogNormalize"} (default) or \code{"SCT"}.
#' @param doublet_vars_to_regress Passed to \code{calldoublet} as
#'   \code{vars.to.regress}. Default \code{"percent.mt"} because percent.mt
#'   has already been computed above; set to \code{NULL} to skip regression.
#' @param doublet_cluster_resolution Passed to \code{calldoublet} as
#'   \code{cluster_resolution}. Default \code{0.1}.
#' @param filter_doublets Logical; if TRUE, subset each object to
#'   \code{doublet_finder == "Singlet"} after doublet calling. Default
#'   \code{FALSE} so the doublet labels are preserved for downstream review.
#' @return A list of Seurat objects
#' @importFrom Seurat Read10X Read10X_h5 CreateSeuratObject PercentageFeatureSet
#' @import Seurat
#' @import stringr
#' @import dplyr
#' @import readr
#' @import ggplot2
#' @importFrom rlang is_empty
#' @importFrom stringr str_detect
#' @importFrom ggplot2 ggplot geom_boxplot labs theme element_text aes
#' @export

CreateRNAObjects <- function(data_dirs, cells = 3, features = 200,
                             mt_pattern = '^mt-',
                             treatment = NULL,
                             object_names = NULL,
                             run_doublet_finder = TRUE,
                             doublet_normalization = c("LogNormalize", "SCT"),
                             doublet_vars_to_regress = "percent.mt",
                             doublet_cluster_resolution = 0.1,
                             filter_doublets = FALSE) {
  doublet_normalization <- match.arg(doublet_normalization)

  message(sprintf('--- Reading data and creating Seurat objects (%d directories) ---',
                  length(data_dirs)))
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
    } else if (sum(stringr::str_detect(list.files(dir), '.h5')) > 0) {
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

  message('--- Calculating percent mitochondrial reads ---')
  seurat_objects <- lapply(seurat_objects, function(obj) {
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = mt_pattern)
    return(obj)
  })

  # Add a to_regress to metadata to specify treatment
  if (is.null(treatment)==FALSE){
    message('--- Adding Treatment metadata column ---')
    seurat_objects <- setNames(lapply(seq_along(seurat_objects), function(i) {
      seurat_obj <- seurat_objects[[i]]
      seurat_obj[["Treatment"]] <- treatment[i]
      return(seurat_obj)
    }), names(seurat_objects))
  }

  seurat_objects <- setNames(lapply(seq_along(seurat_objects), function(i) {
    seurat_obj[["RNA"]] <- methods::as(seurat_obj[["RNA"]], Class = "Assay5")
    return(seurat_obj)
  }), names(seurat_objects))

  # ---- Doublet detection --------------------------------------------------
  if (isTRUE(run_doublet_finder)) {
    message(sprintf('--- Calling doublets with DoubletFinder (%s) ---',
                    doublet_normalization))
    seurat_objects <- setNames(lapply(seq_along(seurat_objects), function(i) {
      lab <- names(seurat_objects)[i]
      message(sprintf('  [%d/%d] %s', i, length(seurat_objects), lab))
      out <- calldoublet(seurat_objects[[i]],
                         samplenameIndex    = i,
                         normalization      = doublet_normalization,
                         vars.to.regress    = doublet_vars_to_regress,
                         cluster_resolution = doublet_cluster_resolution)
      if (isTRUE(filter_doublets)) {
        n_before <- ncol(out)
        out      <- subset(out, doublet_finder == "Singlet")
        message(sprintf('    %s: dropped %d doublets (%d singlets remaining)',
                        lab, n_before - ncol(out), ncol(out)))
      }
      out
    }), names(seurat_objects))
  }

  message('--- Generating unfiltered QC plots ---')
  obj <- merge(seurat_objects[[1]], seurat_objects[-1])
  orig.ident <- nFeature_RNA <- percent.mt <- NULL  # silence R CMD check NSE notes
  gene.plot <- ggplot2::ggplot(obj@meta.data, ggplot2::aes(orig.ident, nFeature_RNA)) +
    ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered') +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  mt.plot <- ggplot2::ggplot(obj@meta.data, ggplot2::aes(orig.ident, percent.mt)) +
    ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered') +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  print(gene.plot + mt.plot)

  return(seurat_objects)
}
