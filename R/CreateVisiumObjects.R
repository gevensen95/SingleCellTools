#' Create Visium Spatial Objects
#'
#' This function takes creates Seurat objects of Visium data, of various
#' formats to make it easier if you download data from GEO
#'
#' @param data_dirs Path to the data directories
#' @param treatment Treatment variable for each sample
#' @param object_names Names of each object
#' @param file_type File type (e.g., h5 or directory)
#' @param mt_pattern Pattern for calculating percent mtDNA (\code{percent.mt}).
#'   Default \code{"^mt-"} (mouse gene symbol convention, e.g.
#'   \code{"mt-Nd1"}); pass \code{"^MT-"} for human tissue.
#' @param rb_pattern Pattern for calculating percent ribosomal-protein reads
#'   (\code{percent.rb}). Default \code{"^(Rp[sl]|RP[SL])"} matches both
#'   mouse and human gene symbol conventions.
#' @param hb_pattern Pattern for calculating percent hemoglobin reads
#'   (\code{percent.hb}) -- relevant for tissue with residual blood
#'   contamination. Default \code{"^(Hb[^p]|HB[^P])"} excludes the unrelated
#'   \code{"Hbp1"}/\code{"HBP1"} gene, a well-known false positive for naive
#'   \code{"^Hb"} patterns.
#' @return A list of Seurat Spatial objects
#' @export

CreateVisiumObjects <- function(data_dirs, treatment = NULL,
                                object_names = NULL, file_type = 'h5',
                                mt_pattern = '^mt-',
                                rb_pattern = '^(Rp[sl]|RP[SL])',
                                hb_pattern = '^(Hb[^p]|HB[^P])') {
  message(sprintf('--- Reading Visium data and creating Seurat objects (%d directories, file_type = %s) ---',
                  length(data_dirs), file_type))
  seurat_objects <- lapply(data_dirs, function(dir) {
    if (file_type == 'h5') {
      seurat_object <- Seurat::Read10X_h5(
        paste(dir,list.files(dir)[sapply(list.files(dir),
                                         function(x) all(c(grepl("filtered", x),
                                                           grepl(".h5", x))))], sep = '/'))
    } else if (file_type == 'directory') {
      seurat_object <- Seurat::Read10X(
        list.dirs(dir)[str_detect(list.dirs(dir), 'filtered')]
      )
    } else {stop('Choose file_type')}

    seurat_object <- CreateSeuratObject(counts = seurat_object, assay="Spatial",
                                        project = basename(dir))
    seurat_object[['RNA']] <- as(object = seurat_object[["Spatial"]],
                                 Class = "Assay5")
    DefaultAssay(seurat_object) <- 'RNA'
    seurat_object[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = mt_pattern)
    seurat_object[["percent.rb"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = rb_pattern)
    seurat_object[["percent.hb"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = hb_pattern)
    return(seurat_object)
  })

  # Name the list elements with the base names of the directories
  if (is.null(object_names) == TRUE) {
    names(seurat_objects) <- basename(data_dirs)
  } else {names(seurat_objects) <- object_names}
  # Add to metadata to specify treatment
  if (is.null(treatment)==FALSE){
    message('--- Adding Treatment metadata column ---')
    seurat_objects <- setNames(lapply(seq_along(seurat_objects), function(i) {
      seurat_obj <- seurat_objects[[i]]
      seurat_obj[["Treatment"]] <- treatment[i]
      return(seurat_obj)
    }), names(seurat_objects))
  }

  message('--- Generating unfiltered QC plots ---')
  obj <- merge(seurat_objects[[1]], seurat_objects[-1])
  gene.plot <- ggplot2::ggplot(obj@meta.data, aes(orig.ident, nFeature_Spatial)) +
    ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered') + Ol_Reliable()
  count.plot <- ggplot2::ggplot(obj@meta.data, aes(orig.ident, nCount_Spatial)) +
    ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered') + Ol_Reliable()
  mt.plot <- ggplot2::ggplot(obj@meta.data, aes(orig.ident, percent.mt)) +
    ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered') + Ol_Reliable()
  print(gene.plot + count.plot + mt.plot)

  message('--- Running edge detection and attaching tissue images ---')
  seurat_objects <- setNames(lapply(seq_along(seurat_objects), function(i) {
    message(sprintf('  Processing object %d of %d: %s',
                    i, length(seurat_objects), names(seurat_objects)[i]))
    obj <- seurat_objects[[i]]
    obj[["barcode"]] <- colnames(obj)
    path_seurat <- paste(names(seurat_objects[i]),
                         list.dirs(names(seurat_objects[i]), full.names = F)[stringr::str_detect(list.dirs(names(seurat_objects[i]), full.names = F),
                                                                                                 pattern = 'spatial')], sep = '/')
    detected <- EdgeDetectionVisium(path_seurat, obj)
    obj$Filter <- detected$Filter
    obj$Filter2 <- detected$Filter2
    #read image
    vis.image <- Read10X_Image(
      image.dir = path_seurat,
      filter.matrix = TRUE
    )
    vis.image@assay <- 'Spatial'
    vis.image@key <- 'slice1'
    obj@images$image <- vis.image
    return(obj)
  }), names(seurat_objects))

  return(seurat_objects)
}
