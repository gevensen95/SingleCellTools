#' Read and Load 10X Genomics Xenium Updated
#'
#' This function is an updated version of LoadXenium from Seurat, as the
#' original function cannot use the new output of Xenium.
#'
#' @param data_dir Path to the directory containing all Xenium data
#' @param outs Molecular outputs to read
#'  - "matrix": counts matrix
#'  - "microns": molecule coordinates
#' @param type Cell spatial coordinate matrices to read
#'  - "centroids": cell centroids
#'  - "segmentations": cell segmentations
#' @param mols.qv.threshold Remvoe transcript molecules wiht a QV less than
#' specified threshold (20 is recommended)
#' @return A list of filtered Seurat objects
#' @export

LoadXenium2 <- function(data_dir, outs = c("matrix", "microns"),
                        type = c("centroids", "segmentations"),
                        mols.qv.threshold = 20)
  {
  type <- match.arg(arg = type, choices = c("centroids", "segmentations"),
                    several.ok = TRUE)
  outs <- match.arg(arg = outs, choices = c("matrix", "microns"),
                    several.ok = TRUE)
  outs <- c(outs, type)
  has_dt <- requireNamespace("data.table", quietly = TRUE) &&
    requireNamespace("R.utils", quietly = TRUE)
  data <- sapply(outs, function(otype) {
    switch(EXPR = otype, matrix = {
      matrix <- suppressWarnings(Read10X(data.dir = file.path(data.dir,
                                                              "cell_feature_matrix/")))
      matrix
    }, centroids = {
      if (has_dt) {
        cell_info <- as.data.frame(data.table::fread(file.path(data.dir,
                                                               "cells.csv.gz")))
      } else {
        cell_info <- read.csv(file.path(data.dir, "cells.csv.gz"))
      }
      cell_centroid_df <- data.frame(x = cell_info$x_centroid,
                                     y = cell_info$y_centroid, cell = cell_info$cell_id,
                                     stringsAsFactors = FALSE)
      cell_centroid_df
    }, segmentations = {
      if (has_dt) {
        cell_boundaries_df <- as.data.frame(data.table::fread(file.path(data.dir,
                                                                        "cell_boundaries.csv.gz")))
      } else {
        cell_boundaries_df <- read.csv(file.path(data.dir,
                                                 "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
      }
      names(cell_boundaries_df) <- c("cell", "x", "y")
      cell_boundaries_df
    }, microns = {

      transcripts <- arrow::read_parquet(file.path(data.dir, "transcripts.parquet"))
      transcripts <- subset(transcripts, qv >= mols.qv.threshold)

      df <- data.frame(x = transcripts$x_location, y = transcripts$y_location,
                       gene = transcripts$feature_name, stringsAsFactors = FALSE)
      df
    }, stop("Unknown Xenium input type: ", otype))
  }, USE.NAMES = TRUE)

  segmentations.data <- list(
    centroids = CreateCentroids(data$centroids),
    segmentation = CreateSegmentation(data$segmentations))
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = data$microns,
    assay = "Xenium")
  xenium.obj <- CreateSeuratObject(
    counts = data$matrix[["Gene Expression"]],
    assay = "Xenium")
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  xenium.obj[["fov"]] <- coords

  return(xenium.obj)
}
