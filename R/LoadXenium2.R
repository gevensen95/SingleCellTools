#' Read and Load 10X Genomics Xenium Updated
#'
#' This function is an updated version of LoadXenium from Seurat, as the
#' original function cannot use the new output of Xenium.
#'
#' @param data_dir Path to the directory containing all Xenium data
#' @param sample_name Name of Sample
#' @param outs Molecular outputs to read
#'  - "matrix": counts matrix
#'  - "microns": molecule coordinates -- reads transcripts.parquet via the
#'    \code{arrow} package (Suggests, not a hard dependency); errors with a
#'    clear message up front if \code{arrow} isn't installed and this is
#'    requested
#' @param type Cell spatial coordinate matrices to read
#'  - "centroids": cell centroids
#'  - "segmentations": cell segmentations
#' @param mols.qv.threshold Remvoe transcript molecules wiht a QV less than
#' specified threshold (20 is recommended)
#' @param microns_lazy Logical; if \code{FALSE} (default), \code{"microns"}
#'   reads the entire \code{transcripts.parquet} into memory before
#'   filtering by \code{mols.qv.threshold} -- unchanged from before. If
#'   \code{TRUE}, the read goes through \code{arrow::open_dataset()} with
#'   the QV filter and column selection pushed down to Arrow's query
#'   engine instead of materializing the full table first, and the
#'   unfiltered dataset connection is additionally attached at
#'   \code{obj@misc$molecules_lazy} so \code{\link{QueryXeniumMolecules}}
#'   can later pull a gene- or region-restricted subset without re-reading
#'   the file. For a whole-slide Xenium run (10^8+ transcript rows) this is
#'   the difference between one big in-memory read and only ever touching
#'   the rows you actually need.
#' @return A list of filtered Seurat objects
#' @export

LoadXenium2 <- function(data_dir, sample_name,
                        outs = c("matrix", "microns"),
                        type = c("centroids", "segmentations"),
                        mols.qv.threshold = 20,
                        microns_lazy = FALSE)
  {
  type <- match.arg(arg = type, choices = c("centroids", "segmentations"),
                    several.ok = TRUE)
  outs <- match.arg(arg = outs, choices = c("matrix", "microns"),
                    several.ok = TRUE)
  outs <- c(outs, type)
  has_dt <- requireNamespace("data.table", quietly = TRUE) &&
    requireNamespace("R.utils", quietly = TRUE)

  # `arrow` is Suggests, not a hard dependency -- only actually needed when
  # "microns" (transcripts.parquet) was requested. Fail fast with a clear
  # message rather than letting arrow::read_parquet() below throw R's raw
  # "there is no package called 'arrow'" partway through the read.
  if ("microns" %in% outs && !requireNamespace("arrow", quietly = TRUE)) {
    stop("Package 'arrow' is required to read transcripts.parquet ",
         "(outs = 'microns'). install.packages('arrow')")
  }

  message(sprintf('--- Loading Xenium sample "%s" from %s ---', sample_name, data_dir))
  message(sprintf('  Outputs requested: %s', paste(outs, collapse = ', ')))

  # Populated by the "microns" branch below when microns_lazy = TRUE, so it
  # can be attached to the returned Seurat object after construction. Left
  # NULL otherwise (eager mode, or "microns" not requested at all).
  molecules_lazy_ds <- NULL

  data <- sapply(outs, function(otype) {
    switch(EXPR = otype, matrix = {
      message('  Reading counts matrix (cell_feature_matrix/)')
      matrix <- suppressWarnings(Read10X(data.dir = file.path(data_dir,
                                                              "cell_feature_matrix/")))
      matrix
    }, centroids = {
      message('  Reading cell centroids (cells.csv.gz)')
      if (has_dt) {
        cell_info <- as.data.frame(data.table::fread(file.path(data_dir,
                                                               "cells.csv.gz")))
      } else {
        cell_info <- read.csv(file.path(data_dir, "cells.csv.gz"))
      }
      cell_centroid_df <- data.frame(x = cell_info$x_centroid,
                                     y = cell_info$y_centroid, cell = cell_info$cell_id,
                                     stringsAsFactors = FALSE)
      cell_centroid_df
    }, segmentations = {
      message('  Reading cell segmentations (cell_boundaries.csv.gz)')
      if (has_dt) {
        cell_boundaries_df <- as.data.frame(data.table::fread(file.path(data_dir,
                                                                        "cell_boundaries.csv.gz")))
      } else {
        cell_boundaries_df <- read.csv(file.path(data_dir,
                                                 "cell_boundaries.csv.gz"), stringsAsFactors = FALSE)
      }
      names(cell_boundaries_df) <- c("cell", "x", "y")
      cell_boundaries_df
    }, microns = {
      parquet_path <- file.path(data_dir, "transcripts.parquet")
      if (isTRUE(microns_lazy)) {
        message(sprintf(
          '  Reading transcripts (transcripts.parquet, qv >= %g) -- lazy via arrow::open_dataset()',
          mols.qv.threshold))
        # Query pushdown: the filter/select happen in Arrow's engine, so only
        # the matching rows/columns are ever materialized into R. The
        # unfiltered dataset itself is kept (via the <<- below) for later
        # windowed/gene-subset queries -- see QueryXeniumMolecules().
        ds <- arrow::open_dataset(parquet_path, format = "parquet")
        molecules_lazy_ds <<- ds
        filtered <- dplyr::filter(ds, qv >= mols.qv.threshold)
        filtered <- dplyr::select(filtered, x = x_location, y = y_location,
                                  gene = feature_name)
        df <- as.data.frame(dplyr::collect(filtered))
      } else {
        message(sprintf('  Reading transcripts (transcripts.parquet, qv >= %g)',
                        mols.qv.threshold))
        transcripts <- arrow::read_parquet(parquet_path)
        transcripts <- subset(transcripts, qv >= mols.qv.threshold)

        df <- data.frame(x = transcripts$x_location, y = transcripts$y_location,
                         gene = transcripts$feature_name, stringsAsFactors = FALSE)
      }
      df
    }, stop("Unknown Xenium input type: ", otype))
  }, USE.NAMES = TRUE)

  message('--- Building FOV (centroids + segmentations + molecules) ---')
  segmentations.data <- list(
    centroids = CreateCentroids(data$centroids),
    segmentation = CreateSegmentation(data$segmentations))
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = data$microns,
    assay = "Xenium")

  message('--- Building Seurat object and attaching control assays ---')
  xenium.obj <- CreateSeuratObject(
    counts = data$matrix[["Gene Expression"]],
    assay = "Xenium",
    project = sample_name)
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  xenium.obj[["fov"]] <- coords

  if (!is.null(molecules_lazy_ds)) {
    message('  Attaching lazy arrow dataset connection at `obj@misc$molecules_lazy`')
    xenium.obj@misc$molecules_lazy <- molecules_lazy_ds
  }

  return(xenium.obj)
}
