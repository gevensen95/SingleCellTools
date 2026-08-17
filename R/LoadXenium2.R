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
#' @param type Cell spatial coordinate matrices to read -- the resulting FOV
#'   only gets the boundary type(s) actually requested here (requesting just
#'   one no longer errors, previously the FOV build unconditionally required
#'   both regardless of this argument).
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
#' @param on_disk Logical; if \code{TRUE}, move the returned object's Xenium
#'   counts layer to an on-disk BPCells matrix via
#'   \code{\link{ConvertToBPCells}} as the very last step. Requires the
#'   \code{BPCells} package (Suggests, not a hard dependency -- checked up
#'   front so this fails fast rather than after the full read). Default
#'   \code{FALSE}.
#' @param bpcells_dir Directory to write the on-disk matrix to when
#'   \code{on_disk = TRUE}. Default \code{file.path(getwd(), "bpcells",
#'   sample_name)}.
#' @return A list of filtered Seurat objects
#' @export

LoadXenium2 <- function(data_dir, sample_name,
                        outs = c("matrix", "microns"),
                        type = c("centroids", "segmentations"),
                        mols.qv.threshold = 20,
                        microns_lazy = FALSE,
                        on_disk = FALSE,
                        bpcells_dir = NULL)
  {
  type <- match.arg(arg = type, choices = c("centroids", "segmentations"),
                    several.ok = TRUE)
  outs <- match.arg(arg = outs, choices = c("matrix", "microns"),
                    several.ok = TRUE)
  outs <- c(outs, type)
  has_dt <- requireNamespace("data.table", quietly = TRUE) &&
    requireNamespace("R.utils", quietly = TRUE)

  if (isTRUE(on_disk) && !requireNamespace("BPCells", quietly = TRUE)) {
    stop("Package 'BPCells' is required for on_disk = TRUE. Install with: ",
         "remotes::install_github('bnprks/BPCells/r')")
  }
  if (is.null(bpcells_dir)) bpcells_dir <- file.path(getwd(), "bpcells", sample_name)

  # `arrow` is Suggests, not a hard dependency -- only actually needed when
  # "microns" (transcripts.parquet) was requested. Fail fast with a clear
  # message rather than letting arrow::read_parquet() below throw R's raw
  # "there is no package called 'arrow'" partway through the read.
  if ("microns" %in% outs && !requireNamespace("arrow", quietly = TRUE)) {
    stop("Package 'arrow' is required to read transcripts.parquet ",
         "(outs = 'microns'). install.packages('arrow')")
  }
  # CreateSeuratObject()/the ControlCodeword/ControlProbe assays below all
  # read from data$matrix unconditionally -- without this check, requesting
  # outs = "microns" alone (a reasonable reading of the docs, which describe
  # "matrix" and "microns" as independently requestable) would silently hand
  # CreateSeuratObject() counts = NULL and fail with a confusing error deep
  # inside Seurat instead of naming the actual problem here.
  if (!"matrix" %in% outs) {
    stop("`outs` must include \"matrix\" -- the counts matrix is required ",
         "to build the Seurat object (got outs = ", paste(outs, collapse = ", "), ").")
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
      matrix <- suppressWarnings(Seurat::Read10X(data.dir = file.path(data_dir,
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
      if (ncol(cell_boundaries_df) != 3) {
        stop("Expected cell_boundaries.csv.gz to have exactly 3 columns ",
             "(cell, x, y), got ", ncol(cell_boundaries_df),
             ". The 10x export format may have changed.")
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

  # Only build the boundary types actually requested via `type` -- `outs`
  # (and therefore `data`) only has entries for those, so unconditionally
  # building both here (as this used to) would pass NULL to
  # CreateCentroids()/CreateSegmentation() whenever the caller requested
  # just one of the two.
  boundaries <- list()
  fov_types  <- character(0)
  # CreateCentroids()/CreateSegmentation()/CreateFOV() are exported by
  # SeuratObject, not Seurat -- Seurat::CreateFOV errors with "'CreateFOV'
  # is not an exported object from 'namespace:Seurat'" (confirmed directly);
  # every fixture in this package's own test suite already builds these via
  # SeuratObject:: for the same reason.
  if ("centroids" %in% type) {
    boundaries$centroids <- SeuratObject::CreateCentroids(data$centroids)
    fov_types <- c(fov_types, "centroids")
  }
  if ("segmentations" %in% type) {
    boundaries$segmentation <- SeuratObject::CreateSegmentation(data$segmentations)
    fov_types <- c(fov_types, "segmentation")
  }
  message(sprintf('--- Building FOV (%s%s) ---',
                  paste(fov_types, collapse = " + "),
                  if ("microns" %in% outs) " + molecules" else ""))
  coords <- SeuratObject::CreateFOV(
    coords = boundaries,
    type = fov_types,
    molecules = data$microns,
    assay = "Xenium")

  message('--- Building Seurat object and attaching control assays ---')
  xenium.obj <- Seurat::CreateSeuratObject(
    counts = data$matrix[["Gene Expression"]],
    assay = "Xenium",
    project = sample_name)
  xenium.obj[["ControlCodeword"]] <- Seurat::CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium.obj[["ControlProbe"]] <- Seurat::CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  xenium.obj[["fov"]] <- coords

  if (!is.null(molecules_lazy_ds)) {
    message('  Attaching lazy arrow dataset connection at `obj@misc$molecules_lazy`')
    xenium.obj@misc$molecules_lazy <- molecules_lazy_ds
  }

  if (isTRUE(on_disk)) {
    xenium.obj <- ConvertToBPCells(xenium.obj, assay = "Xenium",
                                   layers = "counts", path = bpcells_dir)
  }

  return(xenium.obj)
}
