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
#' @param image_backend One of \code{"eager"} (default) or \code{"deferred"}.
#'   \code{"eager"} attaches \code{tissue_hires_image.png} to every sample,
#'   exactly as before (a hires PNG decodes to roughly 100MB in memory --
#'   fine for one sample, but a list of N samples means N times that before
#'   you've done anything with them). \code{"deferred"} instead attaches the
#'   much smaller \code{tissue_lowres_image.png} (~1MB decoded) to every
#'   sample -- coordinates, scale factors, and everything
#'   \code{\link{EdgeDetectionVisium}}/\code{\link{GenerateQCReport}} use are
#'   identical either way, since those only read spot coordinates, never
#'   pixels -- and stashes the hires image's path at
#'   \code{obj@misc$hires_image_path} so \code{\link{GetHiresVisiumImage}}
#'   can decode it later for just the sample(s) you actually need full
#'   detail on. Every sample (regardless of \code{image_backend}) also
#'   gets \code{obj@misc$visium_image_dir} stashed, which
#'   \code{\link{DropSpatialImage}} and \code{\link{SpatialObjectInfo}} use.
#' @param hd_bin_size For Visium HD samples only (auto-detected per sample
#'   by the presence of a \code{binned_outputs/} subdirectory -- regular
#'   Visium samples are unaffected). One of \code{"008um"} (default -- 10x's
#'   commonly recommended bin size, roughly cell-scale), \code{"002um"}
#'   (finest resolution, largest data), or \code{"016um"} (coarsest,
#'   smallest). Selects \code{binned_outputs/square_<hd_bin_size>/} as the
#'   effective sample directory for everything (matrix, image, tissue
#'   positions). \code{binned_outputs.tar.gz}/\code{spatial.tar.gz} etc.
#'   must already be extracted -- this function doesn't extract archives
#'   itself, matching how it already expects regular Visium output
#'   pre-extracted.
#' @param workers Number of parallel workers to use (via \code{future.apply})
#'   for building each sample -- reading the matrix, computing QC metrics,
#'   running \code{\link{EdgeDetectionVisium}}, and attaching the tissue
#'   image, all fully independent across samples. Default \code{1} runs
#'   sequentially exactly as before (with per-sample progress messages);
#'   \code{workers > 1} spins up that many background R sessions via
#'   \code{future::plan(multisession)}, restored on exit. Note each worker
#'   holds its own copy of that sample's data, so peak memory scales with
#'   \code{workers}.
#' @return A list of Seurat Spatial objects
#' @export

CreateVisiumObjects <- function(data_dirs, treatment = NULL,
                                object_names = NULL, file_type = 'h5',
                                mt_pattern = '^mt-',
                                rb_pattern = '^(Rp[sl]|RP[SL])',
                                hb_pattern = '^(Hb[^p]|HB[^P])',
                                image_backend = c("eager", "deferred"),
                                hd_bin_size = c("008um", "002um", "016um"),
                                workers = 1) {
  image_backend <- match.arg(image_backend)
  hd_bin_size <- match.arg(hd_bin_size)
  if (!file_type %in% c('h5', 'directory')) {
    stop("Choose file_type: must be one of 'h5' or 'directory' (got: '", file_type, "').")
  }

  if (workers > 1) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Package 'future.apply' is required for workers > 1. ",
           "install.packages('future.apply')")
    }
    old_plan <- future::plan(future::multisession, workers = workers)
    on.exit(future::plan(old_plan), add = TRUE)
  }

  # A Visium HD sample nests everything (matrix, image, tissue positions)
  # under binned_outputs/square_<hd_bin_size>um/ instead of directly in the
  # sample directory the way regular Visium does. Detected per sample (via
  # a binned_outputs/ subdirectory) so a single call can mix HD and regular
  # Visium samples. Returns `dir` unchanged for regular Visium samples.
  .resolve_visium_sample_dir <- function(dir) {
    hd_root <- file.path(dir, "binned_outputs")
    if (!dir.exists(hd_root)) {
      return(dir)
    }
    hd_dir <- file.path(hd_root, paste0("square_", hd_bin_size))
    if (!dir.exists(hd_dir)) {
      available <- list.dirs(hd_root, full.names = FALSE, recursive = FALSE)
      stop("Detected a Visium HD directory ('", hd_root, "' exists) but '",
           basename(hd_dir), "' was not found there. Available bin size(s): ",
           paste(available, collapse = ", "), ". Set `hd_bin_size` to one of these.")
    }
    hd_dir
  }

  sample_names <- if (is.null(object_names)) basename(data_dirs) else object_names
  n_samples <- length(data_dirs)
  if (!is.null(treatment)) message('--- Treatment metadata column will be added per sample ---')

  message(sprintf('--- Reading Visium data and creating Seurat objects (%d directories, file_type = %s)%s ---',
                  n_samples, file_type,
                  if (workers > 1) sprintf(', %d parallel workers', workers) else ''))
  # Everything for one sample -- read, QC metrics, treatment, edge detection,
  # tissue image -- lives in a single per-sample pass now (previously two
  # separate lapply()s over the same samples: one ending in the QC-metrics
  # step, a second doing edge-detection/image-attachment). Combining them
  # removes a full extra pass over the sample list and, combined with
  # switching the QC plots below to a metadata-only bind_rows instead of a
  # full Seurat::merge(), means nothing here ever needs every sample's full
  # object in memory at once except the returned list itself.
  .build_one <- function(dir, i, nm) {
    if (workers == 1) {
      message(sprintf('  [%d/%d] %s', i, n_samples, nm))
    }
    eff_dir <- .resolve_visium_sample_dir(dir)
    if (file_type == 'h5') {
      seurat_data <- Seurat::Read10X_h5(
        paste(eff_dir, list.files(eff_dir)[sapply(list.files(eff_dir),
                                         function(x) all(c(grepl("filtered", x),
                                                           grepl(".h5", x))))], sep = '/'))
    } else {
      seurat_data <- Seurat::Read10X(
        list.dirs(eff_dir)[str_detect(list.dirs(eff_dir), 'filtered')]
      )
    }

    seurat_object <- CreateSeuratObject(counts = seurat_data, assay = "Spatial",
                                        project = basename(dir))
    seurat_object[['RNA']] <- as(object = seurat_object[["Spatial"]],
                                 Class = "Assay5")
    DefaultAssay(seurat_object) <- 'RNA'
    seurat_object[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = mt_pattern)
    seurat_object[["percent.rb"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = rb_pattern)
    seurat_object[["percent.hb"]] <- Seurat::PercentageFeatureSet(seurat_object, pattern = hb_pattern)

    if (!is.null(treatment)) {
      seurat_object[["Treatment"]] <- treatment[i]
    }

    seurat_object[["barcode"]] <- colnames(seurat_object)
    # Resolve against the actual dir path (not sample_names, which is
    # basename(data_dirs)/object_names and isn't necessarily a valid path
    # on its own -- e.g. when data_dirs holds full/relative paths whose
    # basename differs from a path relative to the working directory) and
    # through the same Visium HD redirection used for the matrix read above.
    path_seurat <- file.path(eff_dir,
                             list.dirs(eff_dir, full.names = FALSE)[stringr::str_detect(
                               list.dirs(eff_dir, full.names = FALSE), pattern = 'spatial')])
    # NB: seurat.obj/coord_path here, matching EdgeDetectionVisium()'s actual
    # signature -- `seurat_object` has no @images entries yet at this point
    # (attached further down), so this always falls through to the
    # coord_path/file read.
    detected <- EdgeDetectionVisium(seurat.obj = seurat_object, coord_path = path_seurat)
    seurat_object$Filter <- detected$Filter
    seurat_object$Filter2 <- detected$Filter2

    image_name <- if (image_backend == "deferred") {
      "tissue_lowres_image.png"
    } else {
      "tissue_hires_image.png"
    }
    vis.image <- Read10X_Image(
      image.dir  = path_seurat,
      image.name = image_name,
      filter.matrix = TRUE
    )
    vis.image@assay <- 'Spatial'
    vis.image@key <- 'slice1'
    seurat_object@images$image <- vis.image
    # Stashed unconditionally (not just under image_backend = "deferred") so
    # DropSpatialImage(mode = "downgrade") can rebuild a lowres image later
    # for an object that was originally loaded eager, and SpatialObjectInfo()
    # can report where each image actually came from.
    seurat_object@misc$visium_image_dir <- path_seurat
    if (image_backend == "deferred") {
      seurat_object@misc$hires_image_path <- file.path(path_seurat, "tissue_hires_image.png")
    }

    seurat_object
  }

  # NB: mapply/future_mapply zips directly over data_dirs/indices/names
  # rather than indexing into a captured list from inside the worker
  # function, so each worker only ever receives the one sample it's
  # processing (same reasoning as CreateRNAObjects()'s workers support).
  seurat_objects <- if (workers > 1) {
    future.apply::future_mapply(
      .build_one, data_dirs, seq_len(n_samples), sample_names,
      SIMPLIFY = FALSE, future.seed = TRUE
    )
  } else {
    mapply(.build_one, data_dirs, seq_len(n_samples), sample_names, SIMPLIFY = FALSE)
  }
  seurat_objects <- setNames(seurat_objects, sample_names)

  message('--- Generating unfiltered QC plots ---')
  # Row-bind just the metadata rather than Seurat::merge()-ing every sample's
  # full object -- merge() would also combine the count matrices and images,
  # which these QC plots never touch, so it's wasted work (and memory) for
  # large/many-sample lists. Same fix already applied in CreateRNAObjects().
  meta <- dplyr::bind_rows(lapply(seurat_objects, function(x) x@meta.data))
  orig.ident <- nFeature_Spatial <- nCount_Spatial <- percent.mt <- NULL  # silence R CMD check NSE notes
  gene.plot <- ggplot2::ggplot(meta, aes(orig.ident, nFeature_Spatial)) +
    ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered') + Ol_Reliable()
  count.plot <- ggplot2::ggplot(meta, aes(orig.ident, nCount_Spatial)) +
    ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered') + Ol_Reliable()
  mt.plot <- ggplot2::ggplot(meta, aes(orig.ident, percent.mt)) +
    ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered') + Ol_Reliable()
  print(gene.plot + count.plot + mt.plot)

  return(seurat_objects)
}
