# Locate a barcodes/features(or genes)/matrix triplet in `dir_path` and read
# it with Seurat::ReadMtx(), which (unlike Read10X()) takes explicit file
# paths rather than assuming the canonical un-prefixed CellRanger names.
# This lets each file be named e.g. "GSM9352918_con11_barcodes.tsv.gz" --
# any sample-name prefix in front of the standard barcodes/features/genes/
# matrix suffix is matched -- which is how GEO commonly redistributes 10X
# Supplementary Files. Returns NULL (so callers can fall through to the
# next candidate location) if `dir_path` doesn't exist or doesn't contain
# a complete triplet; errors if a suffix matches more than one file, since
# silently picking one would risk reading the wrong sample's data.
#' @noRd
.read_10x_triplet <- function(dir_path) {
  if (!dir.exists(dir_path)) return(NULL)
  files <- list.files(dir_path)
  if (length(files) == 0) return(NULL)

  find_one <- function(pattern, label) {
    hits <- grep(pattern, files, value = TRUE, ignore.case = TRUE)
    if (length(hits) > 1) {
      stop("Found more than one candidate ", label, " file in '", dir_path,
           "': ", paste(hits, collapse = ", "),
           ". Expected exactly one -- please split these into separate ",
           "sample directories.")
    }
    if (length(hits) == 1) hits else NA_character_
  }

  barcodes_file <- find_one('barcodes\\.tsv(\\.gz)?$', 'barcodes')
  features_file <- find_one('(features|genes)\\.tsv(\\.gz)?$', 'features/genes')
  matrix_file   <- find_one('matrix\\.mtx(\\.gz)?$', 'matrix')

  if (anyNA(c(barcodes_file, features_file, matrix_file))) return(NULL)

  Seurat::ReadMtx(
    mtx      = file.path(dir_path, matrix_file),
    cells    = file.path(dir_path, barcodes_file),
    features = file.path(dir_path, features_file)
  )
}

#' Create Seurat RNA Objects
#'
#' This function creates multiple Seurat objects. It takes a list of directories
#' as input. In each directory, there should be the contents of the filtered
#' feature matrix folder, a filtered feature matrix folder, or a .h5 file. You
#' can directly give the output folder from cellranger into this function.
#' It will preferentially choose the filtered_feature_matrix folder.
#'
#' Barcodes/features(or genes)/matrix files may use the standard CellRanger
#' names (\code{barcodes.tsv.gz}, \code{features.tsv.gz}, \code{matrix.mtx.gz})
#' or be prefixed with a sample name, e.g.
#' \code{GSM9352918_con11_barcodes.tsv.gz} -- common when downloading 10X
#' Supplementary Files from GEO. Any prefix in front of the standard suffix
#' is matched automatically; \code{genes.tsv.gz} (CellRanger v2 naming) is
#' also recognized alongside \code{features.tsv.gz}.
#'
#' @param data_dirs Path to directories containing matrix.mtx, features.tsv, and
#'  barcodes.tsv or .h5 files.
#' @param cells Features must be expressed in at least this many cells
#' @param features Cells must have at least this many features
#' @param mt_pattern Pattern for calculating percent mtDNA. Default
#'   \code{"^mt-"} (mouse gene symbol convention, e.g. \code{"mt-Nd1"}); pass
#'   \code{"^MT-"} for human data.
#' @param rb_pattern Pattern for calculating percent ribosomal-protein reads
#'   (\code{percent.rb}). Default \code{"^(Rp[sl]|RP[SL])"} matches both
#'   mouse (\code{"Rps"}/\code{"Rpl"}) and human (\code{"RPS"}/\code{"RPL"})
#'   gene symbol conventions.
#' @param hb_pattern Pattern for calculating percent hemoglobin reads
#'   (\code{percent.hb}). Default \code{"^(Hb[^p]|HB[^P])"} matches mouse
#'   (\code{"Hba"}/\code{"Hbb"}) and human (\code{"HBA"}/\code{"HBB"})
#'   hemoglobin genes while excluding the unrelated \code{"Hbp1"}/\code{"HBP1"}
#'   gene, a well-known false positive for naive \code{"^Hb"} patterns.
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
#' @param workers Number of parallel workers to use (via
#'   \code{future.apply}) for reading/creating each sample's Seurat object
#'   and for the per-sample \code{calldoublet} calls -- the two most
#'   expensive steps, and both fully independent across samples. Defaults
#'   to \code{length(data_dirs)} (one worker per sample); errors up front
#'   if that (or an explicit value) exceeds \code{parallel::detectCores()},
#'   naming the number of cores actually available. Pass \code{workers = 1}
#'   to run sequentially instead. \code{workers > 1} spins up that many
#'   background R sessions via \code{future::plan(multisession)}, restored
#'   on exit. Note each worker holds its own copy of that sample's data, so
#'   peak memory scales with \code{workers}.
#' @param on_disk Logical; if \code{TRUE}, move each returned object's RNA
#'   counts layer to an on-disk BPCells matrix via \code{\link{ConvertToBPCells}}
#'   as the very last step, after doublet calling/QC. Requires the
#'   \code{BPCells} package (Suggests, not a hard dependency -- checked up
#'   front so this fails fast rather than after building every object).
#'   Note this doesn't reduce peak memory \emph{during} construction (each
#'   object is still built and doublet-called fully in memory first, since
#'   DoubletFinder isn't written to expect a BPCells-backed matrix) -- it
#'   reduces the memory footprint of the \emph{returned} list, which for a
#'   many-sample GEO-scale run is the piece you'd otherwise be stuck
#'   holding for the rest of the session. Default \code{FALSE}.
#' @param bpcells_dir Directory to write the on-disk matrices to when
#'   \code{on_disk = TRUE}. Default \code{file.path(getwd(), "bpcells")}
#'   (one subdirectory per sample underneath it).
#' @return A list of Seurat objects
#' @export

CreateRNAObjects <- function(data_dirs, cells = 3, features = 200,
                             mt_pattern = '^mt-',
                             rb_pattern = '^(Rp[sl]|RP[SL])',
                             hb_pattern = '^(Hb[^p]|HB[^P])',
                             treatment = NULL,
                             object_names = NULL,
                             run_doublet_finder = TRUE,
                             doublet_normalization = c("LogNormalize", "SCT"),
                             doublet_vars_to_regress = "percent.mt",
                             doublet_cluster_resolution = 0.1,
                             filter_doublets = FALSE,
                             workers = length(data_dirs),
                             on_disk = FALSE,
                             bpcells_dir = NULL) {
  workers <- .resolve_workers(workers, n_samples = length(data_dirs),
                              was_default = missing(workers))
  doublet_normalization <- match.arg(doublet_normalization)

  if (workers > 1) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Package 'future.apply' is required for workers > 1. ",
           "install.packages('future.apply')")
    }
    old_plan <- future::plan(future::multisession, workers = workers)
    on.exit(future::plan(old_plan), add = TRUE)
  }

  if (isTRUE(on_disk) && !requireNamespace("BPCells", quietly = TRUE)) {
    stop("Package 'BPCells' is required for on_disk = TRUE. Install with: ",
         "remotes::install_github('bnprks/BPCells/r')")
  }
  if (is.null(bpcells_dir)) bpcells_dir <- file.path(getwd(), "bpcells")

  message(sprintf('--- Reading data and creating Seurat objects (%d directories)%s ---',
                  length(data_dirs),
                  if (workers > 1) sprintf(', %d parallel workers', workers) else ''))
  # Reading + CreateSeuratObject is fully independent per directory, so this
  # runs in parallel when workers > 1 (see .read_one below); otherwise a
  # plain sequential lapply, unchanged from before.

  .read_one <- function(dir) {
    # Look for a (possibly sample-prefixed) barcodes/features(or genes)/
    # matrix triplet directly in `dir`, then in the two conventional
    # CellRanger subdirectories, before falling back to a .h5 file.
    seurat_data <- .read_10x_triplet(dir)

    if (is.null(seurat_data)) {
      seurat_data <- .read_10x_triplet(paste(dir, 'filtered_feature_bc_matrix', sep = '/'))
    }
    if (is.null(seurat_data)) {
      seurat_data <- .read_10x_triplet(paste(paste(dir, 'outs', sep = '/'),
                                             'filtered_feature_bc_matrix', sep = '/'))
    }

    if (!is.null(seurat_data)) {
      return(Seurat::CreateSeuratObject(counts = seurat_data,
                                        min.cells = cells,
                                        min.features = features,
                                        project = basename(dir)))
    }

    if (sum(stringr::str_detect(list.files(dir), '.h5')) > 0) {
      seurat_data <- Seurat::Read10X_h5(
        paste(dir,
              list.files(dir)[sapply(list.files(dir),
                                     function(x) all(c(grepl("filtered", x),
                                                       grepl(".h5", x))))], sep = '/'))

      return(Seurat::CreateSeuratObject(counts = seurat_data,
                                        min.cells = cells,
                                        min.features = features,
                                        project = dirname(dir)))
    }

    stop("Could not find a barcodes/features/matrix triplet or a .h5 file ",
        "in '", dir, "' (or its filtered_feature_bc_matrix subdirectories).")
  }

  seurat_objects <- if (workers > 1) {
    future.apply::future_lapply(data_dirs, .read_one, future.seed = TRUE)
  } else {
    lapply(data_dirs, .read_one)
  }
  # Name the list elements with the base names of the directories
  if (is.null(object_names) == TRUE) {
    names(seurat_objects) <- basename(data_dirs)
  } else {names(seurat_objects) <- object_names}

  message('--- Calculating percent mitochondrial / ribosomal / hemoglobin reads ---')
  seurat_objects <- lapply(seurat_objects, function(obj) {
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = mt_pattern)
    obj[["percent.rb"]] <- Seurat::PercentageFeatureSet(obj, pattern = rb_pattern)
    obj[["percent.hb"]] <- Seurat::PercentageFeatureSet(obj, pattern = hb_pattern)
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
    seurat_obj <- seurat_objects[[i]]
    seurat_obj[["RNA"]] <- methods::as(seurat_obj[["RNA"]], Class = "Assay5")
    return(seurat_obj)
  }), names(seurat_objects))

  # ---- Doublet detection --------------------------------------------------
  # This is normally the dominant cost of the whole function (DoubletFinder's
  # pK parameter sweep fits many models per sample), and each sample's call
  # is fully independent of every other's, so it parallelizes cleanly.
  if (isTRUE(run_doublet_finder)) {
    message(sprintf('--- Calling doublets with DoubletFinder (%s)%s ---',
                    doublet_normalization,
                    if (workers > 1) sprintf(', %d parallel workers', workers) else ''))
    sample_labels <- names(seurat_objects)
    n_samples     <- length(seurat_objects)

    # NB: mapply/future_mapply over `seurat_objects` (rather than lapply-ing
    # over an index into a captured `seurat_objects` list) so each worker
    # only ever receives the one sample it's processing. Indexing into a
    # shared list from inside the worker function would instead export the
    # *entire* list of objects to *every* worker, multiplying peak memory
    # by `workers` for no reason.
    .call_one <- function(obj, i, lab) {
      # Per-sample progress messages only make sense when running
      # sequentially -- with workers > 1 these run in background sessions
      # and wouldn't surface here in order anyway.
      if (workers == 1) {
        message(sprintf('  [%d/%d] %s', i, n_samples, lab))
      }
      out <- calldoublet(obj,
                         samplenameIndex    = i,
                         normalization      = doublet_normalization,
                         vars.to.regress    = doublet_vars_to_regress,
                         cluster_resolution = doublet_cluster_resolution)
      if (isTRUE(filter_doublets)) {
        n_before <- ncol(out)
        out      <- subset(out, doublet_finder == "Singlet")
        if (workers == 1) {
          message(sprintf('    %s: dropped %d doublets (%d singlets remaining)',
                          lab, n_before - ncol(out), ncol(out)))
        }
      }
      out
    }

    results <- if (workers > 1) {
      future.apply::future_mapply(
        .call_one, seurat_objects, seq_len(n_samples), sample_labels,
        SIMPLIFY = FALSE, future.seed = TRUE
      )
    } else {
      mapply(.call_one, seurat_objects, seq_len(n_samples), sample_labels,
             SIMPLIFY = FALSE)
    }
    seurat_objects <- setNames(results, sample_labels)
  }

  message('--- Generating unfiltered QC plots ---')
  # Row-bind just the metadata rather than Seurat::merge()-ing the objects --
  # merge() would also combine the count matrices, which the QC plots below
  # never touch, so it's wasted work (and memory) for large/many-sample runs.
  meta <- dplyr::bind_rows(lapply(seurat_objects, function(x) x@meta.data))
  orig.ident <- nFeature_RNA <- percent.mt <- Freq <- pct <- label <- NULL  # silence R CMD check NSE notes
  gene.plot <- ggplot2::ggplot(meta, ggplot2::aes(orig.ident, nFeature_RNA)) +
    ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered') +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    Ol_Reliable()
  mt.plot <- ggplot2::ggplot(meta, ggplot2::aes(orig.ident, percent.mt)) +
    ggplot2::geom_boxplot() + ggplot2::labs(title = 'Unfiltered') +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    Ol_Reliable()

  qc_plot <- gene.plot + mt.plot

  # Third panel: doublet composition per sample. Only added when doublet
  # calling actually ran (run_doublet_finder = TRUE), since otherwise
  # there's no doublet_finder column to plot.
  #
  # Raw doublet counts aren't comparable across samples with different
  # cell numbers -- a sample with more cells shows more doublets even at
  # the same underlying doublet rate, which is what actually reflects
  # loading concentration and is what QC review cares about. Instead,
  # plot the Singlet/Doublet proportion per sample (stacked to 100%) and
  # label each bar with the doublet rate and the underlying n, so the
  # rate is visible at a glance without losing the counts.
  if ("doublet_finder" %in% colnames(meta)) {
    doublet_counts <- as.data.frame(table(
      orig.ident      = meta$orig.ident,
      doublet_finder  = meta$doublet_finder
    ))
    doublet_counts <- doublet_counts |>
      dplyr::group_by(orig.ident) |>
      dplyr::mutate(pct = Freq / sum(Freq)) |>
      dplyr::ungroup()

    doublet_labels <- doublet_counts[doublet_counts$doublet_finder == "Doublet", ]
    doublet_labels$label <- sprintf('%.1f%%\n(n=%d)', 100 * doublet_labels$pct, doublet_labels$Freq)

    doublet.plot <- ggplot2::ggplot(doublet_counts,
                                    ggplot2::aes(orig.ident, pct, fill = doublet_finder)) +
      ggplot2::geom_col(position = ggplot2::position_stack()) +
      ggplot2::geom_text(data = doublet_labels,
                         ggplot2::aes(orig.ident, 1, label = label),
                         vjust = -0.3, size = 3, lineheight = 0.9) +
      ggplot2::scale_fill_manual(values = c(Singlet = "grey80", Doublet = "firebrick")) +
      ggplot2::labs(title = 'Doublet rate', y = 'Proportion of cells', fill = NULL) +
      ggplot2::coord_cartesian(clip = 'off') +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      Ol_Reliable()

    qc_plot <- qc_plot + doublet.plot
  }

  print(qc_plot)

  if (isTRUE(on_disk)) {
    seurat_objects <- ConvertToBPCells(seurat_objects, assay = "RNA",
                                       layers = "counts", path = bpcells_dir)
  }

  return(seurat_objects)
}
