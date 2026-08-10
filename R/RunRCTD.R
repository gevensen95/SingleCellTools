# Internal: fetch a counts layer, joining split per-sample layers first if
# needed. Merging multiple Seurat objects without SeuratObject::JoinLayers()
# leaves separate layers per sample (e.g. "counts.sample_A", "counts.sample_B"
# instead of one "counts"). SeuratObject::LayerData(obj, layer = "counts")
# matches all of those by substring and silently returns only the first one
# with a warning ("only the first layer is used") -- which for RunRCTD()
# means every spot/cell outside that one sample would get silently wrong
# counts. Detect that case and join first instead of letting it happen
# silently.
.rctd_get_counts <- function(obj, assay, layer = "counts") {
  avail   <- SeuratObject::Layers(obj[[assay]])
  matches <- avail[avail == layer | startsWith(avail, paste0(layer, "."))]
  if (length(matches) > 1) {
    message(sprintf(
      "  Assay '%s' has %d unjoined '%s' layers (%s) -- joining them first ",
      assay, length(matches), layer, paste(matches, collapse = ", ")),
      "so all samples/cells are included (not just the first layer).")
    obj[[assay]] <- SeuratObject::JoinLayers(obj[[assay]], layers = layer)
  }
  SeuratObject::LayerData(obj, assay = assay, layer = layer)
}

#' Visium spot deconvolution with RCTD (spacexr)
#'
#' Wraps \code{spacexr::create.RCTD} + \code{spacexr::run.RCTD} to
#' estimate per-spot cell-type proportions from a reference single-cell
#' Seurat object. RCTD is the recommended primary annotation strategy for
#' Visium because each spot is a mixture of 1-10 cells of different types;
#' a winner-takes-all classifier (like \code{AnnotateClusters}) loses
#' minority populations by construction.
#'
#' Three modes are supported (RCTD's \code{doublet_mode}):
#' \describe{
#'   \item{\code{"full"} (default)}{Assumes each spot is a mixture of many
#'     cell types. Returns full proportion matrix. Best for standard
#'     Visium.}
#'   \item{\code{"doublet"}}{Assumes each spot has one or two cell types.
#'     Faster, cleaner, but less appropriate for high-density tissue.}
#'   \item{\code{"multi"}}{Iteratively fits up to \code{max_cores_multi}
#'     types per spot. Middle ground.}
#' }
#'
#' Adds the resulting weights matrix (per-cell-type proportions) to the
#' Visium object under \code{obj@misc$rctd_weights}, and — when
#' \code{write_metadata = TRUE} — writes each cell type's per-spot
#' proportion as a metadata column named \code{rctd_<celltype>}, plus three
#' derived columns: \code{rctd_dominant} (the highest-weight cell type),
#' \code{rctd_max_weight} (that type's proportion -- our own simple
#' confidence heuristic, low values mean a genuinely ambiguous mixture, not
#' spacexr-computed), and, in \code{"doublet"}/\code{"multi"} mode only,
#' \code{rctd_spot_class} -- spacexr's own per-spot call ("singlet",
#' "doublet_certain", "doublet_uncertain", or "reject"; \code{"full"} mode
#' has no equivalent built-in field, so this column is absent there rather
#' than faked).
#'
#' \code{obj} may be a single Visium Seurat object or a named list of them
#' (samples that haven't been merged into one object yet). The RCTD
#' reference is built once from \code{reference} and reused across every
#' sample -- reference construction is the expensive, sample-independent
#' part, so looping \code{RunRCTD()} by hand over unmerged samples would
#' otherwise redo that work every time.
#'
#' @param obj A Visium Seurat object, or a (optionally named) list of them.
#' @param reference A single-cell Seurat object with a cell-type column
#'   for the reference.
#' @param celltype_col Reference metadata column with cell-type labels.
#' @param assay_query Assay to read counts from on the Visium object.
#'   Default \code{"Spatial"}.
#' @param assay_ref Assay on the reference. Default \code{"RNA"}.
#' @param mode RCTD \code{doublet_mode}. See Details. Default \code{"full"}.
#' @param max_cells_per_ref_celltype Downsample the reference to this
#'   many cells per cell type before building the RCTD reference. Default
#'   10000. Set to \code{Inf} to disable.
#' @param CELL_MIN_INSTANCE Minimum cells per reference cell type required
#'   to keep it (passed through to \code{create.RCTD}). Default 25.
#' @param write_metadata Logical; if TRUE (default) writes proportion
#'   columns and the derived columns described above into obj@meta.data.
#' @param n_cores Number of cores for RCTD. Default 4.
#' @return The Visium object (or list of them, matching the shape of
#'   \code{obj}) with \code{obj@misc$rctd_weights} (spot x cell-type
#'   matrix) and the associated metadata columns (if \code{write_metadata}).
#' @examples
#' \dontrun{
#' visium <- RunRCTD(visium,
#'                   reference    = pbmc_ref,
#'                   celltype_col = "cell_type",
#'                   mode         = "full",
#'                   n_cores      = 8)
#' colnames(visium@misc$rctd_weights)
#' # Downstream: SpatialFeaturePlot on rctd_<celltype> columns
#' SpatialFeaturePlot(visium, features = "rctd_T_cell")
#'
#' # Multiple unmerged samples -- reference built once, reused for each
#' visium_list <- RunRCTD(list(sample1 = visium1, sample2 = visium2),
#'                        reference    = pbmc_ref,
#'                        celltype_col = "cell_type")
#' }
#' @importFrom Seurat DefaultAssay GetTissueCoordinates
#' @importFrom SeuratObject LayerData
#' @importFrom parallel detectCores
#' @importFrom utils assignInNamespace
#' @export
RunRCTD <- function(obj,
                    reference,
                    celltype_col               = NULL,
                    assay_query                = "Spatial",
                    assay_ref                  = "RNA",
                    mode                       = c("full", "doublet", "multi"),
                    max_cells_per_ref_celltype = 10000,
                    CELL_MIN_INSTANCE          = 25,
                    write_metadata             = TRUE,
                    n_cores                    = 4) {

  mode <- match.arg(mode)
  # Cheap argument validation before the (possibly not installed) spacexr
  # dependency check, so bad input fails fast and clearly regardless of
  # whether spacexr happens to be installed.
  if (!inherits(reference, "Seurat")) {
    stop("`reference` must be a Seurat object.")
  }
  if (is.null(celltype_col) ||
      !celltype_col %in% colnames(reference@meta.data)) {
    stop("`celltype_col` must name a column of reference@meta.data.")
  }

  parsed     <- .as_seurat_list(obj)
  objs       <- parsed$objs
  was_single <- parsed$was_single
  orig_names <- names(objs)

  if (!requireNamespace("spacexr", quietly = TRUE)) {
    stop("'spacexr' is required. Install with ",
         "remotes::install_github('dmcable/spacexr').")
  }

  # ---- Build RCTD reference once, shared across every sample -------------
  ref_counts <- .rctd_get_counts(reference, assay = assay_ref, layer = "counts")
  ref_types  <- factor(as.character(reference@meta.data[[celltype_col]]))
  names(ref_types) <- colnames(reference)

  # Optional downsample per cell type
  if (is.finite(max_cells_per_ref_celltype)) {
    keep <- unlist(lapply(levels(ref_types), function(ct) {
      cells <- names(ref_types)[ref_types == ct]
      if (length(cells) > max_cells_per_ref_celltype) {
        cells <- sample(cells, max_cells_per_ref_celltype)
      }
      cells
    }))
    ref_counts <- ref_counts[, keep, drop = FALSE]
    ref_types  <- ref_types[keep]
  }

  message(sprintf("--- Building RCTD reference (%d cells, %d types) ---",
                  length(ref_types), nlevels(ref_types)))
  ref_nUMI        <- Matrix::colSums(ref_counts)
  names(ref_nUMI) <- colnames(ref_counts)
  rctd_ref <- spacexr::Reference(
    counts     = ref_counts,
    cell_types = ref_types,
    nUMI       = ref_nUMI
  )

  # ---- Per-sample puck + RCTD run -----------------------------------------
  # Every element is already guaranteed to be a Seurat object here --
  # .as_seurat_list() validated that above.
  objs <- lapply(seq_along(objs), function(i) {
    o <- objs[[i]]
    tag <- if (length(objs) > 1) paste0(" ('", orig_names[i], "')") else ""
    .run_rctd_one(o, rctd_ref = rctd_ref, assay_query = assay_query, mode = mode,
                 CELL_MIN_INSTANCE = CELL_MIN_INSTANCE,
                 write_metadata = write_metadata, n_cores = n_cores, tag = tag)
  })
  names(objs) <- orig_names

  if (was_single) return(objs[[1]])
  objs
}


# ============================================================================
# Internal: run RCTD for one Visium object against an already-built RCTD
# reference (spacexr::Reference). Building that reference is the expensive,
# sample-independent step, so RunRCTD() builds it once and calls this once
# per sample -- see the "reference built once" note on RunRCTD() itself.
# `tag` is an optional " ('sample_name')" suffix for messages/errors when
# processing more than one sample, empty string for a single object.
# ============================================================================
#' @keywords internal
#' @noRd
.run_rctd_one <- function(obj, rctd_ref, assay_query, mode, CELL_MIN_INSTANCE,
                          write_metadata, n_cores, tag = "") {

  # ---- Build spatial "puck" ------------------------------------------------
  # RCTD wants a Puck object: coords + counts + nUMI. Visium objects merged
  # from multiple samples/capture areas carry one image per sample under
  # obj@images -- pulling coordinates from only the first image would
  # silently drop every spot belonging to every other sample. Gather
  # coordinates from every image and combine them.
  if (length(obj@images) == 0) {
    stop("`obj` has no images in obj@images", tag,
         " -- RunRCTD() needs tissue coordinates.")
  }
  coords_list <- lapply(names(obj@images), function(img_name) {
    ic <- Seurat::GetTissueCoordinates(obj[[img_name]])
    # Normalize the coords columns depending on Seurat version
    if ("cell" %in% colnames(ic)) {
      rownames(ic) <- ic$cell
      ic <- ic[, c("x", "y")]
    } else if (all(c("imagerow", "imagecol") %in% colnames(ic))) {
      ic <- ic[, c("imagerow", "imagecol")]
      colnames(ic) <- c("x", "y")
    } else {
      ic <- ic[, 1:2]
      colnames(ic) <- c("x", "y")
    }
    ic
  })
  coords <- do.call(rbind, coords_list)
  coords <- coords[!duplicated(rownames(coords)), , drop = FALSE]

  cells_in_both <- intersect(rownames(coords), colnames(obj))
  coords <- coords[cells_in_both, , drop = FALSE]

  q_counts <- .rctd_get_counts(obj, assay = assay_query,
                               layer = "counts")[, cells_in_both, drop = FALSE]
  message(sprintf("--- Building query puck%s (%d spots) ---", tag, ncol(q_counts)))
  q_nUMI        <- Matrix::colSums(q_counts)
  names(q_nUMI) <- colnames(q_counts)
  puck <- spacexr::SpatialRNA(
    coords = coords,
    counts = q_counts,
    nUMI   = q_nUMI
  )

  # ---- Run RCTD -----------------------------------------------------------
  # Defensive guard: spacexr's internals call `parallel::detectCores()` in
  # several places (e.g. inside run.RCTD()'s fitBulk/chooseSigma step) with
  # no NA guard, as in `if (parallel::detectCores() > max_cores) ...`. On
  # minimal HPC/container shells missing core-counting tools (`wc`, `nproc`),
  # detectCores() returns NA and that crashes with "missing value where
  # TRUE/FALSE needed" -- deep inside spacexr, not anything under our
  # control. If that's the environment we're in, temporarily patch
  # detectCores() to report `n_cores` for the duration of the RCTD call, and
  # restore it afterward regardless of how the call finishes.
  if (isTRUE(is.na(suppressWarnings(parallel::detectCores())))) {
    message(sprintf(
      "  parallel::detectCores() returned NA in this environment (likely missing 'wc'/'nproc' on a minimal shell) -- reporting %d core(s) to spacexr for the duration of this call.",
      n_cores))
    ns <- asNamespace("parallel")
    orig_detectCores <- get("detectCores", envir = ns)
    utils::assignInNamespace("detectCores", function(...) n_cores, ns = "parallel")
    on.exit(utils::assignInNamespace("detectCores", orig_detectCores, ns = "parallel"),
           add = TRUE)
  }

  message(sprintf("--- Running RCTD%s (mode = %s, cores = %d) ---",
                  tag, mode, n_cores))
  rctd <- spacexr::create.RCTD(
    spatialRNA        = puck,
    reference         = rctd_ref,
    max_cores         = n_cores,
    CELL_MIN_INSTANCE = CELL_MIN_INSTANCE
  )
  rctd <- spacexr::run.RCTD(rctd, doublet_mode = mode)

  # ---- Extract weights and normalize --------------------------------------
  spot_class <- NULL
  if (mode == "full") {
    weights <- as.matrix(rctd@results$weights)
    # Row-normalize to proportions (they can be counts otherwise)
    weights <- sweep(weights, 1, pmax(rowSums(weights), 1e-8), "/")
  } else {
    # For doublet / multi mode, results is a data frame -- but results_df
    # does NOT carry weight_first/weight_second columns (a latent bug here
    # predating this session assumed it did, and crashed on a real user run
    # with "number of items to replace is not a multiple of replacement
    # length" once it hit a row the wrong assumption didn't silently
    # tolerate). Per spacexr's actual documented output, the doublet-mode
    # weights live in the SEPARATE rctd@results$weights_doublet matrix
    # instead: https://github.com/dmcable/spacexr/issues/45
    df        <- rctd@results$results_df
    w_doublet <- rctd@results$weights_doublet
    if (is.null(w_doublet) || is.null(dim(w_doublet)) || ncol(w_doublet) < 2) {
      stop("Expected rctd@results$weights_doublet with 2 columns (first/",
           "second type weights) for mode = '", mode, "', but got ",
           if (is.null(w_doublet)) "NULL" else paste0(ncol(w_doublet), " column(s)"),
           ". This may indicate an incompatible spacexr version.")
    }
    # Column names for weights_doublet aren't consistently documented across
    # spacexr versions -- match by name if there's an obvious "first"/
    # "second" column, otherwise fall back to position (there are only ever
    # 2 columns: first-type weight, second-type weight).
    w1_col <- grep("first",  colnames(w_doublet), ignore.case = TRUE)[1]
    w2_col <- grep("second", colnames(w_doublet), ignore.case = TRUE)[1]
    if (is.na(w1_col)) w1_col <- 1
    if (is.na(w2_col)) w2_col <- 2

    ct <- levels(factor(c(as.character(df$first_type),
                          as.character(df$second_type))))
    weights <- matrix(0, nrow = nrow(df), ncol = length(ct),
                      dimnames = list(rownames(df), ct))
    for (i in seq_len(nrow(df))) {
      # "reject" spots get no confident type call at all (first_type is NA
      # for them) -- leave that row all-zero rather than erroring on an NA
      # column name.
      if (is.na(df$first_type[i])) next
      spot <- rownames(df)[i]
      if (!spot %in% rownames(w_doublet)) next
      w1 <- w_doublet[spot, w1_col]
      w2 <- w_doublet[spot, w2_col]
      weights[i, as.character(df$first_type[i])]  <- w1
      if (!is.na(df$second_type[i])) {
        weights[i, as.character(df$second_type[i])] <- w2
      }
    }
    weights <- sweep(weights, 1, pmax(rowSums(weights), 1e-8), "/")
    # spacexr's own per-spot classification -- a real diagnostic it computes,
    # not something we're deriving: "singlet" (one confident type),
    # "doublet_certain"/"doublet_uncertain" (two types, confidence differs),
    # "reject" (no confident call). Only available in doublet/multi mode --
    # "full" mode has no equivalent built-in field.
    spot_class <- setNames(as.character(df$spot_class), rownames(df))
  }

  obj@misc$rctd_weights <- weights

  # ---- Optional metadata columns -----------------------------------------
  if (isTRUE(write_metadata)) {
    ct_names <- colnames(weights)
    ct_safe  <- make.names(ct_names)

    # Build every new column (per-cell-type weights + derived columns) as
    # one data frame, row-matched to the full meta.data cell set, then write
    # it into meta.data with a single assignment -- instead of one
    # `obj@meta.data[[col]] <- ...` per cell type, each of which copies the
    # whole meta.data data frame (RCTD reference panels can have dozens of
    # cell types).
    all_cells <- rownames(obj@meta.data)
    match_idx <- match(all_cells, rownames(weights))

    weight_block           <- weights[match_idx, , drop = FALSE]
    rownames(weight_block) <- all_cells
    colnames(weight_block) <- paste0("rctd_", ct_safe)
    weight_df              <- as.data.frame(weight_block)

    dom <- ct_names[apply(weights, 1, which.max)]
    weight_df$rctd_dominant <- dom[match_idx]

    # Our own simple confidence heuristic (available in every mode, unlike
    # rctd_spot_class below): the dominant type's proportion at each spot.
    # Low values mean a spot's mixture is genuinely ambiguous, not a
    # confident call.
    max_weight <- apply(weights, 1, max)
    weight_df$rctd_max_weight <- max_weight[match_idx]

    if (!is.null(spot_class)) {
      weight_df$rctd_spot_class <- spot_class[match_idx]
    }

    obj@meta.data[, colnames(weight_df)] <- weight_df

    message(sprintf("  Wrote %d per-cell-type columns + rctd_dominant/rctd_max_weight%s.",
                    length(ct_names), if (!is.null(spot_class)) "/rctd_spot_class" else ""))
  }

  obj
}
