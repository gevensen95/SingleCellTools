#' Split each spatial FOV into its tissue sections
#'
#' For slides carrying several serial sections per field of view (e.g. a
#' MERSCOPE/Xenium region drawn around a column of spinal cord sections),
#' \code{SplitSections()} finds the separate pieces of tissue in every FOV,
#' optionally applies a table of manual fixes, and writes the result to the
#' object's metadata so each section can be treated as its own image /
#' replicate.
#'
#' Pieces are found by DBSCAN on cell centroids, with a radius that adapts to
#' each FOV's cell density (\code{eps_mult} x the median distance to the
#' 5th-nearest cell). Pieces smaller than \code{min_cells} are treated as
#' debris. A piece much larger than the median piece (\code{split_factor})
#' is assumed to be touching sections and is split at the gaps between them.
#'
#' \strong{Gap search.} Splits (automatic and manual \code{"a;b"}) cut at the
#' lowest dips in cell density along y, i.e. the gaps between sections. The
#' density smoothing is scaled to each FOV's cell spacing (half the median
#' distance to the 5th-nearest cell), so it behaves the same in microns or
#' pixels. A dip only counts as a gap if its density is at most half the
#' piece's typical density, and every resulting part must hold at least
#' \eqn{1/(4k)} of the piece's cells. If fewer clear gaps are found than
#' needed (sections pressed together with no visible boundary), the cut
#' falls halfway between k-means centres instead.
#'
#' \strong{Row axis.} Sections are assumed to be laid out in a row. The row
#' axis (\code{axis = "auto"}, default) is the FOV's longer axis -- whichever
#' of x or y the cells spread further along -- and can be fixed with
#' \code{axis = "x"} or \code{"y"}. Splits cut along it, and automatic labels
#' (\code{S01}, \code{S02}, ...) follow it in increasing order.
#'
#' \strong{Orientation in plots.} Where the row appears depends on the Seurat
#' version: Seurat 5's \code{ImageDimPlot()} draws FOVs with y horizontal, so
#' a row along y reads left to right; Seurat 6 draws x horizontal with y
#' increasing downward, so the same row reads top to bottom. Stored
#' coordinates, and therefore the labels, are the same either way. The check
#' plots below are drawn the way the installed Seurat's \code{ImageDimPlot()}
#' draws each FOV.
#'
#' Torn sections, sections split into halves and touching sections are
#' common, so the intended workflow is two passes:
#' \enumerate{
#'   \item Run with \code{manual = NULL}; inspect the check plots in
#'     \code{plot_dir}.
#'   \item Write a fixes table and run again with \code{manual = <table>}.
#' }
#'
#' The fixes table has columns \code{fov}, \code{from}, \code{to} (extra
#' columns, e.g. \code{comment}, are ignored). \code{from} is an automatic
#' label. \code{to} can be:
#' \itemize{
#'   \item a section name (\code{"sec3"}) -- rename; give several pieces the
#'     same name to merge them (e.g. two halves of one section);
#'   \item \code{NA} / blank -- drop the piece (\code{section_id} = \code{NA});
#'   \item names separated by \code{;} (\code{"sec2;sec1"}) -- split one piece
#'     holding touching sections at the gaps between them; names are assigned
#'     in increasing order along the row axis, i.e. in the same order as the
#'     automatic labels.
#' }
#' Each row refers to the automatic labels, so rows never interact. Pieces of
#' a listed FOV that no row mentions keep their automatic label. Only FOVs
#' that appear in the table get a final \code{section_id}; others get
#' \code{NA}, so every FOV must be checked before it is used.
#'
#' Coordinates are read with the same helper the other spatial functions use,
#' so FOV-based images (Xenium, CosMx, MERSCOPE) and older \code{VisiumV1}
#' images (\code{imagerow}/\code{imagecol}) both work.
#'
#' Two check plots are saved per FOV in \code{plot_dir}:
#' \code{<fov>.png} (pieces with labels and cell counts, drawn like
#' \code{ImageDimPlot()}) and \code{<fov>_pieces.png} (each piece with
#' dorsal-horn marker cells in black and \code{ventral_gene} cells in red; a
#' half section shows one dense dorsal-horn band, a whole section two). The
#' marker defaults are for spinal cord; change \code{dorsal_genes} /
#' \code{ventral_gene} for other tissues, or set \code{plot_dir = NULL}.
#' When \code{manual} is given, \code{<fov>_final.png} shows the result.
#'
#' @param object A Seurat object with one coordinate-based image per region
#'   (\code{Images(object)}).
#' @param manual Optional fixes table: a data frame, or a path to a CSV,
#'   with columns \code{fov}, \code{from}, \code{to}. \code{NULL} (default)
#'   keeps the automatic labels.
#' @param plot_dir Folder for check plots, or \code{NULL} to skip plotting.
#' @param eps_mult DBSCAN radius as a multiple of the median 5th-nearest-cell
#'   distance. Raise (5-6) if torn sections fall into many pieces; lower (3)
#'   if neighbouring sections are joined.
#' @param min_cells Pieces with fewer cells are debris.
#' @param axis Row axis the sections are laid out along: \code{"auto"}
#'   (default; each FOV's longer axis), \code{"x"} or \code{"y"}. Labels are
#'   ordered, and touching sections split, along this axis.
#' @param split_factor Pieces with more than \code{split_factor} x the median
#'   piece's cell count are split automatically. \code{Inf} disables.
#' @param dorsal_genes,ventral_gene Marker genes for the
#'   \code{<fov>_pieces.png} check plot.
#' @param seed Random seed for the k-means fallback (used only when no clear
#'   gap is found). The caller's random-number state is restored on exit.
#' @return \code{object} with metadata columns \code{section_fov} (image
#'   name), \code{section_auto} (automatic label, \code{NA} for debris) and
#'   \code{section_id} (\code{"<fov>_<section>"}; \code{NA} for debris,
#'   dropped pieces, and -- when \code{manual} is given -- FOVs not in it).
#'   An existing \code{fov} column (e.g. CosMx FOV numbers) is left
#'   untouched. The fixes applied and the row axis used for each FOV are
#'   stored in \code{object@misc$SplitSections}.
#' @examples
#' \dontrun{
#' # Pass 1: automatic labels + check plots
#' obj <- SplitSections(obj, plot_dir = "section_checks")
#'
#' # Pass 2: apply fixes
#' fixes <- data.frame(
#'   fov  = "M303agedfemale1",
#'   from = c("S07", "S08", "S05", "S06", "S03", "S04", "S02", "S01"),
#'   to   = c("sec1", "sec1", "sec2", "sec2", "sec3", "sec3", "sec4", "sec5"))
#' obj <- SplitSections(obj, manual = fixes, plot_dir = "section_checks")
#'
#' # split a piece holding two touching sections, drop a fragment
#' fixes <- rbind(fixes,
#'   data.frame(fov = "M319youngmale2", from = "S04", to = "sec2;sec1"),
#'   data.frame(fov = "M320atoh1female1", from = "S02", to = NA))
#'
#' obj <- SubsetSpatial(obj, subset = !is.na(section_id))
#' ImageDimPlot(obj, fov = "M303agedfemale1", group.by = "section_id")
#' }
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @export
SplitSections <- function(object, manual = NULL, plot_dir = NULL,
                          eps_mult = 4, min_cells = 300, split_factor = 1.7,
                          axis = c("auto", "y", "x"),
                          dorsal_genes = c("Tac1", "Calb2", "Prkcg", "Tac2", "Grp", "Cck"),
                          ventral_gene = "Chat", seed = 1) {
  .assert_seurat(object, "object")
  axis <- match.arg(axis)
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("SplitSections() needs the 'dbscan' package: install.packages('dbscan')")
  }
  if (!length(SeuratObject::Images(object))) stop("`object` has no FOV images.")

  # Save/restore the caller's global RNG state (k-means below needs a seed,
  # but the function shouldn't leak it into the caller's session).
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  # ---- coordinates ------------------------------------------------------------
  # .get_fov_coords() handles both FOV images (x/y/cell) and VisiumV1 images
  # (imagerow/imagecol, barcodes as rownames).
  coords <- do.call(rbind, lapply(SeuratObject::Images(object), function(img) {
    xy <- .get_fov_coords(object, img, cells = colnames(object))
    data.frame(cell = rownames(xy), x = xy[, "x"], y = xy[, "y"], fov = img,
               stringsAsFactors = FALSE, row.names = NULL)
  }))
  if (is.null(coords) || !nrow(coords)) stop("No cell coordinates found in `object`'s images.")

  # ---- row axis per FOV -------------------------------------------------------
  by_fov   <- split(coords, coords$fov)
  row_axis <- vapply(by_fov, function(d) {
    if (axis != "auto") return(axis)
    if (diff(range(d$x)) > diff(range(d$y))) "x" else "y"
  }, character(1))

  # ---- automatic pieces -------------------------------------------------------
  auto <- do.call(rbind, Map(.ss_find_pieces, by_fov, row_axis[names(by_fov)],
                             MoreArgs = list(eps_mult = eps_mult, min_cells = min_cells,
                                             split_factor = split_factor)))
  rownames(auto) <- NULL
  sec <- auto

  # ---- manual fixes -----------------------------------------------------------
  if (!is.null(manual)) {
    if (is.character(manual)) {
      manual <- utils::read.csv(manual, stringsAsFactors = FALSE,
                                colClasses = "character", na.strings = c("", "NA"))
    }
    miss <- setdiff(c("fov", "from", "to"), names(manual))
    if (length(miss)) stop("`manual` is missing column(s): ", paste(miss, collapse = ", "))
    manual <- as.data.frame(lapply(manual[, c("fov", "from", "to")],
                                   function(v) trimws(as.character(v))),
                            stringsAsFactors = FALSE)
    manual$to[manual$to %in% ""] <- NA
    unknown <- setdiff(manual$fov, unique(auto$fov))
    if (length(unknown)) warning("`manual` FOVs not in object: ", paste(unknown, collapse = ", "))
    manual <- manual[manual$fov %in% auto$fov, , drop = FALSE]
    sec <- .ss_apply_manual(auto, manual, row_axis)
    sec$section[!sec$fov %in% manual$fov] <- NA
  }

  # ---- plots ------------------------------------------------------------------
  if (!is.null(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    for (f in sort(unique(auto$fov))) {
      ax <- .ss_plot_axes(object, f)
      # Row drawn across the page -> wide plot; drawn down the page -> tall plot.
      dims <- if (row_axis[[f]] == ax$h) c(16, 5.5) else c(7, 16)
      ggplot2::ggsave(file.path(plot_dir, paste0(f, ".png")),
                      .ss_plot_sections(auto, f, ax, row_axis[[f]]),
                      width = dims[1], height = dims[2], dpi = 200, bg = "white")
      ggplot2::ggsave(file.path(plot_dir, paste0(f, "_pieces.png")),
                      .ss_plot_pieces(object, auto, f, dorsal_genes, ventral_gene, ax),
                      width = 16, height = 8, dpi = 200, bg = "white")
      if (!is.null(manual) && f %in% manual$fov) {
        ggplot2::ggsave(file.path(plot_dir, paste0(f, "_final.png")),
                        .ss_plot_sections(sec, f, ax, row_axis[[f]]),
                        width = dims[1], height = dims[2], dpi = 200, bg = "white")
      }
    }
  }

  # ---- write back -------------------------------------------------------------
  # Written to `section_fov`, not `fov`: imaging exports (e.g. CosMx) often
  # already carry a numeric `fov` column, which must not be overwritten.
  section_fov <- coords$fov[match(colnames(object), coords$cell)]
  object$section_fov  <- section_fov
  object$section_auto <- auto$section[match(colnames(object), auto$cell)]
  s <- sec$section[match(colnames(object), sec$cell)]
  object$section_id   <- ifelse(is.na(s), NA_character_, paste(section_fov, s, sep = "_"))
  object@misc$SplitSections <- list(manual = manual, eps_mult = eps_mult,
                                    min_cells = min_cells, split_factor = split_factor,
                                    axis = axis, row_axis = row_axis, seed = seed)

  n_auto <- tapply(auto$section, auto$fov, function(v) length(unique(stats::na.omit(v))))
  message(sprintf("SplitSections: %d FOVs, %s automatic pieces%s",
                  length(n_auto), sum(n_auto),
                  if (is.null(manual)) " (no manual fixes applied)"
                  else sprintf("; manual fixes for %d FOVs", length(unique(manual$fov)))))
  object
}

# ---- internal helpers ---------------------------------------------------------

# Typical cell spacing of a FOV: median distance to the 5th-nearest cell.
#' @keywords internal
#' @noRd
.ss_spacing <- function(xy) {
  xy <- as.matrix(xy)
  k  <- min(6, nrow(xy))
  if (k < 2) return(NA_real_)
  stats::median(RANN::nn2(xy, k = k)$nn.dists[, k])
}

# Split a row of k touching sections at the gaps between them: the k - 1
# lowest dips in cell density along y. k-means midpoints cut into the larger
# section when neighbours differ in size; searching only between k-means
# centres fails when both centres land on the same side of a gap. Instead:
#   - candidate cuts are local minima of the density inside the data range
#     whose density is <= max_rel x the typical (median) density;
#   - candidates are taken lowest-density first, keeping a cut only if every
#     resulting part has >= n / (4k) cells (stops two cuts in one gap);
#   - with fewer than k - 1 clear gaps, fall back to k-means midpoints.
# Returns 1..k, 1 = smallest coordinate along the row axis.
#' @keywords internal
#' @noRd
.ss_split_valley <- function(y, k, bw, max_rel = 0.5) {
  if (k < 2) return(rep(1L, length(y)))
  if (!is.finite(bw) || bw <= 0) bw <- "nrd0"
  d  <- stats::density(y, bw = bw, n = 2048)
  dy <- d$y
  n  <- length(dy)
  inside <- d$x > min(y) & d$x < max(y)
  is_min <- c(FALSE, dy[2:(n - 1)] < dy[1:(n - 2)] & dy[2:(n - 1)] <= dy[3:n], FALSE)
  typical <- stats::median(dy[inside])
  cand <- which(is_min & inside & dy <= max_rel * typical)
  cand <- cand[order(dy[cand])]

  min_n <- length(y) / (4 * k)
  cuts  <- numeric(0)
  for (i in cand) {
    trial <- sort(c(cuts, d$x[i]))
    sizes <- tabulate(findInterval(y, trial) + 1L, nbins = length(trial) + 1L)
    if (all(sizes >= min_n)) cuts <- trial
    if (length(cuts) == k - 1) break
  }
  if (length(cuts) < k - 1) {
    cen  <- sort(stats::kmeans(y, centers = k, nstart = 20)$centers[, 1])
    cuts <- (cen[-1] + cen[-k]) / 2
  }
  findInterval(y, sort(cuts)) + 1L
}

#' @keywords internal
#' @noRd
.ss_find_pieces <- function(d, row_axis, eps_mult, min_cells, split_factor) {
  if (nrow(d) < max(6, min_cells)) {        # too few cells for any piece to count
    d$section <- NA_character_
    return(d)
  }
  xy      <- as.matrix(d[, c("x", "y")])
  spacing <- .ss_spacing(xy)
  eps     <- eps_mult * spacing
  d$piece <- dbscan::dbscan(xy, eps = eps, minPts = 5)$cluster   # 0 = noise
  sz   <- table(d$piece[d$piece > 0])
  keep <- as.integer(names(sz)[sz >= min_cells])
  d$piece[!d$piece %in% keep] <- NA
  if (length(keep)) {
    typical <- stats::median(as.numeric(sz[as.character(keep)]))
    for (p in keep) {
      idx <- which(d$piece == p)
      n_split <- round(length(idx) / typical)
      if (length(idx) > split_factor * typical && n_split >= 2) {
        d$piece[idx] <- p * 100 + .ss_split_valley(d[[row_axis]][idx], n_split, bw = spacing / 2)
      }
    }
  }
  ok  <- !is.na(d$piece)
  cy  <- tapply(d[[row_axis]][ok], d$piece[ok], stats::median)
  ord <- names(sort(cy))
  d$section <- ifelse(ok, sprintf("S%02d", match(as.character(d$piece), ord)), NA_character_)
  d$piece <- NULL
  d
}

#' @keywords internal
#' @noRd
.ss_apply_manual <- function(sec, manual, row_axis) {
  orig <- sec$section                      # rows always refer to automatic labels
  spacing <- list()                         # per-FOV cell spacing, computed on demand
  for (i in seq_len(nrow(manual))) {
    hit <- which(sec$fov == manual$fov[i] & orig %in% manual$from[i])
    if (!length(hit)) {
      warning("manual row ", i, " matched nothing: ", manual$fov[i], " ", manual$from[i])
      next
    }
    to <- manual$to[i]
    if (!is.na(to) && grepl(";", to)) {
      labs <- trimws(strsplit(to, ";")[[1]])
      f    <- manual$fov[i]
      if (is.null(spacing[[f]])) spacing[[f]] <- .ss_spacing(sec[sec$fov == f, c("x", "y")])
      ra   <- row_axis[[f]]
      sec$section[hit] <- labs[.ss_split_valley(sec[[ra]][hit], length(labs), bw = spacing[[f]] / 2)]
    } else {
      sec$section[hit] <- to
    }
  }
  sec
}

# How the installed Seurat's ImageDimPlot() draws a FOV. Seurat 5 always
# put y horizontal and x vertical (flip_xy = TRUE). Seurat 6 puts x
# horizontal with y increasing downward, unless the FOV is flagged
# coords_x_orientation = "vertical" (deprecated LoadXenium(flip.xy = TRUE)).
# Returns the coordinate drawn horizontally (h), vertically (v), and whether
# the vertical axis is reversed.
#' @keywords internal
#' @noRd
.ss_plot_axes <- function(object, fov,
                          seurat_version = utils::packageVersion("Seurat")) {
  img <- object[[fov]]
  flipped <- if (seurat_version >= "6.0.0") {
    methods::.hasSlot(img, "coords_x_orientation") &&
      identical(methods::slot(img, "coords_x_orientation"), "vertical")
  } else {
    TRUE
  }
  if (flipped) list(h = "y", v = "x", rev_v = FALSE) else list(h = "x", v = "y", rev_v = TRUE)
}

# Pieces of one FOV with labels and cell counts, drawn like ImageDimPlot().
# Labels go above the sections when the row runs across the plot, and to
# the right when it runs down the plot.
#' @keywords internal
#' @noRd
.ss_plot_sections <- function(sec, fov, ax, row_axis) {
  d <- sec[sec$fov == fov, ]
  d$ph <- d[[ax$h]]; d$pv <- d[[ax$v]]
  lab <- d[!is.na(d$section), ] %>%
    dplyr::group_by(.data$section) %>%
    dplyr::summarise(h_med = stats::median(.data$ph), v_med = stats::median(.data$pv),
                     h_max = max(.data$ph), v_min = min(.data$pv), v_max = max(.data$pv),
                     n = dplyr::n(), .groups = "drop")
  stagger <- function(k) ifelse(seq_len(k) %% 2 == 1, 0.15, 0.45)
  if (row_axis == ax$h) {                     # row across the plot: labels above
    lab  <- dplyr::arrange(lab, .data$h_med)
    span <- diff(range(d$pv))
    if (ax$rev_v) {                           # reversed axis: top = smallest value
      lab$ly <- min(d$pv) - span * stagger(nrow(lab)); lab$sy <- lab$v_min
    } else {
      lab$ly <- max(d$pv) + span * stagger(nrow(lab)); lab$sy <- lab$v_max
    }
    lab$lx <- lab$h_med; lab$sx <- lab$h_med
  } else {                                    # row down the plot: labels to the right
    lab  <- dplyr::arrange(lab, .data$v_med)
    span <- diff(range(d$ph))
    lab$lx <- max(d$ph) + span * stagger(nrow(lab)); lab$sx <- lab$h_max
    lab$ly <- lab$v_med; lab$sy <- lab$v_med
  }
  lab$txt <- paste0(lab$section, "\n", lab$n, " cells")
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$ph, .data$pv)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$section), size = 0.05, show.legend = FALSE) +
    ggplot2::geom_segment(data = lab, ggplot2::aes(x = .data$lx, xend = .data$sx,
                                                   y = .data$ly, yend = .data$sy),
                          colour = "grey50", linewidth = 0.3) +
    ggplot2::geom_label(data = lab, ggplot2::aes(.data$lx, .data$ly, label = .data$txt,
                                                 fill = .data$section),
                        size = 2.8, show.legend = FALSE) +
    ggplot2::scale_colour_discrete(na.value = "grey80") +
    ggplot2::coord_equal(clip = "off") + ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.margin = ggplot2::margin(20, 60, 10, 10)) +
    ggplot2::labs(title = fov, subtitle = sprintf(
      "Drawn like ImageDimPlot() (Seurat %s). Labels follow %s. Grey = debris / dropped.",
      utils::packageVersion("Seurat"), row_axis))
  if (ax$rev_v) p <- p + ggplot2::scale_y_reverse()
  p
}

# Each piece of one FOV on its own panel, with marker cells highlighted.
#' @keywords internal
#' @noRd
.ss_plot_pieces <- function(object, sec, fov, dorsal_genes, ventral_gene, ax) {
  d <- sec[sec$fov == fov & !is.na(sec$section), ]
  genes <- intersect(c(dorsal_genes, ventral_gene), rownames(object))
  d$dorsal <- FALSE; d$ventral <- FALSE
  if (length(genes)) {
    ex <- Seurat::FetchData(object, vars = genes, cells = d$cell, layer = "counts")
    dg <- intersect(dorsal_genes, genes)
    if (length(dg)) d$dorsal <- rowSums(ex[d$cell, dg, drop = FALSE]) >= 2
    if (ventral_gene %in% genes) d$ventral <- ex[d$cell, ventral_gene] >= 1
  }
  d$ph <- d[[ax$h]]; d$pv <- d[[ax$v]]
  d <- d %>% dplyr::group_by(.data$section) %>%
    dplyr::mutate(px = .data$ph - stats::median(.data$ph), py = .data$pv - stats::median(.data$pv),
                  n = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(panel = paste0(.data$section, " (", .data$n, " cells)"))
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$px, .data$py)) +
    ggplot2::geom_point(colour = "grey85", size = 0.05) +
    ggplot2::geom_point(data = d[d$dorsal, ], colour = "black", size = 0.15) +
    ggplot2::geom_point(data = d[d$ventral, ], colour = "red", size = 0.8) +
    ggplot2::facet_wrap(~panel, nrow = 2) + ggplot2::coord_equal() +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::labs(title = fov, subtitle = paste0(
      "black = dorsal-horn markers (one band = half section, two = whole); red = ", ventral_gene))
  if (ax$rev_v) p <- p + ggplot2::scale_y_reverse()
  p
}
