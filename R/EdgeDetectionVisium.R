#' Identify Spots at Edges in Visium Data
#'
#' Iteratively flags Visium spots that sit on the boundary of the capture
#' area, the tissue itself, or tears within the tissue. Each iteration uses
#' a k-nearest-neighbor search with a radius cutoff: spots with fewer than
#' \code{neighbors} neighbors within the cutoff are flagged as edges, peeled
#' off, and the next iteration runs on the survivors. After four passes the
#' input table has columns \code{Filter} ... \code{Filter4} recording which
#' iteration first flagged each spot (or \code{"Keep"} if it survived all
#' four).
#'
#' \strong{Coordinate source.} The function tries
#' \code{Seurat::GetTissueCoordinates(seurat.obj)} first, which works in
#' Seurat v4/v5 and avoids depending on the spaceranger output directory
#' being present. If \code{seurat.obj} is \code{NULL} or no coordinates can
#' be extracted from it, it falls back to reading
#' \code{tissue_positions(_list).csv} from \code{coord_path}.
#'
#' \strong{Radius cutoff.} The original implementation hard-coded
#' \code{radius = 2}, which is correct for Visium \code{array_row} /
#' \code{array_col} integer coordinates (hex grid spacing). When using
#' \code{GetTissueCoordinates} the coordinates may be in pixels instead,
#' which breaks that fixed value. The function now auto-computes the
#' radius as \code{radius_factor * median nearest-neighbor distance} so
#' the same logic works whether you're in array or pixel units. Pass
#' \code{radius = <number>} to override.
#'
#' @param seurat.obj A Visium Seurat object. If provided, coordinates are
#'   pulled with \code{GetTissueCoordinates} and barcode order is matched
#'   to \code{colnames(seurat.obj)}.
#' @param coord_path Optional path to a spaceranger \code{spatial/}
#'   directory containing \code{tissue_positions(_list).csv}. Used only
#'   when \code{seurat.obj} is \code{NULL} or coordinate extraction from
#'   it fails.
#' @param image Image / FOV name to pull coordinates from when
#'   \code{seurat.obj} has multiple images. \code{NULL} (default) uses
#'   the first one.
#' @param search Search type passed to \code{\link[RANN]{nn2}}. Default
#'   \code{"radius"}.
#' @param neighbors Number of neighbors required. Spots with fewer than
#'   \code{neighbors} neighbors inside the radius cutoff are flagged.
#'   Default 7 (a Visium spot has 6 immediate neighbors on the hex grid;
#'   requiring 7 = self + 6 demands a complete neighborhood).
#' @param radius Absolute radius cutoff in coordinate units. \code{NULL}
#'   (default) auto-computes from \code{radius_factor * median nearest-
#'   neighbor distance}.
#' @param radius_factor Multiplier on the median nearest-neighbor distance
#'   used to auto-compute \code{radius}. Default 1.5.
#' @return A data frame with columns \code{barcode}, the original
#'   coordinate columns, and \code{Filter}, \code{Filter2}, \code{Filter3},
#'   \code{Filter4} (each \code{"Keep"} or \code{"Filter"}).
#' @importFrom Seurat GetTissueCoordinates
#' @importFrom RANN nn2
#' @importFrom stats median
#' @importFrom stringr str_detect
#' @importFrom dplyr filter
#' @importFrom tidyr drop_na
#' @importFrom magrittr %>%
#' @export
EdgeDetectionVisium <- function(seurat.obj    = NULL,
                                coord_path    = NULL,
                                image         = NULL,
                                search        = "radius",
                                neighbors     = 7,
                                radius        = NULL,
                                radius_factor = 1.5) {

  # ---- 1. Pull coordinates ------------------------------------------------
  coords <- NULL

  if (!is.null(seurat.obj)) {
    message("--- Pulling tissue coordinates from Seurat object ---")
    coords <- tryCatch({
      img_names <- names(seurat.obj@images)
      if (length(img_names) == 0L) {
        stop("Seurat object has no @images entries.")
      }
      use_img <- if (is.null(image)) img_names[1] else image
      if (!use_img %in% img_names) {
        stop("`image` '", use_img, "' not found in seurat.obj@images: ",
             paste(img_names, collapse = ", "))
      }
      raw <- Seurat::GetTissueCoordinates(seurat.obj[[use_img]])

      # Normalize into a uniform data frame with a 'barcode' column and
      # two numeric coordinate columns named 'coord1' / 'coord2'. The
      # exact column names from GetTissueCoordinates vary by Seurat
      # version and assay type:
      #   - Visium (v4):  imagerow, imagecol  (rownames = barcodes)
      #   - Visium (v5):  x, y                (column 'cell' = barcodes)
      #   - Imaging:      x, y, cell
      df <- as.data.frame(raw)
      if ("cell" %in% colnames(df)) {
        bc <- as.character(df$cell)
      } else {
        bc <- rownames(df)
      }

      coord_cols <- intersect(c("x", "y", "imagerow", "imagecol"),
                              colnames(df))
      if (length(coord_cols) < 2L) {
        # Fallback: first two numeric columns
        num_cols <- colnames(df)[vapply(df, is.numeric, logical(1))]
        if (length(num_cols) < 2L) {
          stop("Could not find two numeric coordinate columns in ",
               "GetTissueCoordinates output (got: ",
               paste(colnames(df), collapse = ", "), ").")
        }
        coord_cols <- num_cols[1:2]
      } else {
        coord_cols <- coord_cols[1:2]
      }

      out <- data.frame(barcode = bc,
                        coord1  = as.numeric(df[[coord_cols[1]]]),
                        coord2  = as.numeric(df[[coord_cols[2]]]),
                        stringsAsFactors = FALSE)
      # Match order to colnames(seurat.obj) and drop barcodes not in obj
      keep <- intersect(out$barcode, colnames(seurat.obj))
      out <- out[match(keep, out$barcode), , drop = FALSE]
      rownames(out) <- out$barcode
      out
    }, error = function(e) {
      message("  GetTissueCoordinates failed: ", conditionMessage(e))
      message("  Falling back to coord_path.")
      NULL
    })
  }

  if (is.null(coords)) {
    if (is.null(coord_path)) {
      stop("Could not extract coordinates from `seurat.obj` and no ",
           "`coord_path` was provided. Pass one or both.")
    }
    message("--- Reading tissue coordinates from coord_path ---")
    files <- list.files(coord_path)
    hit <- files[stringr::str_detect(files, "tissue_position")]
    if (length(hit) == 0L) {
      stop("No 'tissue_position*' file found in: ", coord_path)
    }
    raw <- read.delim(file.path(coord_path, hit[1]),
                      header = FALSE, sep = ",")
    # Drop a header row if spaceranger v2 wrote one (first column "barcode")
    if (identical(as.character(raw[1, 1]), "barcode")) {
      colnames(raw) <- as.character(unlist(raw[1, ]))
      raw <- raw[-1, , drop = FALSE]
      for (j in 2:6) raw[[j]] <- as.numeric(raw[[j]])
    } else {
      colnames(raw) <- c("barcode", "in_tissue", "array_row", "array_col",
                         "pxl_row_in_fullres", "pxl_col_in_fullres")
    }
    if (!is.null(seurat.obj)) {
      raw <- raw[match(colnames(seurat.obj), raw$barcode), , drop = FALSE]
    }
    coords <- data.frame(barcode = raw$barcode,
                         coord1  = as.numeric(raw$array_row),
                         coord2  = as.numeric(raw$array_col),
                         stringsAsFactors = FALSE)
    rownames(coords) <- coords$barcode
    # Preserve the rest of the columns for the returned table
    coords <- cbind(coords, raw[, setdiff(colnames(raw), "barcode"),
                                drop = FALSE])
  }

  if (nrow(coords) < neighbors + 1L) {
    stop("Only ", nrow(coords), " spots found; need at least ",
         neighbors + 1L, " for k=", neighbors, " nearest-neighbor search.")
  }

  # ---- 2. Decide radius cutoff -------------------------------------------
  if (is.null(radius)) {
    nn_seed <- RANN::nn2(data = coords[, c("coord1", "coord2")],
                         k = 2)$nn.dists[, 2]
    radius <- radius_factor * stats::median(nn_seed)
    message(sprintf(
      "  Auto-computed radius = %.4f (= %.2f * median NN distance %.4f)",
      radius, radius_factor, radius / radius_factor))
  }

  # ---- 3. Helper: one iteration of edge peeling --------------------------
  .one_pass <- function(active_coords) {
    nn <- as.data.frame(
      RANN::nn2(data       = active_coords[, c("coord1", "coord2")],
                k          = neighbors,
                searchtype = search,
                radius     = radius)[[1]])
    nn[nn == 0] <- NA
    rownames(nn) <- active_coords$barcode
    nn$barcode <- rownames(nn)
    kept    <- nn %>% tidyr::drop_na()
    flagged <- nn %>% dplyr::filter(!barcode %in% kept$barcode)
    flagged$barcode
  }

  # ---- 4. Four iterations -------------------------------------------------
  coords$Filter <- "Keep"

  message("--- Edge detection iteration 1 of 4 ---")
  flagged1 <- .one_pass(coords)
  coords[flagged1, "Filter"] <- "Filter"

  message("--- Edge detection iteration 2 of 4 ---")
  coords$Filter2 <- coords$Filter
  flagged2 <- .one_pass(coords[coords$Filter2 == "Keep", , drop = FALSE])
  coords[flagged2, "Filter2"] <- "Filter"

  message("--- Edge detection iteration 3 of 4 ---")
  coords$Filter3 <- coords$Filter2
  flagged3 <- .one_pass(coords[coords$Filter3 == "Keep", , drop = FALSE])
  coords[flagged3, "Filter3"] <- "Filter"

  message("--- Edge detection iteration 4 of 4 ---")
  coords$Filter4 <- coords$Filter3
  flagged4 <- .one_pass(coords[coords$Filter4 == "Keep", , drop = FALSE])
  coords[flagged4, "Filter4"] <- "Filter"

  coords
}
