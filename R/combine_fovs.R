#' Combine FOVs
#'
#' This function combines multiple FOVs into one using a predetermined offset
#'
#' @param obj Seurat object (spatial)
#' @param assay GEX Assay
#' @param n_cols Number of columns for combined FOVs
#' @param offset Number of pixels for offset
#' @param fov_name Name of combined FOV
#' @return a Seurat object
#' @export
combine_fovs = function(obj,
                        assay = "RNA",
                        n_cols = 2,
                        offset = 5000,
                        fov_name = "combined",
                        append = TRUE) {
  if (!("Seurat" %in% .packages())) library(Seurat)

  all_molecules = rownames(obj)

  n_fovs     = length(obj@images)
  n_rows     = ceiling(n_fovs / n_cols)
  starting_x = 0
  starting_y = 0
  final_centroids = NULL
  final_molecules = NULL

  message(sprintf('--- Combining %d FOVs into a %d-column grid (offset = %d) ---',
                  n_fovs, n_cols, offset))

  for (image in 1:length(obj@images)) {
    message(sprintf('  Stitching FOV %d of %d', image, n_fovs))

    # --- Centroids ---
    if (image == 1) {
      final_centroids   = obj@images[[image]]$centroids
      updated_centroids = final_centroids
    } else {
      updated_centroids = obj@images[[image]]$centroids
      updated_centroids@coords[, "x"] = updated_centroids@coords[, "x"] + starting_x
      updated_centroids@coords[, "y"] = updated_centroids@coords[, "y"] + starting_y
      final_centroids@coords          = rbind(final_centroids@coords, updated_centroids@coords)
      final_centroids@cells           = c(final_centroids@cells, updated_centroids@cells)
      final_centroids@bbox["x", "max"] = max(final_centroids@coords[, "x"])
      final_centroids@bbox["y", "max"] = max(final_centroids@coords[, "y"])
    }

    # --- Molecules ---
    if (image == 1) {
      final_molecules = obj@images[[image]]$molecules
    } else {
      updated_molecules = obj@images[[image]]$molecules

      for (molecule in all_molecules) {

        # Skip if this image doesn't have this molecule
        if (is.null(updated_molecules[[molecule]])) next

        # Offset coordinates
        updated_molecules[[molecule]]@coords[, "x"] = updated_molecules[[molecule]]@coords[, "x"] + starting_x
        updated_molecules[[molecule]]@coords[, "y"] = updated_molecules[[molecule]]@coords[, "y"] + starting_y
        updated_molecules[[molecule]]@bbox["x", ]   = updated_molecules[[molecule]]@bbox["x", ] + starting_x
        updated_molecules[[molecule]]@bbox["y", ]   = updated_molecules[[molecule]]@bbox["y", ] + starting_y

        if (is.null(final_molecules[[molecule]])) {
          # First time we've seen this molecule — add it fresh
          final_molecules[[molecule]] = updated_molecules[[molecule]]
        } else {
          # Append coords and expand bbox
          final_molecules[[molecule]]@coords          = rbind(final_molecules[[molecule]]@coords,
                                                              updated_molecules[[molecule]]@coords)
          final_molecules[[molecule]]@bbox["x", "max"] = max(final_molecules[[molecule]]@bbox["x", "max"],
                                                             updated_molecules[[molecule]]@bbox["x", "max"])
          final_molecules[[molecule]]@bbox["y", "max"] = max(final_molecules[[molecule]]@bbox["y", "max"],
                                                             updated_molecules[[molecule]]@bbox["y", "max"])
        }
      }
    }

    # --- Update offsets ---
    if ((image + 1) %% n_cols == 1) {
      starting_x = 0
      starting_y = final_centroids@bbox["y", "max"] + offset
    } else {
      starting_x = max(updated_centroids@coords[, "x"]) + offset
    }
  }

  message('--- Building combined FOV ---')
  combined_fov = CreateFOV(
    coords    = final_centroids,
    molecules = final_molecules,
    assay     = assay,
    key       = paste0(assay, "_")
  )

  if (append == FALSE) {
    message('  Removing original FOVs (append = FALSE)')
    obj@images[1:length(obj@images)] = NULL
  }
  obj@images[[fov_name]] = combined_fov
  return(obj)
}
