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

  message(sprintf('--- Combining %d FOVs into a %d-column grid (offset = %d) ---',
                  n_fovs, n_cols, offset))

  # ---- Pass 1: compute each FOV's grid placement offset -------------------
  # Row/column placement is inherently sequential (each new row's y-offset
  # depends on the tallest FOV placed in prior rows), but only needs each
  # FOV's raw coordinate range to work out -- not the full coordinate or
  # molecule tables -- so this pass stays cheap (O(1) per FOV beyond the
  # unavoidable min/max scan) regardless of FOV/molecule count.
  x_offsets  = numeric(n_fovs)
  y_offsets  = numeric(n_fovs)
  starting_x = 0
  starting_y = 0
  cum_y_max  = -Inf

  for (image in seq_len(n_fovs)) {
    coords_i = obj@images[[image]]$centroids@coords
    x_max_i  = max(coords_i[, "x"])
    y_max_i  = max(coords_i[, "y"])

    x_offsets[image] = starting_x
    y_offsets[image] = starting_y
    cum_y_max        = max(cum_y_max, y_max_i + starting_y)

    if ((image + 1) %% n_cols == 1) {
      starting_x = 0
      starting_y = cum_y_max + offset
    } else {
      starting_x = x_max_i + starting_x + offset
    }
  }

  # ---- Pass 2: apply offsets, accumulate per-FOV pieces in lists ----------
  # Coordinates are collected in lists and combined with a single rbind()
  # per structure at the end, instead of rbind()-ing into an accumulator
  # once per FOV (which is O(n_fovs^2) total copying for many FOVs/genes).
  centroid_coord_list  = vector("list", n_fovs)
  centroid_cells_list  = vector("list", n_fovs)
  molecule_coord_lists = list()  # molecule name -> list of per-FOV coord matrices
  final_molecules      = NULL    # seeded from the first FOV's raw molecules container

  for (image in seq_len(n_fovs)) {
    message(sprintf('  Stitching FOV %d of %d', image, n_fovs))
    dx = x_offsets[image]
    dy = y_offsets[image]

    # --- Centroids ---
    centroids_i = obj@images[[image]]$centroids
    coords_i    = centroids_i@coords
    coords_i[, "x"] = coords_i[, "x"] + dx
    coords_i[, "y"] = coords_i[, "y"] + dy
    centroid_coord_list[[image]] = coords_i
    centroid_cells_list[[image]] = centroids_i@cells
    if (image == 1) final_centroids = centroids_i  # scaffold: bbox min etc. carried over untouched, as before

    # --- Molecules ---
    molecules_i = obj@images[[image]]$molecules
    if (image == 1) final_molecules = molecules_i  # preserves any container entries outside all_molecules, as before

    if (!is.null(molecules_i)) {
      for (molecule in all_molecules) {
        if (is.null(molecules_i[[molecule]])) next

        mol_coords = molecules_i[[molecule]]@coords
        mol_coords[, "x"] = mol_coords[, "x"] + dx
        mol_coords[, "y"] = mol_coords[, "y"] + dy
        molecule_coord_lists[[molecule]] = c(molecule_coord_lists[[molecule]], list(mol_coords))

        if (is.null(final_molecules) ) final_molecules = list()
        if (is.null(final_molecules[[molecule]])) {
          # First time we've seen this molecule -- keep it as the scaffold
          # (class/other slots); @coords/@bbox get overwritten below with
          # the full combined-across-FOVs result.
          final_molecules[[molecule]] = molecules_i[[molecule]]
        }
      }
    }
  }

  # ---- Combine once per structure ------------------------------------------
  final_centroids@coords = do.call(rbind, centroid_coord_list)
  final_centroids@cells  = unlist(centroid_cells_list, use.names = FALSE)
  final_centroids@bbox["x", "max"] = max(final_centroids@coords[, "x"])
  final_centroids@bbox["y", "max"] = max(final_centroids@coords[, "y"])

  for (molecule in names(molecule_coord_lists)) {
    combined_coords = do.call(rbind, molecule_coord_lists[[molecule]])
    mol_obj = final_molecules[[molecule]]
    mol_obj@coords = combined_coords
    mol_obj@bbox["x", "max"] = max(combined_coords[, "x"])
    mol_obj@bbox["y", "max"] = max(combined_coords[, "y"])
    final_molecules[[molecule]] = mol_obj
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
