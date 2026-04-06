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
                            append = TRUE){
  if(!("Seurat" %in% .packages()))
  {
    library(Seurat)
  }

  #get shape of final image
  n_fovs = length(obj@images)
  n_rows = ceiling(n_fovs / n_cols)

  #make new coordinates for combined object
  starting_x = 0 #where to place leftmost point
  starting_y = 0 #where to place topmost point
  final_centroids = NULL #concatenated version
  final_molecules = NULL

  for(image in 1:length(obj@images))
  {
    #update centroid x and y based on previously defined offsets
    if(image == 1)
    {
      #uses metadata (nsides, theta, radius) from first fov only
      final_centroids = obj@images[[image]]$centroids
      updated_centroids = final_centroids
    } else
    {
      updated_centroids = obj@images[[image]]$centroids
      updated_centroids@coords[,"x"] = updated_centroids@coords[,"x"] + starting_x
      updated_centroids@coords[,"y"] = updated_centroids@coords[,"y"] + starting_y
      final_centroids@coords = rbind(final_centroids@coords, updated_centroids@coords)
      final_centroids@cells = c(final_centroids@cells, updated_centroids@cells)
      #expand out limits
      #note that min is inherited from image 1
      final_centroids@bbox["x", "max"] = max(final_centroids@coords[, "x"])
      final_centroids@bbox["y", "max"] = max(final_centroids@coords[, "y"])
    }

    #update molecule locations
    if(image == 1)
    {
      final_molecules = obj@images[[image]]$molecules
    } else
    {
      updated_molecules = obj@images[[image]]$molecules
      for(molecule in names(updated_molecules)){
        updated_molecules[[molecule]]@coords[,"x"] = updated_molecules[[molecule]]@coords[,"x"] + starting_x
        updated_molecules[[molecule]]@coords[,"y"] = updated_molecules[[molecule]]@coords[,"y"] + starting_y
        updated_molecules[[molecule]]@bbox["x",] = updated_molecules[[molecule]]@bbox["x",] + starting_x
        updated_molecules[[molecule]]@bbox["y",] = updated_molecules[[molecule]]@bbox["y",] + starting_y
      }

      for(molecule in names(final_molecules)){
        final_molecules[[molecule]]@coords = rbind(final_molecules[[molecule]]@coords, updated_molecules[[molecule]]@coords)
        #update boundaries
        #note that min is inherited from image 1
        final_molecules[[molecule]]@bbox["x", "max"] = max(final_molecules[[molecule]]@bbox["x", "max"],
                                                           updated_molecules[[molecule]]@bbox["x", "max"])
        final_molecules[[molecule]]@bbox["y", "max"] = max(final_molecules[[molecule]]@bbox["y", "max"],
                                                           updated_molecules[[molecule]]@bbox["y", "max"])
      }
    }

    #define new offsets based on location in grid for next sample
    #go by columns, then rows
    #if modulus of next section / ncol != 1, then we are in the same row, so increase x
    #if modulus == 1, then we are in a new row, so increase y and reset x
    if((image + 1) %% n_cols == 1)
    {
      starting_x = 0
      #use boundary box so we get global maxima (don't assume similar shapes, sizes)
      starting_y = final_centroids@bbox["y", "max"] + offset
    } else
    {
      #here we want local, not global max
      starting_x = max(updated_centroids@coords[, "x"]) + offset
    }
  }

  combined_fov = CreateFOV(coords = final_centroids,
                           molecules = final_molecules,
                           assay = assay,
                           key = paste0(assay, "_"))

  if(append == FALSE)
  {
    obj@images[1:length(obj@images)] = NULL
  }
  obj@images[[fov_name]] = combined_fov

  return(obj)
}
