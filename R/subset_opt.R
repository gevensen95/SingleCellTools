#' Clean Up Molecules Slot
#'
#' This function is removes any molecules that are not found in the FOV, essentially reducing the size of the Seurat object.
#'
#' @param obj Seurat object (spatial)
#' @return a Seurat object
#' @export
CleanMolSlot <- function(obj){

  message(sprintf('--- Cleaning molecules slot across %d FOVs ---',
                  length(names(obj@images))))

  # for each image...
  for(img in names(obj@images)){
    message(paste("  Fixing", img, "FOV..."))


    # convert the molecules slot to a df:
    lapply(obj@images[[img]]$molecules, function(x){
      data.frame(x@coords)
    }) -> test
    for(i in 1:length(test)){
      test[[i]]$gene <- names(test)[i]
    }
    test <- data.table::rbindlist(test)

    # get the FOV edges from
    xmin <- obj@images[[img]]$centroids@bbox[1,1]
    xmax <- obj@images[[img]]$centroids@bbox[1,2]
    ymin <- obj@images[[img]]$centroids@bbox[2,1]
    ymax <- obj@images[[img]]$centroids@bbox[2,2]

    orig <- nrow(test)
    test <- subset(test, x >= xmin & x <= xmax & y >= ymin & y <= ymax)
    new <- nrow(test)
    message(paste0("    Removed ", formatC((orig-new)/orig*100, digits=3), "% of molecules from image slot (", orig-new, ")"))


    fov <- obj@images[[img]]
    # CreateMolecules() is exported by SeuratObject, not attached by this
    # package's NAMESPACE -- an unprefixed call only worked when a caller's
    # session happened to already have SeuratObject attached via library().
    fov[["molecules"]] <- SeuratObject::CreateMolecules(coords = test)
    obj@images[[img]] <- fov
  }

  return(obj)
}

#' Subset A Spatial Seurat Object
#'
#' This function is an updated version of subset to properly subset Cosmx
#' Seurat objects. It has the same parameters as subset.
#'
#' @param object Seurat object
#' @param cleanMolecules Clean molecules slot
#' @return A dataframe of with 4 iterations of filtering (Keep vs. Filter)
#' @export

subset_opt <- function(
    object = NULL,
    subset = NULL,
    cells = NULL,
    idents = NULL,
    features = NULL,
    Update.slots = TRUE,
    Update.object = TRUE,
    cleanMolecules = TRUE,
    ...)
{

  if (Update.slots) {
    message('--- Updating object slots ---')
    object <- SeuratObject::UpdateSlots(object)
  }

  message('--- Cloning object ---')
  obj_subset <- object

  # sanity check - use only cell ids (no indices)
  if (all(is.integer(cells))) {
    cells <- SeuratObject::Cells(obj_subset)[cells]
  }

  if (!missing(subset) || !is.null(idents)) {
    message('--- Extracting cells matched to `subset` and/or `idents` ---')
  }

  if (class(obj_subset) == "FOV") {
    message("  Object class is `FOV`")
    cells <- SeuratObject::Cells(obj_subset)
  } else if (!class(obj_subset) == "FOV" && !missing(subset)) {
    subset <- rlang::enquo(arg = subset)
    # cells to keep in the object
    cells <-
      Seurat::WhichCells(object = obj_subset,
                         cells = cells,
                         idents = idents,
                         expression = subset,
                         return.null = TRUE, ...)
  } else if (!class(obj_subset) == "FOV" && !is.null(idents)) {
    cells <-
      Seurat::WhichCells(object = obj_subset,
                         cells = cells,
                         idents = idents,
                         return.null = TRUE, ...)
  } else if (is.null(cells)) {
    cells <- SeuratObject::Cells(obj_subset)
  }

  # added support for object class `FOV`
  message('--- Matching cells in FOVs ---')
  if (class(obj_subset) == "FOV") {
    message("  Matching cells for object class `FOV`")
    cells_check <- any(SeuratObject::Cells(obj_subset) %in% cells)
  } else {
    # check if cells are present in all FOV
    cells_check <-
      unlist(lapply(seq_along(Seurat::Images(obj_subset)),
             function(i) {
               any(SeuratObject::Cells(obj_subset[[Seurat::Images(obj_subset)[i]]][["centroids"]]) %in% cells)
             }))
  }

  if (all(cells_check)) {
    message('--- Subsetting object (cells found in all FOVs) ---')
    obj_subset <- base::subset(obj_subset,
                               cells = cells,
                               idents = idents,
                               features = features,
                               ...)
    # subset FOVs
    message('--- Subsetting FOVs ---')
    fovs <-
      lapply(seq_along(Seurat::Images(obj_subset)), function(i) {
        base::subset(x = obj_subset[[Seurat::Images(obj_subset)[i]]],
                     cells = cells,
                     idents = idents,
                     features = features,
                     ...)
      })
    # replace subsetted FOVs -- build the updated @images list and assign it
    # once instead of `obj_subset[[name]] <- fovs[[i]]` per FOV (each call
    # goes through Seurat's double-bracket assign validity machinery and a
    # fresh object copy; CosMx/Xenium panels can have 100+ FOVs).
    new_images <- obj_subset@images
    for (i in seq_along(fovs)) { new_images[[Seurat::Images(object)[i]]] <- fovs[[i]] }
    obj_subset@images <- new_images

  } else {
    message('--- Subsetting FOVs (cells found in only some FOVs) ---')
    # if cells are present only in one or several FOVs:
    # subset FOVs
    fovs <-
      lapply(seq_along(Seurat::Images(obj_subset)), function(i) {
        if (any(SeuratObject::Cells(obj_subset[[Seurat::Images(obj_subset)[i]]][["centroids"]]) %in% cells)) {
          message("  Cell subsets are found only in FOV: ", Seurat::Images(obj_subset)[i])
          message("  Subsetting Centroids")
          base::subset(x = obj_subset[[Seurat::Images(obj_subset)[i]]],
                       cells = cells,
                       idents = idents,
                       features = features,
                       ...)
        }
      })
    # remove FOVs with no matching cells
    message("  Removing FOVs where cells are NOT found: ",
            paste0(Seurat::Images(object)[which(!cells_check == TRUE)], collapse = ', '))
    # replace subsetted FOVs -- same single-assignment batching as above.
    # A NULL entry in fovs[[i]] removes that FOV from the list, same as it
    # would removing it from @images one double-bracket assignment at a time.
    new_images <- obj_subset@images
    for (i in seq_along(fovs)) { new_images[[Seurat::Images(object)[i]]] <- fovs[[i]] }
    obj_subset@images <- new_images

    # subset final object
    message('--- Subsetting final object ---')
    obj_subset <-
      base::subset(obj_subset,
                   cells = cells,
                   idents = idents,
                   features = features,
                   ...)
  }

  if (Update.object && !class(obj_subset) == "FOV") {
    message('--- Updating Seurat object ---')
    obj_subset <- Seurat::UpdateSeuratObject(obj_subset) }

  if (cleanMolecules == TRUE) {
    obj_subset <- CleanMolSlot(obj_subset)
    message('--- Object is ready ---')
    return(obj_subset)
  } else {
    message('--- Object is ready ---')
    return(obj_subset)
  }
}
