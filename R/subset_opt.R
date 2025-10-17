#' Clean Up Molecules Slot
#'
#' This function is removes any molecules that are not found in the FOV, essentially reducing the size of the Seurat object.
#'
#' @param obj Seurat object (spatial)
#' @return a Seurat object
#' @export
CleanMolSlot <- function(obj){

  # for each image...
  for(img in names(obj@images)){
    message(paste("Fixing", img, "FOV..."))


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
    message(paste0("Removed ", formatC((orig-new)/orig*100, digits=3), "% of molecules from image slot (", orig-new, ")!\n"))


    fov <- obj@images[[img]]
    fov[["molecules"]] <- CreateMolecules(coords = test)
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
    ...)
{

  if (Update.slots) {
    message("Updating object slots..")
    object %<>% UpdateSlots()
  }

  message("Cloing object..")
  obj_subset <- object

  # sanity check - use only cell ids (no indices)
  if (all(is.integer(cells))) {
    cells <- Cells(obj_subset)[cells]
  }

  if (!missing(subset) || !is.null(idents)) {
    message("Extracting cells matched to `subset` and/or `idents`")
  }

  if (class(obj_subset) == "FOV") {
    message("object class is `FOV` ")
    cells <- Cells(obj_subset)
  } else if (!class(obj_subset) == "FOV" && !missing(subset)) {
    subset <- enquo(arg = subset)
    # cells to keep in the object
    cells <-
      WhichCells(object = obj_subset,
                 cells = cells,
                 idents = idents,
                 expression = subset,
                 return.null = TRUE, ...)
  } else if (!class(obj_subset) == "FOV" && !is.null(idents)) {
    cells <-
      WhichCells(object = obj_subset,
                 cells = cells,
                 idents = idents,
                 return.null = TRUE, ...)
  } else if (is.null(cells)) {
    cells <- Cells(obj_subset)
  }

  # added support for object class `FOV`
  if (class(obj_subset) == "FOV") {
    message("Matching cells for object class `FOV`..")
    cells_check <- any(obj_subset %>% Cells %in% cells)
  } else {
    # check if cells are present in all FOV
    message("Matching cells in FOVs..")
    cells_check <-
      lapply(Images(obj_subset) %>% seq,
             function(i) {
               any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells)
             }) %>% unlist
  }

  if (all(cells_check)) {
    message("Cell subsets are found in all FOVs!", "\n",
            "Subsetting object..")
    obj_subset %<>% base::subset(cells = cells,
                                 idents = idents,
                                 features = features,
                                 ...)
    # subset FOVs
    message("Subsetting FOVs..")
    fovs <-
      lapply(Images(obj_subset) %>% seq, function(i) {
        base::subset(x = obj_subset[[Images(obj_subset)[i]]],
                     cells = cells,
                     idents = idents,
                     features = features,
                     ...)
      })
    # replace subsetted FOVs
    for (i in fovs %>% seq) { obj_subset[[Images(object)[i]]] <- fovs[[i]] }

  } else {
    # if cells are present only in one or several FOVs:
    # subset FOVs
    fovs <-
      lapply(Images(obj_subset) %>% seq, function(i) {
        if (any(obj_subset[[Images(obj_subset)[i]]][["centroids"]] %>% Cells %in% cells)) {
          message("Cell subsets are found only in FOV: ", "\n", Images(obj_subset)[i])
          message("Subsetting Centroids..")
          base::subset(x = obj_subset[[Images(obj_subset)[i]]],
                       cells = cells,
                       idents = idents,
                       features = features,
                       ...)
        }
      })
    # remove FOVs with no matching cells
    message("Removing FOVs where cells are NOT found: ", "\n",
            paste0(Images(object)[which(!cells_check == TRUE)], "\n"))
    # replace subsetted FOVs
    for (i in fovs %>% seq) { obj_subset[[Images(object)[i]]] <- fovs[[i]] }

    # subset final object
    message("..subset final object")
    obj_subset %<>%
      base::subset(cells = cells,
                   idents = idents,
                   features = features,
                   ...)
  }

  if (Update.object && !class(obj_subset) == "FOV") {
    message("Updating object..")
    obj_subset %<>% UpdateSeuratObject() }

  message("Cleaning Molecule Slot")
  obj_subset <- CleanMolSlot(obj_subset)
  message("Object is ready!")
  return(obj_subset)

}
