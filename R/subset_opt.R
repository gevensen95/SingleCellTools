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

  message("Object is ready!")
  return(obj_subset)

}
