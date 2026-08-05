#' Construct an assay for spatial niche analysis
#'
#' This function will construct a new assay where each feature is a
#' cell label The values represents the sum of a particular cell label
#' neighboring a given cell.
#'
#' @param list.object list of Seurat objects to do clustering on
#' @param list.fov list of fov names to use for grabbing cells to cluster from list.object.  Should be the same length as list.object
#' @param group.by Cell classifications to count in spatial neighborhood
#' @param assay Name for spatial neighborhoods assay
#' @param cluster.name Name of output clusters
#' @param neighbors.k Number of neighbors to consider for each cell
#' @param niches.k.range Number of clusters to return based on the niche assay.  provide a range
#' @param batch_size Number of mini-batches for ClusterR::MiniBatchKmeans
#' @param num_init  Number of times the algorithm will be run with different centroid seeds for ClusterR::MiniBatchKmeans
#' @param type Spatial Technology (specify "visium" if not true single, otherwise NULL is sufficient).
#'   Kept for backward compatibility; coordinate columns from
#'   \code{GetTissueCoordinates()} are now auto-detected regardless of this
#'   value (handles both the older Visium \code{imagecol}/\code{imagerow}
#'   layout and the current \code{x}/\code{y} layout used by newer Visium
#'   objects and all FOV-based platforms).
#' @return Seurat object containing a new assay
#' @export
#'
BuildMultipleNicheAssays <- function(
    list.object,
    list.fov,
    group.by,
    assay = "niche",
    cluster.name = "niches",
    neighbors.k = 20,
    niches.k.range = 2:30,
    batch_size = 20,
    num_init = 20,
    type = NULL
) {
  message('--- Validating FOVs across input objects ---')
  # Remove objects if the fov is not found in the images slot
  remove_indices <- c()  # initialize list of indices to remove
  for (i in seq_along(list.object)) {
    object <- list.object[[i]]
    fov <- list.fov[[i]]
    if (!fov %in% names(object@images)) {
      warning("fov is not found in the ", i, "-th object. Removing the object from list.object and list.fov. i = ", i)
      remove_indices <- c(remove_indices, i)
    }
  }
  if (length(remove_indices) > 0) {
    for (i in rev(remove_indices)) {
      list.object[[i]] <- NULL
      list.fov[[i]] <- NULL
    }
  }

  message(sprintf('--- Building niche assay for %d objects ---', length(list.object)))
  # Process each object to create a niche assay
  for (i in seq_along(list.object)) {
    message(sprintf('  Processing object %d of %d', i, length(list.object)))
    object <- list.object[[i]]
    fov <- list.fov[[i]]

    # Initialize a binary matrix (cells x groups) based on group.by annotation
    cells <- SeuratObject::Cells(object[[fov]])
    group.labels <- unlist(object[[group.by]][cells, ])
    groups <- sort(unique(group.labels))
    cell.type.mtx <- matrix(0, nrow = length(cells), ncol = length(groups))
    rownames(cell.type.mtx) <- cells
    colnames(cell.type.mtx) <- groups

    # Populate the matrix (each row gets a 1 in the column for its group)
    cells.idx <- seq_along(cells)
    group.idx <- match(group.labels, groups)
    cell.type.mtx[cbind(cells.idx, group.idx)] <- 1

    # Retrieve tissue coordinates. GetTissueCoordinates() column names vary by
    # Seurat/image-class version rather than strictly by `type`:
    #   - Visium, older Seurat (VisiumV1 image class): imagecol, imagerow
    #     (pixel coordinates; rownames = barcodes)
    #   - Visium, current Seurat/SeuratData (VisiumV2 image class) and all
    #     FOV-based platforms (Xenium/CosMx): x, y (column 'cell' = barcodes)
    # `type` is kept for backward compatibility but coordinate columns are
    # now auto-detected so this doesn't break when a Visium object returns
    # the newer x/y layout (see also EdgeDetectionVisium(), which normalizes
    # the same way).
    coords_df <- as.data.frame(
      Seurat::GetTissueCoordinates(object[[fov]], which = "centroids")
    )
    if (all(c("imagecol", "imagerow") %in% colnames(coords_df))) {
      coord_cols <- c("imagecol", "imagerow")
    } else if (all(c("x", "y") %in% colnames(coords_df))) {
      coord_cols <- c("x", "y")
    } else {
      stop("Could not find coordinate columns in GetTissueCoordinates() ",
           "output for FOV '", fov, "' (got: ",
           paste(colnames(coords_df), collapse = ", "), ").")
    }
    if ("cell" %in% colnames(coords_df)) {
      rownames(coords_df) <- coords_df[["cell"]]
    }
    coords <- as.matrix(coords_df[, coord_cols])
    colnames(coords) <- c("x", "y")

    # Find neighbors using tissue coordinates
    neighbors <- Seurat::FindNeighbors(object = coords,
                                       k.param = neighbors.k,
                                       compute.SNN = FALSE)

    # Create the niche assay using the neighbors information
    sum.mtx <- as.matrix(neighbors[["nn"]] %*% cell.type.mtx)
    niche.assay <- SeuratObject::CreateAssayObject(counts = t(sum.mtx))
    object[[assay]] <- niche.assay
    SeuratObject::DefaultAssay(object) <- assay

    # Scale the data in the niche assay
    object <- Seurat::ScaleData(object)

    # Replace the object in the list with the modified object
    list.object[[i]] <- object
  }

  # ----------------------------------------
  # Aggregate scaled data for MiniBatchKmeans
  # ----------------------------------------
  message('--- Aggregating scaled niche data across objects ---')
  # Each element: rows = cells, columns = features
  DAT <- lapply(seq_along(list.object), function(i) {
    t(list.object[[i]][[assay]]@scale.data)
  })

  # Compute the union of all features (columns) across objects
  all_features <- sort(unique(unlist(lapply(DAT, colnames))))
  message(sprintf('  Union of niche features across all objects: %d', length(all_features)))

  # For each matrix, add any missing features as columns of zeros and reorder columns
  DAT <- lapply(DAT, function(mat) {
    missing_features <- setdiff(all_features, colnames(mat))
    if (length(missing_features) > 0) {
      add <- matrix(0, nrow = nrow(mat), ncol = length(missing_features))
      colnames(add) <- missing_features
      mat <- cbind(mat, add)
    }
    mat <- mat[, all_features, drop = FALSE]
    return(mat)
  })

  DAT <- do.call("rbind", DAT)
  message(sprintf('  Combined matrix: %d cells x %d features', nrow(DAT), ncol(DAT)))

  # ----------------------------
  # Run MiniBatchKmeans clustering
  # ----------------------------
  message(sprintf('--- Running MiniBatchKmeans across k = %d:%d ---',
                  min(niches.k.range), max(niches.k.range)))
  res.clusters <- data.frame(row.names = rownames(DAT))
  for (k in niches.k.range) {
    message(sprintf('  k = %d', k))
    newCol <- paste0("kmeans_", k)
    km_mb <- ClusterR::MiniBatchKmeans(
      data = DAT,
      clusters = k,
      batch_size = batch_size,
      num_init = num_init,
      max_iters = 100,
      init_fraction = 0.2,
      initializer = "kmeans++",
      early_stop_iter = 10,
      verbose = FALSE,
      CENTROIDS = NULL,
      tol = 1e-04,
      tol_optimal_init = 0.3,
      seed = 1
    )

    # Predict clusters using the centroids
    res.clusters[, newCol] <- ClusterR::predict_MBatchKMeans(
      data = DAT,
      CENTROIDS = km_mb$centroids
    )
    res.clusters[, newCol] <- as.factor(res.clusters[, newCol])
  }

  colnames(res.clusters) <- paste0(cluster.name, ".", colnames(res.clusters))

  # Assign the clusters back to each object (in the cell metadata)
  message('--- Writing cluster assignments back to each object ---')
  for (i in seq_along(list.object)) {
    message(sprintf('  Assigning clusters to object %d of %d', i, length(list.object)))
    object <- list.object[[i]]
    object[[]] <- res.clusters[rownames(object[[]]), ]
    list.object[[i]] <- object
  }

  return(list.object)
}
