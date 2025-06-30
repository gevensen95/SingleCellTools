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
#' @param type Spatial Techology (specify "visium" if not true single, otherwise NULL is sufficient)
#' @return Seurat object containing a new assay
#' @export
#'
BuildNicheAssay.multiple_FOVs.MiniBatchKmeans <- function(
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

  # Process each object to create a niche assay
  for (i in seq_along(list.object)) {
    message("Processing object ", i)
    object <- list.object[[i]]
    fov <- list.fov[[i]]

    # Initialize a binary matrix (cells x groups) based on group.by annotation
    cells <- Cells(object[[fov]])
    group.labels <- unlist(object[[group.by]][cells, ])
    groups <- sort(unique(group.labels))
    cell.type.mtx <- matrix(0, nrow = length(cells), ncol = length(groups))
    rownames(cell.type.mtx) <- cells
    colnames(cell.type.mtx) <- groups

    # Populate the matrix (each row gets a 1 in the column for its group)
    cells.idx <- seq_along(cells)
    group.idx <- match(group.labels, groups)
    cell.type.mtx[cbind(cells.idx, group.idx)] <- 1

    # Retrieve tissue coordinates based on type
    if (!is.null(type) && type == "visium") {
      coords <- Seurat::GetTissueCoordinates(object[[fov]], which = "centroids")
      coords <- as.matrix(coords[, c("imagecol", "imagerow")])
      colnames(coords) <- c("x", "y")
    } else {
      coords <- Seurat::GetTissueCoordinates(object[[fov]], which = "centroids")
      rownames(coords) <- coords[["cell"]]
      coords <- as.matrix(coords[, c("x", "y")])
    }

    # Find neighbors using tissue coordinates
    neighbors <- Seurat::FindNeighbors(object = coords,
                                       k.param = neighbors.k,
                                       compute.SNN = FALSE)

    # Create the niche assay using the neighbors information
    sum.mtx <- as.matrix(neighbors[["nn"]] %*% cell.type.mtx)
    niche.assay <- CreateAssayObject(counts = t(sum.mtx))
    object[[assay]] <- niche.assay
    DefaultAssay(object) <- assay

    # Scale the data in the niche assay
    object <- ScaleData(object)

    # Replace the object in the list with the modified object
    list.object[[i]] <- object
  }

  # ----------------------------------------
  # Aggregate scaled data for MiniBatchKmeans
  # ----------------------------------------
  # Each element: rows = cells, columns = features
  DAT <- lapply(seq_along(list.object), function(i) {
    t(list.object[[i]][[assay]]@scale.data)
  })

  # Compute the union of all features (columns) across objects
  all_features <- sort(unique(unlist(lapply(DAT, colnames))))

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

  # ----------------------------
  # Run MiniBatchKmeans clustering
  # ----------------------------
  res.clusters <- data.frame(row.names = rownames(DAT))
  for (k in niches.k.range) {
    message("k=", k)
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
  for (i in seq_along(list.object)) {
    message("Assigning clusters to object ", i)
    object <- list.object[[i]]
    object[[]] <- res.clusters[rownames(object[[]]), ]
    list.object[[i]] <- object
  }

  return(list.object)
}

