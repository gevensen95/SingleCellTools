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
    niches.k.range = 2:30 ,
    batch_size = 20,
    num_init = 20
) {
  # check for fov in sample set
  # remove if not found in object
  remove = NULL # init list of indices to remove
  for ( i in seq_along(list.object) ){ # message(i)
    # get object and fov for each object
    object = list.object[[i]]
    fov = list.fov[[i]]

    if( !fov %in%  names(object@images) ){
      warning( "fov is not found in the i-th object.  Removing the object from the list.object and list.fov.  i =", i)
      remove = c(remove, i)
    }
  }
  for (i in rev(remove) ){
    list.object[[i]] = NULL
    list.fov[[i]] = NULL
  }


  for ( i in seq_along(list.object) ){ message(i)
    # get object and fov for each object
    object = list.object[[i]]
    fov = list.fov[[i]]

    # initialize an empty cells x groups binary matrix
    cells <- Cells( object[[fov]] )
    group.labels <- unlist(object[[group.by]][cells, ] )
    groups <- sort( unique(group.labels) )
    cell.type.mtx <- matrix(
      "data" = 0
      , "nrow" = length(cells)
      , "ncol" = length(groups)
    )
    rownames(cell.type.mtx) <- cells
    colnames(cell.type.mtx) <- groups

    # populate the binary matrix
    cells.idx <- seq_along(cells)
    group.idx <- match(group.labels, groups)
    cell.type.mtx[cbind(cells.idx, group.idx)] <- 1

    # find neighbors based on tissue position
    coords <- Seurat::GetTissueCoordinates( object[[fov]], "which" = "centroids" )
    rownames(coords) <- coords[["cell"]]
    coords <- as.matrix(coords[ , c("x", "y")])
    neighbors <- Seurat::FindNeighbors(
      "object" = coords
      , "k.param" = neighbors.k # Defines k for the k-nearest neighbor algorithm
      , "compute.SNN" = F
    )

    # create niche assay
    sum.mtx <- as.matrix( neighbors[["nn"]] %*% cell.type.mtx )
    niche.assay <- CreateAssayObject( "counts" = t(sum.mtx) )
    object[[assay]] <- niche.assay
    DefaultAssay(object) <- assay

    # scale data
    object <- ScaleData(object)

    # return edited object to list
    list.object[[i]] = object


  }

  # get aggregate data for ClusterR::MiniBatchKmeans
  # columns = features
  # rows = cells
  # cells = values
  DAT = lapply( seq_along(list.object), function(i){
    t( list.object[[i]][[assay]]@scale.data )
  }  )
  DAT <- do.call("rbind", DAT)



  res.clusters = data.frame(row.names = rownames(DAT))

  for ( k in niches.k.range ){ message("k=", k)
    # new column name
    newCol = paste0("kmeans_", k)
    # get centroids
    km_mb = ClusterR::MiniBatchKmeans(
      "data" = DAT
      , "clusters" = k # the number of clusters
      , "batch_size" = batch_size # the size of the mini batches
      , "num_init" = num_init # number of times the algorithm will be run with different centroid seeds
      , "max_iters" = 100 # the maximum number of clustering iterations.
      , "init_fraction" = 0.2 # percentage of data to use for the initialization centroids (applies if initializer is kmeans++ or optimal_init). Should be a float number between 0.0 and 1.0.
      , "initializer" = "kmeans++" # the method of initialization. One of, optimal_init, quantile_init, kmeans++ and random. See details for more information
      , "early_stop_iter" = 10 # continue that many iterations after calculation of the best within-cluster-sum-of-squared-error
      , "verbose" = F
      , "CENTROIDS" = NULL
      , "tol" = 1e-04
      , "tol_optimal_init" = 0.3
      , "seed" = 1
    )

    # use centroids to get clusters

    res.clusters[,newCol] = ClusterR::predict_MBatchKMeans( # This function takes the data and the output centroids and returns the clusters.
      "data" = DAT
      , "CENTROIDS" = km_mb$centroids
    )
    res.clusters[,newCol] = as.factor( res.clusters[,newCol] ) # change clusters to factors

  }

  # get clusters back onto the objects
  colnames(res.clusters) = paste0(cluster.name,".", colnames(res.clusters))
  for ( i in seq_along(list.object) ){ message(i)
    # get object and fov for each object
    object = list.object[[i]]

    # get clusters in correct cell row order into metadata of object
    object[[]] = res.clusters[rownames(object[[]]),]

    # return edited object to list
    list.object[[i]] = object
  }


  return(list.object)
}



