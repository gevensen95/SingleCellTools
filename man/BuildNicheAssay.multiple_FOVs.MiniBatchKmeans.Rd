% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BuildNicheAssay.multiple_FOVs.MiniBatchKmeans.R
\name{BuildNicheAssay.multiple_FOVs.MiniBatchKmeans}
\alias{BuildNicheAssay.multiple_FOVs.MiniBatchKmeans}
\title{Construct an assay for spatial niche analysis}
\usage{
BuildNicheAssay.multiple_FOVs.MiniBatchKmeans(
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
)
}
\arguments{
\item{list.object}{list of Seurat objects to do clustering on}

\item{list.fov}{list of fov names to use for grabbing cells to cluster from list.object.  Should be the same length as list.object}

\item{group.by}{Cell classifications to count in spatial neighborhood}

\item{assay}{Name for spatial neighborhoods assay}

\item{cluster.name}{Name of output clusters}

\item{neighbors.k}{Number of neighbors to consider for each cell}

\item{niches.k.range}{Number of clusters to return based on the niche assay.  provide a range}

\item{batch_size}{Number of mini-batches for ClusterR::MiniBatchKmeans}

\item{num_init}{Number of times the algorithm will be run with different centroid seeds for ClusterR::MiniBatchKmeans}

\item{type}{Spatial Techology (specify "visium" if not true single, otherwise NULL is sufficient)}
}
\value{
Seurat object containing a new assay
}
\description{
This function will construct a new assay where each feature is a
cell label The values represents the sum of a particular cell label
neighboring a given cell.
}
