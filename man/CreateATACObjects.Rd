% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreateATACObjects.R
\name{CreateATACObjects}
\alias{CreateATACObjects}
\title{Create Seurat ATAC Objects}
\usage{
CreateATACObjects(
  data_dirs,
  add_treatment = FALSE,
  treatment = NULL,
  peakwidths_max = 10000,
  peakwidths_min = 20,
  passed_filters_value = 500
)
}
\arguments{
\item{data_dirs}{Path to the directories with ATAC data}

\item{add_treatment}{Whether to add a treatment column to metadata}

\item{treatment}{Treatment metadata value (vector)}

\item{peakwidths_max}{Max peak width for finding combined peaks}

\item{peakwidths_min}{Min peak width for finding combined peaks}

\item{passed_filters_value}{Min value for filtering cells based on
passed_filters column}
}
\value{
A list of Seurat objects
}
\description{
This function creates multiple Seurat objects. It takes a list of directories
as input. In each directory, there should be at least the following files:
peaks.bed, singlecell.csv and fragments.tsv.gz. You can directly give the
output folder from cellranger into this function. While reading in the data,
the function creates a common peak list from the samples. Following this,
peaks on scaffolds and not the major chromosomes are removed.
}
