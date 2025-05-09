% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CreateRNAObjectsFilter.R
\name{CreateRNAObjectsFilter}
\alias{CreateRNAObjectsFilter}
\title{Create and Filter Seurat RNA Objects}
\usage{
CreateRNAObjectsFilter(
  data_dirs,
  cells = 3,
  features = 200,
  treatment = NULL,
  use_quantile = TRUE,
  quantile_value_min = 0.15,
  feature_min = NA,
  feature_max = NA,
  percent_mt_max = NA,
  interactive = FALSE,
  object_names = NULL
)
}
\arguments{
\item{data_dirs}{Path to directories containing matrix.mtx, features.tsv, and
barcodes.tsv or .h5 files.}

\item{cells}{Features must be expressed in at least this many cells}

\item{features}{Cells must have at least this many features}

\item{treatment}{Treatment metadata column (e.g., Age, chemical, etc.)}

\item{use_quantile}{Use quantile filtering method for nFeature_RNA and
percent.mt}

\item{quantile_value_min}{Minimum percentile for filtering (max will be 1-min)}

\item{feature_min}{Minimum cutoff for features}

\item{feature_max}{Maximum cutoff for features}

\item{percent_mt_max}{Maximum cutoff for percent.mt}

\item{interactive}{Run filtering interactively to pick cutoffs. A figure will
be shown of the the features and the percent.mt, or you can open the
automatically saved file}

\item{object_names}{Names for the Seurat objects}
}
\value{
A list of filtered Seurat objects
}
\description{
This function creates multiple Seurat objects and filters them. It takes a
list of directories as input. In each directory, there should be the contents
of the filtered feature matrix folder, a filtered feature matrix folder, or
a .h5 file. You can directly give the output folder from cellranger into this
function. It will preferentially choose the filtered_feature_matrix folder.
The objects are then automatically filtered based on the 15th
and 85th percentile for nFeature_RNA and percent.mt. This can function can
be run interactively to choose cutoffs without having to re-run it.
}
