# SingleCellTools

## Dependencies

  - Seurat
  - Signac
  - EnsDb.Mmusculus.v79
  - dplyr
  - readr
  - ggplot2
  - stringr
  - tidyr
  - RANN
  - rlang
  - purrr

## Purpose of SingleCellTools

This package acts as a wrapper for other single cell functions to reduce the number of functions one has to use. Most of the functions are built around creating, filtering, and integrating single cell data that can be used for downstream analyses. Two of these functions can be run interactively to decide the best possible cutoffs for filtering. Finally, a new function has been added to detect spots at the edges of the Visium capture area, tissue, and tears in tissue. This is vital because these spots tend to have abnormal counts and should be filtered out. 

Be sure to see all options available for each function. They are designed to be flexible for the user. For co-expression of genes, I highly recommend looking at https://github.com/katlande/scCoExpress. 

## Author

K. Garrett Evensen, PhD
Bioinformatics Analyst II, The Salk Institute for Biological Studies
