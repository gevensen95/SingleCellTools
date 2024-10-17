#' Identify Spots at Edges in Visium Data
#'
#' This function takes coordinate data for all spots in the Visium data and
#' determines which are at the edges of the capture area, tissue boundary,
#' and tears within the tissue.
#'
#' @param coord_path Path to the directory with tissue coordinate file
#' @param seurat.obj Seurat object that corresponds to coord_path to ensure
#' barcodes are the same
#' @param search Search type (see nn2)
#' @param neighbors Number of neighbors to look for
#' @return A dataframe of with 4 iterations of filtering (Keep vs. Filter)
#' @export
EdgeDetectionVisium <- function(coord_path, seurat.obj = NULL,
                                search = 'radius', neighbors = 7) {

  #finds and reads in the data
  coords <- read.delim(paste(coord_path,
                             list.files(coord_path)[which(str_detect(
                               list.files(coord_path),
                               'tissue_position')==TRUE)], sep = '/'),
                       header = F, sep = ',')
  colnames(coords) <- c('barcode', 'in_tissue', 'array_row', 'array_col',
                        'pxl_row_in_fullres', 'pxl_col_in_fullres')

  if(is.null(seurat.obj) == FALSE) {
    #double checks everything is in the right order
    coords <- coords[match(colnames(seurat.obj), coords$barcode),]
  }

  rownames(coords) <- coords$barcode
  closest <- as.data.frame(RANN::nn2(data=coords[,c(3,4)], k=7, searchtype = 'radius',
                               radius = 2)[[1]])
  closest[closest == 0] <- NA

  rownames(closest) <- coords$barcode
  closest$barcode <- rownames(closest)
  closest2 <- closest %>% tidyr::drop_na()
  closest3 <- closest %>% dplyr::filter(!barcode %in% closest2$barcode)
  coords$Filter <- coords$barcode
  coords[rownames(closest2),]$Filter <- 'Keep'
  coords[closest3$barcode,]$Filter <- 'Filter'

  #Now do it again
  coords$Filter2 <- coords$Filter
  coords_red <- coords %>% dplyr::filter(Filter == 'Keep')

  closest <- as.data.frame(RANN::nn2(data=coords_red[,c(3,4)], k=neighbors,
                               searchtype = search,
                               radius = 2)[[1]])
  closest[closest == 0] <- NA

  rownames(closest) <- coords_red$barcode
  closest$barcode <- rownames(closest)
  closest2 <- closest %>% tidyr::drop_na()
  closest3 <- closest %>% dplyr::filter(!barcode %in% closest2$barcode)
  coords[closest3$barcode,]$Filter2 <- 'Filter'

  # Now Again
  coords$Filter3 <- coords$Filter2
  coords_red <- coords %>% dplyr::filter(Filter2 == 'Keep')

  closest <- as.data.frame(RANN::nn2(data=coords_red[,c(3,4)], k=7, searchtype = 'radius',
                               radius = 2)[[1]])
  closest[closest == 0] <- NA

  rownames(closest) <- coords_red$barcode
  closest$barcode <- rownames(closest)
  closest2 <- closest %>% tidyr::drop_na()
  closest3 <- closest %>% dplyr::filter(!barcode %in% closest2$barcode)
  coords[closest3$barcode,]$Filter3 <- 'Filter'

  # Last Time
  coords$Filter4 <- coords$Filter3
  coords_red <- coords %>% dplyr::filter(Filter3 == 'Keep')

  closest <- as.data.frame(RANN::nn2(data=coords_red[,c(3,4)], k=7,
                               searchtype = 'radius',
                               radius = 2)[[1]])
  closest[closest == 0] <- NA

  rownames(closest) <- coords_red$barcode
  closest$barcode <- rownames(closest)
  closest2 <- closest %>% tidyr::drop_na()
  closest3 <- closest %>% dplyr::filter(!barcode %in% closest2$barcode)
  coords[closest3$barcode,]$Filter4 <- 'Filter'

  return(coords)
}
