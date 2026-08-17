#' Parse Polygons
#'
#' This function will parse a data frame where the geometry column is formated like this: POLYGON ((x1 y1, x2 y2)). It will out put a list of data frames.
#'
#' @param df Data frame with polygone data
#' @return List of data frames
#' @export
#'
parse_polygons <- function(df) {
  if (!"geometry" %in% names(df)) {
    stop("`df` must have a 'geometry' column (e.g. 'POLYGON ((x1 y1, x2 y2))').")
  }
  bad <- !grepl("^POLYGON \\(\\(.*\\)\\)$", df$geometry)
  if (any(bad)) {
    stop("`df$geometry` has ", sum(bad), " row(s) not matching the expected ",
        "'POLYGON ((x1 y1, x2 y2))' format (checked row(s): ",
        paste(head(which(bad), 5), collapse = ", "),
        if (sum(bad) > 5) ", ..." else "", "). Parsing them as-is would ",
        "silently produce garbage/NA coordinates instead of an error.")
  }

  message(sprintf('--- Parsing %d polygon geometries ---', nrow(df)))

  # Clean geometry column
  df$geometry <- substr(df$geometry, 11, nchar(df$geometry))
  df$geometry <- substr(df$geometry, 1, nchar(df$geometry) - 2)

  # Helper to parse one geometry string into a data frame
  parse_geometry <- function(geometry_str) {
    coords <- strsplit(geometry_str, ",")[[1]]
    xy_list <- strsplit(trimws(coords), " ")
    geometry_df <- as.data.frame(do.call(rbind, xy_list),
                                 stringsAsFactors = FALSE)
    colnames(geometry_df) <- c("x", "y")
    geometry_df$x <- as.numeric(geometry_df$x)
    geometry_df$y <- as.numeric(geometry_df$y)
    return(geometry_df)
  }

  # Apply to every row of geometry column
  geometry_list <- lapply(df$geometry, parse_geometry)

  # Name list entries by EntityID (if present)
  if ("EntityID" %in% names(df)) {
    names(geometry_list) <- df$EntityID
  }

  return(geometry_list)
}
