#' Create a boundary box polygon
#' 
#' Create a boundary box rectable polygon from an sf object or bbox object. To be included in akgfmaps
#' 
#' @param x An sf or bbox object
#' @param crs coordinate reference system for output. Can be automatically detected if sf is a polygon.
#' @export


bbox_to_polygon <- function(x, crs = NA) {
  
  if(!("bbox" %in% class(x))) {
    if(is.na(crs)) {
      crs <- sf::st_crs(x)
    }
    x <- sf::st_bbox(x)
  } 
  
  stopifnot("bbox_to_polygon: x must be either bbox object or an object that can be converted to a polygon using st_bbox()" = "bbox" %in% class(x))  
  
  bbox_poly <- cbind(c(bbox['xmin'], bbox['xmax'], bbox['xmax'], bbox['xmin'], bbox['xmin']), 
                     c(bbox['ymin'], bbox['ymin'], bbox['ymax'], bbox['ymax'], bbox['ymin'])
  ) |>
    list() |>
    sf::st_polygon() |>
    sf::st_sfc() |>
    sf::st_as_sf()
  
  if(!is.na(crs)) {
    sf::st_crs(bbox_poly) <- "EPSG:3338"
  }
  
  return(bbox_poly)
  
}
