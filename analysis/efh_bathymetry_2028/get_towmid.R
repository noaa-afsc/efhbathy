library(akgfmaps)


# Get tow midpoints
st_line_midpoints <- function(sf_lines = NULL) {
  
  g <- sf::st_geometry(sf_lines)
  
  g_mids <- lapply(g, function(x) {
    
    coords <- as.matrix(x)
    
    # this is just a copypaste of View(maptools:::getMidpoints):
    get_mids <- function (coords) {
      dist <- sqrt((diff(coords[, 1])^2 + (diff(coords[, 2]))^2))
      dist_mid <- sum(dist)/2
      dist_cum <- c(0, cumsum(dist))
      end_index <- which(dist_cum > dist_mid)[1]
      start_index <- end_index - 1
      start <- coords[start_index, ]
      end <- coords[end_index, ]
      dist_remaining <- dist_mid - dist_cum[start_index]
      mid <- start + (end - start) * (dist_remaining/dist[start_index])
      return(mid)
    }
    
    mids <- st_point(get_mids(coords))
  })
  
  geometry <- sf::st_sfc(g_mids, 
                    crs = sf::st_crs(sf_lines))
  
  out <- sf::st_drop_geometry(sf_lines)
  
  out <- cbind(out[, which(names(out) != "geometry")], sf::st_sf(geometry))
  
  return(out)
}

channel <- navmaps::get_connected(schema = "AFSC")

tow_region <- c("BS", "AI", "GOA")
label_region <- c("ebs", "ai", "goa")

for(ii in 1:length(tow_region)) {
  
  haul_data <- RODBC::sqlQuery(
    channel = channel,
    query = paste0(
      "select 
        vessel, 
        cruise, 
        haul, 
        performance, 
        start_latitude, 
        end_latitude, 
        start_longitude, 
        end_longitude, 
        bottom_depth 
      from 
        racebase.haul 
      where 
        region = '", tow_region[ii], "' 
        and haul_type = 3 
        and performance >= 0 
        and abundance_haul = 'Y' 
        and bottom_depth > 0
        and cruise > 200000
        and gear in (30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 172)")) |>
      dplyr::filter(!is.na(START_LONGITUDE), !is.na(START_LATITUDE), !is.na(END_LONGITUDE), !is.na(END_LATITUDE))
    
  towpaths <- dplyr::bind_rows(
    haul_data |>
      dplyr::select(-END_LONGITUDE, -END_LATITUDE) |>
      tidyr::pivot_longer(cols = c(START_LONGITUDE, START_LATITUDE)) |>
      dplyr::mutate(name = if_else(name == "START_LATITUDE", "LATITUDE", "LONGITUDE")) |>
      tidyr::pivot_wider(names_from = name, values_from = value),
    haul_data |>
      dplyr::select(-START_LONGITUDE, -START_LATITUDE) |>
      tidyr::pivot_longer(cols = c(END_LONGITUDE, END_LATITUDE)) |>
      dplyr::mutate(name = if_else(name == "END_LATITUDE", "LATITUDE", "LONGITUDE")) |>
      tidyr::pivot_wider(names_from = name, values_from = value)
  ) |>
    dplyr::filter(!is.na(LATITUDE), !is.na(LONGITUDE)) |>
    sf::st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = "WGS84") |>
    dplyr::arrange(VESSEL, CRUISE, HAUL) |>
    dplyr::group_by(VESSEL, CRUISE, HAUL, BOTTOM_DEPTH) |>
    dplyr::summarise(do_union = FALSE) |>
    sf::st_cast(to = "LINESTRING") |>
    st_line_midpoints() |>
    dplyr::rename(DEPTH = BOTTOM_DEPTH)
  
  sf::st_write(obj = towpaths, 
               dsn = here::here("analysis", 
                                "efh_bathymetry_2028", 
                                "data", 
                                "bts_tows", 
                                paste0(
                                  label_region[ii], 
                                  "_", 
                                  "towmid.shp")
               ),
               append = FALSE
  )
  
}
