library(akgfmaps)
library(cowplot)
library(RODBC)

source(here::here("R", "get_connected.R"))

channel <- gapctd::get_connected(schema = "AFSC")

haul_dat <- RODBC::sqlQuery(channel = channel,
                            query = "select e.gear_temperature,
e.surface_temperature,
e.start_latitude,
e.start_longitude,
e.end_latitude,
e.end_longitude,
e.stationid,
e.start_time,
e.stratum,
e.haul_type,
e.performance,
e.cruise, 
e.bottom_depth,
a.survey_definition_id
from racebase.haul e,
race_data.survey_definitions a, 
    race_data.surveys b,
    race_data.cruises c
where b.survey_definition_id = 78
and e.cruise > 198200 
and e.cruise = c.cruise 
and a.survey_definition_id = b.survey_definition_id
and b.survey_id = c.survey_id
and e.performance >= 0 
and c.vessel_id = e.vessel
and (e.haul_type in (3, 13))") |>
  dplyr::mutate(LATITUDE = (START_LATITUDE + END_LATITUDE) / 2,
                LONGITUDE = (START_LONGITUDE + END_LONGITUDE) / 2)


# Load Alaska Region Digital Elevation Model (ARDEM v2.0) - Seth Danielson--------------------------
ARDEM <- raster::raster(here::here("data", "ARDEMv2.0.nc"), varname = "z")

set_extent <- raster::extent(ARDEM)
set_extent@ymin <- 52
set_extent@ymax <- 62
set_extent@xmin <- 180
set_extent@xmax <- 200
ARDEM <- raster::crop(ARDEM, set_extent)

# Reformate coordinates to a standard projection--------------------------------------------------
ARDEM_coords <- raster::coordinates(ARDEM) %>% 
  as.data.frame()
ARDEM_coords$x[ARDEM_coords$x >= 180] <- -360 + ARDEM_coords$x[ARDEM_coords$x >= 180]
ARDEM <- ARDEM %>% 
  as.data.frame() %>% 
  dplyr::bind_cols(ARDEM_coords) %>%
  dplyr::filter(z < 0)

# Back-transform ARDEM to raster------------------------------------------------------------------
coordinates(ARDEM) <- ~ x + y
ARDEM <- raster::rasterFromXYZ(ARDEM, crs = "EPSG:4326")

# Mask to survey area
map_npac <- akgfmaps::get_base_layers(select.region = "ebs.slope", set.crs = "EPSG:4326")
ARDEM <- raster::mask(ARDEM, map_npac$survey.area)

# Rescale by a factor of 10 ------------------------------------------------------------------------
# ARDEM <- raster::aggregate(ARDEM, fact = 10)

# Reproject---------------------------------------------------------------------------------------
map_npac <- akgfmaps::get_base_layers(select.region = "ebs.slope", set.crs = "EPSG:3338")

ARDEM_reproj <- raster::projectRaster(ARDEM, crs = "EPSG:3338")
ARDEM_reproj <- raster::projectRaster(from = ARDEM, to = ARDEM_reproj, crs = "EPSG:3338")
ARDEM_reproj@data@values[which(ARDEM_reproj@data@values >= 0)] <- NA
ARDEM_reproj@data@values <- -1*ARDEM_reproj@data@values

ARDEM_reproj_df <- as.data.frame(raster::rasterToPoints(ARDEM_reproj))

# Load EBS bathymetric raster
Slope_bath_raster <- raster::raster(here::here("data", "EBS", "Bathy.grd"))

Slope_bath_raster_df <- Slope_bath_raster |>
  raster::mask(map_npac$survey.area) |>
  rasterToPoints() |>
  as.data.frame()

png(here::here("plots", "SLOPE_ardem_vs_EBS_bathy.png"), width = 18, height = 10, units = "in", res = 300)
print(
  cowplot::plot_grid(
    ggplot() +
      geom_contour_filled(data = ARDEM_reproj_df,
                          mapping = aes(x = x, y = y, z = layer),
                          breaks = c(seq(200,1200,100), Inf)) +
      scale_fill_viridis_d(name = "Depth (m)") +
      ggtitle("ARDEM v2") +
      theme(legend.position = "bottom"),
    ggplot() +
      geom_contour_filled(data = Slope_bath_raster_df,
                          mapping = aes(x = x, y = y, z = ebs_bath5hac),
                          breaks = c(seq(200,1200,100), Inf)) +
      scale_fill_viridis_d(name = "Depth (m)") +
      ggtitle("EFH_bathy") +
      theme(legend.position = "bottom")
  )
)
dev.off()

# Plot for subareas 1-6
for(ii in 1:6) {
  
  bssa <- akgfmaps::get_base_layers(select.region = paste0("bssa", ii), set.crs = "EPSG:3338")
  
  png(here::here("plots", paste0("SLOPE_bssa_", ii, "_ardem_vs_EBS_bathy.png")), 
      width = 18, height = 10, units = "in", res = 300)
  print(
    cowplot::plot_grid(
      ggplot() +
        geom_contour_filled(data = ARDEM_reproj_df,
                            mapping = aes(x = x, y = y, z = layer),
                            breaks = c(seq(0,1200,100), Inf)) +
        geom_sf(data = bssa$survey.strata, fill = NA, color = "white") +
        coord_sf(xlim = bssa$plot.boundary$x,
                 ylim = bssa$plot.boundary$y) +
        scale_fill_viridis_d(name = "Depth (m)") +
        ggtitle("ARDEM v2") +
        theme(legend.position = "bottom"),
      ggplot() +
        geom_contour_filled(data = Slope_bath_raster_df,
                            mapping = aes(x = x, y = y, z = ebs_bath5hac),
                            breaks = c(seq(0,1200,100), Inf)) +
        geom_sf(data = bssa$survey.strata, fill = NA, color = "white") +
        coord_sf(xlim = bssa$plot.boundary$x,
                 ylim = bssa$plot.boundary$y) +
        scale_fill_viridis_d(name = "Depth (m)") +
        ggtitle("EBS_bathy") +
        theme(legend.position = "bottom")
    )
  )
  dev.off()
  
}


# Examining the correlation between 2022 survey depths and nearest raster cell depths from ARDEM and Slope_bath
haul_loc <- haul_dat |>
  sf::st_as_sf(crs = "EPSG:4326", coords = c("LONGITUDE", "LATITUDE")) |>
  sf::st_transform(crs = "EPSG:3338")

# Convert points to sf
Slope_bath_raster_sf <- sf::st_as_sf(Slope_bath_raster_df, 
                                    coords = c("x", "y"), 
                                    crs = "EPSG:3338")

ARDEM_sf <- sf::st_as_sf(ARDEM_reproj_df, 
                         coords = c("x", "y"), 
                         crs = "EPSG:3338")

haul_loc$EBSBATHY_RASTER_INDEX <- st_nearest_feature(haul_loc, Slope_bath_raster_sf)
haul_loc$EBSBATHY_RASTER_DEPTH <- Slope_bath_raster_sf$ebs_bath5hac[haul_loc$EBSBATHY_RASTER_INDEX ]

haul_loc$ARDEM_RASTER_INDEX <- st_nearest_feature(haul_loc, ARDEM_sf)
haul_loc$ARDEM_RASTER_DEPTH <- ARDEM_sf$layer[haul_loc$ARDEM_RASTER_INDEX]

png(here::here("plots", "SLOPE_ardem_vs_EBS_bathy_points.png"), width = 10, height = 5, units = "in", res = 300)
print(
  ggplot() +
    geom_point(data = haul_loc |>
                 tidyr::pivot_longer(cols = c("EBSBATHY_RASTER_DEPTH", "ARDEM_RASTER_DEPTH")),
               mapping = aes(x = BOTTOM_DEPTH, y = value)) +
    scale_x_continuous(name = "Survey Depth (m)") +
    scale_y_continuous(name = "Raster Depth (m)") +
    facet_grid(~name),
)
dev.off()

png(here::here("plots", "SLOPE_survey_vs_bathy_map.png"), width = 8, height = 8, units = "in", res = 300)
print(
  ggplot() +
    geom_contour_filled(data = Slope_bath_raster_df,
                        mapping = aes(x = x, 
                                      y = y, 
                                      z = ebs_bath5hac),
                        breaks = c(seq(200,1200,100), Inf)) +
    geom_sf(data = haul_loc,
            mapping = aes(color = cut(BOTTOM_DEPTH,
                                      breaks = c(seq(200,1200,100), Inf)))) +
    scale_fill_viridis_d(name = "Depth (m)", drop = FALSE) +
    scale_color_viridis_d(name = "Depth (m)", drop = FALSE, guide = "none") +
    ggtitle("EBS_bathy and EBS slope survey bottom depths") +
    theme(legend.position = "bottom")
)
dev.off()


png(here::here("plots", "SLOPE_survey_vs_ARDEM_map.png"), width = 8, height = 8, units = "in", res = 300)
print(
  ggplot() +
    geom_contour_filled(data = ARDEM_reproj_df,
                        mapping = aes(x = x, y = y, z = layer),
                        breaks = c(seq(200,1200,100), Inf)) +
    geom_sf(data = haul_loc,
            mapping = aes(color = cut(BOTTOM_DEPTH,
                                      breaks = c(seq(200,1200,100), Inf)))) +
    scale_fill_viridis_d(name = "Depth (m)", drop = FALSE) +
    scale_color_viridis_d(name = "Depth (m)", drop = FALSE, guide = "none") +
    ggtitle("ARDEM and 2022 EBS Slope survey bottom depths") +
    theme(legend.position = "bottom")
)
dev.off()




# Effects on slope and aspect ----------------------------------------------------------------------
Slope_bath_slope_df <- 
  raster::terrain(
    Slope_bath_raster |>
      raster::mask(map_npac$survey.area), 
    opt = c('slope', 'aspect'), 
    unit = 'degrees', 
    neighbors = 8 # Queen case (use eight neighboring cells)
  ) |>
  raster::rasterToPoints() |>
  as.data.frame()

ARDEM_slope_df <- raster::terrain(
  ARDEM_reproj |>
    raster::mask(map_npac$survey.area), 
  opt = c('slope', 'aspect'), 
  unit = 'degrees', 
  neighbors = 8 # Queen case (use eight neighboring cells)
) |>
  raster::rasterToPoints() |>
  as.data.frame()

png(here::here("plots", "SLOPE_ardem_vs_EBS_bathy_slope.png"), width = 18, height = 10, units = "in", res = 300)
print(
  cowplot::plot_grid(
    ggplot() +
      geom_tile(data = Slope_bath_slope_df,
                mapping = aes(x = x, 
                              y = y, 
                              fill = cut(slope, c(0, 0.1, 0.2, seq(0.4, 2, 0.4), 10, 40), right = FALSE))) +
      scale_fill_viridis_d(name = expression(Slope~(degree)), drop = FALSE, option = "C") +
      ggtitle("EFH_bathy") +
      theme(legend.position = "bottom"),
    ggplot() +
      geom_tile(data = ARDEM_slope_df,
                mapping = aes(x = x, 
                              y = y, 
                              fill = cut(slope, c(0, 0.1, 0.2, seq(0.4, 2, 0.4), 10, 40), right = FALSE))) +
      scale_fill_viridis_d(name = expression(Slope~(degree)), drop = FALSE, option = "C") +
      ggtitle("ARDEM v2") +
      theme(legend.position = "bottom")
  )
)
dev.off()



png(here::here("plots", "SLOPE_ardem_vs_Slope_bath_aspect.png"), width = 18, height = 10, units = "in", res = 300)
print(
  cowplot::plot_grid(
    ggplot() +
      geom_tile(data = Slope_bath_slope_df,
                mapping = aes(x = x, 
                              y = y, 
                              fill = cut(aspect, breaks = seq(0,360,45)))) +
      scale_fill_viridis_d(name = expression(Aspect~(degree)), drop = FALSE, option = "A") +
      ggtitle("EFH_bathy") +
      theme(legend.position = "bottom"),
    ggplot() +
      geom_tile(data = ARDEM_slope_df,
                mapping = aes(x = x, 
                              y = y, 
                              fill = cut(aspect, breaks = seq(0,360,45)))) +
      scale_fill_viridis_d(name = expression(Aspect~(degree)), drop = FALSE, option = "A") +
      ggtitle("ARDEM v2") +
      theme(legend.position = "bottom")
  )
)
dev.off()