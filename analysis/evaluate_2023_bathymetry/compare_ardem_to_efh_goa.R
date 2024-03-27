library(akgfmaps)
library(cowplot)
library(RODBC)

source(here::here("R", "get_connected.R"))

channel <- gapctd::get_connected(schema = "AFSC")

haul_dat <- RODBC::sqlQuery(channel = channel,
                            query = "select * from racebase.haul where cruise = 202101 and vessel in (148, 176)")

# Load Alaska Region Digital Elevation Model (ARDEM v2.0) - Seth Danielson--------------------------
ARDEM <- raster::raster(here::here("data", "ARDEMv2.0.nc"), varname = "z")

set_extent <- raster::extent(ARDEM)
set_extent@ymin <- 50
set_extent@ymax <- 62
set_extent@xmin <- 180
set_extent@xmax <- 235
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
map_npac <- akgfmaps::get_base_layers(select.region = "goa", set.crs = "EPSG:4326")
ARDEM <- raster::mask(ARDEM, map_npac$survey.area)

# Rescale by a factor of 10 ------------------------------------------------------------------------
# ARDEM <- raster::aggregate(ARDEM, fact = 10)

# Reproject---------------------------------------------------------------------------------------
map_npac <- akgfmaps::get_base_layers(select.region = "goa", set.crs = "EPSG:3338")

ARDEM_reproj <- raster::projectRaster(ARDEM, crs = "EPSG:3338")
ARDEM_reproj <- raster::projectRaster(from = ARDEM, to = ARDEM_reproj, crs = "EPSG:3338")
ARDEM_reproj@data@values[which(ARDEM_reproj@data@values >= 0)] <- NA
ARDEM_reproj@data@values <- -1*ARDEM_reproj@data@values

ARDEM_reproj_df <- as.data.frame(raster::rasterToPoints(ARDEM_reproj))

# Load GOA bathy raster
goa_raster <- EFHSDM::GOA_bathy

goa_raster_df <- goa_raster |>
  raster::mask(map_npac$survey.area) |>
  rasterToPoints() |>
  as.data.frame()

png(here::here("plots", "GOA_ardem_vs_GOA_bathy.png"), width = 10, height = 12, units = "in", res = 300)
print(
  cowplot::plot_grid(
    ggplot() +
      geom_contour_filled(data = ARDEM_reproj_df,
                          mapping = aes(x = x, y = y, z = layer),
                          breaks = c(seq(0,600,50), Inf)) +
      scale_fill_viridis_d(name = "Depth (m)") +
      ggtitle("ARDEM v2") +
      theme(legend.position = "right"),
    ggplot() +
      geom_contour_filled(data = goa_raster_df,
                          mapping = aes(x = x, y = y, z = goa_bathp1c),
                          breaks = c(seq(0,600,50), Inf)) +
      scale_fill_viridis_d(name = "Depth (m)") +
      ggtitle("GOA_bathy") +
      theme(legend.position = "right"),
    nrow = 2
  )
)
dev.off()

# Examining the correlation between 2021 survey depths and nearest raster cell depths from ARDEM and GOA_bathy
haul_loc <- haul_dat |>
  sf::st_as_sf(crs = "EPSG:4326", coords = c("START_LONGITUDE", "START_LATITUDE")) |>
  sf::st_transform(crs = "EPSG:3338")

# Convert points to sf
goa_raster_sf <- sf::st_as_sf(goa_raster_df, 
                                   coords = c("x", "y"), 
                                   crs = "EPSG:3338")

ARDEM_sf <- sf::st_as_sf(ARDEM_reproj_df, 
                         coords = c("x", "y"), 
                         crs = "EPSG:3338")

haul_loc$GOA_BATHY_RASTER_INDEX <- st_nearest_feature(haul_loc, goa_raster_sf)
haul_loc$GOA_BATHY_RASTER_DEPTH <- goa_raster_sf$goa_bathp1c[haul_loc$GOA_BATHY_RASTER_INDEX]

haul_loc$ARDEM_RASTER_INDEX <- st_nearest_feature(haul_loc, ARDEM_sf)
haul_loc$ARDEM_RASTER_DEPTH <- ARDEM_sf$layer[haul_loc$ARDEM_RASTER_INDEX]

png(here::here("plots", "GOA_ardem_vs_goa_bathy_points.png"), width = 10, height = 5, units = "in", res = 300)
print(
  cowplot::plot_grid(
    ggplot() +
      geom_point(data = haul_loc,
                 mapping = aes(x = BOTTOM_DEPTH, GOA_BATHY_RASTER_DEPTH)) +
      scale_x_continuous(name = "GOA_bathy Depth (m)") +
      scale_y_continuous(name = "2021 GOA Survey Depth (m)"),
    ggplot() +
      geom_point(data = haul_loc,
                 mapping = aes(x = BOTTOM_DEPTH, ARDEM_RASTER_DEPTH)) +
      scale_x_continuous(name = "ARDEM Depth (m)") +
      scale_y_continuous(name = "2021 GOA Survey Depth (m)")
  )
)
dev.off()

png(here::here("plots", "GOA_survey_vs_bathy_map.png"), width = 8, height = 8, units = "in", res = 300)
print(
  ggplot() +
    geom_contour_filled(data = goa_raster_df,
                        mapping = aes(x = x, 
                                      y = y, 
                                      z = goa_bathp1c),
                        breaks = c(seq(0,600,50), Inf)) +
    geom_sf(data = haul_loc,
            mapping = aes(color = cut(BOTTOM_DEPTH,
                                      breaks = c(seq(0,600,50), Inf)))) +
    scale_fill_viridis_d(name = "Depth (m)", drop = FALSE) +
    scale_color_viridis_d(name = "Depth (m)", drop = FALSE, guide = "none") +
    ggtitle("GOA_bathy and 2022 GOA survey bottom depths") +
    theme(legend.position = "bottom")
)
dev.off()


png(here::here("plots", "GOA_survey_vs_ARDEM_map.png"), width = 8, height = 8, units = "in", res = 300)
print(
  ggplot() +
    geom_contour_filled(data = ARDEM_reproj_df,
                        mapping = aes(x = x, y = y, z = layer),
                        breaks = c(seq(0,600,50), Inf)) +
    geom_sf(data = haul_loc,
            mapping = aes(color = cut(BOTTOM_DEPTH,
                                      breaks = c(seq(0,600,50), Inf)))) +
    scale_fill_viridis_d(name = "Depth (m)", drop = FALSE) +
    scale_color_viridis_d(name = "Depth (m)", drop = FALSE, guide = "none") +
    ggtitle("ARDEM and 2022 GOA survey bottom depths") +
    theme(legend.position = "bottom")
)
dev.off()
