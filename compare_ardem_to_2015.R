library(akgfmaps)
library(cowplot)
library(RODBC)

source(here::here("R", "get_connected.R"))

channel <- gapctd::get_connected(schema = "AFSC")

haul_dat <- RODBC::sqlQuery(channel = channel,
                            query = "select * from racebase.haul where cruise in (202201, 202202) and vessel in (94, 162)")


# Load Alaska Region Digital Elevation Model (ARDEM v2.0) - Seth Danielson--------------------------
ARDEM <- raster::raster(here::here("data", "ARDEMv2.0.nc"), varname = "z")

set_extent <- raster::extent(ARDEM)
set_extent@ymin <- 52
set_extent@ymax <- 76
set_extent@xmin <- 165
set_extent@xmax <- 220
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
map_npac <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:4326")
ARDEM <- raster::mask(ARDEM, map_npac$survey.area)

# Rescale by a factor of 10 ------------------------------------------------------------------------
# ARDEM <- raster::aggregate(ARDEM, fact = 10)

# Reproject---------------------------------------------------------------------------------------
map_npac <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = "EPSG:3338")

ARDEM_reproj <- raster::projectRaster(ARDEM, crs = "EPSG:3338")
ARDEM_reproj <- raster::projectRaster(from = ARDEM, to = ARDEM_reproj, crs = "EPSG:3338")
ARDEM_reproj@data@values[which(ARDEM_reproj@data@values >= 0)] <- NA
ARDEM_reproj@data@values <- -1*ARDEM_reproj@data@values 


ARDEM_reproj_df <- as.data.frame(raster::rasterToPoints(ARDEM_reproj))

# Load EBS Raster 2015
ebs_raster_2015 <- raster::raster(here::here("data", "bathy_2015_1km_2020mod.grd"))

ebs_raster_2015_df <- ebs_raster_2015 |>
  raster::mask(map_npac$survey.area) |>
  rasterToPoints() |>
  as.data.frame()

png(here::here("plots", "ardem_vs_bathy2015.png"), width = 18, height = 10, units = "in", res = 300)
print(
  cowplot::plot_grid(
    ggplot() +
      geom_contour_filled(data = ARDEM_reproj_df,
                          mapping = aes(x = x, y = y, z = layer*-1),
                          breaks = c(seq(0,200,10), Inf)) +
      scale_fill_viridis_d(name = "Depth (m)") +
      ggtitle("ARDEM v2") +
      theme(legend.position = "bottom"),
    ggplot() +
      geom_contour_filled(data = ebs_raster_2015_df,
                          mapping = aes(x = x, y = y, z = ebs_bath5hac),
                          breaks = c(seq(0,200,10), Inf)) +
      scale_fill_viridis_d(name = "Depth (m)") +
      ggtitle("bathy_2015_1km_2020mod") +
      theme(legend.position = "bottom")
  )
)
dev.off()




png(here::here("plots", "ardem_vs_bathy2015_fermenter.png"), width = 18, height = 10, units = "in", res = 300)
print(
  cowplot::plot_grid(
    ggplot() +
      geom_tile(data = ARDEM_reproj_df,
                          mapping = aes(x = x, y = y, fill = layer*-1)) +
      scale_fill_fermenter(name = "Depth (m)", breaks = c(seq(0,140, 20), 200), direction = -1) +
      ggtitle("ARDEM v2") +
      theme(legend.position = "bottom"),
    ggplot() +
      geom_tile(data = ebs_raster_2015_df,
                          mapping = aes(x = x, y = y, fill = ebs_bath5hac)) +
      scale_fill_fermenter(name = "Depth (m)", breaks = c(seq(0,140,20), 200), direction = -1) +
      ggtitle("bathy_2015_1km_2020mod") +
      theme(legend.position = "bottom")
  )
)
dev.off()



# Examining the correlation between 2022 survey depths and nearest raster cell depths from ARDEM and bathy_2015_1km_2020mod


haul_loc <- haul_dat |>
  sf::st_as_sf(crs = "EPSG:4326", coords = c("START_LONGITUDE", "START_LATITUDE")) |>
  sf::st_transform(crs = "EPSG:3338")

# Convert points to sf
ebs_raster_2015_sf <- sf::st_as_sf(ebs_raster_2015_df, 
             coords = c("x", "y"), 
             crs = "EPSG:3338")

ARDEM_sf <- sf::st_as_sf(ARDEM_reproj_df, 
                                   coords = c("x", "y"), 
                                   crs = "EPSG:3338")

haul_loc$EBS2015_RASTER_INDEX <- st_nearest_feature(haul_loc, ebs_raster_2015_sf)
haul_loc$EBS2015_RASTER_DEPTH <- ebs_raster_2015_sf$ebs_bath5hac[haul_loc$EBS2015_RASTER_INDEX]

haul_loc$ARDEM_RASTER_INDEX <- st_nearest_feature(haul_loc, ARDEM_sf)
haul_loc$ARDEM_RASTER_DEPTH <- ARDEM_sf$layer[haul_loc$ARDEM_RASTER_INDEX]

png(here::here("plots", "ardem_vs_bathy2015_points.png"), width = 10, height = 5, units = "in", res = 300)
print(
  cowplot::plot_grid(
ggplot() +
  geom_point(data = haul_loc,
             mapping = aes(x = BOTTOM_DEPTH, EBS2015_RASTER_DEPTH)) +
  scale_x_continuous(name = "bathy_2015_1km_2020mod Depth (m)") +
  scale_y_continuous(name = "2022 EBS Shelf Survey Depth (m)"),
ggplot() +
  geom_point(data = haul_loc,
             mapping = aes(x = BOTTOM_DEPTH, -1*ARDEM_RASTER_DEPTH)) +
  scale_x_continuous(name = "ARDEM Depth (m)") +
  scale_y_continuous(name = "2022 EBS Shelf Survey Depth (m)")
)
)
dev.off()

png(here::here("plots", "survey_vs_bathy2015_map.png"), width = 8, height = 8, units = "in", res = 300)
print(
  ggplot() +
    geom_contour_filled(data = ebs_raster_2015_df,
                        mapping = aes(x = x, y = y, z = ebs_bath5hac),
                        breaks = c(seq(0,200,10), Inf)) +
    geom_sf(data = haul_loc,
            mapping = aes(fill = cut(BOTTOM_DEPTH, breaks =  c(seq(0,200,10), Inf))),
            shape = 21) +
    scale_fill_viridis_d(name = "Depth (m)") +
    ggtitle("bathy_2015_1km_2020mod") +
    theme(legend.position = "bottom")
)
dev.off()
  
