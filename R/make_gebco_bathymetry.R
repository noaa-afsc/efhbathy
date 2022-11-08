library(akgfmaps)
library(cowplot)
library(RODBC)

dir.create(path = here::here("output", "EBS"))

source(here::here("R", "get_connected.R"))

# Load GEBCO netCTD file (subset of full 2022 GEBCO bathymetry- downloaded on October 26, 2022)
gebco_2022 <- raster::raster(here::here("data", "gebco_2022_n75.0_s50.0_w-179.99_e-130.0.nc"), 
                             varname = "elevation", # Elevation (up = positive)
                             crs = "EPSG:4269") # WGS84

# Load EFH bathymetry
ebs_bathy_raster <- raster::raster(here::here("data", "EBS", "Bathy.grd"))

# Raster of -1 and NA where -1 indicates locations with data in the EFH bathymety. Multiplying this by the GEBCO will (1) transform GEBCO elevation (up = positive) to same scale as EBS bathymetry (down = positive) and (2) mask the transformed GEBCO raster to locations with data without re-masking using the survey area shapefile
ebs_bathy_locs <- round(ebs_bathy_raster/1e8) - 1

# Project GEBCO to the extent and resolution of the EBS bathymetry raster in the CRS of the bathymetry raster
gebco_epsg3338 <- raster::projectRaster(from = gebco_2022, to = ebs_bathy_raster, method = "bilinear")

# Multiply by binary raster to only include data at the same locations as the EBS bathymetry raster
gebco_epsg3338 <- gebco_epsg3338*ebs_bathy_locs

names(gebco_epsg3338) <- "seafloor_depth"

# Write bathymetry raster to .grd file
writeRaster(x = gebco_epsg3338, filename = here::here("output", "EBS", "gebco_bathymetry.grd"))

# Compare EBS EFH 2022 bathymetry raster to GEBCO
ebs_bathy_raster_df <- ebs_bathy_raster |>
  rasterToPoints() |>
  as.data.frame()

# Plot GEBCO and EBS shelf bathymetry rasters
gebco_df <- as.data.frame(raster::rasterToPoints(gebco_epsg3338))

png(here::here("plots", "EBS_GEBCO_vs_EBS_bathy.png"), width = 18, height = 10, units = "in", res = 300)
print(
  cowplot::plot_grid(
    ggplot() +
      geom_contour_filled(data = gebco_df,
                          mapping = aes(x = x, y = y, z = seafloor_depth),
                          breaks = c(seq(0,200,10), Inf)) +
      scale_fill_viridis_d(name = "Depth (m)") +
      ggtitle("ARDEM v2") +
      theme(legend.position = "bottom"),
    ggplot() +
      geom_contour_filled(data = ebs_bathy_raster_df,
                          mapping = aes(x = x, y = y, z = ebs_bath5hac),
                          breaks = c(seq(0,200,10), Inf)) +
      scale_fill_viridis_d(name = "Depth (m)") +
      ggtitle("EBS_bathy") +
      theme(legend.position = "bottom")
  )
)
dev.off()


# Compare GEBCO 2022 to observed bottom-trawl survey depths in 2022
channel <- gapctd::get_connected(schema = "AFSC")

haul_dat <- RODBC::sqlQuery(channel = channel,
                            query = "select * from racebase.haul where cruise in (202201, 202202) and vessel in (94, 162)")

haul_loc <- haul_dat |>
  sf::st_as_sf(crs = "EPSG:4326", coords = c("START_LONGITUDE", "START_LATITUDE")) |>
  sf::st_transform(crs = "EPSG:3338")

gebco_sf <- sf::st_as_sf(gebco_df, 
                         coords = c("x", "y"), 
                         crs = "EPSG:3338")

haul_loc$GEBCO_RASTER_INDEX <- st_nearest_feature(haul_loc, gebco_sf)
haul_loc$GEBCO_RASTER_DEPTH <- gebco_sf$layer[haul_loc$GEBCO_RASTER_INDEX]


png(here::here("plots", "EBS_GEBCO_vs_EBS_bathy_points.png"), width = 5, height = 5, units = "in", res = 300)
print(
  ggplot() +
    geom_point(data = haul_loc,
               mapping = aes(x = BOTTOM_DEPTH, y = GEBCO_RASTER_DEPTH)) +
    scale_x_continuous(name = "2022 EBS Shelf Survey Depth (m)") +
    scale_y_continuous(name = "GEBCO Depth (m)")
)
dev.off()
