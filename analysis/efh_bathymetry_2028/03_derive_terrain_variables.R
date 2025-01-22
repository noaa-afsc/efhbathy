# Generate terrain variables for 2028 EFH review
# Sean Rohan <sean.rohan@noaa.gov>
# January 7, 2025

library(akgfmaps)
library(MultiscaleDTM) # last run: version 0.8.3

dir.create(here::here("analysis", "efh_bathymetry_2028", "output", "terrain"))

# 2023 Review: Queen case (3x3 neighborhood) for deriving seafloor terrain
neighbors <- 8
bpi_w = c(3, 64) # 2023 Review: AI and GOA Radius 3 x 65; EBS radius 3 x 64


# Calculate terrain variables for 100-m resolution
bathy_100m <- terra::rast(here::here("analysis", "efh_bathymetry_2028", "output", "efh_bathy", "efh_bathy_100m.tif"))

# SLOPE using terra ----
start1 <- Sys.time()
slope_t_100m <- terra::terrain(x = bathy_100m, 
                          v = "slope", 
                          unit = "degrees", 
                          neighbors = neighbors)

# MultiscaleDTM: SLOPE (degrees), EASTNESS, AND NORTHNESS
# Queen case with window = 5
slope_aspect_dtm_100m <- MultiscaleDTM::SlpAsp(r = bathy_100m,
                                          unit = "degrees", # Output units
                                          w = c(5, 5),
                                          method = "queen",
                                          metrics = c("slope", "eastness", "northness"),
                                          na.rm = TRUE)

# MultiscaleDTM: MEAN CURVATURE, PROFILE CURVATURE, SLOPE
# Derived from quadratic fit using OLS
curvature_slope_dtm_100m <- MultiscaleDTM::Qfit(r = bathy_100m,
                                           unit = "degrees", # Output units
                                           w = c(5, 5),
                                           metrics = c("meanc", "profc",  "qslope"),
                                           na.rm = TRUE)

# Write 100-m rasters to GeoTIFF files ----

writeRaster(curvature_slope_dtm_100m$meanc, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_meancurv_100m.tif"), 
            overwrite = TRUE)
writeRaster(curvature_slope_dtm_100m$profc, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_profcurv_100m.tif"), 
            overwrite = TRUE)
writeRaster(curvature_slope_dtm_100m$qslope, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_qslope_100m.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_100m$northness, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_northness_100m.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_100m$eastness, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_eastness_100m.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_100m$slope, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_slope_100m.tif"), 
            overwrite = TRUE)
writeRaster(slope_t_100m,
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_t_slope_100m.tif"),
            overwrite = TRUE)


# Resample 100-m rasters to 1-km resolution ----
slope_t_1km <- terra::aggregate(slope_t_100m, fact = 10, fun = "mean")

slope_aspect_dtm_1km <- terra::aggregate(slope_aspect_dtm_100m, fact = 10, fun = "mean")

curvature_slope_dtm_1km <- terra::aggregate(curvature_slope_dtm_100m, fact = 10, fun = "mean")


# Write 1-km rasters to GeoTIFF files ----

writeRaster(curvature_slope_dtm_1km$meanc, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_meancurv_1km.tif"), 
            overwrite = TRUE)
writeRaster(curvature_slope_dtm_1km$profc, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_profcurv_1km.tif"), 
            overwrite = TRUE)
writeRaster(curvature_slope_dtm_1km$qslope, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_qslope_1km.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_1km$northness, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_northness_1km.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_1km$eastness, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_eastness_1km.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_1km$slope, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_slope_1km.tif"), 
            overwrite = TRUE)
writeRaster(slope_t_1km,
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_t_slope_1km.tif"),
            overwrite = TRUE)