# Generate terrain variables for 2028 EFH review
# Sean Rohan <sean.rohan@noaa.gov>
# February 8, 2025

library(akgfmaps)
library(MultiscaleDTM) # last run: version 0.8.3

dir.create(here::here("analysis", "efh_bathymetry_2028", "output", "terrain"))

# Calculate terrain variables for 100-m resolution
bathy_100m <- terra::rast(here::here("analysis", "efh_bathymetry_2028", "output", "efh_bathy", "efh_bathy_100m.tif"))

# SLOPE using terra ----
# 2023 Review: Queen case (3x3 neighborhood) for deriving seafloor terrain
start1 <- Sys.time()
slope_t_100m <- terra::terrain(x = bathy_100m, 
                               v = "slope", 
                               unit = "degrees", 
                               neighbors = 8)

writeRaster(slope_t_100m,
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_t_slope_100m.tif"),
            overwrite = TRUE)

slope_t_1km <- terra::aggregate(slope_t_100m, fact = 10, fun = "mean")

writeRaster(slope_t_1km,
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_t_slope_1km.tif"),
            overwrite = TRUE)

# MultiscaleDTM: SLOPE (degrees), EASTNESS, AND NORTHNESS
# Queen case with windows of 3 and 5
slope_aspect_dtm_100m_3_3 <- MultiscaleDTM::SlpAsp(r = bathy_100m,
                                          unit = "degrees", # Output units
                                          w = c(3, 3),
                                          method = "queen",
                                          metrics = c("slope", "eastness", "northness"),
                                          na.rm = TRUE)

slope_aspect_dtm_100m_5_5 <- MultiscaleDTM::SlpAsp(r = bathy_100m,
                                                   unit = "degrees", # Output units
                                                   w = c(5, 5),
                                                   method = "queen",
                                                   metrics = c("slope", "eastness", "northness"),
                                                   na.rm = TRUE)
start2 <- Sys.time()
difftime(start2, start1)

# MultiscaleDTM: MEAN CURVATURE, PROFILE CURVATURE, SLOPE
# Derived from ordnary least squares quadratic fit
# Queen case with windows of 3 and 5
curvature_slope_dtm_100m_3_3 <- MultiscaleDTM::Qfit(r = bathy_100m,
                                           unit = "degrees", # Output units
                                           w = c(3, 3),
                                           metrics = c("meanc", "profc",  "qslope"),
                                           na.rm = TRUE)

curvature_slope_dtm_100m_3_3 <- MultiscaleDTM::Qfit(r = bathy_100m,
                                                    unit = "degrees", # Output units
                                                    w = c(3, 3),
                                                    metrics = c("meanc", "profc",  "qslope"),
                                                    na.rm = TRUE)

curvature_slope_dtm_100m_5_5 <- MultiscaleDTM::Qfit(r = bathy_100m,
                                                    unit = "degrees", # Output units
                                                    w = c(5, 5),
                                                    metrics = c("meanc", "profc",  "qslope"),
                                                    na.rm = TRUE)

curvature_slope_dtm_100m_5_5 <- MultiscaleDTM::Qfit(r = bathy_100m,
                                                    unit = "degrees", # Output units
                                                    w = c(5, 5),
                                                    metrics = c("meanc", "profc",  "qslope"),
                                                    na.rm = TRUE)

start3 <- Sys.time()
difftime(start3, start2)

# Write 3x3 window 100-m rasters to GeoTIFF files ----

writeRaster(curvature_slope_dtm_100m_3_3$meanc, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_meancurv_100m_3_3.tif"), 
            overwrite = TRUE)
writeRaster(curvature_slope_dtm_100m_3_3$profc, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_profcurv_100m_3_3.tif"), 
            overwrite = TRUE)
writeRaster(curvature_slope_dtm_100m_3_3$qslope, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_qslope_100m_3_3.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_100m_3_3$northness, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_northness_100m_3_3.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_100m_3_3$eastness, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_eastness_100m_3_3.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_100m_3_3$slope, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_slope_100m_3_3.tif"), 
            overwrite = TRUE)

# Write 5x5 window 100-m rasters to GeoTIFF files ----

writeRaster(curvature_slope_dtm_100m_5_5$meanc, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_meancurv_100m_5_5.tif"), 
            overwrite = TRUE)
writeRaster(curvature_slope_dtm_100m_5_5$profc, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_profcurv_100m_5_5.tif"), 
            overwrite = TRUE)
writeRaster(curvature_slope_dtm_100m_5_5$qslope, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_qslope_100m_5_5.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_100m_5_5$northness, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_northness_100m_5_5.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_100m_5_5$eastness, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_eastness_100m_5_5.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_100m_5_5$slope, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_slope_100m_5_5.tif"), 
            overwrite = TRUE)


# Resample 100-m rasters to 1-km resolution ----
slope_aspect_dtm_1km_3_3 <- terra::aggregate(slope_aspect_dtm_100m_3_3, fact = 10, fun = "mean")

curvature_slope_dtm_1km_3_3 <- terra::aggregate(curvature_slope_dtm_100m_3_3, fact = 10, fun = "mean")

slope_aspect_dtm_1km_5_5 <- terra::aggregate(slope_aspect_dtm_100m_5_5, fact = 10, fun = "mean")

curvature_slope_dtm_1km_5_5 <- terra::aggregate(curvature_slope_dtm_100m_5_5, fact = 10, fun = "mean")


# Write 3x3 1-km rasters to GeoTIFF files ----

writeRaster(curvature_slope_dtm_1km_3_3$meanc, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_meancurv_1km_3_3.tif"), 
            overwrite = TRUE)
writeRaster(curvature_slope_dtm_1km_3_3$profc, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_profcurv_1km_3_3.tif"), 
            overwrite = TRUE)
writeRaster(curvature_slope_dtm_1km_3_3$qslope, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_qslope_1km_3_3.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_1km_3_3$northness, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_northness_1km_3_3.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_1km_3_3$eastness, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_eastness_1km_3_3.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_1km_3_3$slope, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_slope_1km_3_3.tif"), 
            overwrite = TRUE)

# Write 5x5 1-km rasters to GeoTIFF files ----
writeRaster(curvature_slope_dtm_1km_5_5$meanc, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_meancurv_1km_5_5.tif"), 
            overwrite = TRUE)
writeRaster(curvature_slope_dtm_1km_5_5$profc, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_profcurv_1km_5_5.tif"), 
            overwrite = TRUE)
writeRaster(curvature_slope_dtm_1km_5_5$qslope, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_qslope_1km_5_5.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_1km_5_5$northness, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_northness_1km_5_5.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_1km_5_5$eastness, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_eastness_1km_5_5.tif"), 
            overwrite = TRUE)
writeRaster(slope_aspect_dtm_1km_5_5$slope, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_dtm_slope_1km_5_5.tif"), 
            overwrite = TRUE)
