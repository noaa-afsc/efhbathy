# Generate BPI for 2028 EFH review
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

# MultiscaleDTM: BATHYMETRIC POSITION INDEX (BPI) ----
start1 <- Sys.time()

bpi_dtm_15_64_100m <- MultiscaleDTM::BPI(r = bathy_100m,
                                         w = c(15, 64),
                                         unit = "cell", # Units for window
                                         stand = "none",
                                         na.rm = TRUE
)

start2 <- Sys.time()
difftime(start2, start1)

start2 <- Sys.time()
writeRaster(bpi_dtm_15_64_100m, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_bpi_15_64_100m.tif"), 
            overwrite = TRUE)

bpi_dtm_30_124_100m <- MultiscaleDTM::BPI(r = bathy_100m,
                                          w = c(30, 124),
                                          unit = "cell", # Units for window
                                          stand = "none",
                                          na.rm = TRUE
)

start3 <- Sys.time()
difftime(start3, start2)

writeRaster(bpi_dtm_30_124_100m, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_bpi_30_124_100m.tif"), 
            overwrite = TRUE)

bpi_dtm_3_64_100m <- MultiscaleDTM::BPI(r = bathy_100m,
                                        w = bpi_w,
                                        unit = "cell", # Units for window
                                        stand = "none",
                                        na.rm = TRUE
)

writeRaster(bpi_dtm_3_64_100m, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_bpi_3_64_100m.tif"), 
            overwrite = TRUE)

start4 <- Sys.time()
difftime(start4, start3)

# Resample 100-m rasters to 1-km resolution ----
bpi_dtm_3_64_1km <- terra::aggregate(bpi_dtm_3_64_100m, fact = 10, fun = "mean")

bpi_dtm_15_64_1km <- terra::aggregate(bpi_dtm_15_64_100m, fact = 10, fun = "mean")

bpi_dtm_30_124_1km <- terra::aggregate(bpi_dtm_30_124_100m, fact = 10, fun = "mean")

# Write 1-km rasters to GeoTIFF files ----

writeRaster(bpi_dtm_3_64_1km, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_bpi_3_64_1km.tif"), 
            overwrite = TRUE)
writeRaster(bpi_dtm_15_64_1km, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_bpi_15_64_1km.tif"), 
            overwrite = TRUE)
writeRaster(bpi_dtm_30_124_1km, 
            here::here("analysis", "efh_bathymetry_2028", "output", "terrain", "efh_bpi_30_124_1km.tif"), 
            overwrite = TRUE)