library(akgfmaps)
library(MultiscaleDTM)

# 2023 Review: Queen case (3x3 neighborhood) for deriving seafloor terrain
neighbors <- 8
bpi_w = c(3, 64) # 2023 Review: AI and GOA Radius 3 x 65; EBS radius 3 x 64

x <- terra::rast(here::here("output", "efh_bathy", "efh_bathy_1km.tif"))

# SLOPE using terra ----
start1 <- Sys.time()
slope_t <- terra::terrain(x = x, 
                          v = "slope", 
                          unit = "degrees", 
                          neighbors = neighbors)

# MultiscaleDTM: SLOPE (degrees), EASTNESS, AND NORTHNESS
# Queen case with window = 5
slope_aspect_dtm <- MultiscaleDTM::SlpAsp(r = x,
                                          unit = "degrees", # Output units
                                          w = c(5, 5),
                                          method = "queen",
                                          metrics = c("slope", "eastness", "northness"),
                                          na.rm = TRUE)

# MultiscaleDTM: MEAN CURVATURE, PROFILE CURVATURE, SLOPE
# Derived from quadratic fit using OLS
curvature_slope_dtm <- MultiscaleDTM::Qfit(r = x,
                                           unit = "degrees", # Output units
                                           w = c(5, 5),
                                           metrics = c("meanc", "profc",  "qslope"),
                                           na.rm = TRUE)

# MultiscaleDTM: BATHYMETRIC POSITION INDEX (BPI) ----
start1 <- Sys.time()
bpi_dtm <- MultiscaleDTM::BPI(r = x,
                              w = bpi_w,
                              unit = "cell", # Units for window
                              stand = "none",
                              na.rm = TRUE
)

end1 <- Sys.time()
difftime(end1, start1)


# Test 100-m resolution

x_100 <- terra::rast(here::here("output", "efh_bathy", "efh_bathy_100m.tif"))

# SLOPE using terra ----
start2 <- Sys.time()
slope_t <- terra::terrain(x = x_100, 
                          v = "slope", 
                          unit = "degrees", 
                          neighbors = neighbors)

# MultiscaleDTM: SLOPE (degrees), EASTNESS, AND NORTHNESS
# Queen case with window = 5
slope_aspect_dtm <- MultiscaleDTM::SlpAsp(r = x_100,
                                          unit = "degrees", # Output units
                                          w = c(5, 5),
                                          method = "queen",
                                          metrics = c("slope", "eastness", "northness"),
                                          na.rm = TRUE)

# MultiscaleDTM: MEAN CURVATURE, PROFILE CURVATURE, SLOPE
# Derived from quadratic fit using OLS
curvature_slope_dtm <- MultiscaleDTM::Qfit(r = x_100,
                                           unit = "degrees", # Output units
                                           w = c(5, 5),
                                           metrics = c("meanc", "profc",  "qslope"),
                                           na.rm = TRUE)
diff(Sys.time(), start2)

# MultiscaleDTM: BATHYMETRIC POSITION INDEX (BPI) ----
start1 <- Sys.time()
bpi_dtm <- MultiscaleDTM::BPI(r = x_100,
                              w = bpi_w,
                              unit = "cell", # Units for window
                              stand = "none",
                              na.rm = TRUE
)

end2 <- Sys.time()
difftime(end2, start2)