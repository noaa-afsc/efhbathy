library(akgfmaps)
library(reticulate)


dir.create(here::here("analysis", "efh_bathymetry_2028", "data", "efh_bathy"), recursive = TRUE)
dir.create(here::here("analysis", "efh_bathymetry_2028", "output", "efh_bathy"), recursive = TRUE)

# Transform GEBCO Z coordinates to match Alaska bathymetry (negative is up)
gebco <- terra::rast(here::here("analysis", "efh_bathymetry_2028", "data", "gebco_2023_n67.0_s54.0_w-175.0_e-154.0.tif"))
gebco <- gebco * -1

map_layers <- akgfmaps::get_base_layers(select.region = "sebs", set.crs = crs(gebco, proj = TRUE))

gebco <- terra::mask(gebco, map_layers$akland[2,], inverse = TRUE)

terra::writeRaster(gebco, filename = here::here("analysis", "efh_bathymetry_2028", "data", "gebco_2023_inverse.tif"))


# ArcPy
# Python 3.9 with ArcGIS pro version 3.1

reticulate::use_python("C:/Program Files/Python39/python.exe")
reticulate::source_python(file = here::here("analysis", "efh_bathymetry_2028", "make_efh_bathy.py"))


## Trim blank areas around the edge and rename the depth band

efh_bathy_trim <- terra::rast(here::here("analysis", "efh_bathymetry_2028", "data", "efh_raw_100m.tif")) |>
  terra::trim()

names(efh_bathy_trim) <- "depth"

terra::writeRaster(efh_bathy_trim, here::here("analysis", "efh_bathymetry_2028", "output", "efh_bathy", "efh_bathy_100m.tif"), overwrite = TRUE)


# Resample the 100 m resolution grid to 1 km
efh_100m <- terra::rast(here::here("analysis", "efh_bathymetry_2028", "output", "efh_bathy", "efh_bathy_100m.tif"))

efh_1km <- terra::aggregate(efh_100m, fact = 10, fun = "mean")

terra::writeRaster(efh_1km, here::here("analysis", "efh_bathymetry_2028", "output", "efh_bathy", "efh_bathy_1km.tif"), overwrite = TRUE)








