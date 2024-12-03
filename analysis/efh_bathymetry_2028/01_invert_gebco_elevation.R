library(akgfmaps)

dir.create(here::here("analysis", "efh_bathymetry_2028", "data", "efh_bathy"), recursive = TRUE)
dir.create(here::here("analysis", "efh_bathymetry_2028", "output", "efh_bathy"), recursive = TRUE)

# Transform GEBCO Z coordinates to match Alaska bathymetry (negative is up)
gebco <- terra::rast(here::here("analysis", "efh_bathymetry_2028", "data", "gebco_2024_n67.0_s54.0_w-175.0_e-154.0.tif"))
gebco <- gebco * -1

map_layers <- akgfmaps::get_base_layers(select.region = "sebs", set.crs = crs(gebco, proj = TRUE))

gebco <- terra::mask(gebco, map_layers$akland[2,], inverse = TRUE)

terra::writeRaster(gebco, filename = here::here("analysis", "efh_bathymetry_2028", "data", "gebco_2024_inverse.tif"))
