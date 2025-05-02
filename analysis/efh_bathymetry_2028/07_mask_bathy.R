# Mask bathymetry and terrain variable rasters to the extent of survey regions
# Created by Sean Rohan

library(akgfmaps) # akgfmaps v4.0.5
library(terra)
library(here)

# Specify which grid files to use and the names for the associated covariate ----
# 3x3 window 1 km terrain variables and 1.5/6.4 km inner/outer annulus BPI
layer_files <-
  c(
    "depth" = "efh_bathy_1km.tif",
    "slope" = "efh_dtm_slope_1km_3_3.tif",
    "meancurv" = "efh_dtm_meancurv_1km_3_3.tif",
    "northness" = "efh_dtm_northness_1km_3_3.tif",
    "eastness" = "efh_dtm_eastness_1km_3_3.tif",
    "bpi" = "efh_bpi_15_64_1km.tif"
  )

layer_files[] <- paste0(layer_files, "$")

layer_files <- 
  sapply(
  X = layer_files, 
      FUN = list.files, 
      path = here::here("analysis", "efh_bathymetry_2028", "output"), 
      recursive = TRUE, 
      full.names = TRUE
  )

# Load files into a raster stack ----
raster_stack <- terra::rast(layer_files)

# Load survey area boundaries from akgfmaps ----
survey_layers <- 
  akgfmaps::get_base_layers(
  select.region = c("goa", "ebs", "ai"), 
  set.crs = "EPSG:3338",
  high.resolution.coast = TRUE # High resolution coastline for masking
)

# Combine EBS and NBS polygons into a single region for masking ----
region_boundaries <- survey_layers$survey.area
region_boundaries$region[region_boundaries$SURVEY_DEFINITION_ID == 52] <- "AI"
region_boundaries$region[region_boundaries$SURVEY_DEFINITION_ID == 47] <- "GOA"
region_boundaries$region[region_boundaries$SURVEY_DEFINITION_ID %in% c(98, 143)] <- "EBS"

region_boundaries <-
  region_boundaries |>
  dplyr::select(region) |>
  dplyr::group_by(region) |>
  dplyr::summarise(do_union = TRUE)

# Create an .rda file containing terrain variables for each region ----
for(ii in c("AI", "GOA", "EBS")) {
  
  # Mask terrain variables by region and land
  region_vars <- 
    raster_stack |>
    terra::mask(
      region_boundaries[region_boundaries$region == ii, ],
      touches = TRUE
    ) |>
    terra::mask(
      survey_layers$akland,
      inverse = TRUE,
      touches = FALSE) |> # Inverse mask w/ land
    terra::trim() # Remove empty area around the raster
  
  names(region_vars) <- names(raster_stack)
  
  # Pack SpatRaster so it can be written to a .rda file
  assign(
    x = paste0("static_variables_", tolower(ii)), 
    value = terra::wrap(region_vars)
    )
  
  # Save wrapped PackedSpatRaster to .rda
  save(
    list = paste0("static_variables_", tolower(ii)),
    file = here::here(
      "analysis", 
      "efh_bathymetry_2028", 
      "output", 
      paste0("efh_static_variables_", ii, ".rda")
    )
  )
  
}

# Test SpatRaster
load(
  here::here(
  "analysis",
  "efh_bathymetry_2028",
  "output",
  paste0("efh_static_variables_", ii, ".rda")
))
