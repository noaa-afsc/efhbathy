# Generate bathymetry raster for EFH

# Software versions ----
# R version: 4.4.1
# Python 3.9
# ArcGIS Pro Version 13.1.0.41833
# arcgisbinding version 1.0.1.311
# Reticulate version 1.35.0

library(akgfmaps) # Version 4.0
library(arcgisbinding) # Version 1.0.1.311: Installed for ArcGIS Pro 13.1.0.41833
library(reticulate) # Version 1.35.0
library(sf)
library(terra)

# Alaska Albers Equal Area
aea_wkt <- 'PROJCS["NAD_1983_Albers",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Albers"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-154.0],PARAMETER["Standard_Parallel_1",55.0],PARAMETER["Standard_Parallel_2",65.0],PARAMETER["Latitude_Of_Origin",50.0],UNIT["Meter",1.0]]'

# Get akgfmaps layers for masking
map_layers <- 
  akgfmaps::get_base_layers(
  select.region = "ebs", 
  set.crs = aea_wkt, 
  high.resolution.coast = TRUE # Use the highest resolution coastline
)

# Setup directories for bathymetry data and output
dir.create(
  here::here("analysis", "efh_bathymetry_2028", "data", "efh_bathy"), 
  recursive = TRUE, 
  showWarnings = FALSE
)

dir.create(
  here::here("analysis", "efh_bathymetry_2028", "output", "efh_bathy"), 
  recursive = TRUE, 
  showWarnings = FALSE
)

# Define ArcGIS Pro license ----
arcgisbinding::arc.check_product()

# Set Python environment - Python 3.9 with arcpy already installed
# Command line install: conda install arcpy=3.1 -c esri
reticulate::use_condaenv("C:\\Users\\sean.rohan\\.conda\\envs\\py39")

# Load arcpy library
ARCPY <- reticulate::import("arcpy")

# Remove eastern portion of AI grid that extends beyond interpolated area
ai_grid <- 
  terra::rast(
    here::here("analysis", "efh_bathymetry_2028", "data", "ai_grid_100m")
  )

unimak_bounds <- 
  sf::st_polygon(
    list(
      matrix(
        c(-167, 40,
          -167, 54.2,
          -177, 53.5,
          -179, 60,
          -150, 60,
          -150, 40, 
          -167, 40),
        ncol = 2,
        byrow = TRUE
      )
    )
  ) |>
  sf::st_sfc(crs = "NAD83") |>
  sf::st_as_sf() |>
  sf::st_transform(crs = aea_wkt)

ai_grid <- 
  terra::mask(
    x = ai_grid, 
    mask = unimak_bounds, 
    inverse = TRUE
    ) |>
  terra::trim()

terra::writeRaster(
  x = ai_grid,
  filename = here::here("analysis", "efh_bathymetry_2028", "data", "ai_100m_masked.tif"),
  overwrite = TRUE
)

# Fix major depth errors in 2024 GEBCO bathymetry (unreasonably large depths in the EBS inner domain)

gebco_inverse <- terra::rast(here::here("analysis", "efh_bathymetry_2028", "data", "gebco_2024_inverse.tif"))

gebco_mask <- gebco_inverse
gebco_mask[!is.na(values(gebco_mask))] <- 1

gebco_ebs_shelf_10_20 <- 
  terra::mask(
    x = gebco_inverse, 
    mask = sf::st_transform(
      dplyr::filter(
        map_layers$survey.strata, 
        Stratum %in% c(10, 20)),
      crs = terra::crs(gebco_inverse)
    )
  )

gebco_ebs_shelf_10_20[values(gebco_ebs_shelf_10_20) < 100] <- NA

gebco_ebs_shelf_10_20[!is.na(values(gebco_ebs_shelf_10_20))] <- -999

gebco_ebs_shelf_10_20 <- 
  as.polygons(gebco_ebs_shelf_10_20) |>
  sf::st_as_sf() |>
  sf::st_buffer(dist = 5000)

gebco_inverse1 <- 
  terra::mask(
    x = gebco_inverse, 
    mask = gebco_ebs_shelf_10_20, 
    inverse = TRUE
  )

gebco_inverse2 <- 
  terra::focal(
    gebco_inverse1, 
    w = 7, 
    fun = mean, 
    na.policy = "only", 
    na.rm = TRUE
  )

for(ii in 1:5) {
  gebco_inverse2 <- 
    terra::focal(
      gebco_inverse2, 
      w = 7, 
      fun = mean, 
      na.policy = "only", 
      na.rm = TRUE
    )
}

plot(gebco_inverse2 * gebco_mask)

terra::writeRaster(
  gebco_inverse2,
  filename = here::here("analysis", "efh_bathymetry_2028", "data", "gebco_2024_inverse_filled.tif")
)


# Set input file paths for bathymetry rasters files. These should be entered in order of quality, 
# with the highest quality raster first, second highest quality second, etc. The order of these
# rasters determines which layer gets used in the mosaic process.
in_paths <- c(
  here::here("analysis", "efh_bathymetry_2028", "data", "goa_bathy"),
  here::here("analysis", "efh_bathymetry_2028", "data", "ai_100m_masked.tif"),
  here::here("analysis", "efh_bathymetry_2028", "data", "Norton_Bathymetry", "norton_final2"),
  here::here("analysis", "efh_bathymetry_2028", "data", "bs_bathy"),
  here::here("analysis", "efh_bathymetry_2028", "data", "gebco_2024_inverse_filled.tif")
)

# Set output locations 
out_dir_path <- here::here("analysis", "efh_bathymetry_2028", "data")

filename_mosaic <- "efh_mosaic_100m.tif"
filename_mosaic_resample <- "efh_mosaic_resample_100m.tif"
filename_tin <- "efh_tin_100m"
filename_masked <- "efh_masked_100m.tif"
filename_raw <- "efh_raw_100m.tif"
filename_output_100m <- "efh_bathy_100m.tif"
filename_output_1km <- "efh_bathy_1km.tif"

# Set valid depth range for pixels to be used for interpolation, in meters
bathy_depth_range <- c(0, 12000)

# Write land mask shapefile so it can be used by arcpy
sf::st_write(
  obj = map_layers$akland, 
  dsn = here::here("analysis", 
                   "efh_bathymetry_2028", 
                   "data", 
                   "akland_temp_mask.shp"),
  append = FALSE
)

# Make 32-bit float raster mosaic, where values are set in the order of paths in in_path
# This creates a 100 x 100 m resolution raster from the input files using a mosaic method where 
# all of the rasters are used to fill in the coverage. In areas where layers overlap, values from the
# first layer in in_paths are used ahead of the second layer, the second layer ahead of the third,
# etc.
ARCPY$MosaicToNewRaster_management(
  input_rasters = paste(in_paths, collapse = ";"),
  output_location = out_dir_path,
  raster_dataset_name_with_extension = filename_mosaic,
  coordinate_system_for_the_raster = aea_wkt,
  pixel_type = "32_BIT_FLOAT",
  cellsize = NULL,
  number_of_bands = 1,
  mosaic_method = "FIRST",
  mosaic_colormap_mode = "FIRST"
  )

# Mask with coastline to remove land pixels
out_raster <- ARCPY$sa$ExtractByMask(
  in_raster = here::here(out_dir_path, filename_mosaic),
  in_mask_data = here::here("analysis", 
                            "efh_bathymetry_2028", 
                            "data", 
                            "akland_temp_mask.shp"),
  extraction_area = "OUTSIDE",
  analysis_extent = here::here(out_dir_path, filename_mosaic)
)

if(file.exists(here::here("analysis", "efh_bathymetry_2028", "data", filename_masked))) {
  ARCPY$Delete_management(
    here::here("analysis", "efh_bathymetry_2028", "data", filename_masked)
  )
}

out_raster$save(here::here("analysis", "efh_bathymetry_2028", "data", filename_masked))

# Remove pixels that are outside of a specified depth range
out_raster <- 
  ARCPY$sa$SetNull(
    in_conditional_raster = out_raster,
    in_false_raster_or_constant = here::here("analysis", "efh_bathymetry_2028", "data", filename_masked),
    where_clause = paste0("VALUE < ", bathy_depth_range[1], " Or VALUE >", bathy_depth_range[2]) 
  )

# Remove duplicate file
if(file.exists(here::here("analysis", "efh_bathymetry_2028", "data", filename_raw ))) {
  ARCPY$Delete_management(
    here::here("analysis", "efh_bathymetry_2028", "data", filename_raw)
    )
}

out_raster$save(
  here::here("analysis", "efh_bathymetry_2028", "data", filename_raw)
)

# Trim blank space outside of the polygons and rename the depth band in the raster
efh_bathy_trim <- terra::rast(
  here::here("analysis", "efh_bathymetry_2028", "data", filename_raw)
  ) |>
  terra::trim()

names(efh_bathy_trim) <- "depth"

terra::writeRaster(
  x = efh_bathy_trim, 
  filename = here::here("analysis", "efh_bathymetry_2028", "output", "efh_bathy", filename_output_100m), 
  overwrite = TRUE
)

# Resample the 100 m resolution grid to 1 km
efh_1km <- 
  terra::aggregate(
    efh_bathy_trim, 
    fact = 10, 
    fun = "mean"
  )

terra::writeRaster(
  x = efh_1km, 
  filename = here::here("analysis", "efh_bathymetry_2028", "output", "efh_bathy", filename_output_1km), 
  overwrite = TRUE
)

# Remove intermediate files
ARCPY$Delete_management(here::here("analysis", "efh_bathymetry_2028", "data", filename_raw))
ARCPY$Delete_management(here::here("analysis", "efh_bathymetry_2028", "data", filename_masked))
ARCPY$Delete_management(here::here("analysis", "efh_bathymetry_2028", "data", filename_mosaic))
ARCPY$Delete_management(here::here("analysis", "efh_bathymetry_2028", "data", "akland_temp_mask"))
ARCPY$Delete_management(here::here("analysis", "efh_bathymetry_2028", "data", here::here("analysis", "efh_bathymetry_2028", "data", "ai_100m_masked.tif")))
ARCPY$Delete_management(here::here("analysis", "efh_bathymetry_2028", "data", filename_tin))
ARCPY$Delete_management(here::here("analysis", "efh_bathymetry_2028", "data", "ai_100m_masked.tif"))
ARCPY$Delete_management(here::here("analysis", "efh_bathymetry_2028", "data", "gebco_2024_inverse_filled.tif"))

