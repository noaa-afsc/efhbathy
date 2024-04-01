# Generate bathymetry raster for EFH

library(akgfmaps) # Version 4.0: install_github("afsc-gap-products/akgfmaps@shp2025")
library(arcgisbinding) # Version 1.0.1.306: Install for a specific R version through ArcGIS Pro
library(reticulate) # Version 1.35.0

# Setup directory
dir.create(here::here("analysis", "efh_bathymetry_2028", "data", "efh_bathy"), recursive = TRUE)
dir.create(here::here("analysis", "efh_bathymetry_2028", "output", "efh_bathy"), recursive = TRUE)

# Define ArcGIS license
arcgisbinding::arc.check_product()

# Set environment - Python 3.9 with arcpy already installed
# Command line install: conda install arcpy=3.1 -c esri
reticulate::use_condaenv("C:\\Users\\sean.rohan\\.conda\\envs\\py39")

# Load arcpy library
ARCPY <- reticulate::import("arcpy")

# Set input file paths - these should be entered in order of quality, with the best or most recent raster first.
in_paths <- c(
  here::here("analysis", "efh_bathymetry_2028", "data", "goa_bathy"),
  here::here("analysis", "efh_bathymetry_2028", "data", "ai_grid_100m", "ai_grid_100m"),
  here::here("analysis", "efh_bathymetry_2028", "data", "Norton_Bathymetry", "norton_final2"),
  here::here("analysis", "efh_bathymetry_2028", "data", "bs_bathy"),
  here::here("analysis", "efh_bathymetry_2028", "data", "gebco_2023_inverse.tif")
)

# Set output locations 
out_dir_path <- here::here("analysis", "efh_bathymetry_2028", "data")

filename_mosaic <- "efh_mosaic_100m.tif"
filename_masked <- "efh_masked_100m.tif"
filename_raw <- "efh_raw_100m.tif"
filename_output_100m <- "efh_bathy_100m.tif"
filename_output_1km <- "efh_bathy_1km.tif"

bathy_depth_range <- c(0, 2000)

aea_wkt <- 'PROJCS["NAD_1983_Albers",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Albers"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-154.0],PARAMETER["Standard_Parallel_1",55.0],PARAMETER["Standard_Parallel_2",65.0],PARAMETER["Latitude_Of_Origin",50.0],UNIT["Meter",1.0]]'


# Get akgfmaps layers for masking
map_layers <- akgfmaps::get_base_layers(select.region = "ebs", set.crs = aea_wkt)

# Write land mask shapefile so it can be used by arcpy
sf::st_write(obj = map_layers$akland, 
             dsn = here::here("analysis", 
                              "efh_bathymetry_2028", 
                              "data", 
                              "akland_temp_mask.shp"),
             append = FALSE
)

# Make 32-bit float raster mosaic, where values are set in the order of paths in in_path 
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
  analysis_extent = NULL
)

if(file.exists(here::here("analysis", "efh_bathymetry_2028", "data", filename_masked))) {
  ARCPY$Delete_management(here::here("analysis", "efh_bathymetry_2028", "data", filename_masked))
}

out_raster$save(here::here("analysis", "efh_bathymetry_2028", "data", filename_masked))

# Remove pixels that are outside of a specified depth range
out_raster <- ARCPY$sa$SetNull(
  in_conditional_raster = out_raster,
  in_false_raster_or_constant = here::here("analysis", "efh_bathymetry_2028", "data", filename_masked),
  where_clause = paste0("VALUE < ", bathy_depth_range[1], " Or VALUE >", bathy_depth_range[2]) 
)

# Remove duplicate file
if(file.exists(here::here("analysis", "efh_bathymetry_2028", "data", filename_raw ))) {
  ARCPY$Delete_management(here::here("analysis", "efh_bathymetry_2028", "data", filename_raw))
}

out_raster$save(here::here("analysis", "efh_bathymetry_2028", "data", filename_raw))

# Trim blank space and rename the depth band
efh_bathy_trim <- terra::rast(here::here("analysis", "efh_bathymetry_2028", "data", filename_raw)) |>
  terra::trim()

names(efh_bathy_trim) <- "depth"

terra::writeRaster(x = efh_bathy_trim, 
                   filename = here::here("analysis", "efh_bathymetry_2028", "output", "efh_bathy", filename_output_100m), 
                   overwrite = TRUE)

# Resample the 100 m resolution grid to 1 km
efh_1km <- terra::aggregate(efh_bathy_trim, fact = 10, fun = "mean")

terra::writeRaster(x = efh_1km, 
                   filename = here::here("analysis", "efh_bathymetry_2028", "output", "efh_bathy", filename_output_1km), 
                   overwrite = TRUE)

# Remove intermediate files
ARCPY$Delete_management(here::here("analysis", "efh_bathymetry_2028", "data", filename_raw))
ARCPY$Delete_management(here::here("analysis", "efh_bathymetry_2028", "data", filename_masked))
ARCPY$Delete_management(here::here("analysis", "efh_bathymetry_2028", "data", filename_mosaic))
ARCPY$Delete_management(here::here("analysis", "efh_bathymetry_2028", "data", "akland_temp_mask"))
