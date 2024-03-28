import arcpy
from arcpy import env

env.workspace = "C:/Users/sean.rohan/Work/afsc/EFH/efh_bathy/data"

arcpy.management.MosaicToNewRaster(
    input_rasters="goa_bathy;ai_grid_100m;norton_final2;bs_bathy;gebco_2023_inverse.tif",
    output_location=r"C:\Users\sean.rohan\Work\afsc\EFH\efh_bathy\data",
    raster_dataset_name_with_extension="efh_mosaic_100m.tif",
    coordinate_system_for_the_raster='PROJCS["NAD_1983_Albers",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Albers"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-154.0],PARAMETER["Standard_Parallel_1",55.0],PARAMETER["Standard_Parallel_2",65.0],PARAMETER["Latitude_Of_Origin",50.0],UNIT["Meter",1.0]]',
    pixel_type="32_BIT_FLOAT",
    cellsize=None,
    number_of_bands=1,
    mosaic_method="FIRST",
    mosaic_colormap_mode="FIRST"
)

out_raster = arcpy.sa.ExtractByMask(
    in_raster="efh_mosaic_100m.tif",
    in_mask_data="Alaska_Coastline",
    extraction_area="OUTSIDE",
    analysis_extent='-2625085.4081599 -217033.970228436 1614614.5918401 2037966.02977156 PROJCS["NAD_1983_Albers",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Albers"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["central_meridian",-154.0],PARAMETER["Standard_Parallel_1",55.0],PARAMETER["Standard_Parallel_2",65.0],PARAMETER["latitude_of_origin",50.0],UNIT["Meter",1.0]]'
)

out_raster = arcpy.sa.ExtractByMask(
    in_raster=out_raster,
    in_mask_data="e_russia",
    extraction_area="OUTSIDE",
    analysis_extent='-2625085.4081599 -217033.970228436 1614614.5918401 2037966.02977156 PROJCS["NAD_1983_Albers",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Albers"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["central_meridian",-154.0],PARAMETER["Standard_Parallel_1",55.0],PARAMETER["Standard_Parallel_2",65.0],PARAMETER["latitude_of_origin",50.0],UNIT["Meter",1.0]]'
)

out_raster = arcpy.sa.SetNull(
    in_conditional_raster=out_raster,
    in_false_raster_or_constant="efh_mosaic_russmask_100m.tif",
    where_clause="VALUE < 0 Or VALUE > 2000"
)
out_raster.save(r"C:\Users\sean.rohan\Work\afsc\EFH\efh_bathy\data\efh_raw_100m.tif")