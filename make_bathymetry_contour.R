library(akgfmaps)


#' Make bathymetry linestring from a raster

raster_to_bathy_contour <- function(x, 
                                    mask_sf = NULL, 
                                    mask_buffer_m = 0,
                                    mask_inverse = FALSE,
                                    scaling_factor = c(1, 1),
                                    scaling_method = "bilinear",
                                    depth_contours = c(seq(-1200, -100, 100), -50, -20), 
                                    min_length_m = 0,
                                    rms_keep = NULL, 
                                    rms_weighting = NULL,
                                    focal_mean_window = NULL,
                                    smooth_method = "ksmooth",
                                    smoothness = 2,
                                    output_crs = "WGS84") {
  
  # Mask the raster - this is extremely important if the raster includes land to avoid artifacts caused by grid cells that are above sea level
  if(!is.null(mask_sf)) {
    
    mask_sf <- sf::st_buffer(mask_sf, 
                             dist = units::as_units(mask_buffer_m))
    
    # Inverse masking is typically for masking with a land object. Using marine boundaries, inverse should be FALSE.
    x <- terra::mask(x, mask_sf, inverse = mask_inverse)
    
  }
  
  if(all(scaling_factor > 1) & !is.null(scaling_factor)) {
    message("Aggregating by a factor of ", paste(scaling_factor, collapse = " by "))
    x <- terra::aggregate(x, 
                          fact = scaling_factor, 
                          method = scaling_method, 
                          fun = "mean", 
                          dissolve = TRUE)
  }
  
  if(all(scaling_factor < 1) & !is.null(scaling_factor)) {
    x <- terra::disagg(x, 
                       fact = 1/scaling_factor, 
                       method = scaling_method)
  }
  
  if(any(scaling_factor) > 1 & any(scaling_factor) < 1 & !is.null(scaling_factor)) {
    stop("make_bathymetry_lines: Scaling factors cannot be greater and less than one." )
  }
  
  if(is.numeric(focal_mean_window)) {
    x <- terra::focal(x = x, w = focal_mean_window, fun = "mean")
  }
  
  run_smooth <- function(x, smooth_method, smoothness) {
    
    if(is.character(smooth_method) & is.numeric(smoothness)) {
      
      return(
        smoothr::smooth(x , 
                        method = "ksmooth", 
                        smoothness = 3)
        )
      
    }
    
    return(x)
    
  }
  
  # Downscaling to a higher resolution (smaller grid cells) or upscaling to a lower resolution (larger grid cells) can reduce the 'jaggedness' of lines.
  output <- terra::as.contour(x, 
                              levels = depth_contours) |>
    sf::st_as_sf(crs = output_crs) |> 
    run_smooth(smooth_method = "ksmooth", smoothness = smoothness) |>
    sf::st_cast(to = "MULTILINESTRING") |>
    sf::st_cast(to = "LINESTRING", 
                do_split = TRUE, 
                group_or_split = TRUE)
  
  # Remove short lines based on minimum length threshold
  output <- output[sf::st_length(output) > 
                     units::as_units(min_length_m, value = "m"), ]
  
  if(!is.null(rms_keep)) {
    output <- rmapshaper::ms_simplify(output, 
                                      keep = rms_keep, 
                                      weighting = rms_weighting)
  }
  
  names(output)[which(names(output) == "level")] <- "depth"
  
  return(output)
  
}


depth_breaks <- c(seq(10,100,10), seq(150, 500, 50), seq(600, 2000, 100))

bathy_options <- expand.grid(min_length_m = 5000,
                             focal_window = c(0, 3),
                             smooth_method = c("ksmooth", "none", "chaikin"),
                             smoothness = 1:3)

bathy_options$smooth_method[bathy_options$smooth_method == "none"] <- NULL

efh_100m <- terra::rast(here::here("output", "efh_bathy", "efh_bathy_100m.tif"))

test_nosmooth_fw3 <- raster_to_bathy_contour(x = efh_100m,
                      scaling_factor = NULL,
                      scaling_method = NULL,
                      depth_contours = depth_breaks,
                      min_length_m = 5000,
                      # rms_keep = 0.5,
                      # rms_weighting = 0.5,
                      focal_mean_window = 3,
                      output_crs = "EPSG:3338")

test_nosmooth <- raster_to_bathy_contour(x = efh_100m, 
                                 scaling_factor = c(1, 1),
                                 scaling_method = "near",
                                 depth_contours = depth_breaks, 
                                 min_length_m = 5000,
                                 # rms_keep = 0.5,
                                 # rms_weighting = 0.5,
                                 focal_mean_window = 3,
                                 output_crs = "EPSG:3338")

test_smooth2 <- raster_to_bathy_contour(x = efh_100m, 
                                 scaling_factor = c(1, 1),
                                 scaling_method = "near",
                                 depth_contours = depth_breaks, 
                                 min_length_m = 1e4,
                                 # rms_keep = 0.5,
                                 # rms_weighting = 0.5,
                                 smooth_method = "ksmooth",
                                 smoothness = 1,
                                 focal_mean_window = 3,
                                 output_crs = "EPSG:3338")

test_smooth3 <- raster_to_bathy_contour(x = efh_100m, 
                                        scaling_factor = c(1, 1),
                                        scaling_method = "near",
                                        depth_contours = depth_breaks, 
                                        min_length_m = 1e4,
                                        # rms_keep = 0.5,
                                        # rms_weighting = 0.5,
                                        smooth_method = "ksmooth",
                                        smoothness = 3,
                                        focal_mean_window = 3,
                                        output_crs = "EPSG:3338")

test2 <- raster_to_bathy_contour(x = efh_100m, 
                                scaling_factor = c(1, 1),
                                scaling_method = "bilinear",
                                depth_contours = depth_breaks, 
                                # min_length_m = 7000,
                                rms_keep = 0.07,
                                rms_weighting = 0.7,
                                focal_mean_window = 3,
                                output_crs = "EPSG:3338")

par(mfrow = c(1,2))
plot(test1)
plot(test2)

plot(test_smooth2)
plot(test_smooth3)

test <- raster_to_bathy_contour(x = efh_100m, 
                              mask_sf = NULL, 
                              mask_buffer_m = 0,
                              mask_inverse = FALSE,
                              scaling_factor = c(10, 10),
                              scaling_method = "bilinear",
                              depth_contours = c(10, 20, 50, seq(100, 1200, 100)), 
                              min_length_m = 5e4,
                              # rms_keep = 0.5, 
                              # rms_weighting = 0.5,
                              output_crs = "EPSG:3338")


crs(efh_100m, parse = TRUE)

sf::st_crs(test)$wkt
