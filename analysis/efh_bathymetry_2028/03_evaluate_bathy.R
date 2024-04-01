library(akgfmaps)

# Setup
dir.create(here::here("analysis", "efh_bathymetry_2028", "plots", "bathy"), recursive = TRUE)

region <- c("goa", "ai", "ebs")

theme_2 <- function() {
  theme_bw() %+replace%
    theme(legend.title = element_text(size = 7, color = "black"),
          panel.border = element_rect(color = "black", fill = NA),
          legend.text = element_text(size = 7, color = "black"),
          legend.key = element_rect(fill = NA, color = NA),
          axis.title = element_blank(),
          axis.text = element_text(size = 6, color = "black"),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA))
}

# Get bottom depth at tow midpoints
source(here::here("analysis", "efh_bathymetry_2028", "get_towmid.R"))


# Comparison metrics
mean_relative_error <- function(x, y) {
  mean(abs((x-y))/y, na.rm = TRUE)
}

mean_absolute_error <- function(x, y) {
  mean(abs(log10(x)-log10(y)), na.rm = TRUE)
}

mean_bias <- function(x, y) {
  10^mean(log10(x)-log10(y), na.rm = TRUE)
}


# Comparison code

for(ii in 1:length(region)) {
  
  old_bathy <- terra::rast(here::here("analysis", "efh_bathymetry_2028", "data", "2022_efh_bathy", 
                                      paste0("Variables_", region[ii], "_1km"), "Bathy.grd"))
  
  new_bathy <- terra::rast(here::here("analysis", "efh_bathymetry_2028", "output", "efh_bathy", "efh_bathy_1km.tif")) 
  
  new_bathy_resampled <- terra::resample(new_bathy, old_bathy, method = "bilinear")
  
  old_bathy_reproj <- terra::project(old_bathy, new_bathy_resampled)
  
  diff_bathy <- terra::trim(new_bathy_resampled - old_bathy_reproj)
  
  map_layers <- akgfmaps::get_base_layers(select.region = region[ii], set.crs = terra::crs(diff_bathy))
  
  plot_cont_bathy_diff <- ggplot() +
    geom_stars(data = stars::st_as_stars(diff_bathy)) +
    geom_sf(data = map_layers$akland) +
    coord_sf(xlim = map_layers$plot.boundary$x, 
             ylim = map_layers$plot.boundary$y) +
    scale_fill_distiller(name = expression(Z[new]-Z[2022]), 
                         palette = "BrBG",
                         na.value = NA) +
    scale_x_continuous(breaks = map_layers$lon.breaks) +
    scale_y_continuous(breaks = map_layers$lat.breaks) +
    theme_2()
  
  plot_disc_bathy_diff <- ggplot() +
    geom_stars(data = stars::st_as_stars(diff_bathy)) +
    geom_sf(data = map_layers$akland) +
    coord_sf(xlim = map_layers$plot.boundary$x, 
             ylim = map_layers$plot.boundary$y) +
    scale_fill_fermenter(name = expression(Z[new]-Z[2022]), 
                         palette = "BrBG", 
                         breaks = c(-Inf, -1000, -100, -10, -1, 1, 10, 100, 1000, Inf),
                         na.value = NA) +
    scale_x_continuous(breaks = map_layers$lon.breaks) +
    scale_y_continuous(breaks = map_layers$lat.breaks) +
    theme_2()
  
  y_dim_mult <- diff(map_layers$plot.boundary$y)/diff(map_layers$plot.boundary$x)
  
  ragg::agg_png(filename = here::here("analysis", "efh_bathymetry_2028", "plots", "bathy", paste0("bathy_diff_2022_current_d_", region[ii], ".png")), 
                width = 6.5, 
                height = y_dim_mult*6.5, 
                units = "in", 
                res = 180)
  print(plot_disc_bathy_diff)
  dev.off()
  
  ragg::agg_png(filename = here::here("analysis", "efh_bathymetry_2028", "plots", "bathy", paste0("bathy_diff_2022_current_c_", region[ii], ".png")), 
                width = 6.5, 
                height = y_dim_mult*6.5, 
                units = "in", 
                res = 180)
  print(plot_cont_bathy_diff)
  dev.off()
  
  
  # Plot observed tow depth versus raster tow depth
  
  depth_bins <- c(0, 50, seq(100, 1000, 100), 1200, 1400, Inf)
  
  towmid <- sf::st_read(here::here("analysis", "efh_bathymetry_2028", "data", "bts_tows", paste0(region[ii], "_towmid.shp"))) |>
    sf::st_transform(crs = map_layers$crs) |>
    dplyr::mutate(DEPTH_BIN = cut(DEPTH, 
                                  breaks = depth_bins,
                                  labels = paste(depth_bins[1:length(depth_bins)-1], 
                                                 depth_bins[2:length(depth_bins)], 
                                                 sep = "-")))
  
  old_bathy_points <- terra::as.points(old_bathy_reproj) |>
    sf::st_as_sf()
  
  new_bathy_points  <- terra::as.points(new_bathy) |>
    sf::st_as_sf()
  
  towmid$OLD_BATHY_DEPTH <- old_bathy_points[[1]][sf::st_nearest_feature(towmid, old_bathy_points)]
  towmid$NEW_BATHY_DEPTH <- new_bathy_points[[1]][sf::st_nearest_feature(towmid, new_bathy_points)]  
  
  xy_range <- range(c(towmid$OLD_BATHY_DEPTH, towmid$NEW_BATHY_DEPTH, towmid$DEPTH), na.rm = TRUE)

  plot_depth_scatter <- cowplot::plot_grid(
    ggplot() +
      geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2) + 
      geom_point(data = towmid, 
                 mapping = aes(x = OLD_BATHY_DEPTH, y = DEPTH), color = "grey50", alpha = 0.3) +
      scale_x_continuous(name = "2022 EFH Depth (m)", 
                         limits = xy_range) +
      scale_y_continuous(name = "Observed tow midpoint depth (m)\n2000-2023", 
                         limits = xy_range) +
      theme_bw(),
    ggplot() +
      geom_abline(intercept = 0, slope = 1, color = "black", linetype = 2) +
      geom_point(data = towmid, 
                 mapping = aes(x = NEW_BATHY_DEPTH, y = DEPTH), color = "grey50", alpha = 0.3) +
      scale_x_continuous(name = "New Compilation Depth (m)", 
                         limits = xy_range) +
      scale_y_continuous(name = "Observed tow midpoint depth (m)\n2000-2023", 
                         limits = xy_range) +
      theme_bw()
  )
  
  
  ragg::agg_png(filename = here::here("analysis", "efh_bathymetry_2028", "plots", "bathy", paste0("bathy_raster_vs_observed_", region[ii], ".png")), 
                width = 6.5, 
                height = 3.25, 
                units = "in", 
                res = 300)
  print(plot_depth_scatter)
  dev.off()
  
  bias_df <- towmid |>
    sf::st_drop_geometry() |>
    dplyr::group_by(DEPTH_BIN) |>
    dplyr::summarise(NEW_MRE = mean_relative_error(DEPTH, NEW_BATHY_DEPTH),
                     NEW_MAE = mean_absolute_error(DEPTH, NEW_BATHY_DEPTH),
                     NEW_BIAS = mean_bias(DEPTH, NEW_BATHY_DEPTH),
                     OLD_MRE = mean_relative_error(DEPTH, OLD_BATHY_DEPTH),
                     OLD_MAE = mean_absolute_error(DEPTH, OLD_BATHY_DEPTH),
                     OLD_BIAS = mean_bias(DEPTH, OLD_BATHY_DEPTH)) |>
    dplyr::mutate(REGION = toupper(region[ii]))
  
  write.csv(bias_df, 
            file = here::here("analysis", "efh_bathymetry_2028", "plots", "bathy", paste0("bathy_raster_accuracy_", region[ii], ".csv")),
            row.names = FALSE)
  
}


