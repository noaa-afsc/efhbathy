library(akgfmaps)
library(MultiscaleDTM)
source(here::here("R", "bbox_to_polygon.R")) # Function to be included in an upcoming akgfmaps release

# zone = "Nazan Bay"
zone = "RKC Savings Area"

locations <- read.csv(file = here::here("analysis", "efh_bathymetry_2028", "data", "terrain_example_zones.csv"))

set_loc <- dplyr::filter(locations, 
                         label == zone)

map_layers <- akgfmaps::get_base_layers(select.region = set_loc$region, set.crs = "EPSG:3338")

# Load 100 m resolution raster and tow paths ----
efh_100m <- terra::rast(here::here("analysis", "efh_bathymetry_2028", "output", "efh_bathy", "efh_bathy_100m.tif"))

region_towpath <- sf::st_read(here::here("analysis", "efh_bathymetry_2028", "data", "bts_tows", paste0(set_loc$region, "_towpath.shp"))) |>
  sf::st_transform(crs = map_layers$crs$input) |>
  dplyr::filter(CRUISE > 199000) |>
  dplyr::mutate(PERFORMANCE_LABEL = factor(sign(PERFORM), levels = c(-1, 0, 1), labels = c("Bad", "Good", "Acceptable")))

region_towpath <- region_towpath[which(!is.na(sf::st_is_valid(region_towpath))), ]

performance_color <- c("Bad" = "red", "Acceptable" = "purple4", "Good" = "blue")


# Setup plot extent and bounding box for focal zone ----

focal_zone <- set_loc |>
  sf::st_as_sf(coords = c("x", "y"), crs = "WGS84") |>
  sf::st_transform(crs = map_layers$crs$input)

plot_bbox <- sf::st_buffer(focal_zone, dist = set_loc$focal_buffer) |>
  sf::st_bbox()

focal_bbox <- sf::st_buffer(focal_zone, dist = set_loc$focal_buffer) |>
  bbox_to_polygon()

expanded_bbox <- sf::st_buffer(focal_zone, dist = set_loc$extended_buffer) |>
  bbox_to_polygon()

expanded_bathy <- terra::crop(efh_100m, 
                              expanded_bbox)

focal_bathy <- terra::mask(expanded_bathy, 
                           focal_bbox)

plot(focal_bathy)

# Calculate BPI following the 2023 method from Lundblad et al. (2006) ----

bpi_focal_2023_method <- MultiscaleDTM::BPI(r = focal_bathy, 
                                            w = c(3, 64), 
                                            unit = "cell", 
                                            na.rm = TRUE) |>
  terra::crop(focal_bathy)

bpi_focal_2023_method <- round(bpi_focal_2023_method + 0.5)

bpi_expanded_2023_method <- MultiscaleDTM::BPI(r = expanded_bathy, 
                                               w = c(3, 64), 
                                               unit = "cell", 
                                               na.rm = TRUE) |>
  terra::crop(focal_bathy)

bpi_expanded_2023_method <- round(bpi_expanded_2023_method + 0.5)

plot_focal_bathy <- ggplot() +
  geom_stars(data = stars::st_as_stars(focal_bathy)) +
  geom_sf(data = map_layers$akland, fill = "grey30", color = NA) +
  geom_sf(data = region_towpath,
          mapping = aes(color = PERFORMANCE_LABEL)) +
  scale_color_manual(name = "Performance", values = performance_color) +
  scale_fill_distiller(name = "Depth (m)", na.value = "#BED2FF", direction = 1) +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']) + c(-3000, 3000),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']) + c(-3000, 3000),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

plot_focal_bpi_2023_method <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_focal_2023_method)) +
  geom_sf(data = map_layers$akland, fill = "grey30", color = NA) +
  geom_sf(data = region_towpath,
          mapping = aes(color = PERFORMANCE_LABEL)) +
  scale_color_manual(name = "Performance", values = performance_color) +
  scale_fill_distiller(name = "Focal BPI", 
                       na.value = "#BED2FF", 
                       direction = 1, 
                       # limits = bpi_range, 
                       palette = "BrBG") +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']) + c(-3000, 3000),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']) + c(-3000, 3000),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())


plot_expanded_bpi_2023_method <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_expanded_2023_method)) +
  geom_sf(data = map_layers$akland, fill = "grey30", color = NA) +
  geom_sf(data = region_towpath,
          mapping = aes(color = PERFORMANCE_LABEL)) +
  scale_color_manual(name = "Performance", values = performance_color) +
  scale_fill_distiller(name = "Extended BPI", 
                       na.value = "#BED2FF", 
                       direction = 1, 
                       # limits = bpi_range, 
                       palette = "BrBG") +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']) + c(-3000, 3000),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']) + c(-3000, 3000),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

plot_expanded_bpi_2023_diff <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_expanded_2023_method)-stars::st_as_stars(bpi_focal_2023_method)) +
  geom_sf(data = map_layers$akland, fill = "grey30", color = NA) +
  geom_sf(data = region_towpath,
          mapping = aes(color = PERFORMANCE_LABEL)) +
  scale_color_manual(name = "Performance", values = performance_color) +
  scale_fill_distiller(name = "Extended-Focal", 
                       na.value = "#BED2FF", 
                       direction = 1, 
                       palette = "PuOr") +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']) + c(-3000, 3000),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']) + c(-3000, 3000),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

cowplot::plot_grid(plot_focal_bathy, 
                   plot_focal_bpi_2023_method, 
                   plot_expanded_bpi_2023_method, 
                   plot_expanded_bpi_2023_diff,
                   nrow = 2)

# Comparing BPI with different annulus size ----
# Keep in mind that the typical tow paths is ~1.5 km

bpi_15_64 <- MultiscaleDTM::BPI(r = expanded_bathy, 
                                w = c(15, 64), 
                                unit = "cell", 
                                na.rm = TRUE) |>
  terra::crop(focal_bbox)

bpi_15_30 <- MultiscaleDTM::BPI(r = expanded_bathy, 
                                w = c(15, 30), 
                                unit = "cell", 
                                na.rm = TRUE) |>
  terra::crop(focal_bbox)

bpi_3_15 <- MultiscaleDTM::BPI(r = expanded_bathy, 
                               w = c(3, 15), 
                               unit = "cell", 
                               na.rm = TRUE) |>
  terra::crop(focal_bbox)

bpi_3_10 <- MultiscaleDTM::BPI(r = expanded_bathy, 
                               w = c(3, 10), 
                               unit = "cell", 
                               na.rm = TRUE) |>
  terra::crop(focal_bbox)

plot_bpi_15_64 <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_15_64)) +
  geom_sf(data = map_layers$akland, fill = "grey30", color = NA) +
  geom_sf(data = region_towpath,
          mapping = aes(color = PERFORMANCE_LABEL)) +
  scale_color_manual(name = "Performance", values = performance_color) +
  scale_fill_distiller(name = "BPI (15, 64)", 
                       na.value = "#BED2FF", 
                       direction = 1, 
                       # limits = bpi_range, 
                       palette = "BrBG") +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

plot_bpi_15_30 <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_15_30)) +
  geom_sf(data = map_layers$akland, fill = "grey30", color = NA) +
  geom_sf(data = region_towpath,
          mapping = aes(color = PERFORMANCE_LABEL)) +
  scale_color_manual(name = "Performance", values = performance_color) +
  scale_fill_distiller(name = "BPI (15, 30)", 
                       na.value = "#BED2FF", 
                       direction = 1, 
                       # limits = bpi_range, 
                       palette = "BrBG") +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

plot_bpi_3_15 <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_3_15)) +
  geom_sf(data = map_layers$akland, fill = "grey30", color = NA) +
  geom_sf(data = region_towpath,
          mapping = aes(color = PERFORMANCE_LABEL)) +
  scale_color_manual(name = "Performance", values = performance_color) +
  scale_fill_distiller(name = "BPI (3, 15)", 
                       na.value = "#BED2FF", 
                       direction = 1, 
                       # limits = bpi_range, 
                       palette = "BrBG") +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

plot_bpi_3_10 <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_3_10)) +
  geom_sf(data = map_layers$akland, fill = "grey30", color = NA) +
  geom_sf(data = region_towpath,
          mapping = aes(color = factor(sign(PERFORM), 
                                       levels = c(-1, 0, 1), 
                                       labels = c("Bad", "Good", "Acceptable"))),
          linewidth = rel(1.1)) +
  scale_color_manual(name = "Performance", values = performance_color) +
  scale_fill_distiller(name = "BPI (3, 10)", 
                       na.value = "#BED2FF", 
                       direction = 1, 
                       # limits = bpi_range, 
                       palette = "BrBG") +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

cowplot::plot_grid(plot_bpi_15_64, 
                   plot_bpi_15_30, 
                   plot_bpi_3_15, 
                   plot_bpi_3_10, 
                   nrow = 2)


# Aggregated to 1 km 

# Comparing BPI at different scales
# Keep in mind that the typical tow paths is ~1.5 km

bpi_15_64_1km <- bpi_15_64 |>
  terra::aggregate(fact = 10, fun = "mean")

bpi_15_30_1km <- bpi_15_30 |>
  terra::aggregate(fact = 10, fun = "mean")

bpi_3_10_1km <- bpi_3_10 |>
  terra::aggregate(fact = 10, fun = "mean")

bpi_3_15_1km <- bpi_3_15 |>
  terra::aggregate(fact = 10, fun = "mean")

plot_bpi_15_64_1km <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_15_64_1km)) +
  geom_sf(data = map_layers$akland, fill = "grey30", color = NA) +
  geom_sf(data = region_towpath,
          mapping = aes(color = PERFORMANCE_LABEL)) +
  scale_color_manual(name = "Performance", values = performance_color) +
  scale_fill_distiller(name = "BPI (15, 64)", 
                       na.value = "#BED2FF", 
                       direction = 1, 
                       # limits = bpi_range, 
                       palette = "BrBG") +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

plot_bpi_15_30_1km <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_15_30_1km)) +
  geom_sf(data = map_layers$akland, fill = "grey30", color = NA) +
  geom_sf(data = region_towpath,
          mapping = aes(color = PERFORMANCE_LABEL)) +
  scale_color_manual(name = "Performance", values = performance_color) +
  scale_fill_distiller(name = "BPI (15, 30)", 
                       na.value = "#BED2FF", 
                       direction = 1, 
                       # limits = bpi_range, 
                       palette = "BrBG") +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

plot_bpi_3_10_1km <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_3_10_1km)) +
  geom_sf(data = map_layers$akland, fill = "grey30", color = NA) +
  geom_sf(data = region_towpath,
          mapping = aes(color = PERFORMANCE_LABEL)) +
  scale_color_manual(name = "Performance", values = performance_color) +
  scale_fill_distiller(name = "BPI (3, 10)", 
                       na.value = "#BED2FF", 
                       direction = 1, 
                       # limits = bpi_range, 
                       palette = "BrBG") +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

plot_bpi_3_15_1km <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_3_15_1km)) +
  geom_sf(data = map_layers$akland, fill = "grey30", color = NA) +
  geom_sf(data = region_towpath,
          mapping = aes(color = PERFORMANCE_LABEL)) +
  scale_color_manual(name = "Performance", values = performance_color) +
  scale_fill_distiller(name = "BPI (3, 15)", 
                       na.value = "#BED2FF", 
                       direction = 1, 
                       # limits = bpi_range, 
                       palette = "BrBG") +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

cowplot::plot_grid(plot_bpi_15_64_1km, 
                   plot_bpi_15_30_1km, 
                   plot_bpi_3_15_1km, 
                   plot_bpi_15_30_1km, 
                   nrow = 2)
