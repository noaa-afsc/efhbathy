library(akgfmaps)
library(MultiscaleDTM)
source(here::here("R", "bbox_to_polygon.R")) # Function to be included in an upcoming akgfmaps release

ai_layers <- akgfmaps::get_base_layers(select.region = "ai", set.crs = "EPSG:3338")

# 100 m resolution raster
efh_100m <- terra::rast(here::here("analysis", "efh_bathymetry_2028", "output", "efh_bathy", "efh_bathy_100m.tif"))

ai_towpath <- sf::st_read(here::here("analysis", "efh_bathymetry_2028", "data", "bts_tows", "ai_towpath.shp")) |>
  sf::st_transform(crs = ai_layers$crs$input)

ai_towpath <- ai_towpath[which(!is.na(sf::st_is_valid(ai_towpath))), ]

nazan_bay <- data.frame(x = -173.9335093, y = 52.1756401, label = "Nazan Bay") |>
  sf::st_as_sf(coords = c("x", "y"), crs = "WGS84") |>
  sf::st_transform(crs = ai_layers$crs$input)

plot_bbox <- sf::st_buffer(nazan_bay, dist = 15000) |>
  sf::st_bbox()

nazan_bbox <- sf::st_buffer(nazan_bay, dist = 15000) |>
  bbox_to_polygon()

expanded_bbox <- sf::st_buffer(nazan_bay, dist = 100000) |>
  bbox_to_polygon()


nazan_bathy <- terra::crop(efh_100m, sf::st_as_sf(nazan_bbox))
expanded_bathy <- terra::crop(efh_100m, sf::st_as_sf(expanded_bbox))

bpi_nazan_2023_method <- MultiscaleDTM::BPI(r = nazan_bathy, w = c(3, 64), unit = "cell", na.rm = TRUE)
bpi_nazan_2023_method <- round(bpi_nazan_2023_method + 0.5)

bpi_expanded_2023_method <- MultiscaleDTM::BPI(r = expanded_bathy, w = c(3, 64), unit = "cell", na.rm = TRUE)
bpi_expanded_2023_method <- round(bpi_expanded_2023_method + 0.5)

plot_nazan_bathy <- ggplot() +
  geom_stars(data = stars::st_as_stars(nazan_bathy)) +
  geom_sf(data = ai_layers$akland) +
  geom_sf(data = ai_towpath,
          mapping = aes(color = factor(sign(PERFORM), levels = c(-1, 0, 1), labels = c("Bad", "Good", "Acceptable")))) +
  scale_color_manual(name = "Performance", values = c("Bad" = "red", "Acceptable" = "purple4", "Good" = "lightgreen")) +
  scale_fill_distiller(na.value = NA, direction = 1) +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())


bpi_range <- range(terra::global(c(bpi_nazan_2023_method, bpi_expanded_2023_method, bpi_expanded_15_64), range, na.rm = TRUE))


plot_nazan_bpi_2023_method <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_nazan_2023_method)) +
  geom_sf(data = ai_layers$akland) +
  geom_sf(data = ai_towpath,
          mapping = aes(color = factor(sign(PERFORM), levels = c(-1, 0, 1), labels = c("Bad", "Good", "Acceptable")))) +
  scale_color_manual(name = "Performance", values = c("Bad" = "red", "Acceptable" = "purple4", "Good" = "lightgreen")) +
  scale_fill_distiller(name = "Focal BPI", na.value = NA, direction = 1, palette = "BrBG", limits = bpi_range) +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())


plot_expanded_bpi_2023_method <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_expanded_2023_method)) +
  geom_sf(data = ai_layers$akland) +
  geom_sf(data = ai_towpath,
          mapping = aes(color = factor(sign(PERFORM), levels = c(-1, 0, 1), labels = c("Bad", "Good", "Acceptable")))) +
  scale_color_manual(name = "Performance", values = c("Bad" = "red", "Acceptable" = "purple4", "Good" = "lightgreen")) +
  scale_fill_distiller(name = "Extended BPI", na.value = NA, direction = 1, palette = "BrBG", limits = bpi_range) +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

plot_expanded_bpi_2023_diff <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_expanded_2023_method)-stars::st_as_stars(bpi_nazan_2023_method)) +
  geom_sf(data = ai_layers$akland) +
  geom_sf(data = ai_towpath,
          mapping = aes(color = factor(sign(PERFORM), levels = c(-1, 0, 1), labels = c("Bad", "Good", "Acceptable")))) +
  scale_color_manual(name = "Performance", values = c("Bad" = "red", "Acceptable" = "purple4", "Good" = "lightgreen")) +
  scale_fill_distiller(name = "Extended-Focal", na.value = NA, direction = 1, palette = "PuOr") +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

cowplot::plot_grid(plot_nazan_bathy, 
                   plot_nazan_bpi_2023_method, 
                   plot_expanded_bpi_2023_method, 
                   plot_expanded_bpi_2023_diff, 
                   nrow = 2)

# Comparing BPI at different scales
# Keep in mind that the typical tow paths is ~1.5 km

bpi_expanded_15_64 <- MultiscaleDTM::BPI(r = expanded_bathy, w = c(15, 64), unit = "cell", na.rm = TRUE)

bpi_expanded_15_30 <- MultiscaleDTM::BPI(r = expanded_bathy, w = c(15, 30), unit = "cell", na.rm = TRUE)

bpi_expanded_3_10 <- MultiscaleDTM::BPI(r = expanded_bathy, w = c(3, 10), unit = "cell", na.rm = TRUE)

bpi_expanded_3_15 <- MultiscaleDTM::BPI(r = expanded_bathy, w = c(3, 15), unit = "cell", na.rm = TRUE)

plot_expanded_bpi_15_64 <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_expanded_15_64)) +
  geom_sf(data = ai_layers$akland) +
  geom_sf(data = ai_towpath,
          mapping = aes(color = factor(sign(PERFORM), levels = c(-1, 0, 1), labels = c("Bad", "Good", "Acceptable")))) +
  scale_color_manual(name = "Performance", values = c("Bad" = "red", "Acceptable" = "purple4", "Good" = "lightgreen")) +
  scale_fill_distiller(name = "BPI (15, 64)", na.value = NA, direction = 1, palette = "BrBG", limits = bpi_range) +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

plot_expanded_bpi_15_30 <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_expanded_15_30)) +
  geom_sf(data = ai_layers$akland) +
  geom_sf(data = ai_towpath,
          mapping = aes(color = factor(sign(PERFORM), levels = c(-1, 0, 1), labels = c("Bad", "Good", "Acceptable")))) +
  scale_color_manual(name = "Performance", values = c("Bad" = "red", "Acceptable" = "purple4", "Good" = "lightgreen")) +
  scale_fill_distiller(name = "BPI (15, 30)", na.value = NA, direction = 1, palette = "BrBG", limits = bpi_range) +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

plot_expanded_bpi_3_10 <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_expanded_3_10)) +
  geom_sf(data = ai_layers$akland) +
  geom_sf(data = ai_towpath,
          mapping = aes(color = factor(sign(PERFORM), levels = c(-1, 0, 1), labels = c("Bad", "Good", "Acceptable")))) +
  scale_color_manual(name = "Performance", values = c("Bad" = "red", "Acceptable" = "purple4", "Good" = "lightgreen")) +
  scale_fill_distiller(name = "BPI (3, 10)", na.value = NA, direction = 1, palette = "BrBG", limits = bpi_range) +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

plot_expanded_bpi_3_15 <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_expanded_3_15)) +
  geom_sf(data = ai_layers$akland) +
  geom_sf(data = ai_towpath,
          mapping = aes(color = factor(sign(PERFORM), levels = c(-1, 0, 1), labels = c("Bad", "Good", "Acceptable")))) +
  scale_color_manual(name = "Performance", values = c("Bad" = "red", "Acceptable" = "purple4", "Good" = "lightgreen")) +
  scale_fill_distiller(name = "BPI (3, 15)", na.value = NA, direction = 1, palette = "BrBG", limits = bpi_range) +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())


cowplot::plot_grid(plot_expanded_bpi_15_64, 
                   plot_expanded_bpi_15_30, 
                   plot_expanded_bpi_3_15, 
                   plot_expanded_bpi_15_30, 
                   nrow = 2)


# Aggregated to 1 km 

# Comparing BPI at different scales
# Keep in mind that the typical tow paths is ~1.5 km

bpi_expanded_15_64_1km <- MultiscaleDTM::BPI(r = expanded_bathy, w = c(15, 64), unit = "cell", na.rm = TRUE) |>
  terra::aggregate(fact = 10, fun = "mean")

bpi_expanded_15_30_1km <- MultiscaleDTM::BPI(r = expanded_bathy, w = c(15, 30), unit = "cell", na.rm = TRUE) |>
  terra::aggregate(fact = 10, fun = "mean")

bpi_expanded_3_10_1km <- MultiscaleDTM::BPI(r = expanded_bathy, w = c(3, 10), unit = "cell", na.rm = TRUE) |>
  terra::aggregate(fact = 10, fun = "mean")

bpi_expanded_3_15_1km <- MultiscaleDTM::BPI(r = expanded_bathy, w = c(3, 15), unit = "cell", na.rm = TRUE) |>
  terra::aggregate(fact = 10, fun = "mean")

plot_expanded_bpi_15_64_1km <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_expanded_15_64_1km)) +
  geom_sf(data = ai_layers$akland) +
  geom_sf(data = ai_towpath,
          mapping = aes(color = factor(sign(PERFORM), levels = c(-1, 0, 1), labels = c("Bad", "Good", "Acceptable")))) +
  scale_color_manual(name = "Performance", values = c("Bad" = "red", "Acceptable" = "purple4", "Good" = "lightgreen")) +
  scale_fill_distiller(name = "BPI (15, 64)", na.value = NA, direction = 1, palette = "BrBG", limits = bpi_range) +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

plot_expanded_bpi_15_30_1km <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_expanded_15_30_1km)) +
  geom_sf(data = ai_layers$akland) +
  geom_sf(data = ai_towpath,
          mapping = aes(color = factor(sign(PERFORM), levels = c(-1, 0, 1), labels = c("Bad", "Good", "Acceptable")))) +
  scale_color_manual(name = "Performance", values = c("Bad" = "red", "Acceptable" = "purple4", "Good" = "lightgreen")) +
  scale_fill_distiller(name = "BPI (15, 30)", na.value = NA, direction = 1, palette = "BrBG", limits = bpi_range) +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

plot_expanded_bpi_3_10_1km <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_expanded_3_10_1km)) +
  geom_sf(data = ai_layers$akland) +
  geom_sf(data = ai_towpath,
          mapping = aes(color = factor(sign(PERFORM), levels = c(-1, 0, 1), labels = c("Bad", "Good", "Acceptable")))) +
  scale_color_manual(name = "Performance", values = c("Bad" = "red", "Acceptable" = "purple4", "Good" = "lightgreen")) +
  scale_fill_distiller(name = "BPI (3, 10)", na.value = NA, direction = 1, palette = "BrBG", limits = bpi_range) +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())

plot_expanded_bpi_3_15_1km <- ggplot() +
  geom_stars(data = stars::st_as_stars(bpi_expanded_3_15_1km)) +
  geom_sf(data = ai_layers$akland) +
  geom_sf(data = ai_towpath,
          mapping = aes(color = factor(sign(PERFORM), levels = c(-1, 0, 1), labels = c("Bad", "Good", "Acceptable")))) +
  scale_color_manual(name = "Performance", values = c("Bad" = "red", "Acceptable" = "purple4", "Good" = "lightgreen")) +
  scale_fill_distiller(name = "BPI (3, 15)", na.value = NA, direction = 1, palette = "BrBG", limits = bpi_range) +
  coord_sf(xlim = c(plot_bbox['xmin'], plot_bbox['xmax']),
           ylim = c(plot_bbox['ymin'], plot_bbox['ymax']),
           expand = 0) +
  theme_bw() +
  theme(axis.title = element_blank())


cowplot::plot_grid(plot_expanded_bpi_15_64_1km, 
                   plot_expanded_bpi_15_30_1km, 
                   plot_expanded_bpi_3_15_1km, 
                   plot_expanded_bpi_15_30_1km, 
                   nrow = 2)
