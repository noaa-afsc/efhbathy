library(akgfmaps)
library(MultiscaleDTM)

# Example of BPI at multiple scales using the 10 x 10 m resolution Maungawhau volcano raster

r <- terra::rast(volcano, 
                 extent= ext(2667400,
                             2667400 + ncol(volcano) * 10,
                             6478700,
                             6478700 + nrow(volcano)*10),
                 crs = "EPSG:27200")

terra::res(r)

bpi_2_4 <- MultiscaleDTM::BPI(r, w = c(2, 4), stand= "none", unit = "cell", na.rm = TRUE)                                
bpi_4_10<- MultiscaleDTM::BPI(r, w = c(4, 10), stand= "none", unit = "cell", na.rm = TRUE)
bpi_10_20 <- MultiscaleDTM::BPI(r, w = c(10, 20), stand= "none", unit = "cell", na.rm = TRUE)
                                
par(mfrow = c(2,2))
plot(r, main = "Elevation")
plot(bpi_2_4, main = "BPI (2, 4)")
plot(bpi_4_10, main = "BPI (4, 10)")
plot(bpi_10_20, main = "BPI (10, 20)")

# Constructing BPI from scratch ----

# Construct annulus
inner_radius <- 2
outer_radius <- 4
resolution <- terra::res(r)

annulus <- MultiscaleDTM::annulus_window(radius = c(inner_radius, outer_radius), 
                                         unit = "cell", 
                                         resolution = resolution)

par(mfrow = c(1,1))
plot(terra::rast(annulus))

# Calculate focal mean
focal_mean <- terra::focal(x = r, w = annulus, fun = mean, na.rm = TRUE)

bpi <- r - focal_mean

# This is identical to the MultiscaleDTM output for bpi_2_4
plot(bpi - bpi_2_4)
all(values(bpi - bpi_2_4) == 0)

# Lundblad et al. (2006) integer BPI and scale factor scheme
bpi_int <- round(bpi + 0.5) 

scale_factor <- resolution * outer_radius
  
par(mfrow = c(2,2))
plot(r, main = "Elevation")
plot(focal_mean, main = "Focal Mean")
plot(bpi_int, main = paste0("Continuous BPI(", scale_factor, ")"))
plot(bpi_int, main = paste0("Discrete BPI(", scale_factor, ")"))




par(mfrow = c(2,2))
plot(r, main = "Elevation")
plot(round(bpi_2_4 + 0.5) , main = "BPI(40)")
plot(round(bpi_4_10 + 0.5) , main = "BPI(100)")
plot(round(bpi_10_20 + 0.5) , main = "BPI(200)")
