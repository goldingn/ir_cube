# The projection shared by the latent smooths of the dynamical model
# (R/latent_smooth.R) and the two-stage correction (R/two_stage_correction.R).
# Functions only; needs sf and terra.

# Albers equal-area conic for Africa, in km, so that the Matern range means
# the same distance everywhere and the range priors can be stated in km
africa_equal_area_crs <- paste(
  "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25",
  "+x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"
)

project_km <- function(lon, lat, crs = africa_equal_area_crs) {
  xy <- sf::sf_project(from = "+proj=longlat +datum=WGS84 +no_defs", to = crs,
                       pts = cbind(lon, lat))
  colnames(xy) <- c("x_km", "y_km")
  xy
}

# Projected centres of the prediction mask's cells, on a coarse grid (every
# `fact`-th cell), for the outer boundary of the meshes. Callers pass them to
# build_correction_meshes()
prediction_mask_coords <- function(file = "data/clean/raster_mask.tif",
                                   fact = 10) {
  mask <- terra::aggregate(terra::rast(file), fact = fact, fun = "min",
                           na.rm = TRUE)
  xy <- terra::xyFromCell(mask, terra::cells(mask))
  project_km(xy[, 1], xy[, 2])
}
