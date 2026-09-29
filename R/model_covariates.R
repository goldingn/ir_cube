# Covariates of the dynamical model that are built the same way for the fits,
# the folds and the maps. Functions only; needs terra. Source from the repo
# root.


# initial-state covariates (#19) -----------------------------------------------

# Static predictors of the logit initial fraction susceptible, which let the
# initial state vary within a country. Each is standardised over every cell of
# the mask, so the values at a cell do not depend on which cells have data, and
# the fits, folds and maps share them:
#   log_pop_2000  log population per cell in 2000, the earliest layer of the
#                 population cube (there is none before 2000)
#   all_crops     agricultural productivity: the "all crops" yield layer of
#                 crop_group_scaled.tif, as log(x + 1e-4). The layer is 0-1
#                 scaled and very skewed (56% zeros, mean 0.003, sd 0.012), so
#                 standardised untransformed it reaches 47 sd at the data
#                 cells; 1e-4 is about the 10th percentile of its non-zero
#                 values, and on this scale the data cells span -0.7 to 4.1 sd
init_covariate_names <- c("log_pop_2000", "all_crops")

init_covariate_layers <- function() {
  mask <- rast("data/clean/raster_mask.tif")
  log_pop <- log(rast("data/clean/pop_cube.tif")[["pop_2000"]])
  crops <- log(rast("data/clean/crop_group_scaled.tif")[["all crops"]] + 1e-4)
  layers <- terra::mask(c(log_pop, crops), mask)
  names(layers) <- init_covariate_names
  moments <- terra::global(layers, c("mean", "sd"), na.rm = TRUE)
  (layers - moments$mean) / moments$sd
}

# The standardised initial-state covariates at mask cells `cells`, as a
# cells x covariates matrix with named columns.
init_covariate_matrix <- function(cells, layers = init_covariate_layers()) {
  x <- as.matrix(terra::extract(layers, cells))
  colnames(x) <- names(layers)
  x
}
