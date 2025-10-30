# Load all the prediction layers from Hancock et al. 2020 for 2006-2017, and
# convert to a consistently named terra spatraster object
# https://doi.org/10.6084/m9.figshare.9912623

# Note: there is also a 2005 layer listed, but Figshare was unable to serve it
# for download

library(raster)
library(terra)

# where the files are, unzipped
hancock_dir <- "data/raw/hancock_2020"

files <- list.files(hancock_dir, pattern = "*.grd")

# extract years as numeric
years <- files %>%
  str_remove("IR rasters ") %>%
  str_remove(".grd") %>%
  as.numeric()

# load in each year (each is a 5 layer raster for 5 insecticides)
year_bricks <- files %>%
  file.path(hancock_dir, .) %>%
  lapply(brick)

# scrape the end off each layer name, replace with the year (for indexing) and
# rename each one
for (i in seq_along(year_bricks)) {
  
  new_layer_names <- year_bricks[[i]] %>%
    names() %>%
    str_remove("_mortality") %>%
    sprintf("%s_%s", ., years[i])
  
  names(year_bricks[[i]]) <- new_layer_names
  
}

# collapse all the layers and convert to a terra raster object
hancock_preds <- year_bricks %>%
  do.call(stack, .) %>%
  terra::rast()

# save to disk
terra::writeRaster(hancock_preds,
                   "data/clean/hancock_2020_predictions.tif")
