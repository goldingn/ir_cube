# prepare rasters
files <- list.files("data/raw/nets_per_capita/",
                    pattern = ".tif",
                    full.names = TRUE)
years <- files %>%
  basename() %>%
  str_remove("ITN_") %>%
  str_remove("_percapita_nets_mean.tif")
nets_per_capita <- rast(files)
names(nets_per_capita) <- paste0("nets_", years)
plot(nets_per_capita)

terra::writeRaster(nets_per_capita, "data/clean/nets_per_capita_cube.tif")
