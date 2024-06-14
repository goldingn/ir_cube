# prepare rasters

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# ITN cube
itn_files <- list.files("data/raw/nets_per_capita/",
                        pattern = ".tif",
                        full.names = TRUE)
itn_years <- itn_files %>%
  basename() %>%
  str_remove("ITN_") %>%
  str_remove("_percapita_nets_mean.tif")
nets_per_capita <- rast(itn_files)
names(nets_per_capita) <- paste0("nets_", itn_years)

terra::writeRaster(nets_per_capita,
                   "data/clean/nets_per_capita_cube.tif",
                   overwrite = TRUE)

# other static ones

# NB I made a symlink from Gerry's dropbox into my into data raw with the
# following mac terminal command (from inside data/raw):
# ln -s ~/Dropbox/sharing/ir_pred_rasters ./

# load crops, standardised to 0-1
crops_std <- rast("data/raw/ir_pred_rasters/crop_2010_std.tif")

# subset to main staple crops (more likely to be used in subsistence farming and
# therefore influence rural malaria vectors): cereals, roots, and bananas/plantains. Exclude 'other' categories, and any with comparative low coverage
crop_type_keep <- c(
  # cereals
  # "wheat",  # too sparse 
  "rice",
  "maize",
  # "barley",  # too sparse 
  "pearl millet",
  # "small millet",  # too sparse 
  "sorghum",
  "potato", # roots
  "sweet potato",
  "yams",
  "cassava",
  "banana",
  "plantain"
)

crop_codes <- geodata::spamCrops() %>%
  as_tibble() %>%
  mutate(
    code = toupper(code)
  ) %>%
  filter(
    crop %in% crop_type_keep
  ) %>%
  pull(code)

crops_std_sub <- crops_std[[crop_codes]]

# do a pca on these to 

# add an epsilon and logit transform these layers
logit_eps <- function(x, epsilon = 1e-6) {
  scale <- 1 - 2 * epsilon
  z <- epsilon + x * scale
  # remove 0s and 1s
  qlogis(z)
}

logit_crops <- app(crops_std_sub, logit_eps)
scale_logit_crops <- scale(logit_crops)


# now PCA
pca <- terra::princomp(scale_logit_crops, maxcell = 1e6)
pca$loadings
crop_pc <- predict(scale_logit_crops, pca, index = 1:5)

# rescale these for modelling
crop_pc_scale <- scale(crop_pc)

# give them names
names(crop_pc_scale) <- paste0("crop_pc_",
                               seq_len(nlyr(crop_pc_scale)))

# # plot them
# ggplot() +
#   geom_spatraster(
#     data = crop_pc_scale
#   ) +
#   scale_fill_viridis_c(
#     na.value = "transparent"
#   ) +
#   facet_wrap(~lyr) +
#   theme_minimal()


# combine these into a multiband geotiff of flat (not temporally static)
# covariates

flat_covs <- c(crop_pc_scale)

# save these out
terra::writeRaster(flat_covs,
                   "data/clean/flat_covariates.tif",
                   overwrite = TRUE)