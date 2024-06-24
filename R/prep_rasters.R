# prepare rasters

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# ITN cubes

# nets per capita
itn_pc_files <- list.files("data/raw/itn/nets_per_capita/",
                        pattern = ".tif",
                        full.names = TRUE)
itn_pc_years <- itn_pc_files %>%
  basename() %>%
  str_remove("ITN_") %>%
  str_remove("_percapita_nets_mean.tif")
nets_per_capita <- rast(itn_pc_files)
names(nets_per_capita) <- paste0("nets_", itn_pc_years)

terra::writeRaster(nets_per_capita,
                   "data/clean/nets_per_capita_cube.tif",
                   overwrite = TRUE)

# and scale it to 0-1
nets_per_capita_max <- max(global(nets_per_capita,
                                  "max",
                                  na.rm = TRUE))
nets_per_capita_scaled <- nets_per_capita / nets_per_capita_max
terra::writeRaster(nets_per_capita_scaled,
                   "data/clean/nets_per_capita_scaled_cube.tif",
                   overwrite = TRUE)

# net use
itn_use_files <- list.files("data/raw/itn/net_use/",
                        pattern = ".tif",
                        full.names = TRUE)
itn_use_years <- itn_use_files %>%
  basename() %>%
  str_remove("ITN_") %>%
  str_remove("_use_mean.tif")
net_use <- rast(itn_use_files)
names(net_use) <- paste0("nets_", itn_use_years)

terra::writeRaster(net_use,
                   "data/clean/net_use_cube.tif",
                   overwrite = TRUE)

# and scale it to 0-1
net_use_max <- max(global(net_use,
                          "max",
                          na.rm = TRUE))
net_use_scaled <- net_use / net_use_max
terra::writeRaster(net_use_scaled,
                   "data/clean/net_use_scaled_cube.tif",
                   overwrite = TRUE)

# make a mask for later use
mask <- nets_per_capita[[1]] * 0
names(mask) <- "mask" 
terra::writeRaster(mask,
                   "data/clean/raster_mask.tif",
                   overwrite = TRUE)

# IRS cube
# these are admin-level IRS coverage copied from
# GBD2023/Processing/Stages/14_IRS_Coverage_Africa/Intermediary_Outputs/admin_cov/20231107
irs_files <- list.files("data/raw/irs/admin_coverage",
                        pattern = ".tif",
                        full.names = TRUE)
irs_years <- irs_files %>%
  basename() %>%
  str_remove("admin_coverage_") %>%
  str_remove(".tif")
irs_coverage <- rast(irs_files)

# extend these to the ITN mask
irs_coverage <- terra::extend(irs_coverage, mask)

# pad them with zeros (where distribution data missing)
irs_coverage[is.na(irs_coverage)] <- 0

# mask them
irs_coverage <- mask(irs_coverage, mask)

# set their names
names(irs_coverage) <- paste0("irs_", irs_years)

terra::writeRaster(irs_coverage,
                   "data/clean/irs_coverage_cube.tif",
                   overwrite = TRUE)

# scale this
irs_coverage_max <- max(global(irs_coverage,
                               "max",
                               na.rm = TRUE))
irs_coverage_scaled <- irs_coverage / irs_coverage_max
terra::writeRaster(irs_coverage_scaled,
                   "data/clean/irs_coverage_scaled_cube.tif",
                   overwrite = TRUE)

# other static ones

# NB I made a symlink from Gerry's dropbox into my into data raw with the
# following mac terminal command (from inside data/raw):
# ln -s ~/Dropbox/sharing/ir_pred_rasters ./

# load crops, both as-is and standardised to 0-1
crops <- rast("data/raw/ir_pred_rasters/crop_2010.tif")

# # subset to main staple crops (more likely to be used in subsistence farming and
# # therefore influence rural malaria vectors): cereals, roots, and bananas/plantains. Exclude 'other' categories, and any with comparative low coverage
# crop_type_staple <- c(
#   # cereals
#   # "wheat",  # too sparse 
#   "rice",
#   "maize",
#   # "barley",  # too sparse 
#   "pearl millet",
#   # "small millet",  # too sparse 
#   "sorghum",
#   "potato", # roots
#   "sweet potato",
#   "yams",
#   "cassava",
#   "banana",
#   "plantain"
# )
# 
# # subset to staples
# crop_codes_staples <- crop_info %>%
#   filter(
#     crop %in% crop_type_staple
#   ) %>%
#   pull(code)
# 
# crops_staples <- crops[[crop_codes_staples]]
# crops_std_staples <- crops_std[[crop_codes_staples]]

cereals <- c("wheat", "rice", "maize", "barley", "pearl millet", "small millet",
             "sorghum", "other cereals")

roots <- c("potato", "sweet potato", "yams", "cassava", "other roots")

pulses <- c("bean", "chickpea",  "cowpea", "pigeonpea", "lentil",
            "other pulses")

oilcrops <- c("soybean", "groundnut", "coconut", "oilpalm", "sunflower",
              "rapeseed", "sesameseed", "other oil crops")

fibrecrops <- c("sugarcane", "sugarbeet", "cotton", "other fibre crops")

othercrops <- c("arabica coffee", "robusta coffee", "cocoa", "tea", "tobacco",
                "banana", "plantain", "tropical fruit", "temperate fruit",
                "vegetables", "rest of crops")

crop_info <- geodata::spamCrops() %>%
  as_tibble() %>%
  mutate(
    code = toupper(code),
    type = case_when(
      crop %in% cereals ~ "cereal",
      crop %in% roots ~ "root",
      crop %in% pulses ~ "pulse",
      crop %in% oilcrops ~ "oil crop",
      crop %in% fibrecrops ~ "fibre crop",
      crop %in% othercrops ~ "other"
    )
  )

# compute and standardise the yield of all crops
crops_group <- app(crops, sum)
crops_group <- crops_group / global(crops_group, "max", na.rm = TRUE)[1, 1]
names(crops_group) <- "all"

# do the same for groups of these and append them
for (crop_type in unique(crop_info$type)) {
  
  crop_codes_type <- crop_info %>%
    filter(
      type == crop_type
    ) %>%
    pull(code)
  
  crops_type <- app(crops[[crop_codes_type]], sum)
  crops_type <- crops_type / global(crops_type, "max", na.rm = TRUE)[1, 1]
  names(crops_type) <- crop_type
  crops_group <- c(crops_group, crops_type)
  
}

names(crops_group) <- paste0("crops_", names(crops_group))

# 
# 
# 
# 
# 
# # do a pca on these to 
# 
# # add an epsilon and logit transform these layers
# logit_eps <- function(x, epsilon = 1e-6) {
#   scale <- 1 - 2 * epsilon
#   z <- epsilon + x * scale
#   # remove 0s and 1s
#   qlogis(z)
# }
# 
# logit_crops <- app(crops_std, logit_eps)
# scale_logit_crops <- scale(logit_crops)
# 
# 
# # now PCA
# pca <- terra::princomp(scale_logit_crops, maxcell = 1e6)
# crop_pc <- predict(scale_logit_crops, pca, index = 1:5)
# 
# # rescale these for modelling
# crop_pc_scale <- scale(crop_pc)
# 
# # give them names
# names(crop_pc_scale) <- paste0("crop_pc_",
#                                seq_len(nlyr(crop_pc_scale)))
# 
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
# ggsave("figures/cov_crop_pcs.png",
#        bg = "white",
#        width = 8,
#        height = 6,
#        dpi = 300)

# combine these into a multiband geotiff of flat (not temporally static)
# covariates

flat_covs <- c(crops_group)

# save these out
terra::writeRaster(flat_covs,
                   "data/clean/flat_covariates.tif",
                   overwrite = TRUE)