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

# pad them with zeros (where distribution data missing, including -9999 coding)
irs_coverage[irs_coverage < 0] <- 0
irs_coverage[is.na(irs_coverage)] <- 0

# mask them
irs_coverage <- mask(irs_coverage, mask)

# set their names
names(irs_coverage) <- paste0("irs_", irs_years)

# clip to the analysis years (not before 2000)
irs_coverage <- irs_coverage[[paste0("irs_", 2000:2022)]]
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

# Human population, combining MAP's WorldPop observed and projected global
# rasters to get 2000-2022

# past population, 2000-2020
pop_past_files <- list.files("data/raw/pop/",
                             pattern = ".tif",
                             full.names = TRUE)
pop_past_years <- pop_past_files %>%
  basename() %>%
  str_remove("WorldPop_UNAdj_v3_DRC_fix.") %>%
  str_remove(".Annual.Data.5km.sum.tif")
pop_past <- rast(pop_past_files)
names(pop_past) <- paste0("pop_", pop_past_years)

# projected future population, 2021-2049
pop_future_files <- list.files("data/raw/pop_future/",
                               pattern = ".tif",
                               full.names = TRUE)
pop_future_years <- pop_future_files %>%
  basename() %>%
  str_remove("WorldPop_UNAdj_v3_DRC_fix_projected.") %>%
  str_remove(".Annual.Data.5km.sum.tif")
pop_future <- rast(pop_future_files)
names(pop_future) <- paste0("pop_", pop_future_years)

# combine, crop, and mask to mastergrids
pop_all <- c(pop_past, pop_future)
pop_all <- crop(pop_all, mask)
pop_all <- mask(pop_all, mask)

# fill in zeros on land with an epsilon, fill in NAs (outside mastergrids) with
# 0, and relevel total populations, just in case I end up using these for population
# stats by mistake.

# Set the epsilon to 1 person per 5km, approximately the 11th percentile of
# non-zero populations, and for Egypt (with most missing population) this
# matches the neighbouring low-population Sahara desert areas well

# annual populations for Africa within the mastergrids file, for relevelling
target_africa_pops <- global(pop_all, "sum", na.rm = TRUE)

# check approximate percentile of the epsilon
epsilon <- 1
epsilon_quantile <- mean(na.omit(c(pop_all[])) <= epsilon)
round(100 * epsilon_quantile)

# fill in the missing populations (and remask the true NA areas)
pop_all[is.na(pop_all)] <- 0
pop_all[pop_all == 0] <- epsilon
pop_all <- mask(pop_all, mask)

# relevel the total continent populations
new_africa_pops <- global(pop_all, "sum", na.rm = TRUE)
ratios <- target_africa_pops / new_africa_pops
for (i in seq_len(terra::nlyr(pop_all))) {
  pop_all[[i]] <- pop_all[[i]] * ratios[[1]][i]
}

# check there's less than a person different in the total populations
final_africa_pops <- global(pop_all, "sum", na.rm = TRUE)
max(abs(final_africa_pops - target_africa_pops)) < 1

# transform and scale to 0-1
# trans_pop_all <- log(pop_all)
trans_pop_all <- identity(pop_all)
trans_pop_all_min <- min(global(trans_pop_all,
                            "min",
                            na.rm = TRUE))
trans_pop_all <- trans_pop_all - trans_pop_all_min
trans_pop_all_max <- max(global(trans_pop_all,
                            "max",
                            na.rm = TRUE))
trans_pop_all <- trans_pop_all / trans_pop_all_max

# clip to modelling years, and save
pop <- pop_all[[paste0("pop_", 2000:2022)]]
trans_pop <- trans_pop_all[[paste0("pop_", 2000:2022)]]

terra::writeRaster(pop,
                   "data/clean/pop_cube.tif",
                   overwrite = TRUE)

terra::writeRaster(trans_pop,
                   "data/clean/pop_scaled_cube.tif",
                   overwrite = TRUE)

# save future projection years, too
pop_future <- pop_all[[paste0("pop_", 2023:2030)]]
trans_pop_future <- trans_pop_all[[paste0("pop_", 2023:2030)]]

terra::writeRaster(pop_future,
                   "data/clean/pop_cube_future.tif",
                   overwrite = TRUE)

terra::writeRaster(trans_pop_future,
                   "data/clean/pop_scaled_cube_future.tif",
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

# # drop 'all' before saving, since it's not useful as a covariate after all
# crops_group <- crops_group[[-1]]

# combine these into a multiband geotiff of flat (not temporally static)
# covariates

flat_covs <- c(crops_group)

# save these out
terra::writeRaster(flat_covs,
                   "data/clean/flat_covariates.tif",
                   overwrite = TRUE)
