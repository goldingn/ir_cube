# make predictions from the fitted greta model

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# set a global RNG seed for prediction
set.seed(1)

# load the fitted model objects here, to set up predictions
load(file = "temporary/fitted_model.RData")

# load covariate rasters (need to reload these even if restoring workspace
# because pointers)

# load the mask
mask <- rast("data/clean/raster_mask.tif")

# load time-varying net use, IRS coverage, population data
nets_cube <- rast("data/clean/net_use_cube.tif")
irs_cube <- rast("data/clean/irs_coverage_scaled_cube.tif")
pop_cube <- rast("data/clean/pop_scaled_cube.tif")

# for each one pad back to the baseline year, repeating the first
nets_cube <- pre_pad_cube(nets_cube, baseline_year)
irs_cube <- pre_pad_cube(irs_cube, baseline_year)
pop_cube <- pre_pad_cube(pop_cube, baseline_year)

# load the non-temporal crop covariate layers

# collated total yields of crop types
crops_group <- rast("data/clean/crop_group_scaled.tif")

# yields of individual crops
crops_all <- rast("data/clean/crop_scaled.tif")

# Pull out crop types implicated in risk for IR. refer to this review, crop type
# section:
# https://malariajournal.biomedcentral.com/articles/10.1186/s12936-016-1162-4
crops_implicated <- c(
  # "increased resistance at cotton growing sites, a finding subsequently
  # supported in eight other papers from five different African countries", "the
  # cash crop with the highest intensity insecticide use of any crop"
  crops_all$cotton,
  # "In eight studies, vegetable cultivation strongly related to
  # insecticide-resistant field collections", "Vegetable production requires
  # significantly higher quantities and/or more frequent application of
  # pesticides than other food crops"
  crops_all$vegetables,
  # "Seven of the studies reviewed here examined the insecticide susceptibility
  # of vector populations at rice-growing sites, and found low-to-moderate
  # resistance levels in these mosquito populations."
  crops_all$rice)

# combine all temporally-static covariates
covs_flat <- c(crops_group, crops_implicated)

# load the country raster for initial condition lookup
country_raster <- rast("data/clean/country_raster.tif")

# make prediction rasters
# (split this into a separate function, to be run on multiple scenarios, in
# parallel, without refitting)

# years to predict to
years_predict <- baseline_year:2030
n_times_predict <- length(years_predict)

# cells to predict to
cells_predict <- terra::cells(mask)
n_cells_predict <- length(cells_predict)

# pad out the nets cube, repeating the final year into the future
end_year <- max(years_predict)
nets_cube <- post_pad_cube(nets_cube, end_year)
irs_cube <- post_pad_cube(irs_cube, end_year)
pop_cube <- post_pad_cube(pop_cube, end_year)

# pull out temporally-static covariates for all cells
flat_extract_predict <- covs_flat %>%
  extract(cells_predict) %>%
  mutate(
    cell = cells_predict,
    .before = everything()
  )

# extract spatiotemporal covariates from the cube
all_extract_predict <- bind_cols(
  terra::extract(nets_cube, cells_predict),
  terra::extract(irs_cube, cells_predict),
  terra::extract(pop_cube, cells_predict)
) %>%
  mutate(
    cell = cells_predict,
    .before = everything()
  ) %>%
  # this stacks all the different cubes in long format, but we want wide on the
  # variable but long on year, so pivot_wider immediately after
  pivot_longer(
    cols = -one_of("cell"),
    names_sep = "_",
    names_to = c("variable", "year"),
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = "variable",
    values_from = "value"
  ) %>%
  mutate(
    year = as.numeric(year)
  ) %>%
  left_join(
    flat_extract_predict,
    by = "cell"
  ) %>%
  mutate(
    cell_id = match(cell, cells_predict),
    year_id = year - baseline_year + 1,
    .before = everything()
  ) %>%
  filter(
    year >= baseline_year
  ) %>%
  select(
    -cell,
    -year
  )

# pull out index to cells and years
cell_years_predict_index <- all_extract_predict %>%
  select(cell_id, year_id)

# get covariates for these cell-years as a matrix
x_cell_years_predict <- all_extract_predict %>%
  select(-cell_id,
         -year_id) %>%
  as.matrix()

# get the country for each cell
cell_country_predict <- extract(country_raster, cells_predict)$country_name

# create batches of cells for processing

# NOTE due to a weird-as-hell greta bug, this object must not be called
# batch_size: https://github.com/greta-dev/greta/issues/634
batch_bigness <- 1e4
batch_idx <- seq_along(cells_predict) %/% batch_bigness
# index to raster, for setting cell values
cell_batches <- split(cells_predict, batch_idx)
# index to cells vector, for getting covariate values
cell_id_batches <- split(seq_along(cells_predict), batch_idx)
n_batches <- length(cell_batches)

# which to predict - just effective LLIN susceptibility for now
types_save <- c(
  "Permethrin",
  "Alpha-cypermethrin",
  "Deltamethrin",
  "llin_effective",
  "DDT",
  "Pirimiphos-methyl", 
  "Fenitrothion",
  "Bendiocarb",
  "Malathion", 
  "Lambda-cyhalothrin"
)

# loop through these batches of cells/years and insecticide outputs, and save a
# file of the batch predictions to disk

predict_batch <- function(
    # number to iterate over
    batch_index,
    
    # vectors/dfs of data in all the batches to later subset the necessary parts
    
    # raster cell numbers to loop over
    cell_batches,
    
    # index to the cells
    cell_id_batches,
    
    # combinations of all cells and years
    cell_years_predict_index,
    
    # data extracted for all cells and years
    x_cell_years_predict,
    
    # country to which each cell belongs
    cell_country_predict,
    
    # the insecticide type to predict for
    insecticide_type,
    
    # an RNG seed to make sure all batches use the same posterior samples of
    # parameters
    rng_seed,
    
    # path to the file containing all the objects from model fitting
    fitted_model_image_file = "temporary/fitted_model.RData"
) {
  
  # load all the objects from the fitted model image here
  load(fitted_model_image_file)
  
  # reload the packages and functions
  source("R/packages.R")
  source("R/functions.R")
  
  # if the output is the llin_effective susceptibility, get IDs for the multiple
  # active ingredients used, and the weights to combine them
  if (insecticide_type == "llin_effective") {
    
    # load the proportional weights of different insecticides in single-AI LLINs
    ingredient_weights <- readRDS("temporary/ingredient_weights.RDS")
    
    # do predictions for these insecticides
    ingredient_ids <- match(names(ingredient_weights), types)
    
  } else {
    
    # otherwise, pull out the single ingredient ID we are modelling
    ingredient_ids <- match(insecticide_type, types)
    
  }

  # how many to process
  n_ingredients <- length(ingredient_ids)
  
  # find cells to write to    
  cell_batch <- cell_batches[[batch_index]]
  batch_n <- length(cell_batch)
  
  # pull out the index to rows of x_cell_years_predict that correspond to this
  # batch of cells
  cell_id_batch <- cell_id_batches[[batch_index]]
  
  # need to find the elements that match the cell_id, but return them in the
  # correct year order
  x_rows_batch <- which(cell_years_predict_index$cell_id %in% cell_id_batch)
  
  # compute selection coefficients for cell-years in this batch and convert to
  # relative fitness
  x_cell_years_batch <- x_cell_years_predict[x_rows_batch, ]
  selection_cell_years_batch <- x_cell_years_batch %*% effect_type[, ingredient_ids]
  fitness_cell_years_batch <- 1 + selection_cell_years_batch
  
  # reformat this in to a 3D array with dimensions:
  #   n_times x n_unique_cells x (n_types = 1) x 1
  # to solve dynamics with time-varying fitness (time must be first, then other
  # two must match state variable, which has a trailing dimension of size 1)
  fitness_array_batch <- fitness_cell_years_batch
  dim(fitness_array_batch) <- c(n_times_predict, batch_n, n_ingredients, 1)
  
  # # check this is the right orientation!
  # which_cell_id <- 3
  # which_type_id <- 1
  # array_subset <- fitness_array_batch[, which_cell_id, which_type_id, ]
  # cell_years_index_batch <- cell_years_predict_index[x_rows_batch, ]
  # cell_idx <- cell_years_index_batch$cell_id == which_cell_id
  # matrix_subset <- fitness_cell_years_batch[cell_idx, 1]
  # sims <- calculate(
  #   array_subset,
  #   matrix_subset,
  #   nsim = 1)
  # identical(sims$array_subset[1, , 1, 1, 1],
  #           sims$matrix_subset[1, , 1])
  
  # find the initial fraction susceptible for all cells and insecticides, using
  # lookup against country-level model, for all countries in the mask
  
  # for each cell in this batch, pull out the country, the region, and therefore
  # the index to the random variables for the initial fraction
  lookup <- country_region_lookup()

  all_countries <- unique(lookup$country_name)
  all_regions <- unique(lookup$region)
  n_all_countries <- length(all_countries)
  n_all_regions <- length(all_regions)
  
  unobserved_countries <- setdiff(all_countries, countries)
  n_countries_unobserved <- length(unobserved_countries)
  unobserved_regions <- setdiff(all_regions, regions)
  n_regions_unobserved <- length(unobserved_regions)
    
  # get the new vectors of region (rearrange) and country (pad and rearrange)
  # parameters, and recombine into initial conditions here
  countries_observed_index <- match(countries, all_countries)
  regions_observed_idx <- match(regions, all_regions)
  
  # create empty matrices for the initial coniditoion parameters
  init_region_raw_pred <- zeros(n_all_regions, n_types)
  init_country_raw_pred <- zeros(n_all_countries, n_types)

  # add in the estimated parameters for observed countries
  init_country_raw_pred[countries_observed_index, ] <- init_country_raw
  init_region_raw_pred[regions_observed_idx, ] <- init_region_raw
  
  # add any the new country or region parameters for countries/regions not in
  # the data, extrapolating with hierarchical model
  if (n_countries_unobserved > 0) {
    countries_unobserved_index <- match(unobserved_countries, all_countries)
    init_country_raw_unobserved <- normal(0, 1, dim = c(n_countries_unobserved, n_types))
    init_country_raw_pred[countries_unobserved_index, ] <- init_country_raw_unobserved
  }
  
  if (n_regions_unobserved > 0) {
    regions_unobserved_index <- match(unobserved_regions, all_regions)
    init_region_raw_unobserved <- normal(0, 1, dim = c(n_regions_unobserved, n_types))
    init_region_raw_pred[regions_unobserved_index, ] <- init_region_raw_unobserved
  }

  # combine to get the regional and country-level deviation from the prior logit
  # mean
  init_region_effect_pred <- sweep(init_region_raw_pred, 2, init_region_sd, FUN = "*")
  init_country_effect_pred <- sweep(init_country_raw_pred, 2, init_country_sd, FUN = "*")
  
  # get a lookup from countries to regions, in the new order
  # for each country in all_countries (ensuring to use the new order), find the regions
  regions_for_all_countries <- lookup$region[match(all_countries, lookup$country_name)]
  # pull out the index of these countries to the full set of regions, in the new order
  all_country_region_index <- match(regions_for_all_countries, all_regions)
  
  # for each country in the full dataset, combine the national and the regional
  # effects and get the logit initial fraction susceptibile for that country
  init_country_overall_effect_pred <- init_country_effect_pred +
    init_region_effect_pred[all_country_region_index, ]
  
  logit_init_country_pred <- sweep(init_country_overall_effect_pred,
                              2,
                              logit_init_mean,
                              FUN = "+")
  
  # convert from relative (0-1) to the constrained scale (above init_frac_min)
  init_country_relative_pred <- ilogit(logit_init_country_pred)
  init_country_magnitude_pred <- sweep(init_country_relative_pred, 2, init_range, FUN = "*") 
  init_country_pred <- sweep(init_country_magnitude_pred, 2, init_frac_min, FUN = "+")
  
  # expand out to all prediction cells
  batch_countries <- cell_country_predict[cell_id_batch]
  batch_country_id <- match(batch_countries, all_countries)
  init_array_batch <- init_country_pred[batch_country_id, ingredient_ids]
  
  # add a trailing dimension to match greta.dynamics interface
  dim(init_array_batch) <- c(dim(init_array_batch), 1)
  
  # iterate through time to get fraction susceptible for all years at all
  # prediction cells
  dynamic_cells_batch <- iterate_dynamic_function(
    transition_function = haploid_next,
    initial_state = init_array_batch,
    niter = n_times_predict,
    w = fitness_array_batch,
    parameter_is_time_varying = c("w"),
    tol = 0)
  
  fraction_susceptible_batch <- dynamic_cells_batch$all_states
  population_mortality_batch <- fraction_susceptible_batch
  
  # if we want the effective susceptibility against LLINs, compute a weighted
  # sum of the active ingredients
  if (insecticide_type == "llin_effective") {
  
    # compute a weighted sum of these mortalities to get effective susceptibility
    # to LLIN insecticides
    effective_susc_batch <- zeros(batch_n, 1, n_times_predict)
    for (i in seq_along(ingredient_weights)) {
      ingredient_susc <- population_mortality_batch[, i, ]
      effective_susc_batch <- effective_susc_batch + ingredient_susc * ingredient_weights[[i]]
    }
    
    # overwrite this with the weighted sum
    population_mortality_batch <- effective_susc_batch
  
  }  
  
  # get posterior draws of these, fixing the RNG seed so it's the same
  # collection of posterior samples for all batches
  pred_batch <- calculate(population_mortality_batch,
                          values = draws,
                          seed = rng_seed,
                          trace_batch_size = 25,
                          nsim = 500)[[1]]
  
  # compute posterior mean over the draws, for cells and years
  batch_pred_mean <- apply(pred_batch, 2:4, mean)
  
  # write this to a tibble (cell, year,value), and save as a csv.
  # then write another script to 
  batch_pred_mean_mat <- batch_pred_mean[, 1, ]
  colnames(batch_pred_mean_mat) <- years_predict
  batch_pred_mean_mat %>%
    as_tibble() %>%
    mutate(cell = cell_batch,
           .before = everything()) %>%
    pivot_longer(
      cols = -any_of(c("cell")),
      names_to = "year",
      values_to = "mean"
    ) %>%
    write.csv(
      sprintf("temporary/prediction_files/pred_mean_%s_batch_%i.csv",
              insecticide_type, batch_index),
      row.names = FALSE
    )
  
  invisible()
  
}

# How many parallel workers to use? n_workers = 1 is sequential mode, note that
# prediction is memory intensive, so pick the number fo workers to fit within
# memory when all running simultaneously
n_workers <- 4

# future workers are re-used within a single plan, which leads to a memory leak
# if workers are used more than once. So define some 'super batches', so that
# each worker only gets one job before it is reset.
# for picking up broken runs
# batches <- 93:n_batches
batches <- seq_len(n_batches)
super_batches <- split(batches, batches %/% n_workers)

# loop through insecticides
for (this_insecticide in types_save) {
  
  # grab the current seed, so we can do calculate in batches but not shuffle
  # parameters over space
  seed <- get_seed()
  
  # do a sub-batch
  for (super_batch in super_batches) {
    
    # set up the future plan, to run each prediction batch in a new process. This is
    # to limit the overall memory use on the machine because of memory growth if all
    # run in a single process, not (necessarily) to parallelise computation.
    
    # destroy any previous parallel workers by switching to sequential model,
    # and then recreate the parallel workers
    plan(sequential)
    plan(multisession,
         workers = n_workers)
    
    # # if the process crashes, any zombie processes can be killed off by doing
    # # resetWorkers(), using an identical call to the parallel call
    # future::resetWorkers(plan(multisession, workers = n_workers))
    
    # loop through these batches, running predictions in new processes
    future_lapply(X = super_batch,
                  FUN = predict_batch,
                  # raster cell numbers to loop over
                  cell_batches = cell_batches,
                  # index to the cells
                  cell_id_batches = cell_id_batches,
                  # combinations of all cells and years
                  cell_years_predict_index = cell_years_predict_index,
                  # data extracted for all cells and years
                  x_cell_years_predict = x_cell_years_predict,
                  # vector of countries to which each cell (in the full
                  # prediction set) belongs
                  cell_country_predict = cell_country_predict,
                  # the insecticide to compute for
                  insecticide_type = this_insecticide,
                  # an RNG seed to make sure all batches use the same posterior samples of
                  # parameters
                  rng_seed = seed,
                  # this probably isn't needed since we fix the seed at the
                  # level of TF, but it stops a warning from future
                  future.seed = TRUE)
    
  }
    
}

for (this_insecticide in types_save) {

  # write code to reassemble these and save rasters to disk
  # load all the prediction files for this insecticide
  this_pattern <- sprintf("pred_mean_%s_batch_",
                          this_insecticide)
  prediction_files <- list.files("temporary/prediction_files",
                                 pattern = this_pattern,
                                 full.names = TRUE)

  for (this_year in years_predict) {
    
    print(this_year)
    
    # loop through these, loading, subsetting to this year, and returning (to
    # minimise memory usage)
    load_year <- function(file, this_year) {
      read_csv(file,
               show_col_types = FALSE) %>%
        filter(year == this_year) %>%
        select(cell, mean)
    }
    
    prediction_data_this_year <- lapply(prediction_files,
                                        load_year,
                                        this_year)
    
    # combine into a single tibble of all the predictions for this insecticide
    # and year
    prediction_this_year <- do.call(bind_rows,
                                    prediction_data_this_year)
    
    
    # create a raster for this year
    this_ir_raster <- mask
    
    # insert all the values in it
    this_ir_raster[prediction_this_year$cell] <- prediction_this_year$mean
    
    # save to disk in the appropriate place
    write_path <- file.path("outputs/ir_maps",
                            this_insecticide,
                            sprintf("ir_%s_susceptibility.tif",
                                    this_year))
    
    dir.create(dirname(write_path))
    writeRaster(this_ir_raster,
                write_path,
                overwrite = TRUE)
    
  }
  
}

# # copy files over to Tas
# this_insecticide <- "llin_effective"
# source_string <- sprintf("outputs/ir_maps/%s",
#                          this_insecticide)
# tifs <- list.files(source_string, pattern = "*.tif", full.names = TRUE)
# new_dest <- "/mnt/efs/transition/gfatm_scenarios/data/IR_rasters"
# new_dir <- file.path(new_dest, this_insecticide)
# dir.create(new_dir, showWarnings = FALSE)  
# lapply(tifs,
#        file.copy,
#        to = new_dir)
