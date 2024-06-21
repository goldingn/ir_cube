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

# load time-varying net coverage data
nets_cube <- rast("data/clean/nets_per_capita_cube.tif")

# load the other layers
covs_flat <- rast("data/clean/flat_covariates.tif")




# make prediction rasters
# (split this into a separate function, to be run on multiple scenarios, in
# parallel, without refitting)

# years to predict to
years_predict <- 2000:2030
n_times_predict <- length(years_predict)

# cells to predict to
cells_predict <- terra::cells(mask)
n_cells_predict <- length(cells_predict)

# pad out the nets cube, repeating the final year into the future
n_future <- n_times_predict - n_times
nets_cube_future <- nets_cube[[n_times]] %>%
  replicate(n_future,
            .,
            simplify = FALSE) %>%
  do.call(c, .)
years_future <- baseline_year + n_times + seq_len(n_future) - 1
names(nets_cube_future) <- paste0("nets_", years_future)
nets_cube_predict <- c(nets_cube, nets_cube_future)

# pull out temporally-static covariates for all cells
flat_extract_predict <- covs_flat %>%
  extract(cells_predict) %>%
  mutate(
    cell = cells_predict,
    .before = everything()
  )

# extract spatiotemporal covariates from the cube
all_extract_predict <- nets_cube_predict %>%
  terra::extract(cells_predict) %>%
  mutate(
    cell = cells_predict,
    .before = everything()
  ) %>%
  pivot_longer(
    cols = starts_with("nets_"),
    names_prefix = "nets_",
    names_to = "year",
    values_to = "net_coverage"
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

# insecticide types to save
types_save <- c("Deltamethrin",
                "Permethrin",
                "Alpha-cypermethrin")


# loop through these batches of cells/years and insecticides, and save a file of
# the batch predictions to disk

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
    
    # the insecticide to compute for
    insecticide,
    
    # an RNG seed to make sure all batches use the same posterior samples of
    # parameters
    rng_seed,
    
    # path to the file containing all the objects from model fitting
    fitted_model_image_file = "temporary/fitted_model.RData"
) {
  
  # load all the objects from the fitted model image here
  load(fitted_model_image_file)
  
  # index to this insecticide
  type_id_predict <- match(insecticide, types)
  
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
  selection_cell_years_batch <- x_cell_years_batch %*% effect_type[, type_id_predict]
  fitness_cell_years_batch <- 1 + selection_cell_years_batch
  
  # reformat this in to a 3D array with dimensions:
  #   n_times x n_unique_cells x (n_types = 1) x 1
  # to solve dynamics with time-varying fitness (time must be first, then other
  # two must match state variable, which has a trailing dimension of size 1)
  fitness_array_batch <- fitness_cell_years_batch
  dim(fitness_array_batch) <- c(n_times_predict, batch_n, 1, 1)
  
  # # check this is the right orientation!
  # which_cell_id <- 3
  # array_subset <- fitness_array_batch[, which_cell_id, 1, ]
  # cell_years_index_batch <- cell_years_predict_index[x_rows_batch, ]
  # cell_idx <- cell_years_index_batch$cell_id == which_cell_id
  # matrix_subset <- fitness_cell_years_batch[cell_idx, 1]
  # sims <- calculate(
  #   array_subset,
  #   matrix_subset,
  #   nsim = 1)
  # identical(sims$array_subset[1, , 1, 1, 1],
  #           sims$matrix_subset[1, , 1])
  
  # expand initial conditions out to all cells with data (plus a trailing
  # dimension to match greta.dynamics interface)
  init_array_batch <- init_fraction_susceptible[type_id_predict] * ones(batch_n)
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
  
  # get posterior draws of these, fixing the RNG seed so it's the same
  # collection of posterior samples for all batches
  pred_batch <- calculate(population_mortality_batch,
                          values = draws,
                          seed = rng_seed,
                          trace_batch_size = 25,
                          nsim = 50)[[1]]
  
  # compute posterior mean over the draws, for cells and years
  batch_pred_mean <- apply(pred_batch, 2:4, mean)
  
  # write this to a tibble (cell, year, insecticide, value), and save as a csv.
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
              insecticide, batch_index),
      row.names = FALSE
    )
  
  invisible()
  
}

# set up the future plan, to run each prediction batch in a new process. This is
# to limit the overall memory use on the machine because of memory growth if all
# run in a single process, not to parallelise computation. So set the number of
# workers to do that.
plan(multisession,
     workers = 1)

# loop through insecticides
for (this_insecticide in types_save) {
  
  # # grab the current seed, so we can do calculate in batches but not shuffle
  # # parameters over space
  # this_seed <- greta::.internals$utils$misc$get_seed()
  seed <- get_seed()
    
    # run predictions in new processes
    future_lapply(seq_len(n_batches),
                  predict_batch,
                  # raster cell numbers to loop over
                  cell_batches = cell_batches,
                  # index to the cells
                  cell_id_batches = cell_id_batches,
                  # combinations of all cells and years
                  cell_years_predict_index = cell_years_predict_index,
                  # data extracted for all cells and years
                  x_cell_years_predict = x_cell_years_predict,
                  # the insecticide to compute for
                  insecticide = this_insecticide,
                  # an RNG seed to make sure all batches use the same posterior samples of
                  # parameters
                  rng_seed = seed,
                  # this probably isn't needed since we fix the seed at the
                  # level of TF, but it stops a warning from future
                  future.seed = TRUE)

    
  # }

}

# write code to reassemble these and save rasters to disk

for (this_insecticide in types_save) {
  
  # load all the prediction files for this insecticide
  prediction_files <- list.files("temporary/prediction_files",
                                 pattern = this_insecticide,
                                 full.names = TRUE)
  
  for (this_year in years_predict) {
    
    print(year)
    
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

# copy files over to Tas
this_insecticide <- "Deltamethrin"
source_string <- sprintf("outputs/ir_maps/%s",
                        this_insecticide)
tifs <- list.files(source_string, pattern = "*.tif", full.names = TRUE)
new_dest <- "/mnt/Z/gfatm_scenarios/data/IR_rasters"
new_dir <- file.path(new_dest, this_insecticide)
dir.create(new_dir, showWarnings = FALSE)  
lapply(tifs,
       file.copy,
       to = new_dir)
