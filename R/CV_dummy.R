

# build a filter flag for in Hancock et al or not to compare with our model
optimal_nn_preds <- optimal_nn_preds %>%
  # pre-emptively subset to Hancock prediction years and insecticides
  mutate(in_hancock =
           (year_start %in% 2006:2017 &
              insecticide_type %in% c("Alpha-cypermethrin",
                                      "DDT",
                                      "Deltamethrin",
                                      "Lambda-cyhalothrin",
                                      "Permethrin") )
  ) 

# build covariate rasters for the proper model

nets_cube <- pre_pad_cube(nets_cube, baseline_year)
irs_cube <- pre_pad_cube(irs_cube, baseline_year)
pop_cube <- pre_pad_cube(pop_cube, baseline_year)


crops_group <- rast("data/clean/crop_group_scaled.tif")
crops_all <- rast("data/clean/crop_scaled.tif")
crops_implicated <- c(
  crops_all$cotton,
  crops_all$vegetables,
  crops_all$rice)
covs_flat <- c(crops_group, crops_implicated)


# index to the classes for each type
classes_index <- df %>%
  distinct(type_id, class_id) %>%
  arrange(type_id) %>%
  pull(class_id)

# pull out concentrations for different types
type_concentrations <- df %>%
  select(type_id, concentration) %>%
  group_by(type_id) %>%
  filter(row_number() == 1) %>%
  arrange(type_id) %>%
  pull(concentration)


# create design matrix at all unique cells and for all years

# pull out temporally-static covariates for all cells
flat_extract <- covs_flat %>%
  extract(unique_cells) %>%
  mutate(
    cell = unique_cells,
    .before = everything()
  )

# extract spatiotemporal covariates from the cube
all_extract <- bind_cols(
  terra::extract(nets_cube, unique_cells),
  terra::extract(irs_cube, unique_cells),
  terra::extract(pop_cube, unique_cells)
) %>%
  mutate(
    cell = unique_cells,
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
    flat_extract,
    by = "cell"
  ) %>%
  mutate(
    cell_id = match(cell, unique_cells),
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
cell_years_index <- all_extract %>%
  select(cell_id, year_id)

# get covariates for these cell-years as a matrix
x_cell_years <- all_extract %>%
  select(-cell_id,
         -year_id) %>%
  as.matrix()

# dimensions of things in the fitting stage
n_covs <- ncol(x_cell_years)
n_obs <- nrow(df)
n_unique_cells <- length(unique_cells)
n_times <- max(df$year_start) - min(df$year_start) + 1
n_classes <- length(classes)
n_types <- length(types)
n_regions <- length(regions)
n_countries <- length(countries)


# greta arrays being defined below:

# this is the full model with no modification to definition (i.e., the model
# dynamically iterate through the whole area * time), the only difference from
# main model is simply that likelihood is defined over the training fold data.
# The model doesn't need to be redefined in each fold in principle, but we do so
# (and delete all greta arrays at the end of the fold loop) to avoid restarts
# and leaks of the HMC sampler

# utility function to force rebuild the model at the end of the loop
purge_greta_model <- function(){
  rm(beta_overall,
     sigma_overall,
     sigma_class,
     beta_class_raw,
     beta_class_sigma,
     beta_class,
     beta_type_raw,
     beta_type_sigma,
     beta_type,
     effect_type,
     selection_cell_years,
     fitness_cell_years,
     fitness_array,
     init_frac_relative_prior,
     logit_init_mean,
     init_region_sd,
     init_country_sd,
     init_region_raw,
     init_country_raw,
     init_region_effect,
     init_country_effect,
     init_country_overall_effect,
     logit_init_country,
     init_country_relative,
     init_range,
     init_country_magnitude,
     init_country,
     init_array,
     dynamic_cells,
     fraction_susceptible_vec,
     population_mortality_vec,
     rho_classes,
     m, 
     draws,
     envir = globalenv()
  )
}



# dynamic iteration function - no need to redefine each time
haploid_next <- function(state, iter, w) {
  # fraction susceptible
  q <- state
  # fraction resistant
  p <- 1 - q
  # number resistant grows relative to resistant by ratio w, re-normalise to get
  # fraction remaining susceptible
  q / (q + p * w)
}

#####do spatial extrapolation 

# initialise result df
dynamic_extrap_result <- optimal_nn_preds %>%
  filter(
    experiment == "spatial_extrapolation"
  ) %>% 
  group_by(country_name) %>% 
  summarise(
    pred_error_nn = betabinom_dev(died = died,
                                  mosquito_number = mosquito_number,
                                  predicted = predicted),
    bias_nn = mean(predicted - observed)
  ) %>% 
  mutate(pred_error_dynamic = NA,,
         bias_dynamic = NA,
         experiment = "dynmaic_spatial_extrapolation") %>% 
 left_join(spatial_extrapolation_intercept_error,by = join_by(country_name)) 

# initialise result df for vs hancock
dynamic_extrap_result_vs_hancock <- dynamic_extrap_result %>%
  select(country_name,pred_error_dynamic,bias_dynamic) 

for (country_idx in seq_along(countries_to_validate)) {
  
  
this_country <- countries_to_validate[country_idx]
# hierarchical regression parameters

# doubly hierarchical version

# between-classes
beta_overall <- normal(0, 1, dim = n_covs)
sigma_overall <- normal(0, 1, dim = n_covs, truncation = c(0, Inf))
sigma_class <- normal(0, 1, dim = n_covs, truncation = c(0, Inf))

# between-types
beta_class_raw <- normal(0, 1, dim = c(n_covs, n_classes))
beta_class_sigma <- sweep(beta_class_raw, 1, sigma_overall, FUN = "*")
beta_class <- sweep(beta_class_sigma, 1, beta_overall, FUN = "+")

# betas for types
beta_type_raw <- normal(0, 1, dim = c(n_covs, n_types))
beta_type_sigma <- sweep(beta_type_raw, 1, sigma_class, FUN = "*")
beta_type <- beta_class[, classes_index] + beta_type_sigma

# multiply through to get relative fitness of resistance for each insecticide
# type at each cell

effect_type <- exp(beta_type)


selection_cell_years <- x_cell_years %*% effect_type
fitness_cell_years <- 1 + selection_cell_years

# reformat this in to a 3D array with dimensions:
#   n_times x n_unique_cells x n_types x 1
# to solve dynamics with time-varying fitness (time must be first, then other
# two must match state variable, which has a trailing dimension of size 1)
fitness_array <- fitness_cell_years
dim(fitness_array) <- c(n_times, n_unique_cells, n_types, 1)

# prior and minimum values for the initial fractions susceptible
init_frac_prior <- ifelse(types == "DDT", 0.9, 0.95)
init_frac_min <- ifelse(types == "DDT", 0.75, 0.9)

# mean logit proportion of the distance from the minimum to 1, for each type
init_frac_relative_prior <- (init_frac_prior - init_frac_min) / (1 - init_frac_min)
logit_init_mean <- normal(qlogis(init_frac_relative_prior), 1, dim = n_types)

# variance at each level, for each insecticide
init_region_sd <- normal(0, 1, truncation = c(0, Inf), dim = n_types)
init_country_sd <- normal(0, 1, truncation = c(0, Inf), dim = n_types)

# unscaled deviation parameters (hierarchical decentring)
init_region_raw <- normal(0, 1, dim = c(n_regions, dim = n_types))
init_country_raw <- normal(0, 1, dim = c(n_countries, dim = n_types))

# combine to get the regional and country-level deviation from the prior logit
# mean
init_region_effect <- sweep(init_region_raw, 2, init_region_sd, FUN = "*")
init_country_effect <- sweep(init_country_raw, 2, init_country_sd, FUN = "*")

country_region_index <- df %>%
  group_by(country_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(country_id, region_id) %>%
  arrange(country_id) %>%
  pull(region_id)

# combine these together into the country-level logit-mean initial value
# logit_init_mean <- qlogis(init_frac_mean)
init_country_overall_effect <- init_country_effect + init_region_effect[country_region_index, ]
logit_init_country <- sweep(init_country_overall_effect,
                            2,
                            logit_init_mean,
                            FUN = "+")

# convert from relative (0-1) to the constrained scale (above init_frac_min)
init_country_relative <- ilogit(logit_init_country)
init_range <- 1 - init_frac_min
init_country_magnitude <- sweep(init_country_relative, 2, init_range, FUN = "*") 
init_country <- sweep(init_country_magnitude, 2, init_frac_min, FUN = "+")

cell_country_lookup <- df %>%
  group_by(cell_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(cell_id, country_id) %>%
  arrange(cell_id) %>%
  pull(country_id)

# expand out to all observations
init_array <- init_country[cell_country_lookup, ]

# add a trailing dimension to match greta.dynamics interface
dim(init_array) <- c(dim(init_array), 1)

# iterate through time to get fraction without resistant allele for all years at all cells
# with data
dynamic_cells <- iterate_dynamic_function(
  transition_function = haploid_next,
  initial_state = init_array,
  niter = n_times,
  w = fitness_array,
  parameter_is_time_varying = c("w"),
  tol = 0)

# pull out the values at the cells, insecticides, and years corresponding to the
# data

# pull out this training data
fold_train_df <- spatial_extrapolation$training[[country_idx]]

index <- cbind(fold_train_df$cell_id, fold_train_df$type_id, fold_train_df$year_id)
fraction_susceptible_vec <- dynamic_cells$all_states[index]

# resistant phenotype expression rate = expected mortality in bioassays
population_mortality_vec <- fraction_susceptible_vec

# define observation model
rho_classes <- normal(0, 0.025,
                      truncation = c(0, 1),
                      dim = n_classes)

distribution(fold_train_df$died) <- betabinomial_p_rho(
  N = fold_train_df$mosquito_number,
  p = population_mortality_vec,
  rho = rho_classes[fold_train_df$class_id])

m <- model(
  beta_overall,
  sigma_overall,
  beta_type_raw,
  rho_classes
)

n_chains <- 4

system.time(
  draws <- mcmc(m,
                chains = n_chains,
                n_samples = 1e3,
                warmup = 1e3)
)


fold_test_df <- spatial_extrapolation$test[[country_idx]]

index_test <- cbind(fold_test_df$cell_id, fold_test_df$type_id, fold_test_df$year_id)
fraction_susceptible_vec_test <- dynamic_cells$all_states[index_test]
population_mortality_vec_test <- fraction_susceptible_vec_test

predicted_test <- calculate(population_mortality_vec_test,
                            values = draws,
                            nsim = 1e3)

rho_classes_test <- calculate(rho_classes,
                              values = draws,
                              nsim = 1e3)

predicted_test_mean <- apply(predicted_test[[1]],2:3,mean) %>% as.numeric()

# add predicted value back to test data for bookkeeping
spatial_extrapolation$test[[country_idx]] <- tibble(spatial_extrapolation$test[[country_idx]], predicted_dynamic = predicted_test_mean)

# calculate posterior rho
rho_classes_test_mean <- apply(rho_classes_test[[1]],2:3,mean) %>% as.numeric()

rho_classes_test_indexed <- rho_classes_test_mean[fold_test_df$class_id]

# put back into the test df
spatial_extrapolation$test[[country_idx]] <- tibble(spatial_extrapolation$test[[country_idx]], rho_classes_test_indexed = rho_classes_test_indexed)


# calculate deviance and bias on test fold results
test_outcome_df <- optimal_nn_preds %>%
  filter(
    experiment == "spatial_extrapolation",
    country_name == this_country
  ) %>% 
  select(died,mosquito_number,observed,in_hancock)

pred_error_dynamic <- betabinom_dev(died = test_outcome_df$died,
                                    mosquito_number = test_outcome_df$mosquito_number,
                                    predicted = predicted_test_mean,
                                    rho = rho_classes_test_indexed)

bias_dynamic  <- mean(predicted_test_mean - test_outcome_df$observed)

# iteratively insert the results into the result df
dynamic_extrap_result[dynamic_extrap_result$country_name == this_country,'pred_error_dynamic'] <- pred_error_dynamic
dynamic_extrap_result[dynamic_extrap_result$country_name == this_country,'bias_dynamic'] <- bias_dynamic


# vs hancock

pred_error_dynamic_vs_hancock <- betabinom_dev(died = test_outcome_df$died[test_outcome_df$in_hancock],
                                               mosquito_number = test_outcome_df$mosquito_number[test_outcome_df$in_hancock],
                                               predicted = predicted_test_mean[test_outcome_df$in_hancock],
                                               rho = rho_classes_test_indexed[test_outcome_df$in_hancock])

bias_dynamic_vs_hancock  <- mean(predicted_test_mean[test_outcome_df$in_hancock] - test_outcome_df$observed[test_outcome_df$in_hancock])

# iteratively insert the results into the result df
dynamic_extrap_result_vs_hancock[dynamic_extrap_result_vs_hancock$country_name == this_country,'pred_error_dynamic'] <- pred_error_dynamic_vs_hancock
dynamic_extrap_result_vs_hancock[dynamic_extrap_result_vs_hancock$country_name == this_country,'bias_dynamic'] <- bias_dynamic_vs_hancock



purge_greta_model()
gc()
}

dynamic_extrap_result
write_csv(dynamic_extrap_result,"outputs/dynamic_extrap_result.csv")

dynamic_extrap_result_vs_hancock_joined <- hancock_extrap_result %>% 
  left_join(dynamic_extrap_result_vs_hancock,by = join_by(country_name)) %>% 
  mutate(experiment = "dynmaic_spatial_extrapolation_vs_hancock") 

dynamic_extrap_result_vs_hancock_joined
write_csv(dynamic_extrap_result_vs_hancock_joined,"outputs/dynamic_extrap_result_vs_hancock.csv")


#####do spatial interpolation


  # hierarchical regression parameters
  
  # doubly hierarchical version
  
  # between-classes
  beta_overall <- normal(0, 1, dim = n_covs)
  sigma_overall <- normal(0, 1, dim = n_covs, truncation = c(0, Inf))
  sigma_class <- normal(0, 1, dim = n_covs, truncation = c(0, Inf))
  
  # between-types
  beta_class_raw <- normal(0, 1, dim = c(n_covs, n_classes))
  beta_class_sigma <- sweep(beta_class_raw, 1, sigma_overall, FUN = "*")
  beta_class <- sweep(beta_class_sigma, 1, beta_overall, FUN = "+")
  
  # betas for types
  beta_type_raw <- normal(0, 1, dim = c(n_covs, n_types))
  beta_type_sigma <- sweep(beta_type_raw, 1, sigma_class, FUN = "*")
  beta_type <- beta_class[, classes_index] + beta_type_sigma
  
  # multiply through to get relative fitness of resistance for each insecticide
  # type at each cell
  
  effect_type <- exp(beta_type)
  
  
  selection_cell_years <- x_cell_years %*% effect_type
  fitness_cell_years <- 1 + selection_cell_years
  
  # reformat this in to a 3D array with dimensions:
  #   n_times x n_unique_cells x n_types x 1
  # to solve dynamics with time-varying fitness (time must be first, then other
  # two must match state variable, which has a trailing dimension of size 1)
  fitness_array <- fitness_cell_years
  dim(fitness_array) <- c(n_times, n_unique_cells, n_types, 1)
  
  # prior and minimum values for the initial fractions susceptible
  init_frac_prior <- ifelse(types == "DDT", 0.9, 0.95)
  init_frac_min <- ifelse(types == "DDT", 0.75, 0.9)
  
  # mean logit proportion of the distance from the minimum to 1, for each type
  init_frac_relative_prior <- (init_frac_prior - init_frac_min) / (1 - init_frac_min)
  logit_init_mean <- normal(qlogis(init_frac_relative_prior), 1, dim = n_types)
  
  # variance at each level, for each insecticide
  init_region_sd <- normal(0, 1, truncation = c(0, Inf), dim = n_types)
  init_country_sd <- normal(0, 1, truncation = c(0, Inf), dim = n_types)
  
  # unscaled deviation parameters (hierarchical decentring)
  init_region_raw <- normal(0, 1, dim = c(n_regions, dim = n_types))
  init_country_raw <- normal(0, 1, dim = c(n_countries, dim = n_types))
  
  # combine to get the regional and country-level deviation from the prior logit
  # mean
  init_region_effect <- sweep(init_region_raw, 2, init_region_sd, FUN = "*")
  init_country_effect <- sweep(init_country_raw, 2, init_country_sd, FUN = "*")
  
  country_region_index <- df %>%
    group_by(country_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(country_id, region_id) %>%
    arrange(country_id) %>%
    pull(region_id)
  
  # combine these together into the country-level logit-mean initial value
  # logit_init_mean <- qlogis(init_frac_mean)
  init_country_overall_effect <- init_country_effect + init_region_effect[country_region_index, ]
  logit_init_country <- sweep(init_country_overall_effect,
                              2,
                              logit_init_mean,
                              FUN = "+")
  
  # convert from relative (0-1) to the constrained scale (above init_frac_min)
  init_country_relative <- ilogit(logit_init_country)
  init_range <- 1 - init_frac_min
  init_country_magnitude <- sweep(init_country_relative, 2, init_range, FUN = "*") 
  init_country <- sweep(init_country_magnitude, 2, init_frac_min, FUN = "+")
  
  cell_country_lookup <- df %>%
    group_by(cell_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(cell_id, country_id) %>%
    arrange(cell_id) %>%
    pull(country_id)
  
  # expand out to all observations
  init_array <- init_country[cell_country_lookup, ]
  
  # add a trailing dimension to match greta.dynamics interface
  dim(init_array) <- c(dim(init_array), 1)
  
  # iterate through time to get fraction without resistant allele for all years at all cells
  # with data
  dynamic_cells <- iterate_dynamic_function(
    transition_function = haploid_next,
    initial_state = init_array,
    niter = n_times,
    w = fitness_array,
    parameter_is_time_varying = c("w"),
    tol = 0)
  
  # pull out the values at the cells, insecticides, and years corresponding to the
  # data
  
  # pull out this training data
  fold_train_df <- spatial_interpolation$training
  
  index <- cbind(fold_train_df$cell_id, fold_train_df$type_id, fold_train_df$year_id)
  fraction_susceptible_vec <- dynamic_cells$all_states[index]
  
  # resistant phenotype expression rate = expected mortality in bioassays
  population_mortality_vec <- fraction_susceptible_vec
  
  # define observation model
  rho_classes <- normal(0, 0.025,
                        truncation = c(0, 1),
                        dim = n_classes)
  
  distribution(fold_train_df$died) <- betabinomial_p_rho(
    N = fold_train_df$mosquito_number,
    p = population_mortality_vec,
    rho = rho_classes[fold_train_df$class_id])
  
  m <- model(
    beta_overall,
    sigma_overall,
    beta_type_raw,
    rho_classes
  )
  
  n_chains <- 4
  
  system.time(
    draws <- mcmc(m,
                  chains = n_chains,
                  n_samples = 1e3,
                  warmup = 1e3)
  )
  
  
  fold_test_df <- spatial_interpolation$test
  
  index_test <- cbind(fold_test_df$cell_id, fold_test_df$type_id, fold_test_df$year_id)
  fraction_susceptible_vec_test <- dynamic_cells$all_states[index_test]
  population_mortality_vec_test <- fraction_susceptible_vec_test
  
  predicted_test <- calculate(population_mortality_vec_test,
                              values = draws,
                              nsim = 1e3)
  
  rho_classes_test <- calculate(rho_classes,
                                values = draws,
                                nsim = 1e3)
  
  predicted_test_mean <- apply(predicted_test[[1]],2:3,mean) %>% as.numeric()
  
  # add predicted value back to test data for bookkeeping
  
  spatial_interpolation$test <- tibble(spatial_interpolation$test, predicted_dynamic = predicted_test_mean)
  
  rho_classes_test_mean <- apply(rho_classes_test[[1]],2:3,mean) %>% as.numeric()
  
  rho_classes_test_indexed <- rho_classes_test_mean[fold_test_df$class_id]
  
  test_outcome_df <- optimal_nn_preds %>%
    filter(
      experiment == "spatial_interpolation"
    ) %>% 
    select(died,mosquito_number,observed,year_start) %>% 
    mutate(predicted_test_mean = predicted_test_mean,
           rho_classes_test_indexed = rho_classes_test_indexed) %>%
    group_by(
      year_start
    ) %>% 
    summarise(
      pred_error_dynamic = betabinom_dev(died = died,
                                       mosquito_number = mosquito_number,
                                       predicted = predicted_test_mean,
                                       rho = rho_classes_test_indexed),
      bias_dynamic = mean(predicted_test_mean - observed),
      .groups = "drop"
    ) 
  
  
  # result df
  dynamic_interp_result <- optimal_nn_preds %>%
    filter(
      experiment == "spatial_interpolation"
    ) %>% 
    group_by(year_start) %>% 
    summarise(
      pred_error_nn = betabinom_dev(died = died,
                                    mosquito_number = mosquito_number,
                                    predicted = predicted),
      bias_nn = mean(predicted - observed)
    ) %>% 
    left_join(test_outcome_df,by = join_by(year_start)) %>% 
    left_join(spatial_interpolation_intercept_error,by = join_by(year_start)) %>% 
    mutate(experiment = "dynmaic_spatial_interpolation")
  
  # save to csv
  dynamic_interp_result
  write_csv(dynamic_interp_result,"outputs/dynamic_interp_result.csv")
  
# vs hancock
  test_outcome_df_vs_hancock <- optimal_nn_preds %>%
    filter(
      experiment == "spatial_interpolation"
    ) %>% 
    select(died,mosquito_number,observed,year_start,in_hancock) %>% 
    mutate(predicted_test_mean = predicted_test_mean,
           rho_classes_test_indexed = rho_classes_test_indexed) %>%
    filter(in_hancock) %>% 
    group_by(
      year_start
    ) %>% 
    summarise(
      pred_error_dynamic = betabinom_dev(died = died,
                                         mosquito_number = mosquito_number,
                                         predicted = predicted_test_mean,
                                         rho = rho_classes_test_indexed),
      bias_dynamic = mean(predicted_test_mean - observed),
      .groups = "drop"
    ) 
  
  dynamic_interp_result_vs_hancock <- hancock_interp_result %>% 
    left_join(test_outcome_df_vs_hancock,by = join_by(year_start)) %>% 
    mutate(experiment = "dynmaic_spatial_interpolation_vs_hancock") 
  
  dynamic_interp_result_vs_hancock
  write_csv(dynamic_interp_result_vs_hancock,"outputs/dynamic_interp_result_vs_hancock.csv")
  
  purge_greta_model()
  gc()


  #####do temporal forecast
  
  
  # hierarchical regression parameters
  
  # doubly hierarchical version
  
  # between-classes
  beta_overall <- normal(0, 1, dim = n_covs)
  sigma_overall <- normal(0, 1, dim = n_covs, truncation = c(0, Inf))
  sigma_class <- normal(0, 1, dim = n_covs, truncation = c(0, Inf))
  
  # between-types
  beta_class_raw <- normal(0, 1, dim = c(n_covs, n_classes))
  beta_class_sigma <- sweep(beta_class_raw, 1, sigma_overall, FUN = "*")
  beta_class <- sweep(beta_class_sigma, 1, beta_overall, FUN = "+")
  
  # betas for types
  beta_type_raw <- normal(0, 1, dim = c(n_covs, n_types))
  beta_type_sigma <- sweep(beta_type_raw, 1, sigma_class, FUN = "*")
  beta_type <- beta_class[, classes_index] + beta_type_sigma
  
  # multiply through to get relative fitness of resistance for each insecticide
  # type at each cell
  
  effect_type <- exp(beta_type)
  
  
  selection_cell_years <- x_cell_years %*% effect_type
  fitness_cell_years <- 1 + selection_cell_years
  
  # reformat this in to a 3D array with dimensions:
  #   n_times x n_unique_cells x n_types x 1
  # to solve dynamics with time-varying fitness (time must be first, then other
  # two must match state variable, which has a trailing dimension of size 1)
  fitness_array <- fitness_cell_years
  dim(fitness_array) <- c(n_times, n_unique_cells, n_types, 1)
  
  # prior and minimum values for the initial fractions susceptible
  init_frac_prior <- ifelse(types == "DDT", 0.9, 0.95)
  init_frac_min <- ifelse(types == "DDT", 0.75, 0.9)
  
  # mean logit proportion of the distance from the minimum to 1, for each type
  init_frac_relative_prior <- (init_frac_prior - init_frac_min) / (1 - init_frac_min)
  logit_init_mean <- normal(qlogis(init_frac_relative_prior), 1, dim = n_types)
  
  # variance at each level, for each insecticide
  init_region_sd <- normal(0, 1, truncation = c(0, Inf), dim = n_types)
  init_country_sd <- normal(0, 1, truncation = c(0, Inf), dim = n_types)
  
  # unscaled deviation parameters (hierarchical decentring)
  init_region_raw <- normal(0, 1, dim = c(n_regions, dim = n_types))
  init_country_raw <- normal(0, 1, dim = c(n_countries, dim = n_types))
  
  # combine to get the regional and country-level deviation from the prior logit
  # mean
  init_region_effect <- sweep(init_region_raw, 2, init_region_sd, FUN = "*")
  init_country_effect <- sweep(init_country_raw, 2, init_country_sd, FUN = "*")
  
  country_region_index <- df %>%
    group_by(country_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(country_id, region_id) %>%
    arrange(country_id) %>%
    pull(region_id)
  
  # combine these together into the country-level logit-mean initial value
  # logit_init_mean <- qlogis(init_frac_mean)
  init_country_overall_effect <- init_country_effect + init_region_effect[country_region_index, ]
  logit_init_country <- sweep(init_country_overall_effect,
                              2,
                              logit_init_mean,
                              FUN = "+")
  
  # convert from relative (0-1) to the constrained scale (above init_frac_min)
  init_country_relative <- ilogit(logit_init_country)
  init_range <- 1 - init_frac_min
  init_country_magnitude <- sweep(init_country_relative, 2, init_range, FUN = "*") 
  init_country <- sweep(init_country_magnitude, 2, init_frac_min, FUN = "+")
  
  cell_country_lookup <- df %>%
    group_by(cell_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(cell_id, country_id) %>%
    arrange(cell_id) %>%
    pull(country_id)
  
  # expand out to all observations
  init_array <- init_country[cell_country_lookup, ]
  
  # add a trailing dimension to match greta.dynamics interface
  dim(init_array) <- c(dim(init_array), 1)
  
  # iterate through time to get fraction without resistant allele for all years at all cells
  # with data
  dynamic_cells <- iterate_dynamic_function(
    transition_function = haploid_next,
    initial_state = init_array,
    niter = n_times,
    w = fitness_array,
    parameter_is_time_varying = c("w"),
    tol = 0)
  
  # pull out the values at the cells, insecticides, and years corresponding to the
  # data
  
  # pull out this training data
  fold_train_df <- temporal_forecasting$training
  
  index <- cbind(fold_train_df$cell_id, fold_train_df$type_id, fold_train_df$year_id)
  fraction_susceptible_vec <- dynamic_cells$all_states[index]
  
  # resistant phenotype expression rate = expected mortality in bioassays
  population_mortality_vec <- fraction_susceptible_vec
  
  # define observation model
  rho_classes <- normal(0, 0.025,
                        truncation = c(0, 1),
                        dim = n_classes)
  
  distribution(fold_train_df$died) <- betabinomial_p_rho(
    N = fold_train_df$mosquito_number,
    p = population_mortality_vec,
    rho = rho_classes[fold_train_df$class_id])
  
  m <- model(
    beta_overall,
    sigma_overall,
    beta_type_raw,
    rho_classes
  )
  
  n_chains <- 4
  
  system.time(
    draws <- mcmc(m,
                  chains = n_chains,
                  n_samples = 1e3,
                  warmup = 1e3)
  )
  
  
  fold_test_df <- temporal_forecasting$test
  
  index_test <- cbind(fold_test_df$cell_id, fold_test_df$type_id, fold_test_df$year_id)
  fraction_susceptible_vec_test <- dynamic_cells$all_states[index_test]
  population_mortality_vec_test <- fraction_susceptible_vec_test
  
  predicted_test <- calculate(population_mortality_vec_test,
                              values = draws,
                              nsim = 1e3)
  
  rho_classes_test <- calculate(rho_classes,
                                values = draws,
                                nsim = 1e3)
  
  predicted_test_mean <- apply(predicted_test[[1]],2:3,mean) %>% as.numeric()
  
  # add predicted value back to test data for bookkeeping
  
  temporal_forecasting$test <- tibble(temporal_forecasting$test, predicted_dynamic = predicted_test_mean)
  
  rho_classes_test_mean <- apply(rho_classes_test[[1]],2:3,mean) %>% as.numeric()
  
  rho_classes_test_indexed <- rho_classes_test_mean[fold_test_df$class_id]
  
  test_outcome_df <- optimal_nn_preds %>%
    filter(
      experiment == "temporal_forecasting"
    ) %>% 
    select(died,mosquito_number,observed,year_start) %>% 
    mutate(predicted_test_mean = predicted_test_mean,
           rho_classes_test_indexed = rho_classes_test_indexed) %>%
    group_by(
      year_start
    ) %>% 
    summarise(
      pred_error_dynamic = betabinom_dev(died = died,
                                         mosquito_number = mosquito_number,
                                         predicted = predicted_test_mean,
                                         rho = rho_classes_test_indexed),
      bias_dynamic = mean(predicted_test_mean - observed),
      .groups = "drop"
    ) 
  
  
  # result df
  dynamic_temporal_forecast_result <- optimal_nn_preds %>%
    filter(
      experiment == "temporal_forecasting"
    ) %>% 
    group_by(year_start) %>% 
    summarise(
      pred_error_nn = betabinom_dev(died = died,
                                    mosquito_number = mosquito_number,
                                    predicted = predicted),
      bias_nn = mean(predicted - observed)
    ) %>% 
    left_join(test_outcome_df,by = join_by(year_start))%>% 
    left_join(temporal_forecasting_intercept_error,by = join_by(year_start)) %>% 
    mutate(experiment = "dynmaic_temporal_forecasting")
  
  # quick_regression <- lm(bias ~ model + country_name + insecticide_type, 
  #    data = dynamic_temporal_forecast_result %>% 
  #      pivot_longer(cols = starts_with("bias"),names_to = "model",values_to = "bias"))
  # 
  # anova(quick_regression)
  
  # save to csv
  View(dynamic_temporal_forecast_result)
  write_csv(dynamic_temporal_forecast_result,"outputs/dynamic_temporal_forecast_result.csv")
  
  purge_greta_model()
  gc()

  # plot
  
  dynamic_temporal_forecast_result %>% 
    pivot_longer(starts_with("pred_error"), 
                 values_to = "deviance", 
                 names_to = "model") %>% 
    mutate(model = stringr::word(model,3,3, sep = "_"),
           model = if_else(model == "nn", "nearest neighbour", model)) %>% 
    select(year_start,model, deviance) %>% 
    ggplot(aes(x = year_start, y = deviance, col = model)) + 
    theme_minimal() +
    scale_x_continuous(n.breaks = 3, name = "year") + 
    geom_point(pch = 1, size = 2) + 
    ggtitle("Predictive deviance in temporal forecast validation experiment")
  
  ggsave("figures/temporal_forecast_CV_deviance.png",
         bg = "white",
         scale = 0.8,
         width = 8,
         height = 8)
  
  dynamic_temporal_forecast_result %>% 
    pivot_longer(starts_with("bias"), 
                 values_to = "bias", 
                 names_to = "model") %>% 
    mutate(model = stringr::word(model,2,2, sep = "_"),
           model = if_else(model == "nn", "nearest neighbour", model)) %>% 
    select(year_start,model, bias) %>% 
    ggplot(aes(x = year_start, y = bias, col = model)) + 
    theme_minimal() +
    scale_x_continuous(n.breaks = 3, name = "year") + 
    geom_point(pch = 2, size = 2) + 
    ggtitle("Predictive bias in temporal forecast validation experiment")
  
  ggsave("figures/temporal_forecast_CV_bias.png",
         bg = "white",
         scale = 0.8,
         width = 8,
         height = 8)

  dynamic_interp_result %>% 
    pivot_longer(starts_with("pred_error"), 
                 values_to = "deviance", 
                 names_to = "model") %>% 
    mutate(model = stringr::word(model,3,3, sep = "_"),
           model = if_else(model == "nn", "nearest neighbour", model)) %>% 
    select(year_start,model, deviance) %>% 
    ggplot(aes(x = year_start, y = deviance, col = model)) + 
    theme_minimal() +
    scale_x_continuous(breaks = scales::breaks_pretty(),name = "year") + 
    geom_point(pch = 1, size = 2) + 
    ggtitle("Predictive deviance in spatial interpolation validation experiment")
  
  ggsave("figures/spatial_interpolation_CV_deviance.png",
         bg = "white",
         scale = 0.8,
         width = 8,
         height = 8)
  
  dynamic_interp_result %>% 
    pivot_longer(starts_with("bias"), 
                 values_to = "bias", 
                 names_to = "model") %>% 
    mutate(model = stringr::word(model,2,2, sep = "_"),
           model = if_else(model == "nn", "nearest neighbour", model)) %>% 
    select(year_start,model, bias) %>% 
    ggplot(aes(x = year_start, y = bias, col = model)) + 
    theme_minimal() +
    scale_x_continuous(breaks = scales::breaks_pretty(), name = "year") + 
    geom_point(pch = 2, size = 2) + 
    ggtitle("Predictive bias in spatial interpolation validation experiment")
  
  ggsave("figures/spatial_interpolation_CV_bias.png",
         bg = "white",
         scale = 0.8,
         width = 8,
         height = 8)
  
  dynamic_extrap_result %>% 
    pivot_longer(starts_with("pred_error"), 
                 values_to = "deviance", 
                 names_to = "model") %>% 
    mutate(model = stringr::word(model,3,3, sep = "_"),
           model = if_else(model == "nn", "nearest neighbour", model)) %>% 
    select(country_name,model, deviance) %>% 
    ggplot(aes(x = country_name, y = deviance, col = model)) + 
    theme_minimal() +
    scale_x_discrete(name = "country") + 
    geom_point(pch = 1, size = 2) + 
    ggtitle("Predictive deviance in spatial extrapolation validation experiment")
  
  ggsave("figures/spatial_extrapolation_CV_deviance.png",
         bg = "white",
         scale = 0.8,
         width = 8,
         height = 8)
  
  
  dynamic_extrap_result %>% 
    pivot_longer(starts_with("bias"), 
                 values_to = "bias", 
                 names_to = "model") %>% 
    mutate(model = stringr::word(model,2,2, sep = "_"),
           model = if_else(model == "nn", "nearest neighbour", model)) %>% 
    select(country_name,model, bias) %>% 
    ggplot(aes(x = country_name, y = bias, col = model)) + 
    theme_minimal() +
    scale_x_discrete( name = "country") + 
    geom_point(pch = 2, size = 2) + 
    ggtitle("Predictive bias in spatial extrapolation validation experiment")
  
  ggsave("figures/spatial_extrapolation_CV_bias.png",
         bg = "white",
         scale = 0.8,
         width = 8,
         height = 8)
  