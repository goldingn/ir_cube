# fit model

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load covariate rasters

# load, temporally average, and scale net coverage data
nets_cube <- rast("data/clean/nets_per_capita_cube.tif")
nets_flat <- terra::app(nets_cube, mean)
names(nets_flat) <- "itn_percap"
nets_flat_std <- nets_flat / global(nets_flat, "max", na.rm = TRUE)[1, 1]

# load the other layers
covs_flat_orig <- rast("data/clean/flat_covariates.tif")

# add on the scaled, flat nets layer
covs_flat <- c(nets_flat_std, covs_flat_orig)

mask <-  covs_flat[[1]] * 0

# load data
ir_mtm_africa <- readRDS(file = "data/clean/mtm_data.RDS")

# set the start of the timeseries considered
baseline_year <- 2000

df <- ir_mtm_africa %>%
  filter(
    # subset to An. gambiae (s.l./s.s.)
    species %in% c("An. gambiae s.l.", "An. gambiae s.s."),
    # subset to WHO tube tests (most of data)
    test_type == "WHO tube test",
    # drop the minor classes: pyrroles and neonicotinoids
    insecticide_class %in% c("Pyrethroids",
                             "Carbamates",
                             "Organochlorines",
                             "Organophosphates"),
    # drop Dieldrin as it has two concentrations in the dataset, and those with
    # fewer than 200 observations (manually as I'm being lazy)
    !(insecticide_type %in% c("Dieldrin",
                              "Carbosulfan",
                              "Cyfluthrin",
                              "Etofenprox",
                              "Propoxur")),
    # drop any from before when we have data on net coverage
    year_start >= baseline_year
  ) %>%
  mutate(
    # convert concentrations into numeric values
    concentration = as.numeric(str_remove(insecticide_conc, "%")),
    # create an index to the simulation year (in 1-indexed integers)
    year_id = year_start - baseline_year + 1,
    # add on cell ids corresponding to these observations,
    cell = cellFromXY(mask, 
                      as.matrix(select(., longitude, latitude)))
  ) %>%
  # drop a handful of datapoints missing covariates
  filter(
    !is.na(extract(mask, cell)[, 1])
  ) %>%
  # add an index to the vector of unique cells (now they have been subsetted)
  mutate(
    cell_id = match(cell, unique(cell))
  )

# pull covariates for all unique cells
unique_cells <- unique(df$cell)
n_unique_cells <- length(unique_cells)

# create design matrix at all unique cells
x_cells <- extract(covs_flat, unique_cells)

# # extract spatiotemporal covariates from the cube
# nets_lookup <- terra::extract(nets_cube, coords) %>%
#   mutate(
#     row = row_number(),
#     .before = everything()
#   ) %>%
#   pivot_longer(cols = starts_with("nets_"),
#                names_prefix = "nets_",
#                names_to = "year",
#                values_to = "net_coverage") %>%
#   mutate(
#     year = as.numeric(year)
#   )

# map the data to the insecticide classes and types
classes <- unique(df$insecticide_class)
types <- unique(df$insecticide_type)

# add in indices to the dataframe to the types and insecticides
df$class_id <- match(df$insecticide_class, classes)
df$type_id <- match(df$insecticide_type, types)

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
  
# dimensions
n_obs <- nrow(df)
n_covs <- ncol(x_cells)
n_classes <- length(classes)
n_types <- length(types)
n_times <- max(df$year_start) - min(df$year_start) + 1

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


# # simpler single-level hierarchy
# beta_overall <- normal(0, 1, dim = n_covs) * 10
# sigma_overall <- normal(0, 1, dim = n_covs, truncation = c(0, Inf))
# 
# # hierarchical decentring implementation
# beta_type_raw <- normal(0, 1, dim = c(n_covs, n_types))
# beta_type_sigma <- sweep(beta_type_raw, 1, sigma_overall, FUN = "*")
# beta_type <- sweep(beta_type_sigma, 1, beta_overall, FUN = "+")


# remove hierarchy for debugging
# beta_type <- beta_type_raw


# multiply through to get relative fitness of resistance for each insecticide
# type at each cell

# convert beta to positive effect sizes, them matrix multiply with indicator
# variables to constrain to positive effects of the covariates on selection
# coefficient/relative fitness
effect_type <- exp(beta_type)

# # compute selection coefficients, convert to relative fitness, and extract at
# # data locations
selection_cells <- x_cells %*% effect_type
fitness_cells <- 1 + selection_cells

# model fractions susceptible prior to 2000
init_fraction_susceptible <- normal(1, 0.01,
                                    dim = n_types,
                                    truncation = c(0, 1))

# compute the fraction susceptible over time against each insecticide via
# haploid trait evolution (equivalent to diploid evolution with heterozygote
# dominance = 0.5)

# solve through time to get fraction susceptible over time
haploid_next <- function(state, iter, w) {
  p <- state
  q <- 1 - p
  p / (p + q * w)
}

# expand initial conditions out to all cells with data (plus a trailing
# dimension to match greta.dynamics interface)
init_cells <- sweep(ones(n_unique_cells, n_types),
                    2,
                    init_fraction_susceptible,
                    FUN = "*")
dim(init_cells) <- c(dim(init_cells), 1)

# pad out the dimension of the fitness to match
dim(fitness_cells) <- c(dim(fitness_cells), 1)

# iterate through time to get fraction susceptible for all years at all cells
# with data
dynamic_cells <- iterate_dynamic_function(
  transition_function = haploid_next,
  initial_state = init_cells,
  niter = n_times,
  w = fitness_cells,
  tol = 0)

# pull out the values at the cells, insecticides, and years corresponding to the
# data
index <- cbind(df$cell_id, df$type_id, df$year_id)
fraction_susceptible_vec <- dynamic_cells$all_states[index]

# model the susceptible & resistant LD50s of the different insecticides
# (skip this until we have multiple concentrations)

# # for single concentration data, fix the variance to 1 and fix the LD50 for
# # susceptibles to correspond to mortality rate of 99% (98% is resistant)
# population_LD50_sd <- 1
# LD50_susceptible <- type_concentrations - population_LD50_sd * qnorm(0.99) 
# 
# # define a reasonable prior for the difference in LD50s: resistants have mortality ~1% (95% CI goes up to 25% mortality))
# # LD50_resistant_lower <- type_concentrations - population_LD50_sd * qnorm(0.25)
# LD50_resistant_mean <- type_concentrations - population_LD50_sd * qnorm(0.01)
# # LD50_resistant_sd <- LD50_resistant_mean - LD50_resistant_lower / 1.96
# # LD50_difference_mean <- LD50_resistant_mean - LD50_susceptible
# 
# # this is a chore
# # LD50_difference_raw <- normal(0, 1, dim = n_types)
# # LD50_difference <- LD50_difference_mean + LD50_difference_raw * LD50_resistant_sd
# 
# # compute LD50 for resistant genotype
# # LD50_resistant <- LD50_susceptible + LD50_difference
# LD50_resistant <- LD50_resistant_mean
# 
# # get the population-level LD50s for the observations
# LD50_vec <- fraction_susceptible_vec * LD50_susceptible[df$type_id] +
#   (1 - fraction_susceptible_vec) * LD50_resistant[df$type_id]
# 
# probit_vec <- (df$concentration - LD50_vec) / population_LD50_sd
# population_mortality_vec <- iprobit(probit_vec)

# for now, map straight from genotype to mortality
population_mortality_vec <- fraction_susceptible_vec

# define observation model

# compute the expected sd of a binomial distribution (ie. assuming independent
# samples)
binomial_sd <- sqrt(df$mosquito_number *
                      population_mortality_vec *
                      (1 - population_mortality_vec))

# model the observation (betabinomial) sd as a multiplier on the binomial sd,
# accounting for additional error due to nonindependent sampling of individuals
# from the population
observation_sd_multiplier_raw <- normal(0, 0.1, truncation = c(0, Inf))
observation_sd_multiplier <- 1 + observation_sd_multiplier_raw
observation_sd <- binomial_sd * observation_sd_multiplier

# compute parameters of the beta-binomial that match this mean and sd
a <- ((1 - population_mortality_vec) / observation_sd ^ 2 - 1 / population_mortality_vec) * population_mortality_vec ^ 2
b <- a * (1 / population_mortality_vec - 1)

# define the distribution over the observations
distribution(df$died) <- beta_binomial(df$mosquito_number, a, b)

m <- model(
  beta_overall,
  sigma_overall,
  beta_type_raw,
  init_fraction_susceptible,
  observation_sd_multiplier
)

# set the inits for observation_sd_multiplier to be large (more params will fit
# OK to data to start with, then chains can move towards better parts of
# parameter space)
n_chains <- 4
inits <- replicate(n_chains,
                   initials(
                     observation_sd_multiplier_raw = 15
                   ),
                   simplify = FALSE)

system.time(
  draws <- mcmc(m,
                initial_values = inits,
                chains = n_chains)
)
# user   system  elapsed 
# 4243.033 2195.452 1795.518 

# check convergence
coda::gelman.diag(draws,
                  autoburnin = FALSE,
                  multivariate = FALSE)

# compute fitted population mortality rates, and plot with bioassay data
sims <- calculate(population_mortality_vec,
                  fraction_susceptible_vec,
                  values = draws, nsim = 100)
fitted_post_mean <- colMeans(sims$population_mortality_vec[, , 1]) * 100
frac_susc_post_mean <- colMeans(sims$fraction_susceptible_vec[, , 1]) * 100

par(mfrow = c(1, 1))
plot(frac_susc_post_mean ~ jitter(df$year_start),
     cex = 0.1)

# calculate expected bioassay mortalities over time for each insecticide, with
# fixed assumption about covariate values, and overplot the data

x_fixed <- t(c(0.5, rep(0, n_covs - 1)))
selection_plot <- t(x_fixed %*% exp(beta_type))
fitness_plot <- 1 + selection_plot
times_plot <- baseline_year + seq_len(n_times) - 1

dynamic_outputs_plot <- iterate_dynamic_function(
  transition_function = haploid_next,
  initial_state = init_fraction_susceptible,
  niter = n_times,
  w = fitness_plot,
  tol = 0)

# pull out the data locations
fraction_susceptible_plot <- t(dynamic_outputs_plot$all_states)
population_mortality_plot <- fraction_susceptible_plot

post_plot_sims <- calculate(
  population_mortality_plot,
  fitness_plot,
  values = draws,
  nsim = 1000)

# summarise these for plotting
pop_mort_post_mean <- apply(post_plot_sims$population_mortality_plot,
                            2:3,
                            mean)
pop_mort_post_ci <- apply(post_plot_sims$population_mortality_plot,
                          2:3,
                          quantile,
                          c(0.025, 0.975))

par(mfrow = n2mfrow(n_types),
    mar = c(2, 2, 4, 2))
for (i in seq_len(n_types)) {
  plot(pop_mort_post_mean[, i] ~ times_plot,
       ylim = c(0, 1),
       type = "l")
  lines(pop_mort_post_ci[1, , i] ~ times_plot,
        lty = 2)
  lines(pop_mort_post_ci[2, , i] ~ times_plot,
        lty = 2)
  title(main = sprintf("%s\n(%s)",
                       types[i],
                       classes[classes_index[i]]))
  
  # overplot all the data for this type
  df %>%
    filter(type_id == i) %>%
    mutate(mortality = died / mosquito_number) %>%
    select(year_start, mortality) %>%
    points(cex = 0.5, lwd = 0.1)
  
}

# plot posteriors
posteriors <- calculate(init_fraction_susceptible,
                        effect_type,
                        observation_sd_multiplier,
                        values = draws,
                        nsim = 1000)

par(mfrow = c(1, 1),
    mar = c(5, 4, 4, 2) + 0.1)
hist(posteriors$observation_sd_multiplier,
     breaks = 100)

par(mfrow = n2mfrow(n_types))
for (i in seq_len(n_types)) {
  hist(posteriors$init_fraction_susceptible[, i, ],
       xlim = c(0.7, 1),
       breaks = 100,
       main = types[i])
  
}

# only shown for ITNs
par(mfrow = n2mfrow(n_types))
for (i in seq_len(n_types)) {
  hist(posteriors$effect_type[, 1, i],
       breaks = 100,
       main = types[i])
  
}

# make prediction rasters

# get all cells in mastergrids mask, and extract the values there
cells <- terra::cells(mask)

# nets_cube_dataframe <- terra::extract(nets_cube, cells) %>%
#   mutate(
#     cell = cells,
#     .before = everything()
#   ) %>%
#   pivot_longer(cols = starts_with("nets_"),
#                names_prefix = "nets_",
#                names_to = "year",
#                values_to = "net_coverage") %>%
#   mutate(
#     year = as.numeric(year)
#   )

x_predict <- terra::extract(covs_flat, cells)

# years to predict to
years_predict <- 2000:2030

n_times_predict <- length(years_predict)

# loop through insecticides making predictions in batches of cells for multiple years simultaneously
template_raster <- nets_cube$nets_2000 * 0
template_cube <- template_raster %>%
  replicate(length(years_predict),
            .,
            simplify = FALSE) %>%
  do.call(c, .)

names(template_cube) <- years_predict

# create batches of cells for processing

# NOTE due to a weird-as-hell greta bug, this object must not be called
# batch_size: https://github.com/greta-dev/greta/issues/634
batch_bigness <- 1e5
batch_idx <- seq_along(cells) %/% batch_bigness
# index to raster, for setting cell values
cell_batches <- split(cells, batch_idx)
# index to x_predict
x_row_batches <- split(seq_along(cells), batch_idx)
n_batches <- length(cell_batches)

# insecticide types to save
types_save <- c("Deltamethrin",
                "Permethrin",
                "Alpha-cypermethrin")

# loop through insecticides
for (this_insecticide in types_save) {
  
  # make a raster cube to stick predictions into
  this_ir_cube <- template_cube
  
  # grab the current seed, so we can do calculate in batches but not shuffle
  # parameters over space
  this_seed <- greta::.internals$utils$misc$get_seed()
  
  # loop through batches of cells
  for (batch_index in seq_len(n_batches)) {
    
    print(sprintf("processing batch %i of %i",
                  batch_index,
                  n_batches))
    
    cell_batch <- cell_batches[[batch_index]]
    x_row_batch <- x_row_batches[[batch_index]]
    batch_n <- length(cell_batch)
    
    # for this insecticide, tile the intial conditions for each cell
    type_id_predict <- match(this_insecticide, types)
    init_batch <- init_fraction_susceptible[type_id_predict] * ones(batch_n)
    
    # get log fitness for population resistant to this insecticide, at these
    # cells
    x_batch <- x_predict[x_row_batch, ]
    selection_batch <- x_batch %*% effect_type[, type_id_predict]
    fitness_batch <- 1 + selection_batch
    dynamics_batch <- iterate_dynamic_function(
      transition_function = haploid_next,
      initial_state = init_batch,
      niter = n_times_predict,
      w = fitness_batch,
      tol = 0)
    
    fraction_susceptible_batch <- dynamics_batch$all_states
    population_mortality_batch <- fraction_susceptible_batch
    
    # get posterior draws of these, fixing the RNG seed so it's the same
    # collection of posterior samples for all batches
    # set.seed(this_seed)
    pred_batch <- calculate(population_mortality_batch,
                          values = draws,
                          seed = this_seed,
                          nsim = 1e2)[[1]]
    
    # compute posterior mean over the draws, for cells and years
    batch_pred_mean <- apply(pred_batch, 2:3, mean)
    
    # stick this batch of predictions in the raster cube, one year at a time
    for (y in seq_along(years_predict)) {
      this_ir_cube[[y]][cell_batch] <- batch_pred_mean[, y]
    }

  }
  
  # now write out the raster
  write_path <- file.path("outputs/ir_maps",
                          this_insecticide,
                          sprintf("ir_%s_susceptibility.tif",
                                  this_year))
  
  # make sure the directory exists
  dir.create(dirname(write_path))
  writeRaster(this_ir_raster,
              write_path)  
}


# things to do next:
# - add in time-varying nets
# - set up code for posterior predictive checking
# - set up code for out of sample (future timesteps, spatial blocks) evaluation of model fit.
# - add in more discriminating bioassay data
# - add in intensity bioassay data (with LD50 model, including possibility of complete but imperfect resistance)
