# fit model

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load covariate rasters

# load and temporally average raster data
nets_cube <- rast("data/clean/nets_per_capita_cube.tif")
nets_flat <- terra::app(nets_cube, mean)
names(nets_flat) <- "itn_percap"
nets_flat_std <- nets_flat / global(nets_flat, "max", na.rm = TRUE)[1, 1]

# load the other layers
covs_flat_orig <- rast("data/clean/flat_covariates.tif")

# add on the scaled, flat nets layer
covs_flat <- c(nets_flat_std, covs_flat_orig)

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
                             "Organophosphates")
  ) %>%
  mutate(
    # convert concentrations into numeric values
    concentration = as.numeric(str_remove(insecticide_conc, "%")),
    # create a time index, relative to the baseline year
    time = year_start - baseline_year
  ) %>%
  filter(
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
  # extract non-temporal covariates at coordinates
  bind_cols(
    terra::extract(covs_flat,
                   select(., longitude, latitude) %>%
                     as.matrix())
  ) %>%
  # drop a handful of points mossing covariates
  filter(
    !is.na(itn_percap)
  )

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

# # summarise the data
# df %>%
#   select(insecticide_class,
#          insecticide_type,
#          concentration,
#          test_type) %>%
#   group_by(insecticide_class,
#            insecticide_type,
#            concentration,
#            test_type) %>%
#   summarise(
#     count = n()
#   ) %>%
#   arrange(desc(count),
#           insecticide_class,
#           insecticide_type)

# which covariates to use (all by default)
covariate_names <- names(covs_flat)

# create design matrix
x <- df %>%
  select(any_of(covariate_names)) %>%
  as.matrix()

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
n_covs <- ncol(x)
n_classes <- length(classes)
n_types <- length(types)

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


# multiply through to get log relative fitness of resistance for each
# insecticide type

# convert beta to positive effect sizes
effect_type <- exp(beta_type)

# this way is marginally faster and probably quite a lot less memory intensive
# than the matrix multiply
effect_mat <- t(effect_type)[df$type_id,]
eta_mat <- x * effect_mat
eta_vec <- greta:::rowSums(eta_mat)

# pull out the etas (observations, by insecticides) corresponding to data
# eta <- x %*% effect_type
# data_index <- cbind(seq_len(n_obs), df$type_id)
# eta_vec <- eta[data_index]

# model fractions susceptible in 2000

# # hierarchical prior (with hierarchical decentring) on the 
# logit_init_overall <- normal(qlogis(0.999), 1)
# logit_init_sd <- normal(0, 0.5, truncation = c(0, Inf))
# logit_init_raw <- normal(0, 1, dim = n_types)
# logit_init <- logit_init_overall + logit_init_raw * logit_init_sd
# init_fraction_susceptible <- ilogit(logit_init)

# more straightforward prior with a harder limit
init_fraction_susceptible <- normal(1, 0.01,
                                    dim = n_types,
                                    truncation = c(0, 1))
# sims <- calculate(init_fraction_susceptible[1], nsim = 1e4)[[1]]
# hist(sims, xlim = c(0.8, 1))
# range(sims)
# compute the fraction susceptible over time against each insecticide

# more stable version of this:
#   fitness = exp(eta)
#   cumulative_fitness = fitness ^ df$time
cumulative_fitnesses_vec <- exp(eta_vec * df$time)

init_vec <- init_fraction_susceptible[df$type_id]
fraction_susceptible_vec <- init_vec / 
  (init_vec + (1 - init_vec) * cumulative_fitnesses_vec)

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

# define betabinomial observation model

# need to constrain this error for now, until there are covariates to explain it
# better
# observation_extra_error <- normal(0, 0.1, truncation = c(0, Inf))

# compute the expected sd of the binomial (assuming independent samples)
binomial_sd <- sqrt(df$mosquito_number *
                      population_mortality_vec *
                      (1 - population_mortality_vec))

# model the observation (betabinomial) sd as a multiplier on the binomial sd
observation_sd_multiplier_raw <- normal(0, 0.1, truncation = c(0, Inf))
observation_sd_multiplier <- 1 + observation_sd_multiplier_raw

# summary(calculate(observation_sd_multiplier, nsim = 1e6)[[1]])
# observation_sd_multiplier <- 1
observation_sd <- binomial_sd * observation_sd_multiplier

a <- ((1 - population_mortality_vec) / observation_sd ^ 2 - 1 / population_mortality_vec) * population_mortality_vec ^ 2
b <- a * (1 / population_mortality_vec - 1)
distribution(df$died) <- beta_binomial(df$mosquito_number, a, b)

# or just use binomial
# distribution(df$died) <- binomial(df$mosquito_number, population_mortality_vec)

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
# user  system elapsed 
# 605.339  71.699 162.177 

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

n_plot <- 100
type_idx <- 1
times_plot <- seq(1990, max(df$year_start), length.out = n_plot)
time_diff_plot <- sweep(matrix(times_plot, n_plot, n_types), 2, 2000, FUN = "-")

x_fixed <- t(c(0.5, rep(0, n_covs - 1)))
eta_plot <- x_fixed %*% exp(beta_type)
cumulative_fitnesses_plot <- exp(sweep(time_diff_plot, 2, t(eta_plot), FUN = "*"))

# compute the fraction susceptible over time against each insecticide
init_mat <- sweep(zeros(n_plot, n_types), 2, init_fraction_susceptible, FUN = "+")
fraction_susceptible_plot <- init_mat /
  (init_mat + (1 - init_mat) * cumulative_fitnesses_plot)
population_mortality_plot <- fraction_susceptible_plot

post_plot_sims <- calculate(
  population_mortality_plot,
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
                        beta_type,
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

par(mfrow = n2mfrow(n_types))
for (i in seq_len(n_types)) {
  hist(posteriors$beta_type[, , i],
       breaks = 100,
       main = types[i])
  
}

# make prediction rasters

# get all cells in mastergrids, and extract the values there
cells <- terra::cells(covs_flat[[1]])

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

x_predict <- terra::extract(covs_flat, cells) %>%
  as_tibble() %>%
  select(any_of(covariate_names)) %>%
  as.matrix()

# years to predict to
years_predict <- 2000:2030
times_predict <- years_predict - baseline_year

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
# batch_bigness <- length(cells) + 1
batch_bigness <- 1e5
batch_idx <- seq_along(cells) %/% batch_bigness
# index to raster, for setting cell values
cell_batches <- split(cells, batch_idx)
# index to x_predict
x_batches <- split(seq_along(cells), batch_idx)
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
    x_batch <- x_batches[[batch_index]]
    
    type_id_predict <- match(this_insecticide, types)
    
    # get log fitness for population resistant to this insecticide, at these
    # cells
    x_predict_batch <- x_predict[x_batch, ]
    effect_type <- exp(beta_type)
    eta_predict <- x_predict_batch %*% effect_type[, type_id_predict]
    
    # accumulate over time and convert to fitness
    times_predict_mat <- outer(rep(1, length(cell_batch)),
                               times_predict,
                               FUN = "*")
    cumulative_eta_predict <- sweep(as_data(times_predict_mat),
                                    1,
                                    eta_predict,
                                    FUN = "*")
    cumulative_fitnesses_predict <- exp(cumulative_eta_predict)
    
    # compute the fraction susceptible to this insecticide over time
    init_predict <- init_fraction_susceptible[type_id_predict]
    fraction_susceptible_predict <- init_predict / 
      (init_predict + (1 - init_predict) * cumulative_fitnesses_predict)
    population_mortality_predict <- fraction_susceptible_predict
    
    # get posterior draws of these, fixing the RNG seed so it's the same
    # collection of posterior samples for all batches
    # set.seed(this_seed)
    pred_mat <- calculate(population_mortality_predict,
                          values = draws,
                          seed = this_seed,
                          nsim = 1e2)[[1]]
    
    # compute posterior mean
    batch_pred_mean <- apply(pred_mat, 2:3, mean)
    
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
