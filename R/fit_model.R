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

# model overdispersion in the data via an overdispersion parameter rho. This
# prior makes rho approximately uniform, but fairly nicely behaved
logit_rho <- normal(0, 1.6)
rho <- ilogit(logit_rho)

distribution(df$died) <- betabinomial_p_rho(N = df$mosquito_number,
                                            p = population_mortality_vec,
                                            rho = rho)

m <- model(
  beta_overall,
  sigma_overall,
  beta_type_raw,
  init_fraction_susceptible,
  logit_rho
)

# set the inits for obs_var_multiplier to be large (more params will fit OK to
# data to start with, then chains can move towards better parts of parameter
# space)
n_chains <- 4
inits <- replicate(n_chains,
                   initials(
                     logit_rho = 10
                   ),
                   simplify = FALSE)

system.time(
  draws <- mcmc(m,
                initial_values = inits,
                chains = n_chains)
)
# user   system  elapsed 
# 4121.944 2225.062 1789.773 

# check convergence
coda::gelman.diag(draws,
                  autoburnin = FALSE,
                  multivariate = FALSE)

# find some locations with lots of years of data and plot predictions and data
# for these
locations_plot <- df %>%
  group_by(latitude, longitude) %>%
  filter(
    n_distinct(year_start) >= 8
  ) %>%
  select(
    country_name,
    latitude,
    longitude,
    cell_id
  ) %>%
  distinct() %>%
  # geocode then find more interpretable names (google to see if this is how
  # they are referred to in IR papers)
  reverse_geocode(
    lat = latitude,
    long = longitude,
    method = 'osm',
    full_results = TRUE
  ) %>%
  mutate(
    place = case_when(
      address == "Soumousso, Houet, Hauts-Bassins, Burkina Faso" ~ "Soumousso, Burkina Faso",
      address == "Katito-Kendu Bay-Homa Bay road, Oriang, Central ward, Karachuonyo, Homa Bay, Nyanza, Kenya" ~ "Homa Bay, Kenya",
      address == "RNIE 7, Kandi 3, Kandi, Alibori, Bénin" ~ "Kandi, Benin",
      address == "Aménagement de périmètre maraîcher, Djougou, Donga, Bénin" ~ "Djougou, Benin"
    )
  ) %>%
  select(
    place,
    country_name,
    latitude,
    longitude,
    cell_id
  )

# pull these out for plotting
preds_plot_setup <- expand_grid(
  cell_id = locations_plot$cell_id,
  year_start = baseline_year:max(df$year_start),
  insecticide_type = types
) %>%
  mutate(
    year_id = year_start - baseline_year + 1,
    type_id = match(insecticide_type, types)
  )

index_plot <- preds_plot_setup %>%
  select(cell_id,
         type_id,
         year_id) %>%
  as.matrix()

fraction_susceptible_plot <- dynamic_cells$all_states[index_plot]
population_mortality_plot <- fraction_susceptible_plot

# simulate mortalities under binomial sampling
sample_size_plot <- 100
binomial_died <- binomial(sample_size_plot, population_mortality_plot)
binomial_mortality <- binomial_died / sample_size_plot

# simulate mortalities under betabinomial sampling
betabinomial_died <- betabinomial_p_rho(N = sample_size_plot,
                                              p = population_mortality_plot,
                                              rho = rho)
betabinomial_mortality <- betabinomial_died / sample_size_plot

sims <- calculate(population_mortality_plot,
                  binomial_mortality,
                  betabinomial_mortality,
                  values = draws,
                  nsim = 1000)

# posterior mean mortality rate, and intervals for the population value
# (posterior uncertainty), and posterior predictive intervals (posterior
# uncertainty and sampling error) for observed mortality from binomial sampling,
# and from betabinomial sampling
pop_mort_mean <- colMeans(sims$population_mortality_plot[, , 1])
pop_mort_ci <- apply(sims$population_mortality_plot[, , 1],
                     2,
                     quantile,
                     c(0.025, 0.975))
binomial_mort_ci <- apply(sims$binomial_mortality[, , 1],
                          2,
                          quantile,
                          c(0.025, 0.975))
betabinomial_mort_ci <- apply(sims$betabinomial_mortality[, , 1],
                              2,
                              quantile,
                              c(0.025, 0.975))

insecticides_plot <- c("Deltamethrin", "Permethrin", "Bendiocarb")

preds_plot <- preds_plot_setup %>%
  mutate(
    mortality = pop_mort_mean,
    pop_lower = pop_mort_ci[1, ],
    pop_upper = pop_mort_ci[2, ],
    binomial_lower = binomial_mort_ci[1, ],
    binomial_upper = binomial_mort_ci[2, ],
    betabinomial_lower = betabinomial_mort_ci[1, ],
    betabinomial_upper = betabinomial_mort_ci[2, ],
  ) %>%
  filter(
    insecticide_type %in% insecticides_plot
  ) %>%
  left_join(
    locations_plot,
    by = "cell_id"
  )

points_plot <- df %>%
  filter(
    cell_id %in% locations_plot$cell_id,
    insecticide_type %in% insecticides_plot
  ) %>%
  mutate(
    mortality = died / mosquito_number
  ) %>%
  left_join(
    locations_plot,
    by = "cell_id"
  )
  

# plot these, then add cell data over the top
preds_plot %>%
  ggplot(
    aes(
      x = year_start,
      y = mortality,
      group = place
    )
  ) +
  # betabinomial posterior predictive intervals on bioassay data (captures
  # sample size and non-independence effect)
  geom_ribbon(
    aes(
      ymax = betabinomial_lower,
      ymin = betabinomial_upper,
    ),
    colour = grey(0.4),
    linewidth = 0.25,
    linetype = 2,
    fill = grey(0.9)
  ) +
  # binomial posterior predictive intervals on bioassay data (captures sample
  # size but assumes independence, so underestimates variance)
  geom_ribbon(
    aes(
      ymax = binomial_lower,
      ymin = binomial_upper,
    ),
    fill = grey(0.6)
  ) +
  # credible intervals on population-level proportion (our best guess at the
  # 'truth')
  geom_ribbon(
    aes(
      ymax = pop_lower,
      ymin = pop_upper,
    ),
    fill = grey(0.4)
  ) +
  scale_y_continuous(labels = scales::percent) +
  facet_grid(place ~ insecticide_type) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  geom_point(
    aes(
      x = year_start,
      y = mortality,
      group = place
    ),
    data = points_plot,
  )

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
# - add in plotting of a few specific data timeseries
# - add in time-varying nets
# - maybe add a resistance cost parameter
#    (positive, as fitness = 1 + selection - cost)
# - set up code for posterior predictive checking
# - set up code for out of sample (future timesteps, spatial blocks) evaluation
#    of model fit.
# - add in more discriminating bioassay data
# - add in intensity bioassay data using LD50 model, including possibility of
#    complete but imperfect resistance (trait fixation, not leading to 100%
#    survival in the test)
