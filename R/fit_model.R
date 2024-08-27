# fit model

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# set the start of the timeseries considered in modelling (the start of
# non-negligible levels of resistance) - assume it's before the mass-rollout of
# nets
baseline_year <- 1995

# load the mask
mask <- rast("data/clean/raster_mask.tif")

# load time-varying rasters

# net use
nets_cube <- rast("data/clean/net_use_cube.tif")

# IRS coverage
irs_cube <- rast("data/clean/irs_coverage_scaled_cube.tif")

# human population
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

# load bioassay data
ir_africa <- readRDS(file = "data/clean/all_gambiae_complex_data.RDS")

df <- ir_africa %>%
  group_by(insecticide_type) %>%
  # subset to the most common concentration for each insecticide
  filter(
   concentration == sample_mode(concentration)
  ) %>%
  ungroup() %>%
  filter(
    # drop any from before the baseline
    year_start >= baseline_year
  ) %>%
  mutate(
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

# pull covariates for all unique cells
unique_cells <- unique(df$cell)
n_unique_cells <- length(unique_cells)

# create design matrix at all unique cells and for all years

# pull out temporally-static covariates for all cells
flat_extract <- covs_flat %>%
  extract(unique_cells) %>%
  mutate(
    cell = unique_cells,
    .before = everything()
  )

# TODO: extract other spatiotemporal covariates here

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
n_obs <- nrow(df)
n_covs <- ncol(x_cell_years)
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

# multiply through to get relative fitness of resistance for each insecticide
# type at each cell

# convert beta to positive effect sizes, them matrix multiply with indicator
# variables to constrain to positive effects of the covariates on selection
# coefficient/relative fitness
effect_type <- exp(beta_type)

# # compute selection coefficients for cell-years, convert to relative fitness,
# and extract at data locations
selection_cell_years <- x_cell_years %*% effect_type
fitness_cell_years <- 1 + selection_cell_years

# reformat this in to a 3D array with dimensions:
#   n_times x n_unique_cells x n_types x 1
# to solve dynamics with time-varying fitness (time must be first, then other
# two must match state variable, which has a trailing dimension of size 1)
fitness_array <- fitness_cell_years
dim(fitness_array) <- c(n_times, n_unique_cells, n_types, 1)

# # check this is the right orientation!
# which_cell_id <- 3
# which_type_id <- 4
# array_subset <- fitness_array[, which_cell_id, which_type_id, ]
# cell_idx <- cell_years_index$cell_id == which_cell_id
# matrix_subset <- fitness_cell_years[cell_idx, which_type_id]
# sims <- calculate(
#   array_subset,
#   matrix_subset,
#   nsim = 1)
# identical(sims$array_subset[1, , 1, 1, 1],
#           sims$matrix_subset[1, , 1])

# model fractions susceptible across the continent prior to baseline. More flex
# for DDT and less for others
init_frac_sd <- ifelse(types == "DDT", 0.01, 0.001)

init_fraction_susceptible <- normal(1, init_frac_sd,
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
init_array <- sweep(ones(n_unique_cells, n_types),
                    2,
                    init_fraction_susceptible,
                    FUN = "*")
dim(init_array) <- c(dim(init_array), 1)

# iterate through time to get fraction susceptible for all years at all cells
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
# normal(0, 1.6) is similar to logit distribution; with hierarcichal prior, set
# overall sd to: sqrt((1.6 ^ 2) - 1)
logit_rho_mean <- normal(0, 1.3)
logit_rho_sd <- normal(0, 1, truncation = c(0, Inf))
logit_rho_raw <- normal(0, 1, dim = n_classes)
logit_rho_classes <- logit_rho_mean + logit_rho_raw * logit_rho_sd
rho_classes <- ilogit(logit_rho_classes)

distribution(df$died) <- betabinomial_p_rho(N = df$mosquito_number,
                                            p = population_mortality_vec,
                                            rho = rho_classes[df$class_id])

m <- model(
  beta_overall,
  sigma_overall,
  beta_type_raw,
  init_fraction_susceptible,
  logit_rho_mean,
  logit_rho_raw
)

# set the inits for obs_var_multiplier to be large (more params will fit OK to
# data to start with, then chains can move towards better parts of parameter
# space)
n_chains <- 4

system.time(
  draws <- mcmc(m,
                chains = n_chains)
)

# new data, from 1995
# user    system   elapsed 
# 22645.760  9301.635  6053.067 

# new data, from 2000
# user    system   elapsed 
# 15255.782  5870.483  3916.685 

# old data, from 2000
# user   system  elapsed
# 7937.520 3529.599 2452.889

save.image(file = "temporary/fitted_model.RData")

# # check convergence
# coda::gelman.diag(draws,
#                   autoburnin = FALSE,
#                   multivariate = FALSE)

# run prediction code across scenario net cubes (bash batching to prevent
# restarts and memory leaks?)

# things to do next:
# - mask plotted outputs by DVS extent mask
# - implement dose-response model and add in intensity bioassay data
# - consider possibility (identifiability) of complete but imperfect resistance
# (trait fixation, not leading to 100% survival in the test)
# - output validation stratified by species
# - visualise and code up functions for validation method
# - run predictions across multiple scenario cubes
# - maybe extend nets to pre-2000 to find a start point with flat resistance
# - maybe use diploid model (hierarchical logit-het.-dominance parameter) to
#    capture temporal variability in development of resistance phenotype
# - compute Africa-average LLIN-effective susceptibility to show along with
#    bioassay data plot. Also plot weighted LLIN usage at these locations, for
#    context.
# - maybe use a very sparse hierarchical GP to model initial susceptibility
# - maybe add a resistance cost parameter
#    (positive, as fitness = 1 + selection - cost)
# - set up code for posterior predictive checking
# - set up code for out of sample (future timesteps, spatial blocks) evaluation
#    of model fit.
