# fit model

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# set the start of the timeseries considered in modelling (the start of
# non-negligible levels of resistance) - assume it's before the mass-rollout of
# nets
baseline_year <- 1995

# set the final year of data (insufficient and spatially biased data for 2025)
final_data_year <- 2024

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

# if necessary, pad forward to the final data year, repeating the last
nets_cube <- post_pad_cube(nets_cube, final_data_year)
irs_cube <- post_pad_cube(irs_cube, final_data_year)
pop_cube <- post_pad_cube(pop_cube, final_data_year)

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

# numbers of unique location/time records per insecticide
record_counts <- ir_africa %>%
  group_by(insecticide_type, latitude, longitude, year_start) %>%
  summarise(
    mosquito_number = sum(mosquito_number),
    .groups = "drop"
  ) %>%
  group_by(insecticide_type) %>%
  summarise(
    n = n(),
    mosquito_number = mean(mosquito_number),
    .groups = "drop"
  ) %>%
  arrange(desc(n))

# keep only a subset of insecticide types

# keep the first 9 of these: those with at least 1000 unique places/times, and
# alpha-cypermethrin (914 unique) because of its use in LLINs.
insecticides_keep <- c("Alpha-cypermethrin",
                       "Deltamethrin",
                       "Lambda-cyhalothrin", 
                       "Permethrin",
                       "Fenitrothion",
                       "Malathion",
                       "Pirimiphos-methyl",
                       "DDT",
                       "Bendiocarb")

# # note: only one study has chlorfenapyr resistance (Benin in 2022)
# ir_africa %>% filter(insecticide_type == "Chlorfenapyr") %>% View()

df <- ir_africa %>%
  filter(insecticide_type %in% insecticides_keep) %>%
  group_by(insecticide_type) %>%
  # subset to the most common concentration for each insecticide
  filter(
   concentration == sample_mode(concentration)
  ) %>%
  ungroup() %>%
  filter(
    # drop any from before the baseline
    year_start >= baseline_year,
    year_start <= final_data_year
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
  )

# create indices to categorical vectors
classes <- unique(df$insecticide_class)
types <- unique(df$insecticide_type)
regions <- unique(df$region)
countries <- unique(df$country_name)
unique_cells <- unique(df$cell)
years <- baseline_year - 1 + sort(unique(df$year_id))

df <- df %>%
  mutate(
    cell_id = match(cell, unique_cells),
    region_id = match(region, regions),
    country_id = match(country_name, countries),
    class_id = match(insecticide_class, classes),
    type_id = match(insecticide_type, types)
  )

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


# use an IID hierarcichal model for the initial fraction susceptible, with a
# prior logit-mean, and variance across regions and countries within regions

# Note: we could make this ICAR for the same computational complexity and possibly
# fewer parameters (we can drop the regional level)
# https://github.com/goldingn/greta.car

# Maybe do that if it's real speckly, and consider adding subnational levels
# too.

# prior assumption of the fractions susceptible across the continent at the
# baseline. More flex for DDT and less for others

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

# hist(calculate(init_array[1, 5, 1], nsim = 1e4)[[1]])

# compute the fraction expressing the susceptible phenotype over time against
# each insecticide with a single-locus population genetics model. Either with
# haploid trait evolution (approximately equivalent to diploid evolution with
# heterozygote dominance = 0.5), or with a diploid model with variable
# heterozygote dominance

# solve through time to get fraction susceptible over time

# haploid version:
# state = current proportion with *susceptible* allele
# iter = iteration number (required by greta.dynamics)
# w = selection pressure (relative fitness) for resistance phenotype
haploid_next <- function(state, iter, w) {
  # fraction susceptible
  q <- state
  # fraction resistant
  p <- 1 - q
  # number resistant grows relative to resistant by ratio w, re-normalise to get
  # fraction remaining susceptible
  q / (q + p * w)
}


# diploid version: state = current frequency of *susceptible* alleles in the
# population iter = iteration number (required by greta.dynamics) w = selection
# pressure (relative fitness) for resistant phenotype h = level of heterozygote
# dominance (0.5 = heterozygotes have intermediate resistance between homozygote
# susceptible and resistant individuals)

# Ref here for a nice example notation and discussion re. mossies:
# https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0001387
diploid_next <- function(state, iter, w, h = 0.5) {
  # frequency of the susceptible allele
  q <- state
  # frequency of the resistant allele
  p <- 1 - q
  # pre-compute things to save a few FLOPs
  s <- w - 1
  pq <- p * q
  p2 <- p ^ 2
  
  # compute the frequency of the resistant allele in the next generation, given
  # the fractions of the population homozygotic resistant (p^2), homozygotic
  # susceptible (q^2) and heterozygotic (pq), and their relative fitnesses (hom.
  # susceptible = 1, hom. resistant = 1+s = w, het. = 1 + hs)
  
  # size of population with resistant allele in next gen, relative to previous gen
  numerator <- p2 * w + pq * (1 + h * s)
  # size of total population in next gen, relative to previous gen, to normalise
  denominator <- 1 + s * (p2 + 2 * h * pq)
  # new proportion with resistant phenotype now
  new_p <- numerator / denominator
  # return proportion susceptible
  new_q <- 1 - new_p
  new_q
}

# # model hierarchical heterozygote dominance
# logit_het_dom_mean <- normal(0, 1)
# logit_het_dom_sd <- normal(0, 1, truncation = c(0, Inf))
# logit_het_dom_raw <- normal(0, 1, dim = n_types)
# logit_het_dom <- logit_het_dom_mean + logit_het_dom_raw * logit_het_dom_sd
# het_dom <- ilogit(logit_het_dom)

# het_dom_array <- sweep(ones(n_unique_cells, n_types),
#                        2,
#                        het_dom,
#                        FUN = "*")
# dim(het_dom_array) <- c(dim(het_dom_array), 1)

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
index <- cbind(df$cell_id, df$type_id, df$year_id)
fraction_susceptible_vec <- dynamic_cells$all_states[index]

# resistant phenotype expression rate = expected mortality in bioassays
population_mortality_vec <- fraction_susceptible_vec

# define observation model

# model overdispersion in the data via an overdispersion parameter rho.

# at p = 0.5 (the highest variance point), with rho = 0.067 the 95% interval of
# a betabinomial is 0.5; a very large error range we wish to treat as unlikely,
# and penalise towards smaller values. Found this by trial and error with
# simulations:

# p <- 0.5
# rho <- 0.067
# size <- 1000
# a <- p * (1 / rho - 1)
# b <- a * (1 - p) / p
# upper <- qbbinom(0.975, size = size, alpha = a, beta = b, nsims = 1e7)
# ci_width <- 2 * ((upper / size) - p)
# ci_width

# So set this ('rho_max' = 0.067) as being at the upper end of the likely
# scenario, say p(rho < rho_max) = 0.99. Compute the corresponding standard
# deviation of a half-normal (ignoring the truncation in calculations as it has
# negligible effect at these ranges

# # just maths:
# # q = 2 * pnorm(rho_max, 0, sd)
# # q / 2 = pnorm(rho_max, 0, sd)
# # qnorm(q / 2, 0, sd) = rho_max
# # qnorm(q / 2, 0, 1) * sd = rho_max
# # sd = rho_max / qnorm(q, 0, 1)

# calculation
# rho_max = 0.067
# q = 0.99
# sd = rho_max / qnorm(q, 0, 1)
# sd
# [1] 0.02880051
# so sd on rho prior to 0.025

rho_classes <- normal(0, 0.025,
                      truncation = c(0, 1),
                      dim = n_classes)

distribution(df$died) <- betabinomial_p_rho(N = df$mosquito_number,
                                            p = population_mortality_vec,
                                            rho = rho_classes[df$class_id])

m <- model(
  # initial fractions susceptible
  init_region_sd,
  init_country_sd,
  init_region_raw,
  init_country_raw,
  # hierarchical regression coefficients
  beta_overall,
  beta_class_raw,
  beta_type_raw,
  sigma_overall,
  sigma_class,
  # dispersion parameters
  rho_classes
)


n_chains <- 4

# used cached posterior means as inits
inits_one <- readRDS("temporary/inits.RDS")
inits <- replicate(n_chains,
                   inits_one,
                   simplify = FALSE)

system.time(
  draws <- mcmc(m,
                chains = n_chains,
                initial_values = inits,
                warmup = 2000,
                n_samples = 1000)
)

# user    system   elapsed 
# 15596.230  6544.659  4148.049 

# check convergence
coda::gelman.diag(draws,
                  autoburnin = FALSE,
                  multivariate = FALSE)

# save fitted model to use for plotting and predictions
save.image(file = "temporary/fitted_model.RData")




# save posterior means as initial values for a future model run

# these have to match the arguments to model(), above, and neet to be greta
# variable nodes (not operation nodes)
posts <- calculate(
  # initial fractions susceptible
  init_region_sd,
  init_country_sd,
  init_region_raw,
  init_country_raw,
  # hierarchical regression coefficients
  beta_overall,
  beta_class_raw,
  beta_type_raw,
  sigma_overall,
  sigma_class,
  # dispersion parameters
  rho_classes,
  values = draws,
  nsim = 100
)

post_means <- lapply(posts,
                     function(x) {
                       apply(x, 2:3, mean)
                     })
inits <- do.call(greta::initials, post_means)
saveRDS(inits, "temporary/inits.RDS")


