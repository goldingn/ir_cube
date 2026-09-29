# fit model

# load packages and functions
source("R/packages.R")
source("R/functions.R")
source("R/dynamical_model.R")

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

# build the model, with the likelihood over all the data (R/dynamical_model.R)
built <- build_dynamical_model(train_df = df,
                               df = df,
                               x_cell_years = x_cell_years,
                               cell_years_index = cell_years_index,
                               classes_index = classes_index,
                               types = types)
m <- built$model

# the variables and derived quantities, as named objects in the saved image,
# which the figure and prediction scripts read
list2env(built$variables, globalenv())
list2env(built$terms, globalenv())
effect_type <- exp(beta_type)
# the states at every cell, type and year, as cells x types x years (created
# after model(), so they are not computed while sampling)
dynamic_cells <- list(all_states = built$all_states())
population_mortality_vec <- built$population_mortality_vec
country_region_index <- built$lookups$country_region_index
cell_country_lookup <- built$lookups$cell_country_lookup
init_frac_min <- init_frac_constants(types)$min
init_range <- 1 - init_frac_min

n_chains <- 8

# used cached posterior means as inits
inits_one <- readRDS("temporary/inits.RDS")
inits <- replicate(n_chains,
                   inits_one,
                   simplify = FALSE)

Lmax <- 30
Lmin <- round(Lmax / 2)
system.time(
  draws <- mcmc(m,
                chains = n_chains,
                initial_values = inits,
                warmup = 2000,
                sampler = hmc(Lmin = Lmin, Lmax = Lmax),
                n_samples = 1000)
)

# user    system   elapsed 
# 15596.230  6544.659  4148.049 

# check convergence
rhats <- coda::gelman.diag(draws,
                           autoburnin = FALSE,
                           multivariate = FALSE)
summary(rhats$psrf)

# save fitted model to use for plotting and predictions
save.image(file = "temporary/fitted_model.RData")


# add more samples and update saved model
draws <- extra_samples(draws, 1000)
rhats <- coda::gelman.diag(draws,
                           autoburnin = FALSE,
                           multivariate = FALSE)
summary(rhats$psrf)

# save fitted model to use for plotting and predictions
save.image(file = "temporary/fitted_model.RData")




# save posterior means as initial values for a future model run

# these have to match the arguments to model(), and need to be greta variable
# nodes (not operation nodes)
posts <- do.call(calculate,
                 c(built$variables, list(values = draws, nsim = 100)))

post_means <- lapply(posts,
                     function(x) {
                       apply(x, 2:3, mean)
                     })
inits <- do.call(greta::initials, post_means)
saveRDS(inits, "temporary/inits.RDS")


