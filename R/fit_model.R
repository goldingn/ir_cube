# fit model

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load covariate rasters

# load time-varying net coverage data
nets_cube <- rast("data/clean/nets_per_capita_cube.tif")
# and scale it to 0-1
nets_cube <- nets_cube / max(global(nets_cube, "max", na.rm = TRUE))
# names(nets_cube) <- str_remove(names(nets_cube), "nets_")

# load the other layers
covs_flat <- rast("data/clean/flat_covariates.tif")

mask <- covs_flat[[1]] * 0

# load bioassay data
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

# extract spatiotemporal covariates from the cube
all_extract <- nets_cube %>%
  terra::extract(unique_cells) %>%
  mutate(
    cell = unique_cells,
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

system.time(
  draws <- mcmc(m,
                chains = n_chains)
)
# user   system  elapsed 
# 7937.520 3529.599 2452.889 

save.image(file = "temporary/fitted_model.RData")
# load(file = "temporary/fitted_model.RData")

# # check convergence
# coda::gelman.diag(draws,
#                   autoburnin = FALSE,
#                   multivariate = FALSE)

## maybe use posterior means of parameters as initial values for next run
# post <- calculate(
#   beta_overall,
#   sigma_overall,
#   sigma_class,
#   beta_class_raw,
#   beta_type_raw,
#   init_fraction_susceptible,
#   logit_rho,
#   values = draws,
#   nsim = 200)
# 
# mn <- function(x) {
#   post_mean <- apply(x, 2:3, mean)
#   if (dim(post_mean)[2] == 1L) {
#     post_mean <- as.vector(post_mean)
#   }
#   post_mean
# }
# 
# post_means <- lapply(post, mn)
# dput(post_means)
# list(beta_overall = c(-1.17145325567021, -2.2902407283197, -1.00374949286669, 
# -1.54847363451857, -1.96474786969832, -1.6178971952705, -2.51359109963177, 
# -2.35087324569348), sigma_overall = c(0.803771981253376, 1.73761854799904, 
# 0.98176653843001, 0.76101902276866, 0.802156224541578, 0.864402469338693, 
# 1.42718957862772, 1.58054956414425), sigma_class = c(0.499882541903096, 
# 0.826596637439702, 0.877049418464544, 1.41693332386631, 1.58518716682321, 
# 1.80891954718765, 0.980687825405729, 0.895005111387372), beta_class_raw = structure(c(0.65743114714174, 
# -1.24157041805265, -0.104361626744828, -0.161663897145048, -0.699437503327882, 
# 0.0247457379826926, -1.07241866200016, -0.933592949375868, -1.33553986374692, 
# -0.945675116001378, -0.99830388401624, -0.252445605955901, -0.112743920589157, 
# -0.688213652138625, -0.770964578682854, -0.88696679728572, 0.132439902288804, 
# -1.00917605117471, 0.382457593127378, 0.128333222621881, -0.449974337571745, 
# 0.138581665072777, -0.960546192516594, -1.05621335957094, -0.198492781299611, 
# -0.961879614212879, -0.845102972984205, -0.691606905556476, -0.159552818300963, 
# -0.281987716453989, -0.658713122320559, -1.00453579263954), dim = c(8L, 
# 4L), dimnames = list(NULL, NULL)), beta_type_raw = structure(c(-0.293551434359401, 
# -0.432301311513869, 0.23973240983031, -0.330909249408167, -0.753805123295738, 
# -1.24706790723403, -0.337376966215799, -0.588581307360449, 1.14123724785774, 
# 0.100031824263388, -0.669743679900074, -0.258000076723278, -0.403330007889042, 
# -1.09751295538928, -0.0116158694144284, 0.306652129455743, -0.90296210767768, 
# 0.0821310284163792, -0.47001614110545, 0.695253602399654, -0.367393984857844, 
# -0.828831285232909, -0.00661011153301158, 0.0950519810073573, 
# -0.23813864457677, -0.649513512994148, 0.316999899288933, 0.54172335466494, 
# -1.23479422070961, 0.571330380825494, -0.38806326640497, -0.387134581399039, 
# -0.297858542786957, -0.1788729283402, -0.799750794294775, 0.514708609985091, 
# -0.311882237986767, 0.00714882736505584, -0.444013595095626, 
# -0.336351487881085, 0.453238153465782, -0.233822151071852, -0.464331694980122, 
# -0.922455876868687, -1.02893703043709, -0.464346837511513, 0.333455192372367, 
# -0.345442155939028, 0.0940335157543912, -0.616469286505075, -1.19817437334276, 
# -1.23562065651018, -0.0882569406928293, -1.26673004720811, -0.22778053669358, 
# -0.298395461281566, -0.35846663959714, -0.252411212846841, 0.142870756639969, 
# -0.677313995525739, 1.4395912928536, -0.468254229621208, -0.61415095198615, 
# 0.0906779679026175, -0.133810202326697, -0.0906130110448876, 
# 0.988292014155293, -0.683954089863114, -0.188389868208588, 1.47103725481245, 
# 0.0341052621797561, -0.368862511979362), dim = 8:9, dimnames = list(
#     NULL, NULL)), init_fraction_susceptible = c(0.876032397175313, 
# 0.923619493444688, 0.975419765723507, 0.838770187761344, 0.857395527137113, 
# 0.953978449782556, 0.935922245295876, 0.960353580692299, 0.910550105709105
# ), logit_rho = -0.299000653322215)
# 
# # find some locations with lots of years of data and plot predictions and data
# # for these
# locations_plot <- df %>%
#   group_by(latitude, longitude) %>%
#   filter(
#     n_distinct(year_start) >= 8
#   ) %>%
#   select(
#     country_name,
#     latitude,
#     longitude,
#     cell_id
#   ) %>%
#   distinct()
# 
# # # geocode then find more interpretable names (google to see if this is how
# # # they are referred to in IR papers)
# # locations_plot_geocoded <- locations_plot %>%
# #   reverse_geocode(
# #     lat = latitude,
# #     long = longitude,
# #     method = 'osm',
# #     full_results = TRUE
# #   ) %>%
# #   mutate(
# #     place = case_when(
# #       address == "Soumousso, Houet, Hauts-Bassins, Burkina Faso" ~ "Soumousso, Burkina Faso",
# #       address == "Katito-Kendu Bay-Homa Bay road, Oriang, Central ward, Karachuonyo, Homa Bay, Nyanza, Kenya" ~ "Homa Bay, Kenya",
# #       address == "RNIE 7, Kandi 3, Kandi, Alibori, Bénin" ~ "Kandi, Benin",
# #       address == "Aménagement de périmètre maraîcher, Djougou, Donga, Bénin" ~ "Djougou, Benin"
# #     )
# #   ) %>%
# #   select(
# #     place,
# #     country_name,
# #     latitude, longitude
# #   )
# 
# locations_plot <- locations_plot %>%
#   mutate(
#     place = case_when(
#       nearish(latitude, 11.011) & nearish(longitude, -4.056) ~ "Soumousso, Burkina Faso",
#       nearish(latitude, -0.392) & nearish(longitude, 34.626) ~ "Homa Bay, Kenya",
#       nearish(latitude, 11.134) & nearish(longitude, 2.938) ~ "Kandi, Benin",
#       nearish(latitude, 9.705) & nearish(longitude, 1.667) ~ "Djougou, Benin"
#     )
#   )
# 
# # pull these out for plotting
# preds_plot_setup <- expand_grid(
#   cell_id = locations_plot$cell_id,
#   year_start = baseline_year:max(df$year_start),
#   insecticide_type = types
# ) %>%
#   mutate(
#     year_id = year_start - baseline_year + 1,
#     type_id = match(insecticide_type, types)
#   )
# 
# index_plot <- preds_plot_setup %>%
#   select(cell_id,
#          type_id,
#          year_id) %>%
#   as.matrix()
# 
# fraction_susceptible_plot <- dynamic_cells$all_states[index_plot]
# population_mortality_plot <- fraction_susceptible_plot
# 
# # simulate mortalities under binomial sampling
# sample_size_plot <- 100
# binomial_died <- binomial(sample_size_plot, population_mortality_plot)
# binomial_mortality <- binomial_died / sample_size_plot
# 
# # simulate mortalities under betabinomial sampling
# betabinomial_died <- betabinomial_p_rho(N = sample_size_plot,
#                                               p = population_mortality_plot,
#                                               rho = rho)
# betabinomial_mortality <- betabinomial_died / sample_size_plot
# 
# sims <- calculate(population_mortality_plot,
#                   binomial_mortality,
#                   betabinomial_mortality,
#                   values = draws,
#                   nsim = 2000)
# 
# # posterior mean mortality rate, and intervals for the population value
# # (posterior uncertainty), and posterior predictive intervals (posterior
# # uncertainty and sampling error) for observed mortality from binomial sampling,
# # and from betabinomial sampling
# pop_mort_mean <- colMeans(sims$population_mortality_plot[, , 1])
# pop_mort_ci <- apply(sims$population_mortality_plot[, , 1],
#                      2,
#                      quantile,
#                      c(0.025, 0.975))
# binomial_mort_ci <- apply(sims$binomial_mortality[, , 1],
#                           2,
#                           quantile,
#                           c(0.025, 0.975))
# betabinomial_mort_ci <- apply(sims$betabinomial_mortality[, , 1],
#                               2,
#                               quantile,
#                               c(0.025, 0.975))
# 
# insecticides_plot <- c("Deltamethrin", "Permethrin", "Bendiocarb")
# 
# preds_plot <- preds_plot_setup %>%
#   mutate(
#     mortality = pop_mort_mean,
#     pop_lower = pop_mort_ci[1, ],
#     pop_upper = pop_mort_ci[2, ],
#     binomial_lower = binomial_mort_ci[1, ],
#     binomial_upper = binomial_mort_ci[2, ],
#     betabinomial_lower = betabinomial_mort_ci[1, ],
#     betabinomial_upper = betabinomial_mort_ci[2, ],
#   ) %>%
#   filter(
#     insecticide_type %in% insecticides_plot
#   ) %>%
#   left_join(
#     locations_plot,
#     by = "cell_id"
#   )
# 
# points_plot <- df %>%
#   filter(
#     cell_id %in% locations_plot$cell_id,
#     insecticide_type %in% insecticides_plot
#   ) %>%
#   mutate(
#     mortality = died / mosquito_number
#   ) %>%
#   left_join(
#     locations_plot,
#     by = "cell_id"
#   )
#   
# 
# # plot these, then add cell data over the top
# preds_plot %>%
#   ggplot(
#     aes(
#       x = year_start,
#       y = mortality,
#       group = place
#     )
#   ) +
#   # betabinomial posterior predictive intervals on bioassay data (captures
#   # sample size and non-independence effect)
#   geom_ribbon(
#     aes(
#       ymax = betabinomial_lower,
#       ymin = betabinomial_upper,
#     ),
#     colour = grey(0.4),
#     linewidth = 0.25,
#     linetype = 2,
#     fill = grey(0.9)
#   ) +
#   # binomial posterior predictive intervals on bioassay data (captures sample
#   # size but assumes independence, so underestimates variance)
#   geom_ribbon(
#     aes(
#       ymax = binomial_lower,
#       ymin = binomial_upper,
#     ),
#     fill = grey(0.6)
#   ) +
#   # credible intervals on population-level proportion (our best guess at the
#   # 'truth')
#   geom_ribbon(
#     aes(
#       ymax = pop_lower,
#       ymin = pop_upper,
#     ),
#     fill = grey(0.4)
#   ) +
#   scale_y_continuous(labels = scales::percent) +
#   facet_grid(place ~ insecticide_type) +
#   coord_cartesian(ylim = c(0, 1)) +
#   theme_minimal() +
#   geom_point(
#     aes(
#       x = year_start,
#       y = mortality,
#       group = place
#     ),
#     data = points_plot,
#   )
# 
# ggsave("figures/fit_subset.png",
#        bg = "white")

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
batch_bigness <- 3e5
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

# loop through insecticides
for (this_insecticide in types_save) {
  
  # make a raster cube to stick predictions into
  this_ir_cube <- template_cube
  
  # grab the current seed, so we can do calculate in batches but not shuffle
  # parameters over space
  this_seed <- greta::.internals$utils$misc$get_seed()
  
  # loop through batches of cells
  for (batch_index in seq_len(n_batches)) {
    
    print(sprintf("processing batch %i of %i for %s",
                  batch_index,
                  n_batches,
                  this_insecticide))
    
    # index to this insecticide
    type_id_predict <- match(this_insecticide, types)
    
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
    # set.seed(this_seed)
    pred_batch <- calculate(population_mortality_batch,
                            values = draws,
                            seed = this_seed,
                            nsim = 50)[[1]]
    
    # compute posterior mean over the draws, for cells and years
    batch_pred_mean <- apply(pred_batch, 2:4, mean)
    
    # stick this batch of predictions in the raster cube, one year at a time
    for (y in seq_along(years_predict)) {
      this_ir_cube[[y]][cell_batch] <- batch_pred_mean[, 1, y]
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
# - save fitted model and run raster predictions separately
# - maybe add a resistance cost parameter
#    (positive, as fitness = 1 + selection - cost)
# - set up code for posterior predictive checking
# - set up code for out of sample (future timesteps, spatial blocks) evaluation
#    of model fit.
# - add in more discriminating bioassay data
# - add in intensity bioassay data using LD50 model, including possibility of
#    complete but imperfect resistance (trait fixation, not leading to 100%
#    survival in the test)
