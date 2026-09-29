# The dynamical model, defined once.
#
# build_dynamical_model() builds the greta model from the data, the covariates
# and a training set. fit_model.R calls it with the full data, fit_fold()
# (R/fit_validation_fold.R) with a cross-validation training fold; the only
# difference between the two is which records enter the likelihood.
#
# The priors and the deterministic transforms from the sampled parameters to the
# selection coefficients and the initial states are in dynamical_variables() and
# dynamical_terms(). dynamical_terms() is written to work both on greta arrays
# and on one posterior draw as plain R arrays, so that the plain-R predictions
# in R/dynamical_predictions.R and R/two_stage_map_functions.R use the same
# definitions. A new model term goes in those two functions (and, if it changes
# the recursion, in the state computation), switched on through
# dynamical_model_options().
#
# Source from the repo root, after R/packages.R and R/functions.R.


# model options ------------------------------------------------------------

# Switches for model terms. The defaults are the current model; the others are
# placeholders for terms still to be added, and build_dynamical_model() refuses
# them until they are implemented:
#   rho               "class": one overdispersion per insecticide class;
#                     "type": one per insecticide type (#20)
#   mortality_floor   an estimated floor on bioassay mortality (#14)
#   init_covariates   names of static covariates for the initial state (#19)
#   selection_columns extra columns of the selection design matrix (#23)
#   reversion         a per-class reversion rate to susceptibility (#24)
dynamical_model_options <- function(rho = c("class", "type"),
                                    mortality_floor = FALSE,
                                    init_covariates = NULL,
                                    selection_columns = NULL,
                                    reversion = FALSE) {
  list(rho = match.arg(rho),
       mortality_floor = mortality_floor,
       init_covariates = init_covariates,
       selection_columns = selection_columns,
       reversion = reversion)
}

check_dynamical_model_options <- function(options) {
  defaults <- dynamical_model_options()
  stopifnot(setequal(names(options), names(defaults)))
  implemented <- list(rho = "class")
  for (name in names(defaults)) {
    allowed <- if (name %in% names(implemented)) implemented[[name]] else
      list(defaults[[name]])
    if (!any(vapply(allowed, identical, logical(1), options[[name]]))) {
      stop("dynamical model option '", name, "' is not implemented yet")
    }
  }
  invisible(options)
}


# fixed quantities ---------------------------------------------------------

# Prior and minimum values for the initial fraction susceptible, per type. The
# initial state is modelled on the logit of its relative position between
# init_frac_min and 1. More flexibility for DDT, less for the others.
init_frac_constants <- function(types) {
  prior <- ifelse(types == "DDT", 0.9, 0.95)
  min <- ifelse(types == "DDT", 0.75, 0.9)
  list(prior = prior,
       min = min,
       # mean logit proportion of the distance from the minimum to 1
       relative_prior = (prior - min) / (1 - min))
}

# Lookups from the full data (never a fold, so that every fold and the full fit
# index the same countries and cells): the region of each country, and the
# country whose initial state each cell takes, i.e. the country of the cell's
# first record.
dynamical_lookups <- function(df) {
  country_region_index <- df %>%
    group_by(country_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(country_id, region_id) %>%
    arrange(country_id) %>%
    pull(region_id)

  cell_country_lookup <- df %>%
    group_by(cell_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(cell_id, country_id) %>%
    arrange(cell_id) %>%
    pull(country_id)

  list(country_region_index = country_region_index,
       cell_country_lookup = cell_country_lookup)
}


# parameters -----------------------------------------------------------------

# The model's variables (the free parameters), with their priors, as a named
# list of greta arrays. Every element is passed to model(), so every one has
# named columns in the draws.
dynamical_variables <- function(n_covs, n_classes, n_types, n_regions,
                                n_countries, types,
                                options = dynamical_model_options()) {

  init <- init_frac_constants(types)

  list(
    # initial fractions susceptible: a prior logit-mean per type, and IID
    # deviations by region and by country within region
    init_region_sd = normal(0, 1, truncation = c(0, Inf), dim = n_types),
    init_country_sd = normal(0, 1, truncation = c(0, Inf), dim = n_types),
    init_region_raw = normal(0, 1, dim = c(n_regions, n_types)),
    init_country_raw = normal(0, 1, dim = c(n_countries, n_types)),
    # hierarchical regression coefficients: overall -> class -> type
    beta_overall = normal(0, 1, dim = n_covs),
    beta_class_raw = normal(0, 1, dim = c(n_covs, n_classes)),
    beta_type_raw = normal(0, 1, dim = c(n_covs, n_types)),
    sigma_overall = normal(0, 1, dim = n_covs, truncation = c(0, Inf)),
    sigma_class = normal(0, 1, dim = n_covs, truncation = c(0, Inf)),
    # observation overdispersion. With rho = 0.067 the 95% interval of a
    # betabinomial at p = 0.5 is 0.5 wide, a range to treat as unlikely: the
    # half-normal sd of 0.025 puts P(rho < 0.067) at about 0.99 (see the
    # history of fit_model.R for the calculation)
    rho_classes = normal(0, 0.025, truncation = c(0, 1), dim = n_classes),
    logit_init_mean = normal(qlogis(init$relative_prior), 1, dim = n_types)
  )
}

# Deterministic transforms from the variables `v` to the quantities the dynamics
# need:
#   beta_type           n_covs x n_types, log selection effect of each covariate
#   logit_init_country  n_countries x n_types, logit of the initial fraction
#                       susceptible q_0 in each country
# and some intermediate quantities.
#
# `v` is a named list of either greta arrays or plain R arrays for a single
# posterior draw (dimensions as in dynamical_variables(), vectors as vectors or
# one-column matrices). `country_region_index` maps the rows of
# v$init_country_raw to the rows of v$init_region_raw; for prediction at
# countries without data, pass raw deviations and an index for all countries.
dynamical_terms <- function(v, classes_index, country_region_index, types,
                            options = dynamical_model_options()) {

  is_greta <- inherits(v$beta_overall, "greta_array")
  if (is_greta) {
    inv_logit <- greta::ilogit
    softplus <- greta::log1pe
  } else {
    v <- lapply(v, function(x) if (is.matrix(x) && ncol(x) == 1) c(x) else x)
    inv_logit <- stats::plogis
    softplus <- function(x) -stats::plogis(-x, log.p = TRUE)
  }

  # selection effects: doubly hierarchical
  beta_class_sigma <- sweep(v$beta_class_raw, 1, v$sigma_overall, FUN = "*")
  beta_class <- sweep(beta_class_sigma, 1, v$beta_overall, FUN = "+")
  beta_type_sigma <- sweep(v$beta_type_raw, 1, v$sigma_class, FUN = "*")
  beta_type <- beta_class[, classes_index] + beta_type_sigma

  # initial state: logit relative position above init_frac_min, the prior mean
  # plus region and country deviations
  init_region_effect <- sweep(v$init_region_raw, 2, v$init_region_sd,
                              FUN = "*")
  init_country_effect <- sweep(v$init_country_raw, 2, v$init_country_sd,
                               FUN = "*")
  init_country_overall_effect <- init_country_effect +
    init_region_effect[country_region_index, ]
  logit_init_relative <- sweep(init_country_overall_effect, 2,
                               v$logit_init_mean, FUN = "+")

  # q_0 = min + range * ilogit(l). Its logit is computed without forming q_0,
  # which rounds to 1 in double precision when l is large: with
  # 1 - q_0 = range * ilogit(-l),
  #   logit q_0 = log(min + range * ilogit(l)) - log(range) + softplus(l)
  init <- init_frac_constants(types)
  init_range <- 1 - init$min
  q0_scaled <- sweep(sweep(inv_logit(logit_init_relative), 2, init_range,
                           FUN = "*"),
                     2, init$min, FUN = "+")
  logit_init_country <- sweep(log(q0_scaled) + softplus(logit_init_relative),
                              2, log(init_range), FUN = "-")

  list(beta_type = beta_type,
       logit_init_country = logit_init_country,
       # intermediate quantities, which the figure scripts read
       beta_class = beta_class,
       init_region_effect = init_region_effect,
       init_country_effect = init_country_effect,
       logit_init_relative = logit_init_relative)
}


# the model ------------------------------------------------------------------

# Build the greta model with the likelihood over `train_df`.
#
#   train_df          the records in the likelihood (the full data, or a fold)
#   df                the full data, for the dimensions and lookups
#   x_cell_years      covariates, one row per (cell_id, year_id), cell-major
#   cell_years_index  data frame of the cell_id and year_id of each row
#   classes_index     class of each type
#   types             type names, in type_id order
#
# Returns a list with the model, its variables (as passed to model()), the
# derived terms, and a function mortality(rows) that returns the greta array of
# predicted bioassay mortality at the (cell_id, type_id, year_id) of `rows`,
# for predictions after sampling.
build_dynamical_model <- function(train_df,
                                  df,
                                  x_cell_years,
                                  cell_years_index,
                                  classes_index,
                                  types,
                                  options = dynamical_model_options()) {

  check_dynamical_model_options(options)

  n_covs <- ncol(x_cell_years)
  n_unique_cells <- max(df$cell_id)
  n_times <- max(cell_years_index$year_id)
  n_classes <- max(df$class_id)
  n_types <- length(types)
  n_regions <- max(df$region_id)
  n_countries <- max(df$country_id)
  stopifnot(
    length(classes_index) == n_types,
    # the layout the state computation relies on
    identical(as.integer(cell_years_index$cell_id),
              rep(seq_len(n_unique_cells), each = n_times)),
    identical(as.integer(cell_years_index$year_id),
              rep(seq_len(n_times), n_unique_cells))
  )

  lookups <- dynamical_lookups(df)

  variables <- dynamical_variables(n_covs = n_covs,
                                   n_classes = n_classes,
                                   n_types = n_types,
                                   n_regions = n_regions,
                                   n_countries = n_countries,
                                   types = types,
                                   options = options)
  terms <- dynamical_terms(variables,
                           classes_index = classes_index,
                           country_region_index = lookups$country_region_index,
                           types = types,
                           options = options)

  # fraction susceptible for every cell, type and year
  effect_type <- exp(terms$beta_type)
  fitness_cell_years <- 1 + x_cell_years %*% effect_type
  fitness_array <- fitness_cell_years
  dim(fitness_array) <- c(n_times, n_unique_cells, n_types, 1)
  init_array <- ilogit(terms$logit_init_country)[lookups$cell_country_lookup, ]
  dim(init_array) <- c(dim(init_array), 1)
  dynamic_cells <- iterate_dynamic_function(
    transition_function = haploid_next,
    initial_state = init_array,
    niter = n_times,
    w = fitness_array,
    parameter_is_time_varying = c("w"),
    tol = 0)

  mortality <- function(rows) {
    dynamic_cells$all_states[cbind(rows$cell_id, rows$type_id, rows$year_id)]
  }

  # likelihood
  population_mortality_vec <- mortality(train_df)
  rho <- variables$rho_classes[train_df$class_id]
  distribution(train_df$died) <- betabinomial_p_rho(
    N = train_df$mosquito_number,
    p = population_mortality_vec,
    rho = rho)

  # model() takes the target names from its call, so it is called with the
  # variables' names as symbols
  target_env <- list2env(variables, parent = environment())
  m <- do.call("model", lapply(names(variables), as.name), envir = target_env)

  list(model = m,
       variables = variables,
       terms = terms,
       effect_type = effect_type,
       population_mortality_vec = population_mortality_vec,
       dynamic_cells = dynamic_cells,
       mortality = mortality,
       lookups = lookups,
       options = options)
}

# haploid selection: state = fraction with the susceptible allele, w = relative
# fitness of the resistant phenotype
haploid_next <- function(state, iter, w) {
  q <- state
  p <- 1 - q
  q / (q + p * w)
}
