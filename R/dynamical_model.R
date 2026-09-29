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

source("R/greta_setup.R")


# model options ------------------------------------------------------------

# Switches for model terms. The defaults are the current model; the others are
# placeholders for terms still to be added, and build_dynamical_model() refuses
# them until they are implemented:
#   rho               "class": one overdispersion per insecticide class;
#                     "type": one per insecticide type, types nested in
#                     class (#20)
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
  implemented <- list(rho = list("class", "type"))
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

  variables <- list(
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
    logit_init_mean = normal(qlogis(init$relative_prior), 1, dim = n_types)
  )

  rho <- switch(
    options$rho,
    # observation overdispersion per class. With rho = 0.067 the 95% interval
    # of a betabinomial at p = 0.5 is 0.5 wide, a range to treat as unlikely:
    # the half-normal sd of 0.025 puts P(rho < 0.067) at about 0.99 (see the
    # history of fit_model.R for the calculation)
    class = list(
      rho_classes = normal(0, 0.025, truncation = c(0, 1), dim = n_classes)
    ),
    # Per type, nested in class, on the logit scale and non-centred, as the
    # replicate-assay estimate in R/fig_illustrate_bioassay_variability.R:
    #   logit rho_type = rho_mu + rho_sigma_class z_class + rho_sigma_type z_type
    # with the same priors. The prior centre is that of the replicate-assay
    # rho (0.15); rho here also absorbs misfit of the model, and the class-level
    # rho above reached posterior means of 0.24-0.38 despite its half-normal
    # prior, so no stronger prior is put on it.
    type = list(
      rho_mu = normal(qlogis(0.15), 1),
      rho_sigma_class = normal(0, 0.5, truncation = c(0, Inf)),
      rho_sigma_type = normal(0, 0.5, truncation = c(0, Inf)),
      rho_class_raw = normal(0, 1, dim = n_classes),
      rho_type_raw = normal(0, 1, dim = n_types)
    )
  )

  c(variables, rho)
}

# Deterministic transforms from the variables `v` to the quantities the dynamics
# need:
#   beta_type           n_covs x n_types, log selection effect of each covariate
#   logit_init_country  n_countries x n_types, logit of the initial fraction
#                       susceptible q_0 in each country
#   rho_types           n_types, the observation overdispersion of each type
#                       (with rho = "class", its class's)
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
    v <- lapply(v, function(x) {
      if (length(dim(x)) <= 1 || (is.matrix(x) && ncol(x) == 1)) c(x) else x
    })
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

  # observation overdispersion per type
  rho_types <- switch(
    options$rho,
    class = v$rho_classes[classes_index],
    type = {
      logit_rho_class <- v$rho_mu + v$rho_sigma_class * v$rho_class_raw
      inv_logit(logit_rho_class[classes_index] +
                  v$rho_sigma_type * v$rho_type_raw)
    }
  )

  list(beta_type = beta_type,
       logit_init_country = logit_init_country,
       rho_types = rho_types,
       # intermediate quantities, which the figure scripts read
       beta_class = beta_class,
       init_region_effect = init_region_effect,
       init_country_effect = init_country_effect,
       logit_init_relative = logit_init_relative)
}


# Initial values for the model's `variables` from a cached set
# (temporary/inits.RDS, posterior means from an earlier fit): the cached values
# of the variables the model has, and, for rho = "type", a start for rho_mu
# from the cached class-level rho. Variables with neither start where greta
# puts them.
dynamical_inits <- function(cached, variables) {
  cached <- unclass(cached)
  out <- cached[intersect(names(cached), names(variables))]
  if ("rho_mu" %in% names(variables) && !is.null(cached$rho_classes)) {
    out$rho_mu <- mean(qlogis(c(cached$rho_classes)))
  }
  do.call(greta::initials, out)
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
# derived terms, and two functions for predictions after sampling:
# mortality(rows) returns the greta array of predicted bioassay mortality at the
# (cell_id, type_id, year_id) of `rows`, and all_states() the states at every
# cell, type and year.
build_dynamical_model <- function(train_df,
                                  df,
                                  x_cell_years,
                                  cell_years_index,
                                  classes_index,
                                  types,
                                  options = dynamical_model_options()) {

  check_dynamical_model_options(options)
  check_greta_fill()

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

  # predicted mortality (the fraction susceptible) at the (cell_id, type_id,
  # year_id) of `rows`, computing the states only for the cell-type pairs there
  mortality <- function(rows) {
    pairs <- distinct(tibble(cell_id = as.integer(rows$cell_id),
                             type_id = as.integer(rows$type_id)))
    states <- closed_form_states(terms, x_cell_years, pairs$cell_id,
                                 pairs$type_id, lookups$cell_country_lookup,
                                 n_times)
    pair_index <- match(paste(rows$cell_id, rows$type_id),
                        paste(pairs$cell_id, pairs$type_id))
    states[cbind(pair_index, rows$year_id)]
  }

  # the states at every cell, type and year, as n_unique_cells x n_types x
  # n_times, for predictions after sampling. (Creating this before model() would
  # put it in the model's graph.)
  all_states <- function() {
    pairs <- expand.grid(cell_id = seq_len(n_unique_cells),
                         type_id = seq_len(n_types))
    states <- closed_form_states(terms, x_cell_years, pairs$cell_id,
                                 pairs$type_id, lookups$cell_country_lookup,
                                 n_times)
    dim(states) <- c(n_unique_cells, n_types, n_times)
    states
  }

  # likelihood
  population_mortality_vec <- mortality(train_df)
  rho <- terms$rho_types[train_df$type_id]
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
       population_mortality_vec = population_mortality_vec,
       mortality = mortality,
       all_states = all_states,
       lookups = lookups,
       options = options)
}


# the selection recursion ---------------------------------------------------

# Haploid selection: with q the fraction susceptible and w the relative fitness
# of the resistant phenotype,
#   q_t = q_{t-1} / (q_{t-1} + (1 - q_{t-1}) w_t),
# which divides the odds of susceptibility by w_t each year, so it is exactly
# additive on the logit scale:
#   logit q_t = logit q_0 - sum_{s <= t} log w_s,  w_s = 1 + x_s' exp(beta).
# The state recorded for year t has had the fitness of years 1..t applied (as
# greta.dynamics recorded it). This is computed by one greta op (#25).
#
# log w is computed from the linear predictor as
#   m + log(exp(-m) + x' exp(beta - m)),  m = max(0, max_k beta_k)
# per type, with the gradient stopped through m. That is exact for any m, and
# cannot overflow in float64. The covariates are all non-negative.

# The TensorFlow side. Arguments are tensors with a leading batch dimension B
# (greta's), then the constants:
#   beta_type            (B, n_covs, n_types)
#   logit_init_country   (B, n_countries, n_types)
#   x_pairs              (J, n_times, n_covs) covariates of each pair's cell
#   pair_type            (J) 0-based type of each pair
#   pair_init            (J) 0-based index of each pair's (country, type) in
#                        the row-major flattened logit_init_country
# Returns (B, J, n_times): the fraction susceptible for each pair and year.
tf_closed_form_states <- function(beta_type, logit_init_country,
                                  x_pairs, pair_type, pair_init) {
  tf <- tensorflow::tf
  dtype <- beta_type$dtype
  # reshaped to vectors, since reticulate passes a length-one R vector as a
  # scalar
  pair_type <- tf$reshape(tf$constant(pair_type, dtype = tf$int32), list(-1L))
  pair_init <- tf$reshape(tf$constant(pair_init, dtype = tf$int32), list(-1L))
  x_pairs <- tf$constant(x_pairs, dtype = dtype)

  # log fitness, (B, J, n_times)
  m <- tf$stop_gradient(tf$maximum(tf$reduce_max(beta_type, axis = 1L,
                                                  keepdims = TRUE),
                                   tf$constant(0, dtype = dtype)))
  effect <- tf$gather(tf$exp(beta_type - m), pair_type, axis = 2L)
  m <- tf$transpose(tf$gather(m, pair_type, axis = 2L), c(0L, 2L, 1L))
  selection <- tf$einsum("jtp,bpj->bjt", x_pairs, effect)
  log_w <- m + tf$math$log(tf$exp(-m) + selection)

  # initial state of each pair, (B, J, 1)
  n_init <- as.integer(prod(dim(logit_init_country)[-1]))
  logit_q0 <- tf$gather(tf$reshape(logit_init_country, c(-1L, n_init)),
                        pair_init, axis = 1L)

  tf$sigmoid(tf$expand_dims(logit_q0, 2L) - tf$cumsum(log_w, axis = 2L))
}

# The greta side: the fraction susceptible for the cell-type pairs
# (pair_cell, pair_type), as a J x n_times greta array. `terms` is the output of
# dynamical_terms(); `x_cell_years` has one row per (cell, year), cell-major;
# each cell takes the initial state of country cell_country_lookup[cell].
closed_form_states <- function(terms, x_cell_years, pair_cell, pair_type,
                               cell_country_lookup, n_times) {
  n_covs <- ncol(x_cell_years)
  n_countries <- nrow(terms$logit_init_country)
  stopifnot(nrow(x_cell_years) %% n_times == 0,
            length(pair_cell) == length(pair_type))

  # (cells, years, covariates) and each pair's slice of it
  x_cells <- aperm(array(x_cell_years,
                         c(n_times, nrow(x_cell_years) / n_times, n_covs)),
                   c(2, 1, 3))
  x_pairs <- x_cells[pair_cell, , , drop = FALSE]
  pair_country <- cell_country_lookup[pair_cell]
  stopifnot(!anyNA(pair_country))

  # the TensorFlow function is found in this small environment, which is saved
  # with the node, so a reloaded draws object can still calculate() through it
  op_env <- new.env(parent = globalenv())
  op_env$tf_closed_form_states <- tf_closed_form_states

  greta:::op("closed_form_states",
             terms$beta_type,
             terms$logit_init_country,
             operation_args = list(
               x_pairs = x_pairs,
               pair_type = as.integer(pair_type - 1),
               pair_init = as.integer((pair_country - 1) *
                                        ncol(terms$logit_init_country) +
                                        pair_type - 1)),
             tf_operation = "tf_closed_form_states",
             tf_function_env = op_env,
             dim = c(length(pair_cell), n_times))
}
