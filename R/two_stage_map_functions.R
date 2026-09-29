# Helpers for mapping the two-stage model (#21) on the full prediction grid:
# the dynamical model's posterior draws at every mask cell, and the stage-A
# correction's latent fields projected from the mesh nodes to those cells.
#
# The dynamical part is in plain R on the logit scale (see
# R/dynamical_predictions.R for why the recursion has a closed form there), so
# that it can be paired draw for draw with the correction's cut-posterior mode
# shift; R/predict.R uses it for the dynamical model's maps. Functions only:
# source from the repo root after R/dynamical_predictions.R and
# R/two_stage_correction.R.


# covariates at the prediction cells ------------------------------------------

# All covariates the dynamical model uses, at the cells `cells` of the mask and
# for the years baseline_year..end_year, built the way predict.R builds
# x_cell_years_predict: the cubes are padded back to the baseline year and
# forward to end_year by repeating their first and last layers.
#
# `design` is the fit's selection design (options$selection_columns, see
# selection_design() in R/model_covariates.R).
#
# Returned as a cells x years x n array of the time-varying covariates (nets,
# irs, pop and any hinge columns, the column order of x_cell_years), a
# cells x 10 matrix of the
# static crop covariates, and a cells x 2 matrix `init` of the initial-state
# covariates (init_covariate_matrix(), R/model_covariates.R), rather than
# predict.R's long (cell, year) matrix: the
# long form is 53M rows x 13 columns (5.5 GB) for 1.48M cells and 36 years,
# whereas this keeps the crops once per cell and lets a chunk of cells be
# assembled one year at a time.
map_covariates <- function(cells, baseline_year = 1995, end_year = 2030,
                           design = selection_design()) {
  list(time_varying = selection_time_varying(cells, baseline_year, end_year,
                                             design),
       flat = selection_flat(cells),
       init = init_covariate_matrix(cells))
}


# initial conditions for every country ----------------------------------------

# Draws of logit q_0, the initial fraction susceptible, for every African
# country in the UNSD lookup (not only those with data), as in predict_batch():
# countries and regions without data get fresh N(0, 1) raw deviations, i.e. the
# hierarchical model's prediction for a new country or region. Returns a
# draws x countries x types array with dimnames[[2]] the country names.
#
# The country -> region mapping for prediction is the UNSD one, as in
# predict.R, where the fit took each country's region from its first record;
# the two agree for every observed country (checked in two_stage_maps.R).
#
# With initial-state covariates (options$init_covariates, #19) the initial
# state varies within a country, so this returns instead the logit relative
# initial state of each country (the same draws x countries x types array),
# with the covariates' coefficients (draws x n_init_covs x types) as attribute
# "init_coef"; map_cell_logit_init() takes either and gives logit q_0 at cells.
#
# `options` are the fit's model options: fold_options(fold), or the full fit's
# model_options.
#
# `logit_init_mean` is needed only for fits whose draws do not name it (see
# logit_init_mean_draws()); classes_index is needed by dynamical_terms(), which
# also computes the selection effects.
map_logit_init <- function(draws_matrix,
                           logit_init_mean,
                           types,
                           classes_index,
                           countries,
                           regions,
                           lookup = country_region_lookup(),
                           seed = 1,
                           options) {

  n_draws <- nrow(draws_matrix)
  n_types <- length(types)

  all_countries <- unique(lookup$country_name)
  all_regions <- unique(lookup$region)

  variables <- variable_draws(draws_matrix, logit_init_mean)
  # the fit's options (fold_options()), which say whether the draws include
  # initial-state covariates
  stopifnot(identical(!is.null(variables$init_coef),
                      !is.null(options$init_covariates)))
  init_region_raw <- variables$init_region_raw
  init_country_raw <- variables$init_country_raw

  # observed countries and regions keep their sampled deviations; the rest are
  # drawn from the prior, once per posterior draw and shared by every cell
  set.seed(seed)
  country_raw <- array(rnorm(n_draws * length(all_countries) * n_types),
                       c(n_draws, length(all_countries), n_types))
  region_raw <- array(rnorm(n_draws * length(all_regions) * n_types),
                      c(n_draws, length(all_regions), n_types))
  country_raw[, match(countries, all_countries), ] <- init_country_raw
  region_raw[, match(regions, all_regions), ] <- init_region_raw
  variables$init_country_raw <- country_raw
  variables$init_region_raw <- region_raw

  country_region <- match(lookup$region[match(all_countries,
                                              lookup$country_name)],
                          all_regions)

  # the model's own transform (R/dynamical_model.R), over all countries
  covariates <- !is.null(options$init_covariates)
  logit_init <- dynamical_terms_draws(
    variables,
    classes_index = classes_index,
    country_region_index = country_region,
    types = types,
    terms = if (covariates) "logit_init_relative" else "logit_init_country",
    options = options)[[1]]
  dimnames(logit_init) <- list(NULL, all_countries, types)
  if (covariates) {
    attr(logit_init, "init_coef") <- array(
      variables$init_coef,
      c(n_draws, length(options$init_covariates), n_types),
      dimnames = list(NULL, options$init_covariates, types))
  }
  logit_init
}

# Logit q_0 at cells, cells x draws, for type k, from map_logit_init()'s
# output: `cell_country` indexes each cell's country in its second dimension,
# and `init` is the cells' initial-state covariates (map_covariates()$init),
# used only when the fit has them.
map_cell_logit_init <- function(logit_init, cell_country, k, init = NULL,
                                types = dimnames(logit_init)[[3]]) {
  n_draws <- dim(logit_init)[1]
  l <- t(matrix(logit_init[, cell_country, k], nrow = n_draws))
  init_coef <- attr(logit_init, "init_coef")
  if (is.null(init_coef)) {
    return(l)
  }
  stopifnot(!is.null(init))
  x <- init[, dimnames(init_coef)[[2]], drop = FALSE]
  l <- l + x %*% t(matrix(init_coef[, , k], nrow = n_draws))
  logit_init_from_relative(l, init_frac_constants(types)$min[k])
}

# map_logit_init()'s output for the draws `draws`, keeping the initial-state
# coefficients with them
subset_logit_init <- function(logit_init, draws) {
  out <- logit_init[draws, , , drop = FALSE]
  init_coef <- attr(logit_init, "init_coef")
  if (!is.null(init_coef)) {
    attr(out, "init_coef") <- init_coef[draws, , , drop = FALSE]
  }
  out
}


# dynamical draws on a chunk of cells ------------------------------------------

# Logit of predicted bioassay mortality at a chunk of cells, for draws of one
# insecticide type, returned as a list over `years_keep` of cells x draws
# matrices.
#
#   effect       draws x n_covs: exp(beta_type) for this type
#   logit_init   cells x draws: logit q_0 at each cell's country
#   time_varying cells x years x 3 and flat cells x n_flat, from
#                map_covariates(), for this chunk
#   floor        draws: the floor on mortality (parameters$mortality_floor
#                from dynamical_parameter_draws()), or NULL for none
#   kappa        draws: the reversion kappa of this type
#                (parameters$kappa_type[, k]), or NULL for none
#
# logit q_t = logit q_0 - sum_{s <= t} log w_s, with the fitness of year 1
# (the baseline year) already applied to the year-1 state, exactly as
# dynamical_predictions() and greta.dynamics do; with reversion, less t kappa
# (reversion_kappa()), t = 1 in the first of `years`, which must be the
# baseline year. Mortality is q_t, or with a floor f, f + (1 - f) q_t.
dynamical_logit_chunk <- function(effect, logit_init, time_varying, flat,
                                  years, years_keep, floor = NULL,
                                  kappa = NULL) {
  stopifnot(all(years_keep %in% years))
  last <- max(match(years_keep, years))
  effect_t <- t(effect)
  cumulative <- 0
  out <- list()
  for (t in seq_len(last)) {
    x_t <- cbind(time_varying[, t, ], flat)
    cumulative <- cumulative + log1p(x_t %*% effect_t)
    if (!is.null(kappa)) {
      cumulative <- cumulative + rep(kappa, each = nrow(logit_init))
    }
    if (years[t] %in% years_keep) {
      out[[as.character(years[t])]] <- floored_logit_mortality(
        logit_init - cumulative,
        if (!is.null(floor)) rep(floor, each = nrow(logit_init)))
    }
  }
  out
}

# dynamical_logit_chunk() for type k at rows `rows` of map_covariates()'s
# output `covariates`, whose cells are in countries `country_index` (into
# dimnames(logit_init)[[2]]), with the fit's floor and reversion. `effect`
# (draws x n_covs x types) and `floor` and `kappa_type` (from
# dynamical_parameter_draws(), NULL when the fit has none) are for the same
# draws as `logit_init` (map_logit_init()).
map_type_logit <- function(k, rows, country_index, effect, logit_init,
                           covariates, years, years_keep, floor = NULL,
                           kappa_type = NULL) {
  dynamical_logit_chunk(
    effect = matrix(effect[, , k], nrow = dim(effect)[1]),
    logit_init = map_cell_logit_init(logit_init, country_index, k,
                                     covariates$init[rows, , drop = FALSE]),
    time_varying = covariates$time_varying[rows, , , drop = FALSE],
    flat = covariates$flat[rows, , drop = FALSE],
    years = years, years_keep = years_keep,
    floor = floor,
    kappa = if (!is.null(kappa_type)) kappa_type[, k])
}


# correction fields at the mesh nodes ----------------------------------------

# Node values of the correction omega + xi for each year in `years`, as
#   omega      n_nodes_omega x n_draws (the same in every year)
#   xi[[y]]    n_nodes_xi x n_draws
# from a matrix of latent vectors `theta` (n_latent x n_draws, in the fit's
# latent order). Years at or before t0 have xi = 0, years in (t0, T] read the
# fitted x, and later years run the AR(1) for eta forward from
# eta_T = x_T - psi x_{T-1}, adding `innovations[[h]]` (n_nodes_xi x n_draws,
# N(0, Q_eta^-1)) at horizon h and accumulating into xi as
# xi_{T+h} = psi xi_{T+h-1} + eta_{T+h}, as predict_correction() does (psi = 1,
# the undamped sum, unless the fit has damped_xi). With innovations = NULL the
# recursion carries the mean forward instead: eta_{T+h} = phi^h eta_T, so
# undamped
#   xi_{T+h} = xi_T + eta_T phi (1 - phi^h) / (1 - phi),
# the plateauing forecast of the issue; damped, the mean decays to 0.
correction_node_fields <- function(fit, theta, years, innovations = NULL) {
  theta <- as.matrix(theta)
  n_nodes_xi <- fit$mesh_xi$n
  omega <- theta[fit$blocks$w_omega, , drop = FALSE]
  out <- list(omega = omega, xi = list())
  if (fit$variant != "omega_xi_u") {
    for (y in years) {
      out$xi[[as.character(y)]] <- matrix(0, n_nodes_xi, ncol(theta))
    }
    return(out)
  }

  x <- theta[fit$blocks$x, , drop = FALSE]
  x_col <- function(j) x[(j - 1) * n_nodes_xi + seq_len(n_nodes_xi), ,
                         drop = FALSE]
  xi_T <- x_col(fit$n_years)
  psi <- correction_psi(fit)
  eta_T <- if (fit$n_years > 1) xi_T - psi * x_col(fit$n_years - 1) else xi_T

  phi <- fit$hyper$phi
  for (y in years) {
    if (y <= fit$t0) {
      xi_y <- matrix(0, n_nodes_xi, ncol(theta))
    } else if (y <= fit$T) {
      xi_y <- x_col(y - fit$t0)
    } else {
      xi_y <- xi_T
      eta <- eta_T
      for (h in seq_len(y - fit$T)) {
        eta <- phi * eta
        if (!is.null(innovations)) {
          eta <- eta + sqrt(1 - phi ^ 2) * innovations[[h]]
        }
        xi_y <- psi * xi_y + eta
      }
    }
    out$xi[[as.character(y)]] <- xi_y
  }
  out
}

# row standard deviations of a matrix
row_sds <- function(x) {
  n <- ncol(x)
  mu <- rowMeans(x)
  sqrt(pmax(rowSums(x ^ 2) - n * mu ^ 2, 0) / (n - 1))
}
