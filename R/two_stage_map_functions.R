# Helpers for running the dynamical model on the full prediction grid, for the
# two-stage maps (R/two_stage_maps.R), R/predict.R and the figure scripts. The
# recursion and the parameters are those of R/dynamical_predictions.R; this
# adds the grid's covariates and the initial states of countries without data.
# Functions only: source from the repo root after R/dynamical_predictions.R.

# All covariates the dynamical model uses, at mask cells `cells` for the years
# baseline_year..end_year, for the fit's selection design `design`
# (options$selection_columns; R/model_covariates.R):
#   time_varying  cells x years x n, the time-varying columns, in the column
#                 order of x_cell_years (selection_time_varying())
#   flat          cells x n_flat, the static crop columns (selection_static())
#   init          cells x 2, the initial-state covariates
# rather than a long (cell, year) matrix (5.5 GB for the 1.48M cells and 36
# years); map_x() assembles a chunk of cells in the column order of
# x_cell_years
map_covariates <- function(cells, baseline_year = 1995, end_year = 2030,
                           design = selection_design()) {
  list(time_varying = selection_time_varying(cells, baseline_year, end_year,
                                             design),
       flat = selection_static(cells, design),
       init = init_covariate_matrix(cells, design))
}

# The covariates of rows `rows` of map_covariates() for its first n_years
# years, as the cells x years x n_covs array dynamical_logit_cells() takes
map_x <- function(covariates, rows, n_years) {
  flat <- covariates$flat[rows, , drop = FALSE]
  n_time_varying <- dim(covariates$time_varying)[3]
  x <- array(NA_real_, c(length(rows), n_years, n_time_varying + ncol(flat)))
  for (t in seq_len(n_years)) {
    x[, t, ] <- cbind(matrix(covariates$time_varying[rows, t, ],
                             nrow = length(rows)), flat)
  }
  x
}

# Draws of the logit relative initial state (as
# parameters$logit_init_relative, see dynamical_parameter_draws()) for every
# African country in the UNSD lookup, not only those with data, as a
# draws x countries x types array with the country names as dimnames[[2]].
# Countries with data keep their fitted states. Countries and regions without
# data get fresh N(0, 1) deviations (the hierarchical model's prediction for a
# new one), one per posterior draw and type, drawn with the caller's RNG and
# shared by every cell of the country: a new country's state is its region's
# level plus init_country_sd times its deviation, and a new region's level
# logit_init_mean plus init_region_sd times its deviation. The region of a country
# is the UNSD one, as in predict.R; the fit took each country's region from
# its first record, and the two must agree for every country with data
map_logit_init <- function(parameters, countries, regions, df,
                           lookup = country_region_lookup()) {
  all_countries <- unique(lookup$country_name)
  unsd_region <- lookup$region[match(all_countries, lookup$country_name)]

  fit_region <- df %>%
    distinct(country_id, .keep_all = TRUE) %>%
    arrange(country_id)
  stopifnot(identical(fit_region$country_id, seq_along(countries)),
            identical(regions[fit_region$region_id],
                      unsd_region[match(countries, all_countries)]))

  v <- parameters$variables
  n_draws <- parameters$n_draws
  n_types <- length(parameters$types)
  stopifnot(dim(v$init_region_raw)[2] == length(regions))
  fresh <- function(n) array(rnorm(n_draws * n * n_types),
                             c(n_draws, n, n_types))

  country <- match(all_countries, countries)
  fitted <- !is.na(country)
  new_countries <- which(!fitted)
  country_z <- fresh(length(new_countries))
  new_regions <- setdiff(unique(unsd_region), regions)
  region_raw <- array(NA_real_, c(n_draws, length(regions) +
                                    length(new_regions), n_types))
  region_raw[, seq_along(regions), ] <- v$init_region_raw
  region_raw[, length(regions) + seq_along(new_regions), ] <-
    fresh(length(new_regions))
  new_region <- match(unsd_region[new_countries], c(regions, new_regions))

  logit_init <- array(NA_real_, c(n_draws, length(all_countries), n_types),
                      dimnames = list(NULL, all_countries, parameters$types))
  logit_init[, fitted, ] <- parameters$logit_init_relative[, country[fitted], ,
                                                           drop = FALSE]
  per_type <- function(name) matrix(v[[name]], n_draws)
  for (k in seq_len(n_types)) {
    region_level <- matrix(region_raw[, new_region, k], n_draws) *
      per_type("init_region_sd")[, k] + per_type("logit_init_mean")[, k]
    logit_init[, new_countries, k] <- region_level +
      matrix(country_z[, , k], n_draws) * per_type("init_country_sd")[, k]
  }
  logit_init
}

# column standard deviations of a matrix, e.g. over the draws
col_sds <- function(x) {
  n <- nrow(x)
  mu <- colMeans(x)
  sqrt(pmax(colSums(x ^ 2) - n * mu ^ 2, 0) / (n - 1))
}

# n* (#40): the number of independent bioassays of N mosquitoes at a
# pixel-year whose mean would estimate population mortality as precisely as
# the prediction. Unlike the posterior SD, it is scaled by a bioassay's noise,
# which also shrinks towards mortality 0 or 1.
#
# A bioassay of type k has variance p_k (1 - p_k) D_k / N about the population
# mortality p_k, with D_k = 1 + (N - 1) rho_k, rho_k the bioassay overdispersion
# from replicate pairs (rho_lookup()). The prediction is q = sum_k w_k p_k (one
# type, or llin_effective). n bioassays split by weight, n w_k of type k, have
# a weighted mean with error variance sum_k w_k p_k (1 - p_k) D_k / (n N).
# Averaging over the posterior (E[p_k (1 - p_k)] = m_k (1 - m_k) - s_k^2) and
# equating it to the posterior variance of q, s^2:
#   n* = sum_k w_k [m_k (1 - m_k) - s_k^2] D_k / (N s^2),
# for one type [m (1 - m) / s^2 - 1] D / N. Averaging rho over the types
# instead would use q (1 - q) >= sum_k w_k p_k (1 - p_k), overstating n*.
#
# s = sum_k w_k s_k: the types' errors perfectly correlated, which gives the
# largest s the marginal SDs allow, and so the smallest n*. The corrections are
# fitted per type, so in the posterior the types' errors are nearly
# independent. At held-out pixel-years with two pyrethroids the map's errors
# correlate at 0.99 (95% CI 0.94-1.03) once the bioassay noise is removed
# (PR #59), but that includes u, p and survey error shared by the
# pyrethroids, which n*_pop leaves out. So s rests on being the largest, not
# on the 0.99.
#
# means, sds: lists of the types' posterior means and SDs (numbers, matrices
# or rasters), of m + omega + xi for n*_pop or with u and p for n*_local;
# rho, weights: one per type
n_star <- function(means, sds, rho, weights = 1, N = 100) {
  noise <- Map(function(m, s_k, rho_k, w) {
    w * (m * (1 - m) - s_k ^ 2) * (1 + (N - 1) * rho_k)
  }, means, sds, rho, weights)
  s <- Reduce(`+`, Map(`*`, sds, weights))
  Reduce(`+`, noise) / (N * s ^ 2)
}
