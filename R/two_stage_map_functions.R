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
