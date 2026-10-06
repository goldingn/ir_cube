# Check that the greta model (build_dynamical_model(), with the closed-form op)
# and the plain-R predictions (R/dynamical_predictions.R) agree, at a random
# free state, for a given set of model options:
#   - predicted mortality at every training assay, and rho per type
#   - the log likelihood: the model's log density over all the data, less that
#     of the same model with one assay in the likelihood (which shares every
#     prior and Jacobian term), against the plain-R betabinomial log likelihood
#     of the other assays
# and print the log density itself, for regression checks between versions.
# With the species model (#47) the predictions are the mixture at each assay's
# arabiensis share, the map path is checked at the arabiensis fraction r(x),
# and the plain-R mixture is checked to reduce to one trajectory when the
# species do not differ. With the kdr covariate (#47), the model with its
# slopes at 0 is checked against the model without it, in greta (the log
# density) and in plain R (the predictions).
#
#   IR_CUBE_MODEL_OPTIONS='<options>' Rscript R/check_dynamical_model.R [seed] [sd]
# (the free state is N(0, sd^2), sd 0.5 by default; a smaller sd avoids states
# where p rounds to 1 at assays with survivors, and the log density is NaN)
# e.g.
#   IR_CUBE_MODEL_OPTIONS='dynamical_model_options(reversion = FALSE)' \
#     Rscript R/check_dynamical_model.R
#   IR_CUBE_MODEL_OPTIONS='dynamical_model_options(species = species_options())' \
#     Rscript R/check_dynamical_model.R
#   IR_CUBE_MODEL_OPTIONS='dynamical_model_options(kdr = kdr_options())' \
#     Rscript R/check_dynamical_model.R
#
# Run with the greta 0.6 environment (doc/cv_run_plan.md, section 1).

arguments <- commandArgs(trailingOnly = TRUE)
seed <- if (length(arguments) >= 1) as.integer(arguments[1]) else 1L
free_sd <- if (length(arguments) >= 2) as.numeric(arguments[2]) else 0.5

source("R/greta_setup.R")
start_greta(threads = 4)
source("R/dynamical_model.R")
suppressMessages({
  sink("/dev/null")
  source("R/validation_folds.R")
  source("R/validation_covariates.R")
  sink()
})
source("R/validation_functions.R")
source("R/dynamical_predictions.R")
source("R/two_stage_map_functions.R")

cat("options:", Sys.getenv("IR_CUBE_MODEL_OPTIONS", "defaults"), "\n")

build <- function(train_df, options = model_options) {
  build_dynamical_model(train_df = train_df,
                        df = df,
                        x_cell_years = x_cell_years,
                        cell_years_index = cell_years_index,
                        classes_index = classes_index,
                        types = types,
                        options = options,
                        x_cells_init = x_cells_init)
}

built <- build(df)
built_one <- build(df[1, ])
# the options as built, with the centre of the initial-state covariates
model_options <- built$options

log_density <- function(model, free) {
  f <- model$dag$generate_log_prob_function(which = "adjusted")
  free_tf <- tensorflow::tf$constant(matrix(free, nrow = 1),
                                     dtype = tensorflow::tf$float64)
  as.numeric(f(free_tf))
}

n_free <- length(unlist(built$model$dag$example_parameters(free = TRUE)))
stopifnot(n_free == length(unlist(
  built_one$model$dag$example_parameters(free = TRUE))))
set.seed(seed)
free <- rnorm(n_free, 0, free_sd)

# the free-state elements of a variable
free_columns <- function(model, name) {
  columns <- free_state_columns(model)
  columns[[attr(columns, "targets")[[name]]]]
}

# A random reversion rate is on the scale of its free state, around 1 per year,
# which drives p to 1 at most assays by the end of the series; so it is set to
# 0.05 per year (its free state is the log of the rate) and checked
if (identical(model_options$reversion, "estimated")) {
  columns <- free_columns(built$model, "reversion_rate")
  free[columns] <- log(0.05)
  trace <- built$model$dag$trace_values(matrix(free, nrow = 1))
  stopifnot(isTRUE(all.equal(
    unname(trace[1, grep("^reversion_rate", colnames(trace))]),
    rep(0.05, length(columns)))))
}

ld_all <- log_density(built$model, free)
ld_one <- log_density(built_one$model, free)

# the variables' values at this free state, as a one-draw fold for the plain-R
# path
trace <- built$model$dag$trace_values(matrix(free, nrow = 1))
stopifnot(identical(unname(trace),
                    unname(built_one$model$dag$trace_values(
                      matrix(free, nrow = 1)))))
fold <- list(draws = coda::mcmc.list(coda::mcmc(trace)),
             options = model_options, x_cells_init = x_cells_init)
parameters <- dynamical_parameter_draws(fold, classes_index, types,
                                        draw_index = 1)

# greta, at these values: calculate() on a one-draw greta_mcmc_list, as on a
# fitted model's draws
values <- greta:::as_greta_mcmc_list(
  coda::mcmc.list(coda::mcmc(trace)),
  list(raw_draws = coda::mcmc.list(coda::mcmc(matrix(free, nrow = 1))),
       model = built$model))
greta_values <- calculate(p = built$population_mortality_vec,
                          rho = built$terms$rho_types,
                          values = values)
greta_values <- as.matrix(greta_values)
p_greta <- greta_values[1, grep("^p\\[", colnames(greta_values))]
rho_greta <- greta_values[1, grep("^rho\\[", colnames(greta_values))]

# plain R
p_r <- plogis(c(dynamical_logit(parameters, df, df, x_cell_years,
                                cell_years_index)))
rho_r <- c(parameters$rho_types)

# the betabinomial log likelihood, parameterised as in betabinomial_p_rho() and
# not clamped (dbetabinom() clamps p away from 0 and 1)
a <- p_r * (1 / rho_r[df$type_id] - 1)
b <- a * (1 - p_r) / p_r
loglik_r <- extraDistr::dbbinom(df$died, df$mosquito_number, alpha = a,
                                beta = b, log = TRUE)

cat(sprintf("free parameters %d, log density %.10g\n", n_free, ld_all))
p_diff <- max(abs(p_greta - p_r))
logit_diff <- max(abs(qlogis(pmin(pmax(p_greta, 1e-12), 1 - 1e-12)) -
                        qlogis(pmin(pmax(p_r, 1e-12), 1 - 1e-12))))
rho_diff <- max(abs(rho_greta - rho_r))
loglik_diff <- ld_all - ld_one - sum(loglik_r[-1])
cat(sprintf("p: max abs diff %.3g, max logit diff %.3g (range %.3g-%.3g)\n",
            p_diff, logit_diff, min(p_r), max(p_r)))
cat(sprintf("rho: max abs diff %.3g\n", rho_diff))
cat(sprintf("log likelihood of assays 2..n: greta %.10g, plain R %.10g, diff %.3g\n",
            ld_all - ld_one, sum(loglik_r[-1]), loglik_diff))
# the tolerances allow for rounding: the logit of p near 1 and the sum of the
# log likelihood over ~27,000 assays lose the most precision
stopifnot(p_diff < 1e-10, logit_diff < 1e-6, rho_diff < 1e-12,
          abs(loglik_diff) < 1e-6)

# the map path (R/two_stage_map_functions.R), at the data cells in 2000, 2012
# and 2024, against dynamical_logit()
map_years <- c(2000, 2012, 2024)
map_rows <- df %>%
  distinct(cell, cell_id) %>%
  mutate(country_name = countries[
    built$lookups$cell_country_lookup[cell_id]])
covariates <- map_covariates(map_rows$cell, baseline_year, max(map_years),
                             model_options$selection_columns)

# the map path's covariates are x_cell_years at the data cells
n_fit_years <- max(cell_years_index$year_id)
x_map <- matrix(aperm(map_x(covariates, seq_len(nrow(map_rows)), n_fit_years),
                      c(2, 1, 3)),
                ncol = ncol(x_cell_years))
x_fit <- x_cell_years[match(paste(rep(map_rows$cell_id, each = n_fit_years),
                                  rep(seq_len(n_fit_years), nrow(map_rows))),
                            paste(cell_years_index$cell_id,
                                  cell_years_index$year_id)), ]
cat(sprintf("map covariates vs x_cell_years, %d cell-years x %d columns: max abs diff %.3g\n",
            nrow(x_fit), ncol(x_fit), max(abs(x_map - x_fit))))
stopifnot(identical(dim(x_map), dim(x_fit)), max(abs(x_map - x_fit)) == 0)

# the trends' product columns are 0 in the baseline year, unless a trend is
# given as a region x year matrix (e.g. g(1995) = 0.37)
trend_columns <- grep(":g_(dom|ag)$", colnames(x_cell_years))
if (length(trend_columns) > 0) {
  x_baseline <- x_cell_years[cell_years_index$year_id == 1, trend_columns]
  cat(sprintf("%d trend product columns in %d: %d non-zero values\n",
              length(trend_columns), baseline_year, sum(x_baseline != 0)))
  if (any(x_baseline != 0)) {
    design <- model_options$selection_columns
    if (!is.matrix(design$trend_pop) && !is.matrix(design$trend_crops)) {
      stop("trend product columns are not 0 in the baseline year")
    }
    warning("trend product columns are not 0 in the baseline year")
  }
}

logit_init_all <- map_logit_init(parameters, countries, regions, df)
cell_country_index <- match(map_rows$country_name,
                            dimnames(logit_init_all)[[2]])
# the share of the complex-wide predictions at each map cell (NULL without the
# species model), and for dynamical_logit(), at each row; and the standardised
# kdr at each map cell (NULL without the kdr covariate)
map_share <- prediction_share(model_options, map_rows$cell)
map_kdr <- prediction_kdr(model_options, map_rows$cell)
x_years <- map_x(covariates, seq_len(nrow(map_rows)),
                 max(map_years) - baseline_year + 1)
clamp <- function(l) pmin(pmax(l, qlogis(1e-12)), qlogis(1 - 1e-12))
map_difference <- 0
for (k in seq_along(types)) {
  dyn <- dynamical_logit_cells(parameters, k,
                               matrix(logit_init_all[, cell_country_index, k],
                                      1),
                               x_years, map_years - baseline_year + 1,
                               x_init = covariates$init, share = map_share,
                               kdr = map_kdr)
  for (y in map_years) {
    rows <- tibble(cell_id = map_rows$cell_id, type_id = k,
                   year_id = y - baseline_year + 1)
    l_rows <- c(dynamical_logit(parameters, rows, df, x_cell_years,
                                cell_years_index, share = map_share))
    l_map <- c(dyn[[as.character(y - baseline_year + 1)]])
    map_difference <- max(map_difference, abs(clamp(l_map) - clamp(l_rows)))
  }
}
cat(sprintf("map path vs plain R, logit, %d cells x %d types x %d years: max abs diff %.3g\n",
            nrow(map_rows), length(types), length(map_years), map_difference))
stopifnot(map_difference < 1e-9)

if (species_on(model_options)) {
  # identified assays (share 0 or 1) are predicted by one trajectory alone
  share <- arabiensis_share(df, model_options)
  cat(sprintf("species: %d assays of arabiensis, %d of other members, %d of the complex (share %.2f-%.2f)\n",
              sum(share == 1), sum(share == 0), sum(share > 0 & share < 1),
              min(share[share > 0 & share < 1]),
              max(share[share > 0 & share < 1])))
  l_assays <- c(dynamical_logit(parameters, df, df, x_cell_years,
                                cell_years_index))
  for (pure in c(0, 1)) {
    rows <- share == pure
    l_pure <- c(dynamical_logit(parameters, df[rows, ], df, x_cell_years,
                                cell_years_index, share = pure))
    stopifnot(max(abs(clamp(l_assays[rows]) - clamp(l_pure))) < 1e-12)
  }
  # with no difference between the species, the mixture is the other
  # members' trajectory at any share
  same <- parameters
  same$gamma_selection[] <- 0
  if (!is.null(same$gamma_cost)) same$gamma_cost[] <- 0
  same$arabiensis_floor <- same$other_floor
  # and no kdr effect, as the species' bands differ
  same$kdr_slopes <- lapply(same$kdr_slopes, function(x) 0 * x)
  rows <- df[seq(1, nrow(df), by = 10), ]
  l_mixed <- dynamical_logit(same, rows, df, x_cell_years, cell_years_index)
  l_other <- dynamical_logit(same, rows, df, x_cell_years, cell_years_index,
                             share = 0)
  species_difference <- max(abs(clamp(l_mixed) - clamp(l_other)))
  cat(sprintf("species: no difference between them, mixture vs one trajectory, logit: max abs diff %.3g\n",
              species_difference))
  stopifnot(species_difference < 1e-9)
}

if (kdr_on(model_options)) {
  cat(sprintf("kdr: standardised by logit kdr (complex) over %d cells, mean %.3f, sd %.3f\n",
              max(df$cell_id), model_options$kdr$centre,
              model_options$kdr$scale))
  # with its slopes at 0, the model is the model without the kdr covariate:
  # in greta, the log density differs by the slopes' priors at 0 only
  options_off <- model_options
  options_off$kdr <- FALSE
  built_off <- build(df, options_off)
  delta_names <- intersect(kdr_slope_names, names(built$variables))
  cat("kdr: slopes", toString(delta_names), "\n")
  free_zero <- free
  for (name in delta_names) {
    free_zero[free_columns(built$model, name)] <- 0
  }
  columns_off <- free_state_columns(built_off$model)
  free_off <- numeric(length(unlist(
    built_off$model$dag$example_parameters(free = TRUE))))
  for (name in names(attr(columns_off, "targets"))) {
    free_off[columns_off[[attr(columns_off, "targets")[[name]]]]] <-
      free_zero[free_columns(built$model, name)]
  }
  ld_difference <- log_density(built$model, free_zero) -
    log_density(built_off$model, free_off) -
    length(delta_names) * dnorm(0, log = TRUE)
  cat(sprintf("kdr: slopes 0 vs no kdr covariate, greta log density: diff %.3g\n",
              ld_difference))
  stopifnot(abs(ld_difference) < 1e-6)
  rm(built_off)
  # and in plain R, the predictions at every assay
  zero <- parameters
  zero$kdr_slopes <- lapply(zero$kdr_slopes, function(x) 0 * x)
  off <- zero
  off$kdr_slopes <- list()
  off$options <- options_off
  l_zero <- dynamical_logit(zero, df, df, x_cell_years, cell_years_index)
  l_off <- dynamical_logit(off, df, df, x_cell_years, cell_years_index)
  kdr_difference <- max(abs(clamp(l_zero) - clamp(l_off)))
  cat(sprintf("kdr: slopes 0 vs no kdr covariate, plain R, logit: max abs diff %.3g\n",
              kdr_difference))
  stopifnot(kdr_difference < 1e-9)
}
cat("all checks passed\n")
