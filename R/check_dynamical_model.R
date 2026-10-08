# Check that the greta model (build_dynamical_model(), with the closed-form op)
# and the plain-R predictions (R/dynamical_predictions.R) agree, at a random
# free state, for a given set of model options:
#   - predicted mortality at every training assay, and rho per type
#   - the log likelihood: the model's log density over all the data, less that
#     of the same model with one assay in the likelihood (which shares every
#     prior and Jacobian term), against the plain-R betabinomial log likelihood
#     of the other assays (with likelihood = "weighted_binomial", the weighted
#     binomial log likelihood at the replicate rho, from the plain-R logit)
# and print the log density itself, for regression checks between versions.
# Whatever the options, the weighted binomial log density of greta
# (weighted_binomial(), R/weighted_binomial.R) is also checked against plain
# R at a few points, from p near 0 to p near 1, with and without a floor and
# as a mixture of two trajectories.
# With the species model (#47) the predictions are the mixture at each assay's
# arabiensis share, the map path is checked at the arabiensis fraction r(x),
# and the plain-R mixture is checked to reduce to one trajectory when the
# species do not differ. With the kdr covariate (#47), the model with its
# slopes at 0 is checked against the model without it, in greta (the log
# density) and in plain R (the predictions), and with the kdr-dependent floor,
# with its slope at 0 against the constant floor. With the latent smooths
# (V5, #47), each is checked to have mean 0 over the modelled cells, and the
# model with every raw weight 0 against the kdr model with the same floor
# intercepts (by class: V4_class) and its slopes at 0, in greta and in plain
# R: both are then the model with a constant floor per class (or one floor,
# or none). With the shear, the model with its loading b at 0 is checked
# against the same smooths without the shear, in greta and in plain R.
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
#   IR_CUBE_MODEL_OPTIONS='dynamical_model_options(mortality_floor = TRUE,
#     floor_prior = c(1, 4), smooth = smooth_options())' \
#     Rscript R/check_dynamical_model.R
#   IR_CUBE_MODEL_OPTIONS='dynamical_model_options(
#     likelihood = "weighted_binomial")' Rscript R/check_dynamical_model.R
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
# (rho too, when it is a parameter; with the weighted binomial it is fixed)
weighted <- !rho_estimated(model_options)
greta_values <- if (weighted) {
  calculate(p = built$population_mortality_vec, values = values)
} else {
  calculate(p = built$population_mortality_vec, rho = built$terms$rho_types,
            values = values)
}
greta_values <- as.matrix(greta_values)
p_greta <- greta_values[1, grep("^p\\[", colnames(greta_values))]
rho_greta <- if (weighted) built$terms$rho_types else
  greta_values[1, grep("^rho\\[", colnames(greta_values))]

# plain R
logit_r <- c(dynamical_logit(parameters, df, df, x_cell_years,
                             cell_years_index))
p_r <- plogis(logit_r)
rho_r <- c(parameters$rho_types)

if (weighted) {
  # the weighted binomial log likelihood at the replicate rho, its log p and
  # log(1 - p) from the plain-R logit
  stopifnot(identical(unname(rho_r), unname(replicate_rho(types))))
  cat(sprintf("weighted binomial: replicate rho %s\n",
              paste(sprintf("%s %.3f", types, rho_r), collapse = ", ")))
  loglik_r <- weighted_binomial_log_lik(
    df$died, df$mosquito_number,
    log_p = plogis(logit_r, log.p = TRUE),
    log_not_p = plogis(logit_r, lower.tail = FALSE, log.p = TRUE),
    weight = design_effect_weight(df$mosquito_number, rho_r[df$type_id]))
} else {
  # the betabinomial log likelihood, parameterised as in betabinomial_p_rho()
  # and not clamped (dbetabinom() clamps p away from 0 and 1)
  a <- p_r * (1 / rho_r[df$type_id] - 1)
  b <- a * (1 - p_r) / p_r
  loglik_r <- extraDistr::dbbinom(df$died, df$mosquito_number, alpha = a,
                                  beta = b, log = TRUE)
}

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
# species model), and for dynamical_logit(), at each row; the standardised
# kdr at each map cell (NULL without the kdr covariate); and the basis of the
# latent smooths (NULL without them)
map_share <- prediction_share(model_options, map_rows$cell)
map_kdr <- prediction_kdr(model_options, map_rows$cell)
map_basis <- prediction_basis(model_options, map_rows$cell)
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
                               kdr = map_kdr, basis = map_basis)
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
  floor_kind <- kdr_floor(model_options)
  # With its slopes at 0, the model is the model without the kdr covariate,
  # and with the kdr-dependent floor's slope at 0 too, the floor is constant:
  # floor_intercept takes the place of mortality_floor (the same free state,
  # qlogis(floor)). In greta, the log density then differs by the priors of
  # the slopes at 0, and of the floor: N(floor_intercept; logit_beta_moments())
  # in place of Beta(floor) times its Jacobian floor (1 - floor). Not checked
  # in greta for one floor per class, which has no single floor to map to
  options_off <- model_options
  options_off$kdr <- FALSE
  zero_names <- intersect(c(kdr_slope_names, "floor_kdr"),
                          names(built$variables))
  cat("kdr: slopes", toString(zero_names), "\n")
  if (!identical(floor_kind, "class")) {
    built_off <- build(df, options_off)
    free_zero <- free
    for (name in zero_names) {
      free_zero[free_columns(built$model, name)] <- 0
    }
    columns_off <- free_state_columns(built_off$model)
    free_off <- numeric(length(unlist(
      built_off$model$dag$example_parameters(free = TRUE))))
    for (name in names(attr(columns_off, "targets"))) {
      source_name <- if (name == "mortality_floor" && isTRUE(floor_kind)) {
        "floor_intercept"
      } else {
        name
      }
      free_off[columns_off[[attr(columns_off, "targets")[[name]]]]] <-
        free_zero[free_columns(built$model, source_name)]
    }
    expected <- length(zero_names) * dnorm(0, log = TRUE)
    if (isTRUE(floor_kind)) {
      intercept <- free_zero[free_columns(built$model, "floor_intercept")]
      f <- plogis(intercept)
      prior <- logit_beta_moments(model_options$floor_prior)
      expected <- expected +
        dnorm(intercept, prior$mean, prior$sd, log = TRUE) -
        (dbeta(f, model_options$floor_prior[1], model_options$floor_prior[2],
               log = TRUE) + log(f) + log1p(-f))
    }
    ld_difference <- log_density(built$model, free_zero) -
      log_density(built_off$model, free_off) - expected
    cat(sprintf("kdr: slopes 0 vs no kdr covariate%s, greta log density: diff %.3g\n",
                if (isTRUE(floor_kind)) " and a constant floor" else "",
                ld_difference))
    stopifnot(abs(ld_difference) < 1e-6)
    rm(built_off)
  }
  # and in plain R, the predictions at every assay, with the floor's slope at
  # 0 and its intercepts equal, as a constant floor
  zero <- parameters
  zero$kdr_slopes <- lapply(zero$kdr_slopes, function(x) 0 * x)
  off <- zero
  off$kdr_slopes <- list()
  off$options <- options_off
  if (!isFALSE(floor_kind)) {
    zero$floor_kdr[] <- 0
    zero$floor_intercept[] <- zero$floor_intercept[, 1]
    off$floor_kdr <- off$floor_intercept <- NULL
    off$mortality_floor <- plogis(zero$floor_intercept[, 1])
  }
  l_zero <- dynamical_logit(zero, df, df, x_cell_years, cell_years_index)
  l_off <- dynamical_logit(off, df, df, x_cell_years, cell_years_index)
  kdr_difference <- max(abs(clamp(l_zero) - clamp(l_off)))
  cat(sprintf("kdr: slopes 0 vs no kdr covariate%s, plain R, logit: max abs diff %.3g\n",
              if (!isFALSE(floor_kind)) " and a constant floor" else "",
              kdr_difference))
  stopifnot(kdr_difference < 1e-9)
  if (!isFALSE(floor_kind)) {
    floors <- plogis(c(parameters$floor_intercept))
    cat(sprintf("kdr floor (%s): floor at the mean kdr %s, floor_kdr %.3f\n",
                if (isTRUE(floor_kind)) "one" else "by class",
                paste(sprintf("%.3f", floors), collapse = ", "),
                parameters$floor_kdr))
  }
}
if (smooth_on(model_options)) {
  smooth <- model_options$smooth
  kinds <- smooth_kinds(model_options)
  cat(sprintf("smooth: %s; kernel %s, m = (%s), %d basis functions; floor intercepts %s\n",
              paste(sprintf("%s (%s)", kinds,
                            vapply(kinds, function(kind) {
                              if (isTRUE(smooth[[kind]])) "every class" else
                                "pyrethroids and DDT"
                            }, "")), collapse = ", "),
              smooth$kernel, toString(smooth$m), nrow(smooth$indices),
              if (smooth_floor_on(model_options)) smooth$floor_intercepts else
                "none"))
  # each smooth has mean 0 over the modelled cells
  basis_cells <- prediction_basis(model_options, map_rows$cell)
  u_cells <- lapply(parameters$smooth_weights, function(w) {
    c(w %*% t(basis_cells))
  })
  for (kind in kinds) {
    cat(sprintf("smooth %s: sd %.3f, range %.0f km%s; at the cells, mean %.2g, range %.3f to %.3f\n",
                kind, trace[1, paste0("smooth_sd_", kind)],
                smooth_range_km(trace, model_options, kind),
                if (smooth_range_fixed(smooth)) " (fixed)" else "",
                mean(u_cells[[kind]]), min(u_cells[[kind]]),
                max(u_cells[[kind]])))
  }
  stopifnot(all(abs(vapply(u_cells, mean, numeric(1))) < 1e-12))
  shear <- isTRUE(smooth$shear)
  if (shear) {
    cat(sprintf("smooth shear: b %.3f\n", trace[1, "smooth_shear"]))
  }

  # With every raw weight 0 the smooths are 0, and the model is the kdr model
  # with the same floor intercepts (kdr_options(floor = "class") for one per
  # class, as V4_class; TRUE for one; no floor without mortality_floor) and
  # its slopes at 0. In greta, the log density then differs by the priors:
  # those of the smooths (the raw weights N(0, 1) at 0, the exponential sd
  # and inverse range (unless the range is fixed) with the Jacobians of their
  # log free states, and the shear loading) less those of the kdr slopes N(0,
  # 1) at 0
  floor_on <- smooth_floor_on(model_options)
  options_base <- model_options
  options_base$smooth <- FALSE
  options_base$kdr <- kdr_options(floor = if (!floor_on) FALSE else
    if (identical(smooth$floor_intercepts, "class")) "class" else TRUE)
  built_base <- build(df, options_base)
  free_zero <- free
  for (kind in kinds) {
    free_zero[free_columns(built$model, smooth_variable_names(kind)[["raw"]])] <- 0
  }
  columns_base <- free_state_columns(built_base$model)
  base_slopes <- intersect(c(kdr_slope_names, "floor_kdr"),
                           names(built_base$variables))
  free_base <- numeric(length(unlist(
    built_base$model$dag$example_parameters(free = TRUE))))
  for (name in setdiff(names(attr(columns_base, "targets")), base_slopes)) {
    free_base[columns_base[[attr(columns_base, "targets")[[name]]]]] <-
      free_zero[free_columns(built$model, name)]
  }
  trace_zero <- built$model$dag$trace_values(matrix(free_zero, nrow = 1))
  rates <- smooth_prior_rates(smooth)
  expected <- -length(base_slopes) * dnorm(0, log = TRUE)
  for (kind in kinds) {
    sd <- trace_zero[1, paste0("smooth_sd_", kind)]
    expected <- expected + nrow(smooth$indices) * dnorm(0, log = TRUE) +
      dexp(sd, rates$sd, log = TRUE) + log(sd)
    if (!smooth_range_fixed(smooth)) {
      inv_range <- trace_zero[1, paste0("smooth_inv_range_", kind)]
      expected <- expected + dexp(inv_range, rates$range, log = TRUE) +
        log(inv_range)
    }
  }
  if (shear) {
    expected <- expected + dnorm(trace_zero[1, "smooth_shear"],
                                 smooth_shear_prior$mean,
                                 smooth_shear_prior$sd, log = TRUE)
  }
  ld_difference <- log_density(built$model, free_zero) -
    log_density(built_base$model, free_base) - expected
  cat(sprintf("smooth: raw weights 0 vs the kdr model (floor %s) with slopes %s at 0, greta log density: diff %.3g\n",
              deparse(options_base$kdr$floor), toString(base_slopes),
              ld_difference))
  stopifnot(abs(ld_difference) < 1e-6)

  # and in plain R, the predictions at every assay
  zero <- parameters
  zero$smooth_weights <- lapply(zero$smooth_weights, function(x) 0 * x)
  base <- zero
  base$smooth_weights <- list()
  base$options <- built_base$options
  base$kdr_slopes <- lapply(
    setNames(nm = intersect(kdr_slope_names, base_slopes)),
    function(name) rep(0, parameters$n_draws))
  if (floor_on) base$floor_kdr <- rep(0, parameters$n_draws)
  l_zero <- dynamical_logit(zero, df, df, x_cell_years, cell_years_index)
  l_base <- dynamical_logit(base, df, df, x_cell_years, cell_years_index)
  smooth_difference <- max(abs(clamp(l_zero) - clamp(l_base)))
  cat(sprintf("smooth: raw weights 0 vs the kdr model with slopes at 0, plain R, logit: max abs diff %.3g\n",
              smooth_difference))
  stopifnot(smooth_difference < 1e-12)
  rm(built_base)

  # With the shear loading b at 0, the selection smooth is v_s alone, and the
  # model is the one without the shear: in greta, the log density differs by
  # b's prior at 0
  if (shear) {
    options_unsheared <- model_options
    options_unsheared$smooth$shear <- FALSE
    built_unsheared <- build(df, options_unsheared)
    free_b0 <- free
    free_b0[free_columns(built$model, "smooth_shear")] <- 0
    columns_unsheared <- free_state_columns(built_unsheared$model)
    free_unsheared <- numeric(length(unlist(
      built_unsheared$model$dag$example_parameters(free = TRUE))))
    for (name in names(attr(columns_unsheared, "targets"))) {
      free_unsheared[
        columns_unsheared[[attr(columns_unsheared, "targets")[[name]]]]] <-
        free_b0[free_columns(built$model, name)]
    }
    shear_difference <- log_density(built$model, free_b0) -
      log_density(built_unsheared$model, free_unsheared) -
      dnorm(0, smooth_shear_prior$mean, smooth_shear_prior$sd, log = TRUE)
    cat(sprintf("smooth shear: b 0 vs no shear, greta log density: diff %.3g\n",
                shear_difference))
    stopifnot(abs(shear_difference) < 1e-6)
    # and in plain R, the predictions at every assay
    b0 <- parameters
    b0$variables$smooth_shear[] <- 0
    b0$smooth_weights <- smooth_weight_terms(b0$variables, model_options)
    unsheared <- parameters
    unsheared$options <- built_unsheared$options
    unsheared$smooth_weights <- smooth_weight_terms(unsheared$variables,
                                                    unsheared$options)
    l_b0 <- dynamical_logit(b0, df, df, x_cell_years, cell_years_index)
    l_unsheared <- dynamical_logit(unsheared, df, df, x_cell_years,
                                   cell_years_index)
    shear_plain <- max(abs(clamp(l_b0) - clamp(l_unsheared)))
    cat(sprintf("smooth shear: b 0 vs no shear, plain R, logit: max abs diff %.3g\n",
                shear_plain))
    stopifnot(shear_plain < 1e-12)
    rm(built_unsheared)
  }

  # With a fixed range, the model is the one with the range estimated, on
  # the same basis, at an inverse range of 1 / range: in greta, the log
  # density differs by the inverse range's prior and the Jacobian of its log
  # free state; in plain R, the weights of the smooths are the same
  if (smooth_range_fixed(smooth)) {
    options_estimated <- model_options
    options_estimated$smooth[["range"]] <- NULL
    built_estimated <- build(df, options_estimated)
    stopifnot(identical(built_estimated$options$smooth$indices,
                        smooth$indices))
    columns_estimated <- free_state_columns(built_estimated$model)
    targets_estimated <- attr(columns_estimated, "targets")
    free_estimated <- numeric(length(unlist(
      built_estimated$model$dag$example_parameters(free = TRUE))))
    inv_range_names <- vapply(kinds, function(kind) {
      smooth_variable_names(kind)[["inv_range"]]
    }, "")
    for (name in setdiff(names(targets_estimated), inv_range_names)) {
      free_estimated[columns_estimated[[targets_estimated[[name]]]]] <-
        free[free_columns(built$model, name)]
    }
    inv_range <- 1 / smooth[["range"]]
    free_estimated[unlist(columns_estimated[
      unlist(targets_estimated[inv_range_names])])] <- log(inv_range)
    range_prior <- length(kinds) *
      (dexp(inv_range, rates$range, log = TRUE) + log(inv_range))
    range_difference <- log_density(built$model, free) -
      log_density(built_estimated$model, free_estimated) + range_prior
    cat(sprintf("smooth range: fixed at %.0f km vs estimated at it, greta log density: diff %.3g\n",
                1000 * smooth[["range"]], range_difference))
    stopifnot(abs(range_difference) < 1e-6)
    trace_estimated <- built_estimated$model$dag$trace_values(
      matrix(free_estimated, nrow = 1))
    v_estimated <- lapply(
      setNames(nm = unique(sub("\\[.*$", "", colnames(trace_estimated)))),
      extract_parameter, draws_matrix = trace_estimated)
    weights_estimated <- smooth_weight_terms(v_estimated,
                                             built_estimated$options)
    weights_plain <- max(abs(unlist(weights_estimated) -
                               unlist(parameters$smooth_weights)))
    cat(sprintf("smooth range: fixed vs estimated at it, plain R, weights: max abs diff %.3g\n",
                weights_plain))
    stopifnot(weights_plain < 1e-12)
    rm(built_estimated)
  }
  if (floor_on) {
    cat(sprintf("smooth floor (%s): floor where u_f is 0 %s\n",
                smooth$floor_intercepts,
                paste(sprintf("%.3f", plogis(c(parameters$floor_intercept))),
                      collapse = ", ")))
  }
}
# The weighted binomial log density of greta against plain R, at points from
# p near 0 to p near 1 (the logit l from -50 to 50), with and without a floor,
# and as a mixture of two trajectories (the species model). The plain-R
# reference takes log p and log(1 - p) from plogis(), not from
# floored_log_probs(), and the weighted binomial log likelihood from
# dbinom() less the binomial coefficient, at the points where p and 1 - p are
# not too near 0 or 1 for that; at l = 50 with survivors, log(1 - p) from p
# itself would be -Inf. And its simulations (sample(), the beta-binomial at
# the replicate rho) against the beta-binomial's mean n p and variance
# n p (1 - p) (1 + (n - 1) rho), where p is in (0.01, 0.99).
check_weighted_binomial_points <- function() {
  l <- c(-50, -4, -0.5, 0, 1.5, 6, 50)
  n <- c(20, 25, 100, 60, 50, 80, 30)
  died <- c(0, 2, 40, 31, 41, 79, 29)
  rho <- c(0.1, 0.25, 0.15, 0.12, 0.18, 0.09, 0.2)
  weight <- design_effect_weight(n, rho)
  floor <- 0.15
  share <- c(0, 0.3, 0.5, 1, 0.8, 0.2, 0.6)
  l_other <- l - 1.3
  floor_other <- 0.05
  reference <- function(l, floor = NULL) {
    log_q <- plogis(l, log.p = TRUE)
    log_not_q <- plogis(l, lower.tail = FALSE, log.p = TRUE)
    if (is.null(floor)) {
      return(list(log_p = log_q, log_not_p = log_not_q))
    }
    list(log_p = log(floor + (1 - floor) * plogis(l)),
         log_not_p = log1p(-floor) + log_not_q)
  }
  mixed <- function(a, b) {
    list(log_p = log(share * exp(a$log_p) + (1 - share) * exp(b$log_p)),
         log_not_p = log(share * exp(a$log_not_p) +
                           (1 - share) * exp(b$log_not_p)))
  }
  cases <- list(
    `no floor` = list(
      greta = function(lv) floored_log_probs(lv),
      r = reference(l)),
    floor = list(
      greta = function(lv) floored_log_probs(lv, as_data(floor)),
      r = reference(l, floor)),
    mixture = list(
      greta = function(lv) mixture_log_probs(
        share, floored_log_probs(lv, as_data(floor)),
        floored_log_probs(lv - 1.3, as_data(floor_other))),
      r = mixed(reference(l, floor), reference(l_other, floor_other))))
  for (name in names(cases)) {
    case <- cases[[name]]
    lv <- variable(dim = length(l))
    probs <- case$greta(lv)
    y <- as_data(died)
    distribution(y) <- weighted_binomial(n, probs$log_p, probs$log_not_p,
                                         weight)
    m <- model(lv)
    density <- m$dag$generate_log_prob_function(which = "unadjusted")
    ld_greta <- as.numeric(density(tensorflow::tf$constant(
      matrix(l, nrow = 1), dtype = tensorflow::tf$float64)))
    log_lik <- weighted_binomial_log_lik(died, n, case$r$log_p,
                                         case$r$log_not_p, weight)
    stopifnot(all(is.finite(log_lik)))
    # the same points' log p and log(1 - p) in greta, by calculate()
    values <- calculate(log_p = probs$log_p, log_not_p = probs$log_not_p,
                        values = list(lv = l))
    probs_diff <- max(abs(c(values$log_p) - case$r$log_p),
                      abs(c(values$log_not_p) - case$r$log_not_p))
    # and dbinom() less the binomial coefficient, where p is not within 1e-6
    # of 0 or 1
    p <- exp(case$r$log_p)
    usable <- p > 1e-6 & p < 1 - 1e-6
    binomial <- weight * (dbinom(died, n, p, log = TRUE) - lchoose(n, died))
    binomial_diff <- max(abs(binomial[usable] - log_lik[usable]))
    cat(sprintf(paste0("weighted binomial, %s, %d points (l %g to %g): ",
                       "greta %.10g, plain R %.10g, diff %.3g; log p and ",
                       "log(1 - p) max diff %.3g; vs dbinom() at %d: max ",
                       "diff %.3g\n"),
                name, length(l), min(l), max(l), ld_greta, sum(log_lik),
                ld_greta - sum(log_lik), probs_diff, sum(usable),
                binomial_diff))
    stopifnot(abs(ld_greta - sum(log_lik)) < 1e-9, probs_diff < 1e-12,
              binomial_diff < 1e-9)
    nsim <- 20000
    simulated <- matrix(calculate(y, values = list(lv = l), nsim = nsim)$y,
                        nsim)
    middle <- p > 0.01 & p < 0.99
    expected_var <- n * p * (1 - p) * (1 + (n - 1) * rho)
    z_mean <- (colMeans(simulated) - n * p) / sqrt(expected_var / nsim)
    var_ratio <- apply(simulated, 2, var) / expected_var
    cat(sprintf(paste0("weighted binomial, %s: %d simulations at %d points, ",
                       "mean z %s, variance ratio %s\n"),
                name, nsim, sum(middle),
                paste(sprintf("%.2f", z_mean[middle]), collapse = " "),
                paste(sprintf("%.3f", var_ratio[middle]), collapse = " ")))
    stopifnot(all(abs(z_mean[middle]) < 5),
              all(abs(var_ratio[middle] - 1) < 0.08))
  }
}
check_weighted_binomial_points()

cat("all checks passed\n")
