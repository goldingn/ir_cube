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
# (V5, #47), each is checked to have mean 0 over the cells it is centred on
# (every modelled cell, or for the floor smooth of the pyrethroids and DDT,
# the cells with their bioassays) and to equal the smooth from the centred
# basis functions, and the model with every raw weight 0 against the kdr
# model with the same floor intercepts (by class: V4_class) and its slopes
# at 0, in greta and in plain R: both are then the model with a constant
# floor per class (or one floor, or none). With the shear, the model with
# its loading b at 0 is checked against the same smooths without the shear,
# in greta and in plain R. The model with the other prior of the smooths'
# sds (PC or half-normal, smooth_options(sd_prior = )), and with the other
# prior of the floor where u_f is 0 (the half-normal on the floor itself or
# the logit-normal on its logit, smooth_options(floor_intercept_prior = )),
# is checked against it at the same free state, in greta: the log density
# differs by the two priors (with the Jacobian of the floor's transform),
# computed here. With
# the smooth of the initial state (smooth_options(init = TRUE)), the
# selection and floor smooths' raw weights at 0 are checked against the
# model without them (the same smooth of the initial state, a constant floor
# per class), in greta and in plain R; and the plain-R initial state at every
# modelled cell and type, with selection, reversion and the floor switched
# off, against logit_init_mean + lambda u_init(x) + the centred covariates'
# effects, computed here from the variables, the basis functions and the
# spectral density, with its mean over the cells logit_init_mean.
#
#   IR_CUBE_MODEL_OPTIONS='<options>' Rscript R/check_dynamical_model.R [seed] [sd]
# (the free state is N(0, sd^2), sd 0.5 by default; a smaller sd avoids states
# where p rounds to 1 at assays with survivors, and the log density is NaN)
# Without IR_CUBE_MODEL_OPTIONS, the defaults: V5h (#47), with the latent
# smooths and the floor per class. Other models, e.g.
#   IR_CUBE_MODEL_OPTIONS='dynamical_model_options(reversion = FALSE)' \
#     Rscript R/check_dynamical_model.R
#   IR_CUBE_MODEL_OPTIONS='dynamical_model_options(mortality_floor = FALSE,
#     smooth = FALSE)' Rscript R/check_dynamical_model.R  # the default before V5h
#   IR_CUBE_MODEL_OPTIONS='dynamical_model_options(mortality_floor = FALSE,
#     smooth = FALSE, species = species_options())' \
#     Rscript R/check_dynamical_model.R
#   IR_CUBE_MODEL_OPTIONS='dynamical_model_options(mortality_floor = FALSE,
#     smooth = FALSE, kdr = kdr_options())' Rscript R/check_dynamical_model.R
#   IR_CUBE_MODEL_OPTIONS='dynamical_model_options(smooth = smooth_options(
#     range = NULL, sd_prior = c(1, 0.05)))' \
#     Rscript R/check_dynamical_model.R  # V5
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

# The log prior density of the floor where u_f is 0, with greta's Jacobian,
# at its free state, the logit floor l, for the smooth options `smooth` (or
# "beta_moments"): with the half-normal prior on f0 = plogis(l)
# (smooth_floor_log_prior()), plus log f0 + log(1 - f0); with the other,
# N(l; the logit moments of Beta(floor_prior)), the free state being l itself
floor_free_log_prior <- function(l, smooth) {
  if (is.list(smooth) && smooth_floor_prior_family(smooth) == "half_normal") {
    f0 <- plogis(l)
    return(smooth_floor_log_prior(f0, smooth) + log(f0) + log1p(-f0))
  }
  prior <- logit_beta_moments(model_options$floor_prior)
  dnorm(l, prior$mean, prior$sd, log = TRUE)
}

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
  # each smooth has mean 0 over the cells of the bioassays it applies to:
  # every modelled cell, or for a smooth of the pyrethroids and DDT, the
  # cells with their bioassays (smooth_centre()); and the smooth from the
  # basis with its column of ones and the weights with their centring term
  # is the smooth from the centred basis functions
  basis_cells <- prediction_basis(model_options, map_rows$cell)
  u_cells <- lapply(parameters$smooth_weights, function(w) {
    c(w %*% t(basis_cells))
  })
  class_cells <- unique(df$cell[df$insecticide_class %in% smooth_classes])
  centre_rows <- function(kind) {
    if (identical(smooth[[kind]], "class") &&
        !is.null(smooth$class_column_means)) {
      which(map_rows$cell %in% class_cells)
    } else {
      seq_len(nrow(map_rows))
    }
  }
  coords_cells <- smooth_cell_coords(map_rows$cell, smooth$crs)
  centring_difference <- max(vapply(kinds, function(kind) {
    w <- parameters$smooth_weights[[kind]]
    centred <- smooth_centred_basis(smooth, coords_cells, kind)
    max(abs(u_cells[[kind]] - c(w[, seq_len(ncol(centred)), drop = FALSE] %*%
                                  t(centred))))
  }, numeric(1)))
  cat(sprintf("smooths: with the centring term vs from the centred basis: max abs diff %.3g; centred over %s\n",
              centring_difference,
              paste(sprintf("%s %d cells", kinds, vapply(kinds, function(kind) {
                length(centre_rows(kind))
              }, integer(1))), collapse = ", ")))
  stopifnot(centring_difference < 1e-12)
  for (kind in kinds) {
    scale <- if (kind == "init") {
      sprintf("sd 1 (fixed), loadings %s", paste(sprintf(
        "%.3f", trace[1, grep("^smooth_loading_init\\[", colnames(trace))]),
        collapse = ", "))
    } else {
      sprintf("sd %.3f", trace[1, paste0("smooth_sd_", kind)])
    }
    cat(sprintf("smooth %s: %s, range %.0f km%s; at the cells, mean %.2g (over the cells it is centred on %.2g), range %.3f to %.3f\n",
                kind, scale,
                smooth_range_km(trace, model_options, kind),
                if (smooth_range_fixed(smooth)) " (fixed)" else "",
                mean(u_cells[[kind]]), mean(u_cells[[kind]][centre_rows(kind)]),
                min(u_cells[[kind]]), max(u_cells[[kind]])))
  }
  stopifnot(all(vapply(kinds, function(kind) {
    abs(mean(u_cells[[kind]][centre_rows(kind)])) < 1e-12
  }, logical(1))))
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
  # 1) at 0. With the smooth of the initial state, which has no counterpart
  # in the kdr model, the raw weights of the other smooths only are set to 0,
  # and the model is the one without them: smooth_options(selection = FALSE,
  # floor = FALSE), with the same smooth of the initial state and floor
  # intercepts
  floor_on <- smooth_floor_on(model_options)
  init_on <- smooth_init_on(model_options)
  zeroed <- setdiff(kinds, "init")
  options_base <- model_options
  if (init_on) {
    options_base$smooth$selection <- FALSE
    options_base$smooth$floor <- FALSE
    options_base$smooth$shear <- FALSE
  } else {
    options_base$smooth <- FALSE
    options_base$kdr <- kdr_options(floor = if (!floor_on) FALSE else
      if (identical(smooth$floor_intercepts, "class")) "class" else TRUE)
  }
  built_base <- build(df, options_base)
  free_zero <- free
  for (kind in zeroed) {
    free_zero[free_columns(built$model, smooth_variable_names(kind)[["raw"]])] <- 0
  }
  columns_base <- free_state_columns(built_base$model)
  base_slopes <- intersect(c(kdr_slope_names, "floor_kdr"),
                           names(built_base$variables))
  free_base <- numeric(length(unlist(
    built_base$model$dag$example_parameters(free = TRUE))))
  # the kdr model's floor intercepts take the free state of the floor at a
  # flat smooth (floor_flat, V5f), which is its logit, the same value
  flat_to_base <- !init_on && "floor_flat" %in% names(built$variables)
  for (name in setdiff(names(attr(columns_base, "targets")), base_slopes)) {
    own <- if (flat_to_base && name == "floor_intercept") "floor_flat" else
      name
    free_base[columns_base[[attr(columns_base, "targets")[[name]]]]] <-
      free_zero[free_columns(built$model, own)]
  }
  trace_zero <- built$model$dag$trace_values(matrix(free_zero, nrow = 1))
  rates <- smooth_prior_rates(smooth)
  expected <- -length(base_slopes) * dnorm(0, log = TRUE)
  if (flat_to_base) {
    logit_f0 <- free_zero[free_columns(built$model, "floor_flat")]
    expected <- expected + sum(floor_free_log_prior(logit_f0, smooth) -
                                 floor_free_log_prior(logit_f0, "beta_moments"))
  }
  for (kind in zeroed) {
    sd <- trace_zero[1, paste0("smooth_sd_", kind)]
    expected <- expected + nrow(smooth$indices) * dnorm(0, log = TRUE) +
      smooth_sd_log_prior(sd, smooth) + log(sd)
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
  base_label <- if (init_on) {
    "the model without them (smooth of the initial state only)"
  } else {
    sprintf("the kdr model (floor %s) with slopes %s at 0",
            deparse(options_base$kdr$floor), toString(base_slopes))
  }
  cat(sprintf("smooth: raw weights of %s 0 vs %s, greta log density: diff %.3g\n",
              toString(zeroed), base_label, ld_difference))
  stopifnot(abs(ld_difference) < 1e-6)

  # and in plain R, the predictions at every assay
  zero <- parameters
  zero$smooth_weights[zeroed] <- lapply(zero$smooth_weights[zeroed],
                                        function(x) 0 * x)
  base <- zero
  base$smooth_weights <- zero$smooth_weights[setdiff(names(
    zero$smooth_weights), zeroed)]
  base$options <- built_base$options
  if (!init_on) {
    base$kdr_slopes <- lapply(
      setNames(nm = intersect(kdr_slope_names, base_slopes)),
      function(name) rep(0, parameters$n_draws))
    if (floor_on) base$floor_kdr <- rep(0, parameters$n_draws)
  }
  l_zero <- dynamical_logit(zero, df, df, x_cell_years, cell_years_index)
  l_base <- dynamical_logit(base, df, df, x_cell_years, cell_years_index)
  smooth_difference <- max(abs(clamp(l_zero) - clamp(l_base)))
  cat(sprintf("smooth: raw weights of %s 0 vs %s, plain R, logit: max abs diff %.3g\n",
              toString(zeroed), base_label, smooth_difference))
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

  # The model with the other prior of the smooths' sds (the PC prior's
  # default, c(1, 0.05), for the half-normal, or the half-normal with scale
  # 0.5 for the PC prior), at the same free state: in greta, the log density
  # differs by the two priors' log densities at each sd (and loading of the
  # smooth of the initial state) only, the free state being log sd in both
  sd_names <- vapply(kinds, function(kind) {
    smooth_variable_names(kind)[[if (kind == "init") "loading" else "sd"]]
  }, "")
  options_other <- model_options
  options_other$smooth$sd_prior <-
    if (smooth_sd_prior_family(smooth) == "half_normal") c(1, 0.05) else
      list(family = "half_normal", scale = 0.5)
  built_other <- build(df, options_other)
  columns_other <- free_state_columns(built_other$model)
  targets_other <- attr(columns_other, "targets")
  stopifnot(setequal(names(targets_other), names(attr(free_state_columns(
    built$model), "targets"))))
  free_other <- numeric(length(free))
  for (name in names(targets_other)) {
    free_other[columns_other[[targets_other[[name]]]]] <-
      free[free_columns(built$model, name)]
  }
  sds <- unlist(lapply(sd_names, function(name) {
    c(extract_parameter(trace, name))
  }))
  prior_difference <- sum(smooth_sd_log_prior(sds, smooth) -
                            smooth_sd_log_prior(sds, options_other$smooth))
  sd_prior_difference <- log_density(built$model, free) -
    log_density(built_other$model, free_other) - prior_difference
  cat(sprintf("smooth sd prior: %s vs %s at sds %s, greta log density: diff %.4f, analytic %.4f, residual %.3g\n",
              smooth_sd_prior_label(smooth),
              smooth_sd_prior_label(options_other$smooth),
              paste(sprintf("%.3f", sds), collapse = ", "),
              prior_difference + sd_prior_difference, prior_difference,
              sd_prior_difference))
  stopifnot(abs(sd_prior_difference) < 1e-6)
  rm(built_other)

  # The model with the other prior of the floor where u_f is 0
  # (smooth_options(floor_intercept_prior = )), at the same free state: the
  # half-normal prior's floor_flat and the other's floor_intercept have the
  # same free state, the logit floor l, so the log density differs by the
  # priors alone: log N(f0; 0, s^2) - log P(0 < f0 < 1) + log f0 + log(1 -
  # f0) (greta's Jacobian of f0 = plogis(l)) against log N(l; the logit
  # moments of Beta(floor_prior)), computed here (floor_free_log_prior())
  if (floor_on) {
    options_floor <- model_options
    half_normal <- smooth_floor_prior_family(smooth) == "half_normal"
    options_floor$smooth$floor_intercept_prior <- if (half_normal) {
      "beta_moments"
    } else {
      list(family = "half_normal", scale = 0.05)
    }
    built_floor <- build(df, options_floor)
    columns_floor <- free_state_columns(built_floor$model)
    targets_floor <- attr(columns_floor, "targets")
    floor_name <- function(half) if (half) "floor_flat" else "floor_intercept"
    free_floor <- numeric(length(free))
    for (name in names(targets_floor)) {
      own <- if (name == floor_name(!half_normal)) floor_name(half_normal) else
        name
      free_floor[columns_floor[[targets_floor[[name]]]]] <-
        free[free_columns(built$model, own)]
    }
    logit_f0 <- free[free_columns(built$model, floor_name(half_normal))]
    prior_difference <- sum(
      floor_free_log_prior(logit_f0, smooth) -
        floor_free_log_prior(logit_f0, options_floor$smooth))
    floor_prior_difference <- log_density(built$model, free) -
      log_density(built_floor$model, free_floor) - prior_difference
    cat(sprintf("smooth floor prior: %s vs %s at f0 %s, greta log density: diff %.4f, analytic %.4f, residual %.3g\n",
                deparse(smooth$floor_intercept_prior),
                deparse(options_floor$smooth$floor_intercept_prior),
                paste(sprintf("%.3f", plogis(logit_f0)), collapse = ", "),
                prior_difference + floor_prior_difference, prior_difference,
                floor_prior_difference))
    stopifnot(abs(floor_prior_difference) < 1e-6)
    rm(built_floor)
  }

  if (floor_on) {
    cat(sprintf("smooth floor (%s): floor where u_f is 0 %s\n",
                smooth$floor_intercepts,
                paste(sprintf("%.3f", plogis(c(parameters$floor_intercept))),
                      collapse = ", ")))
  }

  # The initial state with the smooth of the initial state, at every modelled
  # cell and type, computed here from the variables: the basis functions
  # (hsgp_basis()) centred over the modelled cells, the weights from the
  # squared exponential's spectral density at sd 1, and the covariates
  # centred at their mean over the modelled cells,
  #   logit_init_relative(x, k) = logit_init_mean[k] + lambda[k] u_init(x) +
  #                               (x_init(x) - mean) init_coef[, k]
  #   q_0 = init_frac_min + (1 - init_frac_min) plogis(logit_init_relative)
  # against the plain-R prediction for year index 1 with no selection
  # (exp(beta) 0), no reversion and no floor, which is then logit q_0
  if (init_on) {
    stopifnot(nrow(map_rows) == max(df$cell_id))
    value_of <- function(name) c(extract_parameter(trace, name))
    phi <- hsgp_basis(smooth_cell_coords(map_rows$cell, smooth$crs), smooth)
    phi <- sweep(phi, 2, colMeans(phi))
    omega <- sqrt((pi * smooth$indices[, 1] / (2 * smooth$half_width[1])) ^ 2 +
                    (pi * smooth$indices[, 2] / (2 * smooth$half_width[2])) ^ 2)
    init_range <- if (smooth_range_fixed(smooth)) smooth[["range"]] else
      1 / trace[1, "smooth_inv_range_init"]
    ell <- init_range / 2
    spectral <- if (smooth$kernel == "se") {
      sqrt(2 * pi) * ell * exp(-ell ^ 2 * omega ^ 2 / 4)
    } else {
      c(smooth_sqrt_spectral(omega, 1, 1 / init_range, smooth$kernel))
    }
    u_init <- c(phi %*% (spectral * value_of("smooth_raw_init")))
    lambda <- value_of("smooth_loading_init")
    logit_init_mean <- value_of("logit_init_mean")
    init_min <- init_frac_constants(types)$min
    covariate_effect <- matrix(0, nrow(map_rows), length(types))
    if (!is.null(model_options$init_covariates)) {
      x_init_cells <- x_cells_init[map_rows$cell_id,
                                   model_options$init_covariates,
                                   drop = FALSE]
      x_mean <- colMeans(x_cells_init[seq_len(max(df$cell_id)),
                                      model_options$init_covariates,
                                      drop = FALSE])
      init_coef <- matrix(value_of("init_coef"),
                          length(model_options$init_covariates))
      covariate_effect <- sweep(x_init_cells, 2, x_mean) %*% init_coef
    }
    l_relative <- sweep(outer(u_init, lambda) + covariate_effect, 2,
                        logit_init_mean, FUN = "+")
    q_0 <- sweep(sweep(plogis(l_relative), 2, 1 - init_min, FUN = "*"), 2,
                 init_min, FUN = "+")
    # logit q_0, with 1 - q_0 = (1 - init_frac_min) plogis(-l) computed
    # directly rather than from q_0 near 1
    log_not_q_0 <- sweep(plogis(-l_relative, log.p = TRUE), 2,
                         log1p(-init_min), FUN = "+")
    l_expected <- c(log(q_0) - log_not_q_0)
    off <- parameters
    off$effect_type[] <- 0
    if (!is.null(off$kappa_type)) off$kappa_type[] <- 0
    if (!is.null(off$mortality_floor)) off$mortality_floor[] <- 0
    if (!is.null(off$floor_intercept)) off$floor_intercept[] <- -Inf
    rows <- tibble(cell_id = rep(map_rows$cell_id, length(types)),
                   type_id = rep(seq_along(types), each = nrow(map_rows)),
                   year_id = 1L)
    l_plain <- c(dynamical_logit(off, rows, df, x_cell_years,
                                 cell_years_index))
    init_difference <- max(abs(l_plain - l_expected))
    # the plain-R logit relative initial state, averaged over the cells
    q_plain <- matrix(plogis(l_plain), nrow(map_rows))
    l_relative_plain <- qlogis(sweep(sweep(q_plain, 2, init_min), 2,
                                     1 - init_min, FUN = "/"))
    mean_difference <- max(abs(colMeans(l_relative_plain) - logit_init_mean))
    cat(sprintf(paste0("smooth init: initial state at %d cells x %d types ",
                       "vs computed here, logit: max abs diff %.3g; mean ",
                       "logit relative initial state over the cells vs ",
                       "logit_init_mean: max abs diff %.3g; lambda u_init ",
                       "%.2f to %.2f\n"),
                nrow(map_rows), length(types), init_difference,
                mean_difference, min(outer(u_init, lambda)),
                max(outer(u_init, lambda))))
    stopifnot(init_difference < 1e-9, mean_difference < 1e-8)
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
