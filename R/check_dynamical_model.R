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
# density) and in plain R (the predictions), and with the kdr-dependent floor,
# with its slope at 0 against the constant floor. With the latent smooths
# (V5, #47), each is checked to have mean 0 over the cells it is centred on
# (every modelled cell, or for the floor smooth of the pyrethroids and DDT,
# the cells with their bioassays) and to equal the smooth from the centred
# basis functions; the model with every raw weight 0 against the model
# without the smooths (smooth_options(selection = FALSE, floor = FALSE)),
# which has the same floor intercepts, in greta and in plain R, and in plain
# R against the closed form with a constant floor per class (or one floor,
# or none) put on by hand; with a fixed range, the model against the one with
# the range estimated at it, in greta and in plain R. The model with the
# other prior of the smooths' sds (PC or half-normal, smooth_options(
# sd_prior = )), and with the other
# prior of the floor where u_f is 0 (the half-normal on the floor itself or
# the logit-normal on its logit, smooth_options(floor_intercept_prior = )),
# is checked against it at the same free state, in greta: the log density
# differs by the two priors (with the Jacobian of the floor's transform),
# computed here.
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
greta_values <- as.matrix(calculate(p = built$population_mortality_vec,
                                    rho = built$terms$rho_types,
                                    values = values))
p_greta <- greta_values[1, grep("^p\\[", colnames(greta_values))]
rho_greta <- greta_values[1, grep("^rho\\[", colnames(greta_values))]

# plain R
logit_r <- c(dynamical_logit(parameters, df, df, x_cell_years,
                             cell_years_index))
p_r <- plogis(logit_r)
rho_r <- c(parameters$rho_types)

# the betabinomial log likelihood, parameterised as in betabinomial_p_rho()
# and not clamped (dbetabinom() clamps p away from 0 and 1)
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
    scale <- sprintf("sd %.3f", trace[1, paste0("smooth_sd_", kind)])
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

  # With every raw weight 0 the smooths are 0, and the model is the one
  # without them, smooth_options(selection = FALSE, floor = FALSE), with the
  # same floor intercepts (the floor per class), or none without
  # mortality_floor. In greta, the log density then differs by the priors of
  # the smooths: the raw weights N(0, 1) at 0, and the sd and inverse range
  # (unless the range is fixed) with the Jacobians of their log free states.
  # And in plain R, the model with the smooths at 0 is the closed form
  # (smooth = FALSE) without a floor, with the floor of each type's class
  # put on by hand: f + (1 - f) q, f = plogis(floor_intercept)
  floor_on <- smooth_floor_on(model_options)
  options_base <- model_options
  options_base$smooth$selection <- FALSE
  options_base$smooth$floor <- FALSE
  built_base <- build(df, options_base)
  free_zero <- free
  for (kind in kinds) {
    free_zero[free_columns(built$model, smooth_variable_names(kind)[["raw"]])] <- 0
  }
  columns_base <- free_state_columns(built_base$model)
  free_base <- numeric(length(unlist(
    built_base$model$dag$example_parameters(free = TRUE))))
  for (name in names(attr(columns_base, "targets"))) {
    free_base[columns_base[[attr(columns_base, "targets")[[name]]]]] <-
      free_zero[free_columns(built$model, name)]
  }
  trace_zero <- built$model$dag$trace_values(matrix(free_zero, nrow = 1))
  rates <- smooth_prior_rates(smooth)
  expected <- 0
  for (kind in kinds) {
    sd <- trace_zero[1, paste0("smooth_sd_", kind)]
    expected <- expected + nrow(smooth$indices) * dnorm(0, log = TRUE) +
      smooth_sd_log_prior(sd, smooth) + log(sd)
    if (!smooth_range_fixed(smooth)) {
      inv_range <- trace_zero[1, paste0("smooth_inv_range_", kind)]
      expected <- expected + dexp(inv_range, rates$range, log = TRUE) +
        log(inv_range)
    }
  }
  ld_difference <- log_density(built$model, free_zero) -
    log_density(built_base$model, free_base) - expected
  cat(sprintf("smooth: raw weights of %s 0 vs the model without them, greta log density: diff %.3g\n",
              toString(kinds), ld_difference))
  stopifnot(abs(ld_difference) < 1e-6)

  # and in plain R, the predictions at every assay: against the model
  # without them, and against the closed form with the floors by hand
  zero <- parameters
  zero$smooth_weights <- lapply(zero$smooth_weights, function(x) 0 * x)
  base <- zero
  base$smooth_weights <- list()
  base$options <- built_base$options
  l_zero <- dynamical_logit(zero, df, df, x_cell_years, cell_years_index)
  l_base <- dynamical_logit(base, df, df, x_cell_years, cell_years_index)
  smooth_difference <- max(abs(clamp(l_zero) - clamp(l_base)))
  closed <- zero
  closed$smooth_weights <- list()
  closed$floor_intercept <- NULL
  closed$mortality_floor <- NULL
  closed$options$smooth <- FALSE
  closed$options$mortality_floor <- FALSE
  p_closed <- plogis(c(dynamical_logit(closed, df, df, x_cell_years,
                                       cell_years_index)))
  if (floor_on) {
    f <- plogis(c(zero$floor_intercept[, smooth_intercept_index(
      model_options, classes_index[df$type_id])]))
    p_closed <- f + (1 - f) * p_closed
  }
  closed_difference <- max(abs(plogis(c(l_zero)) - p_closed))
  cat(sprintf("smooth: raw weights of %s 0 vs the model without them, plain R, logit: max abs diff %.3g; vs the closed form with the floor per class by hand, p: max abs diff %.3g\n",
              toString(kinds), smooth_difference, closed_difference))
  stopifnot(smooth_difference < 1e-12, closed_difference < 1e-12)
  rm(built_base)

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
  # differs by the two priors' log densities at each sd only, the free state
  # being log sd in both
  sd_names <- vapply(kinds, function(kind) {
    smooth_variable_names(kind)[["sd"]]
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
}

cat("all checks passed\n")
