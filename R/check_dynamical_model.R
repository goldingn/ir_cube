# Check that the greta model (build_dynamical_model(), with the closed-form op)
# and the plain-R predictions (R/dynamical_predictions.R) agree, at a random
# free state, for a given set of model options:
#   - predicted mortality at every training assay, and rho per type
#   - the log likelihood: the model's log density over all the data, less that
#     of the same model with one assay in the likelihood (which shares every
#     prior and Jacobian term), against the plain-R betabinomial log likelihood
#     of the other assays
# and print the log density itself, for regression checks between versions.
#
#   Rscript R/check_dynamical_model.R '<options>' [seed] [sd]
# (the free state is N(0, sd^2), sd 0.5 by default; a smaller sd avoids states
# where p rounds to 1 at assays with survivors, and the log density is NaN)
# e.g.
#   Rscript R/check_dynamical_model.R 'dynamical_model_options(rho = "type")'
#
# Run with the greta 0.6 environment (doc/cv_run_plan.md, section 1).

arguments <- commandArgs(trailingOnly = TRUE)
options_text <- if (length(arguments) >= 1) arguments[1] else
  "dynamical_model_options()"
seed <- if (length(arguments) >= 2) as.integer(arguments[2]) else 1L
free_sd <- if (length(arguments) >= 3) as.numeric(arguments[3]) else 0.5

source("R/greta_setup.R")
start_greta(threads = 4)
suppressMessages({
  sink("/dev/null")
  source("R/validation_folds.R")
  source("R/validation_covariates.R")
  sink()
})
source("R/validation_functions.R")
source("R/dynamical_predictions.R")
source("R/two_stage_map_functions.R")

model_options <- eval(parse(text = options_text))
cat("options:", options_text, "\n")

build <- function(train_df) {
  build_dynamical_model(train_df = train_df,
                        df = df,
                        x_cell_years = x_cell_years,
                        cell_years_index = cell_years_index,
                        classes_index = classes_index,
                        types = types,
                        options = model_options,
                        x_cells_init = x_cells_init)
}

built <- build(df)
built_one <- build(df[1, ])

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

# the free-state elements of a variable, located as in logit_init_mean_draws()
free_columns <- function(model, name) {
  dag <- model$dag
  free_list <- dag$example_parameters(free = TRUE)
  sizes <- vapply(free_list, length, integer(1))
  ends <- cumsum(sizes)
  node_names <- vapply(dag$node_list, function(node) node$unique_name,
                       character(1))
  node <- greta:::get_node(model$target_greta_arrays[[name]])$unique_name
  i <- match(dag$get_tf_names()[match(node, node_names)], names(free_list))
  (ends[i] - sizes[i] + 1):ends[i]
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
             options = model_options)
# the plain-R path takes the covariates as an argument here, not from the fold

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
p_r <- c(dynamical_predictions(fold, select(df, -country_id), df,
                               x_cell_years, cell_years_index, classes_index,
                               types, draw_index = 1,
                               x_cells_init = x_cells_init))
rho_r <- c(dynamical_parameter_draws(fold, df, classes_index, types,
                                     draw_index = 1)$rho_types)

# the betabinomial log likelihood, parameterised as in betabinomial_p_rho() and
# not clamped (dbetabinom() clamps p away from 0 and 1)
a <- p_r * (1 / rho_r[df$type_id] - 1)
b <- a * (1 - p_r) / p_r
loglik_r <- extraDistr::dbbinom(df$died, df$mosquito_number, alpha = a,
                                beta = b, log = TRUE)

cat(sprintf("free parameters %d, log density %.10g\n", n_free, ld_all))
cat(sprintf("p: max abs diff %.3g, max logit diff %.3g (range %.3g-%.3g)\n",
            max(abs(p_greta - p_r)),
            max(abs(qlogis(pmin(pmax(p_greta, 1e-12), 1 - 1e-12)) -
                      qlogis(pmin(pmax(p_r, 1e-12), 1 - 1e-12)))),
            min(p_r), max(p_r)))
cat(sprintf("rho: max abs diff %.3g\n", max(abs(rho_greta - rho_r))))
cat(sprintf("log likelihood of assays 2..n: greta %.10g, plain R %.10g, diff %.3g\n",
            ld_all - ld_one, sum(loglik_r[-1]),
            ld_all - ld_one - sum(loglik_r[-1])))

# the map path (R/two_stage_map_functions.R), at the data cells in 2000, 2012
# and 2024, against dynamical_predictions()
map_years <- c(2000, 2012, 2024)
map_rows <- df %>%
  distinct(cell, cell_id) %>%
  mutate(country_name = countries[
    built$lookups$cell_country_lookup[cell_id]])
covariates <- map_covariates(map_rows$cell, baseline_year, max(map_years))
parameters <- dynamical_parameter_draws(fold, df, classes_index, types,
                                        draw_index = 1)
logit_init_all <- map_logit_init(trace, NULL, types, classes_index, countries,
                                 regions, country_region_lookup(),
                                 options = model_options)
cell_country_index <- match(map_rows$country_name,
                            dimnames(logit_init_all)[[2]])
map_difference <- 0
for (k in seq_along(types)) {
  dyn <- dynamical_logit_chunk(
    effect = matrix(parameters$effect_type[, , k], nrow = 1),
    logit_init = map_cell_logit_init(logit_init_all, cell_country_index, k,
                                     covariates$init),
    time_varying = covariates$time_varying,
    flat = covariates$flat,
    years = baseline_year:max(map_years), years_keep = map_years,
    floor = parameters$mortality_floor,
    kappa = if (!is.null(parameters$kappa_type)) parameters$kappa_type[, k])
  for (y in map_years) {
    rows <- tibble(cell_id = map_rows$cell_id, type_id = k,
                   year_id = y - baseline_year + 1)
    p_rows <- c(dynamical_predictions(fold, rows, df, x_cell_years,
                                      cell_years_index, classes_index, types,
                                      draw_index = 1,
                                      x_cells_init = x_cells_init))
    l_rows <- qlogis(pmin(pmax(p_rows, 1e-12), 1 - 1e-12))
    l_map <- pmin(pmax(c(dyn[[as.character(y)]]), qlogis(1e-12)),
                  qlogis(1 - 1e-12))
    map_difference <- max(map_difference, abs(l_map - l_rows))
  }
}
cat(sprintf("map path vs plain R, logit, %d cells x %d types x %d years: max abs diff %.3g\n",
            nrow(map_rows), length(types), length(map_years), map_difference))
