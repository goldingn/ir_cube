# Fit the dynamical model to a single cross-validation fold and save the draws.
#
#   Rscript R/run_one_fold.R <experiment> <fold> [n_chains] [threads]
#
# e.g. Rscript R/run_one_fold.R spatial_extrapolation Kenya 4 4
#
# One fold per process, so that greta's warmup and sampling progress goes to
# this fold's own log rather than being buffered out of sight. Dispatched by
# run_validation_folds.R, but can be run directly to redo a single fold.

arguments <- commandArgs(trailingOnly = TRUE)
experiment_name <- arguments[1]
fold_name <- arguments[2]
n_chains <- if (length(arguments) >= 3) as.integer(arguments[3]) else 4L
threads <- if (length(arguments) >= 4) as.integer(arguments[4]) else 4L
# optional overrides, for smoke-testing the path without a real fit
warmup <- if (length(arguments) >= 5) as.integer(arguments[5]) else 2000L
n_samples <- if (length(arguments) >= 6) as.integer(arguments[6]) else 5000L

# Order matters here, and for two separate reasons. TensorFlow refuses to change
# its thread count once initialised, so that has to be set first. And python has
# to be initialised before terra and sf are attached, because those load the
# system XML libraries, against which the conda environment's pyexpat is then
# resolved and tensorflow_probability fails to import. So: load greta, set
# threads, force python up, and only then source anything else.
suppressMessages(library(greta))

tensorflow_module <- reticulate::import("tensorflow")
tensorflow_module$config$threading$set_intra_op_parallelism_threads(
  as.integer(threads))
tensorflow_module$config$threading$set_inter_op_parallelism_threads(
  as.integer(threads))

invisible(calculate(normal(0, 1), nsim = 1))

source("R/validation_functions.R")
source("R/validation_folds.R")
source("R/validation_covariates.R")
source("R/fit_validation_fold.R")

# find the requested fold
before <- NULL
if (experiment_name == "spatial_extrapolation") {
  index <- match(fold_name, countries_to_validate)
  stopifnot(!is.na(index))
  training <- spatial_extrapolation$training[[index]]
  test <- spatial_extrapolation$test[[index]]
} else if (experiment_name == "spatial_interpolation") {
  training <- spatial_interpolation$training
  test <- spatial_interpolation$test
} else if (experiment_name == "spatial_blocks") {
  source("R/validation_blocks.R")
  index <- as.integer(fold_name)
  stopifnot(!is.na(index), index >= 1, index <= length(spatial_blocks))
  training <- spatial_blocks[[index]]$training
  test <- spatial_blocks[[index]]$test
} else if (experiment_name == "temporal_forecasting") {
  training <- temporal_forecasting$training
  test <- temporal_forecasting$test
  # predictions in the before window are needed for the change-based score, and
  # can only be requested at fitting time
  before <- temporal_forecasting$before
} else {
  stop("unknown experiment: ", experiment_name)
}

# CV_DRAWS_DIR lets a smoke test write somewhere harmless, so a five-sample
# wiring check cannot be mistaken for a fold and skipped later
draws_dir <- Sys.getenv("CV_DRAWS_DIR", "outputs/cv_draws")
destination <- file.path(draws_dir,
                         sprintf("dynamical__%s__%s.rds",
                                 experiment_name, fold_name))

# resume by file: a fold already on disk is not refitted, so the dispatcher can
# be restarted without checking anything itself
if (file.exists(destination)) {
  cat(sprintf("%s | %s / %s already fitted, skipping\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              experiment_name, fold_name))
  cat("FOLD COMPLETE\n")
  quit(save = "no", status = 0)
}

cat(sprintf("%s | %s / %s | %i training, %i held out | %i chains, %i threads\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            experiment_name, fold_name, nrow(training), nrow(test),
            n_chains, threads))
flush(stdout())

elapsed <- system.time(
  fit <- fit_fold(
    train_df = training,
    test_df = test,
    before_df = before,
    x_cell_years = x_cell_years,
    df = df,
    classes_index = classes_index,
    types = types,
    n_covs = n_covs,
    n_times = n_times,
    n_unique_cells = n_unique_cells,
    n_classes = n_classes,
    n_types = n_types,
    n_regions = n_regions,
    n_countries = n_countries,
    n_chains = n_chains,
    warmup = warmup,
    n_samples = n_samples
  )
)

cat(sprintf("%s | fit complete in %.1f hours\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            elapsed[["elapsed"]] / 3600))

dir.create(draws_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(
  list(model = "dynamical",
       experiment = experiment_name,
       fold = fold_name,
       # The draws object, and the greta arrays the predictions came from, so
       # that calculate() can be used on a reloaded fold to predict a quantity
       # that was not asked for at fitting time. This is the only supported way
       # to predict from a fitted greta model, so it is kept deliberately even
       # though the scoring path reads only the matrices below; it is also the
       # bulk of each file's size.
       draws = fit$draws,
       prediction_arrays = fit$prediction_arrays,
       p_draws = fit$p_draws,
       rho_class_draws = fit$rho_class_draws,
       class_id = fit$class_id,
       test_df = fit$test_df,
       p_draws_before = fit$p_draws_before,
       before_df = fit$before_df,
       convergence = fit$convergence,
       ess = fit$ess,
       ess_p = fit$ess_p,
       ess_rho = fit$ess_rho,
       n_sampled = fit$n_sampled,
       n_chains = fit$n_chains),
  destination
)

cat(sprintf("%s | saved. p ESS median %.0f, rho ESS median %.0f, worst Rhat %.3f\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            median(fit$ess_p, na.rm = TRUE),
            median(fit$ess_rho, na.rm = TRUE),
            max(fit$convergence[, 1], na.rm = TRUE)))
cat("FOLD COMPLETE\n")
