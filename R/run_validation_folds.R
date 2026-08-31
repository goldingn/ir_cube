# Fit the dynamical model to each cross-validation training fold and save
# posterior predictive draws for the held-out data.
#
# This replaces the point-estimate extraction in dynamic_predictive_validation.R
# (#10). The MCMC is unchanged; what is saved is the posterior draws of the
# predicted population fraction and the observation overdispersion, rather than
# their means, so that held-out data can be scored against the full posterior
# predictive distribution. Scoring and plotting are separate, in
# validation_metrics.R and fig_predictive_validation.R, so metrics can be
# revised without refitting.
#
# The null models are cheap and are run here too, so that every candidate is
# stored in the same format and scored by the same code.

source("R/validation_functions.R")
source("R/validation_folds.R")
source("R/null_models.R")
source("R/validation_covariates.R")
source("R/fit_validation_fold.R")

draws_dir <- "outputs/cv_draws"
dir.create(draws_dir, showWarnings = FALSE, recursive = TRUE)

n_draws <- 1000

# store one fold of one model in a consistent format: the draws, and the
# held-out records they correspond to
save_fold <- function(fit, model, experiment, fold) {
  object <- list(
    model = model,
    experiment = experiment,
    fold = fold,
    p_draws = fit$p_draws,
    rho_draws = fit$rho_draws,
    test_df = fit$test_df
  )
  file <- file.path(draws_dir,
                    sprintf("%s__%s__%s.rds", model, experiment, fold))
  saveRDS(object, file)
  cat(sprintf("saved %s (%i held-out assays)\n", file, nrow(fit$test_df)))
  invisible(file)
}

# the numbers of neighbours already chosen by grid search on internal holdouts
# in predictive_validation.R
optimal_nn <- read.csv("outputs/optimal_nn.csv")
neighbours_for <- function(experiment) {
  optimal_nn$n_neighbours[optimal_nn$experiment == experiment]
}

# the three experiments, as a list of training and test folds
folds <- c(
  lapply(
    seq_along(countries_to_validate),
    function(i) list(experiment = "spatial_extrapolation",
                     fold = countries_to_validate[i],
                     training = spatial_extrapolation$training[[i]],
                     test = spatial_extrapolation$test[[i]],
                     n_years_prior = 1)
  ),
  list(
    list(experiment = "spatial_interpolation",
         fold = "all",
         training = spatial_interpolation$training,
         test = spatial_interpolation$test,
         n_years_prior = 1),
    list(experiment = "temporal_forecasting",
         fold = "all",
         training = temporal_forecasting$training,
         test = temporal_forecasting$test,
         n_years_prior = 3)
  )
)


# null models --------------------------------------------------------------

for (fold in folds) {

  save_fold(
    intercept_null_draws(fold$training, fold$test, n_draws = n_draws),
    model = "intercept",
    experiment = fold$experiment,
    fold = fold$fold
  )

  save_fold(
    nn_null_draws(fold$training,
                  fold$test,
                  n_neighbours = neighbours_for(fold$experiment),
                  n_years_prior = fold$n_years_prior,
                  n_draws = n_draws),
    model = "nearest_neighbour",
    experiment = fold$experiment,
    fold = fold$fold
  )

}


# dynamical model ----------------------------------------------------------

# Each fold is a full HMC run, so this is the expensive step: one fit per
# validation country plus one each for interpolation and forecasting. Folds are
# independent, so the loop can be interrupted and resumed; a fold whose draws
# are already on disk is skipped.
for (fold in folds) {

  file <- file.path(draws_dir,
                    sprintf("dynamical__%s__%s.rds", fold$experiment, fold$fold))

  if (file.exists(file)) {
    cat(sprintf("skipping %s (already fitted)\n", file))
    next
  }

  cat(sprintf("\nfitting %s / %s: %i training, %i held-out assays\n",
              fold$experiment, fold$fold,
              nrow(fold$training), nrow(fold$test)))

  fit <- fit_fold(
    train_df = fold$training,
    test_df = fold$test,
    x_cell_years = x_cell_years,
    df = df,
    classes_index = classes_index,
    types = types,
    n_covs = n_covs,
    n_times = n_times,
    n_classes = n_classes,
    n_types = n_types,
    n_regions = n_regions,
    n_countries = n_countries,
    n_sim = n_draws
  )

  # record the worst convergence diagnostic, to be checked before the draws are
  # used for anything
  cat(sprintf("  worst potential scale reduction factor: %.3f\n",
              max(fit$convergence[, 1], na.rm = TRUE)))

  save_fold(fit,
            model = "dynamical",
            experiment = fold$experiment,
            fold = fold$fold)

  gc()

}
