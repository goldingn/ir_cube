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

# Initialise greta's python session before terra and sf are attached. Those
# load the system XML libraries, which the conda environment's pyexpat is then
# linked against, and tensorflow_probability fails to import as a result.
suppressMessages(library(greta))
invisible(calculate(normal(0, 1), nsim = 1))

library(future)
library(future.apply)

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

  null_file <- function(model) {
    file.path(draws_dir,
              sprintf("%s__%s__%s.rds", model, fold$experiment, fold$fold))
  }

  if (!file.exists(null_file("intercept"))) {
    save_fold(
      intercept_null_draws(fold$training, fold$test, n_draws = n_draws),
      model = "intercept",
      experiment = fold$experiment,
      fold = fold$fold
    )
  }

  if (file.exists(null_file("nearest_neighbour"))) next

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

# Each fold is a full HMC run, and this is the expensive step: one fit per
# validation country plus one each for interpolation and forecasting.
#
# How the work is divided follows from how greta uses the machine. Chains are
# vectorised into a single TensorFlow op rather than run one per core, and that
# op scales poorly beyond about four threads, so the efficient arrangement is
# few chains and several folds at once rather than many chains on one fold.
# Measured on one fold of this model: 8 chains cost 70.6 s per iteration
# against 6.96 s for 2 chains, while effective samples per draw are unchanged,
# and confining a fold to 2 threads costs it only about a fifth of its speed.
#
# Folds are independent and each is skipped if its draws are already on disk,
# so the run can be interrupted and resumed.
n_workers <- 4
threads_per_worker <- 2

pending <- Filter(
  function(fold) {
    !file.exists(file.path(
      draws_dir,
      sprintf("dynamical__%s__%s.rds", fold$experiment, fold$fold)))
  },
  folds
)

cat(sprintf("\n%i folds to fit, %i at a time, %i chains and %i threads each\n",
            length(pending), n_workers, 2, threads_per_worker))

if (length(pending) > 0) {

  plan(future.callr::callr, workers = n_workers)

  fold_diagnostics <- future_lapply(
    pending,
    function(fold) {

      # each worker is a fresh R process, so everything the model definition
      # needs is loaded here rather than inherited: greta.dynamics for the
      # iteration, dplyr for the index lookups, and functions.R for
      # betabinomial_p_rho()
      suppressMessages({
        library(greta)
        library(greta.dynamics)
        library(dplyr)
      })
      source("R/functions.R")
      source("R/validation_functions.R")
      source("R/fit_validation_fold.R")

      fit <- fit_fold(
        train_df = fold$training,
        test_df = fold$test,
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
        n_sim = n_draws,
        threads = threads_per_worker
      )

      object <- list(
        model = "dynamical",
        experiment = fold$experiment,
        fold = fold$fold,
        p_draws = fit$p_draws,
        rho_draws = fit$rho_draws,
        test_df = fit$test_df,
        convergence = fit$convergence,
        ess = fit$ess,
        n_sampled = fit$n_sampled
      )
      saveRDS(object,
              file.path(draws_dir,
                        sprintf("dynamical__%s__%s.rds",
                                fold$experiment, fold$fold)))

      # diagnostics worth seeing before the draws are used for anything
      data.frame(experiment = fold$experiment,
                 fold = fold$fold,
                 n_sampled = fit$n_sampled,
                 min_ess = min(fit$ess, na.rm = TRUE),
                 median_ess = median(fit$ess, na.rm = TRUE),
                 draws_per_ess = fit$draws_per_ess,
                 worst_rhat = max(fit$convergence[, 1], na.rm = TRUE))

    },
    future.seed = TRUE,
    future.globals = c("x_cell_years", "df", "classes_index", "types",
                       "n_covs", "n_times", "n_classes", "n_types",
                       "n_regions", "n_countries", "n_draws", "draws_dir",
                       "threads_per_worker")
  )

  cat("\nper-fold diagnostics:\n")
  print(do.call(rbind, fold_diagnostics))

}
