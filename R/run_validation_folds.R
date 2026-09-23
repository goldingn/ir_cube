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

# This process fits the null models and dispatches the dynamical folds as
# separate processes; it never calls fit_fold() itself. So it needs neither
# greta's python session nor the covariate extraction, both of which used to be
# paid for here for nothing. run_one_fold.R does that setup, in the order it
# has to be done in (#12 review).
source("R/validation_functions.R")
source("R/validation_folds.R")
source("R/null_models.R")

draws_dir <- "outputs/cv_draws"
dir.create(draws_dir, showWarnings = FALSE, recursive = TRUE)

n_draws <- 1000

# store one fold of one model in a consistent format: the draws, and the
# held-out records they correspond to
save_fold <- function(fit, model, experiment, fold, label = experiment) {
  object <- list(
    model = model,
    # the experiment this fold is pooled under when scored, which for the
    # forecasting folds is the origin rather than the bare experiment name; the
    # file keeps the plain name so all the forecasting folds sit together
    experiment = label,
    fold = fold,
    p_draws = fit$p_draws,
    # the overdispersion each null's own residuals imply, a diagnostic only:
    # every model is scored at the external replicate-based estimate
    rho_implied = fit$rho_implied,
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
# Keyed on the experiment name, which is a free-text match: an experiment with
# no row returns numeric(0), and the nearest neighbour null then silently
# reduces to its prior (see predict_null_fixed_nn_counts). Fail here instead.
neighbours_for <- function(experiment) {
  value <- optimal_nn$n_neighbours[optimal_nn$experiment == experiment]
  if (length(value) != 1) {
    stop("outputs/optimal_nn.csv has ", length(value), " rows for experiment '",
         experiment, "'; exactly one is needed. Add it, or point this ",
         "experiment at the tuned value for the experiment it most resembles.")
  }
  value
}

# The folds still to fit. The six leave-one-country-out folds and the
# interpolation fold from the September run are kept as they are: the review of
# #12 found them informative, and refitting them would change nothing. What is
# new is the sub-national block folds, which replace leave-one-country-out as
# the test of spatial skill, and the forecasting fold, whose split was leaking
# post-horizon data. Nulls are rebuilt for every fold either way, since they
# cost minutes.
source("R/validation_blocks.R")

folds <- c(
  lapply(
    seq_along(countries_to_validate),
    function(i) list(experiment = "spatial_extrapolation",
                     fold = countries_to_validate[i],
                     training = spatial_extrapolation$training[[i]],
                     test = spatial_extrapolation$test[[i]],
                     n_years_prior = 1)
  ),
  lapply(
    seq_along(spatial_blocks),
    function(i) list(experiment = "spatial_blocks",
                     fold = as.character(i),
                     training = spatial_blocks[[i]]$training,
                     test = spatial_blocks[[i]]$test,
                     n_years_prior = 1)
  ),
  list(
    list(experiment = "spatial_interpolation",
         fold = "all",
         training = spatial_interpolation$training,
         test = spatial_interpolation$test,
         n_years_prior = 1)
  ),
  # The rolling forecast origins. `n_years_prior` is the window length, not the
  # three it was: the null looks back from each held-out year, so with a
  # five-year window the last lead year needs to reach five years back to see
  # any training data at all. At three it would see none, and
  # predict_null_fixed_nn_counts() now stops rather than quietly pooling the
  # whole training set.
  lapply(
    seq_along(temporal_forecasting_folds),
    function(i) {
      fold <- temporal_forecasting_folds[[i]]
      list(experiment = "temporal_forecasting",
           label = paste0("temporal_forecasting_", fold$cut_year),
           fold = names(temporal_forecasting_folds)[i],
           training = fold$training,
           test = fold$test,
           n_years_prior = fold$window)
    }
  )
)

# the 2020 three-year fold was fitted and scored before the origins were
# labelled, under experiment "temporal_forecasting" and fold "all"
folds <- lapply(folds, function(fold) {
  if (is.null(fold$label)) fold$label <- fold$experiment
  fold
})


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
      fold = fold$fold,
      label = fold$label
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
    fold = fold$fold,
    label = fold$label
  )

}


# dynamical model ----------------------------------------------------------

# Each fold is fitted by its own process, via run_one_fold.R, so that greta's
# warmup and sampling progress goes to that fold's log and can be watched while
# it runs. Two folds at a time, four chains each.
#
# Why this split: TensorFlow vectorises chains into a single op rather than
# running one per core, and that op scales poorly beyond about four threads, so
# a fold confined to four threads loses little. Four chains rather than two
# because greta pools information across chains when adapting during warmup, and
# two chains adapt poorly — the Kenya fold reached Rhat 7.6 that way.
#
# Folds whose draws are already on disk are skipped, so the run resumes.
#
# Two folds at a time, four threads each. Measured on this machine: two
# concurrent at four threads is 62 h per fold, three at three threads is 95 h,
# so with only the two forecast origins left to fit, two concurrent is the whole
# run in about 62 h. Memory is the real constraint — three concurrent took
# available memory down to 3 GB — and the five-year folds ask for four times as
# many predictions as the three-year one did, which is why fit_fold() now thins
# the stored draws and splits the calculate() in two.
n_concurrent <- 2
chains_per_fold <- 4
threads_per_fold <- 4

log_dir <- "outputs/cv_logs"
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

pending <- Filter(
  function(fold) {
    !file.exists(file.path(
      draws_dir,
      sprintf("dynamical__%s__%s.rds", fold$experiment, fold$fold)))
  },
  folds
)

cat(sprintf("\n%i folds to fit, %i at a time, %i chains and %i threads each\n",
            length(pending), n_concurrent, chains_per_fold, threads_per_fold))
cat(sprintf("progress logs: %s/\n\n", log_dir))

running <- list()

launch <- function(fold) {
  log_file <- file.path(log_dir,
                        sprintf("%s__%s.log", fold$experiment, fold$fold))
  cat(sprintf("%s | launching %s / %s -> %s\n",
              format(Sys.time(), "%H:%M:%S"),
              fold$experiment, fold$fold, log_file))
  process <- processx::process$new(
    "Rscript",
    c("R/run_one_fold.R", fold$experiment, fold$fold,
      as.character(chains_per_fold), as.character(threads_per_fold)),
    stdout = log_file,
    stderr = "2>&1"
  )
  list(fold = fold, process = process, log = log_file)
}

queue <- pending

while (length(queue) > 0 || length(running) > 0) {

  # start jobs while there is room
  while (length(running) < n_concurrent && length(queue) > 0) {
    running <- c(running, list(launch(queue[[1]])))
    queue <- queue[-1]
  }

  Sys.sleep(60)

  # reap anything that has finished
  still_running <- list()
  for (job in running) {
    if (job$process$is_alive()) {
      still_running <- c(still_running, list(job))
    } else {
      status <- job$process$get_exit_status()
      cat(sprintf("%s | finished %s / %s (exit %s)\n",
                  format(Sys.time(), "%H:%M:%S"),
                  job$fold$experiment, job$fold$fold, status))
      if (!identical(status, 0L)) {
        cat(sprintf("  NON-ZERO EXIT: see %s\n", job$log))
      }
    }
  }
  running <- still_running

}

cat("\nall folds finished\n")
