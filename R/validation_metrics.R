# Score the saved cross-validation posterior predictive draws.
#
# Reads everything written by run_validation_folds.R and writes tidy tables of
# per-record scores and per-experiment summaries, which the figures then read.
# Keeping scoring separate from fitting means metrics can be revised without
# refitting anything (#10).
#
# Three questions are asked of each model, one measure each:
#
#   is the predictive distribution the right width and shape?
#       coverage of central predictive intervals, with the Cramer-von Mises
#       statistic on the randomised PIT values as a single scalar summary
#   is the whole distribution close to the data?
#       CRPS, in mortality units, and mean squared error against the noise floor
#   is it right on average, and conditional on the prediction?
#       reliability bins, and the mean PIT
#
# Coverage, mean PIT and the Cramer-von Mises statistic are all functionals of
# one PIT distribution, so only these are kept: the Kolmogorov-Smirnov
# statistic added a fourth view of the same object, and a null band for the
# Cramer-von Mises statistic assumed independent PIT values, which does not hold
# for records scored against a shared posterior (#12 review).
#
# Every model is scored at the same overdispersion, the external
# replicate-based estimate, rather than at its own fitted value. Letting each
# model choose made coverage a comparison of dispersion rather than of
# prediction: a model can look calibrated by being vague, and the intercept null
# reached nominal coverage that way. The dispersion each model's own residuals
# imply is reported separately, in cv_rho_comparison.csv.
#
# Scores are also computed on pooled groups of assays. A single bioassay is a
# noisy measurement of the population fraction the model is predicting, so
# aggregation is what brings the comparison to bear on that quantity: the
# reference in every case is the model's own aggregated predictive
# distribution, so heterogeneity in the true fraction within a group is carried
# by the model's predictions rather than assumed away. Only the country-year
# rung is reported; pooling within pixel, year and insecticide averaged 1.3
# assays per group, so it was indistinguishable from the unpooled scores.

source("R/validation_functions.R")

suppressMessages({
  library(dplyr)
  library(tidyr)
})

draws_dir <- "outputs/cv_draws"
n_pit_reps <- 100
coverage_levels <- seq(0.1, 0.95, by = 0.05)

# Folds fitted by MCMC hold up to 20,000 draws (4 chains x 5,000). The scoring
# functions build an n_draws x (died + 1) matrix per observation, so that is
# twenty times the work of the 1,000 draws the null models carry, for no gain:
# at roughly 50 draws per effective sample, every tenth draw retains almost all
# the information. Thinning is applied to the model folds only; the null models
# are analytic and already independent.
max_draws <- 2000

thin_draws <- function(x, maximum = max_draws) {
  if (nrow(x) <= maximum) return(x)
  keep <- round(seq(1, nrow(x), length.out = maximum))
  x[keep, , drop = FALSE]
}

set.seed(2026 - 8 - 31)

# externally estimated overdispersion, from replicate bioassays in the same
# pixel, year and insecticide (estimate_bioassay_rho.R). This sets the noise
# floor, and is independent of any of the models being scored
rho_external <- read.csv("outputs/bioassay_rho.csv")
rho_for_class <- function(insecticide_class) {
  index <- match(insecticide_class, rho_external$insecticide_class)
  # fall back to the pooled estimate for any class not fitted separately
  pooled <- rho_external$rho[rho_external$insecticide_class == "all"]
  ifelse(is.na(index), pooled, rho_external$rho[index])
}

# Every model is scored at that external estimate, so the only thing differing
# between models is the posterior on the population fraction. Set this to FALSE
# to score each model at its own fitted overdispersion instead, which is the
# sensitivity behind the choice: the dynamical model's posterior on p was
# fitted jointly with its own rho, so the swap is not perfectly clean without a
# refit, and it is worth reporting both.
score_at_external_rho <- TRUE

files <- list.files(draws_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0) {
  stop("no draws found in ", draws_dir, "; run R/run_validation_folds.R first")
}

cat(sprintf("scoring %i saved folds\n", length(files)))


# per record ---------------------------------------------------------------

score_fold <- function(file) {

  fold <- readRDS(file)
  test <- fold$test_df

  # predictions come from the saved object. For model folds these were produced
  # by greta's calculate(values = draws) in MCMC order; for the null models they
  # are analytic. Either way they are not recomputed here
  p_draws <- thin_draws(fold$p_draws)

  # the overdispersion each model's own fit implies, kept for the diagnostic
  # table: a posterior for the dynamical model, a single fitted value for the
  # nulls, and absent from folds saved before the two were separated
  rho_fitted <- if (!is.null(fold$rho_class_draws)) {
    # saved compactly as draws by insecticide class, expanded here
    thin_draws(fold$rho_class_draws)[, fold$class_id, drop = FALSE]
  } else if (!is.null(fold$rho_draws)) {
    thin_draws(fold$rho_draws)
  } else {
    matrix(fold$rho_implied, nrow = nrow(p_draws), ncol = nrow(test))
  }

  rho_scoring <- rho_for_class(test$insecticide_class)
  rho_draws <- if (score_at_external_rho) {
    matrix(rho_scoring, nrow = nrow(p_draws), ncol = nrow(test), byrow = TRUE)
  } else {
    rho_fitted
  }

  summary <- ppd_summary(test$died,
                         test$mosquito_number,
                         p_draws,
                         rho_draws)

  pit <- ppd_pit(summary, n_rep = n_pit_reps)
  sims <- ppd_simulate(test$mosquito_number, p_draws, rho_draws)

  scores <- summary %>%
    mutate(
      model = fold$model,
      experiment = fold$experiment,
      fold = fold$fold,
      insecticide_type = test$insecticide_type,
      insecticide_class = test$insecticide_class,
      country_name = test$country_name,
      year_start = test$year_start,
      cell = test$cell,
      # One randomisation replicate, not the mean of them. Averaging over
      # replicates converges to the mid-P value `cdf_below + 0.5 * pmf_at`,
      # which is not uniform under calibration for discrete data: with ~30% of
      # held-out assays at 100% mortality a perfectly calibrated model reads
      # 0.969 at nominal 0.95. The uniformity statistics use the full matrix
      # and were never affected; this column feeds the figures (#12 review)
      pit = pit[, 1],
      crps = ppd_crps(test$died, test$mosquito_number, sims),
      rho_external = rho_scoring,
      rho_fitted = colMeans(rho_fitted),
      .before = everything()
    )

  list(scores = scores,
       pit = pit,
       sims = sims,
       fold = fold)

}

scored <- lapply(files, score_fold)
names(scored) <- basename(files)

all_scores <- bind_rows(lapply(scored, `[[`, "scores"))

write.csv(all_scores, "outputs/cv_scores.csv", row.names = FALSE)


# per experiment -----------------------------------------------------------

# summarise one model in one experiment, pooling its folds
summarise_experiment <- function(scores, pit_list) {

  pit <- do.call(rbind, pit_list)
  n_obs <- nrow(scores)

  coverage <- coverage_curve(pit, levels = c(0.5, 0.95))
  floor_mse <- noise_floor_mse(scores$died,
                               scores$mosquito_number,
                               scores$rho_external)

  data.frame(
    n = n_obs,
    mean_pit = mean(pit),
    coverage_50 = coverage$empirical[1],
    coverage_95 = coverage$empirical[2],
    cvm = pit_statistic(pit, cvm_stat),
    crps = mean(scores$crps),
    elpd = mean(scores$log_score),
    bias = mean(scores$predicted - scores$observed),
    mse = mean((scores$observed - scores$predicted) ^ 2),
    mse_floor = floor_mse
  )

}

keys <- bind_rows(lapply(scored, function(x) {
  data.frame(model = x$fold$model, experiment = x$fold$experiment)
}))

summaries <- lapply(
  split(seq_along(scored), paste(keys$model, keys$experiment)),
  function(index) {
    summarise_experiment(
      bind_rows(lapply(scored[index], `[[`, "scores")),
      lapply(scored[index], `[[`, "pit")
    ) %>%
      mutate(model = keys$model[index[1]],
             experiment = keys$experiment[index[1]],
             .before = everything())
  }
)
summaries <- bind_rows(summaries)

# Variance in the population fraction that each model explains, anchored on the
# intercept null and on the noise floor: 0 is the no-information baseline, 1 is
# as good as bioassay noise allows. The anchor was previously the nearest
# neighbour null, which pinned an informative baseline at zero by construction
# and hid that it is itself worse than a global per-insecticide mean under
# spatial extrapolation (#12 review).
#
# `excess` is reported alongside every ratio: it is mean squared error above the
# noise floor in absolute mortality-squared units, so the conclusion does not
# rest entirely on the floor. `rms_p` is its square root, an error in the
# population fraction itself.
summaries <- summaries %>%
  group_by(experiment) %>%
  mutate(
    excess = mse - mse_floor,
    rms_p = sqrt(pmax(excess, 0)),
    mse_null = mse[model == "intercept"],
    skill = mse_skill(mse, mse_null, mse_floor)
  ) %>%
  ungroup()

write.csv(summaries, "outputs/cv_summary.csv", row.names = FALSE)

cat("\nsummary by experiment and model:\n")
print(summaries %>%
        select(experiment, model, n, coverage_95, mean_pit, crps, mse,
               mse_floor, excess, rms_p, skill, cvm) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
        as.data.frame())


# coverage curves ----------------------------------------------------------

coverage_curves <- lapply(
  split(seq_along(scored), paste(keys$model, keys$experiment)),
  function(index) {
    pit <- do.call(rbind, lapply(scored[index], `[[`, "pit"))
    coverage_curve(pit, levels = coverage_levels) %>%
      mutate(model = keys$model[index[1]],
             experiment = keys$experiment[index[1]],
             .before = everything())
  }
)
coverage_curves <- bind_rows(coverage_curves)
write.csv(coverage_curves, "outputs/cv_coverage.csv", row.names = FALSE)


# reliability --------------------------------------------------------------

reliability <- all_scores %>%
  group_by(model, experiment) %>%
  group_modify(~ {
    bins <- reliability_bins(.x$predicted, .x$observed, n_bins = 10)
    # the scatter a perfect model would still show in each bin, from the
    # external overdispersion, the assay sizes and the number of assays pooled
    bins$envelope <- reliability_envelope(
      p = bins$predicted,
      mosquito_number = median(.x$mosquito_number),
      k = bins$n,
      rho = mean(.x$rho_external)
    )
    bins
  }) %>%
  ungroup()

write.csv(reliability, "outputs/cv_reliability.csv", row.names = FALSE)


# aggregated scores --------------------------------------------------------

# Assays are pooled by country, year and insecticide: 17-18 assays per group,
# so assay noise falls roughly seventeen-fold and the comparison is nearly
# purely about the population fraction, at the cost of testing an aggregate
# rather than any single pixel. Pooling within pixel, year and insecticide was
# also reported, and dropped: it averaged 1.3 assays per group, so it was the
# unpooled comparison under another name (#12 review).
aggregate_fold <- function(entry, grouping) {

  test <- entry$fold$test_df
  group <- switch(
    grouping,
    country_year = paste(test$country_name, test$year_start,
                         test$insecticide_type)
  )

  ppd_aggregate(test$died, test$mosquito_number, group, entry$sims) %>%
    mutate(model = entry$fold$model,
           experiment = entry$fold$experiment,
           grouping = grouping,
           .before = everything())

}

aggregated <- bind_rows(
  unname(lapply(scored, aggregate_fold, grouping = "country_year"))
)

write.csv(aggregated, "outputs/cv_aggregate.csv", row.names = FALSE)

aggregate_summary <- aggregated %>%
  group_by(grouping, experiment, model) %>%
  summarise(
    groups = n(),
    mean_assays = mean(n_assays),
    coverage_95 = mean(observed >= lower & observed <= upper),
    mean_pit = mean(pit),
    bias = mean(predicted - observed),
    rmse = rmse(observed, predicted),
    .groups = "drop"
  )

cat("\naggregated calibration:\n")
print(aggregate_summary %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
        as.data.frame())

write.csv(aggregate_summary, "outputs/cv_aggregate_summary.csv",
          row.names = FALSE)


# the model's overdispersion against the external estimate ------------------

# a fitted overdispersion larger than the replicate-based estimate would mean
# the model is absorbing process misfit into the observation process, which
# would also show as over-coverage
# sampling diagnostics per fold, so that the convergence caveat travels with
# the results. Null model folds are analytic and have none
sampling_diagnostics <- bind_rows(lapply(scored, function(entry) {
  fold <- entry$fold
  if (is.null(fold$ess_p)) return(NULL)
  data.frame(
    model = fold$model,
    experiment = fold$experiment,
    fold = fold$fold,
    n_chains = fold$n_chains,
    n_sampled = fold$n_sampled,
    ess_p_median = median(fold$ess_p, na.rm = TRUE),
    ess_p_min = min(fold$ess_p, na.rm = TRUE),
    ess_rho_median = median(fold$ess_rho, na.rm = TRUE),
    rhat_worst = max(fold$convergence[, 1], na.rm = TRUE),
    rhat_above_1.01 = sum(fold$convergence[, 1] > 1.01, na.rm = TRUE)
  )
}))

if (nrow(sampling_diagnostics) > 0) {
  write.csv(sampling_diagnostics, "outputs/cv_sampling_diagnostics.csv",
            row.names = FALSE)
  cat("\nsampling diagnostics:\n")
  print(sampling_diagnostics %>%
          mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
          as.data.frame())
}


rho_comparison <- bind_rows(lapply(scored, function(entry) {
  entry$scores %>%
    select(model, experiment, fold, insecticide_class, rho_fitted)
})) %>%
  group_by(model, experiment, insecticide_class) %>%
  summarise(rho_fitted = mean(rho_fitted), .groups = "drop") %>%
  mutate(rho_external = rho_for_class(insecticide_class))

write.csv(rho_comparison, "outputs/cv_rho_comparison.csv", row.names = FALSE)

cat("\nfitted against externally estimated overdispersion:\n")
print(rho_comparison %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
        as.data.frame())


# per fold and per year ----------------------------------------------------

# The pooled numbers hide which folds carry the result, and master reported a
# per-country and a per-lead-year breakdown that the first version of this
# pipeline dropped. `excess` is again mean squared error above the noise floor,
# and `explained` is the share of the intercept null's excess that this model
# removes, computed within each group so that groups of differing difficulty are
# not compared on a common denominator (#12 review).
by_group <- function(scores, ...) {
  scores %>%
    group_by(experiment, ..., model) %>%
    summarise(
      n = n(),
      mean_observed = mean(observed),
      mean_predicted = mean(predicted),
      bias = mean(predicted - observed),
      mean_pit = mean(pit),
      coverage_95 = mean(pit > 0.025 & pit < 0.975),
      crps = mean(crps),
      mse = mean((observed - predicted) ^ 2),
      mse_floor = noise_floor_mse(died, mosquito_number, rho_external),
      .groups = "drop"
    ) %>%
    group_by(experiment, ...) %>%
    mutate(excess = mse - mse_floor,
           explained = 1 - excess / excess[model == "intercept"]) %>%
    ungroup()
}

by_fold <- by_group(all_scores, fold)
write.csv(by_fold, "outputs/cv_by_fold.csv", row.names = FALSE)

cat("\nby fold:\n")
print(by_fold %>%
        select(experiment, fold, model, n, bias, coverage_95, excess,
               explained) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
        as.data.frame())

by_year <- by_group(
  all_scores %>% filter(experiment == "temporal_forecasting"),
  year_start
)
write.csv(by_year, "outputs/cv_by_year.csv", row.names = FALSE)

if (nrow(by_year) > 0) {
  cat("\nforecasting, by lead year:\n")
  print(by_year %>%
          select(year_start, model, n, mean_observed, mean_predicted, bias,
                 mean_pit, excess, explained) %>%
          mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
          as.data.frame())
}
