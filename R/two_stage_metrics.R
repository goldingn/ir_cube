# The two-stage model (#21) against the dynamical model: paired differences
# per experiment, and by forecast horizon.
#
#   Rscript R/two_stage_metrics.R
#
# Run after R/validation_metrics.R, which scores every fold in
# outputs/cv_draws, the two-stage model's (two_stage__*.rds) among them, and
# writes the per-record scores this reads (outputs/cv_scores.csv). Nothing is
# rescored here: every model's numbers are #12's, and this adds only the
# comparison with the dynamical model, which is the question for the
# correction, with paired pixel-bootstrap intervals (pixel_bootstrap(), as in
# R/variance_explained.R).
#
# Per group of held-out records and model:
#   elpd       mean beta-binomial log predictive density (higher is better)
#   crps       mean CRPS on the mortality scale (lower is better)
#   mse        mean squared error of the predictive mean (lower is better)
#   cover50, cover95
#              expected coverage of the central 50% and 95% predictive
#              intervals over the PIT randomisation (expected_coverage())
# each with its 95% interval, and its difference from the dynamical model's,
# with the 95% interval of the paired difference and the bootstrap probability
# that the model beats the dynamical model (for coverage: is closer to
# nominal). Variance explained, its ceiling and its paired contrasts, per
# experiment, fold and forecast horizon, are R/variance_explained.R's, on the
# same records.
#
# Groups: each experiment pooled over its folds (spatial_interpolation,
# spatial_blocks, temporal_forecasting_2014, temporal_forecasting_2018), the two
# spatial blocks separately, and the two forecasting origins pooled as
# temporal_forecasting (their test sets overlap; a pixel's records from both
# move together in the bootstrap). By horizon: years ahead of the last training
# year, per origin and pooled.
#
# Writes outputs/two_stage/cv_headline_two_stage.csv and
# outputs/two_stage/cv_horizon_two_stage.csv.

source("R/validation_scoring.R")

set.seed(2026 - 10 - 2)
n_bootstrap <- 2000
output_dir <- "outputs/two_stage"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

scores <- read.csv("outputs/cv_scores.csv", colClasses = c(fold = "character"))
models <- c("dynamical", "two_stage", "nearest_neighbour",
            "nearest_neighbour_oracle", "intercept")
stopifnot(all(models %in% scores$model))

# years ahead of the last training year, for the forecasting folds: the
# training data are year_start < cut, so the cut year itself is one year ahead
scores <- scores %>%
  filter(model %in% models) %>%
  mutate(horizon = ifelse(startsWith(experiment, "temporal_forecasting"),
                          year_start - suppressWarnings(as.integer(fold)) + 1,
                          NA),
         squared_error = (observed - predicted) ^ 2,
         cover50 = expected_coverage(cdf_below, pmf_at, 0.5),
         cover95 = expected_coverage(cdf_below, pmf_at, 0.95))

# one row per held-out record, one column per model and measure. Every model
# scored the fold's records in the same order, which is checked on the record
measures <- c(log_score = "elpd", crps = "crps", squared_error = "mse",
              cover50 = "cover50", cover95 = "cover95")
records <- bind_rows(lapply(split(scores, paste(scores$experiment,
                                                 scores$fold)), function(d) {
  reference <- d %>% filter(model == "dynamical")
  columns <- lapply(models, function(m) {
    x <- d %>% filter(model == m)
    stopifnot(identical(record_key(x), record_key(reference)))
    setNames(x[names(measures)], paste(measures, m, sep = "|"))
  })
  bind_cols(reference %>% select(experiment, fold, cell, horizon),
            columns)
}))


# summaries --------------------------------------------------------------------------

# every measure for every model, from one set of records
measures_for <- function(data) {
  unlist(lapply(models, function(m) {
    value <- vapply(setNames(measures, measures), function(measure) {
      mean(data[[paste(measure, m, sep = "|")]])
    }, numeric(1))
    setNames(value, paste(names(value), m, sep = "|"))
  }))
}

summarise_records <- function(data) {
  point <- measures_for(data)
  replicates <- pixel_bootstrap(data, measures_for, n_bootstrap)
  interval <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
  bind_rows(lapply(models, function(m) {
    row <- tibble(model = m, n = nrow(data), n_pixels = n_distinct(data$cell))
    for (measure in measures) {
      column <- paste(measure, m, sep = "|")
      reference <- paste(measure, "dynamical", sep = "|")
      difference <- replicates[, column] - replicates[, reference]
      nominal <- c(cover50 = 0.5, cover95 = 0.95)[measure]
      better <- switch(
        measure,
        crps = , mse = difference < 0,
        cover50 = , cover95 = abs(replicates[, column] - nominal) <
          abs(replicates[, reference] - nominal),
        difference > 0)
      row[[measure]] <- point[[column]]
      row[[paste0(measure, "_lower")]] <- interval(replicates[, column], 0.025)
      row[[paste0(measure, "_upper")]] <- interval(replicates[, column], 0.975)
      row[[paste0("diff_", measure)]] <- point[[column]] - point[[reference]]
      row[[paste0("diff_", measure, "_lower")]] <- interval(difference, 0.025)
      row[[paste0("diff_", measure, "_upper")]] <- interval(difference, 0.975)
      row[[paste0("prob_better_", measure)]] <- mean(better)
    }
    row
  }))
}

# the four experiments pooled over their folds, the forecasting origins pooled,
# and the two spatial blocks separately
forecasting <- records %>%
  filter(startsWith(experiment, "temporal_forecasting")) %>%
  mutate(experiment = "temporal_forecasting")
blocks <- records %>% filter(experiment == "spatial_blocks")
groups <- c(lapply(split(records, records$experiment), mutate,
                   fold = "pooled"),
            list(mutate(forecasting, fold = "pooled")),
            split(blocks, blocks$fold))
headline <- bind_rows(lapply(groups, function(data) {
  summarise_records(data) %>%
    mutate(experiment = data$experiment[1], fold = data$fold[1], .before = 1)
}))
write.csv(headline, file.path(output_dir, "cv_headline_two_stage.csv"),
          row.names = FALSE)

by_horizon <- bind_rows(records %>% filter(!is.na(horizon)), forecasting)
horizon <- bind_rows(lapply(
  split(by_horizon, paste(by_horizon$experiment, by_horizon$horizon)),
  function(data) {
    if (n_distinct(data$cell) < 2) return(NULL)
    summarise_records(data) %>%
      mutate(experiment = data$experiment[1], horizon = data$horizon[1],
             .before = 1)
  })) %>%
  arrange(experiment, horizon)
write.csv(horizon, file.path(output_dir, "cv_horizon_two_stage.csv"),
          row.names = FALSE)


# report ------------------------------------------------------------------------------

options(width = 200)
show <- function(table, ...) {
  table %>%
    filter(model == "two_stage") %>%
    transmute(..., n,
              elpd = sprintf("%.3f", elpd),
              d_elpd = sprintf("%+.3f [%+.3f, %+.3f]", diff_elpd,
                               diff_elpd_lower, diff_elpd_upper),
              d_crps_e3 = sprintf("%+.1f [%+.1f, %+.1f]", 1e3 * diff_crps,
                                  1e3 * diff_crps_lower, 1e3 * diff_crps_upper),
              cover = sprintf("%.3f / %.3f", cover50, cover95)) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
}
cat("two-stage model, and its difference from the dynamical model:\n")
show(headline, experiment, fold)
cat("\nby forecast horizon:\n")
show(horizon, experiment, horizon)
