# Score the forecasting experiment on the change in mortality, not the level.
#
# With the training set leak fixed, the forecasting experiment is still
# substantially a spatial test: most held-out site-years sit at sites with
# training data one to three years earlier, so a local method gets the level
# nearly free and the comparison is decided by local structure rather than by
# trend. Scoring the change differences the site level out and leaves the local
# slope, which is what the dynamical model claims to know (#12 review 5.1).
#
# For each group g with data in both windows — a before window B of the same
# length as the holdout, immediately preceding the cut, and the holdout window
# H:
#
#   delta_obs(g)  = sum_H died / sum_H tested  -  sum_B died / sum_B tested
#   delta_pred(g) = mean over draws of ( weighted mean of p over H
#                                        - weighted mean of p over B )
#
# The target is a difference of two empirical proportions, with no model
# assumptions in it. The baseline is "no change", which is what the nearest
# neighbour null predicts by construction, so the nulls need no separate
# treatment here.
#
# The floor is the sum of the two windows' irreducible variances, since the
# windows are independent given the fractions; noise_floor_var_pooled() gives
# each one.
#
# Run at two scales. Cell-insecticide is the most local and the thinnest, and
# the floor handles that honestly. Country-insecticide pools many more assays,
# so the floor falls and the comparison is almost purely about the population
# fraction.

source("R/validation_functions.R")

suppressMessages({
  library(dplyr)
  library(tidyr)
})

# Every forecasting origin on disk, scored separately. The design is a rolling
# origin — five-year windows cut at 2014 and 2018, plus the superseded 2020
# three-year fold kept as a supplementary observation — and the origins must not
# be pooled: their holdout windows have true rates of decline differing by a
# factor of two, which is the contrast the test is built on.
fold_files <- sort(list.files(
  "outputs/cv_draws",
  pattern = "^dynamical__temporal_forecasting__.*\\.rds$",
  full.names = TRUE
))
if (length(fold_files) == 0) {
  stop("no forecasting folds in outputs/cv_draws")
}

rho_external <- read.csv("outputs/bioassay_rho.csv")
rho_for_class <- function(insecticide_class) {
  index <- match(insecticide_class, rho_external$insecticide_class)
  pooled <- rho_external$rho[rho_external$insecticide_class == "all"]
  ifelse(is.na(index), pooled, rho_external$rho[index])
}

# thin to a common set of draws; the change is a smooth functional and does not
# need twenty thousand of them
max_draws <- 2000

# the pooled observed proportion and the size-weighted mean of p, per group and
# per window
window_summary <- function(data, p_draws, group) {
  index <- split(seq_len(nrow(data)), group)
  weights <- data$mosquito_number
  list(
    groups = names(index),
    died = vapply(index, function(i) sum(data$died[i]), numeric(1)),
    tested = vapply(index, function(i) sum(weights[i]), numeric(1)),
    assays = vapply(index, length, numeric(1)),
    # n_draws x n_groups, the weighted mean of p in that window
    p = vapply(index, function(i) {
      as.numeric(p_draws[, i, drop = FALSE] %*% weights[i]) / sum(weights[i])
    }, numeric(nrow(p_draws))),
    insecticide_class = vapply(index,
                               function(i) data$insecticide_class[i[1]],
                               character(1)),
    # the assay sizes, kept so the floor can be computed per group
    sizes = lapply(index, function(i) weights[i]),
    counts = lapply(index, function(i) data$died[i])
  )
}

score_scale <- function(scale, holdout, before, p_holdout, p_before) {

  grouping <- function(data) {
    switch(
      scale,
      cell = paste(data$cell, data$insecticide_type),
      country = paste(data$country_name, data$insecticide_type)
    )
  }

  h <- window_summary(holdout, p_holdout, grouping(holdout))
  b <- window_summary(before, p_before, grouping(before))

  shared <- intersect(h$groups, b$groups)
  if (length(shared) == 0) {
    return(NULL)
  }
  hi <- match(shared, h$groups)
  bi <- match(shared, b$groups)

  rho <- rho_for_class(h$insecticide_class[hi])

  delta_observed <- h$died[hi] / h$tested[hi] - b$died[bi] / b$tested[bi]
  delta_draws <- h$p[, hi, drop = FALSE] - b$p[, bi, drop = FALSE]
  delta_predicted <- colMeans(delta_draws)

  # the irreducible variance of the observed change: the two windows are
  # independent given the fractions, so their variances add
  floor_variance <- vapply(seq_along(shared), function(k) {
    noise_floor_var_pooled(h$counts[[hi[k]]], h$sizes[[hi[k]]], rho[k]) +
      noise_floor_var_pooled(b$counts[[bi[k]]], b$sizes[[bi[k]]], rho[k])
  }, numeric(1))

  data.frame(
    scale = scale,
    group = shared,
    stringsAsFactors = FALSE,
    insecticide_class = h$insecticide_class[hi],
    assays_before = b$assays[bi],
    assays_holdout = h$assays[hi],
    delta_observed = delta_observed,
    delta_predicted = delta_predicted,
    delta_sd = apply(delta_draws, 2, sd),
    floor_variance = floor_variance,
    row.names = NULL
  )

}

# score one fitted forecasting fold at both scales
score_fold <- function(fold_file) {

  fold <- readRDS(fold_file)
  if (is.null(fold$p_draws_before)) {
    warning(basename(fold_file), " carries no before-window predictions; it ",
            "predates the change-based score and cannot be used for it ",
            "without refitting. Skipping.")
    return(NULL)
  }

  holdout <- fold$test_df
  before <- fold$before_df
  p_holdout <- fold$p_draws
  p_before <- fold$p_draws_before
  stopifnot(ncol(p_holdout) == nrow(holdout),
            ncol(p_before) == nrow(before),
            nrow(p_holdout) == nrow(p_before))

  if (nrow(p_holdout) > max_draws) {
    keep <- round(seq(1, nrow(p_holdout), length.out = max_draws))
    p_holdout <- p_holdout[keep, , drop = FALSE]
    p_before <- p_before[keep, , drop = FALSE]
  }

  out <- bind_rows(lapply(c("cell", "country"), score_scale,
                          holdout = holdout, before = before,
                          p_holdout = p_holdout, p_before = p_before))
  if (is.null(out) || nrow(out) == 0) {
    return(NULL)
  }

  out %>%
    mutate(experiment = fold$experiment,
           cut_year = as.integer(sub("^temporal_forecasting_", "",
                                     fold$experiment)),
           holdout_years = paste(range(holdout$year_start), collapse = "-"),
           before_years = paste(range(before$year_start), collapse = "-"),
           .before = everything())
}

changes <- bind_rows(lapply(fold_files, score_fold))
if (nrow(changes) == 0) {
  stop("no forecasting fold carries before-window predictions")
}
changes <- changes %>% arrange(cut_year, scale)
write.csv(changes, "outputs/cv_change.csv", row.names = FALSE)

# summarise: mean squared error of the predicted change against "no change",
# both measured above the floor, and the sign test
summary_table <- changes %>%
  filter(!is.na(floor_variance)) %>%
  group_by(cut_year, holdout_years, before_years, scale) %>%
  summarise(
    groups = n(),
    assays = sum(assays_before + assays_holdout),
    mean_observed_change = mean(delta_observed),
    mean_predicted_change = mean(delta_predicted),
    floor = mean(floor_variance),
    mse_model = mean((delta_observed - delta_predicted) ^ 2),
    mse_no_change = mean(delta_observed ^ 2),
    .groups = "drop"
  ) %>%
  mutate(
    excess_model = mse_model - floor,
    excess_no_change = mse_no_change - floor,
    skill = mse_skill(mse_model, mse_no_change, floor)
  )

# direction of change, among groups whose observed change is larger than its
# own noise standard deviation
sign_table <- changes %>%
  filter(!is.na(floor_variance),
         abs(delta_observed) > sqrt(floor_variance)) %>%
  group_by(cut_year, scale) %>%
  summarise(
    groups = n(),
    correct_direction = mean(sign(delta_predicted) == sign(delta_observed)),
    mean_absolute_observed = mean(abs(delta_observed)),
    .groups = "drop"
  )

write.csv(summary_table, "outputs/cv_change_summary.csv", row.names = FALSE)
write.csv(sign_table, "outputs/cv_change_direction.csv", row.names = FALSE)

cat("\nchange in mortality between the window before each cut and the",
    "holdout window:\n")
print(as.data.frame(summary_table %>%
        mutate(across(where(is.numeric), ~ round(.x, 4)))))

cat("\ndirection of change, groups with an observed change above its noise:\n")
print(as.data.frame(sign_table %>%
        mutate(across(where(is.numeric), ~ round(.x, 3)))))
