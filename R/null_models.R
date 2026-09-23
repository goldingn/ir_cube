# Null model definitions used in the cross-validation experiments: the
# insecticide-type intercept model and the nearest neighbour heuristic.
#
# Moved out of predictive_validation.R so that both the existing point-prediction
# scoring and the posterior predictive scoring (#10) use one definition of each
# null model. The lower half of this file adds predictive distributions for the
# nulls, because a proper scoring rule needs a distribution rather than a point
# from every candidate. Every model is scored at the same external,
# replicate-based overdispersion (see validation_metrics.R); what each model's
# own residuals imply is reported separately as a diagnostic.

prop <- function(died, tested) {
  died / tested
}

# approximate non-zero and non-one proportions by applying the empirical logit
# transform to the data and then the inverse logit to provide a proportion
emplog_prop <- function(died, mosquito_number) {
  emplog <- log((died + 0.5) / (mosquito_number - died + 0.5))
  plogis(emplog)
}
# Define the nearest neighbour null model: For each point in the training data,
# average over the X nearest datapoints from the current and previous year

# given vectors of latitude, longitude, year, and insecticide type for test
# data, and a tibble of training data, return a vector of predictions of the
# susceptibility fraction from a weighted average of the `n_nearest_neighbours`
# nearest points in the training data of that insecticide type, and in the same
# year or up to `n_years_prior` earlier years
predict_null_fixed_nn <- function(latitude,
                                  longitude,
                                  year,
                                  insecticide_type,
                                  training_data,
                                  n_nearest_neighbours,
                                  n_years_prior = 1) {

  counts <- predict_null_fixed_nn_counts(
    latitude = latitude,
    longitude = longitude,
    year = year,
    insecticide_type = insecticide_type,
    training_data = training_data,
    n_nearest_neighbours = n_nearest_neighbours,
    n_years_prior = n_years_prior
  )

  emplog_prop(counts$total_died, counts$total_tested)

}

# predictive distributions for the null models ------------------------------

# maximum likelihood estimate of the observation overdispersion implied by a
# set of predictions. This is reported as a diagnostic - a model whose implied
# overdispersion greatly exceeds the replicate-based estimate is absorbing
# process misfit into the observation process - but it is no longer the
# dispersion the model is scored under. Letting each model choose its own rho
# made coverage a comparison of dispersion rather than of prediction: the
# intercept null reached nominal coverage by inflating rho to 0.45 against an
# external estimate of 0.12-0.22 (#12 review)
fit_rho_given_predictions <- function(died, mosquito_number, predicted,
                                      interval = c(1e-4, 0.9)) {
  negative_log_likelihood <- function(rho) {
    -sum(dbetabinom(died, mosquito_number, predicted, rho, log = TRUE))
  }
  optimise(negative_log_likelihood, interval = interval)$minimum
}

# As predict_null_fixed_nn(), but returning the pooled counts over the
# neighbours rather than a single proportion, so that uncertainty in the
# predicted fraction can be represented by a beta posterior on those counts.
# Neighbour selection matches predict_null_fixed_nn() exactly
predict_null_fixed_nn_counts <- function(latitude,
                                         longitude,
                                         year,
                                         insecticide_type,
                                         training_data,
                                         n_nearest_neighbours,
                                         n_years_prior = 1) {

  # A missing row in outputs/optimal_nn.csv returns numeric(0) from the lookup
  # in run_validation_folds.R, and every step below then degrades silently:
  # sort(x)[integer(0)] is numeric(0), `<= numeric(0)` is logical(0),
  # which(logical(0)) is integer(0), and sum() of no elements is 0. The result
  # is a pooled count of 0 died out of 0 tested for every held-out record, whose
  # Jeffreys posterior is Beta(0.5, 0.5) — so the "nearest neighbour null"
  # becomes a random number generator with no data in it, and nothing errors.
  # This is exactly what happened to the two spatial block folds, whose
  # experiment name had no row in that file: their reported score of -0.86 was
  # the score of the prior, not of a nearest neighbour model.
  stopifnot(
    is.numeric(n_nearest_neighbours),
    length(n_nearest_neighbours) == 1,
    is.finite(n_nearest_neighbours),
    n_nearest_neighbours >= 1
  )

  training_coords <- as.matrix(training_data[, c("longitude", "latitude")])
  test_coords <- cbind(longitude, latitude)

  dists <- fields::rdist.earth(test_coords, training_coords,
                               miles = FALSE)

  # The year window is anchored at prediction time, not at the record's own
  # year: the most recent `n_years_prior` + 1 years of data that exist when the
  # prediction is made. For the spatial experiments training spans every year,
  # so the anchor is the record's own year and this is the familiar
  # {year, year - 1}. For a forecast the training data stop at the cut, so the
  # anchor is the last training year and every held-out record draws on the same
  # window — which is the situation a person forecasting from that cut is
  # actually in.
  #
  # Anchoring on the record's own year instead made the null weaker the further
  # ahead it had to predict, for an indexing reason rather than an information
  # one: a record one year past the cut saw the whole window, one five years past
  # saw a single year of it, or none at all.
  anchor <- pmin(year, max(training_data$year_start))
  year_diff <- -1 * seq(0, n_years_prior)

  n_test <- nrow(test_coords)
  total_died <- rep(NA_real_, n_test)
  total_tested <- rep(NA_real_, n_test)

  for (i in seq_len(n_test)) {

    distance_vec <- dists[i, ]
    valid_years <- anchor[i] + year_diff
    valid <- training_data$year_start %in% valid_years &
      training_data$insecticide_type == insecticide_type[i]
    masked_distance_vec <- ifelse(valid, distance_vec, Inf)

    # If no training record matches this record's insecticide and year window,
    # every masked distance is Inf, the threshold is Inf, and `<= threshold`
    # then selects the entire training set — every insecticide, every year. That
    # is not a nearest neighbour prediction, it is a global mean wearing one,
    # and it is silent. It bites when the holdout runs further past the cut than
    # `n_years_prior` reaches back: with a five-year forecast window and
    # `n_years_prior = 3`, the fourth and fifth lead years have no valid
    # training year at all. Fail here instead, so the caller has to set
    # `n_years_prior` to at least the window length.
    if (!any(valid)) {
      stop("no training records in ", paste(range(valid_years), collapse = "-"),
           " for ", insecticide_type[i], ": the nearest neighbour null has ",
           "nothing to predict from")
    }

    threshold_distance <- sort(masked_distance_vec,
                               decreasing = FALSE)[n_nearest_neighbours]
    nearest <- which(masked_distance_vec <= threshold_distance)

    total_died[i] <- sum(training_data$died[nearest])
    total_tested[i] <- sum(training_data$mosquito_number[nearest])

  }

  data.frame(total_died = total_died,
             total_tested = total_tested)

}

# posterior draws of the predicted fraction from pooled counts, under a
# Jeffreys beta prior
pooled_count_draws <- function(total_died, total_tested, n_draws) {
  n_obs <- length(total_died)
  draws <- rbeta(n_draws * n_obs,
                 rep(total_died + 0.5, each = n_draws),
                 rep(total_tested - total_died + 0.5, each = n_draws))
  matrix(draws, nrow = n_draws, ncol = n_obs)
}

# Predictive distribution of the insecticide-type intercept null: the fraction
# for each type has a beta posterior from the pooled training counts for that
# type. `rho_implied` is the overdispersion that best explains this model's own
# training residuals, returned for the diagnostic table only; scoring uses the
# external estimate
intercept_null_draws <- function(training_data, test_data, n_draws = 1000) {

  pooled <- aggregate(
    cbind(died, mosquito_number) ~ insecticide_type,
    data = training_data,
    FUN = sum
  )

  index <- match(test_data$insecticide_type, pooled$insecticide_type)
  p_draws <- pooled_count_draws(pooled$died[index],
                                pooled$mosquito_number[index],
                                n_draws)

  training_index <- match(training_data$insecticide_type,
                          pooled$insecticide_type)
  training_predicted <- pooled$died[training_index] /
    pooled$mosquito_number[training_index]
  rho_implied <- fit_rho_given_predictions(training_data$died,
                                           training_data$mosquito_number,
                                           training_predicted)

  list(p_draws = p_draws,
       rho_implied = rho_implied,
       test_df = test_data)

}

# Predictive distribution of the nearest neighbour null. `n_neighbours` is the
# value already selected by grid search on an internal holdout, and
# `rho_implied` - a diagnostic only, as for the intercept null - is fitted on a
# further internal holdout, so neither quantity is tuned on the test fold
# Predictive distribution of the nearest neighbour null.
#
# This is not a competitor model to be optimised, it is a stand-in for what a
# person would do without a model: read a site's value off the nearby recent
# surveys. Tuning it separately for each experiment would build a series of new
# models to compare against the one model under analysis, so it is not tuned at
# all. Two specifications are reported instead, and neither involves a chosen
# value:
#
#   the practice baseline - `n_neighbours = 1`, `n_years_prior = 1`: the single
#     nearest record in the most recent two years available at prediction time.
#     Short recency is justified a priori, since anyone running this
#     surveillance knows resistance is moving fast, and the k-curves agree.
#
#   the oracle bound - nn_oracle_draws() below: the same method at whichever
#     number of neighbours minimises its own error on the held-out records. It
#     answers a different question, how well any such rule could have done in
#     this test, and a bound is properly per-test.
#
# This replaces a tuned `n_neighbours` read from outputs/optimal_nn.csv. That
# was chosen on an internal test set of 100 records sampled at random from the
# training data, which sat a median of 0 km from their nearest usable neighbour
# - 26% of them at the same site - while the held-out records sit at 36 km
# (interpolation), 139 km (blocks) and 194 km (countries). Optimal k rises with
# the distance that has to be reached, so tuning at 0 km chose it too small and
# handicapped the baseline. The lookup also had no row for the block folds,
# which silently reduced the null there to its prior (#12).
nn_null_draws <- function(training_data, test_data, n_neighbours = 1,
                          n_years_prior = 1, n_draws = 1000,
                          holdout_size = 100, seed = 111) {

  counts <- predict_null_fixed_nn_counts(
    latitude = test_data$latitude,
    longitude = test_data$longitude,
    year = test_data$year_start,
    insecticide_type = test_data$insecticide_type,
    training_data = training_data,
    n_nearest_neighbours = n_neighbours,
    n_years_prior = n_years_prior
  )

  p_draws <- pooled_count_draws(counts$total_died,
                                counts$total_tested,
                                n_draws)

  # Fit the implied overdispersion on an internal holdout from the training
  # data. The seed is fixed so the holdout is reproducible, and the caller's
  # random stream is restored afterwards rather than left reset (#12 review)
  if (exists(".Random.seed", envir = globalenv())) {
    caller_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", caller_seed, envir = globalenv()))
  }
  set.seed(seed)
  holdout <- sample(nrow(training_data), min(holdout_size, nrow(training_data)))
  holdout_data <- training_data[holdout, ]
  remainder <- training_data[-holdout, ]

  holdout_counts <- predict_null_fixed_nn_counts(
    latitude = holdout_data$latitude,
    longitude = holdout_data$longitude,
    year = holdout_data$year_start,
    insecticide_type = holdout_data$insecticide_type,
    training_data = remainder,
    n_nearest_neighbours = n_neighbours,
    n_years_prior = n_years_prior
  )

  holdout_predicted <- emplog_prop(holdout_counts$total_died,
                                   holdout_counts$total_tested)
  rho_implied <- fit_rho_given_predictions(holdout_data$died,
                                           holdout_data$mosquito_number,
                                           holdout_predicted)

  list(p_draws = p_draws,
       rho_implied = rho_implied,
       test_df = test_data)

}


# The oracle bound: the nearest neighbour null at whichever number of
# neighbours minimises its own mean squared error on the held-out records.
#
# This is deliberately given hindsight the dynamical model is not given, so that
# no choice of neighbour count can be said to have handicapped the baseline. The
# selection is one scalar over thousands of records, so the optimism it buys is
# small, and it runs in the conservative direction for any claim that the
# dynamical model beats the null.
#
# Minimised on mean squared error because that is what the comparison is
# reported on - excess MSE and variance explained - so the bound is a bound on
# the statistic actually quoted. The coverage and CRPS of the returned fold are
# therefore at the MSE-optimal k, not at their own optima.
nn_oracle_draws <- function(training_data, test_data,
                            n_years_prior = 1, n_draws = 1000,
                            k_grid = c(1, 2, 3, 5, 8, 12, 20, 30, 50, 80, 120,
                                       200),
                            holdout_size = 100, seed = 111) {

  observed <- test_data$died / test_data$mosquito_number

  mse <- vapply(k_grid, function(k) {
    counts <- predict_null_fixed_nn_counts(
      latitude = test_data$latitude,
      longitude = test_data$longitude,
      year = test_data$year_start,
      insecticide_type = test_data$insecticide_type,
      training_data = training_data,
      n_nearest_neighbours = k,
      n_years_prior = n_years_prior
    )
    mean((observed - emplog_prop(counts$total_died,
                                 counts$total_tested)) ^ 2)
  }, numeric(1))

  best <- k_grid[which.min(mse)]
  if (best == max(k_grid)) {
    warning("the oracle neighbour count is at the top of the grid (", best,
            "); widen k_grid so the minimum is interior")
  }

  fit <- nn_null_draws(training_data, test_data,
                       n_neighbours = best,
                       n_years_prior = n_years_prior,
                       n_draws = n_draws,
                       holdout_size = holdout_size,
                       seed = seed)

  fit$n_neighbours <- best
  fit$k_grid <- k_grid
  fit$k_mse <- mse
  fit

}
