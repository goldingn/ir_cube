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

# given tibbles of test and training data, return the test data tibble augmented
# with observed and predicted values of the susceptibility fraction from a
# weighted average of the nearest points in the training data of that
# insecticide type, and in the same year or up to `n_years_prior` earlier years,
# for multiple values of the number of neighbours, across a grid search between
# the values in `n_nearest_neighbour_range`, evaluating integers
# `n_nearest_neighbour_delta` apart. If plot = TRUE, plot RMSE against the
# numbers of neighbours for visual assessment of convexity. The number of
# neighbours yielding the minimal RMSE is identified by the column
# `nn_is_optimal`; filtering on this column will yield the optimal predictions
predict_null_optimal_nn <- function(test_data,
                                    training_data,
                                    n_nearest_neighbour_range = c(1, 20),
                                    n_nearest_neighbour_delta = 1,
                                    n_years_prior = 1,
                                    plot = TRUE) {
  
  # Do grid search to find the optimal number of nearest neighbours for spatial
  # interpolation
  
  # nearest neighbour values to try
  nn_values <- seq(from = n_nearest_neighbour_range[1],
                   to = n_nearest_neighbour_range[2],
                   by = n_nearest_neighbour_delta)
  
  # test data
  grid_search <- test_data %>%
    mutate(
      observed = prop(died, mosquito_number)
    ) %>%
    # add on nearest neighbour values to try (enforce integers)
    expand_grid(
      n_neighbours = round(nn_values)
    ) %>%
    # batch predictions by each value of the neighbour parameter
    group_by(
      n_neighbours
    ) %>%
    # compute predictions and observed fractions
    mutate(
      predicted = predict_null_fixed_nn(longitude = longitude,
                                        latitude = latitude,
                                        year = year_start,
                                        insecticide_type = insecticide_type,
                                        training = training_data,
                                        n_nearest_neighbours = n_neighbours)
    ) %>%
    ungroup()
  
  # compute rmse for this grid search on numbers of neighbours to find the optimum
  pred_errors <- grid_search %>%
    group_by(
      n_neighbours
    ) %>%
    summarise(
      pred_error = betabinom_dev(died = died,
                                 mosquito_number = mosquito_number,
                                 predicted = predicted)
    )
  
  # maybe plot the relationship to check for convexity
  if (plot) {
    plot(pred_error ~ n_neighbours,
         data = pred_errors,
         type = "b")
  }
  
  # pull out the optimal value
  optimal_nn <- pred_errors %>%
    filter(pred_error == min(pred_error)) %>%
    slice(1) %>%
    pull(n_neighbours)
  
  # add a flag for optimality and return
  grid_search %>%
    mutate(
      nn_is_optimal = n_neighbours == optimal_nn
    )
  
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

  training_coords <- as.matrix(training_data[, c("longitude", "latitude")])
  test_coords <- cbind(longitude, latitude)

  dists <- fields::rdist.earth(test_coords, training_coords,
                               miles = FALSE)
  year_diff <- -1 * seq(0, n_years_prior)

  n_test <- nrow(test_coords)
  total_died <- rep(NA_real_, n_test)
  total_tested <- rep(NA_real_, n_test)

  for (i in seq_len(n_test)) {

    distance_vec <- dists[i, ]
    valid_years <- year[i] + year_diff
    valid <- training_data$year_start %in% valid_years &
      training_data$insecticide_type == insecticide_type[i]
    masked_distance_vec <- ifelse(valid, distance_vec, Inf)

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
nn_null_draws <- function(training_data, test_data, n_neighbours,
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
