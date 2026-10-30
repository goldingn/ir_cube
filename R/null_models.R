# Null model definitions used in the cross-validation experiments: the
# insecticide-type intercept model and the nearest neighbour heuristic.
#
# Separated from the scoring so that both the point-prediction
# scoring and the posterior predictive scoring (#10) use one definition of each
# null model. The lower half of this file adds predictive distributions for the
# nulls, because a proper scoring rule needs a distribution rather than a point
# from every candidate. Every model is scored at the same external,
# replicate-based overdispersion (see validation_metrics.R).
#
# The nulls do not fit an overdispersion of their own. They earn their place on
# point prediction - mean squared error and variance explained - and that needs
# only the predicted fraction, while the dynamical model's calibration is judged
# against held-out data directly. Fitting one rho per null made coverage a
# comparison of vagueness rather than of prediction: the intercept null reached
# nominal coverage by inflating its rho to 0.45 against an external estimate of
# 0.12-0.22 (#12 review).

# the nearest neighbour null, and predictive distributions for the nulls -----

# The year the nearest neighbour null's window ends at, for a record of `year`.
# The window is anchored at prediction time, not at the record's own year: the
# most recent `n_years_prior` + 1 years of data that exist when the prediction
# is made. For the spatial experiments training spans every year, so the anchor
# is the record's own year and the window is the familiar {year, year - 1}.
# For a forecast the training data stop at the cut, so the anchor is the last
# training year and every held-out record draws on the same window — which is
# the situation a person forecasting from that cut is actually in.
#
# Anchoring on the record's own year instead made the null weaker the further
# ahead it had to predict, for an indexing reason rather than an information
# one: a record one year past the cut saw the whole window, one five years past
# saw a single year of it, or none at all.
nn_anchor_year <- function(year, training_years) {
  pmin(year, max(training_years))
}

# For each held-out record, average over the k nearest training records of the
# same insecticide from the most recent two years available.

# Pooled counts over the nearest neighbours of each held-out record - died and
# tested summed over them, rather than a single proportion, so that uncertainty
# in the predicted fraction can be represented by a beta posterior on those
# counts.
#
# Several neighbour counts at once. Neither the distance matrix nor the
# per-record ordering of those distances depends on k, so both are computed once
# and every k is read off the same cumulative sums.
#
# `k_values` is a vector; the returned matrices have one row per held-out record
# and one column per element of it, in that order.
nn_counts_grid <- function(latitude,
                           longitude,
                           year,
                           insecticide_type,
                           training_data,
                           k_values,
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
    is.numeric(k_values),
    length(k_values) >= 1,
    all(is.finite(k_values)),
    all(k_values >= 1)
  )

  training_coords <- as.matrix(training_data[, c("longitude", "latitude")])
  test_coords <- cbind(longitude, latitude)

  dists <- fields::rdist.earth(test_coords, training_coords,
                               miles = FALSE)

  # the year window, nn_anchor_year()
  anchor <- nn_anchor_year(year, training_data$year)
  year_diff <- -1 * seq(0, n_years_prior)

  n_test <- nrow(test_coords)
  total_died <- matrix(NA_real_, nrow = n_test, ncol = length(k_values))
  total_tested <- matrix(NA_real_, nrow = n_test, ncol = length(k_values))

  # which training records a record may draw on depends only on its anchor year
  # and its insecticide, so the masks are built once per distinct pair rather
  # than once per record
  mask_key <- paste(anchor, insecticide_type)
  unique_key <- unique(mask_key)
  key_index <- match(mask_key, unique_key)
  valid_masks <- lapply(match(unique_key, mask_key), function(j) {
    training_data$year %in% (anchor[j] + year_diff) &
      training_data$insecticide_type == insecticide_type[j]
  })

  training_died <- training_data$died
  training_tested <- training_data$mosquito_number

  for (i in seq_len(n_test)) {

    distance_vec <- dists[i, ]
    valid <- valid_masks[[key_index[i]]]
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
      valid_years <- anchor[i] + year_diff
      stop("no training records in ", paste(range(valid_years), collapse = "-"),
           " for ", insecticide_type[i], ": the nearest neighbour null has ",
           "nothing to predict from")
    }

    # the k-th smallest masked distance is the threshold, and every record at or
    # inside it is pooled - so ties at the threshold all count, and a k larger
    # than the number of valid records gives an infinite threshold and pools the
    # whole training set, exactly as the per-k version did. Sorting once and
    # taking cumulative sums along that order makes each k an index lookup
    order_index <- order(masked_distance_vec)
    sorted_distance <- masked_distance_vec[order_index]
    cumulative_died <- cumsum(training_died[order_index])
    cumulative_tested <- cumsum(training_tested[order_index])

    # findInterval() on an ascending vector returns the number of elements at or
    # below the threshold, which is the size of the pooled set
    n_pooled <- findInterval(sorted_distance[k_values], sorted_distance)

    total_died[i, ] <- cumulative_died[n_pooled]
    total_tested[i, ] <- cumulative_tested[n_pooled]

  }

  list(total_died = total_died,
       total_tested = total_tested)

}

# the same for a single neighbour count, as a data frame of one column each
predict_null_fixed_nn_counts <- function(latitude,
                                         longitude,
                                         year,
                                         insecticide_type,
                                         training_data,
                                         n_nearest_neighbours,
                                         n_years_prior = 1) {

  stopifnot(length(n_nearest_neighbours) == 1)

  counts <- nn_counts_grid(
    latitude = latitude,
    longitude = longitude,
    year = year,
    insecticide_type = insecticide_type,
    training_data = training_data,
    k_values = n_nearest_neighbours,
    n_years_prior = n_years_prior
  )

  data.frame(total_died = counts$total_died[, 1],
             total_tested = counts$total_tested[, 1])

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
# type
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

  list(p_draws = p_draws,
       test_df = test_data)

}

# Predictive distribution of the nearest neighbour null.
#
# This is not a competitor model to be optimised, it is a stand-in for what a
# person would do without a model: read a site's value off the nearby recent
# surveys. Tuning it separately for each experiment would build a series of new
# models to compare against the one model under analysis, so it is not tuned at
# all: it is the practice baseline, `n_neighbours = 1`, `n_years_prior = 1`, the
# single nearest record in the most recent two years available at prediction
# time. Short recency is justified a priori, since anyone running this
# surveillance knows resistance is moving fast, and the k-curves agree.
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
                          n_years_prior = 1, n_draws = 1000) {

  counts <- predict_null_fixed_nn_counts(
    latitude = test_data$latitude,
    longitude = test_data$longitude,
    year = test_data$year,
    insecticide_type = test_data$insecticide_type,
    training_data = training_data,
    n_nearest_neighbours = n_neighbours,
    n_years_prior = n_years_prior
  )

  p_draws <- pooled_count_draws(counts$total_died,
                                counts$total_tested,
                                n_draws)

  list(p_draws = p_draws,
       test_df = test_data)

}
