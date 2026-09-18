# How predictive skill depends on how far the held-out record sits from the
# training data.
#
# Each cross-validation experiment reduces to a single number, but the folds
# differ enormously in difficulty: a held-out record two cells from a training
# record is a different problem from one 800 km away. Reporting skill against
# distance and against local data volume turns each fold set into a continuum of
# difficulty rather than one summary, makes the arbitrary geometry of the folds
# matter less, and answers directly the question a user of the map has — at what
# separation, and at what data density, should the mechanistic model be trusted
# over local interpolation (#12 review).
#
# Three axes, for every held-out record:
#
#   distance to the nearest training record of the same insecticide, in km
#   number of training records of the same insecticide within `radius_km`
#   years since the most recent training observation at that same pixel
#
# Distances are between pixel centroids, because the pixel is the unit the model
# predicts at. Sourcing this file needs only the fold definitions and the saved
# scores; nothing is refitted.

source("R/validation_functions.R")
source("R/validation_folds.R")

suppressMessages({
  library(dplyr)
  library(tidyr)
})

radius_km <- 100

# breaks chosen so that no bin is so thin that its excess MSE is noise: the
# distance breaks straddle the ~13 km minimum separation the interpolation fold
# guarantees, and the ~5 km grid resolution
distance_breaks <- c(0, 10, 25, 50, 100, 250, 500, Inf)
volume_breaks <- c(0, 1, 5, 20, 100, Inf)
lag_breaks <- c(-Inf, 0, 1, 2, 3, 5, 10, Inf)

cell_coordinates <- function(cells) {
  xy <- terra::xyFromCell(mask, cells)
  colnames(xy) <- c("longitude", "latitude")
  xy
}

# the three axes for one fold. Computed on distinct (cell, year, insecticide)
# keys rather than per assay, since that is all they depend on
fold_geometry <- function(training, test, experiment, fold) {

  keys <- test %>%
    distinct(cell, year_start, insecticide_type)

  distance <- fields::rdist.earth(cell_coordinates(keys$cell),
                                  cell_coordinates(training$cell),
                                  miles = FALSE)

  same_insecticide <- outer(keys$insecticide_type,
                            training$insecticide_type,
                            FUN = "==")
  distance_same <- distance
  distance_same[!same_insecticide] <- Inf

  # the most recent training year at the same pixel. Negative where the nearest
  # training observation at that pixel is from a later year than the held-out
  # record, which should not happen in a forecasting experiment and is how the
  # training set leak showed itself
  last_year <- tapply(training$year_start, training$cell, max)

  keys %>%
    mutate(experiment = experiment,
           fold = fold,
           distance_any = apply(distance, 1, min),
           distance_same = apply(distance_same, 1, min),
           n_within_radius = rowSums(distance_same <= radius_km),
           years_since_cell = year_start -
             as.numeric(last_year[as.character(cell)]),
           .before = everything())

}

geometry <- bind_rows(
  bind_rows(
    lapply(seq_along(countries_to_validate), function(index) {
      fold_geometry(spatial_extrapolation$training[[index]],
                    spatial_extrapolation$test[[index]],
                    "spatial_extrapolation",
                    countries_to_validate[index])
    })
  ),
  fold_geometry(spatial_interpolation$training,
                spatial_interpolation$test,
                "spatial_interpolation",
                "all"),
  fold_geometry(temporal_forecasting$training,
                temporal_forecasting$test,
                "temporal_forecasting",
                "all")
)

write.csv(geometry, "outputs/cv_geometry.csv", row.names = FALSE)

cat("\ndistance from each held-out pixel to the nearest training record",
    "of the same insecticide (km):\n")
print(geometry %>%
        group_by(experiment) %>%
        summarise(keys = n(),
                  min = min(distance_same),
                  q25 = quantile(distance_same, 0.25),
                  median = median(distance_same),
                  q75 = quantile(distance_same, 0.75),
                  max = max(distance_same),
                  .groups = "drop") %>%
        mutate(across(where(is.numeric), ~ round(.x, 1))) %>%
        as.data.frame())

cat("\nyears since the last training observation at the same pixel",
    "(negative means a later year, which is a training set leak):\n")
print(geometry %>%
        group_by(experiment) %>%
        summarise(keys = n(),
                  pixel_never_sampled = sum(is.na(years_since_cell)),
                  later_year_in_training = sum(years_since_cell < 0,
                                               na.rm = TRUE),
                  median = median(years_since_cell, na.rm = TRUE),
                  max = suppressWarnings(max(years_since_cell, na.rm = TRUE)),
                  .groups = "drop") %>%
        as.data.frame())


# skill against each axis --------------------------------------------------

scores <- read.csv("outputs/cv_scores.csv", encoding = "UTF-8") %>%
  left_join(geometry,
            by = c("experiment", "fold", "cell", "year_start",
                   "insecticide_type"))
stopifnot(!anyNA(scores$distance_same))

# excess mean squared error, and the share of the intercept null's excess that
# each model removes, within bins of one axis
by_axis <- function(scores, axis, breaks, label) {
  scores %>%
    filter(!is.na(.data[[axis]])) %>%
    mutate(bin = cut(.data[[axis]], breaks = breaks, include.lowest = TRUE)) %>%
    group_by(experiment, bin, model) %>%
    summarise(n = n(),
              mse = mean((observed - predicted) ^ 2),
              mse_floor = noise_floor_mse(died, mosquito_number, rho_external),
              .groups = "drop") %>%
    group_by(experiment, bin) %>%
    mutate(axis = label,
           excess = mse - mse_floor,
           explained = 1 - excess / excess[model == "intercept"]) %>%
    ungroup() %>%
    select(axis, experiment, bin, model, n, mse, mse_floor, excess, explained)
}

skill_by_geometry <- bind_rows(
  by_axis(scores, "distance_same", distance_breaks,
          "km to nearest same-insecticide training record"),
  by_axis(scores, "n_within_radius", volume_breaks,
          sprintf("same-insecticide training records within %i km", radius_km)),
  by_axis(scores, "years_since_cell", lag_breaks,
          "years since the last observation at that pixel")
)

write.csv(skill_by_geometry, "outputs/cv_skill_by_geometry.csv",
          row.names = FALSE)

for (this_axis in unique(skill_by_geometry$axis)) {
  cat("\nexcess MSE by", this_axis, ":\n")
  print(skill_by_geometry %>%
          filter(axis == this_axis) %>%
          select(experiment, bin, model, n, excess, explained) %>%
          pivot_wider(names_from = model,
                      values_from = c(n, excess, explained)) %>%
          mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
          as.data.frame())
}
