# Run out-of-sample predictive validation experiments

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load bioassay data
ir_africa <- readRDS(file = "data/clean/all_gambiae_complex_data.RDS")


# Note: there is code in fit_model to subset this to some insecticides. Relevant
# code is copied here for now, but move that into the data preparation scripts
# and save only the model-ready version to load in here

# load the mask
mask <- rast("data/clean/raster_mask.tif")

# an Africa polygon for plotting
gadm_polys <- readRDS("data/clean/gadm_polys.RDS")
africa <- gadm_polys %>%
  # st_combine() %>%
  st_union()


baseline_year <- 1995

insecticides_keep <- c("Alpha-cypermethrin",
                       "Deltamethrin",
                       "Lambda-cyhalothrin", 
                       "Permethrin",
                       "Fenitrothion",
                       "Malathion",
                       "Pirimiphos-methyl",
                       "DDT",
                       "Bendiocarb")

df <- ir_africa %>%
  filter(
    insecticide_type %in% insecticides_keep
  ) %>%
  group_by(
    insecticide_type
  ) %>%
  # subset to the most common concentration for each insecticide
  filter(
    concentration == sample_mode(concentration)
  ) %>%
  ungroup() %>%
  filter(
    # drop any from before the baseline
    year_start >= baseline_year
  ) %>%
  mutate(
    # add on cell ids corresponding to these observations,
    cell = cellFromXY(mask,
                      as.matrix(select(., longitude, latitude)))
  ) %>%
  # drop a handful of datapoints missing covariates
  filter(
    !is.na(extract(mask, cell)[, 1])
  )

# Define the training and test folds for: spatial extrapolation (country
# dropout), spatial interpolation (multi-area dropout), and temporal forecasting
# (last years dropout)

# Find the countries with enough data for model validation (at least 20
# bioassays for each insecticide type) and make a table of counts
country_bioassay_counts <- df %>%
  group_by(
    country_name,
    insecticide_type) %>%
  summarise(
    records = n(),
    .groups = "drop"
  ) %>%
  # find country/insecticide combinations with a reasonable number of records
  filter(
    records > 20
  ) %>%
  # find countries with enough records of all insecticides
  group_by(
    country_name
  ) %>%
  mutate(
    n_insecticides = n()
  ) %>%
  filter(
    n_insecticides == 9
  ) %>%
  pivot_wider(
    names_from = country_name,
    values_from = records
  ) %>%
  select(
    -n_insecticides
  )

# view(country_bioassay_counts)

# pull out the country names
countries_to_validate <- country_bioassay_counts %>%
  select(-insecticide_type) %>%
  colnames()

countries_to_validate

# make a list of training and test datasets for out-of-sample validation

# Subset a dataset based on values in a field. If keep == TRUE (the default), return the
# subset of the dataset where field_name (a character for the column name) is
# one of the elements in field_values (a character vector), if keep = FALSE, return everything except for those records
split_data <- function(field_values, field_name, dataset, keep = TRUE) {
  
  if (keep) {
    
    subsetted <- dataset %>%
      filter_at(
        vars(starts_with(field_name)),
        any_vars(. %in% field_values)
      )
    
  } else {
    
    subsetted <- dataset %>%
      filter_at(
        vars(starts_with(field_name)),
        any_vars(!(. %in% field_values))
      )
    
  }
  
  subsetted
  
}

spatial_extrapolation <- list(
  training = lapply(
    countries_to_validate,
    split_data,
    field_name = "country_name",
    dataset = df,
    keep = FALSE
  ),
  test = lapply(
    countries_to_validate,
    split_data,
    field_name = "country_name",
    dataset = df,
    keep = TRUE
  )
)
names(spatial_extrapolation$training) <- countries_to_validate
names(spatial_extrapolation$test) <- countries_to_validate



# Spatial interpolation

# do k-means clustering to identify region centroids
df_locations <- df %>%
  select(
    longitude,
    latitude
  ) %>%
  as.matrix()

# set the RNG seed, as this is stochastic
set.seed(2025-07-07)
centroids <- df_locations %>%
  unique() %>%
  kmeans(
    centers = 50,
    iter.max = 300,
    nstart = 100
  ) %>%
  `[[`(
    "centers"
  )


plot(df_locations,
     asp = 1,
     pch = 16,
     cex = 0.5,
     col = grey(0.8))

points(centroids,
       pch = 16,
       cex = 0.5,
       col = "red")

# identify all points within some distance of these as, as being test data
dists <- fields::rdist.earth(
  x1 = centroids,
  x2 = df_locations,
  miles = FALSE
)

# for each datapoint, get the minimum distance to a centroid
min_dists <- apply(dists, 2, min)

# get the distance in km such that 5% of records are in the test set
test_threshold <- quantile(min_dists, 0.05)

# set half this to be the buffer distance, and define the outer edge of the
# buffer circle
buffer_distance <- 0.5 * test_threshold
buffer_threshold <- test_threshold + buffer_distance

df_interp <- df %>%
  mutate(
    fold = case_when(
      min_dists < test_threshold ~ "test",
      min_dists > buffer_threshold ~ "training",
      .default = "buffer"
    ),
    fold = factor(fold, levels = c("buffer", "training", "test"))
  )

# check the training and test splits contain all the insecticides
df_interp %>%
  filter(
    fold != "buffer"
  ) %>%
  group_by(
    fold,
    insecticide_type
  ) %>%
  summarise(
    records = n()
  ) %>%
  pivot_wider(
    names_from = fold,
              values_from  = records
  )
  

# plot the training and test split, with circles and coloured points.
mask_poly <- mask %>%
  as.polygons() %>%
  simplifyGeom(tolerance = 0.05)

train_test_col <- RColorBrewer::brewer.pal(3, "Set1")[1:2]

interp_plot <- df_interp %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = crs(mask)
  ) %>%
  arrange(fold) %>%
  ggplot(
    aes(
      colour = fold
    )
  ) +
  geom_spatvector(
    data = mask_poly,
    colour = "transparent",
    fill = grey(0.9)
  ) +
  # geom_spatraster(
  #   data = mask
  # ) +
  # scale_fill_gradient(
  #   low = grey(0.95),
  #   high = grey(0.95),
  #   na.value = "transparent",
  #   guide = "none"
  # ) +
  geom_sf() +
  scale_colour_manual(
    values = c(
      "training" = train_test_col[2],
      "test" = train_test_col[1],
      "buffer" = grey(0.8)
    )
  ) +
  coord_sf(
    ylim = range(df$latitude)
  ) +
  theme_minimal()

interp_plot_small <- interp_plot +
  coord_sf(
    xlim = c(0, 5),
    ylim = c(5, 10)
  )

interp_plot / interp_plot_small

spatial_interpolation <- list(
  training = split_data(
    field_values = "training",
    field_name = "fold",
    dataset = df_interp,
    keep = TRUE
  ),
  test = split_data(
    field_values = "test",
    field_name = "fold",
    dataset = df_interp,
    keep = TRUE
  )
)


# Temporal forecasting

# work out then the latest covariate layers are to get the three test years

nets_cube <- rast("data/clean/net_use_cube.tif")
irs_cube <- rast("data/clean/irs_coverage_scaled_cube.tif")
pop_cube <- rast("data/clean/pop_scaled_cube.tif")
nets_final_year <- nets_cube %>%
  names() %>%
  tail(1) %>%
  str_remove("nets_") %>%
  as.numeric()
irs_final_year <- irs_cube %>%
  names() %>%
  tail(1) %>%
  str_remove("irs_") %>%
  as.numeric()
pop_final_year <- pop_cube %>%
  names() %>%
  tail(1) %>%
  str_remove("pop_") %>%
  as.numeric()
final_year <- min(nets_final_year, irs_final_year, pop_final_year)
validation_years <- final_year + -2:0

# split into a training set (before those years) and a test set
temporal_forecasting <- list(
  training = split_data(
    field_values = validation_years,
    field_name = "year_start",
    dataset = df,
    keep = FALSE
  ),
  test = split_data(
    field_values = validation_years,
    field_name = "year_start",
    dataset = df,
    keep = TRUE
  )
)

# these are the train and test sets
spatial_extrapolation
spatial_interpolation
temporal_forecasting


# Define the null model: For each point in the training data, average over the X
# nearest datapoints from the current and previous year

rmse <- function(observed, predicted) {
  sqrt(mean(observed - predicted) ^ 2)
}

mae <- function(observed, predicted) {
  mean(abs(observed - predicted))
}

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
  
  training_coords <- training_data %>%
    select(longitude, latitude) %>%
    as.matrix()
  
  test_coords <- bind_cols(longitude = longitude,
                           latitude = latitude) %>%
    as.matrix()
  
  # compute the distance matrix between the test and training sets
  dists <- fields::rdist.earth(test_coords,
                               training_coords,
                               miles = FALSE)
  
  year_diff <- -1 * seq(0, n_years_prior)
  
  n_test <- nrow(test_coords)
  preds <- rep(NA, n_test)
  # loop through each test record
  for (i in seq_len(n_test)) {
    # obtain the vector of distances to that record from all training records
    distance_vec <- dists[i, ]
    
    # mask these (with Inf) if they are not in the same or the previous year or
    # not for the same insecticide type
    this_year <- year[i]
    this_insecticide_type <- insecticide_type[i]
    valid_years <- this_year + year_diff
    training_year_valid <- training$year_start %in% valid_years
    insecticide_type_valid <- training$insecticide_type == this_insecticide_type
    valid <- training_year_valid & insecticide_type_valid
    masked_distance_vec <- ifelse(valid, distance_vec, Inf)
    
    # identify the X closest records
    threshold_distance <- sort(masked_distance_vec,
                               decreasing = FALSE)[n_nearest_neighbours]
    nearest_training_idx <- which(masked_distance_vec <= threshold_distance)
    
    # compute a weighted mean of these as the prediction
    total_died <- sum(training_data$died[nearest_training_idx])
    total_tested <- sum(training_data$mosquito_number[nearest_training_idx])
    preds[i] <- total_died / total_tested
    
  }
  
  # return the predictions
  preds
  
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
      observed = died / mosquito_number
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
                                        training = spatial_interpolation$training,
                                        n_nearest_neighbours = n_neighbours[1])
    ) %>%
    ungroup()
  
  # compute rmse for this grid search on numbers of neighbours to find the optimum
  rmses <- grid_search %>%
    group_by(
      n_neighbours
    ) %>%
    summarise(
      rmse = rmse(observed, predicted)
    )
  
  # maybe plot the relationship to check for convexity
  if (plot) {
    plot(rmse ~ n_neighbours,
         data = rmses,
         type = "b")
  }
  
  # pull out the optimal value
  optimal_nn <- rmses %>%
    filter(rmse == min(rmse)) %>%
    slice(1) %>%
    pull(n_neighbours)
  
  # add a flag for optimality and return
  grid_search %>%
    mutate(
      nn_is_optimal = n_neighbours == optimal_nn
    )
  
}


# Do grid search to find the optimal number of nearest neighbours for spatial
# interpolation and temporal forecasting

spatial_interpolation_preds <- spatial_interpolation$test %>%
  mutate(
    mortality = died / mosquito_number
  ) %>%
  predict_null_optimal_nn(
    training_data = spatial_interpolation$training
  )

# for temporal forecasting, use up to 3 years prior to enable prediction to the
# third year into the future
temporal_forecasting_preds <- temporal_forecasting$test %>%
  mutate(
    mortality = died / mosquito_number
  ) %>%
  predict_null_optimal_nn(
    training_data = temporal_forecasting$training,
    n_years_prior = 3
  )

# for spatial extrapolation, need to run it for each country and compute the average rmses to
# RMSEs to identify the optimal number of neighbours
each_country_optimal_nn <- function(test, training) {
  test %>%
    mutate(
      mortality = died / mosquito_number
    ) %>%
    predict_null_optimal_nn(
      training_data = training
    )
}

spatial_extrapolation_preds <- mapply(
  FUN = each_country_optimal_nn,
  spatial_extrapolation$test,
  spatial_extrapolation$training,
  SIMPLIFY = FALSE
) %>%
  do.call(
    bind_rows, .
  )

# get and plot the overall rmses
spatial_extrapolation_rmses <- spatial_extrapolation_preds %>%
  group_by(
    n_neighbours
  ) %>%
  summarise(
    rmse = rmse(observed, predicted)
  )

plot(rmse ~ n_neighbours,
     data = spatial_extrapolation_rmses,
     type = "b")

# find the optimum
spatial_extrapolation_optimal_nn <- spatial_extrapolation_rmses %>%
  filter(rmse == min(rmse)) %>%
  slice(1) %>%
  pull(n_neighbours)

# overwrite the optimal flag
spatial_extrapolation_preds <- spatial_extrapolation_preds %>%
  mutate(
    nn_is_optimal = n_neighbours == spatial_extrapolation_optimal_nn
  )

# subset each of these validation experiments to the optimal numbers of
# neighbours
optimal_nn_preds <- bind_rows(
  spatial_interpolation = spatial_interpolation_preds,
  spatial_extrapolation = spatial_extrapolation_preds,
  temporal_forecasting = temporal_forecasting_preds,
  .id = "experiment"
) %>%
  filter(
    nn_is_optimal
  )

# record the optima
optimal_nn <- optimal_nn_preds %>%
  group_by(
    experiment
  ) %>%
  summarise(
    n_neighbours = n_neighbours[1]
  )
optimal_nn

# compute bias and RMSE for each of these and tabulate
optimal_nn_preds %>%
  group_by(
    experiment
  ) %>%
  summarise(
    rmse = rmse(observed, predicted),
    mae = mae(observed, predicted),
    bias = mean(predicted - observed),
    .groups = "drop"
  )

# split these by country (spatial extrapolation) and years ahead (temporal
# forecasting)
optimal_nn_preds %>%
  filter(
    experiment == "spatial_extrapolation"
  ) %>%
  group_by(
    country_name
  ) %>%
  summarise(
    rmse = rmse(observed, predicted),
    mae = mae(observed, predicted),
    bias = mean(predicted - observed),
    .groups = "drop"
  )

# The biases are positive for some countries and negative for others. They
# cancel out in aggregate, but this is suboptimal


optimal_nn_preds %>%
  filter(
    experiment == "temporal_forecasting"
  ) %>%
  group_by(
    year_start
  ) %>%
  summarise(
    rmse = rmse(observed, predicted),
    mae = mae(observed, predicted),
    bias = mean(predicted - observed),
    .groups = "drop"
  )

# these are substantially more biased further into the future, as expected


# should these validation studies be limited to dropouts in recent years, to
# capture heterogeneity and average levels of resistance?

# should the error and bias measures attempt to account for the binomial
# variance (errors will be bigger closer to 0.5, so scale according to this?)

# Can we use Penny's previously estimated maps (to 2017) in this comparison?
# https://doi.org/10.6084/m9.figshare.9912623 