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

