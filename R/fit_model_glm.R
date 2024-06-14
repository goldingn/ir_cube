# A quick and dirty non-Bayesian model generalised linear mixed effects model of
# IR spread over time.

# The evolution of a single gene over time, with temporally-constant selection
# pressure, can be modelled (with only very slight approximation) as a logistic
# curve on time. This enables us to fit the temporally-static selection pressure
# model via a logistic regression, and capture correlations between insecticide
# types with mixed effect modelling, all of which can be carried out using
# standard, quick, non-Bayesian inference software.

# Note that this trick wil no longer work when the assumption of temporal
# variation in selection pressure is relaxed, since we will need to solve
# through time to compute thcumulative selection pressures and fractions with
# each geontype. We will also need to move to a flexible Bayesian framework to
# capture data from assays with differing concentrations, an incorporating
# imperefect susceptibility and resistance priors.

# In the single-insecticide case, we can use this model to for the evolution of
# resistance:
#   p_mort = plogis(eta)
#   eta = int + beta * time
#   beta = a + b * x

# where int is the logit mortality rate when time = 0 (set to year 2000 for
# interpretability), beta is the log of the relative fitness of susceptible
# mosquitoes over resistant mosquitoes (which we expect to be negative in this
# context), and this in turn is given by a linear model, with intercept a
# (fitness in the absence of covariate x) and slope b (relationship between
# fitness and covariate x)

# we can express the same model in the following way:
#   p_mort = plogis(eta)
#   eta = int + (a + b * x) * time
#       = int + a * time + b * x * time

# and by defining the fixed variables (creating the covariates) x1 = time, and
# x2 = x * time, we can fit the model as a three-parameter logistic regression:
#   p_mort = plogis(eta)
#   eta = alpha + beta1 * x1 + beta2 * x2
# where, if we map back to the original formulation: alpha = int, beta1 = a,
# beta2 = b

# we also wish for all of these parameters to differ between insecticides, so we
# use a random intercepts and random slopes model. In glmer notation, that is
# given as:

# glmer(cbind(died, survived) ~ (1| insecticide_type) +
#         (x1 | insecticide_type) +
#         (x2 | insecticide_type),
#       family = binomial,
#       data = df)

# we should also be able to structure the hierarchical relationship as having a
# level for the class of insecticides, simultaneously consider multiple species,
# and, using INLA or similar, have a betabinomial observation model.



# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load and prep bioassay data
ir_mtm_africa <- readRDS(file = "data/clean/mtm_data.RDS")

# set the start of the timeseries considered
baseline_year <- 1995

df <- ir_mtm_africa %>%
  filter(
    # subset to An. gambiae (s.l./s.s.)
    species %in% c("An. gambiae s.l.", "An. gambiae s.s."),
    # subset to WHO tube tests (most of data)
    test_type == "WHO tube test",
    # drop the minor classes: pyrroles and neonicotinoids, and Dieldrin
    insecticide_class %in% c("Pyrethroids",
                             "Carbamates",
                             "Organochlorines",
                             "Organophosphates"),
    # drop Dieldrin as it has two concentrations in the dataset, and those with
    # fewer than 200 observations (manually as I'm being lazy)
    !(insecticide_type %in% c("Dieldrin",
                              "Carbosulfan",
                              "Cyfluthrin",
                              "Etofenprox",
                              "Propoxur")),
    # drop any from before when we have data on net coverage
    year_start >= baseline_year
  ) %>%
  mutate(
    # convert concentrations into numeric values
    concentration = as.numeric(str_remove(insecticide_conc, "%")),
    # add in time variable, in years since start of timeseries
    time = year_start - baseline_year
  ) # %>%
  # # subset to Burkina Faso for debugging
  # filter(
  #   country_name == "Burkina Faso"
  # )

# map the data to the insecticide classes and types
classes <- unique(df$insecticide_class)
types <- unique(df$insecticide_type)

# load and temporally average raster data
nets_cube <- rast("data/clean/nets_per_capita_cube.tif")
nets_flat <- terra::app(nets_cube, mean)
names(nets_flat) <- "itn_percap"

# load the other layers
covs_flat_orig <- rast("data/clean/flat_covariates.tif")

# add on the flat nets
covs_flat <- c(nets_flat, covs_flat_orig)


# extract at coordinates
coords <- df %>%
  select(longitude, latitude) %>%
  as.matrix()

# # extract from the cube
# nets_lookup <- terra::extract(nets_cube, coords) %>%
#   mutate(
#     row = row_number(),
#     .before = everything()
#   ) %>%
#   pivot_longer(cols = starts_with("nets_"),
#                names_prefix = "nets_",
#                names_to = "year",
#                values_to = "net_coverage") %>%
#   mutate(
#     year = as.numeric(year)
#   )


# create response data and covariates for modelling
df_glmer <- df %>%
  # add on spatial covariates
  bind_cols(
    terra::extract(covs_flat, coords)
  ) %>%
  # # add on spatiotemporal covariates
  # mutate(
  #   row = row_number()
  # ) %>%
  # left_join(nets_lookup,
  #           by = c("row", year_start = "year")) %>%
  # add on the success/failure columns for R's binomial interface
  mutate(
    survived = mosquito_number - died
  ) %>%
  # add an intercept for the fitness sub model
  mutate(
    intercept = 1
  ) %>%
  # interact covariates with time to get fitness submodel covariates factors
  mutate(
    across(
      c(intercept,
        itn_percap,
        crop_pc_1,
        crop_pc_2, 
        crop_pc_3, 
        crop_pc_4, 
        crop_pc_5),
      ~time * .,
      .names = "fitness_{.col}"
    )
  ) %>%
  # drop some observations outside the covariate rasters
  filter(!is.na(itn_percap))

# fit a glm to this
m <- glmer(cbind(died, survived) ~
             (1 + #fitness_intercept +
                fitness_itn_percap | insecticide_type),
                # fitness_crop_pc_1 +
                # fitness_crop_pc_2 +
                # fitness_crop_pc_3 +
                # fitness_crop_pc_4 +
                # fitness_crop_pc_5 + 
           family = stats::binomial,
           data = df_glmer)

print(m)
summary(m)
ranef(m)

# plot fitted vs observed (but noisy AF) bioassay data
plot((I(100 - df_glmer$mortality_adjusted)) ~ I(100 * (1 - fitted(m))),
     xlab = "predicted resistance",
     ylab = "bioassay resistance",
     cex = 0.5,
     lwd = 0.2)

# plot predicted and observed mortality rates over time
df_plot <- expand_grid(
  year = 1990:2024,
  insecticide_type = types,
  itn_percap = c(0.05, 0.5, 0.95)
) %>%
  mutate(
    intercept = 1,
    time = year - baseline_year,
    crop_pc_1 = 0,
    crop_pc_2 = 0,
    crop_pc_3 = 0,
    crop_pc_4 = 0,
    crop_pc_5 = 0
  ) %>%
  mutate(
    across(
      c(intercept,
        itn_percap,
        crop_pc_1,
        crop_pc_2, 
        crop_pc_3, 
        crop_pc_4, 
        crop_pc_5),
      ~time * .,
      .names = "fitness_{.col}"
    )
  ) %>%
  mutate(
    mortality = predict(m,
                   newdata = .,
                   type = "response", ),
    `mortality (%)` = mortality * 100
  )


# plot these
df_plot %>%
  ggplot(
    aes(
      x = year,
      y = `mortality (%)`,
      group = itn_percap,
      colour = itn_percap
    )
  ) +
  geom_line() +
  facet_wrap(~insecticide_type) +
  geom_point(
    aes(
      x = year_start,
      y = mortality_adjusted,
      colour = itn_percap
    ),
    data = df_glmer,
    alpha = 0.1
  ) +
  theme_minimal()

# now predict back to rasters

# get all cells in mastergrids, and extract the values there
cells <- terra::cells(covs_flat)
# nets_cube_dataframe <- terra::extract(nets_cube, cells) %>%
#   mutate(
#     cell = cells,
#     .before = everything()
#   ) %>%
#   pivot_longer(cols = starts_with("nets_"),
#                names_prefix = "nets_",
#                names_to = "year",
#                values_to = "net_coverage") %>%
#   mutate(
#     year = as.numeric(year)
#   )

covs_flat_dataframe <- terra::extract(covs_flat, cells) %>%
  mutate(
    cell = cells,
    .before = everything()
  )

# loop through insecticides and years, making predictions
years <- 2000:2030
template_raster <- nets_cube$nets_2000 * 0
types_plot <- c("Deltamethrin",
                "Permethrin",
                "Alpha-cypermethrin")

mask <- covs_flat$itn_percap > 0
mask[mask == 0] <- NA
for (this_insecticide in types_plot) {
  for (this_year in years) {
    
    # extract the nets data for this year and prepare for prediction
    df_predict <- covs_flat_dataframe %>%
      mutate(
        insecticide_type = this_insecticide,
        year = this_year,
        intercept = 1,
        time = year - baseline_year
      ) %>%
      mutate(
        across(
          c(intercept,
            itn_percap,
            crop_pc_1,
            crop_pc_2, 
            crop_pc_3, 
            crop_pc_4, 
            crop_pc_5),
          ~time * .,
          .names = "fitness_{.col}"
        )
      ) %>%
      mutate(
        mortality = predict(m,
                            newdata = .,
                            type = "response")
      )
    
    # stick this in a raster
    this_ir_raster <- template_raster
    this_ir_raster[cells] <- df_predict$mortality
    this_ir_raster <- this_ir_raster * mask
    
    write_path <- file.path("outputs/ir_maps",
                            this_insecticide,
                            sprintf("ir_%s_susceptibility.tif",
                                    this_year))

    # make sure the directory exists
    dir.create(dirname(write_path))
    writeRaster(this_ir_raster,
                write_path)
    
  }
}

# ggplot() +
#   geom_spatraster(data = this_ir_raster) +
#   scale_fill_gradient(limits = c(0, 1),
#                       na.value = "transparent") +
#   theme_minimal()

# now redo this with more spatial covariates
