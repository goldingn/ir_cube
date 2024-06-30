# figures to validate internal consistency of model predictions, along important
# gradients

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load the fitted model objects here, to set up predictions
load(file = "temporary/fitted_model.RData")

mask <- rast("data/clean/raster_mask.tif")

# get posterior predictive simulations of observations
died_sim <- betabinomial_p_rho(N = df$mosquito_number,
                               p = population_mortality_vec,
                               rho = rho_classes[df$class_id])
mortality_sim <- died_sim / df$mosquito_number

# summarise fit to data
died <- calculate(died_sim, values = draws, nsim = 1e4)[[1]][, , 1]

# create a dharma object to compute randomised quantile residuals and
# corresponding residual z scores
dharma <- DHARMa::createDHARMa(
  simulatedResponse = t(died),
  observedResponse = df$died,
  integerResponse = TRUE
)

# extract RQR, scaled to normal distribution for easier checking
df_validate <- df %>%
  mutate(
    z_resid = qnorm(dharma$scaledResiduals),
    # handle some Infs
    z_resid = pmin(z_resid, max(z_resid[is.finite(z_resid)])),
    z_resid = pmax(z_resid, min(z_resid[is.finite(z_resid)])),
  ) %>%
  # add on covariate values
  left_join(
    all_extract,
    by = c("cell_id", "year_id")
  )

# plot and fit models comparing this with covariate values

df_validate %>%
  ggplot(
    aes(
      x = longitude,
      y = latitude,
      colour = z_resid
    )
  ) +
  geom_spatraster(
    data = mask
  ) +
  geom_point(
    alpha = 0.5,
    size = 2
  ) +
  scale_fill_gradient(
    low = grey(0.9),
    high = grey(0.9),
    na.value = "transparent",
    guide = "none"
  ) +
  scale_colour_viridis_c() +
  theme_minimal() +
  theme_ir_maps() +
  ggtitle("Residuals")

space_fit <- mgcv::gam(
  z_resid ~ s(latitude, longitude, bs = "gp", k = 200),
  method = "REML",
  data = df_validate
)
summary(space_fit)

# make a raster of the GAM residual fit
mask_lores <- mask %>%
  terra::aggregate(10)
coords_pred <- mask_lores %>%
  terra::xyFromCell(cells(.)) %>%
  as_tibble() %>%
  rename(
    latitude = y,
    longitude = x
  ) %>%
  mutate(
    z_resid_pred = predict(space_fit, ., se.fit = TRUE)$fit,
    z_resid_sd = predict(space_fit, ., se.fit = TRUE)$se.fit
  )

z_resid_raster <- c(mask_lores, mask_lores)
names(z_resid_raster) <- c("mean", "sd") 
z_resid_raster$mean[cells(z_resid_raster)] <- as.matrix(coords_pred$z_resid_pred)
z_resid_raster$sd[cells(z_resid_raster)] <- as.matrix(coords_pred$z_resid_sd)

ggplot() +
  geom_spatraster(
    aes(
      fill = mean,
    ),
    data = z_resid_raster
  ) +
  scale_fill_distiller(
    type = "div",
    na.value = "transparent"
  ) +
  theme_ir_maps() +
  ggtitle(
    "Residual spatial correlation in model"
  )

# plot against covariates
cov_names <- colnames(all_extract)[-(1:2)]
df_validate %>%
  select(
    all_of(cov_names),
    z_resid
  ) %>%
  pivot_longer(
    cols = all_of(cov_names),
    names_to = "covariate_name",
    values_to = "covariate_value"
  ) %>%
  ggplot(
    aes(
      x = covariate_value,
      y = z_resid,
      group = covariate_name
    )
  ) +
  geom_point(
    alpha = 0.1
  ) +
  geom_smooth() +
  facet_wrap(~covariate_name) +
  theme_minimal()

# plot against time, in different regions


