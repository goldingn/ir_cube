# summarise model

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load the fitted model objects here, to set up predictions
load(file = "temporary/fitted_model.RData")

# load covariate rasters (need to reload these even if restoring workspace
# because pointers)
mask <- rast("data/clean/raster_mask.tif")
nets_cube <- rast("data/clean/nets_per_capita_cube.tif")
covs_flat <- rast("data/clean/flat_covariates.tif")

# summarise the covariate effect sizes (at insecticide class level)
effect_sizes <- summary(calculate(exp(beta_class[, 1]), values = draws))$statistics[, c("Mean", "SD")]
rownames(effect_sizes) <- c("itn_coverage", names(covs_flat))
round(effect_sizes, 2)


died_sim <- betabinomial_p_rho(N = df$mosquito_number,
                               p = population_mortality_vec,
                               rho = rho)
mortality_sim <- died_sim / df$mosquito_number

# summarise fit to data
died <- calculate(died_sim, values = draws, nsim = 1e3)[[1]][, , 1]

# create a dharma object to compute randomised quantile residuals and
# corresponding residual z scores
dharma <- DHARMa::createDHARMa(
  simulatedResponse = t(died),
  observedResponse = df$died,
  integerResponse = TRUE
)

# looks fairly uniform
length(types)
par(mfrow = n2mfrow(length(types)))
for (type in types) {
  hist(dharma$scaledResiduals[df$insecticide_type == type],
       main = type)
}

# DDT the most obviously skewed
not_ddt_index <- df$insecticide_type != "DDT"
par(mfrow = c(2, 1))
hist(dharma$scaledResiduals,
     breaks = 100,
     main = "all")
hist(dharma$scaledResiduals[not_ddt_index],
     breaks = 100,
     main = "not DDT")

# is there excessive dispersion in DDT and underdispersion in Deltamethrin?
# make rho vary by insecticide class with a hierarchy in logit_rho?

dharma_not_ddt <- DHARMa::createDHARMa(
  simulatedResponse = t(died[, not_ddt_index]),
  observedResponse = df$died[not_ddt_index],
  integerResponse = TRUE
)

par(mfrow = c(1, 1))
plot(dharma_not_ddt)

# no evidence of temporal variation in model misspecification
z <- qnorm(dharma$scaledResiduals)
plot(z ~ jitter(df$year_start),
     cex = 0.5)
abline(h = 0)

plot(dharma)

# evidence of slight misspecification
hist(dharma$scaledResiduals)
plot(z ~ df$)

