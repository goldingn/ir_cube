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

# load data
ir_mtm_africa <- readRDS(file = "data/clean/mtm_data.RDS")

df <- ir_mtm_africa %>%
  filter(
    # subset to An. gambiae (s.l./s.s.)
    species %in% c("An. gambiae s.l.", "An. gambiae s.s."),
    # subset to WHO tube tests (most of data)
    test_type == "WHO tube test",
    # drop DDT (the only organochlorine tested) and the minor classes: pyrroles and
    # neonicotinoids
    insecticide_class %in% c("Pyrethroids",
                             "Carbamates",
                             "Organochlorines",
                             "Organophosphates")
  ) %>%
  # convert concentrations into numeric values
  mutate(
    concentration = as.numeric(str_remove(insecticide_conc, "%"))
  ) %>%
  # subset to Burkina Faso for debugging
  filter(
    country_name == "Burkina Faso"
  )

# fit a model to these data, intercept-only to start

# map the data to the insecticide classes and types
classes <- unique(df$insecticide_class)
types <- unique(df$insecticide_type)


# add in time variable, in years since 2000
baseline_year <- 2000
df$time <- df$year_start - baseline_year

# create response data and covariates for modelling
df_glmer <- df %>%
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
      c(intercept),
      ~time * .,
      .names = "fitness_{.col}"
    )
  )

# fit a glm to this
m <- glmer(cbind(died, survived) ~ (1| insecticide_type) +
        (fitness_intercept | insecticide_type),
      family = stats::binomial,
      data = df_glmer)

# plot predicted and observed mortality rates
df_plot <- expand_grid(
  year = 1990:2024,
  insecticide_type = types
) %>%
  mutate(
    intercept = 1,
    time = year - baseline_year
  ) %>%
  mutate(
    across(
      c(intercept),
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
      group = insecticide_type
    )
  ) +
  geom_line() +
  facet_wrap(~insecticide_type) +
  geom_point(
    aes(
      x = year_start,
      y = mortality_adjusted,
    ),
    data = df,
    alpha = 0.1
  ) +
  theme_minimal()

# now redo this with a spatial covariate of ITN access
