library(tidyverse)
library(readxl)
ir_mtm <- read_xlsx(
  path = "~/Downloads/MTM_DISCRIMINATING_CONCENTRATION_BIOASSAY_20240610.xlsx",
  sheet = "Data") %>%
  mutate(
    across(c(LATITUDE, LONGITUDE, MOSQUITO_NUMBER, MORTALITY_ADJUSTED),
           as.numeric)
  ) %>%
  mutate(DIED = round(MOSQUITO_NUMBER *MORTALITY_ADJUSTED / 100),
         SURVIVED = MOSQUITO_NUMBER - DIED)

africa_countries <- c("Burkina Faso", "Botswana",
                      "Mali", "Kenya",
                      "Malawi", "Cameroon", "Niger", "Benin", "Côte d’Ivoire",
                      "Togo", "Nigeria", "United Republic of Tanzania", "Ethiopia",
                      "Mozambique", "Uganda", "Ghana",
                      "Senegal", "Democratic Republic of the Congo",
                      "Somalia", "South Africa", "Zambia", "Eritrea",
                      "Sudan", "Sierra Leone",
                      "Angola", "Namibia",
                      "Djibouti", "Rwanda",
                      "Chad", "Liberia", "Zimbabwe", "Eswatini",
                      "Guinea", "Equatorial Guinea",
                      "Congo",
                      "Central African Republic",
                      "Madagascar", "Comoros", "Sao Tome and Principe",
                      "Mauritania", "Morocco", "Cabo Verde", "Gabon", "Burundi",
                      "Guinea-Bissau", "Gambia"
)

ir_mtm_africa <- ir_mtm %>%
  filter(
    COUNTRY_NAME %in% africa_countries,
    YEAR_START >= 2010
  ) %>%
  # throw out a dodgy coordinate
  filter(LONGITUDE > -50)

# plot these points
ir_mtm_africa %>%
  arrange(YEAR_START) %>%
  ggplot(
    aes(x = LONGITUDE,
        y = LATITUDE,
        colour = YEAR_START)
  ) +
  geom_point() +
  theme_minimal() +
  coord_fixed()


ir_mtm_africa %>%
  arrange(YEAR_START) %>%
  ggplot(
    aes(x = LONGITUDE,
        y = LATITUDE,
        colour = MORTALITY_ADJUSTED)
  ) +
  geom_point() +
  theme_minimal() +
  coord_fixed()


# plot levels of resistance in different locations over time
# make some timeseries
ir_mtm_africa_ts <- ir_mtm_africa %>%
  group_by(SPECIES, INSECTICIDE_TYPE, INSECTICIDE_CONC, SITE_CODE) %>%
  mutate(
    n_years = n_distinct(YEAR_START),
    .after = ID
  ) %>%
  ungroup() %>%
  filter(n_years > 2) %>%
  arrange(desc(n_years))

# plot the ones with the most data
ir_mtm_africa_ts_plot <- ir_mtm_africa_ts %>% filter(
  INSECTICIDE_TYPE == "Deltamethrin",
  INSECTICIDE_CONC == "0.05%",
  SPECIES %in% c("An. gambiae s.l.", "An. gambiae s.s."))


ir_mtm_africa_ts_plot %>%
  ggplot(
    aes(
      x = YEAR_START,
      y = MORTALITY_ADJUSTED,
      group = SITE_CODE
    )
  ) +
  geom_line() +
  facet_wrap( ~ COUNTRY_NAME)

m <- ir_mtm_africa_ts_plot %>%
  glm(cbind(DIED, SURVIVED) ~ YEAR_START,
      family = stats::binomial,
      data = .)

# tabulate by species, insecticide type, and plot those with multiple observations



# Modelling approach:

# We aim to explain broad-scale spatio-temporal patterns in resistance with a
# Bayesian semi-mechanistic model that captures spatial variation in the
# evolution of resistance as a function of several covariates.

# Model of evolution:

# We use a simple two-genotype (resistant and susceptible) model of evolution,
# and focus our inference on spatial variation in the selection coefficient
# between these genotypes. The selection coefficient is modelled as a function
# of the degree of insecticide exposure, from both malaria control and
# agricultural sources, and also a function of environmental conditions that
# determine the population size, reproduction rate, and survival rates of both
# resistant and susceptible mosquitoes.

# Where the number of mosquitoes in the susceptible genotype for insecticide $I$
# is $S_I$ and the number of the resistant genotype is $R_I$, we focus on
# modelling the evolution over time of the proportion of mosquitoes with the
# susceptible genotype $p_I = S_I / (S_I + R_I)$. The populations of the
# two genotypes reproduce according to fitness parameters $W_I^S$ and $W_I^R$,
# such that after $T$ timesteps (generations or time periods) their
# populations are:
#   $S_{I,T} = \prod_{t=0}{T}(W_I^S S_{I,t})$
#   $R_{I,T} = \prod_{t=0}{T}(W_I^R R_{I,t})$
# so the fraction in the susceptible genotype is:
#   $p_{I,T} = S_{I,T}/(S_{I,T} + R_{I,T})$
# and the susceptible fraction at the next timestep is:
#   $p_{I,t+1} = (W_{I,t}^S S_{I,t}) / (W_{I,t}^S S_{I,t} + W_{I,t}^R R_{I,t})$

# Since this only depends on the ratio of fitness parameters, we remove one
# unidentifiable parameter by defining relative fitnesses of the two genotypes,
# as a ratio to the fitness of the susceptible (wild-type) genotype. Since
# $W_{I,t}^S / W_{I,t}^S = 1$, we focus our attention on the relative fitness of
# the resistance genotype:
#   $w_{I,t} = W_{I,t}^R / W_{I,t}^S$

# Instead of modelling $W_{I,t}^R$, $W_{I,t}^S$ we model $w_{I,t}$ and how it
# varies in space and time as an arbitrary (correlative) function of
# spatio-temporal covariates known to be related to resistance.

# Defining the whole-population average relative fitness at time t as:
#   $\bar{w_{I,t}} = 1 p_{I,t} + w_{I,t} (1 - p_{I,t})
#   $ = p_{I,t} + w_{I,t} (1 - p_{I,t})
# we have the next generation susceptible fraction as:
#   $p_{I,t+1} = \mathfrac{w_{I,t}}{\bar{w_{I,t}}} p_{I,t}$
# and the difference over a single timestep as:
#   $\delta p_{I,t} = p_{I,t+1} - p_{I,t}$
#   $ = \mathfrac{w_{I,t}}{\bar{w_{I,t}}} p_{I,t} -  p_{I,t}$
#   $ = p_{I,t} (\mathfrac{w_{I,t}}{\bar{w_{I,t}}} - 1)$

# with time-varying selection coefficients (given by our correlative model for
# ${w_{I,t}^R}$) we can solve this iteratively for each location to compute the
# fraction susceptible over time, and relate that to resistance bioassay data.

# If we instead assume that the selection pressure for resistance to insecticide
# $I$ is constant, with relative fitness of the resistant genotype $w_{I}$, we
# can model the susceptible fraction in a single spatial location directly as a
# sigmoid-like function of time:
#   $p_{I,T} = p_{I,0} \frac{p_{I,0}}{p_{I,0} + (1 - p_{I,0}) (w_{I}) ^ T}$
# where $p_{I,0}$ is the fraction susceptible at the beginning of the timeseries.

# NOTE: need to confirm I didn't stuff up the maths for this version

# Since we do not know (but have some external information to inform) the time
# at which the original ancestor of the resistant population of each insecticide
# first emerged, we set $t=0$ to the start of 2010 (when resistance data started
# to become more widely available), and learn the time period parameter $tau$
# giving the amount of time from the emergence event (a time at which the
# proportion resistant was at some very low, predetermined fraction $\epsilon =
# 0.0001$) to the start of our timeseries  $t=0$, and hence $p_{I,-\tau} = 1 -
# \epsilon$. Reusing the above equation, we compute $p_{I,0}$ as:
#   $p_{I,0} = p_{I,-\tau} \frac{p_{I,-\tau}}{p_{I,-\tau} + (1 - p_{I,-\tau}) (w_{I}) ^ \tau}$
# We define a mildly informative prior on $\tau$, *a priori* placing it
# uniformly 2000 and 2010.


# Modelling selection for multiple insecticides:

# We model rates of selection for resistance against insecticide type $I$ at
# different times $t$ and locations $i$ as a function of the vector of covariate
# values $x_{i,t}$, which is shared across all insecticide types, multiplied by
# corresponding vector of regression coefficients for that insecticide type
# $\beta_I$:

#   $log(w_{I,t,i}) = x_{i,t} \beta_I$

# We aim to simultaneously infer the rates of spread of resistance to each type
# of insecticide, across several classes by using a three-level hierarchical
# model to share information between insecticide-specific models,
# partially-pooling information both across insecticide classes, and across
# insecticide types within each class. Whilst selection of resistance against
# different insecticides will have different genetic bases, there are likely to
# be common relationships between these selection processes and the spatial
# covaraites we can use as proxies for covariates, resulting in significant
# spatial correlation in resistance across multiple insecticide types.

# Where $I$ is a numerical index denoting the insecticide type (e.g.
# Deltamethrin, Permethrin, Bendiocarb), and C_I is the index to the class to
# which that type of insecticide belongs (e.g. Pyrethroids, Carbamates), we
# model the vector of regression coefficients $beta_I$ hierarchically given
# parameters $\beta_{C_I,k}$ and $\sigma_{C_I,k}^2$ for that class of
# insecticides:
#   $\beta_{I,k} \sim N(\beta_{C_I,k}, \sigma_{C_I,k}^2)$

# where $k$ is the index to the $K$ different covariates for which there is an
# element in $\beta_I$. In turn we model the class-level means from an
# across-class hierarchical distribution:
#   $\beta_{c,k} \sim N(\beta_{k}, \sigma_{k}^2)$

# this doubly-hierarchical structure enables the model to learn common
# relationships between selective pressure and covariates, across all classes
# ($\beta_{k}$) and the degree of similarity across classes $\sigma_{k}^2$, and
# the same again within classes but across types of insecticide. A priori this
# hierarchical structure states that insecticides within the same class are more
# likely to be similar in their relationship with covariates than they are with
# different classes. Hierarchical structures regression coefficientsare also
# regularly found to improve the out-of-sample predictive performance of models,
# our predominant aim in this analysis.


# Combining insecticide concentration data:

# We simultaneously analyse resistance levels using from assay data with
# differing concentrations of each insecticides using a simple dose-response
# relationship. We assume that each genotype has a different LD50 (the dose
# required to kill 50% of the population of the genotype) for each insecticide
# type, with the resistant genotype having a higher LD50 than the susceptible
# genotype for each insecticide $LD50_I^R > LD50_I^S$. For a given
# mixed-genotype population, the population-level LD50 is a weighted average of
# the two:
#   $LD50_{I,t,i} <- LD50_I^S p_{I,t,i} + LD50_I^R (1 - p_{I,t,i})$


# For a given insecticide we model population variability in the dose required
# to kill the mosquito as being normally distributed, with variance $\sigma_I^2$
# of the lethal concentration across the population. We assume that the
# population variance in the lethal dose required is the same for both genotypes
# and therefore constant across all times and places. Given the normal
# distribution described by these two parameters, the proportion of mosquitoes
# killed in a bioassay performed at location $i$ at time $t$, with insecticide concentration $c$ can
# therefore be computed using a probit dose-response model:
#   $q_{I,t,i} = \phi( (c - LD50_{I,t,i}) / \sigma_I)

# Observation model:

# We fit this model to a dataset of resistance bioassay data from discriminating
# concentration bioassays, which expose a known number of mosquitoes to filter
# paper impregnated with insecticide, under standardised conditions and for a
# standardised period of time, and record the number of those exposed that were
# killed. These observations may be modelled with a binomial distribution, given
# the modelled mortality rate described above. However, since the mosquito
# larvae collected form the wild (and then reared) to perform these bioassays
# are collected from a small number of larval habitats, they are very likely to
# be closely related. Consequently, each of the individual mosquito results
# making up a single bioassay result are non-independent, strongly violating the
# binomial assumption. We instead use a betabinomial observation model, which
# assumes that expected mortality rate for the sample (and therefore the sample
# LD50, and sample proportion with the resistance genotype) is itself drawn from
# a population distribution, with mean given by the population mortaliy rate
# $q_{I,t,i}$. This observation model therefore accounts for the inherent
# stochasticity observations due to the small sample size, and the additional
# observation error due to the non-representative larval sampling procedure. The
# observation model on the umber of mosquitoes that died $D_{I, i,t}$, given the
# number exposed $E_{I, i,t}$ is as follows:

#   $a_{I,t,i} = (1 - q_{I,t,i}) / \delta_I^2 - 1 / q_{I,t,i}) * q_{I,t,i} ^ 2$
#   $b_{I,t,i} = a_{I,t,i} (1 / q_{I,t,i} - 1)$
#   $D_{I,i,t} \sim Betabinomial(E_{I, i,t}, a_{I,i,t}, b_{I,i,t})$

# greta implementation

library(greta)

set.seed(2024-06-10)

n_obs <- 20
n_covs <- 5
n_classes <- 3
n_each_type <- sample.int(4, n_classes)
n_types <- sum(n_each_type)

# index to the classes for each type
classes_index <- rep(seq_along(n_each_type), times = n_each_type)

# fake covariates
x <- cbind(1, matrix(rnorm(n_obs * (n_covs - 1)), n_obs, n_covs - 1))

# shrink to scale down predictions
x <- x * 1e-1

# fake observations
df <- data.frame(
  # decimal years
  time = round(runif(n_obs, 0, 20)),
  # aim for 100, but sometime they die early
  n_exposed = 100 - rpois(n_obs, rlnorm(n_obs, -1, 1)),
  # insecticide types
  type = sample.int(n_types, n_obs, replace = TRUE),
  # concentrations
  concentration = c(0.05, 0.1, 0.5)[sample.int(3, n_obs, replace = TRUE)]
)

# hierarchical regression parameters

# between-classes
beta_overall <- normal(0, 1, dim = n_covs)
sigma_overall <- normal(0, 1, dim = n_covs, truncation = c(0, Inf))

# between-types
beta_class_raw <- normal(0, 1, dim = c(n_covs, n_classes))
beta_class_sigma <- sweep(beta_class_raw, 1, sigma_overall, FUN = "*")
beta_class <- sweep(beta_class_sigma, 1, beta_overall, FUN = "+")

# betas for types
sigma_class <- normal(0, 1, dim = n_covs, truncation = c(0, Inf))
beta_type_raw <- normal(0, 1, dim = c(n_covs, n_types))
beta_type_sigma <- sweep(beta_type_raw, 1, sigma_class, FUN = "*")
beta_type <- beta_class[, classes_index] + beta_type_sigma

# multiply through to get log relative fitness of resistance for each
# insecticide type
eta <- x %*% beta_type
# fitnesses <- exp(eta)

# pull out the etas corresponding to data
data_index <- cbind(seq_len(n_obs), df$type)
eta_vec <- eta[data_index]

# model the time of emergence prior to the data timeseries
tau <- uniform(0, 10, dim = n_types)

# get the times from emergence of resistance (tau, prior to data) to each data
# point
# time_pre_data <- sweep(zeros(n_obs, n_types), 2, tau, FUN = "+")
# time <- sweep(time_pre_data, 1, df$time, FUN = "+")
time_vec <- df$time + tau[df$type]

# cumulative_fitnesses <- exp(eta * time)
cumulative_fitnesses_vec <- exp(eta_vec * time_vec)

# compute the fraction susceptible over time against each insecticide
epsilon <- 1e-4
init <- 1 - epsilon
# fraction_susceptible <- init / (init + (1 - init) * cumulative_fitnesses)
fraction_susceptible_vec <- init / (init + (1 - init) * cumulative_fitnesses_vec)

# model the susceptible & resistant LD50s of the different insecticides
LD50_susceptible <- normal(0.05, 0.1, dim = n_types, truncation = c(0, Inf))
LD50_difference <- normal(0.5, 0.1, dim = n_types, truncation = c(0, Inf))
LD50_resistant <- LD50_susceptible + LD50_difference

# get the population-level LD50s for the observations
# LD50_susceptible_weighted <- sweep(fraction_susceptible, 2, LD50_susceptible, FUN = "*")
# LD50_resistant_weighted <- sweep(1 - fraction_susceptible, 2, LD50_resistant, FUN = "*")
# LD50 <- LD50_susceptible_weighted + LD50_resistant_weighted
LD50_susceptible_weighted_vec <- fraction_susceptible_vec *
  LD50_susceptible[df$type]
LD50_resistant_weighted_vec <- (1 - fraction_susceptible_vec) *
  LD50_resistant[df$type]
LD50_vec <- LD50_susceptible_weighted_vec + LD50_resistant_weighted_vec

# and the population mortality rates
population_LD50_sd <- normal(0, 0.1, truncation = c(0, Inf))
# scaled_probit <- sweep(-LD50, 1, df$concentration, FUN = "+") / population_LD50_sd
# population_mortality <- iprobit(scaled_probit)
probit_vec <- (df$concentration - LD50_vec) / population_LD50_sd
population_mortality_vec <- iprobit(probit_vec)

# define betabinomial observation model
observation_extra_error <- normal(0, 0.1, truncation = c(0, Inf))
a <- ((1 - population_mortality_vec) / observation_extra_error ^ 2 - 1 / population_mortality_vec) * population_mortality_vec ^ 2
b <- a * (1 / population_mortality_vec - 1)
died <- beta_binomial(df$n_exposed, a, b)

# now try fitting to data

# check resistance grows and is roughly logistic with sensible values
x <- (1:20) + 10
calculate(1 - fraction_susceptible_vec,
          values = list(
            time_vec = x,
            eta_vec = rep(log(1.6), n_obs)),
          nsim = 1)[[1]][1, , ] -> y
plot(y ~ x, type = "l")

