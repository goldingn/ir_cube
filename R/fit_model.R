# fit model

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
  )

# # summarise the data
# df %>%
#   select(insecticide_class,
#          insecticide_type,
#          concentration,
#          test_type) %>%
#   group_by(insecticide_class,
#            insecticide_type,
#            concentration,
#            test_type) %>%
#   summarise(
#     count = n()
#   ) %>%
#   arrange(desc(count),
#           insecticide_class,
#           insecticide_type)

# fit a model to these data, intercept-only to start

# map the data to the insecticide classes and types
classes <- unique(df$insecticide_class)
types <- unique(df$insecticide_type)

# add in indices to the dataframe to the types and insecticides
df$class_id <- match(df$insecticide_class, classes)
df$type_id <- match(df$insecticide_type, types)

# add in time index
df$time_id <- df$year_start# - 2000

# index to the classes for each type
classes_index <- df %>%
  distinct(type_id, class_id) %>%
  arrange(type_id) %>%
  pull(class_id)

# pull out concentrations for different types
type_concentrations <- df %>%
  select(type_id, concentration) %>%
  group_by(type_id) %>%
  filter(row_number() == 1) %>%
  arrange(type_id) %>%
  pull(concentration)
  
# dimensions
n_obs <- nrow(df)
n_covs <- 1
n_classes <- length(classes)
n_types <- length(types)

# create design matrix
x <- matrix(1, nrow = n_obs, ncol = n_covs)

# hierarchical regression parameters

# between-classes
beta_overall <- normal(0, 1, dim = n_covs) * 10
sigma_overall <- normal(0, 1, dim = n_covs, truncation = c(0, Inf))
sigma_class <- normal(0, 1, dim = n_covs, truncation = c(0, Inf))

# between-types
beta_class_raw <- normal(0, 1, dim = c(n_covs, n_classes))
beta_class_sigma <- sweep(beta_class_raw, 1, sigma_overall, FUN = "*")
beta_class <- sweep(beta_class_sigma, 1, beta_overall, FUN = "+")

# betas for types
beta_type_raw <- normal(0, 1, dim = c(n_covs, n_types))
beta_type_sigma <- sweep(beta_type_raw, 1, sigma_class, FUN = "*")
beta_type <- beta_class[, classes_index] + beta_type_sigma

# multiply through to get log relative fitness of resistance for each
# insecticide type
eta <- x %*% beta_type
# fitnesses <- exp(eta)

# pull out the etas (observations, by insecticides) corresponding to data
data_index <- cbind(seq_len(n_obs), df$type_id)
eta_vec <- eta[data_index]

# model the time of emergence (rescale to enable better priors)
tau_raw <- normal(0, 1, dim = n_types)
tau <- 2000 + tau_raw * 5

# get the times from emergence of resistance (tau) to each data
# point
# time_pre_data <- sweep(zeros(n_obs, n_types), 2, tau, FUN = "+")
# time <- sweep(time_pre_data, 1, df$time, FUN = "+")
time_vec <- df$time_id - tau[df$type_id]

# cumulative_fitnesses <- exp(eta * time)
cumulative_fitnesses_vec <- exp(eta_vec * time_vec)

# compute the fraction susceptible over time against each insecticide
epsilon <- 1e-4
init <- 1 - epsilon
# fraction_susceptible <- init / (init + (1 - init) * cumulative_fitnesses)
fraction_susceptible_vec <- init / (init + (1 - init) * cumulative_fitnesses_vec)

# model the susceptible & resistant LD50s of the different insecticides

# for single concentration data, fix the variance to 1 and fix the LD50 for
# susceptibles to correspond to mortality rate of 99% (98% is resistant)
population_LD50_sd <- 1
LD50_susceptible <- type_concentrations - population_LD50_sd * qnorm(0.99) 

# check mortality prediction
# pnorm( (type_concentrations - LD50_susceptible) / population_LD50_sd)

# define a reasonable prior for the difference in LD50s: resistants have mortality ~1% (95% CI goes up to 25% mortality))
# LD50_resistant_lower <- type_concentrations - population_LD50_sd * qnorm(0.25)
LD50_resistant_mean <- type_concentrations - population_LD50_sd * qnorm(0.01)
# LD50_resistant_sd <- LD50_resistant_mean - LD50_resistant_lower / 1.96
# LD50_difference_mean <- LD50_resistant_mean - LD50_susceptible

# this is a chore
# LD50_difference_raw <- normal(0, 1, dim = n_types)
# LD50_difference <- LD50_difference_mean + LD50_difference_raw * LD50_resistant_sd

# compute LD50 for resistant genotype
# LD50_resistant <- LD50_susceptible + LD50_difference
LD50_resistant <- LD50_resistant_mean

# get the population-level LD50s for the observations
LD50_vec <- fraction_susceptible_vec * LD50_susceptible[df$type_id] +
  (1 - fraction_susceptible_vec) * LD50_resistant[df$type_id]

probit_vec <- (df$concentration - LD50_vec) / population_LD50_sd
population_mortality_vec <- iprobit(probit_vec)

# define betabinomial observation model

# need to. constrain this error for now, until there are covariates to explain
# it better
# observation_extra_error <- normal(0, 0.1, truncation = c(0, Inf))
# observation_extra_error <- normal(0, 0.1, truncation = c(0, 2))
observation_extra_error <- 1
a <- ((1 - population_mortality_vec) / observation_extra_error ^ 2 - 1 / population_mortality_vec) * population_mortality_vec ^ 2
b <- a * (1 / population_mortality_vec - 1)
distribution(df$died) <- beta_binomial(df$mosquito_number, a, b)

# # multiply by 1 to work around a bug in greta for probit probabilities
# distribution(df$died) <- binomial(df$mosquito_number, population_mortality_vec * 1)

# betabinomial is possibly going wacky. Could it be marginalised into the dose-response
# curve? LD50 of the sample has a zero-men Gaussian random variable
# added, implying extra Gaussian noise? That doesn't capture the nonindependece
# though


m <- model(#LD50_difference,
           # observation_extra_error,
           tau,
           beta_overall,
           beta_class_raw,
           beta_type_raw,
           sigma_overall,
           sigma_class)

draws <- mcmc(m)

# check convergence
coda::gelman.diag(draws,
                  autoburnin = FALSE,
                  multivariate = FALSE)

# # compute fitted population mortality rates, and plot with bioassay data
sims <- calculate(population_mortality_vec,
                  fraction_susceptible_vec,
                  values = draws, nsim = 100)
fitted_post_mean <- colMeans(sims$population_mortality_vec[, , 1]) * 100
frac_susc_post_mean <- colMeans(sims$fraction_susceptible_vec[, , 1]) * 100

par(mfrow = c(1, 1))
plot(frac_susc_post_mean ~ df$year_start)


# calculate expected bioassay mortalities over time for each insecticide, and
# overplot the data

n_plot <- 100
times_plot <- seq(1990, max(df$year_start), length.out = n_plot)
time_diff_plot <- sweep(matrix(times_plot, n_plot, n_types), 2, tau, FUN = "-")

eta_plot <- eta[1, ]
cumulative_fitnesses_plot <- exp(sweep(time_diff_plot, 2, t(eta_plot), FUN = "*"))

# compute the fraction susceptible over time against each insecticide
epsilon <- 1e-4
init <- 1 - epsilon
fraction_susceptible_plot <- init / (init + (1 - init) * cumulative_fitnesses_plot)

# get the population-level LD50s for plotting
LD50_susceptible_weighted_plot <- sweep(fraction_susceptible_plot,
                                        2,
                                        LD50_susceptible,
                                        FUN = "*")
LD50_resistant_weighted_plot <- sweep(1 - fraction_susceptible_plot,
                                      2,
                                      LD50_resistant,
                                      FUN = "*")
LD50_plot <- LD50_susceptible_weighted_plot + LD50_resistant_weighted_plot

# convert to proportional mortality expected with random population sampling
probit_plot <- sweep(-LD50_plot, 2, type_concentrations, FUN = "+") / population_LD50_sd
population_mortality_plot <- iprobit(probit_plot)

# # define betabinomial sampling model
# a_plot <- ((1 - population_mortality_plot) / observation_extra_error ^ 2 - 1 / population_mortality_plot) * population_mortality_plot ^ 2
# b_plot <- a_plot * (1 / population_mortality_plot - 1)
# died_plot <- beta_binomial(matrix(1000, n_plot, n_types), a_plot, b_plot)
# assay_mortality_plot <- died_plot / 1000

post_plot_sims <- calculate(
  fraction_susceptible_plot,
  population_mortality_plot,
  # assay_mortality_plot,
  values = draws,
  nsim = 1000)

# summarise these for plotting
frac_susc_post_mean <- apply(post_plot_sims$fraction_susceptible_plot,
                            2:3,
                            mean)
frac_susc_post_ci <- apply(post_plot_sims$fraction_susceptible_plot,
                          2:3,
                          quantile,
                          c(0.025, 0.975))

pop_mort_post_mean <- apply(post_plot_sims$population_mortality_plot,
                            2:3,
                            mean)
pop_mort_post_ci <- apply(post_plot_sims$population_mortality_plot,
                          2:3,
                          quantile,
                          c(0.025, 0.975))

# assay_mort_post_mean <- apply(post_plot_sims$assay_mortality_plot,
#                             2:3,
#                             mean)
# assay_mort_post_ci <- apply(post_plot_sims$assay_mortality_plot,
#                           2:3,
#                           quantile,
#                           c(0.025, 0.975))
# do posterior predictive simulations

par(mfrow = n2mfrow(n_types),
    mar = c(2, 2, 4, 2))
for (i in seq_len(n_types)) {
  plot(frac_susc_post_mean[, i] ~ times_plot,
       ylim = c(0, 1),
       type = "l")
  lines(frac_susc_post_ci[1, , i] ~ times_plot,
        lty = 2)
  lines(frac_susc_post_ci[2, , i] ~ times_plot,
        lty = 2)
  title(main = sprintf("%s\n(%s)",
                       types[i],
                       classes[classes_index[i]]))
  
  # overplot all the data for this type
  df %>%
    filter(type_id == i) %>%
    mutate(mortality = died / mosquito_number) %>%
    select(year_start, mortality) %>%
    points(cex = 0.5, lwd = 0.1)
  
}


par(mfrow = n2mfrow(n_types),
    mar = c(2, 2, 4, 2))
for (i in seq_len(n_types)) {
  plot(pop_mort_post_mean[, i] ~ times_plot,
       ylim = c(0, 1),
       type = "l")
  lines(pop_mort_post_ci[1, , i] ~ times_plot,
        lty = 2)
  lines(pop_mort_post_ci[2, , i] ~ times_plot,
        lty = 2)
  title(main = sprintf("%s\n(%s)",
                       types[i],
                       classes[classes_index[i]]))
  
  # overplot all the data for this type
  df %>%
    filter(type_id == i) %>%
    mutate(mortality = died / mosquito_number) %>%
    select(year_start, mortality) %>%
    points(cex = 0.5, lwd = 0.1)
  
}


# plot posteriors
posteriors <- calculate(tau,
                        # LD50_resistant,
                        eta_plot,
                        values = draws,
                        nsim = 1000)

par(mfrow = n2mfrow(n_types))
for (i in seq_len(n_types)) {
  hist(posteriors$tau[, i, ],
       xlim = c(1980, 2022),
       breaks = 100,
       main = types[i])
  
}

par(mfrow = n2mfrow(n_types))
for (i in seq_len(n_types)) {
  hist(posteriors$eta_plot[, , i],
       breaks = 100,
       xlim = range(posteriors$eta_plot),
       main = types[i])
  
}

# 
# 
# par(mfrow = n2mfrow(n_types))
# for (i in seq_len(n_types)) {
#   vals <- posteriors$LD50_resistant[, i, ]
#   hist(vals,
#        breaks = 100,
#        xlim = range(c(vals, LD50_resistant_mean[i])),
#        main = types[i])
#   abline(v = LD50_resistant_mean[i])
# }

# it seems like there is some non-identifiability in the observation model,
# between the resistant-genotype LD50s, the trend in genotypic resistance
# levels, and the betabinomial residual error. In the diagnostic plots, it looks
# like the model is favouring early resistance, even though the data appear to
# show resistance coming later. Given the chance, it will set the resistant
# LD50s way too low.



# This could also be a coding error between the fitting and secondary prediction
# part, though I can't see anything that could cause it

# It's not obviously due to insufficient variation in beta.

# things to try to fix the model:
# - writing out the model with pen and paper to find nonidentifiabilities
# - automatically finding better initial values
# - removing parts of the model (dropping insecticides, fixing assumptions, etc.) to determine which bits are causing problems


# things to do next:
# - include average coverage as a covariate for resistance
# - set up code for posterior predictive checking
# - set up code for out of sample (future timesteps, spatial blocks) evaluation of model fit.

