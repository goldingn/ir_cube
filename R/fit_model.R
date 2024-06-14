# fit model

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load data
ir_mtm_africa <- readRDS(file = "data/clean/mtm_data.RDS")

# set the start of the timeseries considered
baseline_year <- 2000

df <- ir_mtm_africa %>%
  filter(
    # subset to An. gambiae (s.l./s.s.)
    species %in% c("An. gambiae s.l.", "An. gambiae s.s."),
    # subset to WHO tube tests (most of data)
    test_type == "WHO tube test",
    # drop the minor classes: pyrroles and neonicotinoids
    insecticide_class %in% c("Pyrethroids",
                             "Carbamates",
                             "Organochlorines",
                             "Organophosphates")
  ) %>%
  mutate(
    # convert concentrations into numeric values
    concentration = as.numeric(str_remove(insecticide_conc, "%")),
    # create a time index, relative to the baseline year
    time = year_start - baseline_year
  ) %>%
  filter(
    # drop Dieldrin as it has two concentrations in the dataset, and those with
    # fewer than 200 observations (manually as I'm being lazy)
    !(insecticide_type %in% c("Dieldrin",
                              "Carbosulfan",
                              "Cyfluthrin",
                              "Etofenprox",
                              "Propoxur")),
    # drop any from before when we have data on net coverage
    year_start >= baseline_year
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

# # doubly hierarchical version
# 
# # between-classes
# beta_overall <- normal(0, 1, dim = n_covs) * 10
# sigma_overall <- normal(0, 1, dim = n_covs, truncation = c(0, Inf))
# sigma_class <- normal(0, 1, dim = n_covs, truncation = c(0, Inf))
# 
# # between-types
# beta_class_raw <- normal(0, 1, dim = c(n_covs, n_classes))
# beta_class_sigma <- sweep(beta_class_raw, 1, sigma_overall, FUN = "*")
# beta_class <- sweep(beta_class_sigma, 1, beta_overall, FUN = "+")
# 
# # betas for types
# beta_type_raw <- normal(0, 1, dim = c(n_covs, n_types))
# beta_type_sigma <- sweep(beta_type_raw, 1, sigma_class, FUN = "*")
# beta_type <- beta_class[, classes_index] + beta_type_sigma


# simpler single-level hierarchy
beta_overall <- normal(0, 1, dim = n_covs) * 10
sigma_overall <- normal(0, 1, dim = n_covs, truncation = c(0, Inf))

# hierarchical decentring implementation
beta_type_raw <- normal(0, 1, dim = c(n_covs, n_types))
beta_type_sigma <- sweep(beta_type_raw, 1, sigma_overall, FUN = "*")
beta_type <- sweep(beta_type_sigma, 1, beta_overall, FUN = "+")


# remove hierarchy for debugging
# beta_type <- beta_type_raw


# multiply through to get log relative fitness of resistance for each
# insecticide type

# this way is marginally faster and probably quite a lot less memory intensive
# than the matrix multiply
beta_mat <- t(beta_type)[df$type_id,]
eta_mat <- x * beta_mat
eta_vec <- greta:::rowSums(eta_mat)

# pull out the etas (observations, by insecticides) corresponding to data
# eta <- x %*% beta_type
# data_index <- cbind(seq_len(n_obs), df$type_id)
# eta_vec <- eta[data_index]

# switch this to just modelling the fraction susceptible in 2000

# relax this prior - DDT was pretty resistant by then
init_fraction_susceptible <- normal(1, 0.01,
                                    dim = n_types,
                                    truncation = c(0.9, 1))
# init_sims <- calculate(init_fraction_susceptible[1], nsim = 1e5)[[1]]
# hist(init_sims, xlim = c(0.8, 1))
# mean(init_sims)
# mean(init_sims < 0.9)
# quantile(init_sims, 0.25)


# compute the fraction susceptible over time against each insecticide

# more stable version of this:
#   fitness = exp(eta)
#   cumulative_fitness = fitness ^ df$time
cumulative_fitnesses_vec <- exp(eta_vec * df$time)

init_vec <- init_fraction_susceptible[df$type_id]
fraction_susceptible_vec <- init_vec / 
  (init_vec + (1 - init_vec) * cumulative_fitnesses_vec)

# model the susceptible & resistant LD50s of the different insecticides
# (skip this until we have multiple concentrations)

# # for single concentration data, fix the variance to 1 and fix the LD50 for
# # susceptibles to correspond to mortality rate of 99% (98% is resistant)
# population_LD50_sd <- 1
# LD50_susceptible <- type_concentrations - population_LD50_sd * qnorm(0.99) 
# 
# # define a reasonable prior for the difference in LD50s: resistants have mortality ~1% (95% CI goes up to 25% mortality))
# # LD50_resistant_lower <- type_concentrations - population_LD50_sd * qnorm(0.25)
# LD50_resistant_mean <- type_concentrations - population_LD50_sd * qnorm(0.01)
# # LD50_resistant_sd <- LD50_resistant_mean - LD50_resistant_lower / 1.96
# # LD50_difference_mean <- LD50_resistant_mean - LD50_susceptible
# 
# # this is a chore
# # LD50_difference_raw <- normal(0, 1, dim = n_types)
# # LD50_difference <- LD50_difference_mean + LD50_difference_raw * LD50_resistant_sd
# 
# # compute LD50 for resistant genotype
# # LD50_resistant <- LD50_susceptible + LD50_difference
# LD50_resistant <- LD50_resistant_mean
# 
# # get the population-level LD50s for the observations
# LD50_vec <- fraction_susceptible_vec * LD50_susceptible[df$type_id] +
#   (1 - fraction_susceptible_vec) * LD50_resistant[df$type_id]
# 
# probit_vec <- (df$concentration - LD50_vec) / population_LD50_sd
# population_mortality_vec <- iprobit(probit_vec)

# for now, map straight from genotype to mortality
population_mortality_vec <- fraction_susceptible_vec

# define betabinomial observation model

# need to constrain this error for now, until there are covariates to explain it
# better
# observation_extra_error <- normal(0, 0.1, truncation = c(0, Inf))
# a <- ((1 - population_mortality_vec) / observation_extra_error ^ 2 - 1 / population_mortality_vec) * population_mortality_vec ^ 2
# b <- a * (1 / population_mortality_vec - 1)
# distribution(df$died) <- beta_binomial(df$mosquito_number, a, b)

# just use binomial for now
distribution(df$died) <- binomial(df$mosquito_number, population_mortality_vec)

m <- model(
  beta_overall,
  init_fraction_susceptible,
  beta_type_raw,
  sigma_overall
)

system.time(
  draws <- mcmc(m)
)

# check convergence
coda::gelman.diag(draws,
                  autoburnin = FALSE,
                  multivariate = FALSE)

# compute fitted population mortality rates, and plot with bioassay data
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
type_idx <- 1
times_plot <- seq(1990, max(df$year_start), length.out = n_plot)
time_diff_plot <- sweep(matrix(times_plot, n_plot, n_types), 2, 2000, FUN = "-")

eta_plot <- beta_type
cumulative_fitnesses_plot <- exp(sweep(time_diff_plot, 2, t(eta_plot), FUN = "*"))

# compute the fraction susceptible over time against each insecticide
init_mat <- sweep(zeros(n_plot, n_types), 2, init_fraction_susceptible, FUN = "+")
fraction_susceptible_plot <- init_mat /
  (init_mat + (1 - init_mat) * cumulative_fitnesses_plot)

# # get the population-level LD50s for plotting
# LD50_susceptible_weighted_plot <- sweep(fraction_susceptible_plot,
#                                         2,
#                                         LD50_susceptible,
#                                         FUN = "*")
# LD50_resistant_weighted_plot <- sweep(1 - fraction_susceptible_plot,
#                                       2,
#                                       LD50_resistant,
#                                       FUN = "*")
# LD50_plot <- LD50_susceptible_weighted_plot + LD50_resistant_weighted_plot
# 
# # convert to proportional mortality expected with random population sampling
# probit_plot <- sweep(-LD50_plot, 2, type_concentrations, FUN = "+") / population_LD50_sd
# population_mortality_plot <- iprobit(probit_plot)
population_mortality_plot <- fraction_susceptible_plot

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
posteriors <- calculate(init_fraction_susceptible,
                        beta_type,
                        values = draws,
                        nsim = 1000)

par(mfrow = n2mfrow(n_types))
for (i in seq_len(n_types)) {
  hist(posteriors$init_fraction_susceptible[, i, ],
       xlim = c(0.9, 1),
       breaks = 100,
       main = types[i])
  
}

par(mfrow = n2mfrow(n_types))
for (i in seq_len(n_types)) {
  hist(posteriors$beta_type[, , i],
       breaks = 100,
       main = types[i])
  
}

# things to do next:
# - add in the temporally-static LLIN and crop covariates
# - add in code for spatial predictions
# - add in time-varying nets
# - set up code for posterior predictive checking
# - set up code for out of sample (future timesteps, spatial blocks) evaluation of model fit.
# - add in more discriminating bioassay data
# - add in intensity bioassay data (with LD50 model, including possibility of complete but imperfect resistance)
