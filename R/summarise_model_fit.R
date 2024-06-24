# summarise model

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load the fitted model objects here, to set up predictions
load(file = "temporary/fitted_model.RData")

# load covariate rasters (need to reload these even if restoring workspace
# because pointers)
covs_flat <- rast("data/clean/flat_covariates.tif")
mask <- covs_flat[[1]] * 0

# summarise the covariate effect sizes (at insecticide class level)
effect_sizes <- summary(calculate(exp(beta_class[, 1]), values = draws))$statistics[, c("Mean", "SD")]
rownames(effect_sizes) <- c("itn_coverage", names(covs_flat))
round(effect_sizes, 2)

# get RMSE and MAE for observed data and posterior mean within-sample predictions
pop_susc_sim <- calculate(population_mortality_vec,
                          values = draws,
                          nsim = 1e3)
pop_susc_post_mean <- colMeans(pop_susc_sim$population_mortality_vec[, , 1])
obs <- df$died / df$mosquito_number
rmse <- sqrt(mean((pop_susc_post_mean - obs) ^ 2))
mae <- mean(abs(pop_susc_post_mean - obs))

# get posterior predictive simulations of observations
died_sim <- betabinomial_p_rho(N = df$mosquito_number,
                               p = population_mortality_vec,
                               rho = rho_classes[df$class_id])
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

# some statistical evidence of misspecification, but it's a very large sample
# size and a very small deviation
plot(dharma)

# looks relatively uniform for each insecticide
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


# make posterior predictions of Africa-wide average susceptibility (over
# locations with any bioassays) for all years with the average sample size

# get average total sample size per year, per insecticide, over all samples
average_sample_sizes <- df %>%
  group_by(year_start, type_id) %>%
  summarise(
    mosquito_number = sum(mosquito_number),
    .groups = "drop"
  ) %>%
  group_by(type_id) %>%
  summarise(
    mosquito_number = mean(mosquito_number)
  ) %>%
  arrange(type_id) %>%
  pull(mosquito_number)

# get average annual total sample size in unique cell, per insecticide
cell_sample_sizes <- df %>%
  group_by(cell_id, year_start, type_id) %>%
  # get annual sums 
  summarise(
    mosquito_number = sum(mosquito_number),
    .groups = "drop"
  ) %>%
  # now collapse to cell and insecticide type (average over years to make it
  # constant for plotting)
  group_by(cell_id, type_id) %>%
  summarise(
    # this isn't padded, to to get annual mean we need to divide by all years,
    # not observed years
    mosquito_number = sum(mosquito_number) / n_times,
    .groups = "drop"
  )

# we already have a greta array for susceptibility in all unique cells in the
# dataset, for all insecticides and all years: dynamic_cells$all_states

# for each year and each insecticide, compute the weighted (by annual average of
# the total number of mosquitoes at that cell, per insecticide) mean of the
# susceptibility fractions, and simulate the expected fraction in observed
# summary stats.

sample_size_mat <- expand_grid(
  cell_id = seq_len(n_unique_cells),
  type_id = seq_len(n_types)
) %>%
  left_join(cell_sample_sizes,
            by = c("cell_id", "type_id")) %>%
  replace_na(list(mosquito_number = 0)) %>%
  # make into a matrix and squash names
  pivot_wider(
    names_from = type_id,
    values_from = mosquito_number,
    names_prefix = "type_"
  ) %>%
  select(-cell_id) %>%
  as.matrix() %>%
  `dimnames<-`(NULL)

# get weights
weights_mat <- sweep(sample_size_mat,
                    2,
                    colSums(sample_size_mat),
                    "/")

# expand out into an array, replicating by year
weights_array <- array(rep(c(weights_mat), n_times),
                       dim = c(n_unique_cells, n_types, n_times))

# check dimensions are the right way around
# identical(weights_array[, , 1],  weights_mat)
# all(weights_array[1, 1, ] == weights_array[1, 1, 1])

# multiply through and sum to get average susceptibilities
weighted_susc_array <- dynamic_cells$all_states * weights_array
overall_susc <- apply(weighted_susc_array, 2:3, "sum")

# expand these out into a vector by insecticides and years, to later compile
# sims
indices <- expand_grid(
  type_id = seq_len(n_types),
  time_id = seq_len(n_times)
)

overall_susc_vec <- overall_susc[as.matrix(indices)]

# relabel - this is population-level mortality
overall_pop_mort_vec <- overall_susc_vec

# get posterior samples of these
sims <- calculate(overall_pop_mort_vec,
                  values = draws,
                  nsim = 1e3)

type_years <- indices %>%
  mutate(
    year = baseline_year + time_id - 1,
    insecticide = types[type_id],
    row = row_number()
  )  %>%
  select(row, year, insecticide)

insecticides_plot <- tibble(
  insecticide = types,
  class = classes[classes_index]
) %>%
  arrange(desc(class), insecticide) %>%
  filter(
    !(insecticide %in% c("DDT")) 
  ) %>%
  pull(insecticide)

# make these into tibbles of sims, and summarise
pop_mort_sry <- sims$overall_pop_mort_vec[, , 1] %>%
  t() %>%
  as_tibble() %>%
  mutate(
    row = row_number(),
    .before = everything()
  ) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_prefix = "V",
    names_to = "sim",
    values_to = "susceptibility"
  ) %>%
  group_by(row) %>%
  summarise(
    susc_pop_mean = mean(susceptibility),
    susc_pop_lower = quantile(susceptibility, 0.025),
    susc_pop_upper = quantile(susceptibility, 0.975),
  ) %>%
  left_join(type_years,
            by = "row") %>%
  select(-row) %>%
  filter(
    insecticide %in% insecticides_plot
  ) %>%
  # set the plotting order of the insecticides
  mutate(
    insecticide = factor(insecticide,
                         levels = insecticides_plot)
  )

# now plot these, along with the raw data
df_country_plot <- df %>%
  group_by(country_name, year_start, insecticide_type) %>%
  summarise(
    died = sum(died),
    mosquito_number = sum(mosquito_number),
    .groups = "drop"
  ) %>%
  mutate(
    Susceptibility = died / mosquito_number
  ) %>%
  rename(
    year = year_start,
    country = country_name,
    insecticide = insecticide_type
  ) %>%
  filter(
    insecticide %in% insecticides_plot
  ) %>%
  mutate(
    insecticide = factor(insecticide,
                         levels = insecticides_plot)
  )

df_overall_plot <- df_country_plot %>%
  group_by(year, insecticide) %>%
  summarise(
    died = sum(died),
    mosquito_number = sum(mosquito_number),
    .groups = "drop"
  ) %>%
  mutate(
    Susceptibility = died / mosquito_number
  )

pop_mort_sry %>%
  ggplot(
  aes(
    x = year,
    fill = insecticide
  )
) +
  geom_ribbon(
    aes(
      ymax = susc_pop_upper,
      ymin = susc_pop_lower,
    )
  ) +
  geom_point(
    data = df_overall_plot,
    mapping = aes(
      y = Susceptibility,
      group = "none",
      size = mosquito_number
    ),
    shape = 21,
  ) +
  facet_wrap(~insecticide,
             ncol = 4) +
  guides(
    size = guide_legend(title = "No. tested")
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)) +
  scale_size_continuous(
    labels = scales::number_format(
      accuracy = 1000,
      big.mark = ","
    )
  ) +
  scale_fill_discrete(
    direction = -1,
    guide = FALSE) +
  xlab("") +
  coord_cartesian(xlim = c(1998, 2022)) +
  theme_minimal() +
  ggtitle(
    "Modelled growth of resistance to single insecticides across Africa"
  )

ggsave("figures/all_africa_obs_pred.png",
       bg = "white",
       width = 12,
       height = 5)

# Annual aggregated bioassay results (circles) and predicted population-level
# resistance across the sites where the samples were collected (bold colour
# band; 95%CI)



# plot predicted trends and numbers of nets at a few unique locations

# look up the net coverage for different cells, find those with low, medium, and
# high net coverage
set.seed(1)
exemplar_cells <- all_extract %>%
  group_by(cell_id) %>%
  summarise(
    net_coverage = mean(net_coverage),
    .groups = "drop"
  ) %>%
  # then pull the quartile threshold values
  mutate(
    which = case_when(
      near(net_coverage, quantile(net_coverage, 0.1), tol = 1e-4) ~ "low",
      near(net_coverage, quantile(net_coverage, 0.5), tol = 1e-4) ~ "medium",
      near(net_coverage, quantile(net_coverage, 0.9), tol = 1e-4) ~ "high",
      .default = NA
    )
  ) %>%
  filter(
    !is.na(which)
  ) %>%
  # randomly pick a cell from each group
  arrange(runif(n())) %>%
  filter(!is.na(which)) %>%
  # get the first cell in each level of net coverage
  group_by(which) %>%
  summarise(
    cell_id = cell_id[1],
    .groups = "drop"
  ) %>%
  mutate(
    which = factor(which, levels = c("low", "medium", "high"))
  ) %>%
  # find the coordinates and reverse geocode them
  mutate(
    cell = unique_cells[cell_id],
  ) %>%
  bind_cols(
    xyFromCell(mask, .$cell)
  ) %>%
  reverse_geocode(
    lat = y,
    long = x,
    method = 'osm',
    full_results = TRUE
  )

# find a short name for these places
place_lookup <- exemplar_cells %>%
  mutate(
    precise_place = str_split_i(address, ",", 1),
    place = paste(precise_place, country, sep = ", ")
  ) %>%
  select(
    cell_id,
    which,
    place)

# pull the IR and ITN timeseries for these, and plot

# do predictions of the effective resistance to ITNS
ingredient_weights <- readRDS("temporary/ingredient_weights.RDS")

# now do predictions for these, for deltamethrin
ingredient_ids <- match(names(ingredient_weights), types)
# deltamethrin_id <- match("Deltamethrin", types)

# compute a matrix of effective susceptibilities over these cells and years
effective_susc <- zeros(nrow(place_lookup), 1, n_times)
for (i in seq_along(ingredient_weights)) {
  ingredient_susc <- dynamic_cells$all_states[place_lookup$cell_id,
                                              ingredient_ids[i], ]
  effective_susc <- effective_susc + ingredient_susc * ingredient_weights[[i]]
}

pred_index <- expand_grid(
  cell_id = seq_len(nrow(place_lookup)),
  type_id = 1,
  year_id = seq_len(n_times)
) %>%
  as.matrix()

# subset the susceptibilities to these cells and insecitcides
pred_vec <- effective_susc[pred_index]
pred_vec_sims <- calculate(pred_vec, values = draws, nsim = 1e3)

ir_plot <- all_extract %>%
  filter(
    cell_id %in% place_lookup$cell_id
  ) %>%
  left_join(place_lookup,
            by = "cell_id") %>%
  mutate(
    year = baseline_year + year_id - 1,
    post_mean = colMeans(pred_vec_sims$pred_vec[, , 1]),
    post_lower = apply(pred_vec_sims$pred_vec[, , 1], 2, quantile, 0.025),
    post_upper = apply(pred_vec_sims$pred_vec[, , 1], 2, quantile, 0.975),
  ) %>%
  ggplot(
    aes(
      x = year,
      y = post_mean,
      ymax = post_upper,
      ymin = post_lower,
      group = place
    )
  ) +
  facet_wrap(~place) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent) +
  ylab("Susceptibility to LLIN insecticides") +
  xlab("") +
  geom_ribbon(fill = "#56B1F7") +
  theme_minimal()

itns_plot <- all_extract %>%
  filter(
    cell_id %in% place_lookup$cell_id
  ) %>%
  left_join(place_lookup,
            by = "cell_id") %>%
  mutate(
    year = baseline_year + year_id - 1,
  ) %>%
  ggplot(
    aes(
      x = year,
      y = net_coverage,
      group = place
    )
  ) +
  facet_wrap(~place) +
  scale_y_continuous(limits = c(0, 1),
                     labels = scales::percent) +
  ylab("LLIN use") +
  xlab("") +
  geom_line(
    colour = grey(0.5),
    linewidth = 1.2
  ) + 
  theme_minimal() +
  # suppress facet labels for bottom row
  theme(
    strip.text.x = element_blank()
  )

# # patchwork is bugging out with recent ggplot, use 3.4.4:
# remotes::install_version("ggplot2", version = "3.4.4", repos = "http://cran.us.r-project.org")

ir_plot / itns_plot

ggsave("figures/exemplar_itn_susc.png",
       bg = "white",
       width = 8,
       height = 5)


insecticides_plot_small <- c("Deltamethrin",
                             "Permethrin",
                             "Alpha-cypermethrin")

# find some locations with lots of bioassay data for the pyrethroids and plot
# predictions and data for these
locations_plot <- df %>%
  filter(
    insecticide_type %in% insecticides_plot_small
  ) %>%
  group_by(cell) %>%
  filter(
    n() >= 25,
    n_distinct(year_start) >= 5
  ) %>%
  select(
    country_name,
    cell_id,
    cell
  ) %>%
  distinct() %>%
  bind_cols(
    xyFromCell(mask, .$cell)
  ) %>%
  rename(
    longitude = x,
    latitude = y
  ) %>%
  # geocode then find more interpretable names (google to see if this is how
  # they are referred to in IR papers)
  reverse_geocode(
    lat = latitude,
    long = longitude,
    method = 'osm',
    full_results = TRUE
  ) %>%
  mutate(
    precise_place = str_split_i(address, ",", 1),
    # tidy up some of these where possible
    precise_place = case_when(
      grepl("Soumousso", address) ~ "Soumousso",
      grepl("Busia", address) ~ "Busia",
      grepl("Homa Bay", address) ~ "Homa Bay",
      grepl("Pitoa", address) ~ "Pitoa",
      grepl("Garoua", address) ~ "Garoua",
      .default = precise_place
    ),
    place = paste(precise_place, country_name, sep = ", ")
  ) %>%
  select(
    place,
    latitude,
    longitude,
    cell_id
  )

# pull these out for plotting
preds_plot_setup <- expand_grid(
  cell_id = locations_plot$cell_id,
  year_start = baseline_year:max(df$year_start),
  insecticide_type = types
) %>%
  mutate(
    year_id = year_start - baseline_year + 1,
    type_id = match(insecticide_type, types)
  )

index_plot <- preds_plot_setup %>%
  select(cell_id,
         type_id,
         year_id) %>%
  as.matrix()

fraction_susceptible_plot <- dynamic_cells$all_states[index_plot]
population_mortality_plot <- fraction_susceptible_plot

# simulate mortalities under binomial sampling
sample_size_plot <- 100
binomial_died <- binomial(sample_size_plot, population_mortality_plot)
binomial_mortality <- binomial_died / sample_size_plot

rho_index <- classes_index[index_plot[, "type_id"]]

# simulate mortalities under betabinomial sampling
betabinomial_died <- betabinomial_p_rho(N = sample_size_plot,
                                        p = population_mortality_plot,
                                        rho = rho_classes[rho_index])
betabinomial_mortality <- betabinomial_died / sample_size_plot

sims <- calculate(population_mortality_plot,
                  binomial_mortality,
                  betabinomial_mortality,
                  values = draws,
                  nsim = 2000)

# posterior mean mortality rate, and intervals for the population value
# (posterior uncertainty), and posterior predictive intervals (posterior
# uncertainty and sampling error) for observed mortality from binomial sampling,
# and from betabinomial sampling
pop_mort_mean <- colMeans(sims$population_mortality_plot[, , 1])
pop_mort_ci <- apply(sims$population_mortality_plot[, , 1],
                     2,
                     quantile,
                     c(0.025, 0.975))
binomial_mort_ci <- apply(sims$binomial_mortality[, , 1],
                          2,
                          quantile,
                          c(0.025, 0.975))
betabinomial_mort_ci <- apply(sims$betabinomial_mortality[, , 1],
                              2,
                              quantile,
                              c(0.025, 0.975))

preds_plot <- preds_plot_setup %>%
  mutate(
    Susceptibility = pop_mort_mean,
    pop_lower = pop_mort_ci[1, ],
    pop_upper = pop_mort_ci[2, ],
    binomial_lower = binomial_mort_ci[1, ],
    binomial_upper = binomial_mort_ci[2, ],
    betabinomial_lower = betabinomial_mort_ci[1, ],
    betabinomial_upper = betabinomial_mort_ci[2, ],
  ) %>%
  filter(
    insecticide_type %in% insecticides_plot_small
  ) %>%
  left_join(
    locations_plot,
    by = "cell_id"
  ) %>%
  mutate(
    insecticide_type = factor(insecticide_type,
                              levels = insecticides_plot_small)
  )

points_plot <- df %>%
  filter(
    cell_id %in% locations_plot$cell_id,
    insecticide_type %in% insecticides_plot_small
  ) %>%
  mutate(
    Susceptibility = died / mosquito_number
  ) %>%
  left_join(
    locations_plot,
    by = "cell_id"
  ) %>%
  mutate(
    insecticide_type = factor(insecticide_type,
                              levels = insecticides_plot_small)
  )
  

# set the colours for these insecticides
colour_types <- scales::hue_pal(direction = -1)(8)
types_plot_id <- match(insecticides_plot_small, insecticides_plot)
colours_plot <- colour_types[types_plot_id]

# plot these, then add cell data over the top
preds_plot %>%
  ggplot(
    aes(
      x = year_start,
      y = Susceptibility,
      group = place,
      fill = insecticide_type
    )
  ) +
  # betabinomial posterior predictive intervals on bioassay data (captures
  # sample size and non-independence effect)
  geom_ribbon(
    aes(
      ymax = betabinomial_lower,
      ymin = betabinomial_upper,
    ),
    colour = grey(0.4),
    linewidth = 0.25,
    linetype = 2,
    alpha = 0.1
  ) +
  # binomial posterior predictive intervals on bioassay data (captures sample
  # size but assumes independence, so underestimates variance)
  geom_ribbon(
    aes(
      ymax = binomial_lower,
      ymin = binomial_upper,
    ),
    alpha = 0.4
  ) +
  # credible intervals on population-level proportion (our best guess at the
  # 'truth')
  geom_ribbon(
    aes(
      ymax = pop_lower,
      ymin = pop_upper,
    ),
  ) +
  geom_point(
    aes(
      x = year_start,
      y = Susceptibility,
      group = place,
      size = mosquito_number
    ),
    shape = 21,
    data = points_plot,
  ) +
  scale_fill_manual(
    values = colours_plot,
    guide = FALSE
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  scale_size_binned(
    range = c(1, 4),
    n.breaks = 6
  ) +
  facet_grid(insecticide_type ~ place) +
  guides(
    size = guide_legend(title = "No. tested")
  ) +
  xlab("") +
  coord_cartesian(
    xlim = c(1998, 2022),
    ylim = c(0, 1)
  ) +
  theme_minimal() +
  ggtitle(
    "Modelled population-level susceptibility and bioassay results at heavily sampled locations"
  )

ggsave("figures/fit_subset.png",
       bg = "white",
       scale = 0.8,
       width = 16,
       height = 7)

# these are the only sites with at least 25 pyrethroid bioassay datapoints,
# spanning at least 5 years

