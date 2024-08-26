# make figure of model predictions over time for all of Africa overlaid with raw
# bioassay data, aggregated up at various levels.

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load the fitted model objects here, to set up predictions
load(file = "temporary/fitted_model.RData")

# load time-varying net coverage data and flatten it
nets_cube <- rast("data/clean/nets_per_capita_cube.tif")
years_sub <- 2010:2022
nets_cube_sub <- nets_cube[[paste0("nets_", years_sub)]]
nets_flat <- app(nets_cube_sub, "mean")

# set the start of the timeseries considered
baseline_year_data <- 1990

# set colours
pyrethroid_blue <- "#56B1F7"

# get all the pyrethroids
pyrethroids <- tibble(
  insecticide = types,
  class = classes[classes_index]
) %>%
  arrange(desc(class), insecticide) %>%
  filter(class == "Pyrethroids") %>%
  pull(insecticide)

# subset to pyrethroids
df_pyrethroids <- df %>%
  filter(
    insecticide_class == "Pyrethroids",
  ) %>%
  rename(
    year = year_start,
    country = country_name,
    insecticide = insecticide_type
  )

# data aggregated over all of the pyrethroids
df_pyrethroids_africa_plot <- df_pyrethroids %>%
  group_by(year) %>%
  summarise(
    died = sum(died),
    mosquito_number = sum(mosquito_number),
    .groups = "drop"
  ) %>%
  mutate(
    Susceptibility = died / mosquito_number
  )

# compute predicted population-level resistance to these pyrethroids, weighted
# over the numbers of samples per location, per insecticide type

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
overall_pop_mort_vec <- overall_susc_vec

# now compute a weighted sum over the pyrethroids (weights given by the numbers
# of mosquitoes tested for each pyrethroid, over the whole timeseries)
pyrethroid_weights <- df %>% 
  group_by(
    type_id,
    insecticide_type 
  ) %>%
  summarise(
    mosquito_number = sum(mosquito_number),
    .groups = "drop"
  ) %>%
  mutate(
    is_a_pyrethroid = insecticide_type %in% pyrethroids,
    mosquito_number = mosquito_number * as.numeric(is_a_pyrethroid),
    weight = mosquito_number / sum(mosquito_number)
  ) %>%
  arrange(type_id) %>%
  pull(weight)

pyrethroid_susc <- as_data(t(pyrethroid_weights)) %*% overall_susc

# get posterior samples of these
sims <- calculate(overall_pop_mort_vec,
                  pyrethroid_susc,
                  values = draws,
                  nsim = 1e3)

# posterior summaries of pyrethroid susceptibility by year
pop_mort_sry_pyrethroid <- sims$pyrethroid_susc[,  1, ] %>%
  t() %>%
  as_tibble() %>%
  mutate(
    time_id = row_number(),
    .before = everything()
  ) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_prefix = "V",
    names_to = "sim",
    values_to = "susceptibility"
  ) %>%
  group_by(time_id) %>%
  summarise(
    susc_pop_mean = mean(susceptibility),
    susc_pop_lower = quantile(susceptibility, 0.025),
    susc_pop_upper = quantile(susceptibility, 0.975),
  ) %>%
  mutate(
    year = baseline_year + time_id - 1
  )

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
  # filter(
  #   !(insecticide %in% c("DDT")) 
  # ) %>%
  pull(insecticide)

insecticide_type_labels <- sprintf("%s) %s%s",
                                   LETTERS[1 + seq_along(insecticides_plot)],
                                   insecticides_plot,
                                   ifelse(insecticides_plot %in% pyrethroids,
                                          "*",
                                          ""))

insecticides_plot_lookup <- tibble(
  insecticide = insecticides_plot,
  insecticide_type_label = insecticide_type_labels
)

# summarise posteriors of plots against all insecticides
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
  left_join(
    insecticides_plot_lookup,
    by = "insecticide"
  )

df_overall_plot <- df %>%
  group_by(year_start,
           insecticide_type) %>%
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
    insecticide = insecticide_type
  ) %>%
  filter(
    insecticide %in% insecticides_plot
  ) %>%
  left_join(
    insecticides_plot_lookup,
    by = "insecticide"
  )

# make them share a point size legend
size_limits <- c(min(df_overall_plot$mosquito_number),
                 max(df_pyrethroids_africa_plot$mosquito_number))


# plot all-pyrethroid figure
pyrethroid_fig <- pop_mort_sry_pyrethroid %>%
  ggplot(
    aes(
      x = year
    )
  ) +
  geom_ribbon(
    aes(
      ymax = susc_pop_upper,
      ymin = susc_pop_lower,
    ),
    data = pop_mort_sry_pyrethroid,
    fill = pyrethroid_blue
  ) +
  geom_point(
    aes(
      y = Susceptibility,
      size = mosquito_number
    ),
    shape = 21,
    fill = pyrethroid_blue,
    data = df_pyrethroids_africa_plot
  ) +
  guides(
    size = "none"
  ) +
  facet_wrap(~ "A) All pyrethroids*") +
  scale_y_continuous(
    labels = scales::percent) +
  scale_size_continuous(
    labels = scales::number_format(accuracy = 1000, big.mark = ","),
    limits = size_limits
  ) +
  xlab("") +
  ylab("Susceptibility") +
  coord_cartesian(xlim = c(1998, 2022),
                  ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    strip.text.x = element_text(hjust = 0)
  )

all_insecticides_fig <- pop_mort_sry %>%
  ggplot(
    aes(
      x = year,
      fill = insecticide_type_label
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
  facet_wrap(~insecticide_type_label,
             ncol = 3) +
  guides(
    size = guide_legend(title = "No. tested")
  ) +
  scale_y_continuous(
    labels = scales::percent,
    breaks = c(0, 1),
    limits = c(0, 1)) +
  scale_x_continuous(breaks = c(2000, 2020)) +
  scale_size_continuous(
    labels = scales::number_format(
      accuracy = 1000,
      big.mark = ","
    ),
    limits = size_limits
  ) +
  scale_fill_discrete(
    direction = -1,
    guide = FALSE) +
  xlab("") +
  # suppress ylab, as it is on the other panel
  ylab("") +
  coord_cartesian(xlim = c(1998, 2022)) +
  theme_minimal() +
  theme(
    strip.text.x = element_text(hjust = 0)
  )

# use patchwork to set up the multi-panel plot
pyrethroid_fig + all_insecticides_fig

ggsave("figures/all_africa_pred_data.png",
       bg = "white",
       scale = 0.8,
       width = 14,
       height = 6)

# Annual aggregated bioassay results (circles) and predicted population-level
# resistance across the sites where the samples were collected (bold colour
# band; 95%CI)


# add the low, medium, high LLINs underneath
# find how much data by different coverages

net_coverage_cell_lookup <- df %>%
  group_by(
    cell_id, cell
  ) %>%
  # get one record per observed cell
  slice(1) %>%
  ungroup() %>%
  mutate(
    net_coverage = terra::extract(nets_flat, pull(., cell)),
    net_coverage_class = case_when(
      net_coverage < 0.3 ~ "A) Low use",
      .default = "B) High use")
  ) %>%
  select(
    cell,
    cell_id,
    net_coverage_class
  )

# extract net use over time at all observation locations
df_net_use <- all_extract %>%
  left_join(
    net_coverage_cell_lookup,
    by = "cell_id"
  ) %>%
  group_by(
    year_id,
    net_coverage_class
  ) %>%
  summarise(
    lower = quantile(nets, 0.25),
    upper = quantile(nets, 0.75),
    mean = mean(nets),
    .groups = "drop"
  ) %>%
  mutate(
    year = year_id + baseline_year - 1
  )


# get LLIN-effective resistance at these sites, predictions and points

# do predictions of the effective resistance to ITNS
ingredient_weights <- readRDS("temporary/ingredient_weights.RDS")

ingredient_ids <- match(names(ingredient_weights), types)

# average the LLIN pyrethroids across these to plot against time
df_net_class_points <- df %>%
  filter(
    insecticide_type %in% pyrethroids
  ) %>%
  left_join(
    net_coverage_cell_lookup,
    by = "cell_id"
  ) %>%
  group_by(year_start,
           # insecticide_type,
           net_coverage_class) %>%
  summarise(
    died = sum(died),
    mosquito_number = sum(mosquito_number),
    .groups = "drop"
  ) %>%
  mutate(
    Susceptibility = died / mosquito_number,
  ) %>%
  rename(
    year = year_start
  )

# compute predicted population-level susceptibility to the insecticides
# average net coverage over time in each of those places

# now compute a weighted sum over the pyrethroids (weights given by the numbers
# of mosquitoes tested for each pyrethroid, over the whole timeseries)
pyrethroid_net_class_weights <- df %>%
  left_join(
    net_coverage_cell_lookup,
    by = "cell_id"
  ) %>%
  group_by(
    type_id,
    insecticide_type,
    net_coverage_class
  ) %>%
  summarise(
    mosquito_number = sum(mosquito_number),
    .groups = "drop"
  ) %>%
  group_by(net_coverage_class) %>%
  mutate(
    is_a_pyrethroid = insecticide_type %in% pyrethroids,
    mosquito_number = mosquito_number * as.numeric(is_a_pyrethroid),
    weight = mosquito_number / sum(mosquito_number)
  ) %>%
  ungroup() %>%
  select(
    type_id,
    weight,
    net_coverage_class
  ) %>%
  pivot_wider(
    names_from = net_coverage_class,
    values_from = "weight"
  ) %>%
  arrange(type_id) %>%
  select(-type_id) %>%
  as.matrix()

# AAARGH, this code is horrible, but best I can do for now on this plane.

# compute a weighted average of the insecticides in low, medium, high net
# coverage locations

# just subset the data and repeat what happens in the overall plot

low_idx <- net_coverage_cell_lookup %>%
  filter(
    net_coverage_class == "A) Low use"
  ) %>%
  pull(cell_id)

sample_size_mat_low <- sample_size_mat[low_idx, ]
weights_mat_low <- sweep(sample_size_mat_low,
                         2,
                         colSums(sample_size_mat_low),
                         "/")

# expand out into an array, replicating by year
weights_array_low <- array(rep(c(weights_mat_low), n_times),
                           dim = c(nrow(weights_mat_low), n_types, n_times))

# check dimensions are the right way around
# identical(weights_array_low[, , 1],  weights_mat_low)
# all(weights_array_low[1, 1, ] == weights_array_low[1, 1, 1])

# multiply through and sum to get average susceptibilities
weighted_susc_array_low <- dynamic_cells$all_states[low_idx, , ] * weights_array_low
overall_susc_low <- apply(weighted_susc_array_low, 2:3, "sum")
pyrethroid_net_class_susc_low <- as_data(t(pyrethroid_net_class_weights[, "A) Low use"])) %*% overall_susc_low

# and repeat for high
high_idx <- net_coverage_cell_lookup %>%
  filter(
    net_coverage_class == "B) High use"
  ) %>%
  pull(cell_id)

sample_size_mat_high <- sample_size_mat[high_idx, ]
weights_mat_high <- sweep(sample_size_mat_high,
                         2,
                         colSums(sample_size_mat_high),
                         "/")

# expand out into an array, replicating by year
weights_array_high <- array(rep(c(weights_mat_high), n_times),
                            dim = c(nrow(weights_mat_high), n_types, n_times))

# multiply through and sum to get average susceptibilities
weighted_susc_array_high <- dynamic_cells$all_states[high_idx, , ] * weights_array_high
overall_susc_high <- apply(weighted_susc_array_high, 2:3, "sum")
pyrethroid_net_class_susc_high <- as_data(t(pyrethroid_net_class_weights[, "B) High use"])) %*% overall_susc_high

# get posterior samples of these
sims_net_class <- calculate(pyrethroid_net_class_susc_low,
                            pyrethroid_net_class_susc_high,
                            values = draws,
                            nsim = 1e3)


pop_mort_sry_low <- sims_net_class$pyrethroid_net_class_susc_low[, 1, ] %>%
  t() %>%
  as_tibble() %>%
  mutate(
    time_id = row_number(),
    .before = everything()
  ) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_prefix = "V",
    names_to = "sim",
    values_to = "susceptibility"
  ) %>%
  group_by(time_id) %>%
  summarise(
    susc_pop_mean = mean(susceptibility),
    susc_pop_lower = quantile(susceptibility, 0.025),
    susc_pop_upper = quantile(susceptibility, 0.975),
  ) %>%
  mutate(
    year = baseline_year + time_id - 1
  )

pop_mort_sry_high <- sims_net_class$pyrethroid_net_class_susc_high[, 1, ] %>%
  t() %>%
  as_tibble() %>%
  mutate(
    time_id = row_number(),
    .before = everything()
  ) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_prefix = "V",
    names_to = "sim",
    values_to = "susceptibility"
  ) %>%
  group_by(time_id) %>%
  summarise(
    susc_pop_mean = mean(susceptibility),
    susc_pop_lower = quantile(susceptibility, 0.025),
    susc_pop_upper = quantile(susceptibility, 0.975),
  ) %>%
  mutate(
    year = baseline_year + time_id - 1
  )

pop_mort_sry_net_class <- bind_rows(
  `A) Low use` = pop_mort_sry_low,
  `B) High use` = pop_mort_sry_high,
  .id = "net_coverage_class"
)

net_class_fig <- df_net_class_points %>%
  ggplot(
    aes(
      x = year
    )
  ) +
  geom_ribbon(
    aes(
      x = year,
      ymax = upper,
      ymin = lower,
      group = net_coverage_class
    ),
    data = df_net_use,
    fill = grey(0.9),
    colour = grey(0.7)
  ) +
  geom_ribbon(
    aes(
      ymax = susc_pop_upper,
      ymin = susc_pop_lower,
    ),
    data = pop_mort_sry_net_class,
    fill = pyrethroid_blue
  ) +
  geom_point(
    aes(
      y = Susceptibility,
      size = mosquito_number
    ),
    shape = 21,
    fill = pyrethroid_blue
  ) +
  facet_wrap(~net_coverage_class) +
  guides(
    size = "none"
  ) +
  scale_y_continuous(
    labels = scales::percent,
    sec.axis = sec_axis(
      ~.,
      name = "LLIN use (grey)",
      labels = scales::percent,
    )
  ) +
  scale_size_continuous(
    labels = scales::number_format(accuracy = 1000, big.mark = ","),
    limits = size_limits
  ) +
  xlab("") +
  ylab("Susceptibility to pyrethroids") +
  coord_cartesian(xlim = c(2000, 2022),
                  ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    strip.text.x = element_text(hjust = 0)
  )

net_class_fig

ggsave("figures/all_africa_net_use_pred_data.png",
       bg = "white",
       scale = 0.8,
       width = 8,
       height = 4)

