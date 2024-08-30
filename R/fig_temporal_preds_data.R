# Make figure of model predictions over time for all of Africa overlaid with
# point-estimates of susceptibility bioassay data, aggregated up at various
# levels.

# Compute the point estimates as *weighted average* susceptibilities (weighted
# sum of number that died over weighted sum of number tested) to represent
# continent-wide sampling for that insecticide. Replace legend title with
# 'Effective sample size'

# Repeat the plots, with multiple lines and points for each regions. (Make
# functions to simplify this?)

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

# add on regions and tidy up column names
df_sub <- df %>%
  left_join(
    country_region_lookup(),
    by = "country_name"
  ) %>%
  rename(
    year = year_start,
    country = country_name,
    insecticide = insecticide_type
  )
  
# subset to pyrethroids, and add on UN geoscheme regions for Africa
df_pyrethroids <- df_sub %>%
  filter(
    insecticide_class == "Pyrethroids",
  )

# data aggregated over all of the pyrethroids and each region
df_pyrethroids_region_plot <- df_pyrethroids %>%
  group_by(year, region) %>%
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
average_sample_sizes <- df_sub %>%
  group_by(year, type_id) %>%
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
cell_sample_sizes <- df_sub %>%
  group_by(cell_id, year, type_id) %>%
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
pyrethroid_weights <- df_sub %>% 
  group_by(
    type_id,
    insecticide 
  ) %>%
  summarise(
    mosquito_number = sum(mosquito_number),
    .groups = "drop"
  ) %>%
  mutate(
    is_a_pyrethroid = insecticide %in% pyrethroids,
    mosquito_number = mosquito_number * as.numeric(is_a_pyrethroid),
    weight = mosquito_number / sum(mosquito_number)
  ) %>%
  arrange(type_id) %>%
  pull(weight)

pyrethroid_susc <- as_data(t(pyrethroid_weights)) %*% overall_susc


# get weighted average pyrethroid susceptibility for all cells and times
pyrethroid_weights_array <- array(pyrethroid_weights,
                                  dim = c(n_types, n_times, n_unique_cells))
pyrethroid_weights_array <- aperm(pyrethroid_weights_array, c(3, 1, 2))
pyrethroid_all_cells <- apply(pyrethroid_weights_array * dynamic_cells$all_states, c(1, 3), FUN = "sum")

# vector version of all cells and times for the weighted pyrethroid estimate
indices_all_pyrethroid <- expand_grid(
  cell_id = seq_len(n_unique_cells),
  time_id = seq_len(n_times)
)

pyrethroid_all_cells_vec <- pyrethroid_all_cells[as.matrix(indices_all_pyrethroid)]

# get posterior samples of these
sims <- calculate(overall_pop_mort_vec,
                  pyrethroid_susc,
                  pyrethroid_all_cells_vec,
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

# posterior means for all cells
pop_mort_sry_pyrethroid_all <- sims$pyrethroid_all_cells_vec[,  , 1] %>%
  t() %>%
  as_tibble() %>%
  bind_cols(
    indices_all_pyrethroid,
    .
  ) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_prefix = "V",
    names_to = "sim",
    values_to = "susceptibility"
  ) %>%
  group_by(cell_id, time_id) %>%
  summarise(
    susc_pop_mean = mean(susceptibility),
    .groups = "drop"
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

# now aggregate data for plotting Africa-wide averages
df_overall_plot <- df_sub %>%
  
  # first, group by cells, insecticides, and years
  group_by(cell_id,
           insecticide,
           year) %>%
  summarise(
    died = sum(died),
    mosquito_number = sum(mosquito_number),
    bioassays = n(),
    .groups = "drop"
  ) %>%
  
  # now compute the population fraction for each cell, for each insecticide, in
  # each year and in the full dataset, to compute weights reducing the annual
  # variability in the data estimates
  
  # for this insecticide, how many mosquitoes were collected in total
  group_by(insecticide) %>%
  mutate(
    total_overall_mosquito_number = sum(mosquito_number)
  ) %>%
  
  # for this insecticide, how many mosquitoes were collected in total in
  # this cell
  group_by(insecticide, cell_id) %>%
  mutate(
    cell_overall_mosquito_number = sum(mosquito_number)
  ) %>%
  
  # for this insecticide and this year, how many mosquitoes were collected in
  # total
  group_by(insecticide, year) %>%
  mutate(
    total_year_mosquito_number = sum(mosquito_number)
  ) %>%
  
  # for this insecticide and this year, how many mosquitoes were collected in
  # this cell
  group_by(insecticide, year, cell_id) %>%
  mutate(
    cell_year_mosquito_number = sum(mosquito_number)
  ) %>%
  
  # compute the cell's fraction of the total population in each year and
  # overall, and compute a corresponding weight
  ungroup() %>%
  mutate(
    year_fraction = cell_year_mosquito_number / total_year_mosquito_number,
    overall_fraction = cell_overall_mosquito_number / total_overall_mosquito_number,
    weight = overall_fraction / year_fraction
  ) %>%
  
  # now compute Africa-wide weighted susceptibilities for each year and insecticide
  group_by(
    year,
    insecticide
  ) %>%
  mutate(
    # components of the relative variance of the weights
    relvar_component = (weight - mean(weight)) ^ 2 / (mean(weight) ^ 2)
  ) %>%
  summarise(
    bioassays = sum(bioassays),
    died_weighted = sum(died * weight),
    mosquito_number_weighted = sum(mosquito_number * weight),
    # died = sum(died),
    mosquito_number = sum(mosquito_number),
    relvar = mean(relvar_component),
    .groups = "drop"
  ) %>%
  mutate(
    # Susceptibility_raw = died / mosquito_number,
    Susceptibility = died_weighted / mosquito_number_weighted,
    # design effect and effective sample size due to weighting
    design_effect = 1 + relvar,
    effective_bioassays = bioassays / design_effect,
    effective_samples = mosquito_number / design_effect,
    effective_died = Susceptibility * effective_samples,
    .after = insecticide
  ) %>%
  # add labels for plotting
  left_join(
    insecticides_plot_lookup,
    by = "insecticide"
  )

# now aggregate the pyrethroids, ensuring even weighting across pyrethroid types
# between years
df_pyrethroids_plot <- df_overall_plot %>%
  filter(
    insecticide %in% pyrethroids
  ) %>%

  # compute target weights for different pyrethroids over the whole dataset
  mutate(
    total_overall_samples = sum(effective_samples)
  ) %>%
  
  group_by(insecticide) %>%
  mutate(
    insecticide_overall_samples = sum(effective_samples)
  ) %>%

  group_by(year) %>%
  mutate(
    total_year_samples = sum(effective_samples)
  ) %>%
  
  group_by(insecticide, year) %>%
  mutate(
    insecticide_year_samples = sum(effective_samples)
  ) %>%
  
  ungroup() %>%
  mutate(
    year_fraction = insecticide_year_samples / total_year_samples,
    overall_fraction = insecticide_overall_samples / total_overall_samples,
    weight = overall_fraction / year_fraction
  ) %>%
  
  # collapse down to year
  group_by(year) %>%
  
  mutate(
    relvar_component = (weight - mean(weight)) ^ 2 / (mean(weight) ^ 2)
  ) %>%
  # now compute the year-specific population fraction for this cell in  insecticide
  summarise(
    bioassays_weighted = sum(effective_bioassays * weight),
    died_weighted = sum(effective_died * weight),
    samples_weighted = sum(effective_samples * weight),
    died = sum(effective_died),
    samples = sum(effective_samples),
    bioassays = sum(effective_bioassays),
    relvar = mean(relvar_component),
    .groups = "drop"
  ) %>%
  mutate(
    Susceptibility_raw = died / samples,
    Susceptibility = died_weighted / samples_weighted,
    # design effect and effective sample size due to weighting
    design_effect = 1 + relvar,
    effective_samples = samples / design_effect,
    effective_bioassays = bioassays / design_effect,
    .after = year
  )

# plot(df_pyrethroids_plot$Susceptibility_raw ~ df_pyrethroids_plot$Susceptibility)

# make them share a point size legend
max(c(df_pyrethroids_plot$effective_bioassays,
      df_overall_plot$effective_bioassays))
size_limits <- c(0, 500)

# add lines to prediction ribbon for visibility
line_size <- 0.5

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
    fill = pyrethroid_blue,
    colour = pyrethroid_blue,
    size = line_size,
    lineend = "round"
  ) +
  geom_point(
    aes(
      y = Susceptibility,
      size = effective_bioassays,
    ),
    shape = 21,
    fill = pyrethroid_blue,
    colour = "black",
    data = df_pyrethroids_plot
  ) +
  guides(
    size = "none"
  ) +
  facet_wrap(~ "A) All pyrethroids*") +
  scale_y_continuous(
    labels = scales::percent) +
  scale_size_area(
    limits = size_limits
  ) +
  xlab("") +
  ylab("Susceptibility") +
  coord_cartesian(xlim = c(1995, 2022),
                  ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    strip.text.x = element_text(hjust = 0)
  )

all_insecticides_fig <- pop_mort_sry %>%
  ggplot(
    aes(
      x = year,
      fill = insecticide_type_label,
      colour = insecticide_type_label,
    )
  ) +
  geom_ribbon(
    aes(
      ymax = susc_pop_upper,
      ymin = susc_pop_lower,
    ),
    size = line_size
  ) +
  geom_point(
    data = df_overall_plot,
    mapping = aes(
      y = Susceptibility,
      group = "none",
      size = effective_bioassays
    ),
    shape = 21,
    colour = "black"
  ) +
  facet_wrap(~insecticide_type_label,
             ncol = 3) +
  guides(
    size = guide_legend(title = "Effective\nsamples")
  ) +
  scale_y_continuous(
    labels = scales::percent,
    breaks = c(0, 1),
    limits = c(0, 1)) +
  scale_x_continuous(breaks = c(2000, 2020)) +
  scale_size_area(
    labels = scales::number_format(
      accuracy = 100,
      big.mark = ","
    ),
    breaks = c(1, 3, 5) * 100,
    limits = size_limits
  ) +
  scale_fill_discrete(
    direction = -1,
    guide = FALSE) +
  scale_colour_discrete(
    direction = -1,
    guide = FALSE) +
  xlab("") +
  # suppress ylab, as it is on the other panel
  ylab("") +
  coord_cartesian(xlim = c(1995, 2022)) +
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

