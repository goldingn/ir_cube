# visualise resistance data

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load data
ir_mtm_africa <- readRDS(file = "data/clean/mtm_data.RDS")

# visualise resistance data over time

# only plot An gambiae (s.l. and s.s.) and pyrethroids
ir_mtm_africa %>%
  # fix up some names for plotting
  rename(
    year = year_start,
    mortality = mortality_adjusted,
    country = country_name,
    insecticide = insecticide_type
  ) %>%
  filter(
    species %in% c("An. gambiae s.l.", "An. gambiae s.s."),
    insecticide_class %in% c("Pyrethroids")
  ) %>%
  # for each insecticide, filter to only the most common concentration
  group_by(
    insecticide
  ) %>%
  filter(
    insecticide_conc == sample_mode(insecticide_conc)
  ) %>%
  group_by(
    country
  ) %>%
  filter(
    n_distinct(year) >= 10,
    n() >= 100
  ) %>%
  ungroup() %>%
  ggplot(
    aes(
      x = year,
      y = mortality,
      group = country,
      colour = insecticide
    )
  ) +
  facet_wrap(~country,
             ncol = 3) +
  geom_point(alpha = 0.5) +
  theme_minimal()
