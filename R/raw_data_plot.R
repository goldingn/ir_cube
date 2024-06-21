# fit model

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load time-varying net coverage data and flatten it
nets_cube <- rast("data/clean/nets_per_capita_cube.tif")
years_sub <- 2005:2015
nets_cube_sub <- nets_cube[[paste0("nets_", years_sub)]]
nets_flat <- app(nets_cube_sub, "mean")

# load bioassay data
ir_mtm_africa <- readRDS(file = "data/clean/mtm_data.RDS")

# set the start of the timeseries considered
baseline_year <- 1990

insecticides_plot <- c(
  "Deltamethrin",
  "Permethrin",
  "Alpha-cypermethrin")

df <- ir_mtm_africa %>%
  filter(
    # subset to An. gambiae (s.l./s.s.)
    species %in% c("An. gambiae s.l.",
                   "An. gambiae s.s.",
                   "An. coluzzii",
                   "An. arabiensis"),
    # subset to WHO tube tests (most of data)
    test_type == "WHO tube test",
    # drop the minor classes: pyrroles and neonicotinoids
    insecticide_class %in% c("Pyrethroids",
                             "Carbamates",
                             "Organochlorines",
                             "Organophosphates"),
    # drop Dieldrin as it has two concentrations in the dataset, and those with
    # fewer than 200 observations (manually as I'm being lazy)
    !(insecticide_type %in% c("Dieldrin",
                              "Carbosulfan",
                              "Cyfluthrin",
                              "Etofenprox",
                              "Propoxur")),
    # drop any from before when we have data on net coverage
    year_start >= baseline_year
  ) %>%
  mutate(
    # convert concentrations into numeric values
    concentration = as.numeric(str_remove(insecticide_conc, "%")),
    net_coverage = terra::extract(nets_flat, select(., longitude, latitude))[, 2]
  ) %>%
  filter(
    # drop 4 records with missing values (coastal)
    !is.na(net_coverage)
  )

df_country_plot <- df %>%
  group_by(country_name, year_start, insecticide_type) %>%
  summarise(
    died = sum(died),
    mosquito_number = sum(mosquito_number),
    net_coverage = mean(net_coverage),
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
    net_coverage = mean(net_coverage),
    .groups = "drop"
  ) %>%
  mutate(
    Susceptibility = died / mosquito_number
  )

df_country_plot %>%
  ggplot(
    aes(
      x = year,
      y = Susceptibility,
      fill = insecticide
    )
  ) +
  geom_point(
    aes(
      group = country,
      size = mosquito_number
    ),
    colour = grey(0.9), 
  ) +
  geom_point(
    data = df_overall_plot,
    mapping = aes(
      group = "none",
      size = mosquito_number
    ),
    shape = 21,
  ) +
  guides(
    fill = guide_legend(title = "Insecticide",
                        order = 1,
                        override.aes = list(size = 4)),
    size = guide_legend(title = "No. tested")
  ) +
  scale_y_continuous(
    labels = scales::percent) +
  scale_size_continuous(
    labels = scales::number_format(accuracy = 1000, big.mark = ",")) +
  xlab("") +
  coord_cartesian(xlim = c(1998, 2022)) +
  theme_minimal() +
  ggtitle(
    "Growth of pyrethroid resistance across Africa",
    "Annual aggregate bioassay results across countries (grey) and all of Africa (colour)"
  )

ggsave("figures/bioassay_raw.png",
       bg = "white",
       scale = 0.8,
       width = 8,
       height = 6)

# Pool pyrethroids and plot by level of ITN use
df_country_plot %>%
  group_by(year, country) %>%
  summarise(
    died = sum(died),
    mosquito_number = sum(mosquito_number),
    net_coverage = mean(net_coverage),
    .groups = "drop"
  ) %>%
  mutate(
    Susceptibility = died / mosquito_number
  ) %>%
  ggplot(
    aes(
      x = year,
      y = Susceptibility,
    )
  ) +
  geom_point(
    aes(
      group = country,
      size = mosquito_number,
      colour = net_coverage
    ),
  ) +
  guides(
    colour = guide_legend(title = "LLIN coverage",
                        order = 1,
                        override.aes = list(size = 4)),
    size = guide_legend(title = "No. tested")
  ) +
  scale_colour_steps2(
    midpoint = median(df$net_coverage)
    # low = grey(0.8),
    # high = grey(0.2),
    # breaks = 
  ) +
  scale_y_continuous(
    labels = scales::percent) +
  scale_size_continuous(
    labels = scales::number_format(accuracy = 1000, big.mark = ",")) +
  xlab("") +
  coord_cartesian(xlim = c(1998, 2022)) +
  theme_minimal() +
  ggtitle(
    "Growth of pyrethroid resistance across Africa",
    "Annual aggregate bioassay results across countries (grey) and all of Africa (colour)"
  )



# for model validation, group data by predicted value, and pool all the data
# within those quantiles?