# illustrate the predictive distribution validation metric

source("R/packages.R")
source("R/functions.R")

# load the fitted model objects here, to set up predictions
load(file = "temporary/fitted_model.RData")

# load the mask
mask <- rast("data/clean/raster_mask.tif")

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

insecticides_plot_small <- c("Deltamethrin",
                             "Permethrin",
                             "Alpha-cypermethrin")

# subset df to the cell and year with the most bioassay results
df_most_sampled <- df %>%
  filter(
    insecticide_type %in% insecticides_plot_small
  ) %>%
  group_by(
    cell,
    cell_id,
    insecticide_type,
    year_start
  ) %>%
  summarise(
    count = n(),
    .groups = "drop"
  ) %>%
  arrange(
    desc(count)
  ) %>%
  slice(2) %>%
  bind_cols(
    xyFromCell(mask, .$cell)
  ) %>%
  reverse_geocode(
    lat = y,
    long = x,
    method = 'osm',
    full_results = TRUE
  ) %>%
  mutate(
    precise_place = str_split_i(address, ",", 1),
    precise_place = case_when(
      grepl("Placodji", address) ~ "Placodji, Cotonou",
      grepl("Busia", address) ~ "Busia",
      .default = precise_place
    ),
    place = paste(precise_place, country, sep = ", ")
  ) %>%
  select(
    cell,
    year_start,
    insecticide_type,
    count,
    place
  ) %>%
  left_join(
    mutate(df,
           index = row_number()),
    by = c("cell", "year_start", "insecticide_type")
  ) %>%
  mutate(
    Susceptibility = died / mosquito_number,
    lower = binom::binom.bayes(died, mosquito_number)$lower,
    upper = binom::binom.bayes(died, mosquito_number)$upper
  )

df_most_sampled_overall <- df_most_sampled %>%
  summarise(
    mean = mean(Susceptibility),
    lower_50 = quantile(Susceptibility, 0.25),
    upper_50 = quantile(Susceptibility, 0.75),
    lower_60 = quantile(Susceptibility, 0.2),
    upper_60 = quantile(Susceptibility, 0.8),
    lower_70 = quantile(Susceptibility, 0.15),
    upper_70 = quantile(Susceptibility, 0.85),
    lower_80 = quantile(Susceptibility, 0.1),
    upper_80 = quantile(Susceptibility, 0.9),
    lower_90 = quantile(Susceptibility, 0.05),
    upper_90 = quantile(Susceptibility, 0.95),
    lower_100 = quantile(Susceptibility, 0.00001),
    upper_100 = quantile(Susceptibility, 0.9999)
  )

# colours for different insecticides
colour_types <- scales::hue_pal(direction = -1)(8)
insecticides_plot <- tibble(
  insecticide = types,
  class = classes[classes_index]
) %>%
  arrange(desc(class), insecticide) %>%
  filter(
    !(insecticide %in% c("DDT")) 
  ) %>%
  pull(insecticide)
types_plot_id <- match("Deltamethrin", insecticides_plot)
colours_plot <- colour_types[types_plot_id]

set.seed(1)
df_most_sampled %>%
  # shuffle the x axes so the points don't overlap too much
  mutate(
    x_random = sample(seq(0.2, 1.8, length.out = n()))
  ) %>%
  ggplot(
    aes(
      xmin = 0,
      xmax = 2
    )
  ) +
  geom_rect(
    aes(
      ymax = upper_100,
      ymin = lower_100,
    ),
    data = df_most_sampled_overall,
    fill = colours_plot,
    alpha = 0.15
  ) +
  geom_rect(
    aes(
      ymax = upper_90,
      ymin = lower_90,
    ),
    data = df_most_sampled_overall,
    fill = colours_plot,
    alpha = 0.15
  ) +
  geom_rect(
    aes(
      ymax = upper_80,
      ymin = lower_80,
    ),
    data = df_most_sampled_overall,
    fill = colours_plot,
    alpha = 0.15
  ) +
  geom_rect(
    aes(
      ymax = upper_70,
      ymin = lower_70,
    ),
    data = df_most_sampled_overall,
    fill = colours_plot,
    alpha = 0.15
  ) +
  geom_rect(
    aes(
      ymax = upper_60,
      ymin = lower_60,
    ),
    data = df_most_sampled_overall,
    fill = colours_plot,
    alpha = 0.15
  ) +
  geom_rect(
    aes(
      ymax = upper_50,
      ymin = lower_50,
    ),
    data = df_most_sampled_overall,
    fill = colours_plot,
    alpha = 0.15
  ) +
  geom_hline(
    yintercept = df_most_sampled_overall$mean,
    linewidth = 1,
    colour = grey(0.4)
  ) +
  geom_point(
    aes(
      y = Susceptibility,
      size = mosquito_number,
      x = x_random
    ),
    fill = colours_plot,
    shape = 21
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  scale_x_continuous(
    limits = c(0, 2)
  ) +
  xlab("") +
  ggtitle(
    sprintf("%s bioassays",
            df_most_sampled$insecticide_type[1]),
    sprintf("%s, %s",
            df_most_sampled$place[1],
            df_most_sampled$year_start)
  ) +
  theme_minimal() +
  guides(
    size = guide_legend(title = "No. tested")
  ) +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) 

ggsave(
  "figures/bioassay_variability.png",
  scale = 0.9,
  width = 4,
  height = 5,
  bg = "white"
)

df_most_sampled$count[1]
# 11x discriminating concentration bioassay results against Deltamethrin at a
# single location in 2014, as an example of the inherent variability in this
# type of data. Horizontal line gives the mean of the susceptibility estimates
# and coloured bands give quantiles of the distribution over estimates. Note
# that for this illustration neither the mean or quantiles are weighted by
# sample size.

# get ECDF for each of the 11 data points, and plot the observation against each
# one (do as susceptibiltiy fraction) and report quantile for each one.

# plot histogram of RQRs for this data point, and for a larger sample (all
# Deltamethrin bioassays in Kenya), and QQ plots?

# Display the GoF statistic (Anderson Darling) for each of these histograms,
# versus a uniform distribution

