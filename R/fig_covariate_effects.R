# make a figure showing the fitted relationships between covariates and
# selection for each insecticide

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load the fitted model objects here, to set up predictions
load(file = "temporary/fitted_model.RData")

# get posterior means and standard deviations of these
sims <- calculate(effect_type, values = draws, nsim = 1e4)

# make into tibbles for plotting
mean_tibble <- sims$effect_type %>%
  apply(2:3, mean) %>%
  `colnames<-`(types) %>%
  as_tibble() %>%
  mutate(covariate = colnames(x_cell_years)) %>%
  pivot_longer(
    cols = -covariate,
    names_to = "insecticide",
    values_to = "mean"
  )

sd_tibble <- sims$effect_type %>%
  apply(2:3, sd) %>%
  `colnames<-`(types) %>%
  as_tibble() %>%
  mutate(covariate = colnames(x_cell_years)) %>%
  pivot_longer(
    cols = -covariate,
    names_to = "insecticide",
    values_to = "sd"
  )

# get the mean of the scaled covariates in the data
data_means <- x_cell_years %>%
  as_tibble() %>%
  pivot_longer(
    cols = everything(),
    names_to = "covariate",
    values_to = "value"
  ) %>%
  group_by(covariate) %>%
  summarise(
    value = mean(value),
    .groups = "drop"
  )

stats <- mean_tibble %>%
  left_join(
    sd_tibble,
    by = c("covariate", "insecticide")
  ) %>%
  left_join(
    data_means,
    by = "covariate"
  ) %>%
  mutate(
    mean_scaled = mean * value,
    # is this right foraffine transform of a normal? Can't check as I'm on a
    # plane
    sd_scaled = sd * value,
    precision = 1 / sd ^ 2,
    precision_scaled = 1 / sd_scaled ^ 2,
    inv_coef_variation = mean / sd ^ 2,
    inv_coef_variation_scaled = mean_scaled / sd_scaled ^ 2
  ) 

stats %>%
  ggplot(
    aes(
      x = insecticide,
      y = covariate,
      fill = mean_scaled,
      # alpha = precision_scaled
    )
  ) +
  geom_tile() +
  scale_fill_gradient(
    low = "white",
    high = "red"
  ) +
  scale_x_discrete(
    position = "top"
  ) +
  scale_alpha_continuous(
    range = c(0.2, 1)
  ) +
  theme_minimal() +
  xlab("") +
  guides(
    fill = guide_legend(
      title = "Average selection effect"
    )
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 0)
  )

# The average contribution of each covariate to the selection coefficient for
# resistance to each insecticide. Computed as the posteiror mean of the partial
# selection coefficient, multiplied by the mean over the dataset of the
# indicator variable for the covariate. This captures both the relative strength
# of the covariate in selecting for resistance, and its prevalence.

ggsave(
  "figures/covariate_loadings.png",
  bg = "white",
  scale = 0.8,
  width = 8,
  height = 5
)  

