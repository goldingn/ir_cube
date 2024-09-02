# plot estimated baseline susceptibility

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load the fitted model objects here, to set up predictions
load(file = "temporary/fitted_model.RData")

gadm_polys <- readRDS("data/clean/gadm_polys.RDS")
country_raster <- rast("data/clean/country_raster.tif")
region_raster <- rast("data/clean/region_raster.tif")

# expand out the country-level random effects with the unobserved countries
all_countries <- gadm_polys$country_name
n_countries_all <- length(all_countries)
unobserved_countries <- setdiff(all_countries, countries)
n_countries_unobserved <- length(unobserved_countries)

# combine into one random effect
init_country_raw_all <- zeros(n_countries_all, n_types)
init_country_raw_unobserved <- normal(0, 1, dim = c(n_countries_unobserved, n_types))

countries_observed_index <- match(countries, all_countries)
init_country_raw_all[countries_observed_index, ] <- init_country_raw

countries_unobserved_index <- match(unobserved_countries, all_countries)
init_country_raw_all[countries_unobserved_index, ] <- init_country_raw_unobserved

# combine to get the country-level deviation from the prior logit mean, for all
# countries
init_country_effect_all <- sweep(init_country_raw_all, 2, init_country_sd, FUN = "*")


country_region_index_all <- gadm_polys %>%
  as_tibble() %>%
  mutate(
    all_country_id = match(country_name, all_countries),
    region_id = match(region, regions)
  ) %>%
  group_by(all_country_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(all_country_id, region_id) %>%
  arrange(all_country_id) %>%
  pull(region_id)

# combine these together into the country-level logit-mean initial value
# logit_init_mean <- qlogis(init_frac_mean)
init_country_overall_effect_all <- init_country_effect_all + init_region_effect[country_region_index_all, ]
logit_init_country_all <- sweep(init_country_overall_effect_all,
                                2,
                                logit_init_mean,
                                FUN = "+")

init_country_relative_all <- ilogit(logit_init_country_all)
init_country_magnitude_all <- sweep(init_country_relative_all, 2, init_range, FUN = "*") 
init_country_all <- sweep(init_country_magnitude_all, 2, init_frac_min, FUN = "+")

# rearrange into a vector for summarising
index <- expand_grid(
  all_country_id = seq_len(n_countries_all),
  type_id = seq_len(n_types)
)

init_country_all_vec <- init_country_all[as.matrix(index)]

sim <- calculate(init_country_all_vec, values = draws, nsim = 1e3)

country_init_sims <- sim$init_country_all_vec[, , 1] %>%
  t() %>%
  as_tibble() %>%
  bind_cols(
    index,
    .
  ) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_prefix = "V",
    names_to = "sim",
    values_to = "susceptibility"
  ) %>%
  group_by(all_country_id, type_id) %>%
  summarise(
    susc_mean = mean(susceptibility),
    susc_lower = quantile(susceptibility, 0.025),
    susc_upper = quantile(susceptibility, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    country_name = all_countries[all_country_id],
    insecticide = types[type_id]
  )

init_polys <- gadm_polys %>%
  left_join(
    country_init_sims,
    by = join_by(country_name)
  )

ggplot() +
  geom_sf(
    aes(
      fill = susc_mean
    ),
    data = init_polys
  ) +
  scale_fill_gradient(
    labels = scales::percent,
    high = "palegreen",
    name = "Initial\nsusceptibility",
    limits = c(0.5, 1),
    na.value = "transparent") +
  facet_wrap(~insecticide, ncol = 3) +
  theme_ir_maps()

ggsave("figures/initial_susceptibility_map.png",
       bg = "white",
       width = 8,
       height = 8,
       dpi = 300)
