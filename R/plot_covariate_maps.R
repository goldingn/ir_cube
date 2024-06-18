# plot covariates

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# plot time-varying covariates
itn <- rast("data/clean/nets_per_capita_cube.tif")
itn_plot_years <- c(2000, 2010, 2020, 2022)
itn_plot_layers <- paste0("nets_", itn_plot_years)
itn_plot <- itn[[itn_plot_layers]]
names(itn_plot) <- itn_plot_years

irs <- rast("data/clean/irs_coverage_cube.tif")
irs_plot_years <- c(1997, 2010, 2020, 2022)
irs_plot_layers <- paste0("irs_", irs_plot_years)
irs_plot <- irs[[irs_plot_layers]]
names(irs_plot) <- irs_plot_years

# plot ITN
max_nets <- max(global(itn_plot, "max", na.rm = TRUE))
ggplot() +
  geom_spatraster(
    data = itn_plot
  ) +
  scale_fill_gradient(
    name = "Nets per capita",
    limits = c(0, max_nets),
    low = grey(0.9),
    high = grey(0.2),
    na.value = "transparent") +
  facet_wrap(~lyr, nrow = 2) +
  theme_minimal() +
  ggtitle(
    label = "LLIN coverage"
  )

ggsave("figures/itn_map.png",
       bg = "white",
       width = 8,
       height = 8,
       dpi = 300)

ggplot() +
  geom_spatraster(
    data = irs_plot
  ) +
  scale_fill_gradient(
    labels = scales::percent,
    low = grey(0.9),
    high = grey(0.2),
    name = "Coverage",
    limits = c(0, 1),
    na.value = "transparent") +
  facet_wrap(~lyr, nrow = 2) +
  theme_minimal() +
  ggtitle(
    label = "IRS coverage"
  )

ggsave("figures/irs_map.png",
       bg = "white",
       width = 8,
       height = 8,
       dpi = 300)


crops <- rast("data/clean/flat_covariates.tif")

ggplot() +
  geom_spatraster(
    data = other
  ) +
  scale_fill_gradient(
    # labels = scales::percent,
    transform = "sqrt",
    low = grey(0.9),
    high = grey(0.2),
    name = "Relative yield",
    limits = c(0, 1),
    na.value = "transparent") +
  facet_wrap(~lyr, nrow = 2) +
  theme_minimal() +
  ggtitle(
    label = "Agricultural intensities"
  )

ggsave("figures/crop_map.png",
       bg = "white",
       width = 12,
       height = 8,
       dpi = 300)
