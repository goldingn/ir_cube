# plot covariates

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# plot time-varying covariates
itn <- rast("data/clean/net_use_cube.tif")
itn_plot_years <- c(2000, 2005, 2010, 2015, 2020)
itn_plot_layers <- paste0("nets_", itn_plot_years)
itn_plot <- itn[[itn_plot_layers]]
names(itn_plot) <- itn_plot_years

irs <- rast("data/clean/irs_coverage_cube.tif")
irs_plot_years <- c(2000, 2005, 2010, 2015, 2020)
irs_plot_layers <- paste0("irs_", irs_plot_years)
irs_plot <- irs[[irs_plot_layers]]
names(irs_plot) <- irs_plot_years

pop <- rast("data/clean/pop_cube.tif")
pop_plot_years <- c(2000, 2005, 2010, 2015, 2020)
pop_plot_layers <- paste0("pop_", pop_plot_years)
pop_plot <- pop[[pop_plot_layers]]
names(pop_plot) <- pop_plot_years

# get population density
pop_dens_plot <- pop_plot / terra::cellSize(pop_plot[[1]], unit = "km")

# plot ITN
max_nets <- max(global(itn_plot, "max", na.rm = TRUE))
ggplot() +
  geom_spatraster(
    data = itn_plot
  ) +
  scale_fill_gradient(
    name = "LLIN use",
    labels = scales::percent,
    limits = c(0, 1),
    low = grey(0.9),
    high = grey(0.2),
    na.value = "transparent") +
  facet_wrap(~lyr, nrow = 1) +
  ggtitle(
    label = "LLIN use"
  ) +
  theme_ir_maps()

ggsave("figures/itn_map.png",
       bg = "white",
       width = 18,
       height = 4,
       dpi = 300)

# plot IRS
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
  facet_wrap(~lyr, nrow = 1) +
  ggtitle(
    label = "IRS coverage"
  ) +
  theme_ir_maps()

ggsave("figures/irs_map.png",
       bg = "white",
       width = 18,
       height = 4,
       dpi = 300)


# plot pop
ggplot() +
  geom_spatraster(
    data = pop_dens_plot
  ) +
  scale_fill_continuous(
    trans = "log1p",
    labels = scales::number_format(big.mark = ","),
    low = grey(0.9),
    high = grey(0.2),
    breaks = c(1e1, 1e2, 1e3, 1e4),
    name = "People per km<sup>2</sup>",
    na.value = "transparent") +
  facet_wrap(~lyr, nrow = 1) +
  ggtitle(
    label = "Population density"
  ) +
  theme_ir_maps()

ggsave("figures/pop_map.png",
       bg = "white",
       width = 18,
       height = 4,
       dpi = 300)


# collated total yields of crop types
crops_group <- rast("data/clean/crop_group_scaled.tif")

# yields of individual crops
crops_all <- rast("data/clean/crop_scaled.tif")

# Pull out crop types implicated in risk for IR. refer to this review, crop type
# section:
# https://malariajournal.biomedcentral.com/articles/10.1186/s12936-016-1162-4
crops_implicated <- c(
  # "increased resistance at cotton growing sites, a finding subsequently
  # supported in eight other papers from five different African countries", "the
  # cash crop with the highest intensity insecticide use of any crop"
  crops_all$cotton,
  # "In eight studies, vegetable cultivation strongly related to
  # insecticide-resistant field collections", "Vegetable production requires
  # significantly higher quantities and/or more frequent application of
  # pesticides than other food crops"
  crops_all$vegetables,
  # "Seven of the studies reviewed here examined the insecticide susceptibility
  # of vector populations at rice-growing sites, and found low-to-moderate
  # resistance levels in these mosquito populations."
  crops_all$rice)

# combine all temporally-static covariates
crops <- c(crops_group, crops_implicated)
names(crops) <- str_to_sentence(names(crops))
ggplot() +
  geom_spatraster(
    data = crops
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
  ggtitle(
    label = "Agricultural intensities"
  ) +
  theme_ir_maps()

ggsave("figures/crop_map.png",
       bg = "white",
       width = 14,
       height = 6,
       dpi = 300)
