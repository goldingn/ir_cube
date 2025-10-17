# figure 2: 2025 maps of insecticide resistance 

# load packages and functions
source("R/packages.R")
source("R/functions.R")


ir_yr <- 2025

# load in and prepare rasters

# First, the pyrethroids used in LLINs:

llin_filename <- sprintf(
  "outputs/ir_maps/llin_effective/ir_%s_susceptibility.tif",
  ir_yr
)

llin_rast <- rast(llin_filename)
names(llin_rast) <- "A) LLIN pyrethroids"


# now, all the individual ones
insecticide_names <- c(
  "Alpha-cypermethrin",
  "Deltamethrin",
  "Lambda-cyhalothrin",
  "Permethrin",
  "Fenitrothion", 
  "Malathion",
  "Pirimiphos-methyl", 
  "DDT",
  "Bendiocarb"
)

ir_filenames <- sprintf(
  "outputs/ir_maps/%s/ir_%s_susceptibility.tif",
  insecticide_names,
  ir_yr
)

ir_rasts <- rast(ir_filenames)

insecticide_figure_names <- paste0(
  LETTERS[2:10],
  ") ",
  insecticide_names)

names(ir_rasts) <- insecticide_figure_names


# make a plot of susceptibility to Pyrethroid insecticides in 2025

pyrethroid_fig <- ggplot() +
  geom_spatraster(
    data = llin_rast
  ) +
  facet_wrap(~lyr) +
  scale_fill_gradient(
    labels = scales::percent,
    limits = c(0, 1),
    na.value = "transparent",
    guide = "none"
  ) +
  theme_ir_maps() +
  theme(
    strip.text.x = element_text(hjust = 0)
  )

# make a 9-panel plot of susceptibility to all insecticides in 2025
all_insecticides_fig <- ggplot() +
  geom_spatraster(
    data = ir_rasts
  ) +
  facet_wrap(~lyr,
             ncol = 3) +
  guides(
    size = guide_legend(title = "Susceptibility")
  ) +
  scale_fill_gradient(
    labels = scales::percent,
    name = "Susceptibility",
    limits = c(0, 1),
    na.value = "transparent"
  ) +
  theme_ir_maps() +
  theme(
    strip.text.x = element_text(hjust = 0)
  )
all_insecticides_fig


# combine them to plot
pyrethroid_fig + all_insecticides_fig

# save the plot
ggsave(
  "figures/ir_map_all_insecticides_2025.png",
  bg = "white",
  width = 14,
  height = 6,
  scale = 0.8
)

# need to make colour ramps in all_insecticides_fig match the panels in figure 1

