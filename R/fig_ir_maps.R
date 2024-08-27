# plot IR predictions

# load packages and functions
source("R/packages.R")
source("R/functions.R")

this_insecticide <- "Deltamethrin"
plot_years <- c(2000, 2005, 2010, 2015, 2020, 2025, 2030)

directory <- file.path("outputs/ir_maps", this_insecticide)
files <- file.path(directory,
                   sprintf("ir_%s_susceptibility.tif",
                           plot_years))
this_raster <- rast(files)
names(this_raster) <- plot_years

ggplot() +
  geom_spatraster(
    data = this_raster
  ) +
  scale_fill_gradient(
    labels = scales::percent,
    name = "Susceptibility",
    limits = c(0, 1),
    na.value = "transparent") +
  facet_wrap(~lyr, nrow = 2, ncol = 5) +
  ggtitle(
    label = this_insecticide,
    subtitle = "Susceptibility of An. gambiae (s.l./s.s.) in WHO bioassays"
  ) +
  theme_ir_maps()

figure_name <- file.path(
  "figures",
  sprintf("%s_ir_map.png",
          this_insecticide)
)

ggsave(figure_name,
       bg = "white",
       width = 18,
       height = 8,
       dpi = 300)
