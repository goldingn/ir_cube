# figure 2: 2025 maps of insecticide resistance 

# load packages and functions
source("R/packages.R")
source("R/functions.R")

insecticide_names <- c(
  "Deltamethrin",
  "Alpha-cypermethrin",
  "Pirimiphos-methyl", 
  "Permethrin",
  "DDT",
  "Bendiocarb",
  "Malathion",
  "Fenitrothion", 
  "Lambda-cyhalothrin"#,
  #"llin_effective"
)

ir_yr <- 2025

ir_filenames <- sprintf(
  "outputs/ir_maps/%s/ir_%s_susceptibility.tif",
  insecticide_names,
  ir_yr
)


ir_rasts <- rast(ir_filenames)

names(ir_rasts) <- insecticide_names


ggplot() +
  geom_spatraster(
    data = ir_rasts
  ) +
  scale_fill_gradient(
    labels = scales::percent,
    name = "Susceptibility",
    limits = c(0, 1),
    na.value = "transparent") +
  facet_wrap(~lyr, nrow = 2, ncol = 5) +
  labs(
    #subtitle = "Susceptibility of An. gambiae (s.l./s.s.) in WHO bioassays in 2025"
    subtitle = expression(Susceptibility~of~italic("An. gambiae")~`(s.l./s.s.)`~`in`~WHO~bioassays~`in`~2025)
  ) +
  theme_ir_maps()


ggsave(
  "figures/ir_map_all_insecticides_2025.png",
  bg = "white",
  width = 18,
  height = 8,
  dpi = 300
)
