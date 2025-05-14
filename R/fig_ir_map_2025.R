# figure 2: 2025 maps of insecticide resistance 

# load packages and functions
source("R/packages.R")
source("R/functions.R")


# load in and prepare rasters

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


# standard plot

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


# patchwork plot to get scale into empty facet cell

# function to make plot of each insecticide per above
ir_plot_single <- function(
    r,
    name
){
  
  ggplot() +
    geom_spatraster(
      data = r
    ) +
    scale_fill_gradient(
      labels = scales::percent,
      name = "Susceptibility",
      limits = c(0, 1),
      na.value = "transparent") +
    labs(
      title = name
    ) +
    theme_ir_maps()

}

# check it
ir_plot_single(ir_rasts[[1]], name = insecticide_names[1])

# make list of ggplot objects for each insecticide
ir_map_list <- mapply(
  ir_plot_single,
  ir_rasts,
  insecticide_names,
  SIMPLIFY = FALSE
)


# arrange patchwork and plot
ir_map_list[[1]] +
  ir_map_list[[2]] +
  ir_map_list[[3]] +
  ir_map_list[[4]] +
  ir_map_list[[5]] +
  ir_map_list[[6]] +
  ir_map_list[[7]] +
  ir_map_list[[8]] +
  ir_map_list[[9]] +
  guide_area() + # this directs guide into empty cell in bottom right
  plot_layout(
    guides = "collect",
    ncol = 5
  )

# save
ggsave(
  "figures/ir_map_all_insecticides_2025_fig2.png",
  bg = "white",
  width = 18,
  height = 8,
  dpi = 300
)

# arrange patchwork and plot
ir_map_list[[1]] +
  ir_map_list[[2]] +
  ir_map_list[[3]] +
  ir_map_list[[4]] +
  guide_area() + # this directs guide into empty cell in bottom right
  ir_map_list[[5]] +
  ir_map_list[[6]] +
  ir_map_list[[7]] +
  ir_map_list[[8]] +
  ir_map_list[[9]] +
  plot_layout(
    guides = "collect",
    ncol = 5
  )

# save
ggsave(
  "figures/ir_map_all_insecticides_2025_fig2_alt.png",
  bg = "white",
  width = 18,
  height = 8,
  dpi = 300
)
