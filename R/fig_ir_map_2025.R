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

pyrethroid_names <- c("Alpha-cypermethrin",
                      "Deltamethrin",
                      "Lambda-cyhalothrin",
                      "Permethrin")
pyrethroid_star <- ifelse(insecticide_names %in% pyrethroid_names, "*", "")

# make figure names, with panel letters and asterisks
insecticide_figure_names <- paste0(
  LETTERS[2:10],
  ") ",
  insecticide_names,
  pyrethroid_star)

names(ir_rasts) <- insecticide_figure_names

# change the coordinates for plotting all rasters, to ensure the plots are approximately square
new_ratio <- 0.8

ext <- as.vector(ext(llin_rast))
diffs <- diff(ext)[c(1, 3)]

# keep the y limits the same, and adjust the x limits
xmid <- ext[1] + 0.5 * diffs[1]
xdiff <- diffs[2] / new_ratio

# rebuild the limits
xlim <- xmid + xdiff * c(-0.5, 0.5)


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
    name = "Susceptibility"
    # guide = "none"
  ) +
  # coord_sf(xlim = xlim) +
  theme_ir_maps() +
  theme(
    strip.text.x = element_text(hjust = 0)
  )

# make 9 different plots, for each of the individual insecticides

n_insecticides <- length(insecticide_figure_names)
all_insecticides_fig_list <- list()
all_insecticides_cols <- rev(scales::hue_pal()(n_insecticides))

for (i in seq_len(n_insecticides)) {
  all_insecticides_fig_list[[i]] <- ggplot() +
    geom_spatraster(
      data = ir_rasts[[i]]
    ) +
    facet_wrap(~lyr,
               ncol = 3) +
    scale_fill_gradient(
      labels = scales::percent,
      name = "Susceptibility",
      high = all_insecticides_cols[i],
      guide = "none",
      limits = c(0, 1),
      na.value = "transparent"
    ) +
    coord_sf(xlim = xlim) +
    theme_ir_maps() +
    theme(
      strip.text.x = element_text(hjust = 0)
    )
}

# make a 9-panel plot of susceptibility to all insecticides in 2025
all_insecticides_fig <- patchwork::wrap_plots(all_insecticides_fig_list,
                                              ncol = 3) +
  plot_layout()

# combine them to plot
combined_fig <- pyrethroid_fig + all_insecticides_fig

# save the plot
ggsave(
  filename = "figures/ir_map_all_insecticides_2025.png",
  plot = combined_fig,
  bg = "white",
  width = 14,
  height = 6,
  scale = 0.8
)

# to do:
# - add a legend on the side (wrap_plots(guides = "collect") ?)