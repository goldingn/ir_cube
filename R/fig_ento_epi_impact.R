# Plot the decline (2000-2025 and projected 2025-2030) in the effectiveness of
# ITNs across Africa due to IR. Plot the decline in both direct entomological
# protection, and in epidemiological outcomes, based on Imperial and MAP models.

# A 2-panel plot, x-axis: time, y-axis: relative impact, one each for ento and
# epi outcomes for 2000-2030, with dashed projected lines from the latest ITN
# coverages to 2030. In each panel, one line shows a counterfactual where nets
# retained their full impact (mostly increasing over time due to increasing
# coverage), the other showing the reduced impact due to IR.  0%-100%. Relative
# impact is relative to a hypothetical 100% impact with 100% net use and 0%
# resistance. Spatial average by the infected population (MAP receptive PR times
# annual population).

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# data to load:

# spatial:
# - Country borders for plotting
# - IR layers
# - MAP net use (to manually project to 2030)
# - MAP PR 2000 baseline
# - MAP populations by year (projected to 2030)
# other:
# - Nash/Sherrard-Smith killing effect estimates
# - Symons epi impact estimate

# years to plot for
years_plot <- 2000:2030

# spatial:

# - Country borders for plotting
borders <- readRDS("data/clean/gadm_polys.RDS")

# - IR layers (susceptibility to the pyrethroids used in LLINs)
ir_filenames <- sprintf(
  "outputs/ir_maps/llin_effective/ir_%s_susceptibility.tif",
  years_plot
)
ir <- rast(ir_filenames)
names(ir) <- years_plot

# - MAP net use (to manually project to 2030)
nets <- rast("data/clean/net_use_cube.tif")
names(nets) <- str_remove(names(nets), "^nets_")
latest_nets <- max(as.numeric(names(nets)))
end_year <- max(years_plot)
nets <- post_pad_cube(nets, end_year)


# - MAP PR 2000 baseline
pr2000 <- rast("data/clean/pfpr_2000.tif")

# - MAP populations by year (projected to 2030)
pop <- rast("data/clean/pop_cube.tif")
names(pop) <- str_remove(names(pop), "^pop_")
pop_future <- rast("data/clean/pop_cube_future.tif")
names(pop_future) <- str_remove(names(pop_future), "^pop_")
pop_all <- c(pop, pop_future)

# other:

# - Nash/Sherrard-Smith killing effect estimates.

# load in and prepare curves from Nash/Sherard-Smith, given in supp to Sherard
# Smith paper: https://doi.org/10.1038/s41467-022-30700-1
ss_est <- read_excel("data/raw/41467_2022_30700_MOESM3_ESM.xlsx",
                        sheet = 4,
                        range = "A15:Y117") %>%
  select(
    resistance = "Pyrethroid resistance approximated from susceptibility bioassay",
    kill = "dN0_4"
  ) %>%
  slice(
    -1
  ) %>%
  mutate(
    susceptibility = 1 - resistance,
    relative_killing_effect = kill / max(kill)
  )
# 
# # plot to check against paper figures
# ss_est %>%
#   ggplot(
#     aes(
#       x = resistance,
#       y = kill
#     )
#   ) +
#   geom_line() +
#   theme_minimal()

# produce an interpolating function to get the relative killing effect
killing_function <- splinefun(x = ss_est$susceptibility,
                              y = ss_est$relative_killing_effect)

# check this interpolator behaves appropriately at the edges and throughout
stopifnot(identical(killing_function(0), 0))
stopifnot(identical(killing_function(1), 1))
stopifnot(identical(killing_function(ss_est$susceptibility),
                    ss_est$relative_killing_effect))

# encode Symons et al. PMI function for relative effectiveness of ITNs at
# reducing logit-PfPR (reduction in on the log-odds ratio), compared to in the absence
# of resistance.

# when susceptibility is 100%, the ITN effectiveness parameter has full effect
# (multiplier 1), when susceptibility is 0%, the ITN effectiveness parameter has
# 54% (multiplier 0.54).

epi_function <- function(susceptibility) {
  1 - 0.46 * (1 - susceptibility)
}

# The log-odds ratio for malaria infection under full ITN coverage (relative to
# no ITN coverage) is -1.68, which corresponds to an odds ratio of 0.186374.
# With a fully-resistant vector population, this is would reduce to 54% of its
# previous value, -0.9072, giving an odds ratio of 0.4036529.

# Alternative interpretation: the odds ratio for a population with 100% coverage
# and 100% resistance is equivalent to that of a population with 54% coverage
# and 0% resistance.

# perform calculations:


# make a cube of the killing effect of ITNs, over time and space
relative_killing <- terra::app(ir, killing_function)
names(relative_killing) <- names(ir)

# make a cube of the epi impact of ITNs, over time and space
relative_epi_impact <- terra::app(ir, epi_function)
names(relative_epi_impact) <- names(ir)

# plot(relative_epi_impact[[c("2000", "2010", "2020", "2030")]],
#      range = c(0, 1))

# first, average the killing effects and the epi impact of the nets in use
# across Africa to get the temporal trends in impact of each net

# make a cube of net weights
normalise <- function(x) {
  total <- global(x, "sum", na.rm = TRUE)[1, 1]
  x / total
}
net_weights <- sapp(nets, normalise)

# average the killing effect and epi impact of nets over the nets in use
average_killing <- global(relative_killing * net_weights,
                          fun = "sum",
                          na.rm = TRUE)

average_epi_impact <- global(relative_epi_impact * net_weights,
                          fun = "sum",
                          na.rm = TRUE)


# plot the relative impact of nets on the whole continent by both measures, with
# and without IR

# weight by potential infected population
infection_weights <- sapp(pop_all * pr2000, normalise)

# average net use over these locations
net_impact_no_ir <- global(nets * infection_weights,
                           fun = "sum",
                           na.rm = TRUE)
# apply the impact of IR and recalculate average
net_impact_killing <- global(nets * relative_killing * infection_weights,
                             fun = "sum",
                             na.rm = TRUE)

net_impact_epi_impact <- global(nets * relative_epi_impact * infection_weights,
                             fun = "sum",
                             na.rm = TRUE)

par(mfrow = c(1, 3),
    mar = c(2, 5, 3, 2))

# Effectiveness\n of single-AI nets
no_change <- rep(1, length(years_plot))
# impact of each net
plot(no_change ~ years_plot,
     type = "l",
     col = grey(0.3),
     lwd = 2,
     xlab = "",
     ylab = "Relative effectiveness",
     ylim = c(0, 1))
lines(average_killing$sum ~ years_plot,
     col = "red",
     lwd = 2)
abline(v = latest_nets, lty = 2)
title(main = "Vector killing")

plot(no_change ~ years_plot,
     type = "l",
     col = grey(0.3),
     lwd = 2,
     xlab = "",
     ylab = "Relative effectiveness",
     ylim = c(0, 1))
lines(average_epi_impact$sum ~ years_plot,
      col = "blue",
      lwd = 2)
abline(v = latest_nets, lty = 2)
title(main = "Transmission\nreduction")


plot(net_impact_no_ir$sum ~ years_plot,
     type = "l",
     lwd = 2,
     col = grey(0.3),
     ylab = "Relative impact",
     xlab = "",
     ylim = c(0, 1))
lines(net_impact_epi_impact$sum ~ years_plot,
      col = "blue",
      lwd = 2)
abline(v = latest_nets, lty = 2)
title(main = "Impact of\nnet distribution")

# redo in ggplot and make future lines dotted

# do summaries for a few specific countries, or regions of countries, to make
# the point about stratification (make sure the raw data aligns first!)

# add limits of transmission mask to other maps

