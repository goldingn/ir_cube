# Slices of the log posterior density of a saved full fit through its
# mortality floor(s) (#14, #37, #47).
#
#   Rscript R/floor_slice.R <fitted_model.RData> <label>
#
# For each floor mode of the fit's usable chains (all but those stuck), every
# other free parameter is held at the mean of that mode's chains' free states,
# and the log posterior density is evaluated on a fine grid of each floor (one
# for mortality_floor; with the species model's two floors, each in turn with
# the other at the mode's mean, and on a 41 x 41 grid of both). The density
# is that of the rebuilt model (rebuild_fit_model(), R/species_fit_helpers.R):
# "unadjusted", on the floor's own scale, which is plotted, and "adjusted",
# with the Jacobians of greta's transforms. A slice holds everything else
# fixed, so it shows the shape near each mode, not the marginal posterior of
# the floor (R/floor_profile.R).
#
# Writes outputs/species_runs/floor_slice/<label>.rds and
# figures/species_runs/floor_slice_<label>.png. Run with the greta 0.6
# environment (doc/cv_run_plan.md); FLOOR_PRIOR for fits saved without a
# floor prior (complete_model_options()). About 4 GB and 1-2 minutes.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2)
file <- arguments[1]
label <- arguments[2]

source("R/greta_setup.R")
start_greta(threads = 4)
suppressMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})
source("R/functions.R")
source("R/dynamical_predictions.R")
source("R/species_fit_helpers.R")

output_dir <- "outputs/species_runs/floor_slice"
figure_dir <- "figures/species_runs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

fit <- load_fit(file)
floors <- fit_floor_names(fit$draws)
if (length(floors) == 0) {
  report("%s has no floor; nothing to slice", label)
  quit(save = "no")
}
built <- rebuild_fit_model(fit)
log_density <- log_density_function(built$model)

usable <- usable_chains(fit$draws, label)
mode_of_chain <- chain_floor_mode(fit$draws)
modes <- split(usable, mode_of_chain[usable])
report("%s: floors %s; modes %s", label, toString(floors),
       paste(sprintf("%s (chains %s)", names(modes),
                     vapply(modes, toString, "")), collapse = ", "))

# the fine grid of a floor, dense near 0
floor_grid <- sort(unique(c(exp(seq(log(1e-4), log(0.01), length.out = 30)),
                            seq(0.01, 0.6, by = 0.005))))
coarse_grid <- c(1e-4, 5e-4, 0.001, 0.0025, seq(0.005, 0.6, length.out = 37))

slices <- list()
grids <- list()
for (mode in names(modes)) {
  centre <- colMeans(free_states(fit, built, modes[[mode]]))
  centre_floors <- floor_values(built$model$dag$trace_values(
    matrix(centre, nrow = 1))[, floors, drop = FALSE])[1, ]
  for (name in floors) {
    floor_free_check(built$model, centre, name)
    states <- matrix(centre, length(floor_grid), length(centre), byrow = TRUE)
    states[, free_column(built$model, name)] <- floor_free(floor_grid)
    time <- system.time(values <- log_density(states))[["elapsed"]]
    report("%s, %s: %d points in %.1f s", mode, name, nrow(states), time)
    slices[[length(slices) + 1]] <- tibble(
      mode = mode, floor = name, value = floor_grid,
      adjusted = values[, "adjusted"], unadjusted = values[, "unadjusted"],
      centre = centre_floors[[name]])
  }
  if (length(floors) == 2) {
    both <- expand.grid(first = coarse_grid, second = coarse_grid)
    states <- matrix(centre, nrow(both), length(centre), byrow = TRUE)
    states[, free_column(built$model, floors[1])] <- floor_free(both$first)
    states[, free_column(built$model, floors[2])] <- floor_free(both$second)
    time <- system.time(values <- log_density(states))[["elapsed"]]
    report("%s, both floors: %d points in %.1f s", mode, nrow(states), time)
    grids[[mode]] <- tibble(mode = mode, !!floors[1] := both$first,
                            !!floors[2] := both$second,
                            unadjusted = values[, "unadjusted"],
                            adjusted = values[, "adjusted"])
  }
}
slices <- bind_rows(slices) %>%
  group_by(mode, floor) %>%
  mutate(relative = unadjusted - max(unadjusted)) %>%
  ungroup()
grids <- bind_rows(grids)
saveRDS(list(label = label, file = file, floors = floors, modes = modes,
             slices = slices, grids = grids),
        file.path(output_dir, sprintf("%s.rds", label)))

# the maximum of each slice, and the density there relative to the other modes'
print(as.data.frame(slices %>%
                      group_by(mode, floor) %>%
                      slice_max(unadjusted, n = 1) %>%
                      select(mode, floor, centre, at_maximum = value,
                             unadjusted)), digits = 6)

p <- ggplot(slices, aes(value, pmax(relative, -50), colour = mode)) +
  geom_line() +
  geom_vline(aes(xintercept = centre, colour = mode), linetype = 2,
             data = distinct(slices, mode, floor, centre)) +
  facet_wrap(~ floor) +
  scale_x_sqrt() +
  labs(x = "floor (square-root scale); dashed: the mode's mean",
       y = "log posterior density, relative to each slice's maximum (cut at -50)",
       title = sprintf("%s: slices through the floor%s", label,
                       if (length(floors) > 1) "s" else ""),
       subtitle = "every other parameter at the mean of the mode's chains") +
  theme_bw(base_size = 9)
if (nrow(grids) > 0) {
  grid_plot <- grids %>%
    group_by(mode) %>%
    mutate(relative = pmax(unadjusted - max(unadjusted), -50)) %>%
    ungroup() %>%
    ggplot(aes(.data[[floors[1]]], .data[[floors[2]]], z = relative)) +
    geom_contour_filled(breaks = c(-50, -20, -10, -5, -2, -1, 0)) +
    facet_wrap(~ mode) +
    scale_x_sqrt() +
    scale_y_sqrt() +
    coord_equal() +
    labs(fill = "relative log\ndensity") +
    theme_bw(base_size = 9)
  p <- p / grid_plot
}
ggsave(file.path(figure_dir, sprintf("floor_slice_%s.png", label)), p,
       width = 8, height = if (nrow(grids) > 0) 8 else 4, dpi = 150)
report("saved; peak memory %.1f GB", peak_memory_gb())
