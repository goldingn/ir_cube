# Maps of the latent smooths of a saved full fit (V5, #47; R/latent_smooth.R):
# the posterior mean and sd of u_s(x), the log multiplier of the strength of
# selection, and u_f(x), the shift of the logit mortality floor, over the
# limits of transmission without water bodies
# (data/clean/pfpr_water_mask.tif, aggregated by 3, to about 14 km), from
# about n_draws (500) posterior draws, evenly spaced in each usable chain.
# With the shear (smooth_options(shear = TRUE)), u_s = v_s + b u_f, and v_s,
# the selection smooth's own part, is mapped too, with b in the caption. With
# the smooth of the initial state (smooth_options(init = TRUE)), u_init(x),
# at sd 1, is mapped too, and each type's loading lambda (its sd of the
# field: the type's logit relative initial state is shifted by lambda
# u_init) is in the caption.
#
#   Rscript R/smooth_maps.R <fitted_model.RData> <label>
#
# The means share one RdBu ramp, symmetric about 0: red where selection is
# weaker (u_s < 0, mapped as -u_s and labelled with the multiplier exp(u_s)),
# the floor higher (u_f > 0) or the initial state more susceptible (u_init >
# 0), blue where selection is stronger, the floor lower or the initial state
# more resistant. The sds share one sequential ramp, with the modelled
# bioassay cells as dots. Every smooth is 0 on average over the modelled
# cells. Writes
#   outputs/species_runs/smooth/<label>_smooths.tif  layers <smooth>_mean and
#                                                    <smooth>_sd
#   figures/species_runs/smooth_<label>.png
# Plain R; about 3 GB and 2 minutes.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2)
file <- arguments[1]
label <- arguments[2]
n_draws <- 500

suppressMessages({
  library(greta)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(terra)
  library(tidyterra)
  library(sf)
  library(ggtext)
})
source("R/functions.R")
source("R/dynamical_predictions.R")
source("R/species_fit_helpers.R")

output_dir <- "outputs/species_runs/smooth"
figure_dir <- "figures/species_runs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

fit <- load_fit(file)
if (!smooth_on(fit$options)) {
  report("%s has no latent smooths; nothing to map", label)
  quit(save = "no")
}
draws_used <- even_draws(fit, n_draws, label)
parameters <- fit_parameter_draws(fit, draws_used$index)
weights <- smooth_weight_terms(parameters$variables, parameters$options,
                               own = TRUE)
kinds <- intersect(c("selection", "selection_own", "floor", "init"),
                   names(weights))
shear <- if (isTRUE(fit$options$smooth$shear)) {
  c(parameters$variables$smooth_shear)
}
loadings <- parameters$init_loading

# the posterior mean and sd at the centre of each aggregated cell
grid <- terra::aggregate(rast("data/clean/pfpr_water_mask.tif"), 3,
                         fun = "max", na.rm = TRUE)
grid_cells <- terra::cells(grid)
xy <- terra::xyFromCell(grid, grid_cells)
time <- system.time(
  values <- smooth_posterior_at(parameters,
                                smooth_coords(xy[, 1], xy[, 2],
                                              fit$options$smooth$crs),
                                weights = weights)
)[["elapsed"]]
report("%s: %s at %d cells from %d draws in %.0f s", label, toString(kinds),
       length(grid_cells), parameters$n_draws, time)
layers <- list()
for (kind in kinds) {
  for (statistic in c("mean", "sd")) {
    layer <- rast(grid)
    layer[] <- NA_real_
    layer[grid_cells] <- values[[kind]][, statistic]
    names(layer) <- paste(kind, statistic, sep = "_")
    layers[[names(layer)]] <- layer
  }
}
layers <- rast(layers)
writeRaster(layers, file.path(output_dir, sprintf("%s_smooths.tif", label)),
            overwrite = TRUE, datatype = "FLT4S",
            gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
for (kind in kinds) {
  report("%s: posterior mean %.2f to %.2f, sd %.2f to %.2f", kind,
         min(values[[kind]][, "mean"]), max(values[[kind]][, "mean"]),
         min(values[[kind]][, "sd"]), max(values[[kind]][, "sd"]))
}


# the figure ---------------------------------------------------------------------

borders <- readRDS("data/clean/country_borders.RDS")
cells_xy <- as_tibble(terra::xyFromCell(grid, unique(terra::cellFromXY(
  grid, cbind(fit$df$longitude, fit$df$latitude)))))
# red where selection is weaker, the floor higher or the initial state more
# susceptible: -u_s, u_f and u_init
signed <- c(selection = -1, selection_own = -1, floor = 1, init = 1)
for (kind in kinds) {
  layers[[paste0(kind, "_mean")]] <- signed[[kind]] *
    layers[[paste0(kind, "_mean")]]
}
mean_limit <- max(vapply(kinds, function(kind) {
  max(abs(values[[kind]][, "mean"]))
}, numeric(1)))
sd_limit <- max(vapply(kinds, function(kind) max(values[[kind]][, "sd"]),
                       numeric(1)))
panel <- function(layer, fill_scale, title, dots = FALSE) {
  p <- ggplot() +
    geom_sf(data = borders, fill = grey(0.9), colour = NA) +
    geom_spatraster(data = layers[[layer]]) +
    fill_scale +
    geom_sf(data = borders, fill = NA, colour = grey(0.4), linewidth = 0.1)
  if (dots) {
    p <- p + geom_point(aes(x, y), data = cells_xy, size = 0.05,
                        colour = "black", alpha = 0.5)
  }
  p +
    coord_sf(xlim = c(-18, 52), ylim = c(-35, 25), expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_ir_maps()
}
mean_scale <- function(kind) {
  labels <- if (kind %in% c("floor", "init")) {
    waiver()
  } else {
    function(x) sprintf("%.2g", exp(-x))
  }
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
    limits = c(-mean_limit, mean_limit), oob = scales::squish,
    na.value = "transparent", labels = labels,
    name = switch(kind, selection = "multiplier\nexp(u_s)",
                  selection_own = "multiplier\nexp(v_s)",
                  floor = "logit floor\nshift u_f",
                  init = "initial state\nu_init (sd 1)"))
}
sd_scale <- scale_fill_gradientn(
  colours = RColorBrewer::brewer.pal(9, "YlGnBu"), limits = c(0, sd_limit),
  na.value = "transparent", name = "posterior\nsd")
titles <- c(selection = "selection, u_s: weaker (red) or stronger (blue)",
            selection_own = "selection's own part, v_s",
            floor = "floor, u_f: higher (red) or lower (blue)",
            init = "init, u_init: susceptible (red) or resistant (blue)")
panels <- list()
for (kind in kinds) {
  panels[[length(panels) + 1]] <- panel(paste0(kind, "_mean"),
                                        mean_scale(kind),
                                        sprintf("%s, posterior mean",
                                                titles[[kind]]))
}
for (kind in kinds) {
  panels[[length(panels) + 1]] <- panel(
    paste0(kind, "_sd"), sd_scale,
    sprintf("%s, posterior sd", sub("_own", ", own part", kind)), dots = TRUE)
}
shear_note <- if (!is.null(shear)) {
  sprintf("; u_s = v_s + b u_f, b %.2f (95%%: %.2f to %.2f)", mean(shear),
          quantile(shear, 0.025), quantile(shear, 0.975))
} else {
  ""
}
# the smooths' range: fixed, or the posterior median of each
range_note <- if (smooth_range_fixed(fit$options$smooth)) {
  sprintf("; range fixed at %s km",
          format(1000 * fit$options$smooth[["range"]], big.mark = ","))
} else {
  sprintf("; range (posterior median) %s",
          paste(sprintf("%s %s km", smooth_kinds(fit$options),
                        vapply(smooth_kinds(fit$options), function(kind) {
                          format(round(median(smooth_range_km(
                            parameters$variables, fit$options, kind))),
                            big.mark = ",")
                        }, "")), collapse = ", "))
}
floor_classes <- if (identical(fit$options$smooth$floor, "class")) {
  "; u_f applies to the pyrethroids and DDT only"
} else {
  ""
}
# each type's loading on u_init, the posterior mean
init_note <- if (!is.null(loadings)) {
  sprintf(paste0(";\nu_init shifts each type's logit relative initial ",
                 "state by lambda u_init, lambda (posterior mean) %s"),
          paste(sprintf("%s %.2f", fit$types, colMeans(loadings)),
                collapse = ", "))
} else {
  ""
}
p <- wrap_plots(panels, ncol = length(kinds), byrow = TRUE) +
  plot_annotation(
    title = sprintf("%s: the latent smooths", label),
    caption = sprintf(paste0(
      "posterior over %d draws; each smooth averages 0 over the modelled ",
      "bioassay cells (dots)%s%s%s%s;\none colour scale for the means, on ",
      "the log (selection) and logit (floor, initial state) scales"),
      parameters$n_draws, floor_classes, shear_note, range_note, init_note))
ggsave(file.path(figure_dir, sprintf("smooth_%s.png", label)), p,
       width = 5.2 * length(kinds) + 1, height = 8.2, dpi = 150, bg = "white")
report("saved; peak memory %.1f GB", peak_memory_gb())
