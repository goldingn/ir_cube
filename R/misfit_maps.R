# Maps of a fit's pyrethroid misfit (#37): 2D smooths over longitude and
# latitude, weighted by precision and predicted to the mask (blank more than
# 2.5 degrees from any data), of three quantities, each signed so that
# positive means the fit predicts too much resistance:
#   1. the level residual: empirical minus the fit's predicted logit of
#      mortality of each bioassay (positive: observed mortality higher);
#   2. the change misfit of within-pixel pairs: predicted minus observed fall
#      in logit mortality per year (positive: the fit over-predicts
#      selection; the negative of R/change_misfit_plot.R's);
#   3. the slope of the level residual on cumulative net use, varying in
#      space (residual ~ s(lon, lat) + s(lon, lat, by = cumulative nets)):
#      positive where the fit over-predicts resistance more the more nets
#      have accumulated.
#
#   Rscript R/misfit_maps.R <fit RData> <label> [gap] [logits rds] [name]
#
# With a logits file (from R/net_w_screen.R, a list of each bioassay's fitted
# logit by w) and the name of one of its elements, the predictions are those
# logits rather than the fit's, and the change misfit uses the predicted
# change in mean logit between each pair's pixel-years, not the simulated
# empirical-logit change.
#
# Writes figures/misfit_maps_<label>[_<name>].png.

arguments <- commandArgs(trailingOnly = TRUE)
fit_file <- arguments[1]
label <- arguments[2]
gap <- if (length(arguments) >= 3) as.integer(arguments[3]) else 5L
logits_file <- if (length(arguments) >= 4) arguments[4] else NA
logits_name <- if (length(arguments) >= 5) arguments[5] else NA

source("R/packages.R")
source("R/functions.R")
source("R/dynamical_model.R")
source("R/dynamical_predictions.R")
suppressMessages({
  library(mgcv)
  library(patchwork)
})

f <- new.env()
load(fit_file, envir = f)
df <- f$df
mask <- rast("data/clean/raster_mask.tif")
fold <- list(draws = f$draws, options = f$model_options,
             x_cells_init = f$x_cells_init)
n_total <- sum(vapply(fold$draws, nrow, integer(1)))
par <- dynamical_parameter_draws(fold, f$classes_index, f$types,
                                 draw_index = round(seq(1, n_total,
                                                        length.out = 200)))
predicted <- if (is.na(logits_file)) {
  colMeans(dynamical_logit(par, select(df, cell_id, country_id, type_id,
                                       year_id),
                           df, f$x_cell_years, f$cell_years_index))
} else {
  readRDS(logits_file)[[logits_name]]
}
stopifnot(length(predicted) == nrow(df))
rho <- colMeans(par$rho_types)[df$type_id]
n_years <- max(f$cell_years_index$year_id)
x_row <- matrix(NA_integer_, max(f$cell_years_index$cell_id), n_years)
x_row[cbind(f$cell_years_index$cell_id, f$cell_years_index$year_id)] <-
  seq_len(nrow(f$cell_years_index))
nets <- matrix(f$x_cell_years[x_row, "nets"], nrow(x_row))
cum_nets <- t(apply(nets, 1, cumsum))
cell_xy <- terra::xyFromCell(mask, f$unique_cells)

n <- df$mosquito_number
d <- tibble(class = df$insecticide_class,
            lon = cell_xy[df$cell_id, 1], lat = cell_xy[df$cell_id, 2],
            residual = qlogis((df$died + 0.5) / (n + 1)) - predicted,
            cum_nets = cum_nets[cbind(df$cell_id, df$year_id)],
            w = 1 / ((1 / (df$died + 0.5) + 1 / (n - df$died + 0.5)) *
                       (1 + (n - 1) * rho))) |>
  filter(class == "Pyrethroids") |>
  mutate(w = w / mean(w))

stem <- if (gap > 1) sprintf("%s_gap%d", label, gap) else label
pairs <- read.csv(file.path("outputs/hinge_screen",
                            sprintf("%s_change_pairs.csv", stem))) |>
  filter(insecticide_class == "Pyrethroids")
span <- pairs$year_2 - pairs$year_1
elogit_var <- function(died, n, rho) {
  (1 / (died + 0.5) + 1 / (n - died + 0.5)) * (1 + (n - 1) * rho)
}
if (!is.na(logits_file)) {
  # the predicted change in mean logit between each pair's pixel-years
  key <- paste(df$cell_id, df$type_id, df$year_id)
  group_logit <- tapply(predicted, key, mean)
  pairs$change_pred <-
    group_logit[paste(pairs$cell_id, pairs$type_id, pairs$year_id_2)] -
    group_logit[paste(pairs$cell_id, pairs$type_id, pairs$year_id_1)]
}
pairs <- pairs |>
  mutate(w = span ^ 2 / (elogit_var(died_1, n_1, rho) +
                           elogit_var(died_2, n_2, rho)),
         w = w / mean(w),
         misfit = (change_obs - change_pred) / span,
         lon = cell_xy[cell_id, 1], lat = cell_xy[cell_id, 2])

k <- 60
level_fit <- bam(residual ~ s(lon, lat, k = k), data = d, weights = w,
                 discrete = TRUE)
change_fit <- gam(misfit ~ s(lon, lat, k = 40), data = pairs, weights = w,
                  method = "REML")
slope_fit <- bam(residual ~ s(lon, lat, k = k) +
                   s(lon, lat, by = cum_nets, k = k),
                 data = d, weights = w, discrete = TRUE)

# predictions on the mask (every 4th cell each way, for plotting)
grid <- as.data.frame(terra::aggregate(mask, 4), xy = TRUE, na.rm = TRUE) |>
  transmute(lon = x, lat = y)
grid$level <- predict(level_fit, grid)
grid$change <- predict(change_fit, grid)
slope_terms <- predict(slope_fit, mutate(grid, cum_nets = 1),
                       type = "terms")
grid$slope <- slope_terms[, grep(":cum_nets", colnames(slope_terms),
                                  fixed = TRUE)]
# blank cells more than 2.5 degrees from any data, where the smooths only
# extrapolate
near_data <- function(points) {
  xy <- distinct(points, lon, lat)
  nearest <- vapply(seq_len(nrow(grid)), function(i) {
    min((xy$lon - grid$lon[i]) ^ 2 + (xy$lat - grid$lat[i]) ^ 2)
  }, numeric(1))
  sqrt(nearest) <= 2.5
}
level_near <- near_data(d)
grid$level[!level_near] <- NA
grid$slope[!level_near] <- NA
grid$change[!near_data(pairs)] <- NA

panel <- function(column, points, title, limit, legend) {
  ggplot(grid, aes(lon, lat, fill = .data[[column]])) +
    geom_raster() +
    geom_point(data = distinct(points, lon, lat), aes(lon, lat),
               inherit.aes = FALSE, size = 0.15, alpha = 0.4) +
    scale_fill_distiller(palette = "RdBu", direction = -1,
                         limits = c(-limit, limit), oob = scales::squish,
                         na.value = "grey92") +
    coord_equal() +
    labs(x = NULL, y = NULL, fill = legend, title = title) +
    theme_minimal(base_size = 9)
}
p <- panel("level", d,
           "1. Level residual (observed minus predicted logit mortality)",
           2, "logit") +
  panel("change", pairs,
        sprintf("2. Change misfit, %d-year pairs: predicted minus observed fall per year", gap), 0.5,
        "logit / year") +
  panel("slope", d,
        "3. Slope of level residual on cumulative net use", 1,
        "logit per\nnet-year") +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = sprintf("Pyrethroid misfit of %s%s, smoothed over space", label,
                    if (is.na(logits_name)) "" else
                      sprintf(", net coefficient refitted with w = %s",
                              logits_name)),
    subtitle = paste("Red: the fit predicts too much resistance (1, 3) or",
                     "selection (2); blue: too little. Grey: over 2.5",
                     "degrees from data. Points: data locations"))
ggsave(sprintf("figures/misfit_maps_%s%s.png", label,
               if (is.na(logits_name)) "" else paste0("_w", logits_name)), p, width = 11,
       height = 11, dpi = 130, bg = "white")
cat("edf: level", round(sum(level_fit$edf), 1), "change",
    round(sum(change_fit$edf), 1), "slope model", round(sum(slope_fit$edf), 1),
    "\n")
