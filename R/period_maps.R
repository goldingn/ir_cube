# Maps of pyrethroid mortality in a period (#37): where observed mortality
# levels off and where the fit keeps lowering it. 2D smooths over longitude
# and latitude, weighted by precision, predicted to the mask (blank more than
# 2.5 degrees from the period's data), of the bioassays in [from, to]:
#   1. observed mortality (the smooth of the empirical logit, back-transformed);
#   2. the fit's predicted mortality (the smooth of its predicted logit);
#   3. the fit's level residual (empirical minus predicted logit);
#   4. the level residual under the mock-up's logits (R/net_w_screen.R), if
#      given.
# Residual panels: red where the fit predicts too much resistance.
#
#   Rscript R/period_maps.R <fit RData> <label> <from> <to> [logits rds] [name]
#
# Writes figures/period_maps_<label>_<from>_<to>.png.

arguments <- commandArgs(trailingOnly = TRUE)
fit_file <- arguments[1]
label <- arguments[2]
from <- as.integer(arguments[3])
to <- as.integer(arguments[4])
logits_file <- if (length(arguments) >= 5) arguments[5] else NA
logits_name <- if (length(arguments) >= 6) arguments[6] else NA

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
predicted <- colMeans(dynamical_logit(par, select(df, cell_id, country_id,
                                                  type_id, year_id),
                                      df, f$x_cell_years, f$cell_years_index))
mock <- if (is.na(logits_file)) NULL else readRDS(logits_file)[[logits_name]]
rho <- colMeans(par$rho_types)[df$type_id]
cell_xy <- terra::xyFromCell(mask, f$unique_cells)

n <- df$mosquito_number
keep <- df$insecticide_class == "Pyrethroids" &
  df$year_start >= from & df$year_start <= to
d <- tibble(lon = cell_xy[df$cell_id, 1], lat = cell_xy[df$cell_id, 2],
            empirical = qlogis((df$died + 0.5) / (n + 1)),
            predicted = predicted,
            mock = if (is.null(mock)) NA_real_ else mock,
            w = 1 / ((1 / (df$died + 0.5) + 1 / (n - df$died + 0.5)) *
                       (1 + (n - 1) * rho)))[keep, ] |>
  mutate(w = w / mean(w), residual = empirical - predicted,
         mock_residual = empirical - mock)
cat(sprintf("%d pyrethroid bioassays, %d-%d\n", nrow(d), from, to))

grid <- as.data.frame(terra::aggregate(mask, 4), xy = TRUE, na.rm = TRUE) |>
  transmute(lon = x, lat = y)
smooth_on_grid <- function(response) {
  fit <- bam(as.formula(paste(response, "~ s(lon, lat, k = 60)")), data = d,
             weights = w, discrete = TRUE)
  predict(fit, grid)
}
grid$observed <- plogis(smooth_on_grid("empirical"))
grid$fitted <- plogis(smooth_on_grid("predicted"))
grid$residual <- smooth_on_grid("residual")
if (!is.null(mock)) grid$mock_residual <- smooth_on_grid("mock_residual")
xy <- distinct(d, lon, lat)
near <- vapply(seq_len(nrow(grid)), function(i) {
  min((xy$lon - grid$lon[i]) ^ 2 + (xy$lat - grid$lat[i]) ^ 2)
}, numeric(1)) <= 2.5 ^ 2
grid[!near, setdiff(names(grid), c("lon", "lat"))] <- NA

base <- function(column, title) {
  ggplot(grid, aes(lon, lat, fill = .data[[column]])) +
    geom_raster() +
    geom_point(data = xy, aes(lon, lat), inherit.aes = FALSE, size = 0.15,
               alpha = 0.4) +
    coord_equal() +
    labs(x = NULL, y = NULL, title = title) +
    theme_minimal(base_size = 9)
}
mortality_scale <- scale_fill_viridis_c(limits = c(0, 1), na.value = "grey92",
                                        name = "mortality")
residual_scale <- scale_fill_distiller(palette = "RdBu", direction = -1,
                                       limits = c(-2, 2),
                                       oob = scales::squish,
                                       na.value = "grey92", name = "logit")
panels <- list(
  base("observed", "1. Observed mortality") + mortality_scale,
  base("fitted", "2. Fit's predicted mortality") + mortality_scale,
  base("residual", "3. Level residual, the fit") + residual_scale)
if (!is.null(mock)) {
  panels[[4]] <- base("mock_residual",
                      sprintf("4. Level residual, mock-up with w = %s",
                              logits_name)) + residual_scale
}
p <- wrap_plots(panels, ncol = 2) +
  plot_annotation(
    title = sprintf("Pyrethroid bioassays %d-%d, %s", from, to, label),
    subtitle = paste("Smoothed over space, weighted by precision. Residuals:",
                     "red where the fit predicts too much resistance.",
                     "Grey: over 2.5 degrees from data"))
ggsave(sprintf("figures/period_maps_%s_%d_%d.png", label, from, to), p,
       width = 11, height = 11, dpi = 130, bg = "white")
