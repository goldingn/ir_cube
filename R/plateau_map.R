# Where the dynamical model's plateau is structurally wrong (#37): the change
# in a fit's pyrethroid level residual (empirical minus predicted logit of
# mortality) from an earlier window to a late one, mapped over space. With no
# floor the model keeps lowering mortality wherever exposure continues, so
# where the data have levelled off the residual grows between the windows:
# the model has fallen further than the data. Fitted as
#
#   residual ~ s(lon, lat) + s(lon, lat, by = late)
#
# late = 1 in the late window, 0 in the early one, weighted by precision; the
# map is the `by` smooth, blank where either window has no data within 2.5
# degrees. Red: the model fell further than the data (a missing plateau);
# blue: the data fell further than the model.
#
#   Rscript R/plateau_map.R <fit RData> <label> [early from-to] [late from-to]
#                           [logits rds] [name]
#
# With a logits file (from R/net_w_screen.R) and the name of one of its
# elements, the predictions are those logits rather than the fit's.
#
# Writes figures/plateau_map_<label>[_w<name>].png and the mapped change as
# outputs/hinge_screen/plateau_grid_<label>[_w<name>].csv, and the pyrethroid
# bioassays' residuals as plateau_points_<label>[_w<name>].csv.

arguments <- commandArgs(trailingOnly = TRUE)
fit_file <- arguments[1]
label <- arguments[2]
early <- if (length(arguments) >= 3) as.integer(strsplit(arguments[3], "-")[[1]]) else
  c(2012L, 2016L)
late <- if (length(arguments) >= 4) as.integer(strsplit(arguments[4], "-")[[1]]) else
  c(2019L, 2025L)
logits_file <- if (length(arguments) >= 5) arguments[5] else NA
logits_name <- if (length(arguments) >= 6) arguments[6] else NA

source("R/packages.R")
source("R/functions.R")
source("R/dynamical_model.R")
source("R/dynamical_predictions.R")
suppressMessages(library(mgcv))

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
cell_xy <- terra::xyFromCell(mask, f$unique_cells)
n <- df$mosquito_number
window <- case_when(between(df$year_start, early[1], early[2]) ~ 0,
                    between(df$year_start, late[1], late[2]) ~ 1)
d <- tibble(lon = cell_xy[df$cell_id, 1], lat = cell_xy[df$cell_id, 2],
            residual = qlogis((df$died + 0.5) / (n + 1)) - predicted,
            late = window, year = df$year_start, species = df$species,
            w = 1 / ((1 / (df$died + 0.5) + 1 / (n - df$died + 0.5)) *
                       (1 + (n - 1) * rho)))[
  df$insecticide_class == "Pyrethroids" & !is.na(window), ] |>
  mutate(w = w / mean(w))
cat(sprintf("pyrethroid bioassays: %d in %d-%d, %d in %d-%d\n",
            sum(d$late == 0), early[1], early[2], sum(d$late == 1), late[1],
            late[2]))
fit <- bam(residual ~ s(lon, lat, k = 60) + s(lon, lat, by = late, k = 60),
           data = d, weights = w, discrete = TRUE)

grid <- as.data.frame(terra::aggregate(mask, 4), xy = TRUE, na.rm = TRUE) |>
  transmute(lon = x, lat = y)
terms <- predict(fit, mutate(grid, late = 1), type = "terms")
grid$change <- terms[, grep(":late", colnames(terms), fixed = TRUE)]
near <- function(points) {
  xy <- distinct(points, lon, lat)
  vapply(seq_len(nrow(grid)), function(i) {
    min((xy$lon - grid$lon[i]) ^ 2 + (xy$lat - grid$lat[i]) ^ 2)
  }, numeric(1)) <= 2.5 ^ 2
}
grid$change[!(near(filter(d, late == 0)) & near(filter(d, late == 1)))] <- NA

p <- ggplot(grid, aes(lon, lat, fill = change)) +
  geom_raster() +
  geom_point(data = distinct(d, lon, lat, late), aes(lon, lat,
                                                     shape = factor(late)),
             inherit.aes = FALSE, size = 0.4, alpha = 0.5) +
  scale_shape_manual(values = c(`0` = 1, `1` = 3),
                     labels = c(sprintf("%d-%d", early[1], early[2]),
                                sprintf("%d-%d", late[1], late[2]))) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-2, 2),
                       oob = scales::squish, na.value = "grey92") +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "logit", shape = NULL,
       title = sprintf(paste("Change in the pyrethroid level residual,",
                             "%d-%d to %d-%d, %s"),
                       early[1], early[2], late[1], late[2],
                       if (is.na(logits_name)) label else
                         sprintf("%s, mock-up with w = %s", label,
                                 logits_name)),
       subtitle = paste("Red: the model fell further than the data (a",
                        "missing plateau); blue: the data fell further.",
                        "Grey: no data in both windows")) +
  theme_minimal()
suffix <- if (is.na(logits_name)) "" else paste0("_w", logits_name)
ggsave(sprintf("figures/plateau_map_%s%s.png", label, suffix), p, width = 8,
       height = 8, dpi = 130, bg = "white")
write.csv(grid, sprintf("outputs/hinge_screen/plateau_grid_%s%s.csv", label,
                        suffix), row.names = FALSE)
write.csv(d, sprintf("outputs/hinge_screen/plateau_points_%s%s.csv", label,
                     suffix), row.names = FALSE)
