# The chosen kdr surfaces beside the plateau misfit and the An. arabiensis
# fraction (#46, #47), masked to the project's map mask
# (data/clean/pfpr_water_mask.tif: the P. falciparum limits, less water).
#
# kdr surfaces, as chosen by cross-validation (R/kdr_cv.R,
# R/kdr_cv_extrap.R), fitted to all records and predicted at 2015:
#   995F  thin plate, range-space sp 10^1, null space penalised (sp 10^10)
#   995S  thin plate, sp 10^2
# and total kdr, their sum (capped at 1).
#
#   Rscript R/kdr_species_plateau_maps.R <plateau grid csv> <label>
#
# Writes figures/kdr_species_plateau_<label>.png (plateau, arabiensis, total
# kdr) and figures/kdr_variants_<label>.png (995F, 995S).

arguments <- commandArgs(trailingOnly = TRUE)
plateau_file <- arguments[1]
label <- arguments[2]

source("R/packages.R")
source("R/kdr_functions.R")
source("R/kdr_cv_functions.R")

grid <- read.csv(plateau_file)
points <- as.matrix(grid[, c("lon", "lat")])

in_mask <- !is.na(terra::extract(rast("data/clean/pfpr_water_mask.tif"),
                                  points)[[1]])
grid$arabiensis <- terra::extract(
  terra::aggregate(rast("data/clean/arabiensis_fraction.tif"), 4,
                   na.rm = TRUE), points)[[1]]

records <- kdr_cv_records()
chosen <- list(
  `995F` = candidate(tp_smooth, sp = 10 ^ c(1, 10), select = TRUE,
                     type = "tp, null space penalised"),
  `995S` = candidate(tp_smooth, sp = 10 ^ 2, type = "tp"))
new <- data.frame(longitude = grid$lon, latitude = grid$lat, year_c = 0)
for (v in names(chosen)) {
  model <- fit_candidate(chosen[[v]], allele_data(records, v))
  cat(sprintf("%s: edf %.1f, scale %.1f\n", v, spatial_edf(model),
              model$scale))
  grid[[v]] <- as.numeric(predict(model, new, type = "response"))
}

grid$total <- pmin(grid$`995F` + grid$`995S`, 1)
for (column in c("change", "arabiensis", "995F", "995S", "total")) {
  grid[[column]][!in_mask] <- NA
}

base <- function(title) {
  ggplot(grid, aes(lon, lat)) +
    geom_raster(fill = "grey92") +
    coord_equal() +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.text = element_blank())
}
fraction_panel <- function(column, title) {
  base(title) +
    geom_raster(data = filter(grid, !is.na(.data[[column]])),
                aes(fill = .data[[column]])) +
    scale_fill_viridis_c(limits = c(0, 1), name = NULL)
}
plateau_panel <- base(
  "Missing plateau (red): change in level\nresidual, 2012-16 to 2019-25") +
  geom_raster(data = filter(grid, !is.na(change)), aes(fill = change)) +
  scale_fill_distiller(palette = "RdBu", limits = c(-2, 2),
                       oob = scales::squish, name = "logit")
p <- patchwork::wrap_plots(
  plateau_panel,
  fraction_panel("arabiensis",
                 "An. arabiensis fraction of the complex\n(Vector Atlas)"),
  fraction_panel("total", "Total kdr (995F + 995S), 2015\n(CV-chosen smooths)"),
  nrow = 1)
ggsave(sprintf("figures/kdr_species_plateau_%s.png", label), p, width = 16,
       height = 6, dpi = 110, bg = "white")
p <- patchwork::wrap_plots(
  fraction_panel("995F", "Vgsc-995F allele frequency, 2015\n(CV-chosen smooth)"),
  fraction_panel("995S", "Vgsc-995S allele frequency, 2015\n(CV-chosen smooth)"),
  nrow = 1)
ggsave(sprintf("figures/kdr_variants_%s.png", label), p, width = 11,
       height = 6, dpi = 110, bg = "white")
