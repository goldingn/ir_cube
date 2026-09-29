# How the two changes of #26 alter the net covariate: the use layers (legacy to
# MITN run 06) and the pyrethroid-only share (R/prep_net_use_pyrethroid.R).
#   A  legacy use, all nets (net_use_cube.tif; the fits to date)
#   B  run 06 use, all nets (net_use_run06_cube.tif)
#   C  run 06 use x pyrethroid-only share, w = 0.25 (the new default) and
#      w = 0.47
# Writes figures/net_covariate_maps.png, net_covariate_differences.png,
# net_covariate_series.png and net_covariate_scatter.png.

source("R/packages.R")
source("R/functions.R")
source("R/model_covariates.R")

years <- 2000:2024
map_years <- c(2005, 2012, 2018, 2024)
scatter_years <- c(2010, 2022)
use_label <- "Proportion of people using a net"

covariates <- c(
  "A: legacy use, all nets" = "data/clean/net_use_cube.tif",
  "B: run 06 use, all nets" = "data/clean/net_use_run06_cube.tif",
  "C: run 06 use, pyrethroid-only, w = 0.25" =
    "data/clean/net_use_pyrethroid_cube_w0.25.tif",
  "C: run 06 use, pyrethroid-only, w = 0.47" =
    "data/clean/net_use_pyrethroid_cube_w0.47.tif"
)
cubes <- lapply(covariates, rast)
short <- setNames(c("A", "B", "C (w = 0.25)", "C (w = 0.47)"),
                  names(covariates))

mask <- rast("data/clean/raster_mask.tif")
mask_cells <- which(!is.na(terra::values(mask, mat = FALSE)))

bioassays <- readRDS("data/clean/all_gambiae_complex_data.RDS") %>%
  mutate(cell = terra::cellFromXY(mask, cbind(longitude, latitude)),
         year = pmin(pmax(year_start, min(years)), max(years))) %>%
  filter(!is.na(cell), cell %in% mask_cells) %>%
  distinct(cell, year)

# bioassay pixel-years at the cells where run 06 has no use, set to 0
filled <- rast("data/clean/net_use_run06_filled.tif")
filled_at <- as.matrix(terra::extract(filled, bioassays$cell))[
  cbind(seq_len(nrow(bioassays)), bioassays$year - min(years) + 1)]
cat("bioassay pixel-years:", nrow(bioassays), "; in cells where run 06 use",
    "is set to 0:", sum(filled_at == 1, na.rm = TRUE), "\n")

at_bioassays <- bind_rows(lapply(names(cubes), function(name) {
  x <- as.matrix(terra::extract(cubes[[name]], bioassays$cell))
  tibble(covariate = name, cell = bioassays$cell, year = bioassays$year,
         value = x[cbind(seq_len(nrow(x)), bioassays$year - min(years) + 1)])
}))


# maps -------------------------------------------------------------------------

map_frame <- function(layers) {
  layers <- terra::aggregate(layers, fact = 3, fun = "mean", na.rm = TRUE)
  as.data.frame(layers, xy = TRUE, na.rm = FALSE) %>%
    pivot_longer(-c(x, y), names_to = "layer", values_to = "value") %>%
    filter(!is.na(value))
}
map_theme <- theme_void() +
  theme(strip.text = element_text(size = 9),
        legend.position = "bottom",
        legend.key.width = unit(1.5, "cm"))

levels_panel <- c("A", "B", "C (w = 0.25)")
map_layers <- do.call(c, lapply(names(covariates)[1:3], function(name) {
  layers <- cubes[[name]][[paste0("nets_", map_years)]]
  names(layers) <- paste(short[[name]], map_years, sep = "|")
  layers
}))
maps <- map_frame(map_layers) %>%
  separate(layer, c("panel", "year"), sep = "\\|") %>%
  mutate(panel = factor(panel, levels_panel))
p_maps <- ggplot(maps, aes(x, y, fill = value)) +
  geom_raster() +
  facet_grid(panel ~ year, switch = "y") +
  scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0, 1),
                       name = use_label) +
  coord_equal() +
  labs(caption = paste("A: legacy use, all nets. B: run 06 use, all nets.",
                       "C: run 06 use x pyrethroid-only share (w = 0.25).")) +
  map_theme
ggsave("figures/net_covariate_maps.png", p_maps, width = 11, height = 9,
       bg = "white")

difference_layers <- do.call(c, lapply(map_years, function(year) {
  layer <- paste0("nets_", year)
  out <- c(cubes[[2]][[layer]] - cubes[[1]][[layer]],
           cubes[[3]][[layer]] - cubes[[2]][[layer]])
  names(out) <- paste(c("B - A", "C - B"), year, sep = "|")
  out
}))
differences <- map_frame(difference_layers) %>%
  separate(layer, c("panel", "year"), sep = "\\|")
limit <- max(abs(differences$value))
p_differences <- ggplot(differences, aes(x, y, fill = value)) +
  geom_raster() +
  facet_grid(panel ~ year, switch = "y") +
  scale_fill_distiller(palette = "RdBu", direction = 1,
                       limits = c(-limit, limit),
                       name = paste("Difference in", tolower(use_label))) +
  coord_equal() +
  labs(caption = paste("B - A: run 06 use less legacy use (all nets).",
                       "C - B: pyrethroid-only (w = 0.25) less all nets,",
                       "both run 06 use.")) +
  map_theme
ggsave("figures/net_covariate_differences.png", p_differences, width = 11,
       height = 6.5, bg = "white")


# time series ------------------------------------------------------------------

band_summary <- function(data) {
  data %>%
    group_by(covariate, year) %>%
    summarise(q05 = quantile(value, 0.05, na.rm = TRUE),
              q25 = quantile(value, 0.25, na.rm = TRUE),
              median = median(value, na.rm = TRUE),
              q75 = quantile(value, 0.75, na.rm = TRUE),
              q95 = quantile(value, 0.95, na.rm = TRUE),
              .groups = "drop")
}
cell_bands <- bind_rows(lapply(names(cubes), function(name) {
  v <- terra::values(cubes[[name]], mat = TRUE)[mask_cells, ]
  q <- apply(v, 2, quantile, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
  tibble(covariate = name, year = years, q05 = q[1, ], q25 = q[2, ],
         median = q[3, ], q75 = q[4, ], q95 = q[5, ])
}))
series <- bind_rows(
  band_summary(at_bioassays) %>%
    mutate(where = "At bioassay pixel-years"),
  cell_bands %>% mutate(where = "All mask cells")
) %>%
  mutate(covariate = factor(short[covariate], short),
         where = factor(where, c("At bioassay pixel-years",
                                 "All mask cells")))
p_series <- ggplot(series, aes(year)) +
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = "#2a78d6", alpha = 0.2) +
  geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#2a78d6", alpha = 0.4) +
  geom_line(aes(y = median), colour = "#1c5cab", linewidth = 0.7) +
  facet_grid(where ~ covariate) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = use_label,
       caption = paste("Line: median. Bands: 50% (dark) and 90% (light)",
                       "intervals over cells. Bioassay years before 2000",
                       "are shown as 2000.\nA: legacy use, all nets.",
                       "B: run 06 use, all nets. C: run 06 use x",
                       "pyrethroid-only share.")) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())
ggsave("figures/net_covariate_series.png", p_series, width = 11, height = 6,
       bg = "white")


# scatter at bioassay pixel-years ---------------------------------------------

scatter_data <- at_bioassays %>%
  filter(year %in% scatter_years) %>%
  mutate(covariate = short[covariate]) %>%
  pivot_wider(names_from = covariate, values_from = value) %>%
  mutate(region = cell_regions(cell))
region_colours <- c("Eastern Africa" = "#2a78d6",
                    "Middle Africa" = "#eb6834",
                    "Northern Africa" = "#e87ba4",
                    "Southern Africa" = "#1baf7a",
                    "Western Africa" = "#eda100")
scatter_data$region <- factor(scatter_data$region, names(region_colours))
# the legend on the first panel only (2010 has every region)
scatter_panel <- function(x, y, year_value, x_label, y_label) {
  legend <- if (x == "A" && year_value == scatter_years[1]) "top" else "none"
  ggplot(filter(scatter_data, year == year_value),
         aes(.data[[x]], .data[[y]], colour = region)) +
    geom_abline(colour = "grey60", linetype = "dashed") +
    geom_point(size = 1.6, alpha = 0.8) +
    scale_colour_manual(values = region_colours, drop = FALSE,
                        limits = names(region_colours),
                        name = NULL,
                        guide = guide_legend(nrow = 2,
                                             override.aes = list(size = 3))) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = x_label, y = y_label, title = year_value) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), legend.position = legend)
}
p_scatter <- wrap_plots(
  scatter_panel("A", "B", 2010, "A: legacy use, all nets",
                "B: run 06 use, all nets"),
  scatter_panel("A", "B", 2022, "A: legacy use, all nets",
                "B: run 06 use, all nets"),
  scatter_panel("B", "C (w = 0.25)", 2010, "B: run 06 use, all nets",
                "C: pyrethroid-only, w = 0.25"),
  scatter_panel("B", "C (w = 0.25)", 2022, "B: run 06 use, all nets",
                "C: pyrethroid-only, w = 0.25"),
  ncol = 2
) +
  plot_annotation(caption = paste("Each point a bioassay pixel-year;",
                                  "axes are the proportion of people using",
                                  "a net. Dashed line: equality."))
ggsave("figures/net_covariate_scatter.png", p_scatter, width = 9, height = 10,
       bg = "white")
