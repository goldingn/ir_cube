# Diagnostic maps of the West's LLIN-pyrethroid bioassays, fitted mortality
# and misfit (#47): is the West's plateau in the regional trend figure
# (R/species_compare.R) misplaced because the bioassays move between
# 2010-2015 and 2019 onwards, and can the West be split into subregions that
# behave differently?
#
#   Rscript R/west_figures.R
#
# The bioassays are the modelled LLIN-pyrethroid bioassays (alpha-cypermethrin,
# deltamethrin and permethrin, as the trend figure) in the West
# (analysis_region(), R/species_fit_helpers.R: UNSD Western Africa), in two
# windows, 2010-2015 and 2019 onwards (to 2024, the last year of data). A
# site is a model cell; its observed mortality in a window is its pooled died
# over tested. From outputs/species_runs/misfit/<label>_bioassays.csv
# (R/species_misfit.R) for ref_f0, V3f, V4_class and V5, and the grids of
# R/west_predictions.R (outputs/species_runs/west/grid_<label>.rds). Every map
# shows the limits of transmission without water bodies
# (data/clean/pfpr_water_mask.tif, aggregated by 4 as the grid) and country
# borders (data/clean/country_borders.RDS); mortality is on one RdYlBu scale,
# 0-100%, red where mortality is low (resistant).
#
# Writes, in figures/species_runs/west/:
#   sites.png        sites in each window, coloured by observed mortality and
#                    sized by mosquitoes tested; which sites have data in both
#                    windows; bioassays per country and window
#   fitted_2019.png  posterior mean predicted mortality, 2019-2025, per model,
#                    with the 2019+ sites' observed mortality
#   misfit_2019.png, misfit_2010.png
#                    mean logit misfit per site (empirical minus predicted
#                    logit, as R/species_misfit.R: red where more died than
#                    predicted, the model predicting too much resistance)
#   subregions.png   the proposed subregions (west_subregion()), and per
#                    subregion the observed pooled mortality per year with the
#                    models' site-matched predictions (the mean prediction at
#                    that year's bioassays, weighted by mosquitoes tested)
# and in outputs/species_runs/west/:
#   bioassays_by_country.csv, bioassays_by_subregion.csv
#                    bioassays, sites, mosquitoes tested, observed pooled
#                    mortality, and per model the site-matched predicted
#                    mortality and the mean misfit, per window; the
#                    subregions also with the eastern one split at 9 N
# Plain R, under a minute.

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(terra)
  library(sf)
  library(ggtext)
})
source("R/functions.R")
# analysis_region() and report()
source("R/species_fit_helpers.R")

labels <- c("ref_f0", "V3f", "V4_class", "V5")
misfit_labels <- c("V3f", "V4_class", "V5", "ref_f0")
fit_colours <- c(ref_f0 = grey(0.45), V3f = "#E69F00", V4_class = "#CC79A7",
                 V5 = "#009E73")
llin_pyrethroids <- c("Alpha-cypermethrin", "Deltamethrin", "Permethrin")
windows <- c(early = "2010-2015", late = "2019-2024")
figure_dir <- "figures/species_runs/west"
output_dir <- "outputs/species_runs/west"
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dpi <- 120

# The proposed subregions of the West, by country, from the maps of the
# sites, fitted mortality and misfit: the split is west to east, not north to
# south. Mali (11 bioassays from 2019) joins the Atlantic countries, and Togo
# (25) the eastern ones; either could go with the middle subregion. The
# eastern subregion can be split at 9 N (east_split), into the coast, whose
# mortality fell, and the north, Niger most of it from 2019
west_subregion <- function(country) {
  case_when(
    country %in% c("Senegal", "Gambia", "Guinea-Bissau", "Guinea",
                   "Sierra Leone", "Liberia", "Mauritania",
                   "Mali") ~ "Senegal to Liberia, and Mali",
    country %in% c("Côte d’Ivoire", "Burkina Faso", "Ghana") ~
      "Côte d'Ivoire, Burkina Faso, Ghana",
    country %in% c("Togo", "Benin", "Nigeria", "Niger") ~
      "Togo, Benin, Nigeria, Niger")
}
subregion_colours <- c("Senegal to Liberia, and Mali" = "#56B4E9",
                       "Côte d'Ivoire, Burkina Faso, Ghana" = "#D55E00",
                       "Togo, Benin, Nigeria, Niger" = "#009E73")
east_split <- 9


# the bioassays ---------------------------------------------------------------------

read_misfit <- function(label) {
  read.csv(file.path("outputs/species_runs/misfit",
                     sprintf("%s_bioassays.csv", label)))
}
bioassays <- read_misfit(labels[1]) %>%
  select(longitude, latitude, cell, year_start, country_name, region,
         insecticide_type, species, died, mosquito_number)
for (label in labels) {
  m <- read_misfit(label)
  stopifnot(identical(m$cell, bioassays$cell),
            identical(m$year_start, bioassays$year_start),
            identical(m$insecticide_type, bioassays$insecticide_type))
  bioassays[[paste0("predicted_", label)]] <- m$predicted
  bioassays[[paste0("misfit_", label)]] <- m$misfit
}
mask <- rast("data/clean/raster_mask.tif")
west <- bioassays %>%
  filter(region == "West", insecticide_type %in% llin_pyrethroids) %>%
  mutate(window = case_when(year_start %in% 2010:2015 ~ "early",
                            year_start >= 2019 ~ "late"),
         subregion = west_subregion(country_name))
stopifnot(!anyNA(west$subregion))
cell_xy <- terra::xyFromCell(mask, west$cell)
west$cell_x <- cell_xy[, 1]
west$cell_y <- cell_xy[, 2]
in_windows <- filter(west, !is.na(window))

# the summaries of a group of bioassays: counts, observed pooled mortality,
# and per model the site-matched prediction (the mean prediction at the
# bioassays, weighted by mosquitoes tested) and the mean misfit
summarise_bioassays <- function(data) {
  data %>%
    summarise(bioassays = n(), sites = n_distinct(cell),
              tested = sum(mosquito_number),
              observed = sum(died) / sum(mosquito_number),
              across(starts_with("predicted_"),
                     ~ sum(.x * mosquito_number) / sum(mosquito_number)),
              across(starts_with("misfit_"), mean),
              .groups = "drop")
}
by_country <- in_windows %>%
  group_by(subregion, country = country_name, window = windows[window]) %>%
  summarise_bioassays() %>%
  arrange(subregion, country, window)
east_name <- names(subregion_colours)[3]
by_subregion <- bind_rows(
  in_windows %>% group_by(subregion, window = windows[window]) %>%
    summarise_bioassays(),
  in_windows %>%
    filter(subregion == east_name) %>%
    mutate(subregion = paste0(east_name, ifelse(cell_y < east_split,
                                                 ": south of ", ": north of "),
                              east_split, " N")) %>%
    group_by(subregion, window = windows[window]) %>%
    summarise_bioassays(),
  in_windows %>% group_by(window = windows[window]) %>%
    summarise_bioassays() %>% mutate(subregion = "West")) %>%
  arrange(factor(subregion, unique(c(names(subregion_colours), subregion,
                                     "West"))), window)
write.csv(by_country, file.path(output_dir, "bioassays_by_country.csv"),
          row.names = FALSE)
write.csv(by_subregion, file.path(output_dir, "bioassays_by_subregion.csv"),
          row.names = FALSE)
options(width = 220)
print(as.data.frame(by_country), digits = 2)
print(as.data.frame(by_subregion), digits = 2)

# per site and window (a cell's country and subregion are those of its first
# bioassay)
sites <- in_windows %>%
  group_by(window, cell) %>%
  summarise(x = first(cell_x), y = first(cell_y),
            country_name = first(country_name), subregion = first(subregion),
            bioassays = n(), tested = sum(mosquito_number),
            observed = sum(died) / sum(mosquito_number),
            across(starts_with("misfit_"), mean), .groups = "drop")
coverage <- sites %>%
  group_by(cell, x, y) %>%
  summarise(tested = sum(tested),
            coverage = if (n_distinct(window) == 2) "both windows" else
              paste(windows[first(window)], "only"),
            .groups = "drop") %>%
  mutate(coverage = factor(coverage, c("both windows",
                                       paste(windows, "only"))))
report("sites: %s", paste(names(table(coverage$coverage)),
                          table(coverage$coverage), collapse = ", "))
# the share of each window's mosquitoes tested in each subregion
print(as.data.frame(in_windows %>%
                      group_by(window = windows[window], subregion) %>%
                      summarise(tested = sum(mosquito_number),
                                .groups = "drop_last") %>%
                      mutate(share = tested / sum(tested))), digits = 2)


# the map layers -------------------------------------------------------------------

borders <- readRDS("data/clean/country_borders.RDS")
grids <- lapply(setNames(nm = labels), function(label) {
  readRDS(file.path(output_dir, sprintf("grid_%s.rds", label)))
})
coarse <- terra::aggregate(mask, grids[[1]]$aggregation, fun = "max",
                           na.rm = TRUE)
xlim <- c(-17.9, 15.6)
ylim <- c(4, 20.9)
water_mask <- terra::aggregate(rast("data/clean/pfpr_water_mask.tif"),
                               grids[[1]]$aggregation, fun = "max",
                               na.rm = TRUE)
water_mask <- terra::crop(water_mask, terra::ext(xlim[1] - 1, xlim[2] + 1,
                                                 ylim[1] - 1, ylim[2] + 1))
background <- as.data.frame(water_mask, xy = TRUE, na.rm = TRUE)
tile <- res(water_mask)

# the land in light grey, the limits of transmission a little darker, then
# `layers`, then the borders
base_map <- function(layers = NULL, title = NULL) {
  ggplot() +
    geom_sf(data = borders, fill = grey(0.94), colour = NA) +
    geom_tile(aes(x, y), data = background, fill = grey(0.83),
              width = tile[1], height = tile[2]) +
    layers +
    geom_sf(data = borders, fill = NA, colour = grey(0.25),
            linewidth = 0.25) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_ir_maps() +
    theme(plot.title = element_markdown(size = 11),
          panel.background = element_rect(fill = "#F4F8FB", colour = NA))
}
annotation_theme <- theme(plot.caption = element_text(hjust = 0, size = 9),
                          plot.title = element_text(size = 13))

# one mortality scale for every figure: RdYlBu, red where mortality is low
mortality_colours <- RColorBrewer::brewer.pal(11, "RdYlBu")
mortality_scale <- function(name = "observed or<br>predicted<br>mortality") {
  scale_fill_gradientn(colours = mortality_colours, limits = c(0, 1),
                       labels = scales::percent, oob = scales::squish,
                       name = name)
}
size_scale <- function(range = c(2, 7.5)) {
  scale_size(range = range, limits = c(0, 5200),
             breaks = c(100, 500, 1500, 4000), name = "mosquitoes<br>tested")
}
# sites as filled points with dark outlines, the largest drawn first
site_points <- function(data, fill) {
  geom_point(aes(x, y, fill = .data[[fill]], size = tested),
             data = arrange(data, desc(tested)), shape = 21,
             colour = grey(0.1), stroke = 0.45)
}


# 1. the sites ---------------------------------------------------------------------

site_panel <- function(w) {
  data <- filter(sites, window == w)
  base_map(site_points(data, "observed"),
           sprintf("%s: %d sites, %d bioassays, %s mosquitoes", windows[[w]],
                   nrow(data), sum(data$bioassays),
                   format(sum(data$tested), big.mark = ","))) +
    mortality_scale("observed<br>mortality") + size_scale()
}
coverage_panel <- base_map(
  geom_point(aes(x, y, fill = coverage, size = tested),
             data = arrange(coverage, desc(tested)), shape = 21,
             colour = grey(0.1), stroke = 0.45),
  sprintf("Sites by window: %s", paste(sprintf(
    "%d %s", table(coverage$coverage), levels(coverage$coverage)),
    collapse = ", "))) +
  scale_fill_manual(values = c("#7B3294", "#9ECAE1", "#F4A261"),
                    name = "data in",
                    guide = guide_legend(override.aes = list(size = 4))) +
  size_scale()
country_order <- by_country %>%
  group_by(country) %>%
  summarise(subregion = first(subregion), total = sum(bioassays)) %>%
  arrange(desc(factor(subregion, names(subregion_colours))), total)
country_bars <- by_country %>%
  mutate(country = factor(country, country_order$country)) %>%
  ggplot(aes(bioassays, country, fill = window)) +
  geom_col(position = position_dodge(width = 0.8, preserve = "single"),
           width = 0.75) +
  geom_text(aes(label = sprintf("%d (%.0f%%)", bioassays, 100 * observed)),
            position = position_dodge(width = 0.8, preserve = "single"),
            hjust = -0.08, size = 2.7) +
  scale_fill_manual(values = c("2010-2015" = "#9ECAE1",
                               "2019-2024" = "#F4A261"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = "LLIN-pyrethroid bioassays (observed pooled mortality)", y = NULL,
       title = "Bioassays per country and window, west to east by subregion") +
  theme_minimal() +
  theme(plot.title = element_text(size = 11), legend.position = "top",
        panel.grid.major.y = element_blank())
sites_figure <- (site_panel("early") + site_panel("late") +
                   plot_layout(guides = "collect")) /
  (coverage_panel + country_bars + plot_layout(widths = c(1.45, 1))) +
  plot_annotation(
    title = "West Africa: LLIN-pyrethroid bioassay sites, 2010-2015 and 2019 onwards",
    caption = paste(
      "Sites are model cells (5 km). Colour: observed mortality pooled over",
      "the site's alpha-cypermethrin, deltamethrin and permethrin bioassays in",
      "the window (died / tested), red where mortality is low; size:",
      "mosquitoes tested.\nDarker grey: the limits of transmission without",
      "water bodies (pfpr_water_mask). The data end in 2024."),
    theme = annotation_theme)
ggsave(file.path(figure_dir, "sites.png"), sites_figure, width = 15,
       height = 9, dpi = dpi, bg = "white")


# 2. fitted mortality, 2019-2025 ---------------------------------------------------

late_sites <- filter(sites, window == "late")
west_late <- filter(by_subregion, subregion == "West",
                    window == windows[["late"]])
fitted_panel <- function(label) {
  g <- grids[[label]]
  data <- g$grid %>% filter(window == "2019-2025")
  centre <- terra::xyFromCell(coarse, terra::cellFromXY(
    coarse, cbind(data$longitude, data$latitude)))
  data$x <- centre[, 1]
  data$y <- centre[, 2]
  base_map(list(geom_tile(aes(x, y, fill = mean), data = data,
                          width = tile[1], height = tile[2]),
                site_points(late_sites, "observed")),
           sprintf("%s: map mean %.0f%%; at the 2019+ bioassays %.0f%% (observed %.0f%%)",
                   label, 100 * mean(data$mean),
                   100 * west_late[[paste0("predicted_", label)]],
                   100 * west_late$observed)) +
    mortality_scale() + size_scale(c(1.6, 6))
}
weights <- grids[[1]]$type_weights
fitted_figure <- wrap_plots(lapply(labels, fitted_panel), ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Posterior mean predicted mortality to the LLIN pyrethroids, 2019-2025, with the 2019+ bioassay sites",
    caption = paste0(
      "Tiles: the mean over 2019-2025 and over the three pyrethroids, ",
      "weighted by mosquitoes tested in the West as the trend figure (",
      paste(sprintf("%s %.2f", tolower(weights$insecticide_type),
                    weights$weight), collapse = ", "),
      "), of ", grids[[1]]$n_draws, " posterior draws (V5: chains ",
      toString(grids[["V5"]]$chains), "),\nat every 4th model cell ",
      "(0.17 degrees), with each model's own covariates (V3f: the whole ",
      "complex at the arabiensis share r(x); V4_class: kdr; V5: the smooths). ",
      "Points: observed mortality pooled per site, 2019-2024, same scale.\n",
      "Title: the map's mean over the West's cells; the mean prediction at ",
      "the 2019+ bioassays, weighted by mosquitoes tested (site-matched); ",
      "observed pooled mortality there."),
    theme = annotation_theme)
ggsave(file.path(figure_dir, "fitted_2019.png"), fitted_figure, width = 15,
       height = 8.6, dpi = dpi, bg = "white")


# 3. misfit per site ---------------------------------------------------------------

misfit_limit <- 3
misfit_scale <- scale_fill_gradientn(
  colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
  limits = c(-misfit_limit, misfit_limit), oob = scales::squish,
  name = paste0("mean misfit<br>(logit)<br><br>",
                "red: more died<br>than predicted<br>(model too<br>",
                "resistant)<br><br>blue: fewer died<br>than predicted<br>",
                "(model not<br>resistant enough)<br>"))
misfit_figure <- function(w) {
  data <- filter(sites, window == w)
  panels <- lapply(misfit_labels, function(label) {
    column <- paste0("misfit_", label)
    mean_misfit <- mean(in_windows[[column]][in_windows$window == w])
    base_map(site_points(data, column),
             sprintf("%s%s: mean misfit %.2f over %d bioassays", label,
                     if (label == "ref_f0") " (reference, no floor)" else "",
                     mean_misfit, sum(in_windows$window == w))) +
      misfit_scale + size_scale(c(1.6, 6.5))
  })
  wrap_plots(panels, ncol = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("LLIN-pyrethroid level misfit per site, %s",
                      windows[[w]]),
      caption = paste(
        "Misfit: empirical logit, log((died + 0.5) / (survived + 0.5)),",
        "minus the logit of the posterior mean prediction",
        "(R/species_misfit.R), averaged over the site's bioassays in the",
        "window;\ncolours squished at +-3; size: mosquitoes tested. Title:",
        "the mean over the West's bioassays in the window. V5: chains 1 and",
        "2."),
      theme = annotation_theme)
}
ggsave(file.path(figure_dir, "misfit_2019.png"), misfit_figure("late"),
       width = 15, height = 8.6, dpi = dpi, bg = "white")
ggsave(file.path(figure_dir, "misfit_2010.png"), misfit_figure("early"),
       width = 15, height = 8.6, dpi = dpi, bg = "white")


# 4. the subregions ----------------------------------------------------------------

west_borders <- borders %>%
  filter(region == "Western Africa") %>%
  mutate(subregion = west_subregion(country_name)) %>%
  filter(!is.na(subregion))
counts <- by_subregion %>%
  filter(subregion %in% names(subregion_colours)) %>%
  select(subregion, window, bioassays, sites) %>%
  pivot_wider(names_from = window, values_from = c(bioassays, sites))
subregion_labels <- setNames(
  sprintf("%s<br>%s: %d bioassays, %d sites<br>%s: %d bioassays, %d sites",
          counts$subregion,
          windows[["early"]], counts[[paste0("bioassays_", windows[["early"]])]],
          counts[[paste0("sites_", windows[["early"]])]],
          windows[["late"]], counts[[paste0("bioassays_", windows[["late"]])]],
          counts[[paste0("sites_", windows[["late"]])]]),
  counts$subregion)
subregion_map <- base_map(
  list(geom_sf(aes(fill = subregion), data = west_borders, alpha = 0.45,
               colour = NA, inherit.aes = FALSE),
       geom_hline(yintercept = east_split, linetype = "dashed",
                  colour = grey(0.3), linewidth = 0.3),
       geom_point(aes(x, y, shape = coverage), data = coverage,
                  size = 1.2, colour = grey(0.1))),
  sprintf("Proposed subregions (dashed: the optional split of the east at %d N)",
          east_split)) +
  scale_fill_manual(values = subregion_colours, labels = subregion_labels,
                    breaks = names(subregion_colours), name = NULL) +
  scale_shape_manual(values = c(16, 1, 4), name = "sites with data in") +
  theme(legend.text = element_markdown(size = 9),
        legend.key.spacing.y = unit(6, "pt"))

# per subregion and year: observed pooled mortality and the site-matched
# predictions
yearly <- west %>%
  filter(year_start >= 2008) %>%
  group_by(subregion, year = year_start) %>%
  summarise_bioassays()
yearly_long <- yearly %>%
  select(subregion, year, starts_with("predicted_")) %>%
  pivot_longer(starts_with("predicted_"), names_to = "model",
               names_prefix = "predicted_", values_to = "predicted") %>%
  mutate(model = factor(model, labels))
# the mean misfit of each model in each window, as text
misfit_note <- function(s) {
  rows <- filter(by_subregion, subregion == s)
  paste(vapply(names(windows), function(w) {
    row <- rows[rows$window == windows[[w]], ]
    sprintf("%s misfit: %s", windows[[w]],
            paste(sprintf("%s %.2f", misfit_labels,
                          unlist(row[paste0("misfit_", misfit_labels)])),
                  collapse = ", "))
  }, character(1)), collapse = "\n")
}
trend_panels <- lapply(names(subregion_colours), function(s) {
  ggplot(mapping = aes(year)) +
    annotate("rect", xmin = c(2009.5, 2018.5), xmax = c(2015.5, 2024.5),
             ymin = -Inf, ymax = Inf, fill = grey(0.94)) +
    geom_line(aes(y = predicted, colour = model),
              data = filter(yearly_long, subregion == s), linewidth = 0.7) +
    geom_point(aes(y = observed, size = tested),
               data = filter(yearly, subregion == s), shape = 21,
               fill = subregion_colours[[s]], colour = grey(0.1),
               stroke = 0.4) +
    annotate("text", x = 2007.6, y = 0.01, hjust = 0, vjust = 0, size = 2.5,
             label = misfit_note(s), lineheight = 0.95) +
    scale_colour_manual(values = fit_colours, name = "site-matched\nprediction") +
    scale_size_area(max_size = 6, limits = c(0, 40000),
                    breaks = c(1000, 5000, 20000),
                    name = "mosquitoes tested",
                    guide = guide_legend(
                      override.aes = list(fill = grey(0.6)))) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = NULL, y = "LLIN-pyrethroid mortality", title = s) +
    theme_minimal() +
    theme(plot.title = element_text(size = 11,
                                    colour = subregion_colours[[s]]),
          panel.grid.minor = element_blank())
})
subregion_figure <- subregion_map /
  (wrap_plots(trend_panels, nrow = 1) + plot_layout(guides = "collect")) +
  plot_layout(heights = c(1.25, 1)) +
  plot_annotation(
    title = "Proposed subregions of the West, and their LLIN-pyrethroid mortality by year",
    caption = paste(
      "Points: observed mortality pooled over the subregion's bioassays in",
      "the year (died / tested), sized by mosquitoes tested. Lines: each",
      "model's mean prediction at the same bioassays, weighted by mosquitoes",
      "tested (site-matched, so it\nfollows where the bioassays are; not the",
      "trend figure's fixed-weight average). Shaded: the windows 2010-2015",
      "and 2019-2024. Misfit: empirical minus predicted logit, negative where",
      "fewer died than predicted. V5: chains 1 and 2."),
    theme = annotation_theme)
ggsave(file.path(figure_dir, "subregions.png"), subregion_figure,
       width = 15, height = 10, dpi = dpi, bg = "white")
report("written %s", toString(file.path(figure_dir, c(
  "sites.png", "fitted_2019.png", "misfit_2019.png", "misfit_2010.png",
  "subregions.png"))))
