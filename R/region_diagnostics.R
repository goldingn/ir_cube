# Regional trend diagnostics of the fits of #47 on the diagnostic regions of
# R/diagnostic_regions.R (groups of neighbouring countries with alike current
# LLIN-pyrethroid susceptibility, from the data alone): does each model get
# the broad trend, and above all the current level, right in each part of
# Africa? A plateau at the wrong level is the main failure to look for.
#
#   Rscript R/region_diagnostics.R
#
# The bioassays are the modelled LLIN-pyrethroid bioassays (alpha-cypermethrin,
# deltamethrin and permethrin, as the trend figure, R/species_compare.R) of
# every year, with each fit's posterior mean prediction and logit misfit at
# each, from outputs/species_runs/misfit/<label>_bioassays.csv
# (R/species_misfit.R) for ref_f0, V3f, V4_class and V5 (chains 1 and 2), and
# each country's region from outputs/species_runs/regions/regions.csv (the
# countries without current bioassays in the region they were joined to). The
# windows are current, 2019-2024 (the data end in 2024), and early,
# 2010-2015.
#
# Per region and year (or window): the observed pooled mortality (died over
# tested), and per model the site-matched prediction, the mean prediction at
# the same bioassays weighted by mosquitoes tested (so it follows where and
# when the bioassays are; not the trend figure's fixed-weight average, which
# R/region_trends.R draws), and the mean logit misfit (empirical minus
# predicted logit, as R/species_misfit.R: positive where more died than
# predicted, i.e. the model predicts too much resistance).
#
# Writes
#   figures/species_runs/regions/region_trends.png
#       one panel per region: observed pooled mortality per year, sized by
#       mosquitoes tested, and each model's site-matched prediction for the
#       year as small coloured points (offset a little), joined by thin lines
#   figures/species_runs/regions/current_levels.png
#       per region and window, each model's site-matched prediction minus
#       the observed pooled mortality, and its mean misfit
#   outputs/species_runs/regions/current_levels.csv
#       per region (and Africa) and window: counts, observed pooled
#       mortality, and per model the site-matched prediction, its difference
#       from the observed (percentage points) and the mean misfit
# Plain R, under a minute.

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ggtext)
})
source("R/functions.R")
# report()
source("R/species_fit_helpers.R")

labels <- c("ref_f0", "V3f", "V4_class", "V5")
# or the labels given (Rscript R/region_diagnostics.R <label> ...)
if (length(commandArgs(trailingOnly = TRUE)) > 0) {
  labels <- commandArgs(trailingOnly = TRUE)
}
# Okabe-Ito, as R/west_figures.R; the observed points are light grey
fit_colours <- c(ref_f0 = grey(0.3), V3f = "#E69F00", V4_class = "#CC79A7",
                 V5 = "#009E73")
fit_colours <- c(fit_colours, V4 = "#0072B2")
llin_pyrethroids <- c("Alpha-cypermethrin", "Deltamethrin", "Permethrin")
windows <- c(early = "2010-2015", current = "2019-2024")
years <- 2005:2024
figure_dir <- "figures/species_runs/regions"
output_dir <- "outputs/species_runs/regions"
dpi <- 120
options(width = 220)


# the bioassays ---------------------------------------------------------------------

regions <- read.csv(file.path(output_dir, "regions.csv"))
region_counts <- read.csv(file.path(output_dir, "region_counts.csv"))
read_misfit <- function(label) {
  read.csv(file.path("outputs/species_runs/misfit",
                     sprintf("%s_bioassays.csv", label)))
}
bioassays <- read_misfit(labels[1]) %>%
  select(cell, year_start, country_name, insecticide_type, died,
         mosquito_number)
for (label in labels) {
  m <- read_misfit(label)
  stopifnot(identical(m$cell, bioassays$cell),
            identical(m$year_start, bioassays$year_start),
            identical(m$insecticide_type, bioassays$insecticide_type))
  bioassays[[paste0("predicted_", label)]] <- m$predicted
  bioassays[[paste0("misfit_", label)]] <- m$misfit
}
bioassays <- bioassays %>%
  filter(insecticide_type %in% llin_pyrethroids) %>%
  mutate(region = regions$region[match(country_name, regions$country)],
         site = paste(country_name, cell),
         window = case_when(year_start %in% 2010:2015 ~ "early",
                            year_start %in% 2019:2024 ~ "current"))
stopifnot(!anyNA(bioassays$region))

# counts, observed pooled mortality and, per model, the site-matched
# prediction, its difference from the observed and the mean misfit
summarise_bioassays <- function(data) {
  data %>%
    summarise(bioassays = n(), sites = n_distinct(site),
              tested = sum(mosquito_number),
              observed = sum(died) / sum(mosquito_number),
              across(starts_with("predicted_"),
                     ~ sum(.x * mosquito_number) / sum(mosquito_number)),
              across(starts_with("misfit_"), mean),
              .groups = "drop")
}
long_by_model <- function(summary) {
  summary %>%
    pivot_longer(matches("^(predicted|misfit)_"),
                 names_to = c(".value", "model"),
                 names_pattern = "^(predicted|misfit)_(.*)$") %>%
    mutate(model = factor(model, labels),
           difference = 100 * (predicted - observed))
}


# the levels per region and window ------------------------------------------------

in_windows <- filter(bioassays, !is.na(window))
levels_wide <- bind_rows(
  in_windows %>% group_by(region, window) %>% summarise_bioassays(),
  in_windows %>% group_by(window) %>% summarise_bioassays() %>%
    mutate(region = "Africa")) %>%
  mutate(name = coalesce(region_counts$name[match(region,
                                                  region_counts$region)],
                         "Africa: all regions"),
         window = factor(windows[window], windows)) %>%
  arrange(window, region == "Africa", region)
current_levels <- long_by_model(levels_wide) %>%
  select(region, name, window, bioassays, sites, tested, observed, model,
         predicted, difference, misfit)
write.csv(current_levels, file.path(output_dir, "current_levels.csv"),
          row.names = FALSE)

# printed: per region, observed and each model's prediction (misfit)
printable <- current_levels %>%
  mutate(cell = sprintf("%3.0f%% (%+.2f)", 100 * predicted, misfit)) %>%
  select(window, region, bioassays, sites, observed, model, cell) %>%
  mutate(observed = sprintf("%.0f%%", 100 * observed)) %>%
  pivot_wider(names_from = model, values_from = cell)
print(as.data.frame(printable), right = FALSE)


# 1. the trends per region ------------------------------------------------------

yearly <- bioassays %>%
  filter(year_start %in% years) %>%
  group_by(region, year = year_start) %>%
  summarise_bioassays()
offsets <- setNames(seq(-0.24, 0.24, length.out = length(labels)), labels)
yearly_long <- long_by_model(yearly) %>%
  mutate(x = year + offsets[as.character(model)])
max_tested <- max(yearly$tested)

model_key <- function(values, format) {
  keys <- sprintf("<span style='color:%s'>**%s** %s</span>",
                  fit_colours[labels], labels, sprintf(format, values))
  # three to a line, so that more than four models fit the panel
  lines <- split(keys, ceiling(seq_along(keys) / 3))
  paste(vapply(lines, paste, "", collapse = ", "), collapse = "<br>")
}
trend_panel <- function(r) {
  counts <- region_counts[region_counts$region == r, ]
  current <- filter(levels_wide, region == r, window == windows[["current"]])
  title <- sub("^([A-Z]): ", "**\\1**: ", counts$name)
  if (!is.na(counts$attached) && nzchar(counts$attached)) {
    title <- sprintf("%s (+ %s*)", title, counts$attached)
  }
  title <- paste(strwrap(title, 58), collapse = "<br>")
  subtitle <- sprintf(
    "2019-2024: observed %.0f%% (%d bioassays, %d sites)<br>%s",
    100 * current$observed, current$bioassays, current$sites,
    model_key(100 * unlist(current[paste0("predicted_", labels)]), "%.0f%%"))
  ggplot(mapping = aes(year)) +
    annotate("rect", xmin = c(2009.5, 2018.5), xmax = c(2015.5, 2024.5),
             ymin = -Inf, ymax = Inf, fill = grey(0.94)) +
    geom_point(aes(y = observed, size = tested),
               data = filter(yearly, region == r), shape = 21,
               fill = grey(0.82), colour = grey(0.1), stroke = 0.5) +
    geom_line(aes(x, predicted, colour = model, group = model),
              data = filter(yearly_long, region == r), linewidth = 0.35,
              alpha = 0.8) +
    geom_point(aes(x, predicted, fill = model),
               data = filter(yearly_long, region == r), shape = 21,
               colour = "white", stroke = 0.3, size = 2.4) +
    scale_colour_manual(values = fit_colours, guide = "none") +
    scale_fill_manual(values = fit_colours,
                      name = "site-matched\nprediction",
                      guide = guide_legend(override.aes = list(size = 4))) +
    scale_size_area(max_size = 9, limits = c(0, max_tested),
                    breaks = c(500, 2000, 10000, 30000),
                    labels = scales::comma,
                    name = "observed:\nmosquitoes tested") +
    scale_x_continuous(limits = range(years) + c(-0.5, 0.5),
                       breaks = seq(2005, 2025, by = 5)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.25)) +
    labs(x = NULL, y = "LLIN-pyrethroid mortality", title = title,
         subtitle = subtitle) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_markdown(size = 12, lineheight = 1.1),
          plot.subtitle = element_markdown(size = 10, lineheight = 1.2),
          panel.grid.minor = element_blank())
}
region_letters <- sort(unique(regions$region))
trend_figure <- wrap_plots(lapply(region_letters, trend_panel), ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "LLIN-pyrethroid mortality by year in the diagnostic regions: observed, and each model's site-matched prediction",
    caption = paste(
      "Grey circles: observed mortality pooled over the region's",
      "alpha-cypermethrin, deltamethrin and permethrin bioassays in the year",
      "(died / tested), sized by mosquitoes tested. Coloured points: each",
      "model's posterior mean prediction at the same bioassays,\nweighted by",
      "mosquitoes tested (site-matched: it follows where the bioassays are,",
      "not the trend figure's fixed-weight average; R/region_trends.R),",
      "offset a little to the sides of the year. Shaded: the windows",
      "2010-2015 and 2019-2024.\nRegions: R/diagnostic_regions.R (*: no",
      "2019-2024 bioassays, joined by the longest border). V5: chains 1 and",
      "2."),
    theme = theme(plot.title = element_text(size = 15),
                  plot.caption = element_text(hjust = 0, size = 10)))
ggsave(file.path(figure_dir, "region_trends.png"), trend_figure,
       width = 17, height = 15, dpi = dpi, bg = "white")


# 2. the current level ---------------------------------------------------------------

# a region's short name: its letter and its first two countries
short_names <- c(setNames(vapply(strsplit(sub("^[A-Z]: ", "", region_counts$name),
                                          ", "), function(countries) {
  paste0(paste(head(countries, 2), collapse = ", "),
         if (length(countries) > 2) ", ..." else "")
}, character(1)), region_counts$region), Africa = "all regions")
short_names <- setNames(sprintf("%s: %s", names(short_names), short_names),
                        names(short_names))
level_data <- current_levels %>%
  mutate(region_label = factor(short_names[region], rev(short_names)),
         y = as.numeric(region_label) +
           - offsets[as.character(model)] * 0.9)
level_panel <- function(column, x_label, reference = 0) {
  ggplot(level_data, aes(.data[[column]], y, colour = model)) +
    geom_vline(xintercept = reference, colour = grey(0.4)) +
    geom_point(size = 3.2) +
    facet_wrap(~ window, nrow = 1) +
    scale_y_continuous(breaks = seq_along(levels(level_data$region_label)),
                       labels = levels(level_data$region_label)) +
    scale_colour_manual(values = fit_colours, name = NULL) +
    labs(x = x_label, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(colour = grey(0.9)),
          strip.text = element_text(size = 13, hjust = 0))
}
key_text <- region_counts %>%
  mutate(text = sprintf("%s (%.0f%% observed 2019-2024)", name,
                        100 * observed_current)) %>%
  pull(text)
level_figure <- (level_panel("difference",
                             "site-matched prediction minus observed pooled mortality (percentage points)") /
                   level_panel("misfit",
                               "mean logit misfit (empirical minus predicted; positive: model too resistant)")) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Level misfit per diagnostic region and window",
    caption = paste0(
      paste(strwrap(paste(key_text, collapse = "; "), 190),
            collapse = "\n"),
      "\nTop: the mean prediction at the window's bioassays weighted by ",
      "mosquitoes tested, minus their pooled mortality. Bottom: the mean ",
      "over the bioassays of the empirical logit\nminus the logit of the ",
      "prediction (with many bioassays at 100%, the empirical logit is ",
      "large, so the two can disagree in the far south). V5: chains 1 and ",
      "2."),
    theme = theme(plot.title = element_text(size = 15),
                  plot.caption = element_text(hjust = 0, size = 10)))
ggsave(file.path(figure_dir, "current_levels.png"), level_figure,
       width = 15, height = 11, dpi = dpi, bg = "white")
report("written %s, %s and %s",
       file.path(figure_dir, "region_trends.png"),
       file.path(figure_dir, "current_levels.png"),
       file.path(output_dir, "current_levels.csv"))
