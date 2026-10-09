# Fixed-weight regional trends of LLIN-pyrethroid susceptibility in the
# diagnostic regions of R/diagnostic_regions.R (#47): the trend figure of
# R/species_compare.R, with its regions replaced by the diagnostic ones. The
# site-matched predictions of R/region_diagnostics.R follow where the
# bioassays are each year; these hold the cells fixed.
#
#   USE_CHAINS="V5=1,2" Rscript R/region_trends.R [<label>=<fitted_model.RData> ...]
#
# The fits are ref_f0, V3f, V4_class and V5 (`fits` below), or those given.
# For each, from n_draws posterior draws, the same number from each usable
# chain (even_draws(), R/species_fit_helpers.R; with USE_CHAINS, as
# R/species_misfit.R used for V5), the LLIN-pyrethroid (alpha-cypermethrin,
# deltamethrin, permethrin) predicted mortality in every year at each cell
# and insecticide with bioassays in a region, weighted by its share of the
# region's mosquitoes tested over the whole series (constant over the years),
# as figure 1 (R/fig_temporal_preds_data.R) and R/species_compare.R weight
# them. Each country's
# region is in outputs/species_runs/regions/regions.csv. Cached in
# outputs/species_runs/regions/trends_<label>.rds (recomputed when the fit
# is newer or the regions differ).
#
# The data points are figure 1's within each region (pyrethroid_points(), as
# R/species_compare.R): per year, died over tested with each cell reweighted
# to its share over the whole series, then the insecticides reweighted the
# same way; size the effective number of bioassays.
#
# Writes figures/species_runs/regions/region_trends_fixed.png and
# outputs/species_runs/regions/region_trends_fixed.csv.
# Plain R with greta loaded (the draws are greta objects); one fit loaded at
# a time; about 8 GB at peak (V3f) and 2 minutes, then seconds from the
# caches.

suppressMessages({
  library(greta)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ggtext)
})
source("R/functions.R")
source("R/dynamical_predictions.R")
source("R/species_fit_helpers.R")

n_draws <- 500
scratchpad <- paste0("/tmp/claude-1000/-home-nick-Dropbox-github-ir-cube/",
                     "be75c64a-3bb7-4b3e-a81c-c664fe72f5e2/scratchpad/species")
fits <- c(
  ref_f0 = paste0("../ir_cube_netscreen/outputs/pod_jobs/dh270_lin_f0_full/",
                  "temporary/fitted_model.RData"),
  V3f = "outputs/pod_jobs/sp_v3_floor/temporary/fitted_model.RData",
  V4_class = file.path(scratchpad,
                       "local_sp_v4_class/temporary/fitted_model.RData"),
  V5 = file.path(scratchpad, "local_sp_v5/temporary/fitted_model.RData"))
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) > 0) {
  fits <- character(0)
  for (argument in arguments) {
    parts <- strsplit(argument, "=", fixed = TRUE)[[1]]
    stopifnot(length(parts) == 2)
    fits[[parts[1]]] <- parts[2]
  }
}
stopifnot(all(file.exists(fits)))
fit_colours <- c(ref_f0 = grey(0.3), V3f = "#E69F00", V4_class = "#CC79A7",
                 V5 = "#009E73")
fit_colours <- c(fit_colours, V4 = "#0072B2")
llin_pyrethroids <- c("Alpha-cypermethrin", "Deltamethrin", "Permethrin")
output_dir <- "outputs/species_runs/regions"
figure_dir <- "figures/species_runs/regions"
regions <- read.csv(file.path(output_dir, "regions.csv"))
region_counts <- read.csv(file.path(output_dir, "region_counts.csv"))
region_letters <- sort(unique(regions$region))
dpi <- 120


# per fit: the fixed-weight trends per region ------------------------------------

# one row per region, cell and insecticide with LLIN-pyrethroid bioassays,
# weighted by its share of the region's mosquitoes tested
region_weights <- function(df) {
  df %>%
    filter(insecticide_type %in% llin_pyrethroids) %>%
    mutate(region = regions$region[match(country_name, regions$country)]) %>%
    group_by(region, cell_id, type_id, cell) %>%
    summarise(tested = sum(mosquito_number), .groups = "drop") %>%
    group_by(region) %>%
    mutate(weight = tested / sum(tested)) %>%
    ungroup()
}

summarise_fit <- function(label, file) {
  cache <- file.path(output_dir, sprintf("trends_%s.rds", label))
  assignment <- setNames(regions$region, regions$country)
  if (file.exists(cache) && file.mtime(cache) > file.mtime(file)) {
    cached <- readRDS(cache)
    if (identical(cached$assignment, assignment)) return(cached)
  }
  fit <- load_fit(file)
  df <- fit$df
  chosen <- even_draws(fit, n_draws, label)
  parameters <- fit_parameter_draws(fit, chosen$index)
  weights <- region_weights(df)
  stopifnot(!anyNA(weights$region))
  n_times <- max(fit$cell_years_index$year_id)
  pairs <- distinct(weights, cell_id, type_id, cell)
  rows <- pairs[rep(seq_len(nrow(pairs)), n_times), ]
  rows$year_id <- rep(seq_len(n_times), each = nrow(pairs))
  time <- system.time(
    p <- plogis(dynamical_logit(parameters, rows, df, fit$x_cell_years,
                                fit$cell_years_index))
  )[["elapsed"]]
  report("%s: predictions at %d cell-insecticide-years from %d draws in %.0f s",
         label, nrow(rows), nrow(p), time)
  key <- paste(rows$cell_id, rows$type_id)
  trends <- bind_rows(lapply(region_letters, function(r) {
    w <- weights[weights$region == r, ]
    w_rows <- w$weight[match(key, paste(w$cell_id, w$type_id))]
    w_rows[is.na(w_rows)] <- 0
    by_year <- vapply(seq_len(n_times), function(t) {
      at <- rows$year_id == t
      c(p[, at, drop = FALSE] %*% w_rows[at])
    }, numeric(nrow(p)))
    tibble(region = r,
           year = fit$baseline_year - 1 + seq_len(n_times),
           mean = colMeans(by_year),
           lower = apply(by_year, 2, quantile, 0.025),
           upper = apply(by_year, 2, quantile, 0.975))
  })) %>%
    mutate(label = label, .before = 1)
  out <- list(trends = trends, n_draws = nrow(p),
              chains = sort(unique(chosen$chain)), assignment = assignment)
  saveRDS(out, cache)
  report("%s saved; peak memory %.1f GB", label, peak_memory_gb())
  rm(fit, parameters, p)
  invisible(gc())
  out
}

summaries <- lapply(setNames(nm = names(fits)), function(label) {
  summarise_fit(label, fits[[label]])
})


# the data points ---------------------------------------------------------------

# Figure 1's data points (R/fig_temporal_preds_data.R, df_overall_plot then
# df_pyrethroids_plot), for the records `data` of one region; as
# pyrethroid_points() in R/species_compare.R
pyrethroid_points <- function(data) {
  per_insecticide <- data %>%
    filter(insecticide_type %in% llin_pyrethroids) %>%
    group_by(cell, insecticide = insecticide_type, year = year_start) %>%
    summarise(died = sum(died), mosquito_number = sum(mosquito_number),
              bioassays = n(), .groups = "drop") %>%
    group_by(insecticide) %>%
    mutate(total_overall_mosquito_number = sum(mosquito_number)) %>%
    group_by(insecticide, cell) %>%
    mutate(cell_overall_mosquito_number = sum(mosquito_number)) %>%
    group_by(insecticide, year) %>%
    mutate(total_year_mosquito_number = sum(mosquito_number)) %>%
    group_by(insecticide, year, cell) %>%
    mutate(cell_year_mosquito_number = sum(mosquito_number)) %>%
    ungroup() %>%
    mutate(year_fraction = cell_year_mosquito_number /
             total_year_mosquito_number,
           overall_fraction = cell_overall_mosquito_number /
             total_overall_mosquito_number,
           weight = overall_fraction / year_fraction) %>%
    group_by(year, insecticide) %>%
    mutate(relvar_component = (weight - mean(weight)) ^ 2 /
             (mean(weight) ^ 2)) %>%
    summarise(bioassays = sum(bioassays),
              died_weighted = sum(died * weight),
              mosquito_number_weighted = sum(mosquito_number * weight),
              mosquito_number = sum(mosquito_number),
              relvar = mean(relvar_component),
              .groups = "drop") %>%
    mutate(Susceptibility = died_weighted / mosquito_number_weighted,
           design_effect = 1 + relvar,
           effective_bioassays = bioassays / design_effect,
           effective_samples = mosquito_number / design_effect,
           effective_died = Susceptibility * effective_samples)
  per_insecticide %>%
    mutate(total_overall_samples = sum(effective_samples)) %>%
    group_by(insecticide) %>%
    mutate(insecticide_overall_samples = sum(effective_samples)) %>%
    group_by(year) %>%
    mutate(total_year_samples = sum(effective_samples)) %>%
    group_by(insecticide, year) %>%
    mutate(insecticide_year_samples = sum(effective_samples)) %>%
    ungroup() %>%
    mutate(year_fraction = insecticide_year_samples / total_year_samples,
           overall_fraction = insecticide_overall_samples /
             total_overall_samples,
           weight = overall_fraction / year_fraction) %>%
    group_by(year) %>%
    mutate(relvar_component = (weight - mean(weight)) ^ 2 /
             (mean(weight) ^ 2)) %>%
    summarise(died_weighted = sum(effective_died * weight),
              samples_weighted = sum(effective_samples * weight),
              bioassays = sum(effective_bioassays),
              relvar = mean(relvar_component),
              .groups = "drop") %>%
    mutate(Susceptibility = died_weighted / samples_weighted,
           design_effect = 1 + relvar,
           effective_bioassays = bioassays / design_effect)
}

# the modelled bioassays (every fit's are the same)
records <- read.csv("outputs/species_runs/misfit/ref_f0_bioassays.csv") %>%
  mutate(region = regions$region[match(country_name, regions$country)])
points <- bind_rows(lapply(region_letters, function(r) {
  pyrethroid_points(filter(records, region == r)) %>%
    mutate(region = r, .before = 1)
}))


# the figure --------------------------------------------------------------------

trends <- bind_rows(lapply(summaries, `[[`, "trends")) %>%
  mutate(label = factor(label, names(fits)))
write.csv(trends, file.path(output_dir, "region_trends_fixed.csv"),
          row.names = FALSE)
years <- c(2000, 2024)
current_note <- function(r) {
  at <- trends %>% filter(region == r, year %in% 2019:2024) %>%
    group_by(label) %>% summarise(mean = mean(mean), .groups = "drop")
  observed <- points %>% filter(region == r, year %in% 2019:2024)
  keys <- sprintf("<span style='color:%s'>**%s** %.0f%%</span>",
                  fit_colours[as.character(at$label)], at$label,
                  100 * at$mean)
  # three to a line, so that more than four models fit the panel
  lines <- split(keys, ceiling(seq_along(keys) / 3))
  sprintf("2019-2024 mean of the lines:<br>%s",
          paste(vapply(lines, paste, "", collapse = ", "), collapse = "<br>"))
}
trend_panel <- function(r) {
  counts <- region_counts[region_counts$region == r, ]
  title <- sub("^([A-Z]): ", "**\\1**: ", counts$name)
  if (!is.na(counts$attached) && nzchar(counts$attached)) {
    title <- sprintf("%s (+ %s*)", title, counts$attached)
  }
  title <- paste(strwrap(title, 58), collapse = "<br>")
  ggplot(mapping = aes(x = year)) +
    annotate("rect", xmin = c(2009.5, 2018.5), xmax = c(2015.5, 2024.5),
             ymin = -Inf, ymax = Inf, fill = grey(0.94)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = label),
                data = filter(trends, region == r), alpha = 0.15,
                colour = NA) +
    geom_point(aes(y = Susceptibility, size = effective_bioassays),
               data = filter(points, region == r, year >= years[1]),
               shape = 21, fill = grey(0.82), colour = grey(0.1),
               stroke = 0.5) +
    geom_line(aes(y = mean, colour = label),
              data = filter(trends, region == r), linewidth = 0.8) +
    scale_colour_manual(values = fit_colours, name = "fixed-weight\nprediction") +
    scale_fill_manual(values = fit_colours, guide = "none") +
    scale_size_area(limits = c(0, max(points$effective_bioassays)),
                    max_size = 8, breaks = c(5, 25, 100, 300),
                    name = "observed:\neffective bioassays") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.25)) +
    coord_cartesian(xlim = years + c(-0.5, 0.5)) +
    labs(x = NULL, y = "LLIN-pyrethroid mortality", title = title,
         subtitle = current_note(r)) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_markdown(size = 12, lineheight = 1.1),
          plot.subtitle = element_markdown(size = 10),
          panel.grid.minor = element_blank())
}
trend_figure <- wrap_plots(lapply(region_letters, trend_panel), ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "LLIN pyrethroids in the diagnostic regions: fixed-weight predictions and data, as the trend figure",
    caption = paste0(
      "Lines and 95% bands: posterior predicted mortality averaged over each ",
      "region's bioassay cells and the three pyrethroids, each weighted by ",
      "its mosquitoes tested over the whole series (fixed over the years; ",
      "figure 1 per region), from about ", n_draws,
      " draws\n(V5: chains ", toString(summaries[["V5"]]$chains), "). ",
      "Points: figure 1's weighted annual estimates within the region, size ",
      "the effective number of bioassays. Shaded: the windows 2010-2015 and ",
      "2019-2024. Regions: R/diagnostic_regions.R (*: no 2019-2024 ",
      "bioassays)."),
    theme = theme(plot.title = element_text(size = 15),
                  plot.caption = element_text(hjust = 0, size = 10)))
ggsave(file.path(figure_dir, "region_trends_fixed.png"), trend_figure,
       width = 17, height = 15, dpi = dpi, bg = "white")
report("written %s and %s", file.path(figure_dir, "region_trends_fixed.png"),
       file.path(output_dir, "region_trends_fixed.csv"))
