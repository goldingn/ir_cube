# Compare fits of #47 in depth: regional pyrethroid trends against the data,
# and the in-sample variance explained in within-pixel change.
#
#   Rscript R/species_compare.R [<label>=<fitted_model.RData> ...]
#
# The fits are those of `fits` below (the reference ref_f0, drawn thin and
# grey, and V3f, V4 and V4_class), or those given. For each, from n_draws
# posterior draws, the same number from each usable chain (even_draws(),
# R/species_fit_helpers.R), cached in outputs/species_runs/compare/<label>.rds:
#
# 1. Regional trends of susceptibility to the LLIN pyrethroids (the
#    pyrethroids in all nets recorded in surveys: alpha-cypermethrin,
#    deltamethrin and permethrin), the first panel of figure 1
#    (R/fig_temporal_preds_data.R) per region (analysis_region(); West,
#    Central, East, Southern, Horn) and for Ethiopia. Figure 1 averages the
#    predictions over the cells with bioassays, not over the map: each
#    insecticide's predicted mortality at each cell weighted by the share of
#    that insecticide's mosquitoes tested there over the whole series, and the
#    insecticides weighted by their mosquitoes tested; together, a weight on
#    each cell and insecticide proportional to its mosquitoes tested, constant
#    over the years. This does the same within each region, from the fits'
#    draws (dynamical_logit(), R/dynamical_predictions.R). For the species
#    model, the prediction at a cell is that of the whole complex, the mixture
#    at the arabiensis fraction r(x), as all_states() gives figure 1. The data
#    points are figure 1's, computed within the region: per year, died over
#    tested with each cell reweighted to its share over the whole series, then
#    the insecticides reweighted the same way, the point size the effective
#    number of bioassays after the weighting's design effect.
#    Writes figures/species_runs/regional_trends_pyrethroids.png.
#
# 2. The in-sample variance explained in within-pixel change, as the
#    "temporal change" experiment's change score (R/validation_change.R) with
#    the variance-explained bars of R/variance_explained.R: for each cell and
#    insecticide with bioassays in both the window before a cut and the window
#    from it (five years each, cuts 2014 and 2018, as the forecasting folds),
#    the observed change in pooled mortality, and the predicted change, the
#    difference of the posterior mean predicted mortality at the window's
#    bioassays, weighted by mosquitoes tested; the two cuts pooled. Explained
#    is 100 (1 - MSE / Var(observed change)); the bioassay noise share is 100
#    times the mean irreducible variance of the observed change, the sum of
#    the two windows' noise_floor_var_pooled() at the replicate-based rho of
#    the type, over Var(observed change). Intervals by resampling pixels
#    (pixel_bootstrap()), and for the noise share from the posterior of rho.
#    In sample: every fit has seen both windows, so this is optimistic, and
#    unlike the CV, where the windows after the cut were held out.
#    Writes outputs/species_runs/insample_change_variance.csv,
#    outputs/species_runs/insample_change_bias.csv (mean observed and
#    predicted change, and their correlation) and
#    figures/species_runs/insample_change_variance.png. For reference, the
#    same quantity out of sample from the CV folds of the model before #47
#    (outputs/cv_change.csv, cell scale, both cuts) is -25% (-33% for the
#    pyrethroids), with the same 44% noise share.
#
# Plain R; one fit loaded at a time; about 4 GB for a species fit.

suppressMessages({
  library(greta)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
source("R/functions.R")
source("R/dynamical_predictions.R")
source("R/species_fit_helpers.R")
# pixel_bootstrap(), and through R/validation_functions.R rho_lookup(),
# rho_for_record() and noise_floor_var_pooled()
source("R/validation_scoring.R")
# bar_layers(), region_key() and base_theme
source("R/fig_variance_bars.R")

n_draws <- 500
scratchpad <- paste0("/tmp/claude-1000/-home-nick-Dropbox-github-ir-cube/",
                     "be75c64a-3bb7-4b3e-a81c-c664fe72f5e2/scratchpad/species")
fits <- c(
  ref_f0 = paste0("../ir_cube_netscreen/outputs/pod_jobs/dh270_lin_f0_full/",
                  "temporary/fitted_model.RData"),
  V3f = "outputs/pod_jobs/sp_v3_floor/temporary/fitted_model.RData",
  V4 = file.path(scratchpad, "local_sp_v4/temporary/fitted_model.RData"),
  V4_class = file.path(scratchpad,
                       "local_sp_v4_class/temporary/fitted_model.RData"))
for (argument in commandArgs(trailingOnly = TRUE)) {
  parts <- strsplit(argument, "=", fixed = TRUE)[[1]]
  stopifnot(length(parts) == 2)
  fits[[parts[1]]] <- parts[2]
}
stopifnot(all(file.exists(fits)))
reference <- "ref_f0"
fit_colours <- c(ref_f0 = grey(0.55), V3f = "#E69F00", V4 = "#0072B2",
                 V4_class = "#CC79A7", V5 = "#009E73", V5_shear = "#D55E00")
# the weighted binomial fits (#47) in the colours of the fits they copy
fit_colours <- c(fit_colours, wb_ref = grey(0.3), wb_bf = "#56B4E9",
                 wb_v3f = "#E69F00", wb_v4 = "#0072B2", wb_v4_class = "#CC79A7",
                 wb_v5 = "#009E73")

output_dir <- "outputs/species_runs/compare"
figure_dir <- "figures/species_runs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

llin_pyrethroids <- c("Alpha-cypermethrin", "Deltamethrin", "Permethrin")
areas <- c("West", "Central", "East", "Southern", "Horn", "Ethiopia")
cuts <- c(2014, 2018)
window <- 5


# per fit: posterior mean mortality at the bioassays, and the regional trends

# the regional weights of figure 1: one row per area, cell and insecticide
# with LLIN-pyrethroid bioassays, weighted by its share of the area's
# mosquitoes tested
area_weights <- function(df) {
  rows <- df %>%
    filter(insecticide_type %in% llin_pyrethroids) %>%
    mutate(area = analysis_region(country_name, region))
  bind_rows(rows, rows %>% filter(country_name == "Ethiopia") %>%
              mutate(area = "Ethiopia")) %>%
    filter(area %in% areas) %>%
    group_by(area, cell_id, type_id, cell) %>%
    summarise(tested = sum(mosquito_number), .groups = "drop") %>%
    group_by(area) %>%
    mutate(weight = tested / sum(tested)) %>%
    ungroup()
}

summarise_fit <- function(label, file) {
  cache <- file.path(output_dir, sprintf("%s.rds", label))
  if (file.exists(cache) && file.mtime(cache) > file.mtime(file)) {
    return(readRDS(cache))
  }
  fit <- load_fit(file)
  df <- fit$df
  chosen <- even_draws(fit, n_draws, label)
  parameters <- fit_parameter_draws(fit, chosen$index)

  # posterior mean mortality at every bioassay, as the fit predicts it (with
  # the species model, the mixture at the bioassay's arabiensis share)
  time <- system.time(
    p_mean <- colMeans(plogis(dynamical_logit(parameters, df, df,
                                              fit$x_cell_years,
                                              fit$cell_years_index)))
  )[["elapsed"]]
  report("%s: mean predictions at %d bioassays from %d draws in %.0f s",
         label, nrow(df), nrow(parameters$effect_type), time)

  # every year at each weighted cell and insecticide; with the species model
  # the whole complex, at r(x)
  weights <- area_weights(df)
  n_times <- max(fit$cell_years_index$year_id)
  pairs <- distinct(weights, cell_id, type_id, cell)
  rows <- pairs[rep(seq_len(nrow(pairs)), n_times), ]
  rows$year_id <- rep(seq_len(n_times), each = nrow(pairs))
  share <- prediction_share(fit$options, rows$cell)
  time <- system.time(
    p <- plogis(dynamical_logit(parameters, rows, df, fit$x_cell_years,
                                fit$cell_years_index, share = share))
  )[["elapsed"]]
  report("%s: predictions at %d cell-insecticide-years in %.0f s", label,
         nrow(rows), time)

  # draws x (area, year) by a sparse sum over each area's rows
  key <- paste(rows$cell_id, rows$type_id)
  trends <- bind_rows(lapply(areas, function(area) {
    w <- weights[weights$area == area, ]
    w_rows <- w$weight[match(key, paste(w$cell_id, w$type_id))]
    w_rows[is.na(w_rows)] <- 0
    by_year <- vapply(seq_len(n_times), function(t) {
      at <- rows$year_id == t
      c(p[, at, drop = FALSE] %*% w_rows[at])
    }, numeric(nrow(p)))
    tibble(area = area,
           year = fit$baseline_year - 1 + seq_len(n_times),
           mean = colMeans(by_year),
           lower = apply(by_year, 2, quantile, 0.025),
           upper = apply(by_year, 2, quantile, 0.975))
  })) %>%
    mutate(label = label, .before = 1)

  out <- list(label = label, file = file, p_mean = p_mean,
              key = df %>% select(cell, type_id, year_start, died,
                                  mosquito_number),
              trends = trends, n_draws = nrow(parameters$effect_type))
  saveRDS(out, cache)
  report("%s saved; peak memory %.1f GB", label, peak_memory_gb())
  rm(fit, parameters, p)
  invisible(gc())
  out
}

summaries <- lapply(names(fits), function(label) {
  summarise_fit(label, fits[[label]])
})
names(summaries) <- names(fits)
# every fit is of the same bioassays
for (s in summaries[-1]) {
  stopifnot(identical(s$key, summaries[[1]]$key))
}

# the data, as the fits saw them
fit_env <- new.env()
load(fits[[1]], envir = fit_env)
df <- fit_env$df
rm(fit_env)
invisible(gc())


# 1. regional trends ---------------------------------------------------------------

# Figure 1's data points (R/fig_temporal_preds_data.R, df_overall_plot then
# df_pyrethroids_plot), for the records `data` of one area
pyrethroid_points <- function(data) {
  per_insecticide <- data %>%
    filter(insecticide_type %in% llin_pyrethroids) %>%
    group_by(cell_id, insecticide = insecticide_type, year = year_start) %>%
    summarise(died = sum(died), mosquito_number = sum(mosquito_number),
              bioassays = n(), .groups = "drop") %>%
    group_by(insecticide) %>%
    mutate(total_overall_mosquito_number = sum(mosquito_number)) %>%
    group_by(insecticide, cell_id) %>%
    mutate(cell_overall_mosquito_number = sum(mosquito_number)) %>%
    group_by(insecticide, year) %>%
    mutate(total_year_mosquito_number = sum(mosquito_number)) %>%
    group_by(insecticide, year, cell_id) %>%
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

records <- df %>% mutate(area = analysis_region(country_name, region))
points <- bind_rows(lapply(areas, function(area) {
  data <- if (area == "Ethiopia") {
    filter(records, country_name == "Ethiopia")
  } else {
    filter(records, .data$area == !!area)
  }
  pyrethroid_points(data) %>% mutate(area = area, .before = 1)
}))
area_labels <- setNames(
  sprintf("%s) %s", LETTERS[seq_along(areas)], areas), areas)
points$area_label <- factor(area_labels[points$area], levels = area_labels)

trends <- bind_rows(lapply(summaries, `[[`, "trends")) %>%
  mutate(area_label = factor(area_labels[area], levels = area_labels),
         label = factor(label, levels = names(fits)))
compared <- filter(trends, label != reference)

trend_figure <- ggplot(mapping = aes(x = year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = label),
              data = compared, alpha = 0.22, colour = NA) +
  geom_line(aes(y = mean), data = filter(trends, label == reference),
            colour = fit_colours[[reference]], linewidth = 0.35) +
  geom_line(aes(y = mean, colour = label), data = compared,
            linewidth = 0.7) +
  geom_point(aes(y = Susceptibility, size = effective_bioassays),
             data = points, shape = 21, fill = grey(0.85),
             colour = "black", stroke = 0.3) +
  facet_wrap(~ area_label, nrow = 2) +
  scale_colour_manual(values = fit_colours, name = NULL) +
  scale_fill_manual(values = fit_colours, guide = "none") +
  scale_size_area(limits = c(0, 500), max_size = 6,
                  name = "Effective\nsamples") +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(xlim = c(1995, 2024), ylim = c(0, 1)) +
  labs(x = NULL, y = "Susceptibility",
       title = "LLIN pyrethroids: regional predictions and data",
       caption = paste0(
         "Lines and 95% bands: posterior predicted mortality averaged over ",
         "each region's bioassay cells, weighted by mosquitoes tested (as ",
         "figure 1, per region); ", reference, " thin grey.\n",
         "Points: figure 1's weighted annual estimates within the region, ",
         "size the effective number of bioassays.")) +
  theme_minimal() +
  theme(strip.text.x = element_text(hjust = 0),
        legend.position = "right",
        plot.caption = element_text(hjust = 0, size = 8))
ggsave(file.path(figure_dir, "regional_trends_pyrethroids.png"),
       trend_figure, width = 11, height = 6.5, dpi = 200, bg = "white")
write.csv(trends, file.path(output_dir, "regional_trends.csv"),
          row.names = FALSE)
write.csv(points, file.path(output_dir, "regional_points.csv"),
          row.names = FALSE)


# 2. in-sample variance explained in within-pixel change -------------------------

rho_spec <- rho_lookup()
rho_draws <- readRDS("outputs/bioassay_rho_type_draws.rds")
predicted <- as_tibble(lapply(summaries, `[[`, "p_mean"))
model_columns <- names(fits)

# per cell and insecticide with bioassays in both windows of a cut: the
# observed change in pooled mortality, the predicted change of each model,
# and what the floor needs
change_groups <- bind_rows(lapply(cuts, function(cut) {
  before <- seq(cut - window, cut - 1)
  after <- seq(cut, cut + window - 1)
  data <- bind_cols(df, predicted) %>%
    mutate(when = case_when(year_start %in% before ~ "before",
                            year_start %in% after ~ "after")) %>%
    filter(!is.na(when))
  summary <- data %>%
    group_by(cell, insecticide_type, insecticide_class, when) %>%
    summarise(died = sum(died),
              tested = sum(mosquito_number),
              # for noise_floor_var_pooled(): sum(n) and sum(n (n - 1))
              size_sq = sum(mosquito_number * (mosquito_number - 1)),
              assays = n(),
              across(all_of(model_columns),
                     ~ sum(.x * mosquito_number) / sum(mosquito_number)),
              .groups = "drop")
  summary %>%
    pivot_wider(names_from = when,
                values_from = c(died, tested, size_sq, assays,
                                all_of(model_columns))) %>%
    filter(!is.na(tested_before), !is.na(tested_after)) %>%
    mutate(cut = cut,
           observed = died_after / tested_after -
             died_before / tested_before)
})) %>%
  mutate(rho = rho_for_record(., rho_spec))
for (model in model_columns) {
  change_groups[[paste0("p_", model)]] <-
    change_groups[[paste0(model, "_after")]] -
    change_groups[[paste0(model, "_before")]]
}

# noise_floor_var_pooled() of a window, vectorised over groups and rho
window_floor <- function(died, tested, size_sq, rho) {
  inflation <- (tested + rho * size_sq) / tested ^ 2
  yhat <- died / tested
  out <- pmax(yhat * (1 - yhat) / (1 - inflation), 0) * inflation
  out[!is.finite(inflation) | inflation >= 1] <- NA_real_
  out
}
change_floor <- function(groups, rho) {
  window_floor(groups$died_before, groups$tested_before,
               groups$size_sq_before, rho) +
    window_floor(groups$died_after, groups$tested_after,
                 groups$size_sq_after, rho)
}
change_groups$floor_variance <- change_floor(change_groups,
                                             change_groups$rho)
# the vectorised floor is noise_floor_var_pooled()'s, checked on one group
check_group <- bind_cols(df, predicted) %>%
  filter(cell == change_groups$cell[1],
         insecticide_type == change_groups$insecticide_type[1],
         year_start %in% seq(change_groups$cut[1] - window,
                             change_groups$cut[1] - 1))
stopifnot(abs(noise_floor_var_pooled(check_group$died,
                                     check_group$mosquito_number,
                                     change_groups$rho[1]) -
                window_floor(change_groups$died_before[1],
                             change_groups$tested_before[1],
                             change_groups$size_sq_before[1],
                             change_groups$rho[1])) < 1e-12)
# as R/validation_change.R, groups whose floor is undefined are left out
change_groups <- filter(change_groups, !is.na(floor_variance))

prediction_columns <- paste0("p_", model_columns)
explained_for <- function(data) {
  variance <- mean((data$observed - mean(data$observed)) ^ 2)
  vapply(prediction_columns, function(column) {
    100 * (1 - mean((data$observed - data[[column]]) ^ 2) / variance)
  }, numeric(1))
}
set.seed(47)
summarise_change <- function(data, stratum) {
  variance <- mean((data$observed - mean(data$observed)) ^ 2)
  point <- explained_for(data)
  replicates <- pixel_bootstrap(data, explained_for, 2000)
  # the noise share over draws of rho, per type
  index <- match(data$insecticide_type, colnames(rho_draws))
  stopifnot(!anyNA(index))
  noise <- vapply(sample(nrow(rho_draws), 1000, replace = TRUE), function(i) {
    100 * mean(change_floor(data, rho_draws[i, index]), na.rm = TRUE) /
      variance
  }, numeric(1))
  bind_rows(
    tibble(quantity = model_columns, kind = "model",
           estimate = point[prediction_columns],
           lower = apply(replicates[, prediction_columns, drop = FALSE], 2,
                         quantile, 0.025, na.rm = TRUE),
           upper = apply(replicates[, prediction_columns, drop = FALSE], 2,
                         quantile, 0.975, na.rm = TRUE)),
    tibble(quantity = "bioassay variability", kind = "noise",
           estimate = 100 * mean(data$floor_variance) / variance,
           lower = quantile(noise, 0.025), upper = quantile(noise, 0.975))
  ) %>%
    mutate(stratum = stratum, groups = nrow(data),
           pixels = n_distinct(data$cell),
           assays = sum(data$assays_before + data$assays_after),
           variance = variance,
           mean_observed_change = mean(data$observed), .before = 1)
}
strata <- c("all insecticides", "Pyrethroids", "Organochlorines",
            "Organophosphates", "Carbamates")
change_table <- bind_rows(lapply(strata, function(stratum) {
  data <- if (stratum == "all insecticides") change_groups else
    filter(change_groups, insecticide_class == stratum)
  if (n_distinct(data$cell) < 2) return(NULL)
  summarise_change(data, stratum)
}))
write.csv(change_table, "outputs/species_runs/insample_change_variance.csv",
          row.names = FALSE)

# where the error comes from: the mean observed and predicted change, and
# their correlation, per stratum and model
change_bias <- bind_rows(lapply(strata, function(stratum) {
  data <- if (stratum == "all insecticides") change_groups else
    filter(change_groups, insecticide_class == stratum)
  bind_rows(lapply(model_columns, function(model) {
    predicted_change <- data[[paste0("p_", model)]]
    tibble(stratum = stratum, model = model, groups = nrow(data),
           mean_observed = mean(data$observed),
           mean_predicted = mean(predicted_change),
           correlation = cor(data$observed, predicted_change),
           sd_observed = sd(data$observed),
           sd_predicted = sd(predicted_change))
  }))
}))
write.csv(change_bias, "outputs/species_runs/insample_change_bias.csv",
          row.names = FALSE)
print(as.data.frame(change_bias), digits = 3)
options(width = 160)
print(as.data.frame(change_table %>%
                      mutate(value = sprintf("%5.1f [%5.1f, %5.1f]", estimate,
                                             lower, upper)) %>%
                      select(stratum, groups, pixels, quantity, value) %>%
                      pivot_wider(names_from = quantity,
                                  values_from = value)))

# the bars, in the encoding of R/fig_variance_bars.R
laid_models <- change_table %>%
  filter(kind == "model") %>%
  mutate(quantity = factor(quantity, levels = model_columns),
         position = as.integer(quantity),
         experiment = factor(stratum, levels = strata))
laid_noise <- change_table %>%
  filter(kind == "noise") %>%
  select(stratum, estimate, lower, upper) %>%
  right_join(laid_models %>% select(stratum, position, experiment),
             by = "stratum") %>%
  mutate(kind = "noise")
laid <- bind_rows(laid_models, laid_noise)
reference_bar <- laid_models %>% filter(stratum == strata[1], quantity == "V4")
reference_noise <- laid_noise %>% filter(stratum == strata[1]) %>% slice(1)
key <- lapply(region_key(reference_bar$estimate, reference_noise$estimate,
                         length(model_columns) + 0.6, wrap = TRUE),
              function(layer) {
                layer$data$experiment <- factor(strata[1], levels = strata)
                layer
              })
counts <- change_table %>%
  filter(kind == "noise") %>%
  mutate(experiment = factor(stratum, levels = strata),
         text = sprintf("%d pixel-insecticides\n%d pixels", groups, pixels))
change_figure <- ggplot() +
  bar_layers(laid) +
  key +
  geom_text(aes(x = (length(model_columns) + 1) / 2, y = -9, label = text),
            data = counts, size = 2.4, colour = grey(0.35)) +
  facet_wrap(~ experiment, nrow = 1) +
  scale_fill_manual(values = fit_colours, breaks = model_columns) +
  scale_x_continuous(breaks = seq_along(model_columns),
                     labels = model_columns,
                     expand = expansion(add = 0)) +
  scale_y_continuous(breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0.6, length(model_columns) + 0.4),
                  ylim = c(-14, 100), clip = "off") +
  labs(x = NULL, y = "% of variance in within-pixel change, in sample",
       caption = paste(
         "Change in pooled mortality between the five years before and from",
         "a cut (2014 and 2018, pooled), per pixel and insecticide with",
         "bioassays in both. In sample: every fit has seen both windows.")) +
  base_theme +
  theme(panel.spacing.x = unit(26, "pt"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        plot.caption = element_text(hjust = 0, size = 8),
        plot.margin = margin(6, 90, 6, 6))
ggsave(file.path(figure_dir, "insample_change_variance.png"), change_figure,
       width = 12, height = 4.8, dpi = 200, bg = "white")
report("written %s and %s", file.path(figure_dir,
                                      "regional_trends_pyrethroids.png"),
       file.path(figure_dir, "insample_change_variance.png"))
