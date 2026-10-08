# Do the latent smooths of V5 (#47; R/latent_smooth.R) trade off with its
# selection effects? Across posterior draws, the correlations of the
# selection coefficients of the pyrethroids (beta_overall, beta_class of the
# pyrethroids and beta_type of the four pyrethroid types, for every
# covariate) with the smooths' amplitudes and ranges and the floor
# intercepts:
#   log_sd_selection, log_sd_floor        log marginal sd of u_s and u_f
#   log_range_selection, log_range_floor  log range (km)
#   log_M_cells   log of the mean of exp(u_s(x)) over the modelled bioassay
#                 cells (each once)
#   log_M_pop     log of the mean of exp(u_s(x)) over the limits of
#                 transmission without water bodies
#                 (data/clean/pfpr_water_mask.tif), weighted by population in
#                 2020 (data/clean/pop_cube.tif), on the grid aggregated by 3
#                 (about 14 km)
#   floor_<class> the logit floor intercept of each class
# and of the post-hoc LLIN-use pressure on each pyrethroid type,
# log(exp(beta_type) mean(x_nets exp(u_s(x)))) over the modelled cell-years
# (R/fig_covariate_effects.R). exp(u_s) multiplies the cumulative log fitness
# of every class at a cell, so for weak selection a larger M can stand in
# for larger selection effects everywhere: a negative correlation between
# beta and log M is that trade-off, and it is confounding only if the
# post-hoc pressure moves with M too.
#
#   Rscript R/field_fixed_correlations.R <file>=<label> ...
#
# per fit, all its usable chains (USE_CHAINS, e.g. "V5=1,2"), with
# correlations pooled over them and within each. Writes
#   outputs/species_runs/grids/field_fixed_correlations.csv  every pair
#   outputs/species_runs/grids/field_fixed_draws.rds         the draws
#   figures/species_runs/grids/field_fixed_correlations.png  the key pairs
# Plain R; about 4 GB and 5 minutes for two fits.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) >= 1, all(grepl("=", arguments)))
files <- setNames(sub("=.*$", "", arguments), sub("^.*=", "", arguments))
n_draws <- 4000
chunk <- 4000

suppressMessages({
  library(greta)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(terra)
})
source("R/dynamical_predictions.R")
source("R/species_fit_helpers.R")

output_dir <- "outputs/species_runs/grids"
figure_dir <- "figures/species_runs/grids"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# the population in 2020 in the limits of transmission without water bodies,
# summed over blocks of 3 x 3 cells, at the blocks' centres with any
mask_population <- function() {
  water_mask <- rast("data/clean/pfpr_water_mask.tif")
  population <- rast("data/clean/pop_cube.tif")[["pop_2020"]]
  # the grids differ by 1e-12 degrees in origin
  ext(population) <- ext(water_mask)
  population <- mask(population, water_mask)
  blocks <- aggregate(population, 3, fun = "sum", na.rm = TRUE)
  cells <- cells(blocks)
  weight <- blocks[cells][[1]]
  keep <- !is.na(weight) & weight > 0
  list(xy = xyFromCell(blocks, cells[keep]), weight = weight[keep])
}
population <- mask_population()
report("population weights: %d blocks, %.0f million people",
       length(population$weight), sum(population$weight) / 1e6)

# the mean of exp(u) per draw over points `coords` weighted by `weight`, for
# smooth weights `w` (draws x basis functions), `chunk` points at a time
weighted_mean_exp <- function(w, smooth, coords, weight) {
  total <- numeric(nrow(w))
  for (rows in split(seq_len(nrow(coords)),
                     ceiling(seq_len(nrow(coords)) / chunk))) {
    basis <- smooth_basis_at(smooth, coords[rows, , drop = FALSE])
    total <- total + c(exp(w %*% t(basis)) %*% weight[rows])
  }
  total / sum(weight)
}

draws <- list()
for (label in names(files)) {
  fit <- load_fit(files[[label]])
  stopifnot(smooth_on(fit$options))
  used <- even_draws(fit, n_draws, label)
  parameters <- fit_parameter_draws(fit, used$index)
  v <- parameters$variables
  terms <- dynamical_terms_draws(v, fit$classes_index, fit$types,
                                 terms = c("beta_type", "beta_class"),
                                 options = fit$options)
  covariates <- colnames(fit$x_cell_years)
  short <- str_remove(str_remove(covariates, ":g_ag$"), ":g_dom$")
  short <- str_replace_all(short, " ", "_")
  pyrethroid_types <- fit$types[fit$classes[fit$classes_index] ==
                                  "Pyrethroids"]
  pyrethroid_class <- match("Pyrethroids", fit$classes)

  out <- tibble(fit = label, chain = used$chain, draw = seq_along(used$chain))
  for (k in seq_along(covariates)) {
    out[[paste0("beta_overall_", short[k])]] <- v$beta_overall[, k]
    out[[paste0("beta_class_", short[k])]] <-
      terms$beta_class[, k, pyrethroid_class]
    for (type in pyrethroid_types) {
      out[[paste0("beta_type_", short[k], "_", type)]] <-
        terms$beta_type[, k, match(type, fit$types)]
    }
  }

  # the smooths' hyperparameters and floor intercepts
  for (kind in smooth_kinds(fit$options)) {
    names <- smooth_variable_names(kind)
    out[[paste0("log_sd_", kind)]] <- log(c(v[[names[["sd"]]]]))
    out[[paste0("log_range_", kind)]] <-
      log(1000 / c(v[[names[["inv_range"]]]]))
  }
  intercepts <- matrix(v$floor_intercept, nrow(out))
  for (class in seq_len(ncol(intercepts))) {
    out[[paste0("floor_", fit$classes[class])]] <- intercepts[, class]
  }

  # M at the modelled cells and over the population, and the post-hoc
  # LLIN-use pressure on each pyrethroid type
  smooth <- fit$options$smooth
  w <- parameters$smooth_weights$selection
  n_cells <- max(fit$df$cell_id)
  cells <- fit$df$cell[match(seq_len(n_cells), fit$df$cell_id)]
  u_cells <- w %*% t(prediction_basis(fit$options, cells))
  out$log_M_cells <- log(rowMeans(exp(u_cells)))
  # the realised spread of each smooth over the modelled cells, against its
  # marginal sd
  out$log_spread_selection <- log(apply(u_cells, 1, sd))
  if ("floor" %in% smooth_kinds(fit$options)) {
    out$log_spread_floor <- log(apply(
      parameters$smooth_weights$floor %*%
        t(prediction_basis(fit$options, cells)), 1, sd))
  }
  out$log_M_pop <- log(weighted_mean_exp(
    w, smooth, smooth_coords(population$xy[, 1], population$xy[, 2],
                             smooth$crs),
    population$weight))
  x_nets <- rowsum(fit$x_cell_years[, "nets"], fit$cell_years_index$cell_id,
                   reorder = TRUE)[, 1] / nrow(fit$x_cell_years)
  nets_weighted <- c(exp(u_cells) %*% x_nets)
  for (type in pyrethroid_types) {
    out[[paste0("log_posthoc_nets_", type)]] <-
      terms$beta_type[, match("nets", covariates), match(type, fit$types)] +
      log(nets_weighted)
  }
  report("%s: %d draws of chains %s; median M over cells %.3g, over the population %.3g",
         label, nrow(out), toString(sort(unique(out$chain))),
         median(exp(out$log_M_cells)), median(exp(out$log_M_pop)))
  print(out %>%
          group_by(chain) %>%
          summarise(across(c(log_sd_selection, log_spread_selection,
                             any_of(c("log_sd_floor", "log_spread_floor")),
                             log_M_cells, log_M_pop),
                           ~ sprintf("%.2f (%.2f-%.2f)", median(exp(.x)),
                                     quantile(exp(.x), 0.05),
                                     quantile(exp(.x), 0.95)))) %>%
          rename_with(~ sub("^log_", "", .x)) %>%
          as.data.frame(), row.names = FALSE)
  draws[[label]] <- out
  rm(fit, parameters, v, terms, u_cells)
  gc()
}
draws <- bind_rows(draws)
saveRDS(draws, file.path(output_dir, "field_fixed_draws.rds"))


# correlations ------------------------------------------------------------------

field_names <- c(grep("^log_(sd|range|M|spread)_", names(draws),
                      value = TRUE),
                 grep("^floor_", names(draws), value = TRUE))
beta_names <- c(grep("^beta_", names(draws), value = TRUE),
                grep("^log_posthoc_", names(draws), value = TRUE))
# pooled over each fit's chains, and within each chain
groups <- bind_rows(
  draws %>% mutate(group = fit),
  draws %>% mutate(group = paste0(fit, " ch", chain)))
correlations <- groups %>%
  group_by(fit, group) %>%
  group_modify(function(d, key) {
    r <- cor(as.matrix(d[, beta_names]), as.matrix(d[, field_names]))
    as_tibble(r, rownames = "coefficient") %>%
      pivot_longer(-coefficient, names_to = "field", values_to = "r")
  }) %>%
  ungroup()
write.csv(correlations, file.path(output_dir, "field_fixed_correlations.csv"),
          row.names = FALSE)

options(width = 160)
cat("\nthe strongest correlations of each group (|r|, coefficients against field quantities):\n")
correlations %>%
  filter(!is.na(r)) %>%
  group_by(group) %>%
  slice_max(abs(r), n = 12) %>%
  mutate(r = round(r, 2)) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\nLLIN use on the pyrethroids against the field amplitude and ranges, by group:\n")
correlations %>%
  filter(str_detect(coefficient, "_nets"),
         str_detect(coefficient, "^(beta_overall|beta_class)|Deltamethrin|Permethrin|posthoc"),
         field %in% c("log_sd_selection", "log_M_cells", "log_M_pop",
                      "log_range_selection", "log_range_floor",
                      "floor_Pyrethroids")) %>%
  mutate(r = round(r, 2)) %>%
  pivot_wider(names_from = field, values_from = r) %>%
  select(-fit) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\nthe field quantities with each other (pooled and within chains):\n")
for (g in unique(groups$group)) {
  d <- filter(groups, group == g)
  cat("\n", g, "\n")
  print(round(cor(as.matrix(d[, field_names])), 2))
}

# the trade-off of nets on the pyrethroids with log M: the sd of beta_type,
# of beta_type + log M (the common shift M stands for) and of the post-hoc
# pressure, and the correlation of beta_type with log M
cat("\nnets on the pyrethroid types: posterior sd of beta, of beta + log M_cells and of the log post-hoc pressure, and cor(beta, log M_cells):\n")
pyrethroid_types <- sub("^beta_type_nets_", "",
                        grep("^beta_type_nets_", names(groups), value = TRUE))
bind_rows(lapply(pyrethroid_types, function(type) {
  groups %>%
    transmute(group, type = type,
              beta = .data[[paste0("beta_type_nets_", type)]],
              posthoc = .data[[paste0("log_posthoc_nets_", type)]],
              log_M_cells)
})) %>%
  group_by(group, type) %>%
  summarise(median_beta = median(beta),
            sd_beta = sd(beta),
            sd_beta_plus_log_M = sd(beta + log_M_cells),
            sd_log_posthoc = sd(posthoc),
            cor_beta_log_M = cor(beta, log_M_cells),
            median_log_posthoc = median(posthoc),
            .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  as.data.frame() %>%
  print(row.names = FALSE)


# the figure --------------------------------------------------------------------

# the key pairs: LLIN use on the pyrethroids (overall, class, the mean over
# the four types, and the post-hoc pressure on deltamethrin) against the
# selection smooth's sd and M and both ranges, 600 draws per chain
effect_rows <- c(
  "beta_overall_nets" = "beta_overall\nnets",
  "beta_class_nets" = "beta_class nets\npyrethroids",
  "beta_type_nets_mean" = "beta_type nets\nmean of 4\npyrethroids",
  "log_posthoc_nets_Deltamethrin" = "log post-hoc\nLLIN pressure\ndeltamethrin")
field_columns <- c(
  "log_sd_selection" = "sd u_s",
  "log_M_cells" = "M (cells)",
  "log_M_pop" = "M (population)",
  "log_range_selection" = "range u_s (km)",
  "log_range_floor" = "range u_f (km)")
plot_draws <- draws %>%
  mutate(beta_type_nets_mean = rowMeans(across(matches("^beta_type_nets_"))),
         series = paste0(fit, " ch", chain)) %>%
  group_by(series) %>%
  slice(round(seq(1, n(), length.out = min(n(), 600)))) %>%
  ungroup() %>%
  select(fit, series, all_of(names(effect_rows)), all_of(names(field_columns))) %>%
  pivot_longer(all_of(names(effect_rows)), names_to = "effect",
               values_to = "effect_value") %>%
  pivot_longer(all_of(names(field_columns)), names_to = "field",
               values_to = "field_value") %>%
  mutate(effect = factor(effect_rows[effect], levels = effect_rows),
         field = factor(field_columns[field], levels = field_columns))
# the correlations are of the logs; the axes are in natural units, and the
# draws in random order, so that no chain hides another
set.seed(1)
plot_draws <- plot_draws[sample(nrow(plot_draws)), ]
# the correlation in each panel: within each chain, averaged, by fit
panel_r <- plot_draws %>%
  group_by(fit, series, effect, field) %>%
  summarise(r = cor(effect_value, field_value), .groups = "drop") %>%
  group_by(fit, effect, field) %>%
  summarise(r = mean(r), .groups = "drop") %>%
  group_by(effect, field) %>%
  summarise(text = paste(sprintf("%s %.2f", fit, r), collapse = "\n"),
            .groups = "drop")
series_levels <- sort(unique(plot_draws$series))
series_colours <- setNames(
  c("#2a78d6", "#7fb2ef", "#eb6834", "#1baf7a", "#4a3aa7", "#e87ba4",
    "#eda100", "#008300")[seq_along(series_levels)],
  series_levels)
correlation_plot <- ggplot(plot_draws, aes(x = exp(field_value),
                                           y = effect_value,
                                           colour = series)) +
  geom_point(size = 0.35, alpha = 0.35) +
  geom_text(data = panel_r, aes(x = -Inf, y = Inf, label = text),
            inherit.aes = FALSE, hjust = -0.05, vjust = 1.1, size = 2.4,
            lineheight = 0.9) +
  facet_grid(effect ~ field, scales = "free") +
  scale_colour_manual(values = series_colours, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        panel.spacing.x = unit(14, "pt"),
        strip.text.y = element_text(angle = 0, hjust = 0)) +
  labs(x = NULL, y = NULL,
       title = "LLIN-use selection on the pyrethroids against the latent smooths",
       subtitle = paste("Posterior draws (600 per chain); text: correlation",
                        "of the logs within chains, averaged over each fit's",
                        "chains. V5: chains 1-2; wb_v5 not converged"))
ggsave(file.path(figure_dir, "field_fixed_correlations.png"),
       correlation_plot, bg = "white", width = 11, height = 8.5)
report("wrote %s", file.path(figure_dir, "field_fixed_correlations.png"))
