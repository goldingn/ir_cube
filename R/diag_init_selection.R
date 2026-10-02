# Diagnostics of the trade-off between the initial-state population effect
# (init_coef[pop_2000, type]) and the trended population selection effect
# (beta_type[pop_enc:g_dom, type]) in the full fit:
#   3. posterior correlation of the two, and of init_coef with the net effect
#      of population on logit mortality
#   4. fitted trends against the bioassays, at selected pixels and pooled by
#      population density class
#   5. the logit initial state and the cumulative logit change from each
#      selection term at the bioassay pixels, by density class
#
#   Rscript R/diag_init_selection_draws.R   # once, to extract the draws
#   Rscript R/diag_init_selection.R
#
# Predictions use the shared plain-R functions (map_covariates(),
# map_type_logit(), R/two_stage_map_functions.R) with the fit's parameter
# draws from dynamical_parameter_draws(). Writes figures/diag_*.png and
# outputs/diag_init_selection_correlations.csv.

suppressMessages({
  library(tidyverse)
  library(terra)
  library(patchwork)
})
source("R/functions.R")
source("R/model_covariates.R")
source("R/dynamical_predictions.R")
source("R/two_stage_map_functions.R")

fit <- readRDS("outputs/diag_init_selection_draws.RDS")
par <- fit$parameters
types <- fit$types
design <- fit$options$selection_columns
baseline_year <- fit$baseline_year
end_year <- 2030
years <- baseline_year:end_year
n_draws_all <- par$n_draws
pop_col <- match("pop_enc:g_dom", par$covariate_names)
type_colours <- insecticide_colours()

# 500 of the 2000 draws for the predictions
pred_draws <- round(seq(1, n_draws_all, length.out = 500))

# population density in 2000 at the bioassay pixels, and its classes
mask <- rast("data/clean/raster_mask.tif")
area <- terra::cellSize(mask, unit = "km")
pop_2000 <- pop_zero_filled(rast("data/clean/pop_cube.tif")[["pop_2000"]])
density_cells <- terra::extract(pop_2000 / area, fit$unique_cells)[, 1]
density_class_of <- function(d) {
  cut(d, c(-Inf, 10, 100, 1000, Inf),
      labels = c("<10", "10-100", "100-1,000", ">1,000"), right = FALSE)
}
cell_class <- density_class_of(density_cells)
density_colours <- setNames(c("#c6dbef", "#6baed6", "#2171b5", "#08306b"),
                            levels(cell_class))

# standardised pop_2000 at encounter 0 (empty) and 1 (saturated)
pop_layer <- terra::mask(init_pop_layer(complete_selection_design(design)),
                         mask)
pop_moments <- terra::global(pop_layer, c("mean", "sd"), na.rm = TRUE)
z_empty <- (0 - pop_moments$mean) / pop_moments$sd
z_full <- (1 - pop_moments$mean) / pop_moments$sd
rm(pop_layer, pop_2000)
invisible(gc())

df <- fit$df %>%
  mutate(density_class = cell_class[cell_id],
         insecticide_type = factor(insecticide_type,
                                   levels = insecticides_plot_order))


# figure 3: posterior correlations ------------------------------------------

# net effect of population on logit mortality in 2020: a saturated pixel
# (encounter 1) against an empty one (encounter 0), other selection
# covariates 0 and the same country: the initial-state difference less the
# cumulative selection from pop_enc x g_dom over 1995-2020
g_dom <- selection_trend_matrix(1, baseline_year, 2020, design$trend_pop,
                                design)[1, ]
correlations <- map_dfr(seq_along(types), function(k) {
  coef <- par$init_coef[, 1, k]
  beta <- log(par$effect_type[, pop_col, k])
  init_diff <- coef * (z_full - z_empty)
  selection_2020 <- rowSums(log1p(outer(exp(beta), g_dom)))
  tibble(draw = seq_len(n_draws_all),
         type = types[k],
         init_coef_pop = coef,
         beta_pop = beta,
         init_diff = init_diff,
         selection_2020 = selection_2020,
         net_2020 = init_diff - selection_2020)
})

correlation_table <- correlations %>%
  group_by(type) %>%
  summarise(
    cor_coef_beta = cor(init_coef_pop, beta_pop),
    cor_coef_net = cor(init_coef_pop, net_2020),
    cor_beta_net = cor(beta_pop, net_2020),
    init_coef_pop_mean = mean(init_coef_pop),
    beta_pop_mean = mean(beta_pop),
    init_diff_mean = mean(init_diff),
    selection_2020_mean = mean(selection_2020),
    net_2020_mean = mean(net_2020),
    net_2020_lower = quantile(net_2020, 0.025),
    net_2020_upper = quantile(net_2020, 0.975),
    .groups = "drop") %>%
  mutate(type = factor(type, levels = insecticides_plot_order)) %>%
  arrange(type)
write_csv(correlation_table, "outputs/diag_init_selection_correlations.csv")
print(correlation_table, width = Inf)

key_types <- c("Deltamethrin", "Permethrin", "DDT", "Bendiocarb")
scatter <- correlations %>%
  filter(type %in% key_types) %>%
  mutate(type = factor(type, levels = key_types))
cor_labels <- correlation_table %>%
  filter(type %in% key_types) %>%
  mutate(type = factor(type, levels = key_types))

p3a <- ggplot(scatter, aes(init_coef_pop, beta_pop, colour = type)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_text(aes(x = -Inf, y = Inf,
                label = sprintf("r = %.2f", cor_coef_beta)),
            data = cor_labels, hjust = -0.1, vjust = 1.3, size = 3,
            colour = "black") +
  facet_wrap(~type, scales = "free", nrow = 1) +
  scale_colour_manual(values = type_colours, guide = "none") +
  labs(x = "init_coef[pop_2000] (logit per sd)",
       y = "beta_type[pop_enc:g_dom]\n(log selection coefficient)") +
  theme_minimal(base_size = 9)

p3b <- ggplot(scatter, aes(init_coef_pop, net_2020, colour = type)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_text(aes(x = -Inf, y = Inf,
                label = sprintf("r = %.2f", cor_coef_net)),
            data = cor_labels, hjust = -0.1, vjust = 1.3, size = 3,
            colour = "black") +
  facet_wrap(~type, scales = "free", nrow = 1) +
  scale_colour_manual(values = type_colours, guide = "none") +
  labs(x = "init_coef[pop_2000] (logit per sd)",
       y = "Net population effect on\nlogit mortality in 2020") +
  theme_minimal(base_size = 9)

p3 <- p3a / p3b + plot_annotation(caption = paste0(
  "Points are 2,000 posterior draws; r is their correlation. Net effect: ",
  "logit mortality at a saturated pixel\n(encounter 1, ",
  sprintf("%.2f", z_full), " sd) less an empty one (",
  sprintf("%.2f", z_empty), " sd) in 2020, ",
  "same country, other selection covariates 0:\ninit_coef x ",
  sprintf("%.2f", z_full - z_empty),
  " less the cumulative log(1 + g_dom exp(beta)) over 1995-2020."))
ggsave("figures/diag_init_selection_correlation.png", p3, width = 9,
       height = 6, dpi = 150, bg = "white")
rm(correlations, scatter)


# predictions at the bioassay pixels -----------------------------------------

covariates <- map_covariates(fit$unique_cells, baseline_year, end_year,
                             design)
stopifnot(isTRUE(all.equal(unname(covariates$init[, fit$options$init_covariates]),
                           unname(fit$x_cells_init[, fit$options$init_covariates]),
                           tolerance = 1e-6)))

# logit relative initial state per country, with the coefficients attached,
# as map_logit_init() returns it
logit_init <- par$logit_init_relative[pred_draws, , , drop = FALSE]
dimnames(logit_init) <- list(NULL, fit$countries, types)
attr(logit_init, "init_coef") <- array(
  par$init_coef[pred_draws, , , drop = FALSE],
  c(length(pred_draws), length(fit$options$init_covariates), length(types)),
  dimnames = list(NULL, fit$options$init_covariates, types))
effect <- par$effect_type[pred_draws, , , drop = FALSE]
floor <- par$mortality_floor[pred_draws]
kappa_type <- par$kappa_type[pred_draws, , drop = FALSE]

# each cell's country, from its first record as dynamical_predictions() takes
# it (a few cells have records in two countries)
cell_country <- fit$df %>%
  distinct(cell_id, .keep_all = TRUE) %>%
  arrange(cell_id)
cell_country <- cell_country$country_id[match(seq_along(fit$unique_cells),
                                              cell_country$cell_id)]

# mortality draws (cells x draws per year) of type k at cells `rows`
predict_cells <- function(k, rows) {
  logit <- map_type_logit(k, rows, cell_country[rows], effect, logit_init,
                          covariates, years, years, floor, kappa_type)
  lapply(logit, plogis)
}

# for each type: the posterior mean at every cell with its bioassays and
# year, and the posterior mean and interval of the mean over each density
# class's cells (those with bioassays of the type at any time)
cell_means <- list()
class_trends <- list()
for (k in seq_along(types)) {
  rows <- sort(unique(fit$df$cell_id[fit$df$type_id == k]))
  p <- predict_cells(k, rows)
  cell_means[[k]] <- tibble(
    type = types[k],
    cell_id = rep(rows, length(years)),
    year = rep(years, each = length(rows)),
    predicted = unlist(lapply(p, rowMeans)))
  classes_k <- cell_class[rows]
  class_trends[[k]] <- map_dfr(levels(cell_class), function(cl) {
    in_class <- classes_k == cl
    if (!any(in_class)) return(NULL)
    map_dfr(seq_along(years), function(j) {
      m <- colMeans(p[[j]][in_class, , drop = FALSE])
      tibble(type = types[k], density_class = cl, year = years[j],
             n_cells = sum(in_class), mean = mean(m),
             lower = quantile(m, 0.025), upper = quantile(m, 0.975))
    })
  })
  rm(p)
  invisible(gc())
}
cell_means <- bind_rows(cell_means)
class_trends <- bind_rows(class_trends) %>%
  mutate(density_class = factor(density_class, levels = levels(cell_class)),
         type = factor(type, levels = insecticides_plot_order))


# figure 4a: selected pixels -------------------------------------------------

trend_types <- c("Deltamethrin", "DDT")
# two pixels per density class: those with the most distinct years of data
# over the two types together
pixels <- df %>%
  filter(insecticide_type %in% trend_types) %>%
  group_by(cell_id, density_class, country_name) %>%
  summarise(n_years = n_distinct(year_start),
            first_year = min(year_start), .groups = "drop") %>%
  group_by(density_class) %>%
  slice_max(n_years, n = 2, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(label = sprintf("%s, %s/km² (%s)", country_name,
                         formatC(density_cells[cell_id], format = "d",
                                 big.mark = ","),
                         density_class))
print(pixels)

pixel_draws <- map_dfr(match(trend_types, types), function(k) {
  p <- predict_cells(k, pixels$cell_id)
  map_dfr(seq_along(years), function(j) {
    tibble(type = types[k], cell_id = pixels$cell_id, year = years[j],
           mean = rowMeans(p[[j]]),
           lower = apply(p[[j]], 1, quantile, 0.025),
           upper = apply(p[[j]], 1, quantile, 0.975))
  })
})
pixel_draws <- left_join(pixel_draws, select(pixels, cell_id, label),
                         by = "cell_id")
pixel_data <- df %>%
  filter(insecticide_type %in% trend_types, cell_id %in% pixels$cell_id) %>%
  transmute(type = as.character(insecticide_type), cell_id, year = year_start,
            mortality = died / mosquito_number, mosquito_number) %>%
  left_join(select(pixels, cell_id, label), by = "cell_id")
pixel_levels <- pixels$label
pixel_draws$label <- factor(pixel_draws$label, levels = pixel_levels)
pixel_data$label <- factor(pixel_data$label, levels = pixel_levels)

p4a <- ggplot(pixel_draws, aes(year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = type), alpha = 0.25) +
  geom_line(aes(y = mean, colour = type)) +
  geom_point(aes(y = mortality, size = mosquito_number, colour = type),
             data = pixel_data, alpha = 0.7, shape = 16) +
  facet_wrap(~label, ncol = 4) +
  scale_colour_manual(values = type_colours, name = "Insecticide") +
  scale_fill_manual(values = type_colours, name = "Insecticide") +
  scale_size_area(max_size = 3, name = "Mosquitoes\ntested") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Year", y = "Bioassay mortality",
       caption = paste(
         "Lines: posterior mean predicted bioassay mortality at the pixel",
         "(500 draws); bands: 95% credible interval.",
         "Points: observed bioassays, sized by mosquitoes tested.\n",
         "Panels: the two pixels per 2000 density class with the most",
         "distinct years of deltamethrin and DDT data; title gives the",
         "country and 2000 density.")) +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(size = 7.5))
ggsave("figures/diag_fitted_trends_pixels.png", p4a, width = 10, height = 5.5,
       dpi = 150, bg = "white")


# figure 4b: pooled by density class ----------------------------------------

observed_class <- df %>%
  group_by(insecticide_type, density_class, year = year_start) %>%
  summarise(mortality = sum(died) / sum(mosquito_number),
            n = n(), .groups = "drop") %>%
  rename(type = insecticide_type)
# the predicted means at the same bioassays (same pixels and years)
predicted_class <- df %>%
  transmute(type = as.character(insecticide_type), cell_id,
            year = year_start, density_class, mosquito_number) %>%
  left_join(cell_means, by = c("type", "cell_id", "year")) %>%
  group_by(type, density_class, year) %>%
  summarise(predicted = weighted.mean(predicted, mosquito_number),
            .groups = "drop") %>%
  mutate(type = factor(type, levels = insecticides_plot_order))

p4b <- ggplot(class_trends, aes(year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = density_class),
              alpha = 0.25) +
  geom_line(aes(y = mean, colour = density_class)) +
  geom_point(aes(y = mortality, size = n, colour = density_class),
             data = observed_class, alpha = 0.7, shape = 16) +
  geom_point(aes(y = predicted, colour = density_class),
             data = predicted_class, shape = 4, size = 0.8) +
  facet_grid(type ~ density_class) +
  scale_colour_manual(values = density_colours,
                      name = "2000 density\n(people per km²)") +
  scale_fill_manual(values = density_colours,
                    name = "2000 density\n(people per km²)") +
  scale_size_area(max_size = 3, name = "Bioassays") +
  scale_y_continuous(labels = scales::percent, breaks = c(0, 0.5, 1)) +
  labs(x = "Year", y = "Bioassay mortality",
       caption = paste(
         "Lines and bands: posterior mean and 95% interval of predicted",
         "mortality averaged over the pixels with bioassays of the type in",
         "the density class (fixed set of pixels).\n",
         "Dots: observed pooled mortality (died / tested) in the year,",
         "sized by bioassays. Crosses: posterior mean prediction at the same",
         "bioassays (weighted by mosquitoes tested).")) +
  theme_minimal(base_size = 8) +
  theme(strip.text.y = element_text(size = 6.5))
ggsave("figures/diag_fitted_trends_density_class.png", p4b, width = 9,
       height = 13, dpi = 150, bg = "white")

# early data by class, for the report
early <- df %>%
  filter(year_start < 2005) %>%
  count(insecticide_type, density_class) %>%
  pivot_wider(names_from = density_class, values_from = n, values_fill = 0)
print(early)


# figure 5: decomposition at the deltamethrin bioassay pixels ----------------

decomp_type <- "Deltamethrin"
k <- match(decomp_type, types)
rows <- sort(unique(fit$df$cell_id[fit$df$type_id == k]))
n_d <- length(pred_draws)
effect_k <- matrix(effect[, , k], nrow = n_d)
term_of <- c("Nets", "IRS", "Population x g_dom",
             rep("Crops x g_ag", length(par$covariate_names) - 3))
terms <- c("Nets", "IRS", "Population x g_dom", "Crops x g_ag", "Reversion")
logit_q0 <- map_cell_logit_init(logit_init, cell_country[rows], k,
                                covariates$init[rows, , drop = FALSE])

# log w_t = log1p(sum_j x_j e_j) split over the terms in proportion to
# x_j e_j; the cumulative logit change in mortality from a term is minus its
# cumulative share, and from reversion - t kappa (>= 0)
cumulative <- setNames(lapply(terms, function(x) {
  matrix(0, length(rows), n_d)
}), terms)
decomp <- list()
for (j in seq_along(years)) {
  x_t <- covariates$time_varying[rows, j, ]
  parts <- lapply(seq_along(term_of), function(i) {
    outer(x_t[, i], effect_k[, i])
  })
  total <- Reduce(`+`, parts)
  log_w <- log1p(total)
  share <- ifelse(total > 0, log_w / total, 0)
  for (term in terms[1:4]) {
    part <- Reduce(`+`, parts[term_of == term])
    cumulative[[term]] <- cumulative[[term]] - part * share
  }
  cumulative$Reversion <- matrix(-kappa_type[, k] * j, length(rows), n_d,
                                 byrow = TRUE)
  for (name in terms) {
    decomp[[length(decomp) + 1]] <- tibble(
      year = years[j], term = name, density_class = cell_class[rows],
      value = rowMeans(cumulative[[name]]))
  }
}
decomp <- bind_rows(decomp) %>%
  group_by(year, term, density_class) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  mutate(term = factor(term, levels = terms))
init_class <- tibble(density_class = cell_class[rows],
                     logit_q0 = rowMeans(logit_q0)) %>%
  group_by(density_class) %>%
  summarise(logit_q0 = mean(logit_q0), n_cells = n(), .groups = "drop")
total_class <- decomp %>%
  group_by(year, density_class) %>%
  summarise(change = sum(value), .groups = "drop") %>%
  left_join(init_class, by = "density_class") %>%
  mutate(logit_q = logit_q0 + change)
print(init_class)
print(filter(decomp, year %in% c(2005, 2015, 2025)) %>%
        pivot_wider(names_from = term, values_from = value), width = Inf)

term_colours <- c("Nets" = "#1b9e77", "IRS" = "#d95f02",
                  "Population x g_dom" = "#7570b3",
                  "Crops x g_ag" = "#e6ab02", "Reversion" = "#666666")
class_labels <- with(init_class, setNames(
  sprintf("%s per km² (%i pixels)", density_class, n_cells), density_class))

p5 <- ggplot(decomp, aes(year)) +
  geom_hline(yintercept = 0, colour = grey(0.7)) +
  geom_line(aes(y = value, colour = term), linewidth = 0.7) +
  geom_line(aes(y = logit_q, linetype = "Logit q_t (initial + all terms)"),
            data = total_class) +
  geom_hline(aes(yintercept = logit_q0, linetype = "Initial logit q_0"),
             data = init_class) +
  facet_wrap(~density_class, nrow = 1,
             labeller = labeller(density_class = class_labels)) +
  scale_colour_manual(values = term_colours,
                      name = "Cumulative logit change from") +
  scale_linetype_manual(values = c("Initial logit q_0" = "dashed",
                                   "Logit q_t (initial + all terms)" =
                                     "solid"),
                        name = NULL) +
  labs(x = "Year", y = "Logit fraction susceptible",
       caption = paste(
         "Deltamethrin, averaged over the pixels with deltamethrin bioassays",
         "in each 2000 density class (posterior means, 500 draws).",
         "Coloured lines: cumulative change in logit q from each selection",
         "term,\nlog w_t = log(1 + sum_j x_j exp(beta_j)) split over the",
         "terms in proportion to x_j exp(beta_j); reversion adds -t kappa.",
         "Dashed: initial logit q_0; solid black: logit q_t,",
         "before the mortality floor.")) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", legend.box = "vertical")
ggsave("figures/diag_selection_decomposition.png", p5, width = 10,
       height = 5, dpi = 150, bg = "white")
cat("done\n")
