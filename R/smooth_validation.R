# Post-hoc checks of the latent smooths of a saved full fit (V5, #47;
# R/latent_smooth.R) against the marker layers and the kdr records.
#
#   Rscript R/smooth_validation.R <fitted_model.RData> <label>
#
# From about n_draws (500) posterior draws, evenly spaced in each usable
# chain, the posterior mean of each smooth (u_s, the log multiplier of
# selection; u_f, the shift of the logit floor; with the shear, also v_s,
# the own part of u_s = v_s + b u_f, and b, reported; with the smooth of the
# initial state, u_init, at sd 1, before the types' loadings):
#   layers   at each modelled bioassay cell (each once), regressed on the
#            marker layers, each standardised over the cells: logit total kdr
#            (995F + 995S) in 2015 of the whole complex
#            (data/clean/kdr_total_2015.tif, band "complex") and logit
#            arabiensis fraction r(x) (data/clean/arabiensis_fraction.tif).
#            R^2 of each alone and of both, and the partial correlation of
#            each with the smooth given the other. The kdr work (#46,
#            ../ir_cube_kdr) saves no maps of 995F and 995S apart, so total
#            kdr only
#   records  at each sample of the kdr records (data/clean/kdr_records.csv)
#            since 2010 with one record each of 995F and 995S and a sample
#            size (allele counts consistent, as joint_observations() on branch
#            latent-kdr), against the empirical logit of its total kdr, log((y
#            + 0.5) / (2n - y + 0.5)) for y the 995F and 995S alleles of 2n:
#            correlation and R^2, for all samples and by species group
#   pairs    at the modelled cells, the correlation between u_s and u_f (and
#            v_s and u_f, and u_s and v_s, with the shear; and u_init with
#            u_s and u_f): per draw, as
#            posterior mean and 95% interval, and of the posterior means
#   spread   how much each smooth varies over the modelled cells: the sd,
#            range and central 95% of its posterior mean, and the mean over
#            draws of its sd over the cells; for u_s also as the multiplier
#            exp(u_s)
# Writes outputs/species_runs/smooth/<label>_validation.csv and
# figures/species_runs/smooth_validation_<label>.png. Plain R; about 3 GB
# and 1 minute.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2)
file <- arguments[1]
label <- arguments[2]
n_draws <- 500
records_since <- 2010

suppressMessages({
  library(greta)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})
source("R/dynamical_predictions.R")
source("R/species_fit_helpers.R")

output_dir <- "outputs/species_runs/smooth"
figure_dir <- "figures/species_runs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

fit <- load_fit(file)
if (!smooth_on(fit$options)) {
  report("%s has no latent smooths; nothing to check", label)
  quit(save = "no")
}
draws_used <- even_draws(fit, n_draws, label)
parameters <- fit_parameter_draws(fit, draws_used$index)
weights <- smooth_weight_terms(parameters$variables, parameters$options,
                               own = TRUE)
kinds <- intersect(c("selection", "selection_own", "floor", "init"),
                   names(weights))
shear <- if (isTRUE(fit$options$smooth$shear)) {
  b <- c(parameters$variables$smooth_shear)
  tibble(label = label, smooth = "shear b", n = length(b), mean = mean(b),
         q2.5 = quantile(b, 0.025), q97.5 = quantile(b, 0.975))
}
if (!is.null(shear)) {
  cat(sprintf("\n%s: shear loading b, u_s = v_s + b u_f: %.3f (95%%: %.3f to %.3f)\n",
              label, shear$mean, shear$q2.5, shear$q97.5))
}
crs <- fit$options$smooth$crs


# the marker layers at the modelled cells ---------------------------------------

cells <- fit$df$cell[match(seq_len(max(fit$df$cell_id)), fit$df$cell_id)]
at_cells <- smooth_posterior_at(parameters, smooth_cell_coords(cells, crs),
                                weights = weights)
standardise <- function(x) (x - mean(x)) / sd(x)
cell_data <- tibble(
  kdr = standardise(kdr_logit_at(cells, "complex")),
  arabiensis = standardise(qlogis(clamp_probability(
    arabiensis_fraction_at(cells)))))
for (kind in kinds) {
  cell_data[[kind]] <- at_cells[[kind]][, "mean"]
}

# R^2 of each layer alone and of both, and the partial correlations
partial_correlation <- function(y, x, z) {
  cor(resid(lm(y ~ z)), resid(lm(x ~ z)))
}
layers <- bind_rows(lapply(kinds, function(kind) {
  y <- cell_data[[kind]]
  tibble(
    label = label, smooth = kind, n = length(y),
    r2_kdr = summary(lm(y ~ cell_data$kdr))$r.squared,
    r2_arabiensis = summary(lm(y ~ cell_data$arabiensis))$r.squared,
    r2_both = summary(lm(y ~ cell_data$kdr + cell_data$arabiensis))$r.squared,
    partial_kdr = partial_correlation(y, cell_data$kdr,
                                      cell_data$arabiensis),
    partial_arabiensis = partial_correlation(y, cell_data$arabiensis,
                                             cell_data$kdr),
    slope_kdr = coef(lm(y ~ cell_data$kdr + cell_data$arabiensis))[[2]],
    slope_arabiensis = coef(lm(y ~ cell_data$kdr +
                                 cell_data$arabiensis))[[3]])
}))
options(width = 160)
cat(sprintf("\n%s: posterior mean smooths at %d modelled cells, on the standardised marker layers\n",
            label, length(cells)))
print(as.data.frame(layers), digits = 3)


# the smooths' correlations, and how much each varies, over the cells -----------

cell_basis <- smooth_basis_at(parameters$options$smooth,
                              smooth_cell_coords(cells, crs))
cell_draws <- lapply(weights[kinds], function(w) w %*% t(cell_basis))
for (kind in kinds) {
  # the draws give the posterior means of smooth_posterior_at()
  stopifnot(max(abs(colMeans(cell_draws[[kind]]) - cell_data[[kind]])) < 1e-8)
}
smooth_pairs <- Filter(function(pair) all(pair %in% kinds),
                       list(c("selection", "floor"),
                            c("selection_own", "floor"),
                            c("selection", "selection_own"),
                            c("selection", "init"),
                            c("floor", "init")))
pairs <- bind_rows(lapply(smooth_pairs, function(pair) {
  per_draw <- vapply(seq_len(nrow(cell_draws[[pair[1]]])), function(d) {
    cor(cell_draws[[pair[1]]][d, ], cell_draws[[pair[2]]][d, ])
  }, numeric(1))
  tibble(label = label, smooth = paste(pair, collapse = " ~ "),
         n = length(cells),
         correlation_of_means = cor(cell_data[[pair[1]]],
                                    cell_data[[pair[2]]]),
         mean = mean(per_draw), q2.5 = quantile(per_draw, 0.025),
         q97.5 = quantile(per_draw, 0.975))
}))
cat(sprintf("\n%s: correlations between the smooths over the %d modelled cells (per draw: mean, 95%%; and of the posterior means)\n",
            label, length(cells)))
print(as.data.frame(pairs), digits = 3)
spread <- bind_rows(lapply(kinds, function(kind) {
  u <- cell_data[[kind]]
  tibble(label = label, smooth = kind, n = length(u),
         sd_of_mean = sd(u), min_of_mean = min(u), max_of_mean = max(u),
         q2.5_of_mean = quantile(u, 0.025), q97.5_of_mean = quantile(u, 0.975),
         mean_sd_per_draw = mean(apply(cell_draws[[kind]], 1, sd)))
}))
cat(sprintf("\n%s: spread of each smooth over the %d modelled cells\n", label,
            length(cells)))
print(as.data.frame(spread), digits = 3)
for (kind in intersect(c("selection", "selection_own"), kinds)) {
  row <- spread[spread$smooth == kind, ]
  cat(sprintf(paste0("%s, %s: posterior mean sd %.2f over the cells; ",
                     "multiplier exp() %.2f to %.2f (central 95%% of cells ",
                     "%.2f to %.2f)\n"),
              label, kind, row$sd_of_mean, exp(row$min_of_mean),
              exp(row$max_of_mean), exp(row$q2.5_of_mean),
              exp(row$q97.5_of_mean)))
}


# the kdr records ----------------------------------------------------------------

alleles <- read.csv("data/clean/kdr_records.csv") %>%
  filter(year_start >= records_since, variant %in% c("995F", "995S"),
         !is.na(n_tested), n_tested > 0) %>%
  mutate(m = 2 * n_tested, y = round(frequency * m)) %>%
  group_by(sample_id) %>%
  filter(sum(variant == "995F") == 1, sum(variant == "995S") == 1) %>%
  summarise(longitude = first(longitude), latitude = first(latitude),
            species_group = first(species_group), m = first(m),
            y = sum(y), .groups = "drop") %>%
  # a count over by one is rounding
  mutate(y = if_else(y == m + 1, m, y)) %>%
  filter(y <= m) %>%
  mutate(empirical_logit = log((y + 0.5) / (m - y + 0.5)))
# the samples inside the mask
mask <- terra::rast("data/clean/raster_mask.tif")
sample_cells <- terra::cellFromXY(mask, cbind(alleles$longitude,
                                              alleles$latitude))
inside <- !is.na(sample_cells) &
  !is.na(terra::values(mask, mat = FALSE)[sample_cells])
report("%d samples since %d typed for both 995F and 995S, %d inside the mask",
       nrow(alleles), records_since, sum(inside))
alleles <- alleles[inside, ]
at_samples <- smooth_posterior_at(parameters,
                                  smooth_coords(alleles$longitude,
                                                alleles$latitude, crs),
                                  weights = weights)
for (kind in kinds) {
  alleles[[kind]] <- at_samples[[kind]][, "mean"]
}
groups <- c(list(all = alleles), split(alleles, alleles$species_group))
records <- bind_rows(lapply(names(groups), function(group) {
  data <- groups[[group]]
  bind_rows(lapply(kinds, function(kind) {
    tibble(label = label, smooth = kind, group = group, n = nrow(data),
           correlation = if (nrow(data) > 2) {
             cor(data$empirical_logit, data[[kind]])
           } else NA_real_,
           r2 = correlation ^ 2)
  }))
}))
cat(sprintf("\n%s: posterior mean smooths at the kdr samples since %d, against their empirical logit total kdr\n",
            label, records_since))
print(as.data.frame(records), digits = 3)

write.csv(bind_rows(layers %>% mutate(check = "layers", .before = 1),
                    records %>% mutate(check = "records", .before = 1),
                    pairs %>% mutate(check = "pairs", .before = 1),
                    spread %>% mutate(check = "spread", .before = 1),
                    if (!is.null(shear)) {
                      mutate(shear, check = "shear", .before = 1)
                    }),
          file.path(output_dir, sprintf("%s_validation.csv", label)),
          row.names = FALSE)


# the figure -----------------------------------------------------------------------

points <- bind_rows(
  cell_data %>%
    pivot_longer(all_of(kinds), names_to = "smooth", values_to = "u") %>%
    pivot_longer(c(kdr, arabiensis), names_to = "x_name", values_to = "x") %>%
    mutate(x_name = recode(x_name,
                           kdr = "logit total kdr, 2015 map (standardised)",
                           arabiensis = "logit arabiensis fraction (standardised)")),
  alleles %>%
    pivot_longer(all_of(kinds), names_to = "smooth", values_to = "u") %>%
    transmute(smooth, u, x = empirical_logit,
              x_name = sprintf("empirical logit total kdr, samples since %d",
                               records_since)))
p <- ggplot(points, aes(x, u)) +
  geom_point(size = 0.4, alpha = 0.3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) +
  facet_grid(smooth ~ x_name, scales = "free") +
  labs(x = NULL, y = "posterior mean of the smooth",
       title = sprintf("%s: the latent smooths against the kdr and species markers",
                       label),
       caption = paste("u_s (selection): log multiplier of the cumulative log",
                       "fitness; u_f (floor): shift of the logit floor;",
                       "v_s (selection_own): u_s less b u_f, with the shear;",
                       "u_init (init): the smooth of the logit relative",
                       "initial state, at sd 1")) +
  theme_bw(base_size = 9)
ggsave(file.path(figure_dir, sprintf("smooth_validation_%s.png", label)), p,
       width = 10, height = 2.4 + 2.6 * length(kinds), dpi = 150)
report("saved; peak memory %.1f GB", peak_memory_gb())
