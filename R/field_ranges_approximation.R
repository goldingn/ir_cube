# Are the latent smooths' fitted ranges (V5, #47; R/latent_smooth.R) where
# their Hilbert-space approximation is accurate? The posterior of each
# smooth's range and sd, by chain, against the approximation error of the
# fit's own basis at ranges of 200-4,000 km: as R/check_latent_smooth.R, the
# covariance the centred basis implies at the modelled bioassay cells against
# the exactly centred kernel, over all pairs of cells, as the maximum and root
# mean square absolute error relative to sd^2 (which does not depend on sd),
# and the largest error of a variance (the diagonal). Also the share of the
# kernel's variance at frequencies beyond the basis's highest, omega_max,
# which it cannot represent: for the squared exponential with lengthscale ell =
# range / 2 in two dimensions, exp(-ell^2 omega_max^2 / 2); and the
# half-wavelength of omega_max, pi / omega_max, the shortest feature the basis
# can draw.
#
#   Rscript R/field_ranges_approximation.R <file>=<label> ...
#
# per fit, its usable chains (USE_CHAINS, e.g. "V5=1,2"). Writes
#   outputs/species_runs/grids/field_ranges_approximation.csv  the posterior
#     summaries by fit, chain and smooth, with the error at the median and
#     the 5% quantile of the range
#   outputs/species_runs/grids/field_approximation_error.csv   the error curve
#   figures/species_runs/grids/field_ranges_vs_approximation.png
# Plain R; about 2 GB and 5 minutes.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) >= 1, all(grepl("=", arguments)))
files <- setNames(sub("=.*$", "", arguments), sub("^.*=", "", arguments))
error_ranges <- c(seq(0.2, 1.5, by = 0.05), seq(1.6, 4, by = 0.1))
design_range <- 1

suppressMessages({
  library(greta)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
source("R/dynamical_predictions.R")
source("R/species_fit_helpers.R")

output_dir <- "outputs/species_runs/grids"
figure_dir <- "figures/species_runs/grids"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# the hyperparameter draws of each fit's usable chains, and its basis and
# modelled cells
hyper <- list()
bases <- list()
for (label in names(files)) {
  fit <- load_fit(files[[label]])
  stopifnot(smooth_on(fit$options))
  chains <- usable_chains(fit$draws, label)
  for (chain in chains) {
    m <- as.matrix(fit$draws[[chain]])
    for (kind in smooth_kinds(fit$options)) {
      names <- smooth_variable_names(kind)
      hyper[[length(hyper) + 1]] <- tibble(
        fit = label, chain = chain, smooth = kind,
        sd = m[, names[["sd"]]],
        range_km = 1000 / m[, names[["inv_range"]]])
    }
  }
  n_cells <- max(fit$df$cell_id)
  bases[[label]] <- list(
    smooth = fit$options$smooth,
    cells = fit$df$cell[match(seq_len(n_cells), fit$df$cell_id)])
  rm(fit)
  gc()
}
hyper <- bind_rows(hyper)
# every fit must have the same basis and cells, for one error curve
for (label in names(bases)[-1]) {
  stopifnot(identical(bases[[label]]$cells, bases[[1]]$cells),
            isTRUE(all.equal(bases[[label]]$smooth[c("origin", "half_width",
                                                     "indices", "kernel")],
                             bases[[1]]$smooth[c("origin", "half_width",
                                                 "indices", "kernel")])))
}
smooth <- bases[[1]]$smooth
cells <- bases[[1]]$cells
omega_max <- min(pi * smooth$m / (2 * smooth$half_width))
half_wavelength_km <- 1000 * pi / omega_max
rates <- smooth_prior_rates(smooth)
report("basis: %s kernel, m = (%d, %d), %d functions, half-widths (%.0f, %.0f) km; omega_max %.3f per 1,000 km, half-wavelength pi / omega_max = %.0f km",
       smooth$kernel, smooth$m[1], smooth$m[2], nrow(smooth$indices),
       1000 * smooth$half_width[1], 1000 * smooth$half_width[2], omega_max,
       half_wavelength_km)


# the approximation error against range ----------------------------------------

coords <- smooth_cell_coords(cells, smooth$crs)
distance <- as.matrix(dist(coords))
basis <- smooth_basis_at(smooth, coords)
omega <- hsgp_frequencies(smooth$indices, smooth$half_width)
kernel <- function(d, rho) {
  ell <- rho / 2
  switch(smooth$kernel,
         matern52 = {
           r <- sqrt(5) * d / ell
           (1 + r + r ^ 2 / 3) * exp(-r)
         },
         se = exp(-d ^ 2 / (2 * ell ^ 2)))
}
centre_both <- function(k) {
  k <- sweep(k, 1, rowMeans(k))
  sweep(k, 2, colMeans(k))
}
time <- system.time({
  error_curve <- bind_rows(lapply(error_ranges, function(rho) {
    s <- smooth_sqrt_spectral(omega, 1, 1 / rho, smooth$kernel)
    approximate <- tcrossprod(sweep(basis, 2, c(s), FUN = "*"))
    difference <- approximate - centre_both(kernel(distance, rho))
    tibble(range_km = 1000 * rho,
           max_error = max(abs(difference)),
           rms_error = sqrt(mean(difference ^ 2)),
           max_variance_error = max(abs(diag(difference))),
           truncated_variance = if (smooth$kernel == "se") {
             exp(-(rho / 2) ^ 2 * omega_max ^ 2 / 2)
           } else NA_real_)
  }))
})[["elapsed"]]
report("error curve at %d ranges in %.0f s", length(error_ranges), time)
write.csv(error_curve, file.path(output_dir, "field_approximation_error.csv"),
          row.names = FALSE)
options(width = 160)
print(as.data.frame(error_curve %>%
                      filter(range_km %in% c(200, 300, 400, 500, 550, 600,
                                             650, 700, 800, 900, 1000, 1500,
                                             2000, 4000))),
      digits = 3, row.names = FALSE)

# the error at range r (km), interpolated on log range
error_at <- function(r, column) {
  approx(log(error_curve$range_km), error_curve[[column]], log(r),
         rule = 2)$y
}


# the posteriors against the prior and the error -------------------------------

summarise_hyper <- function(d) {
  d %>%
    summarise(
      n = n(),
      range_median = median(range_km),
      range_q05 = quantile(range_km, 0.05),
      range_q95 = quantile(range_km, 0.95),
      share_below_1000 = mean(range_km < 1000),
      share_below_1500 = mean(range_km < 1500),
      prior_cdf_at_median = exp(-rates$range / (range_median / 1000)),
      max_error_at_median = error_at(range_median, "max_error"),
      rms_error_at_median = error_at(range_median, "rms_error"),
      max_error_at_q05 = error_at(range_q05, "max_error"),
      rms_error_at_q05 = error_at(range_q05, "rms_error"),
      truncated_at_median = error_at(range_median, "truncated_variance"),
      sd_median = median(sd),
      sd_q05 = quantile(sd, 0.05),
      sd_q95 = quantile(sd, 0.95),
      prior_tail_at_sd_median = exp(-rates$sd * sd_median),
      .groups = "drop")
}
posterior <- bind_rows(
  hyper %>% group_by(fit, smooth) %>% summarise_hyper() %>%
    mutate(chain = "all"),
  hyper %>% group_by(fit, chain, smooth) %>% summarise_hyper() %>%
    mutate(chain = as.character(chain))) %>%
  relocate(chain, .after = fit)
write.csv(posterior, file.path(output_dir, "field_ranges_approximation.csv"),
          row.names = FALSE)
cat("\nposterior range (km) and sd of each smooth, the error (share of sd^2) at the range's median and 5% quantile, and the priors' probabilities:\n")
posterior %>%
  mutate(across(c(range_median, range_q05, range_q95), round),
         across(c(share_below_1000, share_below_1500, max_error_at_median,
                  rms_error_at_median, max_error_at_q05, rms_error_at_q05,
                  truncated_at_median, sd_median, sd_q05, sd_q95), ~ round(.x, 3)),
         across(c(prior_cdf_at_median, prior_tail_at_sd_median),
                ~ signif(.x, 2))) %>%
  as.data.frame() %>%
  print(row.names = FALSE)


# the figure --------------------------------------------------------------------

hyper <- hyper %>%
  mutate(series = paste0(fit, " ch", chain),
         smooth = factor(smooth, levels = c("selection", "floor"),
                         labels = c("u_s (selection)", "u_f (floor)")))
series_levels <- sort(unique(hyper$series))
series_colours <- setNames(
  c("#2a78d6", "#7fb2ef", "#eb6834", "#1baf7a", "#4a3aa7", "#e87ba4",
    "#eda100", "#008300")[seq_along(series_levels)],
  series_levels)
range_limits <- c(200, 4000)
range_breaks <- c(200, 300, 500, 700, 1000, 1500, 2000, 3000, 4000)
markers <- tibble(
  range_km = c(half_wavelength_km, 1000 * design_range,
               1000 * smooth$range_prior[1]),
  what = c(sprintf("half-wavelength\npi / omega_max\n(%.0f km)",
                   half_wavelength_km),
           "design limit\n(1,000 km)",
           "prior 5% quantile\n(1,500 km)"))
marker_lines <- function() {
  Map(function(x, type) {
    geom_vline(xintercept = x, linetype = type, colour = grey(0.35),
               linewidth = 0.4)
  }, markers$range_km, c("dotted", "dashed", "longdash"))
}
theme_ranges <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0))

# the largest error over pairs is that of a variance at every range here, so
# the variances' are not drawn apart
stopifnot(isTRUE(all.equal(error_curve$max_error,
                           error_curve$max_variance_error)))
error_long <- error_curve %>%
  pivot_longer(c(max_error, rms_error, truncated_variance),
               names_to = "measure", values_to = "error") %>%
  mutate(measure = factor(measure,
                          levels = c("max_error", "rms_error",
                                     "truncated_variance"),
                          labels = c("max over pairs (a variance)",
                                     "rms over pairs",
                                     "kernel variance beyond omega_max")))
error_panel <- ggplot(error_long, aes(x = range_km, y = 100 * error,
                                      linetype = measure)) +
  marker_lines() +
  geom_text(data = markers, aes(x = range_km, y = Inf, label = what),
            inherit.aes = FALSE, vjust = 1.1, hjust = -0.05, size = 2.6,
            lineheight = 0.9, colour = grey(0.3)) +
  geom_line(linewidth = 0.7) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted"),
                        name = NULL) +
  scale_x_log10(limits = range_limits, breaks = range_breaks) +
  theme_ranges +
  theme(legend.position = c(0.82, 0.55),
        legend.background = element_rect(fill = "white", colour = NA)) +
  labs(x = NULL, y = "Error (% of sd^2)",
       title = "Approximation error of the basis at the modelled cells, against range")

range_panel <- ggplot(hyper, aes(x = range_km, colour = series)) +
  marker_lines() +
  geom_density(linewidth = 0.6, adjust = 1.2) +
  facet_grid(smooth ~ ., scales = "free_y") +
  scale_colour_manual(values = series_colours, name = NULL) +
  scale_x_log10(limits = range_limits, breaks = range_breaks) +
  theme_ranges +
  theme(legend.position = "top") +
  labs(x = "Range (km, log scale)", y = "Posterior density\n(per log10 km)",
       title = "Posterior ranges",
       subtitle = sprintf(
         "Prior: P(range < 1,500 km) = 0.05; P(range < %.0f km) = %.1g",
         600, exp(-rates$range / 0.6)))

sd_panel <- ggplot(hyper, aes(x = sd, colour = series)) +
  geom_vline(xintercept = smooth$sd_prior[1], linetype = "longdash",
             colour = grey(0.35), linewidth = 0.4) +
  geom_density(linewidth = 0.6, adjust = 1.2) +
  facet_grid(smooth ~ ., scales = "free_y") +
  scale_colour_manual(values = series_colours, name = NULL) +
  scale_x_continuous(limits = c(0, NA)) +
  theme_ranges +
  theme(legend.position = "none") +
  labs(x = "Marginal sd", y = "Posterior density",
       title = "Posterior sds",
       subtitle = sprintf(
         "Prior: P(sd > 1) = 0.05 (dashed line); P(sd > 3) = %.1g",
         exp(-rates$sd * 3)))

ranges_plot <- (error_panel / range_panel / sd_panel) +
  plot_layout(heights = c(1.1, 1.6, 1.4))
ggsave(file.path(figure_dir, "field_ranges_vs_approximation.png"),
       ranges_plot, bg = "white", width = 9, height = 11)
report("wrote %s", file.path(figure_dir, "field_ranges_vs_approximation.png"))
