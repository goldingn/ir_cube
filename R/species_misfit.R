# The level misfit of a saved full fit at the bioassays, and its pointwise log
# likelihood for PSIS-LOO (#37, #47).
#
#   Rscript R/species_misfit.R <fitted_model.RData> <label>
#
# From about n_draws (500) posterior draws, evenly spaced in each usable chain
# (all but those stuck), the predicted mortality at every modelled bioassay
# (dynamical_logit(), R/dynamical_predictions.R: with the species model, the
# mixture at the bioassay's arabiensis share; with the kdr covariate, at its
# cell's kdr), and
#   misfit      empirical logit, log((died + 0.5) / (survived + 0.5)), minus
#               the logit of the posterior mean predicted mortality: positive
#               where more died than predicted, i.e. the model predicts too
#               much resistance (red in the maps)
#   log lik     the beta-binomial log likelihood of each bioassay under each
#               draw (betabinomial_p_rho(), with the type's rho; p clamped to
#               [1e-12, 1 - 1e-12]), for PSIS-LOO (loo::loo(), with relative
#               efficiencies by chain). Bioassays are clustered (repeats at a
#               site and year, sites in a study), so leaving one out tests
#               interpolation, not prediction to new places or years. For a
#               fit with the weighted binomial likelihood (#47), its own log
#               likelihood, the binomial weighted by the design effect at
#               the fixed replicate rho (weighted_binomial_log_lik(),
#               R/weighted_binomial.R): not a normalised density of the
#               data, so its elpd_loo is NOT comparable with a beta-binomial
#               fit's, only with other weighted binomial fits' at the same rho
# Writes
#   outputs/species_runs/misfit/<label>_bioassays.csv  misfit per bioassay
#   outputs/species_runs/misfit/<label>_regions.csv    mean pyrethroid misfit
#                                                      by window and region
#   outputs/species_runs/misfit/<label>_loglik.rds     draws x bioassays
#   outputs/species_runs/misfit/<label>_loo.rds        the loo object
#   figures/species_runs/misfit_<label>.png            maps of the mean
#                                                      pyrethroid misfit per
#                                                      site, 2012-16 and
#                                                      2019-25, Africa and
#                                                      Ethiopia
# Plain R; FLOOR_PRIOR for fits saved without a floor prior (not used here,
# but load_fit() completes the options). About 3 GB and 2 minutes.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2)
file <- arguments[1]
label <- arguments[2]
n_draws <- 500

suppressMessages({
  library(greta)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(terra)
  library(tidyterra)
  library(sf)
  library(ggtext)
})
source("R/functions.R")
source("R/dynamical_predictions.R")
source("R/species_fit_helpers.R")

output_dir <- "outputs/species_runs/misfit"
figure_dir <- "figures/species_runs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

fit <- load_fit(file)
df <- fit$df
# about n_draws draws, the same number from each usable chain (even_draws())
draws_used <- even_draws(fit, n_draws, label)
draw_chain_id <- draws_used$chain
parameters <- fit_parameter_draws(fit, draws_used$index)
time <- system.time(
  logit <- dynamical_logit(parameters, df, df, fit$x_cell_years,
                           fit$cell_years_index)
)[["elapsed"]]
report("%s: predictions at %d bioassays x %d draws in %.0f s", label,
       ncol(logit), nrow(logit), time)


# misfit -------------------------------------------------------------------------------

p <- plogis(logit)
p_mean <- colMeans(p)
clamp <- function(x) pmin(pmax(x, 1e-12), 1 - 1e-12)

windows <- list(`2012-16` = 2012:2016, `2019-25` = 2019:2025)
bioassays <- df %>%
  transmute(longitude, latitude, cell, year_start, country_name,
            region = analysis_region(country_name, region),
            insecticide_class, insecticide_type, species, died,
            mosquito_number,
            empirical_logit = log((died + 0.5) /
                                    (mosquito_number - died + 0.5)),
            predicted = p_mean,
            predicted_logit = qlogis(clamp(p_mean)),
            misfit = empirical_logit - predicted_logit,
            window = case_when(year_start %in% windows[[1]] ~ names(windows)[1],
                               year_start %in% windows[[2]] ~ names(windows)[2]))
if (species_on(fit$options)) {
  bioassays$arabiensis_share <- arabiensis_share(df, fit$options)
}
write.csv(bioassays, file.path(output_dir, sprintf("%s_bioassays.csv", label)),
          row.names = FALSE)

pyrethroids <- bioassays %>%
  filter(insecticide_class == "Pyrethroids", !is.na(window))
summarise_misfit <- function(data) {
  summarise(data, bioassays = n(), sites = n_distinct(cell),
            mean_misfit = mean(misfit), median_misfit = median(misfit),
            .groups = "drop")
}
regions <- bind_rows(
  pyrethroids %>% group_by(window, area = region) %>% summarise_misfit(),
  pyrethroids %>% filter(country_name == "Ethiopia") %>%
    group_by(window) %>% summarise_misfit() %>% mutate(area = "Ethiopia"),
  pyrethroids %>% group_by(window) %>% summarise_misfit() %>%
    mutate(area = "all")) %>%
  mutate(label = label, .before = 1) %>%
  arrange(window, area)
write.csv(regions, file.path(output_dir, sprintf("%s_regions.csv", label)),
          row.names = FALSE)
options(width = 160)
cat("\nmean pyrethroid level misfit (empirical - predicted logit; positive:",
    "the model predicts too much resistance)\n")
print(as.data.frame(regions), digits = 3)


# pointwise log likelihood and PSIS-LOO ------------------------------------------------

# rho per type in each draw: the fit's, or with the weighted binomial, the
# fixed replicate rho (fixed_rho_types(), R/dynamical_model.R)
rho <- parameters$rho_types
weighted <- !rho_estimated(fit$options)
loglik <- matrix(NA_real_, nrow(p), ncol(p))
for (d in seq_len(nrow(p))) {
  rd <- rho[d, df$type_id]
  if (weighted) {
    # from the logit, which keeps log p and log(1 - p) precise near 0 and 1
    loglik[d, ] <- weighted_binomial_log_lik(
      df$died, df$mosquito_number,
      log_p = plogis(logit[d, ], log.p = TRUE),
      log_not_p = plogis(logit[d, ], lower.tail = FALSE, log.p = TRUE),
      weight = design_effect_weight(df$mosquito_number, rd))
    next
  }
  pd <- clamp(p[d, ])
  a <- pd * (1 / rd - 1)
  loglik[d, ] <- extraDistr::dbbinom(df$died, df$mosquito_number, alpha = a,
                                     beta = a * (1 - pd) / pd, log = TRUE)
}
stopifnot(all(is.finite(loglik)))
rm(logit, p)
invisible(gc())
saveRDS(loglik, file.path(output_dir, sprintf("%s_loglik.rds", label)))
# relative_eff() needs the chains numbered 1, 2, ..., whichever are used
r_eff <- loo::relative_eff(exp(loglik),
                           chain_id = match(draw_chain_id,
                                            unique(draw_chain_id)))
loo_fit <- loo::loo(loglik, r_eff = r_eff)
saveRDS(loo_fit, file.path(output_dir, sprintf("%s_loo.rds", label)))
cat(sprintf("\n%s: elpd_loo %.1f (se %.1f), p_loo %.1f, over %d bioassays, %d draws\n",
            label, loo_fit$estimates["elpd_loo", "Estimate"],
            loo_fit$estimates["elpd_loo", "SE"],
            loo_fit$estimates["p_loo", "Estimate"], ncol(loglik),
            nrow(loglik)))
print(loo::pareto_k_table(loo_fit))
cat("bioassays are clustered (repeats at a site and year, sites in a study):",
    "leave-one-out tests interpolation, and the elpd differences between fits",
    "are optimistic about their precision\n")
if (weighted) {
  cat(label, "has the weighted binomial likelihood: its elpd_loo is of the",
      "design-effect weighted binomial log likelihood, not comparable with",
      "the beta-binomial fits'\n")
}


# maps --------------------------------------------------------------------------------

borders <- readRDS("data/clean/country_borders.RDS")
# the limits of transmission, without water bodies, as the background
water_mask <- terra::aggregate(rast("data/clean/pfpr_water_mask.tif"), 3,
                               fun = "max", na.rm = TRUE)
sites <- pyrethroids %>%
  group_by(window, cell) %>%
  summarise(longitude = mean(longitude), latitude = mean(latitude),
            bioassays = n(), misfit = mean(misfit), .groups = "drop") %>%
  arrange(abs(misfit))
limit <- ceiling(quantile(abs(sites$misfit), 0.98))
# one diverging ramp: red, positive, where the model predicts too much
# resistance; blue too little
ramp <- scale_colour_gradientn(
  colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
  limits = c(-limit, limit), oob = scales::squish,
  name = "misfit\n(logit)")
map <- function(data, xlim = NULL, ylim = NULL, size = 2.2, title = NULL) {
  ggplot() +
    geom_sf(data = borders, fill = grey(0.9), colour = NA) +
    geom_spatraster(data = water_mask, show.legend = FALSE) +
    scale_fill_gradient(low = grey(0.78), high = grey(0.78),
                        na.value = "transparent", guide = "none") +
    geom_sf(data = borders, fill = NA, colour = grey(0.5), linewidth = 0.1) +
    geom_point(aes(longitude, latitude, colour = misfit), data = data,
               size = size, alpha = 0.9) +
    ramp +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_ir_maps()
}
panels <- list()
for (w in names(windows)) {
  data <- filter(sites, window == w)
  panels[[length(panels) + 1]] <- map(data, title = sprintf("%s, Africa", w))
  panels[[length(panels) + 1]] <- map(data, xlim = c(33, 48), ylim = c(3, 15),
                                      size = 3.5,
                                      title = sprintf("%s, Ethiopia", w))
}
p_maps <- wrap_plots(panels, ncol = 2, widths = c(1, 1)) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = sprintf("%s: pyrethroid level misfit at the bioassay sites", label),
    caption = paste("mean over a site's pyrethroid bioassays in the window of",
                    "the empirical logit minus the logit of the posterior mean",
                    "prediction; red: the model predicts too much resistance"))
ggsave(file.path(figure_dir, sprintf("misfit_%s.png", label)), p_maps,
       width = 10, height = 9, dpi = 150)
report("saved; peak memory %.1f GB", peak_memory_gb())
