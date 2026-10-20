# Checks of candidate smoothers for the kdr and species surfaces (#46): for
# each quantity (Vgsc-995F and 995S allele frequency, and the share of An.
# arabiensis among samples genotyped to species), quasi-binomial GAMs with a
# linear year trend and one of three spatial smooths:
#
#   tp80: thin-plate regression spline, 80 basis functions
#   tp30: thin-plate regression spline, 30 basis functions
#   gp:   Matern (kappa = 2.5) Gaussian process, range 10 degrees, 80 basis
#         functions, no linear trend, so far from data it reverts to the mean
#
# all with REML smoothness selection. Writes, per quantity, a figure of the
# prediction at 2015 over the whole mask (no masking, so the extrapolation is
# visible) and its standard error on the logit scale, with the record
# locations, and a table of convergence, edf, k-index and deviance explained.
#
#   Rscript R/kdr_surface_check.R
#
# Writes figures/kdr_surface_check_<quantity>.png and
# outputs/kdr/kdr_surface_check.csv.

source("R/packages.R")
suppressMessages(library(mgcv))

mask <- rast("data/clean/raster_mask.tif")
grid <- as.data.frame(terra::aggregate(mask, 4), xy = TRUE, na.rm = TRUE) |>
  transmute(longitude = x, latitude = y, year_c = 0)
records <- read.csv("data/clean/kdr_records.csv") |>
  mutate(year_c = year_start - 2015)

allele_data <- function(v) {
  records |>
    filter(variant == v, !is.na(n_tested), n_tested > 0) |>
    mutate(trials = 2 * n_tested, successes = round(frequency * trials))
}
quantities <- list(
  `995F` = allele_data("995F"),
  `995S` = allele_data("995S"),
  arabiensis = records |>
    filter(species_group %in% c("gambiae_ss", "coluzzii", "arabiensis")) |>
    distinct(sample_id, .keep_all = TRUE) |>
    mutate(trials = 1, successes = as.numeric(species_group == "arabiensis")))

smooths <- list(
  tp80 = "s(longitude, latitude, bs = 'tp', k = 80)",
  tp30 = "s(longitude, latitude, bs = 'tp', k = 30)",
  gp = "s(longitude, latitude, bs = 'gp', k = 80, m = c(-4, 10))")

dir.create("outputs/kdr", showWarnings = FALSE, recursive = TRUE)
checks <- list()
for (quantity in names(quantities)) {
  data <- quantities[[quantity]]
  maps <- list()
  for (smooth in names(smooths)) {
    formula <- as.formula(paste("cbind(successes, trials - successes) ~",
                                smooths[[smooth]], "+ year_c"))
    model <- gam(formula, family = stats::quasibinomial(), data = data,
                 method = "REML")
    k_check <- k.check(model)
    checks[[length(checks) + 1]] <- tibble(
      quantity, smooth, records = nrow(data),
      converged = model$converged,
      outer_converged = model$outer.info$conv,
      edf = round(sum(model$edf) - 2, 1),
      k = k_check[1, "k'"], k_index = round(k_check[1, "k-index"], 2),
      k_p = k_check[1, "p-value"],
      deviance_explained = round(summary(model)$dev.expl, 3),
      scale = round(model$scale, 1),
      year_slope = round(coef(model)[["year_c"]], 3))
    p <- predict(model, grid, se.fit = TRUE)
    maps[[smooth]] <- grid |>
      mutate(smooth = smooth, mean = plogis(as.numeric(p$fit)),
             se = as.numeric(p$se.fit))
  }
  maps <- bind_rows(maps) |>
    mutate(smooth = factor(smooth, names(smooths)))
  points <- distinct(data, longitude, latitude)
  base <- function(fill) {
    ggplot(maps, aes(longitude, latitude, fill = .data[[fill]])) +
      geom_raster() +
      geom_point(data = points, aes(longitude, latitude), inherit.aes = FALSE,
                 size = 0.15, colour = "grey50", alpha = 0.5) +
      facet_wrap(~ smooth, nrow = 1) +
      coord_equal() +
      labs(x = NULL, y = NULL) +
      theme_minimal() +
      theme(panel.grid = element_blank(), axis.text = element_blank())
  }
  top <- base("mean") +
    scale_fill_viridis_c(limits = c(0, 1)) +
    labs(fill = "2015", title = sprintf("%s: prediction at 2015", quantity))
  bottom <- base("se") +
    scale_fill_viridis_c(option = "magma", limits = c(0, 3),
                         oob = scales::squish) +
    labs(fill = "se (logit)", title = "standard error, logit scale")
  ggsave(sprintf("figures/kdr_surface_check_%s.png", quantity),
         patchwork::wrap_plots(top, bottom, ncol = 1), width = 12,
         height = 9, dpi = 110, bg = "white")
}
checks <- bind_rows(checks)
write.csv(checks, "outputs/kdr/kdr_surface_check.csv", row.names = FALSE)
print(as.data.frame(checks))
