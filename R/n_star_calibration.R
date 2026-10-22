# EXPERIMENT STUB (archived, not part of the pipeline): calibrating the maps'
# n* (#40, #59) against cross-validation, and a CV-adjusted n* map.
#
# n*_local, mapped from the posterior of ilogit(m + omega + xi + u + p), is the
# model's own estimate of the map's worth in bioassays. n*_CV (R/validation_level.R)
# is the out-of-sample one. At the held-out bioassays, with F_i the bioassay
# noise (the rho-hat floor), s_i^2 the fold's posterior variance and e_i^2 the
# squared error of its posterior mean, calibration means E[e^2] = s^2 + F, so
#
#   n*_pred = mean F / mean s^2        (the map's n*, at the held-out sizes)
#   n*_CV   = mean F / (MSE - mean F)
#   kappa   = n*_CV / n*_pred = mean s^2 / (MSE - mean F)
#
# kappa is 1 / the geostatistical mean squared deviation ratio (MSDR), after
# removing the bioassay noise. N cancels, so kappa x n* is a CV-adjusted map at
# any N (n_star_cv in the figures step of R/two_stage_maps.R).
#
# Findings on the V5pg folds (22-23 Oct 2026):
#   - kappa_local 0.49 [0.46, 0.53] pooled; interpolation 0.56 [0.40, 0.93],
#     extrapolation 0.48 [0.43, 0.54], forecasts 0.50 [0.46, 0.54]; LLIN
#     pyrethroids 0.50 [0.46, 0.55]. The maps' n* is ~2x too high out of sample
#   - flat-ish in distance to the nearest training assay (pooled 0.56 within
#     25 km, 0.41-0.46 beyond 50 km) and in forecast horizon
#   - binned by mapped n*, kappa falls with n*, but the high-n* bins are
#     near-100%-mortality types (noisy excess) and a few gross pyrethroid
#     misses: not robust. Variance forms b s^2 + a [m(1 - m)]^2 do not transfer
#     between experiments
#   - per type it varies (Alpha-cypermethrin 0.18, Bendiocarb 0.71)
#   - coverage of the predictive intervals is near nominal (95% -> 90-93%):
#     insensitive, as the bioassay noise dominates the interval
#   - kappa x n* rather than inflating s^2 by 1 / kappa: inflation also shrinks
#     E[p(1 - p)] (effective ratio 0.32-0.45) and is undefined where
#     s^2 / kappa > m(1 - m)
# To recompute on the #58 folds (sy, CCD hyperparameters) before any use;
# with sy, its variance belongs in s^2 here (it is in the held-out assays) but
# not in the maps.
#
# Reads outputs/cv_scores.csv, outputs/cv_draws/two_stage__*.rds,
# outputs/cv_geometry.csv and outputs/bioassay_rho_type_draws.rds. Prints the
# tables; writes nothing.

suppressMessages({
  library(dplyr)
  library(tidyr)
})

set.seed(2026 - 10 - 22)
n_boot <- 2000
options(width = 200)

col_var <- function(x) {
  n <- nrow(x)
  mu <- colMeans(x)
  (colSums(x ^ 2) - n * mu ^ 2) / (n - 1)
}

# noise_floor_record() (R/validation_functions.R), for a matrix of rho
floor_at <- function(died, n, rho) {
  yhat <- died / n
  inflation <- (1 + (n - 1) * rho) / n
  pmax(yhat * (1 - yhat) / (1 - inflation), 0) * inflation
}


# per held-out assay ----------------------------------------------------------------

scores <- read.csv("outputs/cv_scores.csv", colClasses = c(fold = "character"))
scores <- scores[scores$model == "two_stage", ]
folds <- tibble(
  exp = c("B", "C", "C", "F", "F"),
  experiment = c("spatial_interpolation", "spatial_blocks", "spatial_blocks",
                 "temporal_forecasting_2014", "temporal_forecasting_2018"),
  fold = c("all", "1", "2", "2014", "2018"),
  file = c("spatial_interpolation__all", "spatial_blocks__1",
           "spatial_blocks__2", "temporal_forecasting__2014",
           "temporal_forecasting__2018"))

d <- bind_rows(lapply(seq_len(nrow(folds)), function(i) {
  f <- folds[i, ]
  x <- readRDS(file.path("outputs/cv_draws",
                         paste0("two_stage__", f$file, ".rds")))
  s <- scores[scores$experiment == f$experiment & scores$fold == f$fold, ]
  stopifnot(nrow(s) == ncol(x$map_draws),
            max(abs(s$map - colMeans(x$map_draws))) < 1e-9)
  out <- tibble(exp = f$exp, experiment = f$experiment, fold = f$fold,
                cell = s$cell, year = s$year_start, type = s$insecticide_type,
                died = s$died, n = s$mosquito_number, y = s$observed,
                map = s$map, local = colMeans(x$p_draws),
                s2_pop = col_var(x$map_draws), s2_local = col_var(x$p_draws),
                floor = s$floor, rho = s$rho_external, pit = s$pit,
                pit_pop = s$cdf_below_map + runif(nrow(s)) * s$pmf_at_map)
  rm(x)
  invisible(gc())
  out
})) %>%
  filter(!is.na(floor)) %>%
  mutate(e2_pop = (y - map) ^ 2, e2_local = (y - local) ^ 2,
         # n* at N = 100, as n_star() maps it
         n100_local = (local * (1 - local) - s2_local) * (1 + 99 * rho) /
           (100 * s2_local))

geometry <- read.csv("outputs/cv_geometry.csv",
                     colClasses = c(fold = "character"))
d <- d %>%
  left_join(geometry %>% select(experiment, fold, cell, year_start,
                                insecticide_type, distance_same),
            by = c("experiment", "fold", "cell", "year" = "year_start",
                   "type" = "insecticide_type"),
            relationship = "many-to-one")


# kappa with a pixel bootstrap ------------------------------------------------------

# Means over assays from per-pixel sums; a replicate resamples pixels as
# multinomial counts and takes a posterior draw of rho by type for the floors
rho_all <- readRDS("outputs/bioassay_rho_type_draws.rds")
rho_draws <- rho_all[sample(nrow(rho_all), n_boot), ]

kappa <- function(g, reading = "local") {
  if (n_distinct(g$cell) < 5) return("-")
  s2 <- g[[paste0("s2_", reading)]]
  e2 <- g[[paste0("e2_", reading)]]
  pixel <- factor(paste(g$exp, g$cell))
  s2_pixel <- rowsum(s2, pixel)[, 1]
  e2_pixel <- rowsum(e2, pixel)[, 1]
  floor_pixel <- rowsum(floor_at(g$died, g$n,
                                 t(rho_draws[, g$type, drop = FALSE])), pixel)
  W <- rmultinom(n_boot, nlevels(pixel), rep(1, nlevels(pixel)))
  r <- colSums(W * s2_pixel) /
    (colSums(W * e2_pixel) - colSums(W * floor_pixel))
  sprintf("%.2f [%.2f, %.2f]", sum(s2) / (sum(e2) - sum(g$floor)),
          quantile(r, 0.025), quantile(r, 0.975))
}

by_group <- function(data, ..., reading = "local") {
  data %>%
    group_by(...) %>%
    group_modify(~ tibble(assays = nrow(.x), kappa = kappa(.x, reading)),
                 .keep = TRUE) %>%
    ungroup()
}

cat("\n== kappa by experiment (B interpolation, C extrapolation, F forecasts) ==\n")
print(as.data.frame(
  by_group(d, exp) %>% rename(kappa_local = kappa) %>%
    left_join(by_group(d, exp, reading = "pop") %>%
                select(exp, kappa_pop = kappa), by = "exp")), row.names = FALSE)

cat("\n== kappa_local pooled, LLIN pyrethroids, and per type ==\n")
pyrethroids <- c("Alpha-cypermethrin", "Deltamethrin", "Permethrin")
print(as.data.frame(bind_rows(
  tibble(group = "all", by_group(d, exp) %>% pivot_wider(
    names_from = exp, values_from = c(assays, kappa))),
  tibble(group = "LLIN pyrethroids", by_group(filter(d, type %in% pyrethroids),
                                              exp) %>%
           pivot_wider(names_from = exp, values_from = c(assays, kappa))),
  by_group(d, type, exp) %>% rename(group = type) %>%
    pivot_wider(names_from = exp, values_from = c(assays, kappa))) %>%
  left_join(bind_rows(tibble(group = "all", pooled = kappa(d)),
                      tibble(group = "LLIN pyrethroids",
                             pooled = kappa(filter(d, type %in% pyrethroids))),
                      by_group(d, type) %>% transmute(group = type,
                                                      pooled = kappa)),
            by = "group")), row.names = FALSE)

cat("\n== kappa_local by distance to the nearest training assay (km) ==\n")
print(as.data.frame(
  d %>% mutate(km = cut(distance_same, c(0, 25, 50, 100, Inf),
                        right = FALSE)) %>%
    by_group(km, exp) %>%
    mutate(kappa = sprintf("%s (%i)", kappa, assays)) %>%
    select(-assays) %>%
    pivot_wider(names_from = exp, values_from = kappa)), row.names = FALSE)

cat("\n== kappa_local by mapped n*_local (N = 100), B and C ==\n")
print(as.data.frame(
  d %>% filter(exp %in% c("B", "C")) %>%
    mutate(bin = cut(n100_local, c(0, 0.5, 1, 2, 5, Inf), right = FALSE)) %>%
    by_group(exp, bin)), row.names = FALSE)


# smooth variance forms -------------------------------------------------------------

# E[e^2 - F] = b s^2 + a [m(1 - m)]^2 by least squares (the data floor, as
# n*_CV), fitted on C and checked by binned kappa on every experiment
d <- d %>% mutate(excess = e2_local - floor, s2 = s2_local,
                  v = (local * (1 - local)) ^ 2)
forms <- list(b = "s2", a = "v", ab = c("s2", "v"))
for (form in names(forms)) {
  cf <- coef(lm(reformulate(forms[[form]], "excess", intercept = FALSE),
                data = filter(d, exp == "C")))
  adjusted <- Reduce(`+`, Map(function(term, k) k * d[[term]], names(cf), cf))
  cat(sprintf("\n== form %s fitted on C: %s; kappa raw > adjusted ==\n", form,
              paste(names(cf), signif(cf, 3), sep = " = ", collapse = ", ")))
  print(as.data.frame(
    d %>% mutate(adjusted = adjusted,
                 bin = cut(n100_local, c(0, 0.5, 1, 2, 5, Inf),
                           right = FALSE)) %>%
      group_by(exp, bin) %>%
      summarise(v = sprintf("%.2f > %.2f (%i)", mean(s2) / mean(excess),
                            mean(adjusted) / mean(excess), n()),
                .groups = "drop") %>%
      pivot_wider(names_from = exp, values_from = v)), row.names = FALSE)
}


# coverage --------------------------------------------------------------------------

cat("\n== coverage of central predictive intervals (randomised PIT), % ==\n")
inside <- function(pit, level) 100 * mean(abs(pit - 0.5) < level / 2)
print(as.data.frame(d %>% group_by(exp) %>% summarise(
  assays = n(),
  local_50 = inside(pit, 0.5), local_80 = inside(pit, 0.8),
  local_95 = inside(pit, 0.95),
  pop_50 = inside(pit_pop, 0.5), pop_80 = inside(pit_pop, 0.8),
  pop_95 = inside(pit_pop, 0.95)) %>%
  mutate(across(-c(exp, assays), ~ round(.x, 1)))), row.names = FALSE)
