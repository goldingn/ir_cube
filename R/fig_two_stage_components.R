# Supplementary figures for the two-stage model (#21): how the final model's
# components build a prediction at a pixel, for Deltamethrin.
#
#   figures/two_stage/supp_components_realisations.png (and .pdf)
#     joint random realisations of every component through time, stacked on a
#     shared time axis: net use, the dynamical model's logit m, the annual
#     anomalies eta, their accumulation xi, the static omega and p, the
#     combined logit m + omega + xi + p, and mortality with the data
#   figures/two_stage/supp_components_intervals.png (and .pdf)
#     posterior mean and 95% intervals of mortality: the dynamical model alone
#     against the two-stage final model, with the data
#
#   Rscript R/fig_two_stage_components.R
#
# Run with OpenBLAS (see R/two_stage_maps.R), e.g.
#   LD_PRELOAD=.../libopenblas.so.0 OPENBLAS_NUM_THREADS=4 nice -n 10 Rscript ...
#
# The final model (doc/two_stage_plan.md, "Final model"): on the logit scale
#   lambda(s, t) = m(s, t) + omega(s) + xi(s, t) + p(s) [+ u(s, t)],
# with m the dynamical model's prediction, xi(s, t) = psi xi(s, t - 1) +
# eta(s, t) (psi = 1 in the undamped model: xi is the sum of eta), eta AR(1)
# in time with Matern innovations, and p iid per pixel. Beyond the last data
# year T, eta is simulated forward. u, iid per pixel-year, is observation-level
# noise (with the assay noise), not part of the inferred process, so neither
# figure shows it: both target m + omega + xi + p.
#
# The fit saved by R/two_stage_maps.R (outputs/two_stage/maps/<type>/fit.rds)
# is a light copy without the Hessian factor, so this script refits the final
# model for the one type exactly as two_stage_maps.R does (same data, meshes,
# seed and code), checks the mode against the saved one, and caches the fit
# with the dynamical draws in outputs/two_stage/components/. Delete the cache
# to refit (about 20 minutes).
#
# `two_stage_code_dir` and `correction_template_path` (environment variables
# TWO_STAGE_CODE_DIR, TWO_STAGE_TEMPLATE) let the two-stage code and template
# be taken from a pinned copy rather than R/ and tmb/.

type <- "Deltamethrin"

# the model of R/two_stage_maps.R (model_config there): "omega_xi_u_p_pql", or
# with damped xi "omega_xi_u_p_psi_pql"
model_config <- "omega_xi_u_p_pql"
damped_xi <- grepl("_psi", model_config)

code_dir <- Sys.getenv("TWO_STAGE_CODE_DIR", "R")
template_path <- Sys.getenv("TWO_STAGE_TEMPLATE",
                            "tmb/two_stage_correction.cpp")

cache_dir <- "outputs/two_stage/components"
model_tag <- if (damped_xi) "_psi" else ""
cache_file <- file.path(cache_dir, sprintf("%s%s_fit_cache.rds", type,
                                           model_tag))
draws_file <- file.path(cache_dir, sprintf("%s%s_pixel_draws.rds", type,
                                           model_tag))
figure_dir <- "figures/two_stage"
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

report <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", sprintf(...), "\n")
  flush(stdout())
}

suppressMessages({
  sink("/dev/null")
  source("R/validation_folds.R")
  source("R/validation_covariates.R")
  sink()
})
source(file.path(code_dir, "dynamical_predictions.R"))
source(file.path(code_dir, "two_stage_correction.R"))
source(file.path(code_dir, "two_stage_pql.R"))
source(file.path(code_dir, "two_stage_map_functions.R"))
correction_template <- template_path

# settings of R/two_stage_maps.R
mesh_config <- "omega5000_xi2500"
t0 <- baseline_year
end_year <- 2030
years_all <- baseline_year:end_year
clamp <- 1e-12
safe_logit <- function(p) qlogis(pmin(pmax(p, clamp), 1 - clamp))
logit_max <- qlogis(1 - clamp)

k <- match(type, types)
rows_k <- which(df$type_id == k)


# 1. the fit and the dynamical draws (cached) -----------------------------------

if (!file.exists(cache_file)) {

  fit_env <- new.env()
  load("temporary/fitted_model.RData", envir = fit_env)
  stopifnot(isTRUE(all.equal(fit_env$df, df)),
            identical(fit_env$types, types),
            identical(fit_env$unique_cells, unique_cells))
  fold <- list(draws = fit_env$draws, options = fit_env$model_options,
                x_cells_init = fit_env$x_cells_init)
  rm(fit_env)
  invisible(gc())

  draw_index <- paired_draw_index(fold)
  draws_matrix <- as.matrix(fold$draws)[draw_index, , drop = FALSE]
  logit_init_mean <- logit_init_mean_draws(fold, draw_index)
  parameters <- dynamical_parameter_draws(fold, df = df,
                                          classes_index = classes_index,
                                          types = types,
                                          draw_index = draw_index,
                                          logit_init_mean = logit_init_mean)

  # m_ref over all assays, then this type's columns (as two_stage_maps.R)
  p_train <- dynamical_predictions(fold, select(df, -country_id), df,
                                   x_cell_years, cell_years_index,
                                   classes_index, types,
                                   draw_index = draw_index)
  logit_train_k <- safe_logit(p_train[, rows_k, drop = FALSE])
  rm(p_train)
  m_ref <- colMeans(logit_train_k)

  lookup <- country_region_lookup()
  logit_init_all <- map_logit_init(draws_matrix, logit_init_mean, types,
                                   classes_index, countries, regions, lookup)
  logit_init_k <- logit_init_all[, , k]
  effect_k <- parameters$effect_type[, , k]
  rm(fold, draws_matrix, parameters, logit_init_all)
  invisible(gc())

  rho_table <- read.csv("outputs/bioassay_rho_hierarchical.csv")
  rho <- rho_table$rho[rho_table$insecticide_type == type]
  train_k <- tibble(
    lon = df$longitude[rows_k],
    lat = df$latitude[rows_k],
    year = df$year_start[rows_k],
    cell = df$cell[rows_k],
    died = df$died[rows_k],
    mosquito_number = df$mosquito_number[rows_k],
    m = m_ref,
    rho = rho
  )
  stage_a <- empirical_logit(train_k$died, train_k$mosquito_number, rho)
  train_k$z <- stage_a$z
  train_k$v <- stage_a$v
  T_k <- max(train_k$year)

  coords <- coords_km(train_k)
  meshes <- suppressMessages(build_correction_meshes(coords, mesh_config))

  set.seed(2026 + k)
  time_fit <- system.time({
    fit_a <- fit_correction(train_k, variant = "omega_xi_u", t0 = t0,
                            T = T_k, mesh = meshes$omega,
                            mesh_xi = meshes$xi, pixel_effect = TRUE,
                            damped_xi = damped_xi)
    stopifnot(fit_a$opt$convergence == 0)
    fit <- fit_correction_pql(fit_a, train_k)
  })
  rm(fit_a)
  report("%s refitted in %.0f s", type, time_fit[["elapsed"]])

  # the refit must be the fit behind the maps
  saved <- readRDS(file.path("outputs/two_stage/maps", type, "fit.rds"))
  stopifnot(identical(read.csv(file.path("outputs/two_stage/maps", type,
                                         "hyperparameters.csv"))$model_config,
                      model_config))
  mode_diff <- max(abs(fit$mode - saved$mode))
  report("max |mode - saved mode| = %.2e; phi %.4f vs %.4f", mode_diff,
         fit$hyper$phi, saved$hyper$phi)
  stopifnot(mode_diff < 1e-3)

  # drop the TMB object (external pointers) before caching
  fit$obj <- NULL
  saveRDS(list(fit = fit, train = train_k, logit_train = logit_train_k,
               effect = effect_k, logit_init = logit_init_k, rho = rho),
          cache_file)
  rm(fit, train_k, logit_train_k, effect_k, logit_init_k)
  invisible(gc())
}

cache <- readRDS(cache_file)
fit <- cache$fit
train_k <- cache$train
hyper <- fit$hyper
T_k <- fit$T


# 2. example pixels ---------------------------------------------------------------

# The data per pixel, and the fitted smooth correction at the mode
# (omega + xi at T) at each sampled pixel
coords_train <- coords_km(train_k)
pixel_summary <- train_k %>%
  mutate(x_km = coords_train[, 1], y_km = coords_train[, 2]) %>%
  group_by(cell) %>%
  summarise(lon = mean(lon), lat = mean(lat),
            x_km = mean(x_km), y_km = mean(y_km),
            n_assays = n(), n_years = n_distinct(year),
            first = min(year), last = max(year),
            n_before = sum(year <= 2008), n_after = sum(year >= 2014),
            .groups = "drop")
nodes_T <- correction_node_fields(fit, fit$mode, T_k)
pixel_coords <- as.matrix(pixel_summary[, c("x_km", "y_km")])
pixel_summary$smooth_T <- as.vector(
  mesh_basis(fit$mesh, pixel_coords) %*% nodes_T$omega +
    mesh_basis(fit$mesh_xi, pixel_coords) %*% nodes_T$xi[[1]])

# (i) well sampled: data both before the net scale-up (to 2008) and after it
# (from 2014), and the most distinct years of data among such pixels
well <- pixel_summary %>%
  filter(n_before > 0, n_after > 0) %>%
  arrange(desc(n_years), desc(n_assays)) %>%
  slice(1)

# (iii) large correction: among pixels with at least three years of data, the
# largest |omega + xi| at T at the mode, i.e. where the data pull the
# prediction furthest from the dynamical model
large <- pixel_summary %>%
  filter(n_years >= 3, cell != well$cell) %>%
  arrange(desc(abs(smooth_T))) %>%
  slice(1)

# (ii) no data: a cell of the prediction mask inside the limits of Pf
# transmission, 400-600 km from the nearest Deltamethrin assay, in a country
# with Deltamethrin data: far beyond the range of omega (~43 km), so omega and
# p revert to their priors and the prediction to the dynamical model, but
# within the reach of the large-scale xi (range ~1000 km). The candidate with
# the highest net use in 2020 is taken, so that the dynamical model has a
# trajectory to follow
pf_water_mask <- rast("data/clean/pfpr_water_mask.tif")
country_raster <- rast("data/clean/country_raster.tif")
mask_cells <- terra::cells(terra::mask(mask, pf_water_mask))
set.seed(1)
candidate_cells <- sort(sample(mask_cells, min(length(mask_cells), 40000)))
xy_candidates <- terra::xyFromCell(mask, candidate_cells)
coords_candidates <- project_km(xy_candidates[, 1], xy_candidates[, 2])
nearest_km <- apply(coords_candidates, 1, function(xy) {
  sqrt(min((pixel_coords[, 1] - xy[1]) ^ 2 + (pixel_coords[, 2] - xy[2]) ^ 2))
})
candidate_country <- as.character(
  terra::extract(country_raster, candidate_cells)$country_name)
data_countries <- unique(df$country_name[rows_k])
nets_2020 <- terra::extract(rast("data/clean/net_use_cube.tif")[["nets_2020"]],
                            candidate_cells)[, 1]
ok <- nearest_km > 400 & nearest_km < 600 &
  candidate_country %in% data_countries & !is.na(nets_2020)
none_cell <- candidate_cells[ok][which.max(nets_2020[ok])]
none <- tibble(cell = none_cell,
               lon = xy_candidates[match(none_cell, candidate_cells), 1],
               lat = xy_candidates[match(none_cell, candidate_cells), 2],
               n_assays = 0L, n_years = 0L,
               nearest_km = nearest_km[match(none_cell, candidate_cells)])

pixels <- bind_rows(
  mutate(well, role = "well"),
  mutate(none, role = "none"),
  mutate(large, role = "large")
) %>%
  mutate(country = as.character(
    terra::extract(country_raster, cell)$country_name))
report("pixels: %s", paste(sprintf("%s cell %i (%s, %.2f, %.2f; %i assays, %i years)",
                                  pixels$role, pixels$cell, pixels$country,
                                  pixels$lon, pixels$lat, pixels$n_assays,
                                  pixels$n_years), collapse = "; "))


# 3. draws at the pixels ------------------------------------------------------------

# Every component at each pixel and year 1995-2030, for n_draws joint draws:
# draw d pairs dynamical draw d with a latent draw from N(mode, H^-1) shifted by
# the cut-posterior formula for that dynamical draw (as predict_correction()
# does), the AR(1) forecast of eta beyond T with fresh Matern innovations, and
# p from the latent draw where the pixel has data and fresh from its prior
# elsewhere. (u, observation noise, is not drawn)
if (!file.exists(draws_file) ||
    !identical(readRDS(draws_file)$cells, pixels$cell)) {

  n_draws <- nrow(cache$logit_train)
  n_pix <- nrow(pixels)
  n_years_all <- length(years_all)

  # dynamical logit at the pixels, all years (years x pixels x draws)
  covariates <- map_covariates(pixels$cell, baseline_year, end_year)
  pixel_country <- match(pixels$country, dimnames(cache$logit_init)[[2]])
  stopifnot(!anyNA(pixel_country))
  dyn <- dynamical_logit_chunk(
    effect = cache$effect,
    logit_init = t(cache$logit_init[, pixel_country, drop = FALSE]),
    time_varying = covariates$time_varying,
    flat = covariates$flat,
    years = years_all, years_keep = years_all)
  m <- aperm(simplify2array(dyn), c(3, 1, 2))
  m <- pmin(pmax(m, -logit_max), logit_max)

  # check: at the sampled pixels' data years, the same draws as at the assays
  for (i in which(pixels$n_assays > 0)) {
    rows_i <- which(train_k$cell == pixels$cell[i])
    j <- match(train_k$year[rows_i], years_all)
    m_i <- t(matrix(m[j, i, ], length(j), n_draws))
    difference <- max(abs(m_i - cache$logit_train[, rows_i]))
    report("pixel %i: max |m - m at the assays| = %.1e", pixels$cell[i],
           difference)
    stopifnot(difference < 1e-6)
  }

  pixel_xy <- terra::xyFromCell(mask, pixels$cell)
  coords_pix <- project_km(pixel_xy[, 1], pixel_xy[, 2])
  A_omega <- mesh_basis(fit$mesh, coords_pix)
  A_xi <- mesh_basis(fit$mesh_xi, coords_pix)
  n_nodes_xi <- fit$mesh_xi$n
  n_latent <- length(fit$mode)

  Q_eta <- matern_precision_r(fit$fem_xi, hyper$kappa_eta, hyper$sigma_eta)
  Q_eta_chol <- Matrix::Cholesky(Matrix::forceSymmetric(Q_eta),
                                 perm = TRUE, LDL = FALSE, super = TRUE)

  psi <- correction_psi(fit)
  p_index <- fit$pixels$p_index[match(pixels$cell, fit$pixels$cell)]

  empty <- function() array(NA_real_, c(n_years_all, n_pix, n_draws))
  xi <- empty()
  eta <- empty()
  omega <- matrix(NA_real_, n_pix, n_draws)
  p <- matrix(NA_real_, n_pix, n_draws)

  set.seed(4026 + k)
  batches <- split(seq_len(n_draws), ceiling(seq_len(n_draws) / 250))
  for (batch in batches) {
    nb <- length(batch)
    theta <- sample_latent_deviation(fit$H_chol, n_latent, nb) + fit$mode +
      correction_mode_shift(fit, t(cache$logit_train[batch, , drop = FALSE]))
    omega[, batch] <- as.matrix(A_omega %*% theta[fit$blocks$w_omega, ])
    x <- theta[fit$blocks$x, , drop = FALSE]
    xi_nodes <- matrix(0, n_nodes_xi, nb)
    eta_nodes <- matrix(0, n_nodes_xi, nb)
    for (j in seq_along(years_all)) {
      y <- years_all[j]
      if (y > t0 && y <= T_k) {
        xi_new <- x[(y - t0 - 1) * n_nodes_xi + seq_len(n_nodes_xi), ,
                    drop = FALSE]
        eta_nodes <- xi_new - psi * xi_nodes
        xi_nodes <- xi_new
      } else if (y > T_k) {
        eta_nodes <- hyper$phi * eta_nodes + sqrt(1 - hyper$phi ^ 2) *
          sample_latent_deviation(Q_eta_chol, n_nodes_xi, nb)
        xi_nodes <- psi * xi_nodes + eta_nodes
      }
      xi[j, , batch] <- as.matrix(A_xi %*% xi_nodes)
      eta[j, , batch] <- as.matrix(A_xi %*% eta_nodes)
    }
    p[, batch] <- matrix(rnorm(n_pix * nb, 0, hyper$sigma_p), n_pix, nb)
    seen_p <- !is.na(p_index)
    if (any(seen_p)) {
      p[seen_p, batch] <- theta[fit$blocks$p[p_index[seen_p]], ,
                                drop = FALSE]
    }
  }

  nets <- covariates$time_varying[, , "nets"]
  saveRDS(list(cells = pixels$cell, pixels = pixels, years = years_all,
               m = m, omega = omega, xi = xi, eta = eta, p = p, psi = psi,
               nets = nets, hyper = hyper, T = T_k),
          draws_file)
  report("draws at %i pixels x %i years x %i draws saved", n_pix,
         n_years_all, n_draws)
}


# 4. figures ------------------------------------------------------------------------

draws <- readRDS(draws_file)
pixels <- draws$pixels
years <- draws$years
n_draws <- dim(draws$m)[3]
T_k <- draws$T

# panel titles, in the order of `pixels`
role_title <- c(well = "A) Well sampled",
                none = "B) No data",
                large = "C) Large correction")
pixels <- pixels %>%
  mutate(title = sprintf("%s: %s\n%s; %s",
                         role_title[role], country,
                         sprintf("%.1f°%s, %.1f°%s",
                                 abs(lat), ifelse(lat >= 0, "N", "S"),
                                 abs(lon), ifelse(lon >= 0, "E", "W")),
                         ifelse(n_assays > 0,
                                sprintf("%i assays, %i years", n_assays,
                                        n_years),
                                sprintf("nearest assay %.0f km",
                                        nearest_km))),
         title = factor(title, levels = title))

# the style of the dynamical-model time-series figures (R/summarise_model_fit.R,
# R/fig_temporal_preds_net_use.R): theme_minimal, the pyrethroid blue, net use
# as a thick grey line, percentages, no x label
pyrethroid_blue <- "#56B1F7"
net_grey <- grey(0.5)
realisation_cols <- c("#E69F00", "#009E73", "#CC79A7")
year_breaks <- seq(2000, 2030, by = 10)

projection_shading <- list(
  annotate("rect", xmin = T_k + 0.5, xmax = 2030.5, ymin = -Inf, ymax = Inf,
           fill = grey(0.93)),
  geom_vline(xintercept = T_k + 0.5, colour = grey(0.5), linewidth = 0.3,
             linetype = "dashed")
)
base_theme <- theme_minimal(base_size = 9) +
  theme(strip.text.x = element_text(hjust = 0, size = 8.5),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 8.5),
        legend.position = "none",
        plot.margin = margin(2, 4, 2, 4))
x_scale <- scale_x_continuous(breaks = year_breaks, limits = c(1994.5, 2030.5),
                              expand = c(0, 0))

# long data frames: one row per pixel x year (x draw)
pixel_year <- function(values, name) {
  tibble(title = rep(pixels$title, each = length(years)),
         year = rep(years, nrow(pixels)),
         !!name := as.vector(values))
}
# realisation draws: years x pixels x draws array (or pixels x draws matrix
# for the static terms) at the chosen draws
realisations <- function(values, which) {
  if (length(dim(values)) == 2) {
    values <- array(rep(values, each = length(years)),
                    c(length(years), dim(values)))
  }
  bind_rows(lapply(seq_along(which), function(r) {
    pixel_year(values[, , which[r]], "value") %>%
      mutate(realisation = factor(r))
  }))
}

bioassays <- train_k %>%
  filter(cell %in% pixels$cell) %>%
  mutate(title = pixels$title[match(cell, pixels$cell)],
         mortality = died / mosquito_number)

# 4a. realisations ------------------------------------------------------------------

# three joint draws, spread over the posterior of the dynamical model: draw d
# pairs the dynamical draw d with the correction's latent draw d
set.seed(5)
which_draws <- sort(sample(n_draws, 3))

# the inferred process m + omega + xi + p; u is observation noise and left out
lambda <- draws$m + draws$xi +
  array(rep(draws$omega + draws$p, each = length(years)), dim(draws$m))

row_plot <- function(data, ylab, geoms, y_scale = NULL, strip = FALSE,
                     bottom = FALSE) {
  plot <- ggplot(data, aes(x = year)) +
    projection_shading +
    geoms +
    facet_wrap(~title, nrow = 1) +
    x_scale +
    ylab(ylab) +
    xlab(NULL) +
    base_theme
  if (!is.null(y_scale)) plot <- plot + y_scale
  if (!strip) plot <- plot + theme(strip.text.x = element_blank())
  if (!bottom) plot <- plot + theme(axis.text.x = element_blank())
  plot
}
realisation_colour <- scale_colour_manual(values = realisation_cols)
zero_line <- geom_hline(yintercept = 0, colour = grey(0.6), linewidth = 0.3)

m_mean <- pixel_year(apply(draws$m, 1:2, mean), "value")

p_nets <- row_plot(
  pixel_year(t(draws$nets), "value"), "LLIN use",
  geom_line(aes(y = value), colour = net_grey, linewidth = 1),
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     breaks = c(0, 0.5, 1)),
  strip = TRUE)
p_m <- row_plot(
  realisations(draws$m, which_draws), "dynamical\nlogit m",
  list(geom_line(aes(y = value), data = m_mean, colour = "black",
                 linewidth = 0.6),
       geom_line(aes(y = value, colour = realisation), linewidth = 0.4),
       realisation_colour))
p_eta <- row_plot(
  realisations(draws$eta, which_draws), "annual\nanomaly η",
  list(zero_line,
       geom_line(aes(y = value, colour = realisation), linewidth = 0.3,
                 alpha = 0.6),
       geom_point(aes(y = value, colour = realisation), size = 0.6),
       realisation_colour))
p_xi <- row_plot(
  realisations(draws$xi, which_draws),
  if (draws$psi < 1) "damped sum\nξ = ψξ + η" else "cumulative\nξ = Ση",
  list(zero_line,
       geom_line(aes(y = value, colour = realisation), linewidth = 0.5),
       realisation_colour))
p_omega <- row_plot(
  realisations(draws$omega, which_draws), "static\nfield ω",
  list(zero_line,
       geom_line(aes(y = value, colour = realisation), linewidth = 0.5),
       realisation_colour))
p_p <- row_plot(
  realisations(draws$p, which_draws), "pixel\neffect p",
  list(zero_line,
       geom_line(aes(y = value, colour = realisation), linewidth = 0.5),
       realisation_colour))
p_lambda <- row_plot(
  realisations(lambda, which_draws), "logit\nm+ω+ξ+p",
  list(geom_line(aes(y = value, colour = realisation), linewidth = 0.4),
       realisation_colour))
p_mort <- row_plot(
  realisations(plogis(lambda), which_draws), "mortality",
  list(geom_point(aes(y = mortality, size = mosquito_number),
                  data = bioassays, shape = 21, fill = grey(0.8),
                  colour = grey(0.35), stroke = 0.2, alpha = 0.7),
       geom_line(aes(y = value, colour = realisation), linewidth = 0.4),
       realisation_colour,
       scale_size_area(max_size = 2.5)),
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     breaks = c(0, 0.5, 1)),
  bottom = TRUE)

fig_realisations <- patchwork::wrap_plots(p_nets, p_m, p_eta, p_xi, p_omega, p_p, p_lambda,
                      p_mort, ncol = 1,
                      heights = c(1, 1.1, 1, 1.1, 0.8, 0.8, 1.1, 1.2)) &
  theme(plot.margin = margin(1, 4, 1, 4))
for (ext in c("png", "pdf")) {
  ggsave(file.path(figure_dir, sprintf("supp_components_realisations.%s", ext)),
         plot = fig_realisations,
         bg = "white", width = 7.5, height = 9.2, dpi = 300,
         device = if (ext == "pdf") cairo_pdf else NULL)
}


# 4b. intervals ---------------------------------------------------------------------

# mortality intervals for the population fraction at each pixel-year: the
# dynamical model alone (plogis(m)) and the two-stage model
# (plogis(m + omega + xi + p), with the cut-posterior shift). u is observation
# noise, like the assay-level (beta-binomial) noise, so neither is included
summarise_draws <- function(values, model) {
  pixel_year(apply(values, 1:2, mean), "mean") %>%
    mutate(lower = as.vector(apply(values, 1:2, quantile, 0.025)),
           upper = as.vector(apply(values, 1:2, quantile, 0.975)),
           model = model)
}
intervals <- bind_rows(
  summarise_draws(plogis(draws$m), "dynamical model"),
  summarise_draws(plogis(lambda), "two-stage model")
) %>%
  mutate(model = factor(model, c("dynamical model", "two-stage model")))
model_cols <- c(`dynamical model` = grey(0.45),
                `two-stage model` = pyrethroid_blue)

p_int <- ggplot(intervals, aes(x = year)) +
  projection_shading +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = 0.35) +
  geom_line(aes(y = mean, colour = model), linewidth = 0.6) +
  geom_point(aes(y = mortality, size = mosquito_number), data = bioassays,
             shape = 21, fill = "white", colour = "black", stroke = 0.3,
             alpha = 0.8) +
  facet_wrap(~title, nrow = 1) +
  scale_fill_manual(values = model_cols, name = NULL) +
  scale_colour_manual(values = model_cols, name = NULL) +
  scale_size_area(max_size = 2.5, guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  x_scale +
  xlab(NULL) +
  ylab("Susceptibility to deltamethrin") +
  base_theme +
  theme(legend.position = "top",
        legend.justification = "left",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        axis.text.x = element_blank())
p_int_nets <- row_plot(
  pixel_year(t(draws$nets), "value"), "LLIN use",
  geom_line(aes(y = value), colour = net_grey, linewidth = 1),
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     breaks = c(0, 0.5, 1)),
  bottom = TRUE)

fig_intervals <- patchwork::wrap_plots(p_int, p_int_nets, ncol = 1, heights = c(3, 1))
for (ext in c("png", "pdf")) {
  ggsave(file.path(figure_dir, sprintf("supp_components_intervals.%s", ext)),
         plot = fig_intervals,
         bg = "white", width = 8, height = 4.5, dpi = 300,
         device = if (ext == "pdf") cairo_pdf else NULL)
}
report("figures written")
