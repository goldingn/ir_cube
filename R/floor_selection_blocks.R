# A pseudo-likelihood screen for V5 (#47): do the selection multipliers and
# mortality floors that pyrethroid bioassays need locally lie on one axis (one
# latent spatial smooth could drive both) or vary independently (two)?
#
#   Rscript R/floor_selection_blocks.R [block size, degrees: 5] [fit file]
#
# The base is the reference fit ref_f0 (no floor, kdr or species), every
# parameter held at its posterior mean over the usable chains. At every
# pyrethroid bioassay it gives, in the outer form (outer_mortality(),
# R/dynamical_model.R), the logit initial state l0, the cumulative log fitness
# C_t and the reversion t kappa, checked against the fit's own predictions.
# The bioassays are grouped in square spatial blocks (5 degrees; 7.5 as a
# check), keeping blocks with at least min_bioassays pyrethroid bioassays over
# at least min_years years. In each block b, (s_b, u_b = logit f_b) maximise
#   sum_i log BetaBinomial(died_i | n_i, p_i, rho_type)
#     + log N(s_b; 0, 1) + log N(u_b; -1.83, 1.39^2)
#   p_i = f_b + (1 - f_b) plogis(l0_i - exp(s_b) C_i - t_i kappa_i),
# a selection multiplier exp(s_b) and a floor f_b; the priors keep the
# estimates finite, the floor's that of the floors of the species runs (the
# logit moments of Beta(1, 4)). The inverse Hessian of the negative log
# posterior gives each block's covariance of (s_b, u_b).
#
# Then:
#   correlation  across blocks of (s_b, u_b), corrected for estimation noise
#                by a bivariate random-effects meta-analysis: y_b ~ N(mu,
#                Sigma + S_b), S_b the block's covariance; Sigma, the
#                between-block covariance, by maximum likelihood, its
#                correlation with a Wald interval on atanh() and a bootstrap
#                over blocks
#   one axis     the summed log posterior with (s_b, u_b) free against the
#                constraint that every block lies on one line in (s, u)
#                space, any direction, fitted jointly (one latent value per
#                block): likelihood difference and AIC (2B against B + 2
#                parameters)
#   regions      each block's region (analysis_region() of most of its
#                bioassays; Ethiopia apart), and the mean of each region's
#                estimates
# Writes outputs/species_runs/blocks/blocks_<size>.csv, summary_<size>.csv,
# and figures/species_runs/floor_selection_blocks_<size>.png (the scatter with
# 50% and 95% ellipses) and floor_selection_blocks_map_<size>.png (maps of
# s_b and f_b; red where selection is weaker or the floor higher, the
# direction of a missing plateau). Plain R; about 2 GB.
#
# Caveats: the base fit is held fixed (its initial states, covariate effects
# and reversion absorb what they can before the blocks are fitted); the blocks
# are arbitrary; and in a short or flat series a slower decline and a higher
# floor are hard to tell apart, which the Hessian covariance carries.

arguments <- commandArgs(trailingOnly = TRUE)
block_size <- if (length(arguments) >= 1) as.numeric(arguments[1]) else 5
fit_file <- if (length(arguments) >= 2) arguments[2] else
  paste0("../ir_cube_netscreen/outputs/pod_jobs/dh270_lin_f0_full/",
         "temporary/fitted_model.RData")
min_bioassays <- 40
min_years <- 8
# the floor's prior: the logit moments of Beta(1, 4) (logit_beta_moments())
floor_prior <- list(mean = -1.8333, sd = 1.3888)
selection_prior_sd <- 1

suppressMessages({
  library(greta)
  library(dplyr)
  library(tidyr)
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

output_dir <- "outputs/species_runs/blocks"
figure_dir <- "figures/species_runs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
size_label <- format(block_size)

stopifnot(all.equal(logit_beta_moments(c(1, 4))$mean, floor_prior$mean,
                    tolerance = 1e-4),
          all.equal(logit_beta_moments(c(1, 4))$sd, floor_prior$sd,
                    tolerance = 1e-4))


# the base fit, at its posterior mean -----------------------------------------------

fit <- load_fit(fit_file)
stopifnot(!species_on(fit$options), !kdr_on(fit$options),
          isFALSE(fit$options$mortality_floor))
usable <- usable_chains(fit$draws, "base")
means <- colMeans(as.matrix(fit$draws[usable]))
parameters <- dynamical_parameter_draws(
  list(draws = coda::mcmc.list(coda::mcmc(matrix(means, 1,
                                                 dimnames = list(NULL,
                                                                 names(means))))),
       options = fit$options, x_cells_init = fit$x_cells_init),
  fit$classes_index, fit$types, draw_index = 1, options = fit$options)
df <- fit$df
rho <- c(parameters$rho_types)
kappa <- c(parameters$kappa_type)

pyrethroids <- df %>%
  filter(insecticide_class == "Pyrethroids") %>%
  mutate(row = row_number())

# l0, C and t kappa at each bioassay, as dynamical_trajectories() forms them
cell_country <- dynamical_lookups(df)$cell_country_lookup
n_times <- max(fit$cell_years_index$year_id)
x_row <- matrix(NA_integer_, max(fit$cell_years_index$cell_id), n_times)
x_row[cbind(fit$cell_years_index$cell_id, fit$cell_years_index$year_id)] <-
  seq_len(nrow(fit$cell_years_index))
pyrethroids$l0 <- NA_real_
pyrethroids$C <- NA_real_
for (k in sort(unique(pyrethroids$type_id))) {
  rows_k <- which(pyrethroids$type_id == k)
  cells <- sort(unique(pyrethroids$cell_id[rows_k]))
  x <- array(fit$x_cell_years[as.vector(x_row[cells, ]), , drop = FALSE],
             c(length(cells), n_times, ncol(fit$x_cell_years)))
  l0 <- cell_logit_init(parameters, k,
                        matrix(parameters$logit_init_relative[
                          , cell_country[cells], k], 1),
                        parameters$x_cells_init[cells, , drop = FALSE])
  effect <- matrix(parameters$effect_type[, , k], 1)
  cumulative <- matrix(NA_real_, length(cells), n_times)
  running <- 0
  for (t in seq_len(n_times)) {
    running <- running + log1p(effect %*% t(matrix(x[, t, ],
                                                   nrow = length(cells))))
    cumulative[, t] <- running
  }
  at <- match(pyrethroids$cell_id[rows_k], cells)
  pyrethroids$l0[rows_k] <- l0[1, at]
  pyrethroids$C[rows_k] <- cumulative[cbind(at, pyrethroids$year_id[rows_k])]
}
pyrethroids <- pyrethroids %>%
  mutate(reversion = year_id * kappa[type_id],
         rho = rho[type_id],
         logit_base = l0 - C - reversion)
# the outer form reproduces the fit's own predictions at these parameters
check <- c(dynamical_logit(parameters, pyrethroids, df, fit$x_cell_years,
                           fit$cell_years_index))
difference <- max(abs(pmin(pmax(check, -30), 30) -
                        pmin(pmax(pyrethroids$logit_base, -30), 30)))
report("outer-form quantities at %d pyrethroid bioassays; against the fit's predictions, max |logit diff| %.2g",
       nrow(pyrethroids), difference)
stopifnot(difference < 1e-8)


# blocks ------------------------------------------------------------------------------

pyrethroids <- pyrethroids %>%
  mutate(block_x = floor(longitude / block_size),
         block_y = floor(latitude / block_size),
         block = paste(block_x, block_y),
         area = analysis_region(country_name, region),
         area = ifelse(country_name == "Ethiopia", "Ethiopia", area))
block_table <- pyrethroids %>%
  group_by(block, block_x, block_y) %>%
  summarise(bioassays = n(),
            years = max(year_start) - min(year_start) + 1,
            first_year = min(year_start), last_year = max(year_start),
            cells = n_distinct(cell),
            area = names(which.max(table(area))),
            .groups = "drop")
kept <- block_table %>%
  filter(bioassays >= min_bioassays, years >= min_years)
report("%s-degree blocks: %d with pyrethroid bioassays, %d kept (>= %d bioassays over >= %d years), holding %d of %d bioassays",
       size_label, nrow(block_table), nrow(kept), min_bioassays, min_years,
       sum(kept$bioassays), nrow(pyrethroids))


# the per-block fits -------------------------------------------------------------------

clamp <- function(p) pmin(pmax(p, 1e-12), 1 - 1e-12)

# the beta-binomial log likelihood of a block's bioassays at (s, u)
block_log_lik <- function(data, s, u) {
  f <- plogis(u)
  p <- clamp(f + (1 - f) * plogis(data$l0 - exp(s) * data$C -
                                    data$reversion))
  a <- p * (1 / data$rho - 1)
  sum(extraDistr::dbbinom(data$died, data$mosquito_number, alpha = a,
                          beta = (1 - p) * (1 / data$rho - 1), log = TRUE))
}
log_prior <- function(s, u) {
  dnorm(s, 0, selection_prior_sd, log = TRUE) +
    dnorm(u, floor_prior$mean, floor_prior$sd, log = TRUE)
}
block_log_post <- function(data, theta) {
  block_log_lik(data, theta[1], theta[2]) + log_prior(theta[1], theta[2])
}

# the maximum from a grid of starts, and the covariance from the Hessian
fit_block <- function(data) {
  starts <- expand.grid(s = c(-1.5, 0, 1), u = c(-5, -1.8, 0))
  runs <- lapply(seq_len(nrow(starts)), function(i) {
    optim(unlist(starts[i, ]), function(theta) -block_log_post(data, theta),
          method = "BFGS", control = list(reltol = 1e-12, maxit = 500))
  })
  best <- runs[[which.min(vapply(runs, `[[`, numeric(1), "value"))]]
  hessian <- optimHess(best$par,
                       function(theta) -block_log_post(data, theta))
  covariance <- solve(hessian)
  tibble(s = best$par[1], u = best$par[2],
         se_s = sqrt(covariance[1, 1]), se_u = sqrt(covariance[2, 2]),
         cov_su = covariance[1, 2],
         log_post = -best$value,
         log_lik = block_log_lik(data, best$par[1], best$par[2]),
         log_lik_base = block_log_lik(data, -Inf, -Inf),
         converged = best$convergence == 0)
}

by_block <- split(pyrethroids, pyrethroids$block)
time <- system.time(
  estimates <- bind_rows(lapply(kept$block, function(b) {
    fit_block(by_block[[b]]) %>% mutate(block = b, .before = 1)
  }))
)[["elapsed"]]
blocks <- left_join(kept, estimates, by = "block") %>%
  mutate(f = plogis(u), multiplier = exp(s),
         correlation_within = cov_su / (se_s * se_u))
report("fitted %d blocks in %.0f s; %d converged", nrow(blocks), time,
       sum(blocks$converged))
write.csv(blocks, file.path(output_dir, sprintf("blocks_%s.csv", size_label)),
          row.names = FALSE)


# correlation across blocks, corrected for estimation noise ------------------------------

# y_b ~ N(mu, Sigma + S_b): Sigma from log sds and atanh of the correlation
meta_log_lik <- function(par, y, S) {
  mu <- par[1:2]
  sds <- exp(par[3:4])
  r <- tanh(par[5])
  Sigma <- diag(sds) %*% matrix(c(1, r, r, 1), 2) %*% diag(sds)
  total <- 0
  for (b in seq_len(nrow(y))) {
    R <- tryCatch(chol(Sigma + S[[b]]), error = function(e) NULL)
    if (is.null(R)) return(-1e10)
    z <- backsolve(R, y[b, ] - mu, transpose = TRUE)
    total <- total - sum(log(diag(R))) - 0.5 * sum(z ^ 2) - log(2 * pi)
  }
  total
}
fit_meta <- function(y, S) {
  start <- c(colMeans(y), log(apply(y, 2, sd)), 0)
  optim(start, function(par) -meta_log_lik(par, y, S), method = "BFGS",
        control = list(maxit = 1000, reltol = 1e-12))
}
y <- cbind(blocks$s, blocks$u)
S <- lapply(seq_len(nrow(blocks)), function(b) {
  matrix(c(blocks$se_s[b] ^ 2, blocks$cov_su[b], blocks$cov_su[b],
           blocks$se_u[b] ^ 2), 2)
})
meta <- fit_meta(y, S)
meta_hessian <- optimHess(meta$par, function(par) -meta_log_lik(par, y, S))
z_se <- sqrt(solve(meta_hessian)[5, 5])
set.seed(47)
bootstrap <- replicate(1000, {
  i <- sample(nrow(y), replace = TRUE)
  tanh(fit_meta(y[i, , drop = FALSE], S[i])$par[5])
})
raw_correlation <- cor(blocks$s, blocks$u)
median_se_s <- median(blocks$se_s)
median_se_u <- median(blocks$se_u)
correlation <- tibble(
  block_size = block_size,
  blocks = nrow(blocks),
  raw_correlation = raw_correlation,
  between_correlation = tanh(meta$par[5]),
  wald_lower = tanh(meta$par[5] - 1.96 * z_se),
  wald_upper = tanh(meta$par[5] + 1.96 * z_se),
  bootstrap_lower = quantile(bootstrap, 0.025),
  bootstrap_upper = quantile(bootstrap, 0.975),
  between_sd_s = exp(meta$par[3]),
  between_sd_u = exp(meta$par[4]),
  mean_s = meta$par[1], mean_u = meta$par[2],
  median_within_sd_s = median_se_s,
  median_within_sd_u = median_se_u)


# one axis against two -----------------------------------------------------------------

# every block on the line {d n + z t}, t = (cos theta, sin theta), n = (-sin
# theta, cos theta): for given (theta, d), each block's z by a 1-D maximum
profile_line <- function(theta, d) {
  direction <- c(cos(theta), sin(theta))
  normal <- c(-sin(theta), cos(theta))
  vapply(kept$block, function(b) {
    data <- by_block[[b]]
    -optimize(function(z) {
      point <- d * normal + z * direction
      -block_log_post(data, point)
    }, c(-12, 12), tol = 1e-8)$objective
  }, numeric(1))
}
# start from the principal axis of the free estimates
axis <- prcomp(y)
theta_start <- atan2(axis$rotation[2, 1], axis$rotation[1, 1])
d_start <- sum(c(-sin(theta_start), cos(theta_start)) * axis$center)
line <- optim(c(theta_start, d_start),
              function(par) -sum(profile_line(par[1], par[2])),
              method = "Nelder-Mead",
              control = list(reltol = 1e-10, maxit = 500))
line_log_post <- -line$value
free_log_post <- sum(blocks$log_post)
n_blocks <- nrow(blocks)
one_axis <- tibble(
  block_size = block_size,
  blocks = n_blocks,
  log_post_free = free_log_post,
  log_post_one_axis = line_log_post,
  difference = free_log_post - line_log_post,
  parameters_free = 2 * n_blocks,
  parameters_one_axis = n_blocks + 2,
  aic_free = -2 * free_log_post + 2 * 2 * n_blocks,
  aic_one_axis = -2 * line_log_post + 2 * (n_blocks + 2),
  delta_aic_one_axis_minus_free = aic_one_axis - aic_free,
  lrt_p = pchisq(2 * (free_log_post - line_log_post), n_blocks - 2,
                 lower.tail = FALSE),
  line_theta_degrees = line$par[1] * 180 / pi,
  line_slope_du_ds = tan(line$par[1]))

# regions ----------------------------------------------------------------------------------

areas <- c("West", "Central", "East", "Southern", "Horn", "Ethiopia")
regions <- blocks %>%
  group_by(area) %>%
  summarise(blocks = n(), bioassays = sum(bioassays),
            mean_s = mean(s), se_mean_s = sd(s) / sqrt(n()),
            mean_u = mean(u), se_mean_u = sd(u) / sqrt(n()),
            mean_multiplier = mean(exp(s)), mean_floor = mean(plogis(u)),
            .groups = "drop") %>%
  arrange(match(area, areas))

summary <- list(correlation = correlation, one_axis = one_axis,
                regions = regions,
                counts = tibble(block_size = block_size,
                                blocks_with_data = nrow(block_table),
                                blocks_kept = nrow(blocks),
                                bioassays_kept = sum(blocks$bioassays),
                                bioassays = nrow(pyrethroids)))
saveRDS(summary, file.path(output_dir, sprintf("summary_%s.rds", size_label)))
write.csv(bind_cols(summary$counts, select(correlation, -block_size, -blocks),
                    select(one_axis, -block_size, -blocks)),
          file.path(output_dir, sprintf("summary_%s.csv", size_label)),
          row.names = FALSE)
options(width = 160)
print(as.data.frame(summary$counts))
print(as.data.frame(correlation), digits = 3)
print(as.data.frame(one_axis), digits = 6)
print(as.data.frame(regions), digits = 3)


# figures ----------------------------------------------------------------------------------

area_colours <- c(West = "#0072B2", Central = "#009E73", East = "#E69F00",
                  Southern = "#CC79A7", Horn = "#D55E00", Ethiopia = "#000000",
                  North = grey(0.5))

# the ellipse of a bivariate normal at coverage `level`
ellipse <- function(centre, covariance, level, n = 60) {
  angle <- seq(0, 2 * pi, length.out = n)
  radius <- sqrt(qchisq(level, 2))
  e <- eigen(covariance, symmetric = TRUE)
  points <- t(e$vectors %*% diag(sqrt(pmax(e$values, 0))) %*%
                rbind(cos(angle), sin(angle))) * radius
  tibble(s = centre[1] + points[, 1], u = centre[2] + points[, 2])
}
ellipses <- bind_rows(lapply(seq_len(nrow(blocks)), function(b) {
  bind_rows(
    ellipse(c(blocks$s[b], blocks$u[b]), S[[b]], 0.5) %>%
      mutate(level = "50%"),
    ellipse(c(blocks$s[b], blocks$u[b]), S[[b]], 0.95) %>%
      mutate(level = "95%")) %>%
    mutate(block = blocks$block[b], area = blocks$area[b])
}))
line_points <- {
  direction <- c(cos(line$par[1]), sin(line$par[1]))
  normal <- c(-sin(line$par[1]), cos(line$par[1]))
  z <- seq(-10, 10, length.out = 200)
  tibble(s = line$par[2] * normal[1] + z * direction[1],
         u = line$par[2] * normal[2] + z * direction[2])
}
range_s <- range(c(blocks$s - 2.5 * blocks$se_s, blocks$s + 2.5 * blocks$se_s))
range_u <- range(c(blocks$u - 2.5 * blocks$se_u, blocks$u + 2.5 * blocks$se_u))
scatter <- ggplot(mapping = aes(s, u)) +
  geom_polygon(aes(group = interaction(block, level), fill = area,
                   alpha = level), data = ellipses, colour = NA) +
  geom_line(data = line_points, colour = grey(0.4), linetype = 2) +
  geom_point(aes(colour = area, size = bioassays), data = blocks) +
  scale_alpha_manual(values = c(`50%` = 0.22, `95%` = 0.08),
                     name = "estimate\nellipse") +
  scale_colour_manual(values = area_colours, name = "region") +
  scale_fill_manual(values = area_colours, guide = "none") +
  scale_size_area(max_size = 4, name = "bioassays") +
  scale_y_continuous(sec.axis = sec_axis(~ plogis(.),
                                         breaks = c(0.01, 0.05, 0.1, 0.2,
                                                    0.3, 0.5, 0.7),
                                         name = "floor f")) +
  coord_cartesian(xlim = range_s, ylim = range_u) +
  labs(x = "s, log selection multiplier (block vs the base fit)",
       y = "logit floor",
       title = sprintf("Pyrethroid bioassays in %s-degree blocks: selection multiplier and floor",
                       size_label),
       subtitle = sprintf(paste0(
         "%d blocks; between-block correlation %.2f [%.2f, %.2f] (meta-",
         "analysis, bootstrap interval), raw %.2f; one axis (dashed) vs ",
         "free: delta AIC %+.0f"),
         nrow(blocks), correlation$between_correlation,
         correlation$bootstrap_lower, correlation$bootstrap_upper,
         correlation$raw_correlation,
         one_axis$delta_aic_one_axis_minus_free),
       caption = paste(
         "Base: ref_f0 at its posterior mean, held fixed. Per block,",
         "maximum of the beta-binomial log likelihood with priors",
         "s ~ N(0, 1), logit f ~ N(-1.83, 1.39^2); ellipses from the",
         "Hessian.")) +
  theme_bw(base_size = 9) +
  theme(plot.caption = element_text(hjust = 0))
ggsave(file.path(figure_dir,
                 sprintf("floor_selection_blocks_%s.png", size_label)),
       scatter, width = 8.5, height = 6.5, dpi = 200, bg = "white")

borders <- readRDS("data/clean/country_borders.RDS")
# the limits of transmission, without water bodies, as polygons, so that the
# blocks have the only fill scale
water_mask <- sf::st_as_sf(terra::as.polygons(terra::aggregate(
  rast("data/clean/pfpr_water_mask.tif"), 4, fun = "max", na.rm = TRUE)))
squares <- blocks %>%
  mutate(xmin = block_x * block_size, xmax = xmin + block_size,
         ymin = block_y * block_size, ymax = ymin + block_size)
# both centred on the mean across blocks (the meta-analysis mean, mu), so that
# the colour is the block's departure from the typical block; red where
# selection is weaker or the floor higher
ramp <- function(name, centre, half_width, breaks = waiver(),
                 labels = waiver()) {
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "RdBu")),
                       limits = centre + c(-1, 1) * half_width,
                       oob = scales::squish, name = name, breaks = breaks,
                       labels = labels)
}
block_map <- function(fill, scale, title) {
  ggplot() +
    geom_sf(data = borders, fill = grey(0.92), colour = NA) +
    geom_sf(data = water_mask, fill = grey(0.8), colour = NA) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  fill = .data[[fill]]), data = squares, alpha = 0.85,
              colour = "white", linewidth = 0.3) +
    scale +
    geom_sf(data = borders, fill = NA, colour = grey(0.45), linewidth = 0.1) +
    coord_sf(xlim = c(-18, 52), ylim = c(-35, 25), expand = FALSE) +
    labs(title = title) +
    guides(fill = guide_colourbar(barwidth = 12)) +
    theme_ir_maps()
}
squares$negative_s <- -squares$s
half_s <- max(abs(squares$s - correlation$mean_s))
half_u <- max(abs(squares$u - correlation$mean_u))
floor_breaks <- c(0.02, 0.05, 0.2, 0.5, 0.9)
floor_breaks <- floor_breaks[abs(qlogis(floor_breaks) - correlation$mean_u) <=
                               half_u]
maps <- block_map("negative_s",
                  ramp("multiplier", -correlation$mean_s, half_s,
                       labels = function(x) sprintf("%.2g", exp(-x))),
                  "selection multiplier: weaker (red), stronger (blue)") +
  block_map("u",
            ramp("floor", correlation$mean_u, half_u,
                 breaks = qlogis(floor_breaks), labels = floor_breaks),
            "floor: higher (red), lower (blue)") +
  plot_annotation(caption = sprintf(paste(
    "%s-degree blocks of pyrethroid bioassays; colours centred on the mean",
    "across blocks (multiplier %.2f, floor %.2f, relative to the base fit",
    "ref_f0)"), size_label, exp(correlation$mean_s),
    plogis(correlation$mean_u)))
ggsave(file.path(figure_dir,
                 sprintf("floor_selection_blocks_map_%s.png", size_label)),
       maps & theme(legend.position = "bottom"), width = 11, height = 5.5,
       dpi = 200, bg = "white")
report("written; peak memory %.1f GB", peak_memory_gb())
