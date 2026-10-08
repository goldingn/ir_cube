# The spatial range over which pyrethroid plateaus (mortality floors) vary
# (#47): a justification, from the data, for the range of the latent fields of
# stage 1. V5's two HSGP smooths (R/latent_smooth.R) fitted ranges of about
# 500-650 km, far into the tail of their prior, while stage 2 has pixel and
# pixel-year terms for fine-scale residual structure. So: at what scale do the
# plateaus vary broadly, and does their spatial structure have a long-range
# component separable from a short one?
#
#   Rscript R/plateau_range.R blocks <label> [fit file]
#   Rscript R/plateau_range.R levels
#   Rscript R/plateau_range.R pheno
#   Rscript R/plateau_range.R range [pheno] [cached]
#
# blocks  The block screen of R/floor_selection_blocks.R and
#         R/west_block_check.R, with the weighted binomial likelihood
#         (R/weighted_binomial.R). The base is a weighted-binomial fit
#         (outputs/pod_jobs/<label>/temporary/fitted_model.RData unless given;
#         wb_ref, no floor, so that the blocks carry the plateau; wb_bf, one
#         global floor, as a check), every parameter at its posterior mean over
#         the usable chains. At every pyrethroid bioassay it gives, in the outer
#         form, the logit initial state l0, the cumulative log fitness C_t and
#         the reversion t kappa, checked against the fit's own predictions. The
#         bioassays are grouped in square blocks of block_sizes degrees; a block
#         is kept if it has at least min_bioassays pyrethroid bioassays in at
#         least min_years distinct years, at least min_late of them from
#         late_from. In each kept block b, two variants maximise
#           sum_i w_i [y_i log p_i + (n_i - y_i) log(1 - p_i)] + log prior
#           p_i = f_b + (1 - f_b) plogis(l0_i + d_b - exp(s_b) C_i - t_i kappa_i)
#         with w_i = 1 / (1 + (n_i - 1) rho_type) at the replicate rho:
#           a  (s_b, u_b = logit f_b), d_b = 0, as the screen
#           d  (s_b, u_b, d_b), a block shift of the pyrethroids' logit
#              initial state too (variant d of R/west_block_check.R), since a
#              floor fitted with l0 held fixed is fragile; the main estimates
#         priors s ~ N(0, 1), u ~ N(-1.83, 1.39^2) (the logit moments of
#         Beta(1, 4)), d ~ N(0, 1). The covariance: the inverse Hessian of the
#         negative log posterior, and a sandwich with the scores clustered by
#         site (model cell), since bioassays at a site share more than the
#         replicate rho allows; each parameter's variance is the larger of the
#         two. Also, the likelihood-only (prior removed) Gaussian
#         approximation of each estimate, theta_L = m + L^-1 P0 (m - mu0), L =
#         H - P0, for a check that the prior's shrinkage of weakly identified
#         blocks does not drive the result. Writes, in outputs/species_runs/
#         range/, estimates_<label>_<size>.csv (per block: counts, centroid,
#         eligibility, and per variant the estimates) and bioassays_<label>.csv
#         (the pyrethroid bioassays with l0, C_t, t kappa and w).
# levels  The model-free current LLIN-pyrethroid level of each block: a binomial
#         GLMM (lme4) of the LLIN-pyrethroid bioassays (alpha-cypermethrin,
#         deltamethrin, permethrin) of current_years, logit p = block + type +
#         (1 | site) + (1 | bioassay), as the country GLMM of
#         R/diagnostic_regions.R, for blocks with at least min_current
#         bioassays at min_current_sites sites; each block's level at the
#         window's type mix, with its standard error. Writes
#         levels_<size>.csv.
# pheno   The block floors without the dynamical model: no covariates, no
#         country initial states. The blocks stage's blocks and pyrethroid
#         bioassays (bioassays_wb_ref.csv, its raw columns; the weights
#         recomputed from the replicate rho), each block a phenomenological
#         decline to a floor,
#           p_i = f_b + (1 - f_b) plogis(a_b + gamma_type(i)
#                                        - exp(beta_b) (year_i - pheno_origin))
#         f_b = plogis(u_b), with the weighted binomial likelihood and weak
#         priors a ~ N(0, 5^2), beta ~ N(-1, 1.5^2), u ~ N(0, 2.5^2). The type
#         offsets gamma (sum to zero) come from one joint fit of all the kept
#         blocks of pheno_gamma_size degrees, each with its own (a, beta, u),
#         and are then held fixed (the joint fit at the other size, a check);
#         each block's (a, beta, u) is refitted from a grid of starts, with
#         SEs as in the blocks stage (the larger of the Hessian and the
#         site-clustered sandwich). Writes estimates_pheno_<size>.csv (variant
#         "pheno") and pheno_gamma.csv.
# range   For each base, block size and quantity (u_b and s_b of variant d, the
#         current level; as checks u_b of variant a and u_b with the prior
#         removed), with great-circle distances between block centroids (the
#         mean location of the block's bioassays):
#           variogram  the empirical semivariogram in distance bins of
#                      bin_width km to max_distance km, raw, and less each
#                      pair's estimation noise (SE_i^2 + SE_j^2) / 2
#           single     y_b ~ N(mu, sd^2 K(d; range) + tau^2 I + diag(SE_b^2))
#                      by restricted maximum likelihood (REML; mu the only
#                      fixed effect), K squared exponential, the HSGP's kernel,
#                      with range rho = 2 ell as in R/latent_smooth.R, so
#                      exp(-2 (d / rho)^2), correlation 0.135 at d = rho; and
#                      as a check exponential, exp(-2 d / rho), with the same
#                      correlation at rho. The REML log likelihood profiled over
#                      the range on range_grid (200-6,000 km): the maximum and
#                      the approximate 95% interval {2 (l_max - l) <= 3.84}
#           nested     sd_s^2 K(d; r_s) + sd_l^2 K(d; r_l) in place of
#                      sd^2 K(d; range), squared exponential, with the short
#                      range r_s fixed (short_fixed) or profiled over
#                      short_grid, and the long range r_l profiled over
#                      long_grid; against the short component alone and the
#                      single range by likelihood ratio and AIC
#           fixed      the single and nested fits with the (long) range fixed at
#                      fixed_ranges, against their maxima
#         Writes range_fits.csv, variograms.csv, profiles.csv, and figures in
#         figures/species_runs/range/: variograms_<label>.png,
#         profiles_<label>.png and floor_map_<label>.png. With pheno, the
#         same for u_b of the pheno stage alone (base "pheno"), its files
#         suffixed _pheno, and a table against range_fits.csv.
# Plain R; about 2 GB to load a fit; each stage a few minutes.
#
# Caveats: the base fit is held fixed, so its country initial states and
# covariate effects absorb what they can before the blocks are fitted, and the
# block floors are conditional on it; a block's floor, initial state and
# selection trade off (the Hessian covariance carries this); the blocks are
# arbitrary, and ranges shorter than about the block spacing (278 km at 2.5
# degrees, 556 km at 5) are indistinguishable from the nugget, while ranges
# beyond about half the extent of the data (some 3,000-4,000 km) are hard to
# tell from a trend; and the GP measurement model treats each block's
# posterior mode and SE as an unbiased measurement, which the prior removed
# check tests. The pheno floors have no timing from a model: where a block's
# mortality is flat, its decline is placed before the data and its floor at
# the observed level, so they track the current level more than the
# model-conditioned floors do.

arguments <- commandArgs(trailingOnly = TRUE)
stage <- if (length(arguments) >= 1) arguments[1] else "range"
stopifnot(stage %in% c("blocks", "levels", "pheno", "range"))

block_sizes <- c(2.5, 5)
min_bioassays <- 30
min_years <- 5
late_from <- 2016
min_late <- 5
# the priors of R/floor_selection_blocks.R and R/west_block_check.R
prior_mean <- c(s = 0, u = -1.8333, d = 0)
prior_sd <- c(s = 1, u = 1.3888, d = 1)
starts_grid <- list(s = c(-1.5, 0, 1), u = c(-5, -1.8, 0), d = c(-1, 0, 1))
variant_parameters <- list(a = c("s", "u"), d = c("s", "u", "d"))
llin_pyrethroids <- c("Alpha-cypermethrin", "Deltamethrin", "Permethrin")
current_years <- 2019:2024
min_current <- 10
min_current_sites <- 2
# the pheno stage: the time origin, the weak priors and the starts of each
# block's (a, beta, u), and the block size of the joint fit for gamma
pheno_origin <- 2010
pheno_prior_mean <- c(a = 0, beta = -1, u = 0)
pheno_prior_sd <- c(a = 5, beta = 1.5, u = 2.5)
pheno_starts <- list(a = c(1, 3, 5), beta = c(-3, -1.5, 0), u = c(-4, -1.5, 1))
pheno_gamma_size <- 2.5
labels <- c("wb_ref", "wb_bf")
bin_width <- 250
max_distance <- 4000
range_grid <- exp(seq(log(200), log(6000), length.out = 41))
short_fixed <- c(300, 600)
short_grid <- exp(seq(log(200), log(1000), length.out = 9))
long_grid <- exp(seq(log(1000), log(6000), length.out = 25))
fixed_ranges <- c(1000, 1500, 2000, 2500, 3000)
earth_radius <- 6371

output_dir <- "outputs/species_runs/range"
figure_dir <- "figures/species_runs/range"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# greta (for the blocks stage's draws) before dplyr, which must mask it
suppressMessages({
  if (stage == "blocks") library(greta)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})
# report() and peak_memory_gb()
source("R/two_stage_helpers.R")
options(width = 200)

# the region of each record (analysis_region(), R/species_fit_helpers.R),
# Ethiopia apart as in the block screen
block_area <- function(country, region) {
  area <- dplyr::case_when(
    country %in% c("Djibouti", "Eritrea", "Ethiopia", "Somalia",
                   "Sudan") ~ "Horn",
    country %in% c("Comoros", "Madagascar", "Malawi", "Mauritius",
                   "Mozambique", "Zambia", "Zimbabwe") ~ "Southern",
    region %in% c("West", "Western Africa") ~ "West",
    region %in% c("Central", "Middle Africa") ~ "Central",
    region %in% c("East", "Eastern Africa") ~ "East",
    region %in% c("Southern", "Southern Africa") ~ "Southern",
    TRUE ~ "North")
  ifelse(country == "Ethiopia", "Ethiopia", area)
}

# the square block of `size` degrees of each record, and its key
add_blocks <- function(data, size) {
  data %>%
    mutate(block_x = floor(longitude / size),
           block_y = floor(latitude / size),
           block = paste(block_x, block_y))
}

# the most common value
most_common <- function(x) names(which.max(table(x)))

# each block's counts, centroid and eligibility, from pyrethroid bioassays
# with blocks (add_blocks())
summarise_blocks <- function(data) {
  data %>%
    group_by(block, block_x, block_y) %>%
    summarise(bioassays = n(), years = n_distinct(year_start),
              late = sum(year_start >= late_from),
              first_year = min(year_start), last_year = max(year_start),
              sites = n_distinct(cell),
              llin_current = sum(insecticide_type %in% llin_pyrethroids &
                                   year_start %in% current_years),
              longitude = mean(longitude), latitude = mean(latitude),
              area = most_common(area), country = most_common(country_name),
              observed = sum(died) / sum(mosquito_number),
              .groups = "drop") %>%
    mutate(eligible = bioassays >= min_bioassays & years >= min_years &
             late >= min_late)
}

# The log p and log(1 - p) of p = f + (1 - f) q, q = plogis(l), f =
# plogis(u), without losing precision: log p = log(f + (1 - f) q) by
# log-sum-exp of log f and log(1 - f) + log q, so that it stays finite as f
# underflows to 0 (u -> -Inf, the base without a floor) and as q does; log(1
# - p) = log(1 - f) + log(1 - q), each from plogis(log.p = TRUE), never 1
# minus a number near 1
floored_log_probs_logit <- function(l, u) {
  log_f <- plogis(u, log.p = TRUE)
  log_not_f <- plogis(-u, log.p = TRUE)
  log_q <- plogis(l, log.p = TRUE)
  a <- rep_len(log_f, length(l))
  b <- log_not_f + log_q
  larger <- pmax(a, b)
  list(log_p = larger + log1p(exp(pmin(a, b) - larger)),
       log_not_p = log_not_f + plogis(-l, log.p = TRUE))
}


# blocks ---------------------------------------------------------------------------------

if (stage == "blocks") {
  source("R/functions.R")
  source("R/dynamical_predictions.R")
  source("R/species_fit_helpers.R")

  label <- arguments[2]
  stopifnot(!is.na(label))
  fit_file <- if (length(arguments) >= 3) arguments[3] else
    file.path("outputs/pod_jobs", label, "temporary/fitted_model.RData")
  stopifnot(all.equal(logit_beta_moments(c(1, 4))$mean, prior_mean[["u"]],
                      tolerance = 1e-4),
            all.equal(logit_beta_moments(c(1, 4))$sd, prior_sd[["u"]],
                      tolerance = 1e-4))

  # the base fit, at its posterior mean (as R/floor_selection_blocks.R)
  fit <- load_fit(fit_file)
  stopifnot(model_likelihood(fit$options) == "weighted_binomial",
            !species_on(fit$options), !kdr_on(fit$options),
            !smooth_on(fit$options))
  usable <- usable_chains(fit$draws, label)
  means <- colMeans(as.matrix(fit$draws[usable]))
  parameters <- dynamical_parameter_draws(
    list(draws = coda::mcmc.list(coda::mcmc(matrix(
      means, 1, dimnames = list(NULL, names(means))))),
      options = fit$options, x_cells_init = fit$x_cells_init),
    fit$classes_index, fit$types, draw_index = 1, options = fit$options)
  df <- fit$df
  rho <- c(parameters$rho_types)
  stopifnot(isTRUE(all.equal(unname(rho), unname(replicate_rho(fit$types)))))
  kappa <- c(parameters$kappa_type)
  stopifnot(length(kappa) == length(fit$types))
  base_floor <- parameters$mortality_floor
  report("%s: chains %s; base floor %s", label, toString(usable),
         if (is.null(base_floor)) "none" else sprintf("%.3f", base_floor))

  pyrethroids <- df %>%
    filter(insecticide_class == "Pyrethroids")

  # l0 and C at each bioassay, as dynamical_trajectories() forms them (the
  # loop of R/floor_selection_blocks.R)
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
           weight = design_effect_weight(mosquito_number, rho),
           logit_base = l0 - C - reversion)
  # the outer form, with the base's floor, reproduces the fit's own
  # predictions at these parameters
  check <- c(dynamical_logit(parameters, pyrethroids, df, fit$x_cell_years,
                             fit$cell_years_index))
  base_logit <- floored_logit(pyrethroids$logit_base, base_floor)
  difference <- max(abs(pmin(pmax(check, -30), 30) -
                          pmin(pmax(base_logit, -30), 30)))
  report("outer-form quantities at %d pyrethroid bioassays; against the fit's predictions, max |logit diff| %.2g",
         nrow(pyrethroids), difference)
  stopifnot(difference < 1e-8)
  pyrethroids <- pyrethroids %>%
    mutate(area = block_area(country_name, region))
  write.csv(pyrethroids %>%
              select(longitude, latitude, cell, country_name, area,
                     year_start, insecticide_type, died, mosquito_number,
                     l0, C, reversion, rho, weight, logit_base),
            file.path(output_dir, sprintf("bioassays_%s.csv", label)),
            row.names = FALSE)

  # each bioassay's weighted binomial log likelihood at theta (named, of s,
  # u, d; those missing at the base: s = 0, d = 0, the base's floor)
  bioassay_log_lik <- function(data, theta) {
    s <- if ("s" %in% names(theta)) theta[["s"]] else 0
    d <- if ("d" %in% names(theta)) theta[["d"]] else 0
    u <- if ("u" %in% names(theta)) theta[["u"]] else
      if (is.null(base_floor)) -Inf else qlogis(base_floor)
    probs <- floored_log_probs_logit(data$l0 + d - exp(s) * data$C -
                                       data$reversion, u)
    weighted_binomial_log_lik(data$died, data$mosquito_number, probs$log_p,
                              probs$log_not_p, data$weight)
  }
  log_prior <- function(theta) {
    sum(dnorm(theta, prior_mean[names(theta)], prior_sd[names(theta)],
              log = TRUE))
  }

  # one variant in one block: the maximum from a grid of starts; the Hessian
  # and site-clustered sandwich covariances; the prior-removed estimates; the
  # Pearson dispersion
  fit_variant <- function(data, names) {
    starts <- expand.grid(starts_grid[names])
    objective <- function(theta) {
      theta <- setNames(theta, names)
      -(sum(bioassay_log_lik(data, theta)) + log_prior(theta))
    }
    runs <- lapply(seq_len(nrow(starts)), function(i) {
      optim(unlist(starts[i, , drop = FALSE]), objective, method = "BFGS",
            control = list(reltol = 1e-12, maxit = 1000))
    })
    values <- vapply(runs, `[[`, numeric(1), "value")
    best <- runs[[which.min(values)]]
    m <- setNames(best$par, names)
    hessian <- optimHess(m, objective)
    v_hessian <- solve(hessian)
    # the scores of each site's bioassays, by central differences
    step <- 1e-5
    scores <- vapply(names, function(j) {
      up <- m
      down <- m
      up[[j]] <- up[[j]] + step
      down[[j]] <- down[[j]] - step
      (bioassay_log_lik(data, up) - bioassay_log_lik(data, down)) / (2 * step)
    }, numeric(nrow(data)))
    scores <- matrix(scores, nrow(data), dimnames = list(NULL, names))
    site_scores <- rowsum(scores, data$cell)
    sites <- nrow(site_scores)
    meat <- crossprod(site_scores) * sites / max(sites - 1, 1)
    v_sandwich <- v_hessian %*% meat %*% v_hessian
    variance <- pmax(diag(v_hessian), diag(v_sandwich))
    # the likelihood alone, prior removed
    precision_prior <- diag(1 / prior_sd[names] ^ 2, length(names))
    likelihood_precision <- hessian - precision_prior
    positive <- all(eigen(likelihood_precision, symmetric = TRUE,
                          only.values = TRUE)$values > 1e-8)
    if (positive) {
      v_likelihood <- solve(likelihood_precision)
      m_likelihood <- c(m + v_likelihood %*% precision_prior %*%
                          (m - prior_mean[names]))
    } else {
      v_likelihood <- matrix(NA_real_, length(names), length(names))
      m_likelihood <- rep(NA_real_, length(names))
    }
    probs <- floored_log_probs_logit(
      data$l0 + (if ("d" %in% names) m[["d"]] else 0) -
        exp(m[["s"]]) * data$C - data$reversion, m[["u"]])
    p <- exp(probs$log_p)
    pearson <- sum(data$weight * (data$died - data$mosquito_number * p) ^ 2 /
                     (data$mosquito_number * p * (1 - p))) /
      (nrow(data) - length(names))
    out <- tibble(
      log_post = -best$value,
      log_lik = sum(bioassay_log_lik(data, m)),
      converged = best$convergence == 0,
      start_spread = max(values) - min(values),
      clusters = sites,
      dispersion = pearson)
    for (j in seq_along(names)) {
      name <- names[j]
      out[[name]] <- m[[name]]
      out[[paste0("se_", name)]] <- sqrt(variance[j])
      out[[paste0("se_hessian_", name)]] <- sqrt(v_hessian[j, j])
      out[[paste0("se_sandwich_", name)]] <- sqrt(v_sandwich[j, j])
      out[[paste0(name, "_likelihood")]] <- m_likelihood[j]
      out[[paste0("se_likelihood_", name)]] <- sqrt(v_likelihood[j, j])
      out[[paste0("shrinkage_", name)]] <- v_hessian[j, j] / prior_sd[[name]] ^ 2
    }
    out$cor_su <- v_hessian[1, 2] / sqrt(v_hessian[1, 1] * v_hessian[2, 2])
    if ("d" %in% names) {
      out$cor_ud <- v_hessian[2, 3] / sqrt(v_hessian[2, 2] * v_hessian[3, 3])
      out$cor_sd <- v_hessian[1, 3] / sqrt(v_hessian[1, 1] * v_hessian[3, 3])
    }
    out
  }

  for (size in block_sizes) {
    data <- add_blocks(pyrethroids, size)
    blocks <- summarise_blocks(data)
    kept <- filter(blocks, eligible)
    report("%s, %s-degree blocks: %d with pyrethroid bioassays, %d kept (>= %d bioassays in >= %d years, >= %d from %d), holding %d of %d bioassays",
           label, format(size), nrow(blocks), nrow(kept), min_bioassays,
           min_years, min_late, late_from, sum(kept$bioassays), nrow(data))
    by_block <- split(data, data$block)
    time <- system.time(
      estimates <- bind_rows(lapply(kept$block, function(b) {
        bind_rows(lapply(names(variant_parameters), function(v) {
          fit_variant(by_block[[b]], variant_parameters[[v]]) %>%
            mutate(block = b, variant = v, .before = 1)
        }))
      }))
    )[["elapsed"]]
    base_lik <- vapply(kept$block, function(b) {
      sum(bioassay_log_lik(by_block[[b]], c()))
    }, numeric(1))
    estimates <- estimates %>%
      left_join(tibble(block = kept$block, log_lik_base = base_lik),
                by = "block") %>%
      mutate(f = plogis(u), f_lower = plogis(u - 1.96 * se_u),
             f_upper = plogis(u + 1.96 * se_u), multiplier = exp(s),
             gain = log_lik - log_lik_base)
    report("fitted %d blocks x %d variants in %.0f s; %d not converged; start spread > 0.01 in %d",
           nrow(kept), length(variant_parameters), time,
           sum(!estimates$converged), sum(estimates$start_spread > 0.01))
    out <- blocks %>%
      left_join(estimates, by = "block") %>%
      mutate(block_size = size, base = label, .before = 1)
    write.csv(out, file.path(output_dir, sprintf("estimates_%s_%s.csv", label,
                                                 format(size))),
              row.names = FALSE)
    print(as.data.frame(out %>%
                          filter(eligible, variant == "d") %>%
                          transmute(block, area, country, bioassays, years,
                                    late, sites, s, se_s, u, se_u,
                                    se_hessian_u, se_sandwich_u, f, d, se_d,
                                    u_likelihood, se_likelihood_u,
                                    shrinkage_u, cor_su, cor_ud,
                                    dispersion, gain)),
          digits = 3)
  }
  report("peak memory %.1f GB", peak_memory_gb())
}


# levels -----------------------------------------------------------------------------

if (stage == "levels") {
  suppressMessages(library(lme4))
  bioassays <- read.csv(file.path(output_dir, "bioassays_wb_ref.csv")) %>%
    filter(insecticide_type %in% llin_pyrethroids,
           year_start %in% current_years) %>%
    mutate(survived = mosquito_number - died, site = factor(cell),
           bioassay = factor(seq_len(n())), type_f = factor(insecticide_type))
  contrasts(bioassays$type_f) <- contr.sum(3)
  report("current LLIN-pyrethroid bioassays: %d at %d sites",
         nrow(bioassays), n_distinct(bioassays$cell))
  for (size in block_sizes) {
    data <- add_blocks(bioassays, size)
    counts <- data %>%
      group_by(block, block_x, block_y) %>%
      summarise(bioassays = n(), sites = n_distinct(cell),
                tested = sum(mosquito_number),
                observed = sum(died) / sum(mosquito_number),
                longitude = mean(longitude), latitude = mean(latitude),
                area = most_common(area), country = most_common(country_name),
                .groups = "drop") %>%
      mutate(eligible = bioassays >= min_current & sites >= min_current_sites)
    data <- data %>%
      filter(block %in% counts$block[counts$eligible]) %>%
      mutate(block_f = factor(block), site = factor(cell),
             bioassay = factor(seq_len(n())),
             type_f = factor(insecticide_type))
    contrasts(data$type_f) <- contr.sum(3)
    time <- system.time(
      glmm <- glmer(cbind(died, survived) ~ 0 + block_f + type_f +
                      (1 | site) + (1 | bioassay),
                    family = binomial, data = data,
                    control = glmerControl(optimizer = "bobyqa",
                                           optCtrl = list(maxfun = 1e5)))
    )[["elapsed"]]
    variances <- as.data.frame(VarCorr(glmm))
    beta <- fixef(glmm)
    type_levels <- levels(data$type_f)
    type_share <- tapply(data$mosquito_number, data$type_f, sum)[type_levels]
    type_share <- type_share / sum(type_share)
    block_names <- levels(data$block_f)
    contrast <- matrix(0, length(block_names), length(beta),
                       dimnames = list(block_names, names(beta)))
    for (i in seq_along(block_names)) {
      contrast[i, paste0("block_f", block_names[i])] <- 1
    }
    contrast[, "type_f1"] <- type_share[1] - type_share[3]
    contrast[, "type_f2"] <- type_share[2] - type_share[3]
    levels <- tibble(block = block_names,
                     level = c(contrast %*% beta),
                     se_level = sqrt(diag(contrast %*% as.matrix(vcov(glmm)) %*%
                                            t(contrast))))
    report("%s-degree blocks: %d with current bioassays, %d kept (>= %d bioassays at >= %d sites); GLMM in %.0f s, sd site %.2f, sd bioassay %.2f%s",
           format(size), nrow(counts), sum(counts$eligible), min_current,
           min_current_sites, time,
           variances$sdcor[variances$grp == "site"],
           variances$sdcor[variances$grp == "bioassay"],
           if (is.null(glmm@optinfo$conv$lme4$messages)) "" else
             paste(";", toString(glmm@optinfo$conv$lme4$messages)))
    out <- counts %>%
      left_join(levels, by = "block") %>%
      mutate(block_size = size, level_p = plogis(level), .before = 1)
    write.csv(out, file.path(output_dir, sprintf("levels_%s.csv", format(size))),
              row.names = FALSE)
    print(as.data.frame(filter(out, eligible) %>%
                          select(block, area, country, bioassays, sites,
                                 observed, level_p, level, se_level)),
          digits = 3)
  }
}


# pheno ------------------------------------------------------------------------------

if (stage == "pheno") {
  source("R/weighted_binomial.R")
  pyrethroids <- read.csv(file.path(output_dir, "bioassays_wb_ref.csv")) %>%
    select(longitude, latitude, cell, country_name, area, year_start,
           insecticide_type, died, mosquito_number, weight_file = weight)
  types <- sort(unique(pyrethroids$insecticide_type))
  rho <- replicate_rho(types)
  pyrethroids <- pyrethroids %>%
    mutate(type_id = match(insecticide_type, types),
           time = year_start - pheno_origin,
           weight = design_effect_weight(mosquito_number, rho[type_id]))
  stopifnot(max(abs(pyrethroids$weight - pyrethroids$weight_file)) < 1e-12)
  report("%d pyrethroid bioassays at %d sites, %d-%d; replicate rho %s",
         nrow(pyrethroids), n_distinct(pyrethroids$cell),
         min(pyrethroids$year_start), max(pyrethroids$year_start),
         toString(sprintf("%s %.3f", types, rho)))
  # the type offsets gamma = type_contrast g sum to zero
  type_contrast <- contr.sum(length(types))

  # Each bioassay's weighted binomial log likelihood at logit l and logit
  # floor u (each per bioassay), and its derivatives in l and u,
  #   d log p / dl = (1 - f) q (1 - q) / p,  d log(1 - p) / dl = -q
  #   d log p / du = f (1 - f) (1 - q) / p,  d log(1 - p) / du = -f
  # the ratios formed on the log scale (p >= f, so log p is finite)
  pheno_log_lik <- function(data, l, u) {
    probs <- floored_log_probs_logit(l, u)
    log_f <- plogis(u, log.p = TRUE)
    log_not_f <- plogis(-u, log.p = TRUE)
    log_q <- plogis(l, log.p = TRUE)
    log_not_q <- plogis(-l, log.p = TRUE)
    y <- data$died
    n <- data$mosquito_number
    w <- data$weight
    list(log_lik = weighted_binomial_log_lik(y, n, probs$log_p,
                                             probs$log_not_p, w),
         d_l = w * (y * exp(log_not_f + log_q + log_not_q - probs$log_p) -
                      (n - y) * exp(log_q)),
         d_u = w * (y * exp(log_f + log_not_f + log_not_q - probs$log_p) -
                      (n - y) * exp(log_f)),
         p = exp(probs$log_p), p_not = exp(probs$log_not_p))
  }

  # The negative log posterior of blocks' a, beta and u (vectors over the
  # blocks data$block_id = 1, ..., B, all present) at type offsets gamma,
  # with its gradient in each and in gamma
  pheno_posterior <- function(data, a, beta, u, gamma) {
    slope <- exp(beta)[data$block_id]
    terms <- pheno_log_lik(data, a[data$block_id] + gamma[data$type_id] -
                             slope * data$time, u[data$block_id])
    by_block <- function(x) as.vector(rowsum(x, data$block_id, reorder = TRUE))
    prior <- function(x, name) {
      list(log = sum(dnorm(x, pheno_prior_mean[[name]],
                           pheno_prior_sd[[name]], log = TRUE)),
           d = -(x - pheno_prior_mean[[name]]) / pheno_prior_sd[[name]] ^ 2)
    }
    prior_a <- prior(a, "a")
    prior_beta <- prior(beta, "beta")
    prior_u <- prior(u, "u")
    list(value = -(sum(terms$log_lik) + prior_a$log + prior_beta$log +
                     prior_u$log),
         d_a = -(by_block(terms$d_l) + prior_a$d),
         d_beta = -(by_block(-terms$d_l * slope * data$time) + prior_beta$d),
         d_u = -(by_block(terms$d_u) + prior_u$d),
         d_gamma = -vapply(seq_along(gamma), function(k) {
           sum(terms$d_l[data$type_id == k])
         }, numeric(1)),
         terms = terms, slope = slope)
  }

  # One block at type offsets gamma: the maximum from the grid of starts
  # (and `start`, if given); the Hessian and site-clustered sandwich
  # covariances, as fit_variant() of the blocks stage; the Pearson
  # dispersion; the midpoint year of the decline and the share of it left
  # at the block's last year (at the mean type, gamma 0)
  pheno_block <- function(data, gamma, start = NULL) {
    data$block_id <- 1L
    names <- c("a", "beta", "u")
    objective <- function(theta) {
      pheno_posterior(data, theta[1], theta[2], theta[3], gamma)$value
    }
    gradient <- function(theta) {
      x <- pheno_posterior(data, theta[1], theta[2], theta[3], gamma)
      c(x$d_a, x$d_beta, x$d_u)
    }
    starts <- as.matrix(expand.grid(pheno_starts))
    if (!is.null(start)) starts <- rbind(start, starts)
    runs <- lapply(seq_len(nrow(starts)), function(i) {
      optim(starts[i, ], objective, gradient, method = "BFGS",
            control = list(reltol = 1e-12, maxit = 2000))
    })
    values <- vapply(runs, `[[`, numeric(1), "value")
    best <- runs[[which.min(values)]]
    m <- setNames(best$par, names)
    hessian <- optimHess(m, objective, gradient)
    v_hessian <- solve(hessian)
    at <- pheno_posterior(data, m[["a"]], m[["beta"]], m[["u"]], gamma)
    scores <- cbind(at$terms$d_l, -at$terms$d_l * at$slope * data$time,
                    at$terms$d_u)
    site_scores <- rowsum(scores, data$cell)
    sites <- nrow(site_scores)
    meat <- crossprod(site_scores) * sites / max(sites - 1, 1)
    v_sandwich <- v_hessian %*% meat %*% v_hessian
    variance <- pmax(diag(v_hessian), diag(v_sandwich))
    p <- at$terms$p
    pearson <- sum(data$weight * (data$died - data$mosquito_number * p) ^ 2 /
                     (data$mosquito_number * p * at$terms$p_not)) /
      (nrow(data) - length(names))
    slope <- exp(m[["beta"]])
    out <- tibble(
      log_post = -best$value,
      log_lik = sum(at$terms$log_lik),
      converged = best$convergence == 0,
      start_spread = max(values) - min(values),
      gradient_max = max(abs(gradient(m))),
      clusters = sites,
      dispersion = pearson,
      slope = slope,
      midpoint = pheno_origin + m[["a"]] / slope,
      left_last = plogis(m[["a"]] - slope * (max(data$year_start) -
                                               pheno_origin)))
    for (j in seq_along(names)) {
      name <- names[j]
      out[[name]] <- m[[name]]
      out[[paste0("se_", name)]] <- sqrt(variance[j])
      out[[paste0("se_hessian_", name)]] <- sqrt(v_hessian[j, j])
      out[[paste0("se_sandwich_", name)]] <- sqrt(v_sandwich[j, j])
      out[[paste0("shrinkage_", name)]] <- v_hessian[j, j] /
        pheno_prior_sd[[name]] ^ 2
    }
    correlation <- cov2cor(v_hessian)
    out %>%
      mutate(cor_a_beta = correlation[1, 2], cor_a_u = correlation[1, 3],
             cor_beta_u = correlation[2, 3])
  }

  # The joint fit of all blocks (data$block_id = 1, ..., B), each with its
  # own (a, beta, u), and the type offsets gamma, from `starts` (B x 3) and
  # gamma 0; gamma's covariance from the inverse of the joint Hessian
  pheno_joint <- function(data, starts) {
    B <- nrow(starts)
    K <- length(types)
    unpack <- function(theta) {
      list(a = theta[seq_len(B)], beta = theta[B + seq_len(B)],
           u = theta[2 * B + seq_len(B)],
           gamma = c(type_contrast %*% theta[3 * B + seq_len(K - 1)]))
    }
    objective <- function(theta) {
      x <- unpack(theta)
      pheno_posterior(data, x$a, x$beta, x$u, x$gamma)$value
    }
    gradient <- function(theta) {
      x <- unpack(theta)
      g <- pheno_posterior(data, x$a, x$beta, x$u, x$gamma)
      c(g$d_a, g$d_beta, g$d_u, c(crossprod(type_contrast, g$d_gamma)))
    }
    theta <- c(starts[, 1], starts[, 2], starts[, 3], rep(0, K - 1))
    fit <- optim(theta, objective, gradient, method = "BFGS",
                 control = list(reltol = 1e-14, maxit = 20000))
    hessian <- optimHess(fit$par, objective, gradient)
    index <- 3 * B + seq_len(K - 1)
    v_g <- solve(hessian)[index, index]
    v_gamma <- type_contrast %*% v_g %*% t(type_contrast)
    x <- unpack(fit$par)
    list(gamma = setNames(x$gamma, types),
         se_gamma = setNames(sqrt(diag(v_gamma)), types),
         value = fit$value, convergence = fit$convergence,
         iterations = fit$counts[["function"]],
         gradient_max = max(abs(gradient(fit$par))),
         start_value = objective(theta),
         blocks = cbind(a = x$a, beta = x$beta, u = x$u))
  }

  # the kept blocks at each size (the blocks stage's), each block's starts
  # alone at gamma 0, and the joint fit for gamma; the joint fit's size
  # first, as its gamma is used at both
  sizes <- c(pheno_gamma_size, setdiff(block_sizes, pheno_gamma_size))
  prepared <- list()
  gammas <- list()
  for (size in sizes) {
    data <- add_blocks(pyrethroids, size)
    blocks <- summarise_blocks(data)
    kept <- filter(blocks, eligible)
    base_blocks <- read.csv(file.path(output_dir, sprintf(
      "estimates_wb_ref_%s.csv", format(size)))) %>%
      filter(eligible, variant == "d")
    stopifnot(setequal(kept$block, base_blocks$block))
    data <- data %>%
      filter(block %in% kept$block) %>%
      mutate(block_id = match(block, kept$block))
    by_block <- split(data, data$block_id)
    time <- system.time({
      initial <- t(vapply(by_block, function(d) {
        x <- pheno_block(d, rep(0, length(types)))
        c(a = x$a, beta = x$beta, u = x$u)
      }, numeric(3)))
      joint <- pheno_joint(data, initial)
    })[["elapsed"]]
    report("%s-degree blocks: %d kept (as the blocks stage), %d bioassays; joint fit for gamma in %.0f s (%d evaluations, convergence %d, max |gradient| %.2g, -log posterior %.1f from %.1f): gamma %s",
           format(size), nrow(kept), nrow(data), time, joint$iterations,
           joint$convergence, joint$gradient_max, joint$value,
           joint$start_value,
           toString(sprintf("%s %.3f (%.3f)", types, joint$gamma,
                            joint$se_gamma)))
    prepared[[format(size)]] <- list(data = data, blocks = blocks,
                                     kept = kept, by_block = by_block,
                                     joint = joint, base = base_blocks)
    gammas[[format(size)]] <- tibble(block_size = size, type = types,
                                     rho = unname(rho), gamma = joint$gamma,
                                     se_gamma = joint$se_gamma,
                                     used = size == pheno_gamma_size)
  }
  gamma <- prepared[[format(pheno_gamma_size)]]$joint$gamma
  write.csv(bind_rows(gammas), file.path(output_dir, "pheno_gamma.csv"),
            row.names = FALSE)

  # each block at gamma, refitted from the grid and its joint-fit maximum
  # (at the joint fit's size) or its maximum alone at gamma 0
  for (size in block_sizes) {
    x <- prepared[[format(size)]]
    time <- system.time(
      estimates <- bind_rows(lapply(seq_along(x$by_block), function(i) {
        pheno_block(x$by_block[[i]], gamma, x$joint$blocks[i, ]) %>%
          mutate(block = x$kept$block[i], .before = 1)
      }))
    )[["elapsed"]]
    estimates <- estimates %>%
      mutate(f = plogis(u), f_lower = plogis(u - 1.96 * se_u),
             f_upper = plogis(u + 1.96 * se_u))
    if (size == pheno_gamma_size) {
      # at the joint fit's gamma, the blocks' maxima can only match or
      # improve on the joint one
      stopifnot(-sum(estimates$log_post) <= x$joint$value + 1e-4)
    }
    out <- x$blocks %>%
      left_join(estimates, by = "block") %>%
      mutate(block_size = size, base = "pheno", variant = "pheno",
             .before = 1)
    write.csv(out, file.path(output_dir, sprintf("estimates_pheno_%s.csv",
                                                 format(size))),
              row.names = FALSE)
    report("%s-degree blocks: %d fitted in %.0f s; %d not converged, start spread > 0.01 in %d, max |gradient| %.2g; floors prior-dominated (SE of u > 1) %d; SE of u median %.2f; the sandwich SE the larger in %d; decline more than half left at the last year (left_last > 0.5) %d, under 0.1 %d",
           format(size), nrow(estimates), time, sum(!estimates$converged),
           sum(estimates$start_spread > 0.01),
           max(estimates$gradient_max), sum(estimates$se_u > 1),
           median(estimates$se_u),
           sum(estimates$se_sandwich_u > estimates$se_hessian_u),
           sum(estimates$left_last > 0.5), sum(estimates$left_last < 0.1))
    # against the model-conditioned floors (variant d, wb_ref) and the
    # current level, block by block
    compare <- estimates %>%
      inner_join(x$base %>% select(block, u_model = u, se_u_model = se_u),
                 by = "block")
    well <- compare$se_u < 1 & compare$se_u_model < 1
    levels <- read.csv(file.path(output_dir, sprintf("levels_%s.csv",
                                                     format(size)))) %>%
      filter(eligible) %>%
      inner_join(estimates %>% select(block, u, se_u), by = "block")
    report("u against the model-conditioned floors (wb_ref, d): correlation %.2f over %d blocks, %.2f over the %d with both SEs < 1 (mean difference %.2f); against the current level: correlation %.2f over %d blocks",
           cor(compare$u, compare$u_model), nrow(compare),
           cor(compare$u[well], compare$u_model[well]), sum(well),
           mean(compare$u[well] - compare$u_model[well]),
           cor(levels$u, levels$level), nrow(levels))
    print(as.data.frame(out %>%
                          filter(eligible) %>%
                          transmute(block, area, country, bioassays, years,
                                    first_year, last_year, sites, observed,
                                    a, slope, midpoint, left_last, u, se_u,
                                    se_hessian_u, se_sandwich_u, f,
                                    cor_beta_u, cor_a_u, dispersion)),
          digits = 3)
  }
}


# range ------------------------------------------------------------------------------

if (stage == "range") {
  suppressMessages({
    library(ggplot2)
    library(patchwork)
    library(sf)
    library(terra)
    library(ggtext)
  })
  source("R/functions.R")

  # the set: the model-conditioned floors with the current level, or with
  # pheno the phenomenological floors alone (files suffixed _pheno)
  pheno <- "pheno" %in% arguments[-1]
  suffix <- if (pheno) "_pheno" else ""
  # the quantities: the main ones, then the checks
  if (pheno) {
    quantities <- c(u = "logit floor u (phenomenological)")
    main_quantities <- "u"
    labels <- "pheno"
  } else {
    quantities <- c(u = "logit floor u (d)",
                    s = "log selection s (d)",
                    level = "current level 2019-24",
                    u_a = "check: u (a, l0 fixed)",
                    u_lik = "check: u (d, no prior)")
    main_quantities <- c("u", "s", "level")
  }

  # every (base, block size, quantity): the estimates, their SEs and the
  # block centroids
  dataset <- function(base, size, quantity, data, y, se) {
    keep <- is.finite(y) & is.finite(se)
    tibble(base = base, block_size = size, quantity = quantity,
           block = data$block[keep], longitude = data$longitude[keep],
           latitude = data$latitude[keep], area = data$area[keep],
           y = y[keep], se = se[keep])
  }
  data <- list()
  if (pheno) {
    for (size in block_sizes) {
      estimates <- read.csv(file.path(output_dir, sprintf(
        "estimates_pheno_%s.csv", format(size)))) %>%
        filter(eligible)
      data <- c(data, list(dataset("pheno", size, "u", estimates,
                                   estimates$u, estimates$se_u)))
    }
  }
  for (label in setdiff(labels, "pheno")) {
    for (size in block_sizes) {
      file <- file.path(output_dir, sprintf("estimates_%s_%s.csv", label,
                                            format(size)))
      if (!file.exists(file)) next
      estimates <- read.csv(file) %>% filter(eligible)
      d <- filter(estimates, variant == "d")
      a <- filter(estimates, variant == "a")
      data <- c(data, list(
        dataset(label, size, "u", d, d$u, d$se_u),
        dataset(label, size, "s", d, d$s, d$se_s),
        dataset(label, size, "u_a", a, a$u, a$se_u),
        dataset(label, size, "u_lik", d, d$u_likelihood, d$se_likelihood_u)))
    }
  }
  # the current level, not with pheno
  for (size in if (pheno) NULL else block_sizes) {
    levels <- read.csv(file.path(output_dir, sprintf("levels_%s.csv",
                                                     format(size)))) %>%
      filter(eligible)
    for (label in labels) {
      data <- c(data, list(dataset(label, size, "level", levels,
                                   levels$level, levels$se_level)))
    }
  }
  data <- bind_rows(data)

  # great-circle distances (km) between points
  great_circle <- function(longitude, latitude) {
    phi <- latitude * pi / 180
    lambda <- longitude * pi / 180
    a <- sin(outer(phi, phi, "-") / 2) ^ 2 +
      outer(cos(phi), cos(phi)) * sin(outer(lambda, lambda, "-") / 2) ^ 2
    # pmin() takes its attributes from its first argument: the matrix first
    2 * earth_radius * asin(pmin(sqrt(a), 1))
  }
  # the correlation functions, with range r the distance at which the
  # correlation is exp(-2) = 0.135: the squared exponential of the HSGP
  # (R/latent_smooth.R: lengthscale ell = r / 2, exp(-d^2 / (2 ell^2))), and
  # the exponential
  kernels <- list(se = function(d, r) exp(-2 * (d / r) ^ 2),
                  exp = function(d, r) exp(-2 * d / r))

  # The REML log likelihood (up to a constant) of y ~ N(mu 1, V), mu the only
  # fixed effect: -1/2 [log|V| + log(1' V^-1 1) + r' V^-1 r], r = y - mu_hat,
  # mu_hat the GLS mean; by the Cholesky factor
  reml_log_lik <- function(y, V) {
    R <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(R)) return(-1e10)
    z_y <- backsolve(R, y, transpose = TRUE)
    z_1 <- backsolve(R, rep(1, length(y)), transpose = TRUE)
    a <- sum(z_1 ^ 2)
    mu <- sum(z_1 * z_y) / a
    -sum(log(diag(R))) - 0.5 * log(a) - 0.5 * sum((z_y - mu * z_1) ^ 2)
  }
  gls_mean <- function(y, V) {
    R <- chol(V)
    z_y <- backsolve(R, y, transpose = TRUE)
    z_1 <- backsolve(R, rep(1, length(y)), transpose = TRUE)
    sum(z_1 * z_y) / sum(z_1 ^ 2)
  }

  # the REML maximum over the sds of the structured components (correlation
  # matrices `Ks`) and the nugget tau, with the estimation noise se2 known:
  # V = sum_j sd_j^2 K_j + tau^2 I + diag(se2); log sds bounded to
  # [1e-3, 20]
  fit_variances <- function(y, se2, Ks, start = NULL) {
    k <- length(Ks)
    n <- length(y)
    objective <- function(par) {
      V <- diag(se2 + exp(2 * par[k + 1]), n)
      for (j in seq_len(k)) V <- V + exp(2 * par[j]) * Ks[[j]]
      -reml_log_lik(y, V)
    }
    signal <- max(var(y) - mean(se2), 0.05 * var(y))
    shares <- if (k == 0) list(1) else if (k == 1) {
      list(c(0.8, 0.2), c(0.2, 0.8))
    } else {
      list(c(0.4, 0.4, 0.2), c(0.1, 0.7, 0.2), c(0.7, 0.1, 0.2))
    }
    starts <- lapply(shares, function(share) 0.5 * log(signal * share))
    if (!is.null(start)) starts <- c(list(start), starts)
    runs <- lapply(starts, function(start) {
      optim(start, objective, method = "L-BFGS-B",
            lower = rep(log(1e-3), k + 1), upper = rep(log(20), k + 1))
    })
    best <- runs[[which.min(vapply(runs, `[[`, numeric(1), "value"))]]
    list(log_lik = -best$value, sd = exp(best$par[seq_len(k)]),
         tau = exp(best$par[k + 1]), par = best$par)
  }

  # the approximate 95% interval of a profile (values `r`, log likelihoods
  # `ll`) about its maximum ll_max: where 2 (ll_max - ll) <= 3.84, the ends
  # interpolated on log r; NA at an end that reaches the edge of the grid
  profile_interval <- function(r, ll, ll_max) {
    deviance <- 2 * (ll_max - ll)
    cut <- qchisq(0.95, 1)
    inside <- which(deviance <= cut)
    if (length(inside) == 0) return(c(NA_real_, NA_real_))
    lo <- min(inside)
    hi <- max(inside)
    crossing <- function(i, j) {
      exp(approx(deviance[c(i, j)], log(r[c(i, j)]), xout = cut)$y)
    }
    c(if (lo == 1) NA_real_ else crossing(lo - 1, lo),
      if (hi == length(r)) NA_real_ else crossing(hi, hi + 1))
  }

  # the profile of a one-range model over `grid`, refined at the maximum by
  # a 1-D search on log r between the grid's neighbours of the best point
  profile_range <- function(y, se2, D, grid, components) {
    fits <- list()
    start <- NULL
    for (i in seq_along(grid)) {
      fits[[i]] <- fit_variances(y, se2, components(grid[i]), start)
      start <- fits[[i]]$par
    }
    ll <- vapply(fits, `[[`, numeric(1), "log_lik")
    best <- which.max(ll)
    bounds <- log(grid[c(max(best - 1, 1), min(best + 1, length(grid)))])
    refined <- optimize(function(log_r) {
      fit_variances(y, se2, components(exp(log_r)), fits[[best]]$par)$log_lik
    }, bounds, maximum = TRUE, tol = 1e-3)
    if (refined$objective > ll[best]) {
      r_hat <- exp(refined$maximum)
      ll_max <- refined$objective
    } else {
      r_hat <- grid[best]
      ll_max <- ll[best]
    }
    fit_hat <- fit_variances(y, se2, components(r_hat), fits[[best]]$par)
    interval <- profile_interval(grid, ll, ll_max)
    list(grid = grid, ll = ll, r_hat = r_hat, ll_max = ll_max,
         lower = interval[1], upper = interval[2], fit = fit_hat,
         at_edge = best %in% c(1, length(grid)))
  }

  analyse <- function(set) {
    y <- set$y
    se2 <- set$se ^ 2
    D <- great_circle(set$longitude, set$latitude)
    nugget <- fit_variances(y, se2, list())
    single <- lapply(kernels, function(kernel) {
      profile_range(y, se2, D, range_grid,
                    function(r) list(kernel(D, r)))
    })
    short_only <- lapply(short_fixed, function(r_s) {
      fit_variances(y, se2, list(kernels$se(D, r_s)))
    })
    nested_fixed <- lapply(short_fixed, function(r_s) {
      K_s <- kernels$se(D, r_s)
      profile_range(y, se2, D, long_grid,
                    function(r) list(K_s, kernels$se(D, r)))
    })
    # both ranges: r_s on short_grid, r_l on long_grid; the profile of r_l
    # maximises over r_s
    K_short <- lapply(short_grid, function(r) kernels$se(D, r))
    both <- matrix(NA_real_, length(short_grid), length(long_grid))
    both_fits <- list()
    for (j in seq_along(long_grid)) {
      K_l <- kernels$se(D, long_grid[j])
      start <- NULL
      for (i in seq_along(short_grid)) {
        fit <- fit_variances(y, se2, list(K_short[[i]], K_l), start)
        start <- fit$par
        both[i, j] <- fit$log_lik
        both_fits[[paste(i, j)]] <- fit
      }
    }
    best_both <- which(both == max(both), arr.ind = TRUE)[1, ]
    ll_both <- max(both)
    profile_long <- apply(both, 2, max)
    interval_long <- profile_interval(long_grid, profile_long, ll_both)
    fit_both <- both_fits[[paste(best_both[1], best_both[2])]]
    # the short component alone, its range profiled on short_grid: the
    # null for a long component
    short_profile <- vapply(K_short, function(K_s) {
      fit_variances(y, se2, list(K_s))$log_lik
    }, numeric(1))
    # fixed (long) ranges: the single range, and the long range with r_s
    # profiled on short_grid, with the nested fit's sds
    fixed <- bind_rows(lapply(fixed_ranges, function(r) {
      K_l <- kernels$se(D, r)
      nested_fits <- lapply(K_short, function(K_s) {
        fit_variances(y, se2, list(K_s, K_l))
      })
      ll <- vapply(nested_fits, `[[`, numeric(1), "log_lik")
      best <- nested_fits[[which.max(ll)]]
      single_fit <- fit_variances(y, se2, list(K_l))
      tibble(range = r, ll_single = single_fit$log_lik,
             sd_single = single_fit$sd, tau_single = single_fit$tau,
             ll_nested = max(ll), short_nested = short_grid[which.max(ll)],
             sd_short_nested = best$sd[1], sd_long_nested = best$sd[2],
             tau_nested = best$tau)
    }))
    list(set = set, D = D, nugget = nugget, single = single,
         short_only = short_only, nested_fixed = nested_fixed,
         both = list(ll = ll_both, r_s = short_grid[best_both[1]],
                     r_l = long_grid[best_both[2]], fit = fit_both,
                     lower = interval_long[1], upper = interval_long[2],
                     profile = profile_long,
                     at_edge = best_both[2] %in% c(1, length(long_grid))),
         fixed = fixed,
         short = list(ll = max(short_profile),
                      r_s = short_grid[which.max(short_profile)]),
         mean = gls_mean(y, diag(se2 + nugget$tau ^ 2, length(y))))
  }

  groups <- data %>%
    group_by(base, block_size, quantity) %>%
    group_split()
  # the analyses, cached in results.rds for `range cached` (the figures and
  # tables again without refitting)
  cache <- file.path(output_dir, sprintf("results%s.rds", suffix))
  if ("cached" %in% arguments[-1]) {
    results <- readRDS(cache)
    report("%d analyses from %s", length(results), cache)
  } else {
    time <- system.time(
      results <- lapply(groups, analyse)
    )[["elapsed"]]
    saveRDS(results, cache)
    report("%d analyses in %.0f s", length(results), time)
  }

  # the table of fits
  aic <- function(ll, k) -2 * ll + 2 * k
  fits <- bind_rows(lapply(results, function(x) {
    set <- x$set
    se <- x$single$se
    ex <- x$single$exp
    out <- tibble(
      base = set$base[1], block_size = set$block_size[1],
      quantity = set$quantity[1], blocks = nrow(set),
      sd_y = sd(set$y), median_se = median(set$se),
      ll_nugget = x$nugget$log_lik, tau_nugget = x$nugget$tau,
      ll_se = se$ll_max, range_se = se$r_hat, lower_se = se$lower,
      upper_se = se$upper, sd_se = se$fit$sd, tau_se = se$fit$tau,
      ll_exp = ex$ll_max, range_exp = ex$r_hat, lower_exp = ex$lower,
      upper_exp = ex$upper, sd_exp = ex$fit$sd, tau_exp = ex$fit$tau)
    for (i in seq_along(short_fixed)) {
      tag <- format(short_fixed[i])
      nf <- x$nested_fixed[[i]]
      out[[paste0("ll_short_", tag)]] <- x$short_only[[i]]$log_lik
      out[[paste0("ll_nested_", tag)]] <- nf$ll_max
      out[[paste0("long_nested_", tag)]] <- nf$r_hat
      out[[paste0("lower_nested_", tag)]] <- nf$lower
      out[[paste0("upper_nested_", tag)]] <- nf$upper
      out[[paste0("sd_short_nested_", tag)]] <- nf$fit$sd[1]
      out[[paste0("sd_long_nested_", tag)]] <- nf$fit$sd[2]
      out[[paste0("tau_nested_", tag)]] <- nf$fit$tau
    }
    out <- out %>%
      mutate(ll_both = x$both$ll, short_both = x$both$r_s,
             long_both = x$both$r_l, lower_both = x$both$lower,
             upper_both = x$both$upper, sd_short_both = x$both$fit$sd[1],
             sd_long_both = x$both$fit$sd[2], tau_both = x$both$fit$tau,
             long_share_both = x$both$fit$sd[2] ^ 2 /
               sum(c(x$both$fit$sd ^ 2, x$both$fit$tau ^ 2)),
             ll_short_max = x$short$ll, short_max = x$short$r_s,
             lr_long_vs_short = 2 * (x$both$ll - x$short$ll))
    for (i in seq_len(nrow(x$fixed))) {
      tag <- format(x$fixed$range[i])
      f <- x$fixed[i, ]
      out[[paste0("dll_single_", tag)]] <- f$ll_single - se$ll_max
      out[[paste0("dll_nested_", tag)]] <- f$ll_nested - x$both$ll
      out[[paste0("lr_long_", tag)]] <- 2 * (f$ll_nested - x$short$ll)
      out[[paste0("short_at_", tag)]] <- f$short_nested
      out[[paste0("long_share_", tag)]] <- f$sd_long_nested ^ 2 /
        (f$sd_short_nested ^ 2 + f$sd_long_nested ^ 2 + f$tau_nested ^ 2)
      out[[paste0("sd_long_", tag)]] <- f$sd_long_nested
      out[[paste0("sd_short_", tag)]] <- f$sd_short_nested
      out[[paste0("tau_", tag)]] <- f$tau_nested
    }
    out %>%
      mutate(aic_nugget = aic(ll_nugget, 2), aic_se = aic(ll_se, 4),
             aic_exp = aic(ll_exp, 4),
             aic_short_300 = aic(ll_short_300, 3),
             aic_short_600 = aic(ll_short_600, 3),
             aic_nested_300 = aic(ll_nested_300, 5),
             aic_nested_600 = aic(ll_nested_600, 5),
             aic_both = aic(ll_both, 6),
             lr_long_given_300 = 2 * (ll_nested_300 - ll_short_300),
             lr_long_given_600 = 2 * (ll_nested_600 - ll_short_600),
             lr_two_vs_single = 2 * (ll_both - ll_se),
             lr_single_vs_nugget = 2 * (ll_se - ll_nugget))
  })) %>%
    mutate(quantity = factor(quantity, names(quantities))) %>%
    arrange(base, quantity, block_size)
  write.csv(fits, file.path(output_dir, sprintf("range_fits%s.csv", suffix)),
            row.names = FALSE)

  # the profiles, as deviances from each model's maximum
  profiles <- bind_rows(lapply(results, function(x) {
    head <- tibble(base = x$set$base[1], block_size = x$set$block_size[1],
                   quantity = x$set$quantity[1])
    bind_rows(
      tibble(model = "single, squared exponential", range = range_grid,
             ll = x$single$se$ll, ll_max = x$single$se$ll_max),
      tibble(model = "single, exponential", range = range_grid,
             ll = x$single$exp$ll, ll_max = x$single$exp$ll_max),
      tibble(model = "nested: long range (short profiled)",
             range = long_grid, ll = x$both$profile, ll_max = x$both$ll),
      bind_rows(lapply(seq_along(short_fixed), function(i) {
        tibble(model = sprintf("nested: long range (short %s km)",
                               format(short_fixed[i])),
               range = long_grid, ll = x$nested_fixed[[i]]$ll,
               ll_max = x$nested_fixed[[i]]$ll_max)
      }))) %>%
      mutate(deviance = 2 * (ll_max - ll), base = head$base,
             block_size = head$block_size, quantity = head$quantity,
             .before = 1)
  }))
  write.csv(profiles, file.path(output_dir, sprintf("profiles%s.csv", suffix)),
            row.names = FALSE)

  # the empirical semivariograms, and the fitted ones (the semivariance of
  # the true values, without estimation noise: tau^2 + sum sd_j^2 (1 - K_j)).
  # Each pair's half squared difference has expectation gamma(h) + noise_ij,
  # noise_ij = (SE_i^2 + SE_j^2) / 2, and for normal estimates variance 2
  # (gamma + noise_ij)^2; so the pairs are weighted by 1 / (v + noise_ij)^2,
  # v the true-value variance of the single-range fit (sd^2 + tau^2), an
  # approximate inverse variance that keeps a few blocks with huge SEs from
  # swamping a bin. The weights depend on the SEs, not on the estimates, so
  # the weighted mean less the weighted noise stays unbiased for the
  # (weighted) semivariance
  breaks <- seq(0, max_distance, by = bin_width)
  variograms <- bind_rows(lapply(results, function(x) {
    set <- x$set
    v <- x$single$se$fit$sd ^ 2 + x$single$se$fit$tau ^ 2
    pairs <- which(upper.tri(x$D), arr.ind = TRUE)
    tibble(distance = x$D[pairs],
           half_square = (set$y[pairs[, 1]] - set$y[pairs[, 2]]) ^ 2 / 2,
           noise = (set$se[pairs[, 1]] ^ 2 + set$se[pairs[, 2]] ^ 2) / 2,
           weight = 1 / (v + noise) ^ 2) %>%
      filter(distance < max_distance) %>%
      mutate(bin = cut(distance, breaks, right = FALSE)) %>%
      group_by(bin) %>%
      summarise(pairs = n(), distance = weighted.mean(distance, weight),
                raw = weighted.mean(half_square, weight),
                noise = weighted.mean(noise, weight),
                raw_unweighted = mean(half_square),
                noise_unweighted = mean(noise),
                .groups = "drop") %>%
      mutate(corrected = raw - noise, base = set$base[1],
             block_size = set$block_size[1], quantity = set$quantity[1],
             .before = 1)
  }))
  write.csv(variograms, file.path(output_dir, sprintf("variograms%s.csv",
                                                      suffix)),
            row.names = FALSE)
  fit_models <- c("single SE, range at its maximum",
                  "single exponential, range at its maximum",
                  "nested SE, short + long at their maxima",
                  "single SE, range fixed at 2,500 km",
                  "nested SE, short profiled + long fixed at 2,500 km")
  h <- seq(1, max_distance, length.out = 200)
  curves <- bind_rows(lapply(results, function(x) {
    se <- x$single$se
    ex <- x$single$exp
    both <- x$both
    fixed <- filter(x$fixed, range == 2500)
    bind_rows(
      tibble(model = fit_models[1], distance = h,
             gamma = se$fit$tau ^ 2 + se$fit$sd ^ 2 *
               (1 - kernels$se(h, se$r_hat))),
      tibble(model = fit_models[2], distance = h,
             gamma = ex$fit$tau ^ 2 + ex$fit$sd ^ 2 *
               (1 - kernels$exp(h, ex$r_hat))),
      tibble(model = fit_models[3], distance = h,
             gamma = both$fit$tau ^ 2 +
               both$fit$sd[1] ^ 2 * (1 - kernels$se(h, both$r_s)) +
               both$fit$sd[2] ^ 2 * (1 - kernels$se(h, both$r_l))),
      tibble(model = fit_models[4], distance = h,
             gamma = fixed$tau_single ^ 2 + fixed$sd_single ^ 2 *
               (1 - kernels$se(h, 2500))),
      tibble(model = fit_models[5], distance = h,
             gamma = fixed$tau_nested ^ 2 +
               fixed$sd_short_nested ^ 2 *
               (1 - kernels$se(h, fixed$short_nested)) +
               fixed$sd_long_nested ^ 2 * (1 - kernels$se(h, 2500)))) %>%
      mutate(base = x$set$base[1], block_size = x$set$block_size[1],
             quantity = x$set$quantity[1], .before = 1)
  }))

  # printed summary
  cat("\nsingle-range and nested fits (ranges in km; ll REML; intervals NA at a grid edge, 200 or 6,000 km single, 1,000 or 6,000 km long)\n")
  print(as.data.frame(fits %>%
                        transmute(base, size = block_size, quantity, n = blocks,
                                  sd_y, med_se = median_se,
                                  se = sprintf("%4.0f [%4.0f, %4.0f]", range_se,
                                               lower_se, upper_se),
                                  sd_se, tau_se,
                                  exp = sprintf("%4.0f [%4.0f, %4.0f]",
                                                range_exp, lower_exp,
                                                upper_exp),
                                  long_600 = sprintf("%4.0f [%4.0f, %4.0f]",
                                                     long_nested_600,
                                                     lower_nested_600,
                                                     upper_nested_600),
                                  both = sprintf("%3.0f + %4.0f [%4.0f, %4.0f]",
                                                 short_both, long_both,
                                                 lower_both, upper_both),
                                  long_share = long_share_both,
                                  lr_vs_nug = lr_single_vs_nugget,
                                  lr_long_300 = lr_long_given_300,
                                  lr_long_600 = lr_long_given_600,
                                  lr_two_vs_one = lr_two_vs_single)),
        digits = 3)
  cat("\nAIC less the best of each analysis\n")
  print(as.data.frame(fits %>%
                        select(base, block_size, quantity, starts_with("aic_")) %>%
                        rowwise() %>%
                        mutate(best = min(c_across(starts_with("aic_")))) %>%
                        ungroup() %>%
                        mutate(across(starts_with("aic_"), ~ .x - best)) %>%
                        select(-best)),
        digits = 3)
  cat("\na long component at fixed ranges, the short range profiled over 200-1,000 km: LR against the short component alone, and the long share of the true-value variance (sd_l^2 / (sd_s^2 + sd_l^2 + tau^2))\n")
  print(as.data.frame(fits %>%
                        transmute(base, size = block_size, quantity,
                                  short_alone = short_max,
                                  lr_long_free = lr_long_vs_short,
                                  across(matches("^(lr_long|long_share)_[0-9]+$"))
                                  )),
        digits = 2)
  cat("\nlog likelihood at fixed ranges less the maximum (single SE; nested with the short range profiled)\n")
  print(as.data.frame(fits %>%
                        select(base, block_size, quantity,
                               starts_with("dll_"))),
        digits = 3)
  # the phenomenological floors against the model-conditioned ones and the
  # current level, from the default set's range_fits.csv
  if (pheno) {
    comparison <- read.csv(file.path(output_dir, "range_fits.csv")) %>%
      filter(quantity == "u" | (quantity == "level" & base == "wb_ref")) %>%
      bind_rows(mutate(fits, quantity = as.character(quantity))) %>%
      mutate(quantity = case_when(
        base == "pheno" ~ "floor, phenomenological",
        quantity == "level" ~ "current level (GLMM)",
        TRUE ~ sprintf("floor, model-conditioned (%s)", base))) %>%
      arrange(block_size, quantity)
    cat("\nphenomenological floors against the model-conditioned floors (variant d) and the current level: ranges in km [95% interval]\n")
    print(as.data.frame(comparison %>%
                          transmute(quantity, size = block_size, n = blocks,
                                    med_se = median_se,
                                    se = sprintf("%4.0f [%4.0f, %4.0f]",
                                                 range_se, lower_se, upper_se),
                                    exp = sprintf("%4.0f [%4.0f, %4.0f]",
                                                  range_exp, lower_exp,
                                                  upper_exp),
                                    both = sprintf("%3.0f + %4.0f [%4.0f, %4.0f]",
                                                   short_both, long_both,
                                                   lower_both, upper_both),
                                    long_share = long_share_both,
                                    lr_vs_nug = lr_single_vs_nugget,
                                    lr_two_vs_one = lr_two_vs_single,
                                    dll_1500 = dll_single_1500)),
          digits = 3)
  }


  # figures ---------------------------------------------------------------------------

  base_names <- c(wb_ref = "wb_ref (no floor)", wb_bf = "wb_bf (one global floor)",
                  pheno = "none (phenomenological floors)")
  pheno_formula <- sprintf(paste0(
    "p = f_b + (1 - f_b) plogis(a_b + gamma_type - exp(beta_b) (year - %d));",
    "\ngamma fixed from a joint fit; weighted binomial; no dynamical model."),
    pheno_origin)
  quantity_note <- if (pheno) {
    paste("\nu: logit floor of each block's own decline,", pheno_formula)
  } else {
    paste("u and s: variant d (s, logit floor and an l0 shift per",
          "block);\nlevel: logit current LLIN-pyrethroid level, binomial",
          "GLMM of the 2019-24 bioassays (model-free, the same for both",
          "bases).")
  }
  facet_labels <- labeller(
    quantity = function(x) unname(quantities[x]),
    block_size = function(x) paste0(x, "-degree blocks"))
  # Okabe-Ito, checked with the dataviz palette validator (light surface:
  # all pass), with line types as a second encoding; the fixed-range
  # reference grey and dotted
  model_colours <- setNames(c("#0072B2", "#D55E00", "#009E73", grey(0.55),
                              "#000000"), fit_models)
  model_lines <- setNames(c("solid", "42", "4212", "11", "2212"), fit_models)

  for (label in intersect(labels, unique(data$base))) {
    v <- variograms %>%
      filter(base == label, quantity %in% main_quantities) %>%
      mutate(quantity = factor(quantity, main_quantities))
    if (nrow(v) == 0) next
    v_long <- v %>%
      select(block_size, quantity, distance, pairs, raw, corrected) %>%
      pivot_longer(c(raw, corrected), names_to = "estimate",
                   values_to = "gamma") %>%
      mutate(estimate = factor(estimate, c("corrected", "raw"),
                               c("less estimation noise", "raw")))
    cv <- curves %>%
      filter(base == label, quantity %in% main_quantities) %>%
      mutate(quantity = factor(quantity, main_quantities),
             model = factor(model, fit_models))
    figure <- ggplot(mapping = aes(distance, gamma)) +
      geom_hline(yintercept = 0, colour = grey(0.7), linewidth = 0.3) +
      geom_line(aes(colour = model, linetype = model), data = cv,
                linewidth = 0.6) +
      geom_point(aes(size = pairs, shape = estimate), data = v_long,
                 colour = grey(0.15), fill = grey(0.55), stroke = 0.4) +
      facet_grid(quantity ~ block_size, scales = "free_y",
                 labeller = facet_labels) +
      scale_colour_manual(values = model_colours, name = NULL) +
      scale_linetype_manual(values = model_lines, name = NULL) +
      guides(colour = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2)) +
      scale_shape_manual(values = c(`less estimation noise` = 21, raw = 1),
                         name = "empirical") +
      scale_size_area(max_size = 3, name = "pairs") +
      scale_x_continuous(labels = scales::comma,
                         breaks = seq(0, max_distance, 1000)) +
      labs(x = "great-circle distance between block centroids (km)",
           y = "semivariance",
           title = sprintf("Semivariograms of block estimates, base %s",
                           base_names[[label]]),
           caption = paste(
             "Points: empirical semivariance in 250 km bins, pairs weighted by",
             "1 / (v + noise_ij)^2, raw (open) and less the pairs' estimation",
             "noise (SE_i^2 + SE_j^2) / 2 (filled).\nLines: REML fits of y_b ~ N(mu, sd^2 K(d) + tau^2 I",
             "+ diag(SE_b^2)), drawn without the estimation noise (tau^2 +",
             "sd^2 (1 - K));\nSE and exponential kernels with correlation",
             "0.135 at the range.", quantity_note)) +
      theme_bw(base_size = 9) +
      theme(legend.position = "bottom", legend.box = "vertical",
            plot.caption = element_text(hjust = 0),
            panel.grid.minor = element_blank())
    ggsave(file.path(figure_dir, sprintf("variograms_%s.png", label)), figure,
           width = 9, height = if (pheno) 5.5 else 9, dpi = 150, bg = "white")

    profile_models <- c("single, squared exponential", "single, exponential",
                        "nested: long range (short profiled)")
    profile_colours <- setNames(c("#0072B2", "#D55E00", "#009E73"),
                                profile_models)
    profile_lines <- setNames(c("solid", "42", "4212"), profile_models)
    p <- profiles %>%
      filter(base == label, model %in% profile_models) %>%
      mutate(quantity = factor(quantity, names(quantities)),
             model = factor(model, profile_models),
             deviance = pmin(deviance, 12))
    figure <- ggplot(p, aes(range, deviance, colour = model,
                            linetype = model)) +
      annotate("rect", xmin = 2000, xmax = 3000, ymin = -Inf, ymax = Inf,
               fill = grey(0.92)) +
      geom_hline(yintercept = qchisq(0.95, 1), colour = grey(0.5),
                 linewidth = 0.3) +
      geom_vline(xintercept = 1500, colour = grey(0.5), linewidth = 0.3,
                 linetype = "22") +
      geom_line(linewidth = 0.6) +
      facet_grid(quantity ~ block_size, labeller = facet_labels) +
      scale_x_log10(breaks = c(200, 500, 1000, 2000, 3000, 6000),
                    labels = scales::comma) +
      scale_colour_manual(values = profile_colours, name = NULL) +
      scale_linetype_manual(values = profile_lines, name = NULL) +
      coord_cartesian(ylim = c(0, 12)) +
      labs(x = "range (km; correlation 0.135 at the range)",
           y = "2 x (max log likelihood - profile), REML",
           title = sprintf("Profile likelihoods of the range, base %s",
                           base_names[[label]]),
           caption = paste(
             "Horizontal line: 3.84, the approximate 95% interval. Dashed",
             "vertical: 1,500 km, the 5% point of V5's range prior. Shaded:",
             "2,000-3,000 km.\nNested: short SE component (range profiled",
             "over 200-1,000 km) plus a long one (1,000-6,000 km), the",
             "profile of the long range. Deviances capped at 12.")) +
      theme_bw(base_size = 9) +
      theme(legend.position = "bottom",
            plot.caption = element_text(hjust = 0),
            panel.grid.minor = element_blank())
    ggsave(file.path(figure_dir, sprintf("profiles_%s.png", label)), figure,
           width = 8, height = if (pheno) 4.5 else 11, dpi = 150, bg = "white")
  }

  # maps of the block floors (variant d, or pheno) with their SEs
  borders <- readRDS("data/clean/country_borders.RDS")
  water_mask <- sf::st_as_sf(terra::as.polygons(terra::aggregate(
    rast("data/clean/pfpr_water_mask.tif"), 4, fun = "max", na.rm = TRUE)))
  for (label in intersect(labels, unique(data$base))) {
    map_variant <- if (label == "pheno") "pheno" else "d"
    squares <- bind_rows(lapply(block_sizes, function(size) {
      file <- file.path(output_dir, sprintf("estimates_%s_%s.csv", label,
                                            format(size)))
      read.csv(file) %>%
        filter(eligible, variant == map_variant) %>%
        mutate(xmin = block_x * size, xmax = xmin + size,
               ymin = block_y * size, ymax = ymin + size)
    })) %>%
      mutate(size_label = paste0(block_size, "-degree blocks"))
    extent <- c(range(c(squares$xmin, squares$xmax)) + c(-2, 2),
                range(c(squares$ymin, squares$ymax)) + c(-2, 2))
    block_map <- function(fill, scale, title) {
      ggplot() +
        geom_sf(data = borders, fill = grey(0.95), colour = NA) +
        geom_sf(data = water_mask, fill = grey(0.85), colour = NA) +
        geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                      fill = .data[[fill]]), data = squares,
                  colour = "white", linewidth = 0.3) +
        scale +
        geom_sf(data = borders, fill = NA, colour = grey(0.4),
                linewidth = 0.1) +
        facet_wrap(~ size_label, ncol = 1) +
        coord_sf(xlim = extent[1:2], ylim = extent[3:4], expand = FALSE) +
        labs(title = title) +
        guides(fill = guide_colourbar(barwidth = 12)) +
        theme_ir_maps() +
        theme(legend.position = "bottom",
              strip.text = element_text(size = 9))
    }
    # the phenomenological floors reach 0.98
    floor_breaks <- c(0.05, 0.1, 0.2, 0.4, 0.7, if (label == "pheno") 0.95)
    floor_limits <- c(0.03, if (label == "pheno") 0.98 else 0.85)
    blocks_note <- sprintf(paste(
      "Pyrethroid bioassays in square blocks (>= %d bioassays in >= %d",
      "years, >= %d from %d)."), min_bioassays, min_years, min_late,
      late_from)
    map_caption <- if (label == "pheno") {
      paste0(blocks_note, " Each block's own decline to a floor,\n",
             pheno_formula, " SE: the larger of the Hessian and the",
             " site-clustered sandwich;\nfloor prior N(0, 2.5^2), so an SE",
             " above 1 means the data barely identify the floor.")
    } else {
      sprintf(paste(
        blocks_note, "Variant d: selection multiplier, floor and",
        "an initial-state shift per block,\nwith the weighted binomial, the",
        "base %s held at its posterior mean. SE: the larger of the Hessian",
        "and the site-clustered sandwich;\nfloor prior N(-1.83, 1.39^2), so",
        "an SE near 1.39 means the data barely identify the floor."),
        base_names[[label]])
    }
    maps <- block_map("u", scale_fill_distiller(
      palette = "Blues", direction = 1, name = "floor f_b",
      breaks = qlogis(floor_breaks), labels = floor_breaks,
      limits = qlogis(floor_limits), oob = scales::squish),
      "block floor f_b (logit scale)") +
      block_map("se_u", scale_fill_distiller(
        palette = "Oranges", direction = 1, name = "SE of logit f_b",
        limits = c(0, if (label == "pheno") 2.5 else 1.4),
        oob = scales::squish),
        "its standard error (logit scale)") +
      plot_annotation(caption = map_caption,
        theme = theme(plot.caption = element_text(hjust = 0, size = 8)))
    ggsave(file.path(figure_dir, sprintf("floor_map_%s.png", label)), maps,
           width = 10, height = 9.5, dpi = 150, bg = "white")
  }
  report("written %s and figures in %s; peak memory %.1f GB",
         file.path(output_dir, sprintf("range_fits%s.csv", suffix)),
         figure_dir, peak_memory_gb())
}
