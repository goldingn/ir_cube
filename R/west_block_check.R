# A likelihood check of what stops the models fitting the fall of pyrethroid
# mortality in Côte d'Ivoire, Burkina Faso and Ghana (#47). There the
# LLIN-pyrethroid bioassays fall steadily, from about 0.55 in 2010-2015 to
# 0.12-0.20 from 2019, and V5 meets them with selection 1.65 times stronger
# and a floor of 0.24, which its pyrethroid curve reaches by about 2015: too
# resistant in 2010-2015, flat after. Is that because V5's selection
# multiplier is shared by every insecticide class, so that DDT and
# bendiocarb, which lose susceptibility there faster than predicted, pull it
# up? Can a multiplier and floor for the pyrethroids alone fit the subregion
# early and late together?
#
#   Rscript R/west_block_check.R [fit file]
#
# The base, as in R/floor_selection_blocks.R, is the reference fit ref_f0 (no
# floor, kdr or species) at its posterior mean over the usable chains, in its
# outer form at every bioassay of the West, now of every class:
#   logit q = l0 + d - exp(s) C_t - exp(r) t kappa,  p = f + (1 - f) q,
# the logit initial state l0, cumulative log fitness C_t and reversion
# t kappa checked against the fit's own predictions (s = r = d = 0, f = 0).
# The blocks are the subregions of the West of R/west_figures.R
# (west_subregion()):
#   A  Senegal to Liberia, with Mauritania and Mali
#   B  Côte d'Ivoire, Burkina Faso and Ghana
#   C  Togo, Benin, Nigeria and Niger
# In each, the beta-binomial log posterior, with the priors of
# R/floor_selection_blocks.R (s ~ N(0, 1), u = logit f ~ N(-1.83, 1.39^2)), is
# maximised in four variants:
#   a  pyrethroid only    (s, u) from the pyrethroid bioassays, as the screen;
#                         the other classes stay at the base
#   b  reversion free     as a, with r ~ N(0, 1), a multiplier exp(r) on the
#                         base fit's reversion of the pyrethroids
#   c  shared multiplier  (s, u) from the bioassays of every class, one s for
#                         all, as V5's u_s; the floor on the pyrethroids only,
#                         none on the others, as in ref_f0
#   d  initial state      as a, with d ~ N(0, 1), a shift of the pyrethroids'
#                         logit initial state l0 in the block: whether the
#                         base's initial states, held fixed here (V5 refits
#                         them, per country and type), force the shape
# with standard errors from the Hessian; and for each class alone, s without
# a floor, the selection each would have. The fit of the bioassays under
# each, with the base and V5 (its posterior mean predictions,
# outputs/species_runs/misfit/V5_bioassays.csv, R/species_misfit.R) for
# comparison, is summarised before 2010, in 2010-2015, 2016-2018 and 2019
# onwards: observed and predicted pooled mortality (died / tested, the
# predictions weighted by mosquitoes tested, so site-matched), the mean logit
# misfit (empirical logit, log((died + 0.5) / (survived + 0.5)), minus the
# predicted, as R/species_misfit.R: positive where more died than predicted,
# the model too resistant), and the log likelihood (at the base's rho; none
# for V5).
#
# Writes, in outputs/species_runs/west/:
#   block_check_estimates.csv  per block and variant: the estimates and their
#                              SEs, and log likelihoods and posteriors of the
#                              pyrethroid bioassays, the others and all
#   block_check_classes.csv    per block and class: s of the class alone
#   block_check_selection_by_window.csv
#                              per block and window (with pre-2010): s of
#                              the window's pyrethroid bioassays alone, with
#                              no floor and with a's floor, and their mean
#                              l0, C_t and t kappa
#   block_check_windows.csv    per block, variant, set of bioassays (the LLIN
#                              pyrethroids alpha-cypermethrin, deltamethrin
#                              and permethrin; all pyrethroids; each other
#                              class) and window: the fit
#   block_check_constant.csv   per block and window, the pyrethroid
#                              bioassays: the spread of their mortality, the
#                              best constant mortality under the base's rho
#                              (the level the likelihood wants, one free per
#                              window) and its log likelihood, against each
#                              variant's
#   block_check_floor_profile.csv
#                              per block and fixed floor (0.02-0.3), s (a),
#                              or s and d (d), refitted to the pyrethroid
#                              bioassays: log likelihood and LLIN-pyrethroid
#                              predicted pooled mortality by window
# and figures/species_runs/west/block_check.png: per block, observed
# LLIN-pyrethroid mortality by year with each variant's site-matched
# prediction, and in B the other classes under the base and variant c. Plain
# R; about 2 GB, a minute or two.
#
# Caveats: the base fit is held fixed, so its country initial states (fitted
# with the base's selection and no floor) absorb level differences between
# and within the blocks before s, f, r and d are fitted, and its covariate
# effects set the timing of the cumulative log fitness C_t; the variants
# change only selection, the floor, reversion and the initial state, each
# uniformly within a block.

arguments <- commandArgs(trailingOnly = TRUE)
fit_file <- if (length(arguments) >= 1) arguments[1] else
  paste0("../ir_cube_netscreen/outputs/pod_jobs/dh270_lin_f0_full/",
         "temporary/fitted_model.RData")
# the priors of R/floor_selection_blocks.R, and the reversion multiplier's
# and the initial-state shift's
floor_prior <- list(mean = -1.8333, sd = 1.3888)
selection_prior_sd <- 1
reversion_prior_sd <- 1
shift_prior_sd <- 1
llin_pyrethroids <- c("Alpha-cypermethrin", "Deltamethrin", "Permethrin")
windows <- c("2010-15", "2016-18", "2019+")
window_of <- function(year) {
  dplyr::case_when(year %in% 2010:2015 ~ windows[1],
                   year %in% 2016:2018 ~ windows[2],
                   year >= 2019 ~ windows[3])
}

suppressMessages({
  library(greta)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
source("R/functions.R")
source("R/dynamical_predictions.R")
source("R/species_fit_helpers.R")

output_dir <- "outputs/species_runs/west"
figure_dir <- "figures/species_runs/west"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(all.equal(logit_beta_moments(c(1, 4))$mean, floor_prior$mean,
                    tolerance = 1e-4),
          all.equal(logit_beta_moments(c(1, 4))$sd, floor_prior$sd,
                    tolerance = 1e-4))

# the subregions of west_subregion() (R/west_figures.R, which runs on
# sourcing), as blocks A, B and C
west_block <- function(country) {
  case_when(
    country %in% c("Senegal", "Gambia", "Guinea-Bissau", "Guinea",
                   "Sierra Leone", "Liberia", "Mauritania", "Mali") ~ "A",
    country %in% c("Côte d’Ivoire", "Burkina Faso", "Ghana") ~ "B",
    country %in% c("Togo", "Benin", "Nigeria", "Niger") ~ "C")
}
block_names <- c(A = "A: Senegal to Liberia, and Mali",
                 B = "B: Côte d'Ivoire, Burkina Faso, Ghana",
                 C = "C: Togo, Benin, Nigeria, Niger")


# the base fit, at its posterior mean (as R/floor_selection_blocks.R) ----------------

fit <- load_fit(fit_file)
stopifnot(!smooth_on(fit$options), isFALSE(fit$options$mortality_floor))
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
stopifnot(length(kappa) == length(fit$types))

west <- df %>%
  mutate(row = row_number()) %>%
  filter(analysis_region(country_name, region) == "West") %>%
  mutate(block = west_block(country_name),
         pyrethroid = insecticide_class == "Pyrethroids",
         llin = insecticide_type %in% llin_pyrethroids,
         window = window_of(year_start))
stopifnot(!anyNA(west$block))

# l0 and C at every bioassay, of every class, as dynamical_logit_cells()
# forms them (the loop of R/floor_selection_blocks.R over all types)
cell_country <- dynamical_lookups(df)$cell_country_lookup
n_times <- max(fit$cell_years_index$year_id)
x_row <- matrix(NA_integer_, max(fit$cell_years_index$cell_id), n_times)
x_row[cbind(fit$cell_years_index$cell_id, fit$cell_years_index$year_id)] <-
  seq_len(nrow(fit$cell_years_index))
west$l0 <- NA_real_
west$C <- NA_real_
for (k in sort(unique(west$type_id))) {
  rows_k <- which(west$type_id == k)
  cells <- sort(unique(west$cell_id[rows_k]))
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
  at <- match(west$cell_id[rows_k], cells)
  west$l0[rows_k] <- l0[1, at]
  west$C[rows_k] <- cumulative[cbind(at, west$year_id[rows_k])]
}
west <- west %>%
  mutate(reversion = year_id * kappa[type_id],
         rho = rho[type_id],
         logit_base = l0 - C - reversion)
# the outer form reproduces the fit's own predictions at these parameters,
# class by class
check <- c(dynamical_logit(parameters, west, df, fit$x_cell_years,
                           fit$cell_years_index))
west$check_difference <- abs(pmin(pmax(check, -30), 30) -
                               pmin(pmax(west$logit_base, -30), 30))
checks <- west %>%
  group_by(insecticide_class) %>%
  summarise(bioassays = n(), max_difference = max(check_difference),
            .groups = "drop")
for (i in seq_len(nrow(checks))) {
  report("outer form at %d %s bioassays of the West; against the fit's predictions, max |logit diff| %.2g",
         checks$bioassays[i], checks$insecticide_class[i],
         checks$max_difference[i])
}
stopifnot(max(checks$max_difference) < 1e-8)
report("reversion kappa (logit per year) of the West's types: %s",
       toString(sprintf("%s %.2g", fit$types[sort(unique(west$type_id))],
                        kappa[sort(unique(west$type_id))])))
report("t kappa at the West's bioassays from 2019: %s",
       paste(capture.output(print(summary(
         west$reversion[west$year_start >= 2019]))), collapse = " "))


# the variants -----------------------------------------------------------------------

clamp <- function(p) pmin(pmax(p, 1e-12), 1 - 1e-12)

# predicted mortality at bioassays `data`: selection multiplier exp(s), the
# floor plogis(u) and a shift d of the initial state on the pyrethroids only,
# and a reversion multiplier exp(r); the base fit at s = 0, u = -Inf, r = 0,
# d = 0. With r = d = 0, and data of pyrethroids alone, the p of
# block_log_lik() in R/floor_selection_blocks.R
outer_p <- function(data, s = 0, u = -Inf, r = 0, d = 0) {
  f <- ifelse(data$pyrethroid, plogis(u), 0)
  shift <- ifelse(data$pyrethroid, d, 0)
  clamp(f + (1 - f) * plogis(data$l0 + shift - exp(s) * data$C -
                               exp(r) * data$reversion))
}
# the beta-binomial log likelihood of each bioassay at mortality p, at the
# base's rho of its type
bioassay_log_lik <- function(data, p) {
  extraDistr::dbbinom(data$died, data$mosquito_number,
                      alpha = p * (1 / data$rho - 1),
                      beta = (1 - p) * (1 / data$rho - 1), log = TRUE)
}
log_lik <- function(data, s = 0, u = -Inf, r = 0, d = 0) {
  sum(bioassay_log_lik(data, outer_p(data, s, u, r, d)))
}
# the priors of the parameters given (s always; u, r and d if fitted)
log_prior <- function(s, u = NULL, r = NULL, d = NULL) {
  dnorm(s, 0, selection_prior_sd, log = TRUE) +
    (if (is.null(u)) 0 else dnorm(u, floor_prior$mean, floor_prior$sd,
                                  log = TRUE)) +
    (if (is.null(r)) 0 else dnorm(r, 0, reversion_prior_sd, log = TRUE)) +
    (if (is.null(d)) 0 else dnorm(d, 0, shift_prior_sd, log = TRUE))
}
log_post <- function(data, theta) {
  theta <- as.list(theta)
  do.call(log_lik, c(list(data), theta)) + do.call(log_prior, theta)
}

# the maximum from a grid of starts (the screen's, and for r and d -1, 0,
# 1), and
# the covariance from the Hessian; `names`, the parameters fitted
starts_grid <- list(s = c(-1.5, 0, 1), u = c(-5, -1.8, 0), r = c(-1, 0, 1),
                    d = c(-1, 0, 1))
fit_variant <- function(data, names) {
  starts <- expand.grid(starts_grid[names])
  objective <- function(theta) -log_post(data, setNames(theta, names))
  runs <- lapply(seq_len(nrow(starts)), function(i) {
    optim(unlist(starts[i, , drop = FALSE]), objective, method = "BFGS",
          control = list(reltol = 1e-12, maxit = 500))
  })
  best <- runs[[which.min(vapply(runs, `[[`, numeric(1), "value"))]]
  covariance <- solve(optimHess(best$par, objective))
  list(estimate = setNames(best$par, names),
       se = setNames(sqrt(diag(covariance)), names),
       log_post = -best$value, converged = best$convergence == 0,
       spread = diff(range(vapply(runs, `[[`, numeric(1), "value"))))
}

# s alone, at a fixed logit floor u (the pyrethroids' only), by its 1-D
# maximum, with its SE from the Hessian
fit_s <- function(data, u = -Inf) {
  objective <- function(s) -(log_lik(data, s, u) + log_prior(s))
  best <- optimize(objective, c(-6, 6), tol = 1e-10)
  c(s = best$minimum,
    se_s = sqrt(1 / c(optimHess(best$minimum, objective))))
}

variants <- c(base = "base (ref_f0)", a = "a: pyrethroid only",
              b = "b: reversion free", c = "c: shared multiplier",
              d = "d: a + initial-state shift")
estimates <- list()
classes <- list()
by_window_s <- list()
for (v in names(variants)) west[[paste0("predicted_", v)]] <- NA_real_
for (b in names(block_names)) {
  data <- filter(west, block == b)
  pyrethroids <- filter(data, pyrethroid)
  others <- filter(data, !pyrethroid)
  fits <- list(a = fit_variant(pyrethroids, c("s", "u")),
               b = fit_variant(pyrethroids, c("s", "u", "r")),
               c = fit_variant(data, c("s", "u")),
               d = fit_variant(pyrethroids, c("s", "u", "d")))
  # each variant's parameters, and the others' (the base in a and b, which
  # change the pyrethroids only; s in c)
  values <- list(base = c(s = 0, u = -Inf, r = 0, d = 0))
  for (v in names(fits)) {
    values[[v]] <- c(s = 0, u = -Inf, r = 0, d = 0)
    values[[v]][names(fits[[v]]$estimate)] <- fits[[v]]$estimate
  }
  for (v in names(variants)) {
    theta <- values[[v]]
    other_s <- if (v == "c") theta[["s"]] else 0
    pyr_lik <- log_lik(pyrethroids, theta[["s"]], theta[["u"]], theta[["r"]],
                       theta[["d"]])
    other_lik <- log_lik(others, other_s)
    se <- if (v == "base") c(s = NA, u = NA, r = NA, d = NA) else
      fits[[v]]$se[c("s", "u", "r", "d")]
    fitted <- if (v == "base") character(0) else names(fits[[v]]$estimate)
    prior <- if (v == "base") 0 else
      do.call(log_prior, as.list(theta[fitted]))
    estimates[[length(estimates) + 1]] <- tibble(
      block = b, variant = v,
      s = theta[["s"]], se_s = unname(se["s"]), multiplier = exp(theta[["s"]]),
      u = theta[["u"]], se_u = unname(se["u"]), f = plogis(theta[["u"]]),
      f_lower = plogis(theta[["u"]] - 1.96 * unname(se["u"])),
      f_upper = plogis(theta[["u"]] + 1.96 * unname(se["u"])),
      r = theta[["r"]], se_r = unname(se["r"]),
      reversion_multiplier = exp(theta[["r"]]),
      d = theta[["d"]], se_d = unname(se["d"]),
      pyrethroid_bioassays = nrow(pyrethroids),
      other_bioassays = nrow(others),
      log_lik_pyrethroids = pyr_lik,
      log_post_pyrethroids = pyr_lik + prior,
      log_lik_others = other_lik,
      log_lik_all = pyr_lik + other_lik,
      log_post_all = pyr_lik + other_lik + prior,
      converged = if (v == "base") NA else fits[[v]]$converged,
      start_spread = if (v == "base") NA else fits[[v]]$spread)
    west[west$block == b, paste0("predicted_", v)] <-
      if (v == "c") {
        outer_p(data, theta[["s"]], theta[["u"]])
      } else {
        ifelse(data$pyrethroid,
               outer_p(data, theta[["s"]], theta[["u"]], theta[["r"]],
                       theta[["d"]]),
               outer_p(data))
      }
  }
  # the selection each class would have alone, without a floor
  for (class in sort(unique(data$insecticide_class))) {
    rows <- filter(data, insecticide_class == class)
    alone <- fit_variant(rows, "s")
    classes[[length(classes) + 1]] <- tibble(
      block = b, class = class, bioassays = nrow(rows),
      types = toString(sort(unique(rows$insecticide_type))),
      s = alone$estimate[["s"]], se_s = alone$se[["s"]],
      multiplier = exp(alone$estimate[["s"]]),
      log_lik_gain = log_lik(rows, alone$estimate[["s"]]) - log_lik(rows))
  }
  # the selection each window's pyrethroid bioassays would have alone, with
  # no floor and with a's: if it rises from window to window, the base's
  # cumulative log fitness C_t grows too slowly late for any one multiplier;
  # with the mean l0, C_t and t kappa there
  for (w in c("pre-2010", windows)) {
    rows <- filter(pyrethroids, coalesce(window, "pre-2010") == w)
    no_floor <- fit_s(rows)
    a_floor <- fit_s(rows, values$a[["u"]])
    by_window_s[[length(by_window_s) + 1]] <- tibble(
      block = b, window = w, bioassays = nrow(rows),
      observed = sum(rows$died) / sum(rows$mosquito_number),
      mean_l0 = mean(rows$l0), mean_C = mean(rows$C),
      mean_reversion = mean(rows$reversion),
      s_no_floor = no_floor[["s"]], se_no_floor = no_floor[["se_s"]],
      multiplier_no_floor = exp(no_floor[["s"]]),
      s_a_floor = a_floor[["s"]], se_a_floor = a_floor[["se_s"]],
      multiplier_a_floor = exp(a_floor[["s"]]))
  }
}
estimates <- bind_rows(estimates)
classes <- bind_rows(classes)
by_window_s <- bind_rows(by_window_s)
stopifnot(all(estimates$converged, na.rm = TRUE))
write.csv(estimates, file.path(output_dir, "block_check_estimates.csv"),
          row.names = FALSE)
write.csv(classes, file.path(output_dir, "block_check_classes.csv"),
          row.names = FALSE)
write.csv(by_window_s, file.path(output_dir, "block_check_selection_by_window.csv"),
          row.names = FALSE)


# V5, for comparison ---------------------------------------------------------------------

v5 <- read.csv("outputs/species_runs/misfit/V5_bioassays.csv")
stopifnot(nrow(v5) == nrow(df), all(v5$cell == df$cell),
          all(v5$year_start == df$year_start),
          all(v5$insecticide_type == df$insecticide_type),
          all(v5$died == df$died))
west$predicted_V5 <- v5$predicted[west$row]
variants <- c(variants, V5 = "V5 (posterior mean)")


# fit by window --------------------------------------------------------------------------

# observed and predicted pooled mortality, the mean logit misfit and the
# beta-binomial log likelihood (at the base's rho; none for V5, whose rho
# differs) of each variant, in long form
summarise_fit <- function(data) {
  data %>%
    mutate(empirical = log((died + 0.5) / (mosquito_number - died + 0.5))) %>%
    pivot_longer(starts_with("predicted_"), names_to = "variant",
                 names_prefix = "predicted_", values_to = "predicted") %>%
    mutate(log_lik = bioassay_log_lik(
      list(died = died, mosquito_number = mosquito_number, rho = rho),
      clamp(predicted))) %>%
    group_by(block, set, window, variant) %>%
    summarise(bioassays = n(), tested = sum(mosquito_number),
              observed = sum(died) / sum(mosquito_number),
              predicted_pooled = sum(predicted * mosquito_number) /
                sum(mosquito_number),
              misfit = mean(empirical - qlogis(predicted)),
              log_lik = if (first(variant) == "V5") NA else sum(log_lik),
              .groups = "drop")
}
sets <- bind_rows(
  west %>% filter(llin) %>% mutate(set = "LLIN pyrethroids"),
  west %>% filter(pyrethroid) %>% mutate(set = "all pyrethroids"),
  west %>% filter(!pyrethroid) %>%
    mutate(set = ifelse(insecticide_type %in% c("DDT", "Bendiocarb"),
                        insecticide_type, insecticide_class)))
by_window <- sets %>%
  mutate(window = coalesce(window, "pre-2010")) %>%
  summarise_fit() %>%
  mutate(variant = factor(variant, names(variants)),
         set = factor(set, unique(sets$set))) %>%
  arrange(block, set, variant, window)
write.csv(by_window, file.path(output_dir, "block_check_windows.csv"),
          row.names = FALSE)


# what level the likelihood wants -----------------------------------------------------------

# per block and window, the spread of the pyrethroid bioassays' mortality and
# the best constant mortality under the base's rho (the beta-binomial
# maximum), against each variant's log likelihood there: a reference with
# one free level per window
constant <- west %>%
  filter(pyrethroid) %>%
  mutate(window = coalesce(window, "pre-2010"),
         proportion = died / mosquito_number) %>%
  group_by(block, window) %>%
  group_modify(function(data, key) {
    best <- optimize(function(l) -sum(bioassay_log_lik(data, plogis(l))),
                     c(-8, 8), tol = 1e-10)
    tibble(bioassays = nrow(data),
           pooled = sum(data$died) / sum(data$mosquito_number),
           mean_proportion = mean(data$proportion),
           median_proportion = median(data$proportion),
           share_at_most_5 = mean(data$proportion <= 0.05),
           share_above_30 = mean(data$proportion > 0.3),
           constant_p = plogis(best$minimum),
           log_lik_constant = -best$objective)
  }) %>%
  ungroup() %>%
  left_join(by_window %>%
              filter(set == "all pyrethroids", variant != "V5") %>%
              select(block, window, variant, log_lik) %>%
              pivot_wider(names_from = variant, values_from = log_lik,
                          names_prefix = "log_lik_"),
            by = c("block", "window"))
write.csv(constant, file.path(output_dir, "block_check_constant.csv"),
          row.names = FALSE)

# per block, the pyrethroid fit at fixed floors: s (a), or s and d (d),
# refitted; the log likelihood and the LLIN pyrethroids' predicted pooled
# mortality by window: which windows hold the floor where it is
floor_grid <- c(0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
floor_profile <- bind_rows(lapply(names(block_names), function(b) {
  pyrethroids <- west %>%
    filter(block == b, pyrethroid) %>%
    mutate(window = coalesce(window, "pre-2010"))
  bind_rows(lapply(floor_grid, function(f) {
    bind_rows(lapply(c("a", "d"), function(v) {
      names <- if (v == "a") "s" else c("s", "d")
      objective <- function(theta) {
        -log_post(pyrethroids, c(setNames(theta, names), u = qlogis(f)))
      }
      best <- optim(rep(0, length(names)), objective, method = "BFGS",
                    control = list(reltol = 1e-12, maxit = 500))
      theta <- c(s = 0, d = 0)
      theta[names] <- best$par
      pyrethroids %>%
        mutate(p = outer_p(pyrethroids, theta[["s"]], qlogis(f), 0,
                           theta[["d"]]),
               log_lik = bioassay_log_lik(pyrethroids, p)) %>%
        group_by(window) %>%
        summarise(log_lik = sum(log_lik),
                  llin_predicted = sum((p * mosquito_number)[llin]) /
                    sum(mosquito_number[llin]),
                  .groups = "drop") %>%
        mutate(block = b, variant = v, f = f, s = theta[["s"]],
               d = theta[["d"]], log_lik_total = sum(log_lik), .before = 1)
    }))
  }))
}))
write.csv(floor_profile, file.path(output_dir, "block_check_floor_profile.csv"),
          row.names = FALSE)


# printed summary --------------------------------------------------------------------------

options(width = 200)
cat("\nestimates (s: log selection multiplier against the base; f = plogis(u), 95% from the Hessian; r: log reversion multiplier)\n")
print(as.data.frame(estimates %>%
                      transmute(block, variant, s, se_s, multiplier, u, se_u,
                                f, f_lower, f_upper, r, se_r, d, se_d,
                                lp_pyr = log_post_pyrethroids,
                                ll_pyr = log_lik_pyrethroids,
                                ll_others = log_lik_others,
                                lp_all = log_post_all, ll_all = log_lik_all,
                                start_spread)),
      digits = 4)
cat("\nselection of each class alone, no floor\n")
print(as.data.frame(classes), digits = 3)
cat("\nselection of each window's pyrethroid bioassays alone, no floor and at a's floor\n")
print(as.data.frame(by_window_s), digits = 3)
cat("\nthe pyrethroid bioassays by window: spread, best constant p under the base's rho, and log likelihoods\n")
print(as.data.frame(constant), digits = 4)
cat("\nthe pyrethroid fit at fixed floors: log likelihood (LLIN predicted pooled) by window\n")
print(as.data.frame(floor_profile %>%
                      mutate(cell = sprintf("%6.0f (%.2f)", log_lik,
                                            llin_predicted)) %>%
                      select(block, variant, f, s, d, log_lik_total, window,
                             cell) %>%
                      pivot_wider(names_from = window, values_from = cell)),
      digits = 4, right = FALSE)
cat("\nfit by window: observed / predicted pooled mortality, mean misfit, log likelihood\n")
wide <- by_window %>%
  mutate(cell = sprintf("%.2f/%.2f %+.2f %6.0f", observed, predicted_pooled,
                        misfit, log_lik)) %>%
  select(block, set, variant, window, cell) %>%
  pivot_wider(names_from = window, values_from = cell)
counts <- by_window %>%
  filter(variant == "base") %>%
  select(block, set, window, bioassays) %>%
  pivot_wider(names_from = window, values_from = bioassays,
              names_prefix = "n ")
print(as.data.frame(left_join(wide, counts, by = c("block", "set")) %>%
                      filter(!(set %in% c("LLIN pyrethroids",
                                          "all pyrethroids")) &
                               variant %in% c("base", "c") |
                               set %in% c("LLIN pyrethroids",
                                          "all pyrethroids"))),
      right = FALSE)


# figure ---------------------------------------------------------------------------------

# Okabe-Ito; the base grey and dashed, b dotted, as b and V5 are close in
# deuteranopia, and d dot-dashed
variant_colours <- c(base = grey(0.45), a = "#0072B2", b = "#CC79A7",
                     c = "#D55E00", d = "#56B4E9", V5 = "#009E73")
variant_lines <- c(base = "42", a = "solid", b = "22", c = "solid",
                   d = "4212", V5 = "solid")
yearly <- sets %>%
  filter(year_start >= 2008) %>%
  mutate(window = as.character(year_start)) %>%
  summarise_fit() %>%
  mutate(year = as.integer(window),
         variant = factor(variant, names(variants)))
estimate_note <- function(b) {
  e <- filter(estimates, block == b, variant %in% c("a", "b", "c", "d"))
  paste(sprintf("%s: x%.2f, f %.2f%s", e$variant, e$multiplier, e$f,
                case_when(
                  e$variant == "b" ~ sprintf(", reversion x%.2f",
                                             e$reversion_multiplier),
                  e$variant == "d" ~ sprintf(", l0 %+.2f", e$d),
                  TRUE ~ "")),
        collapse = "\n")
}
trend_panel <- function(b, s, keep = names(variants), title = NULL,
                        note = NULL) {
  points <- filter(yearly, block == b, set == s, variant == "base")
  lines <- filter(yearly, block == b, set == s, variant %in% keep)
  ggplot(mapping = aes(year)) +
    annotate("rect", xmin = c(2009.5, 2018.5), xmax = c(2015.5, 2024.5),
             ymin = -Inf, ymax = Inf, fill = grey(0.94)) +
    geom_line(aes(y = predicted_pooled, colour = variant,
                  linetype = variant), data = lines, linewidth = 0.7) +
    geom_point(aes(y = observed, size = tested), data = points, shape = 21,
               fill = grey(0.85), colour = grey(0.1), stroke = 0.4) +
    {if (!is.null(note)) annotate("text", x = 2008, y = 0.02, hjust = 0,
                                  vjust = 0, size = 2.6, label = note,
                                  lineheight = 0.95)} +
    scale_colour_manual(values = variant_colours, labels = variants,
                        drop = FALSE, name = "prediction") +
    scale_linetype_manual(values = variant_lines, labels = variants,
                          drop = FALSE, name = "prediction") +
    scale_size_area(max_size = 5, limits = c(0, 40000),
                    breaks = c(500, 5000, 20000), name = "mosquitoes tested") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = NULL, y = "mortality", title = if (is.null(title)) s else title) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10))
}
top <- lapply(names(block_names), function(b) {
  trend_panel(b, "LLIN pyrethroids",
              title = sprintf("%s\nLLIN pyrethroids", block_names[[b]]),
              note = estimate_note(b))
})
other_sets <- c("DDT", "Bendiocarb", "Organophosphates")
bottom <- lapply(other_sets, function(s) {
  e <- filter(estimates, block == "B", variant == "c")
  class_s <- unique(west$insecticide_class[west$insecticide_type == s |
                                             west$insecticide_class == s])
  alone <- filter(classes, block == "B", class == class_s)
  trend_panel("B", s, keep = c("base", "c"),
              title = sprintf("B: %s", s),
              note = sprintf("c: x%.2f\nclass alone: x%.2f (s %.2f, SE %.2f)",
                             e$multiplier, alone$multiplier, alone$s,
                             alone$se_s)) +
    guides(colour = "none", linetype = "none")
})
figure <- wrap_plots(c(top, bottom), ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "West Africa subregions: what a block-level selection multiplier, floor, reversion and initial state do to the fit of the pyrethroid decline",
    caption = paste(
      "Points: observed mortality pooled by year (died / tested), sized by",
      "mosquitoes tested. Lines: each variant's prediction at the same",
      "bioassays, weighted by mosquitoes tested. Shaded: 2010-2015 and",
      "2019-2024.\nBase: ref_f0 at its posterior mean, held fixed (its",
      "country initial states absorb level differences). a and b: a",
      "selection multiplier exp(s) and floor f for the pyrethroids alone (b",
      "also a multiplier on the base reversion);\nc: one multiplier for",
      "every class, the floor on the pyrethroids only; d: as a, with a shift",
      "of the pyrethroids' logit initial state l0. Priors s ~ N(0, 1),",
      "logit f ~ N(-1.83, 1.39^2), r ~ N(0, 1), d ~ N(0, 1).\nV5: its",
      "posterior mean prediction (chains 1 and 2). Bottom row: block B's",
      "other classes under the base (also a, b and d) and c, the same line",
      "styles."),
    theme = theme(plot.caption = element_text(hjust = 0, size = 8.5),
                  plot.title = element_text(size = 12)))
ggsave(file.path(figure_dir, "block_check.png"), figure, width = 14,
       height = 9, dpi = 130, bg = "white")
report("written %s and %s; peak memory %.1f GB",
       file.path(output_dir, "block_check_*.csv"),
       file.path(figure_dir, "block_check.png"), peak_memory_gb())
