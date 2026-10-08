# Why V5 does not follow the continued fall of LLIN-pyrethroid mortality in
# diagnostic region A (Burkina Faso and Côte d'Ivoire, with Liberia attached;
# outputs/species_runs/regions/regions.csv, R/diagnostic_regions.R; #47).
# There the bioassays fall from about 85% (2006) to about 30% (2014-2015),
# then more slowly to 11-15% (2021-2024), but both V5 fits, wb_v5 (weighted
# binomial likelihood, all four chains; not converged) and V5 (beta-binomial,
# chains 1 and 2), level off at 23-24% from 2017, although their latent
# smooths (R/latent_smooth.R) could move A's selection multiplier exp(u_s)
# and floor f = plogis(floor_intercept + u_f) on their own.
#
#   USE_CHAINS="V5=1,2" Rscript R/v5_region_a.R
#
# Each fit is taken apart, draw by draw (n_draws, the same number from each
# usable chain), at every bioassay of regions A, B (the Sahel band) and C
# (Benin, Ghana, Togo), in its outer form (outer_logit(), R/dynamical_model.R):
#   logit q = l0 - exp(u_s) C_t - t kappa,  p = f + (1 - f) q,
#   f = plogis(floor_intercept[class] + u_f)
# with l0 the logit initial state of the cell's country and type, C_t the
# cumulative log fitness of the covariates to the bioassay's year, t kappa
# the reversion (kappa <= 0, so - t kappa raises logit q), and u_f 0 for the
# classes without the floor smooth (all but the pyrethroids and DDT); checked
# against the fit's own predictions (dynamical_logit()). Then:
#   1. the smooths at the LLIN-pyrethroid bioassays (alpha-cypermethrin,
#      deltamethrin, permethrin) of A (and of its countries), B and C: per
#      draw, means weighted by mosquitoes tested of u_s, u_f, exp(u_s) and
#      the pyrethroid floor f; their posterior mean, sd and 95% interval,
#      over all usable chains and per chain; and the spread between the
#      bioassays of the posterior mean
#   2. the decomposition by year in A, 2005-2024, at A's LLIN-pyrethroid
#      bioassays of each year (site-matched, weighted by mosquitoes tested):
#      the floor f, the susceptible fraction q, exp(u_s) C_t, C_t, the
#      reversion - t kappa, l0 and logit q, and the predicted p against the
#      observed pooled mortality (died / tested); and the same at fixed
#      weights (each of A's LLIN-pyrethroid cells and insecticides weighted
#      by its mosquitoes tested over the whole series), to follow one set of
#      sites through time
#   3. for wb_v5, a profile of the weighted binomial log likelihood (the
#      model's own, R/weighted_binomial.R, at the replicate rho) of A's
#      pyrethroid bioassays over two shifts constant across A: Ds added to
#      u_s and Df to the logit floor. Every other quantity is held at its
#      posterior mean at each bioassay (l0, exp(u_s) C_t, t kappa and the
#      logit floor; the posterior mean of the parameters themselves is not a
#      coherent point where the chains disagree on the smooths' sds), over
#      all chains and per chain; the plug-in's predictions are checked
#      against the posterior mean predictions. The optimum, its fit by window
#      (before 2010, 2010-15, 2016-18, 2019-24) and what the same shifts
#      would do to A's other classes and to B and C (Ds on every class, as
#      u_s; Df on the pyrethroids and DDT, as u_f)
#   4. selection by window: a shift Ds_w of u_s for each window, with the
#      floor free (one Df), at wb_v5's floor (Df = 0) and without a floor (f
#      = 0, as R/west_block_check.R's "no floor"), fitted to A's pyrethroid
#      bioassays by maximum likelihood, SEs from the Hessian: whether the
#      selection the data want falls relative to C_t from window to window,
#      i.e. whether C_t's timing is wrong
#   5. maps over West Africa of the posterior means of u_s (as the
#      multiplier) and of the pyrethroid floor plogis(floor_intercept + u_f)
#      (which, unlike u_f, compares between chains), per fit and chain set,
#      with regions A, B and C outlined (the continental maps of u_s and u_f
#      are R/smooth_maps.R's)
#
# Writes, in outputs/species_runs/west/:
#   v5_region_a_smooths.csv        1, per fit, chain set and block
#   v5_region_a_decomposition.csv  2, per fit, chain set, weighting and year
#   v5_region_a_profile_grid.csv   3, wb_v5's log likelihood over the grid
#   v5_region_a_profile.csv        3, per chain set: optimum and fit by window
#   v5_region_a_neighbours.csv     3, the change in log likelihood per block
#                                  and set of bioassays at A's optimum
#   v5_region_a_sub_profiles.csv   3, the profile on parts of A: all classes,
#                                  each country, the late bioassays
#   v5_region_a_plug_in.rds        3, wb_v5's plug-in terms and the bioassays
#   v5_region_a_window_selection.csv
#                                  4, per variant and window
# and in figures/species_runs/west/:
#   v5_decomposition_A.png  2, wb_v5 and V5 side by side
#   v5_profile_A.png        3 and 4
#   smooth_west.png         5
# Plain R with greta loaded (the draws are greta objects); one fit at a time;
# about 6.5 GB and 3 minutes.

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

scratchpad <- paste0("/tmp/claude-1000/-home-nick-Dropbox-github-ir-cube/",
                     "be75c64a-3bb7-4b3e-a81c-c664fe72f5e2/scratchpad/species")
fits <- c(wb_v5 = "outputs/pod_jobs/wb_v5/temporary/fitted_model.RData",
          V5 = file.path(scratchpad, "local_sp_v5/temporary/fitted_model.RData"))
# 250 draws per chain: four chains of wb_v5, two of V5 (USE_CHAINS)
n_draws <- c(wb_v5 = 1000, V5 = 500)
profile_fit <- "wb_v5"
llin_pyrethroids <- c("Alpha-cypermethrin", "Deltamethrin", "Permethrin")
windows <- c("pre-2010", "2010-15", "2016-18", "2019-24")
window_of <- function(year) {
  case_when(year < 2010 ~ windows[1], year <= 2015 ~ windows[2],
            year <= 2018 ~ windows[3], TRUE ~ windows[4])
}
years_shown <- 2005:2024
output_dir <- "outputs/species_runs/west"
figure_dir <- "figures/species_runs/west"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
options(width = 200)

regions <- read.csv("outputs/species_runs/regions/regions.csv")
blocks <- c("A", "B", "C")
block_names <- setNames(regions$region_name[match(blocks, regions$region)],
                        blocks)

# the grid of the West Africa maps: the cells of the water mask aggregated by
# 3 (about 14 km, as R/smooth_maps.R) in the maps' window
west_xlim <- c(-17.9, 24.5)
west_ylim <- c(4, 21)
west_grid <- terra::crop(
  terra::aggregate(rast("data/clean/pfpr_water_mask.tif"), 3, fun = "max",
                   na.rm = TRUE),
  terra::ext(west_xlim[1] - 1, west_xlim[2] + 1, west_ylim[1] - 1,
             west_ylim[2] + 1))
west_cells <- terra::cells(west_grid)
west_xy <- terra::xyFromCell(west_grid, west_cells)


# 0. each fit's terms at the bioassays of A, B and C -----------------------------

# The outer form's terms at every bioassay of A, B and C, draw by draw, as
# draws x bioassays matrices: l0, C (C_t), us (u_s, 0 where the class has no
# selection smooth), uf (u_f, 0 where it has no floor smooth), lf (the logit
# floor, floor_intercept + u_f) and reversion (t kappa, <= 0); and C at
# every year of each of A's LLIN-pyrethroid (cell, type) pairs, draws x pairs
# x years. With the bioassays (rows), the pairs, the chain of each draw and
# the fixed rho of each type.
fit_terms <- function(fit, label, n) {
  used <- even_draws(fit, n, label)
  parameters <- fit_parameter_draws(fit, used$index)
  options <- fit$options
  stopifnot(smooth_on(options), smooth_floor_on(options),
            !species_on(options), !kdr_on(options))
  D <- parameters$n_draws
  df <- fit$df
  rows <- df %>%
    mutate(row = row_number()) %>%
    left_join(transmute(regions, country, block = region),
              by = c("country_name" = "country")) %>%
    filter(block %in% blocks) %>%
    mutate(llin = insecticide_type %in% llin_pyrethroids,
           pyrethroid = insecticide_class == "Pyrethroids",
           floor_smooth = options$smooth$term_classes[class_id],
           window = window_of(year_start))
  pairs <- rows %>%
    filter(block == "A", llin) %>%
    group_by(cell_id, type_id) %>%
    summarise(tested = sum(mosquito_number), .groups = "drop") %>%
    mutate(pair = row_number())
  n_times <- max(fit$cell_years_index$year_id)
  mask_cells <- df$cell[match(seq_len(max(df$cell_id)), df$cell_id)]
  basis <- prediction_basis(options, mask_cells)
  cell_country <- dynamical_lookups(df)$cell_country_lookup
  x_row <- matrix(NA_integer_, max(fit$cell_years_index$cell_id), n_times)
  x_row[cbind(fit$cell_years_index$cell_id, fit$cell_years_index$year_id)] <-
    seq_len(nrow(fit$cell_years_index))
  empty <- matrix(NA_real_, D, nrow(rows))
  terms <- list(l0 = empty, C = empty, us = empty, uf = empty, lf = empty,
                reversion = empty)
  pair_C <- array(NA_real_, c(D, nrow(pairs), n_times))
  for (k in sort(unique(rows$type_id))) {
    at_k <- which(rows$type_id == k)
    cells <- sort(unique(rows$cell_id[at_k]))
    class <- fit$classes_index[k]
    x <- array(fit$x_cell_years[as.vector(x_row[cells, ]), , drop = FALSE],
               c(length(cells), n_times, ncol(fit$x_cell_years)))
    l0 <- cell_logit_init(parameters, k,
                          matrix(parameters$logit_init_relative[
                            , cell_country[cells], k], D),
                          parameters$x_cells_init[cells, , drop = FALSE])
    smooth_at <- function(kind) {
      if (smooth_class_weight(options, kind, class) == 1) {
        parameters$smooth_weights[[kind]] %*%
          t(basis[cells, , drop = FALSE])
      } else {
        matrix(0, D, length(cells))
      }
    }
    us <- smooth_at("selection")
    uf <- smooth_at("floor")
    intercept <- parameters$floor_intercept[
      , smooth_intercept_index(options, class)]
    kappa <- if (is.null(parameters$kappa_type)) rep(0, D) else
      parameters$kappa_type[, k]
    effect <- matrix(parameters$effect_type[, , k], D)
    column <- match(rows$cell_id[at_k], cells)
    terms$l0[, at_k] <- l0[, column]
    terms$us[, at_k] <- us[, column]
    terms$uf[, at_k] <- uf[, column]
    terms$lf[, at_k] <- intercept + uf[, column]
    terms$reversion[, at_k] <- outer(kappa, rows$year_id[at_k])
    pairs_k <- which(pairs$type_id == k)
    cumulative <- 0
    for (t in seq_len(n_times)) {
      cumulative <- cumulative +
        log1p(effect %*% t(matrix(x[, t, ], nrow = length(cells))))
      at_t <- at_k[rows$year_id[at_k] == t]
      terms$C[, at_t] <- cumulative[, match(rows$cell_id[at_t], cells)]
      if (length(pairs_k) > 0) {
        pair_C[, pairs_k, t] <- cumulative[, match(pairs$cell_id[pairs_k],
                                                   cells)]
      }
    }
  }
  # the outer form against the fit's own predictions
  logit_q <- terms$l0 - exp(terms$us) * terms$C - terms$reversion
  own <- dynamical_logit(parameters, rows, df, fit$x_cell_years,
                         fit$cell_years_index)
  ours <- floored_logit(logit_q, plogis(terms$lf))
  difference <- max(abs(pmin(pmax(own, -30), 30) - pmin(pmax(ours, -30), 30)))
  report("%s: terms at %d bioassays of A, B and C x %d draws; against the fit's predictions, max |logit diff| %.2g",
         label, nrow(rows), D, difference)
  stopifnot(difference < 1e-8)
  # per (cell, type) pair of A's LLIN pyrethroids, the terms that do not
  # change with time, from any of its bioassays
  first <- match(paste(pairs$cell_id, pairs$type_id),
                 paste(rows$cell_id, rows$type_id))
  # over the West Africa grid, per chain set: the posterior means of u_s and
  # of the pyrethroid floor plogis(floor_intercept + u_f), which, unlike u_f,
  # compares between chains that split the floor differently between the
  # two
  basis_west <- smooth_basis_at(options$smooth,
                                smooth_coords(west_xy[, 1], west_xy[, 2],
                                              options$smooth$crs))
  us_west <- parameters$smooth_weights$selection %*% t(basis_west)
  pyrethroids <- match("Pyrethroids", fit$classes)
  floor_west <- plogis(
    parameters$floor_intercept[, smooth_intercept_index(options, pyrethroids)] +
      parameters$smooth_weights$floor %*% t(basis_west))
  west <- lapply(chain_sets(list(label = label, chain = used$chain)),
                 function(i) {
                   list(u_s = colMeans(us_west[i, , drop = FALSE]),
                        floor = colMeans(floor_west[i, , drop = FALSE]))
                 })
  list(label = label, rows = rows, terms = terms, pairs = pairs, west = west,
       pair_C = pair_C,
       pair_terms = lapply(terms[c("l0", "us", "uf", "lf")],
                           function(m) m[, first, drop = FALSE]),
       pair_kappa = terms$reversion[, first, drop = FALSE] /
         rep(rows$year_id[first], each = D),
       chain = used$chain, baseline_year = fit$baseline_year,
       rho = c(parameters$rho_types[1, ]), types = fit$types,
       floor_intercept = parameters$floor_intercept,
       likelihood = fit$options$likelihood)
}

# the chain sets of a fit: all its usable chains, then each (wb_v5 only)
chain_sets <- function(x) {
  sets <- list(all = seq_along(x$chain))
  if (x$label == profile_fit) {
    for (chain in sort(unique(x$chain))) {
      sets[[paste0("chain ", chain)]] <- which(x$chain == chain)
    }
  }
  sets
}

# per draw, the means of the columns `columns` of draws x bioassays matrix
# `m` weighted by `weight`
weighted_by_draw <- function(m, columns, weight) {
  c(m[, columns, drop = FALSE] %*% weight[columns]) / sum(weight[columns])
}
posterior <- function(x) {
  c(mean = mean(x), sd = sd(x), lower = unname(quantile(x, 0.025)),
    upper = unname(quantile(x, 0.975)))
}

results <- list()
for (label in names(fits)) {
  fit <- load_fit(fits[[label]])
  results[[label]] <- fit_terms(fit, label, n_draws[[label]])
  rm(fit)
  gc(verbose = FALSE)
}


# 1. the smooths at the LLIN-pyrethroid bioassays --------------------------------

smooth_rows <- list()
for (label in names(results)) {
  x <- results[[label]]
  groups <- c(setNames(lapply(blocks, function(b) {
    which(x$rows$block == b & x$rows$llin)
  }), blocks),
  lapply(split(which(x$rows$block == "A" & x$rows$llin),
               x$rows$country_name[x$rows$block == "A" & x$rows$llin]),
         identity))
  names(groups)[-(1:3)] <- paste0("A: ", names(groups)[-(1:3)])
  weight <- x$rows$mosquito_number
  for (chain_set in names(chain_sets(x))) {
    draws <- chain_sets(x)[[chain_set]]
    for (group in names(groups)) {
      columns <- groups[[group]]
      quantities <- list(
        u_s = x$terms$us[draws, , drop = FALSE],
        u_f = x$terms$uf[draws, , drop = FALSE],
        multiplier = exp(x$terms$us[draws, , drop = FALSE]),
        floor = plogis(x$terms$lf[draws, , drop = FALSE]))
      for (quantity in names(quantities)) {
        by_draw <- weighted_by_draw(quantities[[quantity]], columns, weight)
        cell_means <- colMeans(quantities[[quantity]][, columns, drop = FALSE])
        w <- weight[columns] / sum(weight[columns])
        smooth_rows[[length(smooth_rows) + 1]] <- tibble(
          fit = label, chains = chain_set, group = group, quantity = quantity,
          bioassays = length(columns), cells = n_distinct(x$rows$cell[columns]),
          !!!as.list(posterior(by_draw)),
          spread = sqrt(sum(w * (cell_means - sum(w * cell_means)) ^ 2)))
      }
    }
    pyr <- x$floor_intercept[draws, 1]
    smooth_rows[[length(smooth_rows) + 1]] <- tibble(
      fit = label, chains = chain_set, group = "all", quantity = "pyrethroid floor intercept, plogis",
      bioassays = NA, cells = NA, !!!as.list(posterior(plogis(pyr))),
      spread = NA)
  }
}
smooth_table <- bind_rows(smooth_rows)
write.csv(smooth_table, file.path(output_dir, "v5_region_a_smooths.csv"),
          row.names = FALSE)
cat("\n1. the smooths at the LLIN-pyrethroid bioassays, weighted by mosquitoes tested (mean, sd, 95%; spread: sd between bioassays of the posterior mean)\n")
print(as.data.frame(smooth_table %>%
                      mutate(across(c(mean, sd, lower, upper, spread),
                                    ~ round(.x, 3)))), right = FALSE)


# 2. the decomposition by year in A ------------------------------------------------

# per draw, the weighted means of each term at bioassays `columns` (site-
# matched) or at pairs (fixed weights)
term_means <- function(terms, columns, weight) {
  q <- plogis(terms$l0 - exp(terms$us) * terms$C - terms$reversion)
  f <- plogis(terms$lf)
  quantities <- list(
    p = f + (1 - f) * q, f = f, q = q, above_floor = (1 - f) * q,
    selection = exp(terms$us) * terms$C, C = terms$C,
    multiplier = exp(terms$us), cost = -terms$reversion, l0 = terms$l0,
    logit_q = terms$l0 - exp(terms$us) * terms$C - terms$reversion)
  lapply(quantities, weighted_by_draw, columns = columns, weight = weight)
}
summarise_quantities <- function(means) {
  bind_rows(lapply(names(means), function(quantity) {
    tibble(quantity = quantity, !!!as.list(posterior(means[[quantity]])))
  }))
}
decomposition_rows <- list()
for (label in names(results)) {
  x <- results[[label]]
  weight <- x$rows$mosquito_number
  a_llin <- which(x$rows$block == "A" & x$rows$llin)
  for (chain_set in names(chain_sets(x))) {
    draws <- chain_sets(x)[[chain_set]]
    terms <- lapply(x$terms, function(m) m[draws, , drop = FALSE])
    for (year in years_shown) {
      columns <- a_llin[x$rows$year_start[a_llin] == year]
      if (length(columns) > 0) {
        decomposition_rows[[length(decomposition_rows) + 1]] <-
          summarise_quantities(term_means(terms, columns, weight)) %>%
          mutate(fit = label, chains = chain_set, weighting = "site-matched",
                 year = year, bioassays = length(columns),
                 tested = sum(weight[columns]),
                 observed = sum(x$rows$died[columns]) / sum(weight[columns]),
                 .before = 1)
      }
      # fixed weights: every pair at this year
      t <- year - x$baseline_year + 1
      pair_terms <- c(lapply(x$pair_terms, function(m) m[draws, , drop = FALSE]),
                      list(C = x$pair_C[draws, , t],
                           reversion = x$pair_kappa[draws, , drop = FALSE] * t))
      decomposition_rows[[length(decomposition_rows) + 1]] <-
        summarise_quantities(term_means(pair_terms, seq_len(nrow(x$pairs)),
                                        x$pairs$tested)) %>%
        mutate(fit = label, chains = chain_set, weighting = "fixed", year = year,
               bioassays = NA, tested = NA, observed = NA, .before = 1)
    }
  }
}
decomposition <- bind_rows(decomposition_rows)
write.csv(decomposition, file.path(output_dir, "v5_region_a_decomposition.csv"),
          row.names = FALSE)
cat("\n2. decomposition by year in A, LLIN pyrethroids, site-matched, posterior means (all usable chains)\n")
print(as.data.frame(decomposition %>%
                      filter(chains == "all", weighting == "site-matched") %>%
                      select(fit, year, bioassays, observed, quantity, mean) %>%
                      pivot_wider(names_from = quantity, values_from = mean) %>%
                      mutate(across(where(is.double), ~ round(.x, 3)))))
cat("\n   the same, at fixed weights\n")
print(as.data.frame(decomposition %>%
                      filter(chains == "all", weighting == "fixed") %>%
                      select(fit, year, quantity, mean) %>%
                      pivot_wider(names_from = quantity, values_from = mean) %>%
                      mutate(across(where(is.double), ~ round(.x, 3)))))
cat("\n   wb_v5 per chain, site-matched p, f and q by year\n")
print(as.data.frame(decomposition %>%
                      filter(fit == profile_fit, weighting == "site-matched",
                             quantity %in% c("p", "f", "q")) %>%
                      mutate(value = round(mean, 3)) %>%
                      select(year, observed, quantity, chains, value) %>%
                      pivot_wider(names_from = c(quantity, chains),
                                  values_from = value) %>%
                      mutate(observed = round(observed, 3))))


# 3. the profile of A's pyrethroids over shifts of u_s and the logit floor ----------

x <- results[[profile_fit]]
stopifnot(x$likelihood == "weighted_binomial")
rows <- x$rows
rows$weight <- design_effect_weight(rows$mosquito_number,
                                    x$rho[rows$type_id])

# the plug-in terms of a chain set: each term's posterior mean at each
# bioassay, with the selection term exp(u_s) C_t as one
plug_in <- function(draws) {
  m <- function(name) colMeans(x$terms[[name]][draws, , drop = FALSE])
  list(l0 = m("l0"),
       selection = colMeans(exp(x$terms$us[draws, , drop = FALSE]) *
                              x$terms$C[draws, , drop = FALSE]),
       reversion = m("reversion"), lf = m("lf"))
}
# logit q and the floor at bioassays `columns` with shifts ds (of u_s, on
# every class) and df (of the logit floor, on the classes with u_f), each
# one, or one per bioassay; or with the floor `floor` in place of the
# fit's, if given
plug_in_parts <- function(point, columns, ds = 0, df = 0, floor = NULL) {
  logit_q <- point$l0[columns] - exp(ds) * point$selection[columns] -
    point$reversion[columns]
  f <- if (is.null(floor)) {
    plogis(point$lf[columns] + df * rows$floor_smooth[columns])
  } else {
    rep_len(floor, length(columns))
  }
  list(logit_q = logit_q, f = f)
}
# predicted mortality at bioassays `columns` (plug_in_parts())
plug_in_p <- function(point, columns, ds = 0, df = 0, floor = NULL) {
  parts <- plug_in_parts(point, columns, ds, df, floor)
  parts$f + (1 - parts$f) * plogis(parts$logit_q)
}
# the weighted binomial log likelihood of each bioassay in `columns`, as the
# model's (R/weighted_binomial.R; plug_in_parts())
plug_in_log_lik <- function(point, columns, ds = 0, df = 0, floor = NULL) {
  parts <- plug_in_parts(point, columns, ds, df, floor)
  probs <- floored_log_probs(parts$logit_q,
                             if (all(parts$f == 0)) NULL else parts$f)
  weighted_binomial_log_lik(rows$died[columns], rows$mosquito_number[columns],
                            probs$log_p, probs$log_not_p,
                            rows$weight[columns])
}
a_pyr <- which(rows$block == "A" & rows$pyrethroid)
grid <- expand.grid(ds = seq(-0.6, 0.4, by = 0.01), df = seq(-1.2, 0.6, by = 0.02))

# the fit of A's bioassays by window: observed and predicted pooled
# mortality (the prediction weighted by mosquitoes tested) and the log
# likelihood, LLIN pyrethroids and all pyrethroids
window_fit <- function(p, columns, log_lik) {
  tibble(row = columns, p = p, log_lik = log_lik) %>%
    mutate(window = rows$window[columns], llin = rows$llin[columns],
           died = rows$died[columns], tested = rows$mosquito_number[columns]) %>%
    { bind_rows(mutate(., set = "all pyrethroids"),
                filter(., llin) %>% mutate(set = "LLIN pyrethroids")) } %>%
    group_by(set, window) %>%
    summarise(bioassays = n(), observed = sum(died) / sum(tested),
              predicted = sum(p * tested) / sum(tested),
              log_lik = sum(log_lik), .groups = "drop")
}

profiles <- list()
profile_grids <- list()
neighbours <- list()
window_fits <- list()
optima <- list()
for (chain_set in names(chain_sets(x))) {
  draws <- chain_sets(x)[[chain_set]]
  point <- plug_in(draws)
  # the plug-in against the posterior mean prediction at A's LLIN pyrethroids
  q_draws <- plogis(x$terms$l0[draws, a_pyr] -
                      exp(x$terms$us[draws, a_pyr]) * x$terms$C[draws, a_pyr] -
                      x$terms$reversion[draws, a_pyr])
  f_draws <- plogis(x$terms$lf[draws, a_pyr])
  posterior_p <- colMeans(f_draws + (1 - f_draws) * q_draws)
  base_p <- plug_in_p(point, a_pyr)
  objective <- function(theta) {
    -sum(plug_in_log_lik(point, a_pyr, theta[1], theta[2]))
  }
  starts <- expand.grid(ds = c(-1, 0, 1), df = c(-2, 0, 1))
  runs <- lapply(seq_len(nrow(starts)), function(i) {
    optim(unlist(starts[i, ]), objective, method = "BFGS",
          control = list(reltol = 1e-12, maxit = 1000))
  })
  best <- runs[[which.min(vapply(runs, `[[`, numeric(1), "value"))]]
  se <- sqrt(diag(solve(optimHess(best$par, objective))))
  optimum <- best$par
  optima[[chain_set]] <- optimum
  base_ll <- -objective(c(0, 0))
  profile_grids[[chain_set]] <- grid %>%
    mutate(chains = chain_set,
           log_lik = vapply(seq_len(nrow(grid)), function(i) {
             -objective(c(grid$ds[i], grid$df[i]))
           }, numeric(1)) - base_ll)
  window_fits[[chain_set]] <- bind_rows(
    window_fit(posterior_p, a_pyr, rep(NA_real_, length(a_pyr))) %>%
      mutate(at = "posterior mean prediction"),
    window_fit(base_p, a_pyr, plug_in_log_lik(point, a_pyr)) %>%
      mutate(at = "plug-in, no shift"),
    window_fit(plug_in_p(point, a_pyr, optimum[1], optimum[2]), a_pyr,
               plug_in_log_lik(point, a_pyr, optimum[1], optimum[2])) %>%
      mutate(at = "plug-in, optimum")) %>%
    mutate(chains = chain_set, .before = 1)
  profiles[[chain_set]] <- tibble(
    chains = chain_set, ds = optimum[1], se_ds = se[1], df = optimum[2],
    se_df = se[2], multiplier_shift = exp(optimum[1]),
    log_lik_gain = -best$value - base_ll,
    start_spread = diff(range(vapply(runs, `[[`, numeric(1), "value"))),
    plug_in_vs_posterior_max_pp = 100 * max(abs(
      (window_fits[[chain_set]] %>% filter(at == "plug-in, no shift"))$predicted -
        (window_fits[[chain_set]] %>% filter(at == "posterior mean prediction"))$predicted)))
  # the same shifts elsewhere: Ds on every class, Df on the pyrethroids and
  # DDT
  sets <- list(
    "LLIN pyrethroids" = rows$llin, "all pyrethroids" = rows$pyrethroid,
    "DDT" = rows$insecticide_type == "DDT",
    "other classes" = !rows$pyrethroid & rows$insecticide_type != "DDT")
  for (b in blocks) {
    for (s in names(sets)) {
      columns <- which(rows$block == b & sets[[s]])
      neighbours[[length(neighbours) + 1]] <- tibble(
        chains = chain_set, block = b, set = s, bioassays = length(columns),
        log_lik_change = sum(plug_in_log_lik(point, columns, optimum[1],
                                             optimum[2])) -
          sum(plug_in_log_lik(point, columns)),
        observed_2019 = with(rows[columns, ], sum(died[year_start >= 2019]) /
                               sum(mosquito_number[year_start >= 2019])),
        predicted_2019 = {
          late <- columns[rows$year_start[columns] >= 2019]
          sum(plug_in_p(point, late) * rows$mosquito_number[late]) /
            sum(rows$mosquito_number[late])
        },
        shifted_2019 = {
          late <- columns[rows$year_start[columns] >= 2019]
          sum(plug_in_p(point, late, optimum[1], optimum[2]) *
                rows$mosquito_number[late]) / sum(rows$mosquito_number[late])
        })
    }
  }
}
profiles <- bind_rows(profiles)
profile_grid <- bind_rows(profile_grids)
neighbours <- bind_rows(neighbours)
window_fits <- bind_rows(window_fits)
# a reference: the best constant mortality of A's pyrethroid bioassays in
# each window under the weighted binomial (its maximum, sum w y / sum w n),
# the same at every bioassay of the window
constant_p <- tapply(rows$weight[a_pyr] * rows$died[a_pyr],
                     rows$window[a_pyr], sum) /
  tapply(rows$weight[a_pyr] * rows$mosquito_number[a_pyr],
         rows$window[a_pyr], sum)
p_constant <- unname(constant_p[rows$window[a_pyr]])
window_fits <- bind_rows(
  window_fits,
  window_fit(p_constant, a_pyr,
             weighted_binomial_log_lik(rows$died[a_pyr],
                                       rows$mosquito_number[a_pyr],
                                       log(p_constant), log1p(-p_constant),
                                       rows$weight[a_pyr])) %>%
    mutate(at = "best constant per window", chains = "all"))
# wb_v5's plug-in terms over all chains, for further checks
saveRDS(list(rows = rows, point = plug_in(chain_sets(x)$all)),
        file.path(output_dir, "v5_region_a_plug_in.rds"))
# the same profile (all chains) on parts of A's bioassays: all classes (Ds on
# every class, Df on the pyrethroids and DDT, as the smooths), each country's
# pyrethroids, and the late pyrethroids alone
point <- plug_in(chain_sets(x)$all)
parts <- list(
  "A, pyrethroids" = a_pyr,
  "A, all classes" = which(rows$block == "A"),
  "Burkina Faso, pyrethroids" = which(rows$country_name == "Burkina Faso" &
                                        rows$pyrethroid),
  "Côte d'Ivoire, pyrethroids" = which(rows$country_name == "Côte d’Ivoire" &
                                         rows$pyrethroid),
  "A, pyrethroids 2016-24" = which(rows$block == "A" & rows$pyrethroid &
                                     rows$year_start >= 2016),
  "A, pyrethroids 2019-24" = which(rows$block == "A" & rows$pyrethroid &
                                     rows$year_start >= 2019))
sub_profiles <- bind_rows(lapply(names(parts), function(part) {
  columns <- parts[[part]]
  objective <- function(theta) {
    -sum(plug_in_log_lik(point, columns, theta[1], theta[2]))
  }
  runs <- lapply(list(c(0, 0), c(-1, -1), c(1, 0), c(-0.5, -2)), function(s) {
    optim(s, objective, method = "BFGS",
          control = list(reltol = 1e-12, maxit = 2000))
  })
  best <- runs[[which.min(vapply(runs, `[[`, numeric(1), "value"))]]
  se <- sqrt(diag(solve(optimHess(best$par, objective))))
  llin <- columns[rows$llin[columns] & rows$year_start[columns] >= 2019]
  tibble(part = part, bioassays = length(columns), ds = best$par[1],
         se_ds = se[1], df = best$par[2], se_df = se[2],
         log_lik_gain = -best$value + objective(c(0, 0)),
         llin_floor_2019 = if (length(llin) > 0) {
           sum(plogis(point$lf[llin] + best$par[2]) *
                 rows$mosquito_number[llin]) / sum(rows$mosquito_number[llin])
         } else NA_real_)
}))
write.csv(sub_profiles, file.path(output_dir, "v5_region_a_sub_profiles.csv"),
          row.names = FALSE)
write.csv(profile_grid, file.path(output_dir, "v5_region_a_profile_grid.csv"),
          row.names = FALSE)
write.csv(left_join(window_fits, profiles, by = "chains"),
          file.path(output_dir, "v5_region_a_profile.csv"), row.names = FALSE)
write.csv(neighbours, file.path(output_dir, "v5_region_a_neighbours.csv"),
          row.names = FALSE)
cat("\n3. profile of A's pyrethroid bioassays over (Ds, Df), wb_v5 at its posterior mean terms\n")
print(as.data.frame(profiles), digits = 3)
cat("\n   fit by window (observed / predicted pooled, log likelihood)\n")
print(as.data.frame(window_fits %>%
                      mutate(cell = ifelse(is.na(log_lik),
                                           sprintf("%.3f/%.3f", observed, predicted),
                                           sprintf("%.3f/%.3f %7.1f", observed,
                                                   predicted, log_lik))) %>%
                      select(chains, set, at, window, cell) %>%
                      pivot_wider(names_from = window, values_from = cell)),
      right = FALSE)
cat("\n   the same profile on parts of A (all chains; llin_floor_2019: the LLIN-pyrethroid floor at the optimum, 2019-24)\n")
print(as.data.frame(sub_profiles), digits = 3)
cat("\n   the change in log likelihood if A's optimal shifts applied elsewhere (and pooled 2019+ observed / predicted / shifted)\n")
print(as.data.frame(neighbours %>%
                      mutate(across(where(is.double), ~ round(.x, 3)))))


# 4. selection by window ------------------------------------------------------------

# A shift Ds_w of u_s in each window w, fitted to A's pyrethroid bioassays at
# wb_v5's terms (all chains), with the floor
#   floor free     one shift Df of the logit floor for all windows
#   wb_v5 floor    Df = 0
#   floor c        the floor fixed at c (0.15, 0.10, 0.05) at every bioassay
#   no floor       f = 0
#   window floors  a shift Df_w of the logit floor in each window too
# If C_t's timing fits A, the shifts are alike across windows at some floor;
# if they fall from window to window, C_t grows too fast late (or the decline
# came early) for one multiplier.
point <- plug_in(chain_sets(x)$all)
window_index <- match(rows$window[a_pyr], windows)
n_windows <- length(windows)
llin_a <- rows$llin[a_pyr]
variants <- list("floor free" = list(df = "one"),
                 "wb_v5 floor" = list(df = "none"),
                 "floor 0.15" = list(floor = 0.15),
                 "floor 0.10" = list(floor = 0.10),
                 "floor 0.05" = list(floor = 0.05),
                 "no floor" = list(floor = 0),
                 "window floors" = list(df = "window"))
by_window <- list()
by_window_p <- list()
for (variant in names(variants)) {
  setting <- variants[[variant]]
  n_df <- switch(coalesce(setting$df, "none"), none = 0, one = 1,
                 window = n_windows)
  # ds per bioassay, df per bioassay (or 0), and the fixed floor (or NULL)
  shifts <- function(theta) {
    ds <- theta[window_index]
    df <- switch(as.character(n_df), "0" = 0,
                 "1" = theta[n_windows + 1],
                 theta[n_windows + window_index])
    list(ds = ds, df = df)
  }
  objective <- function(theta) {
    s <- shifts(theta)
    -sum(plug_in_log_lik(point, a_pyr, s$ds, s$df, setting$floor))
  }
  starts <- lapply(c(-1, -0.5, 0, 0.5), function(v) {
    c(rep(v, n_windows), rep(0, n_df))
  })
  runs <- lapply(starts, function(start) {
    optim(start, objective, method = "BFGS",
          control = list(reltol = 1e-12, maxit = 5000))
  })
  values <- vapply(runs, `[[`, numeric(1), "value")
  best <- runs[[which.min(values)]]
  hessian <- optimHess(best$par, objective)
  se <- tryCatch(sqrt(pmax(diag(solve(hessian)), 0)),
                 error = function(e) rep(NA_real_, length(best$par)))
  s <- shifts(best$par)
  p <- plug_in_p(point, a_pyr, s$ds, s$df, setting$floor)
  f <- plug_in_parts(point, a_pyr, s$ds, s$df, setting$floor)$f
  ll <- plug_in_log_lik(point, a_pyr, s$ds, s$df, setting$floor)
  base_ll <- plug_in_log_lik(point, a_pyr)
  tested <- rows$mosquito_number[a_pyr]
  by_window[[variant]] <- tibble(
    variant = variant, window = windows,
    ds = best$par[seq_len(n_windows)], se_ds = se[seq_len(n_windows)],
    multiplier_shift = exp(ds),
    # the tested-weighted exp(u_s + Ds_w) at the window's LLIN pyrethroids
    multiplier = vapply(seq_len(n_windows), function(w) {
      columns <- a_pyr[window_index == w & llin_a]
      sum(exp(colMeans(x$terms$us[, columns, drop = FALSE]) +
                best$par[w]) * rows$mosquito_number[columns]) /
        sum(rows$mosquito_number[columns])
    }, numeric(1)),
    df = vapply(seq_len(n_windows), function(w) {
      if (n_df == 0) NA_real_ else if (n_df == 1) best$par[n_windows + 1] else
        best$par[n_windows + w]
    }, numeric(1)),
    se_df = vapply(seq_len(n_windows), function(w) {
      if (n_df == 0) NA_real_ else if (n_df == 1) se[n_windows + 1] else
        se[n_windows + w]
    }, numeric(1)),
    # LLIN pyrethroids, by window: the floor, observed and predicted pooled
    # mortality, and the log likelihood (all pyrethroids) against wb_v5's
    llin_floor = vapply(seq_len(n_windows), function(w) {
      i <- window_index == w & llin_a
      sum(f[i] * tested[i]) / sum(tested[i])
    }, numeric(1)),
    observed = vapply(seq_len(n_windows), function(w) {
      i <- window_index == w & llin_a
      sum(rows$died[a_pyr][i]) / sum(tested[i])
    }, numeric(1)),
    predicted = vapply(seq_len(n_windows), function(w) {
      i <- window_index == w & llin_a
      sum(p[i] * tested[i]) / sum(tested[i])
    }, numeric(1)),
    log_lik_gain_window = vapply(seq_len(n_windows), function(w) {
      sum(ll[window_index == w]) - sum(base_ll[window_index == w])
    }, numeric(1)),
    log_lik_gain = sum(ll) - sum(base_ll),
    start_spread = diff(range(values)))
  by_window_p[[variant]] <- tibble(row = a_pyr, p = p, variant = variant)
}
by_window <- bind_rows(by_window)
write.csv(by_window, file.path(output_dir, "v5_region_a_window_selection.csv"),
          row.names = FALSE)
cat("\n4. selection by window: shifts Ds_w of u_s (multiplier_shift = exp(Ds_w), relative to wb_v5's exp(u_s) C_t), A's pyrethroids; LLIN floor, observed and predicted pooled, log likelihood gain over wb_v5\n")
print(as.data.frame(by_window %>% mutate(across(where(is.double), ~ signif(.x, 3)))))

# figures ---------------------------------------------------------------------------

fit_titles <- c(wb_v5 = "wb_v5 (weighted binomial, all 4 chains)",
                V5 = "V5 (beta-binomial, chains 1 and 2)")

# the decomposition
observed <- decomposition %>%
  filter(chains == "all", weighting == "site-matched", quantity == "p") %>%
  select(fit, year, observed, tested)
panel_theme <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 10))
decomposition_panel <- function(label) {
  d <- filter(decomposition, fit == label, chains == "all")
  site <- filter(d, weighting == "site-matched")
  fixed <- filter(d, weighting == "fixed")
  shade <- annotate("rect", xmin = c(2009.5, 2018.5), xmax = c(2015.5, 2024.5),
                    ymin = -Inf, ymax = Inf, fill = grey(0.95))
  mortality <- ggplot(mapping = aes(year)) + shade +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                data = filter(site, quantity == "p"), fill = "#009E73",
                alpha = 0.2) +
    geom_line(aes(y = mean, colour = "predicted p", linetype = "site-matched"),
              data = filter(site, quantity == "p"), linewidth = 0.8) +
    geom_line(aes(y = mean, colour = "predicted p", linetype = "fixed weights"),
              data = filter(fixed, quantity == "p"), linewidth = 0.5) +
    geom_line(aes(y = mean, colour = "floor f", linetype = "site-matched"),
              data = filter(site, quantity == "f"), linewidth = 0.8) +
    geom_line(aes(y = mean, colour = "(1 - f) q", linetype = "site-matched"),
              data = filter(site, quantity == "above_floor"), linewidth = 0.8) +
    geom_point(aes(y = observed, size = tested),
               data = filter(observed, fit == label), shape = 21,
               fill = grey(0.8), colour = grey(0.1), stroke = 0.4) +
    scale_colour_manual(values = c("predicted p" = "#009E73",
                                   "floor f" = "#D55E00",
                                   "(1 - f) q" = "#0072B2"),
                        name = NULL) +
    scale_linetype_manual(values = c("site-matched" = "solid",
                                     "fixed weights" = "22"), name = NULL) +
    scale_size_area(max_size = 4, name = "mosquitoes tested",
                    breaks = c(2000, 10000, 20000)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = NULL, y = "mortality",
         title = sprintf("%s\nmortality p = f + (1 - f) q", fit_titles[[label]])) +
    panel_theme
  susceptible <- ggplot(mapping = aes(year)) + shade +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                data = filter(site, quantity == "q"), fill = "#0072B2",
                alpha = 0.2) +
    geom_line(aes(y = mean, linetype = "site-matched"),
              data = filter(site, quantity == "q"), colour = "#0072B2",
              linewidth = 0.8) +
    geom_line(aes(y = mean, linetype = "fixed weights"),
              data = filter(fixed, quantity == "q"), colour = "#0072B2",
              linewidth = 0.5) +
    scale_linetype_manual(values = c("site-matched" = "solid",
                                     "fixed weights" = "22"), name = NULL) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = NULL, y = "susceptible fraction q",
         title = "susceptible fraction q (mortality above the floor)") +
    guides(linetype = "none") + panel_theme
  logit_terms <- c(l0 = "l0, logit initial state",
                   selection = "exp(u_s) C_t, selection",
                   C = "C_t alone",
                   cost = "-t kappa, reversion",
                   logit_q = "logit q = l0 - exp(u_s) C_t + (-t kappa)")
  logit_colours <- setNames(c(grey(0.4), "#D55E00", "#E69F00", "#CC79A7",
                              "#0072B2"), logit_terms)
  logit <- ggplot(mapping = aes(year)) + shade +
    geom_hline(yintercept = 0, colour = grey(0.6), linewidth = 0.3) +
    geom_line(aes(y = mean, colour = unname(logit_terms[quantity]),
                  linetype = "site-matched"),
              data = filter(site, quantity %in% names(logit_terms)),
              linewidth = 0.8) +
    geom_line(aes(y = mean, colour = unname(logit_terms[quantity]),
                  linetype = "fixed weights"),
              data = filter(fixed, quantity %in% names(logit_terms)),
              linewidth = 0.5) +
    scale_colour_manual(values = logit_colours, breaks = unname(logit_terms),
                        name = NULL) +
    scale_linetype_manual(values = c("site-matched" = "solid",
                                     "fixed weights" = "22"), name = NULL) +
    labs(x = NULL, y = "logit scale",
         title = "the terms of logit q") +
    guides(linetype = "none") + panel_theme
  list(mortality, susceptible, logit)
}
panels <- lapply(names(results), decomposition_panel)
figure <- wrap_plots(c(panels[[1]], panels[[2]]), ncol = 2, byrow = FALSE) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = sprintf("Region A (%s): what makes V5's LLIN-pyrethroid mortality level off",
                    sub("^A: ", "", block_names[["A"]])),
    caption = paste(
      "Solid: means over each year's LLIN-pyrethroid bioassays in A (with Liberia, attached to A), weighted by mosquitoes tested",
      "(site-matched), of each term's posterior draws; bands 95%.\nDashed: the same at fixed weights (each cell and insecticide",
      "of A by its mosquitoes tested over the whole series), following one set of sites. Points: observed pooled mortality",
      "(died / tested).\nlogit q = l0 - exp(u_s) C_t - t kappa (kappa <= 0; the pyrethroids' is about 0); f = plogis(floor_intercept + u_f).",
      "Shaded: 2010-2015 and 2019-2024."),
    theme = theme(plot.caption = element_text(hjust = 0, size = 8.5)))
ggsave(file.path(figure_dir, "v5_decomposition_A.png"), figure, width = 12,
       height = 11, dpi = 130, bg = "white")

# the profile and selection by window
heat <- profile_grid %>% filter(chains == "all")
heat_panel <- ggplot(heat, aes(ds, df)) +
  geom_raster(aes(fill = pmax(log_lik, -40))) +
  geom_contour(aes(z = log_lik), breaks = c(-30, -20, -10, -5, 0, 5, 10),
               colour = grey(0.2), linewidth = 0.25) +
  geom_point(aes(ds, df, shape = chains), data = profiles, size = 2) +
  annotate("point", x = 0, y = 0, shape = 4, size = 3, colour = "white") +
  scale_fill_viridis_c(name = "log likelihood\nchange (>= -40;\ncontours -30, -20,\n-10, -5, 0, 5, 10)") +
  scale_shape_manual(values = c(all = 16, "chain 1" = 1, "chain 2" = 2,
                                "chain 3" = 0, "chain 4" = 5), name = "optimum") +
  labs(x = "Ds, shift of u_s (log multiplier)", y = "Df, shift of the logit floor",
       title = "a. wb_v5: A's pyrethroid log likelihood over block shifts\n(cross: no shift; points: optimum, all chains and per chain)") +
  theme_minimal(base_size = 10)
yearly <- function(p, columns, name) {
  tibble(row = columns, p = p) %>%
    mutate(year = rows$year_start[columns], tested = rows$mosquito_number[columns],
           llin = rows$llin[columns]) %>%
    filter(llin, year %in% years_shown) %>%
    group_by(year) %>%
    summarise(predicted = sum(p * tested) / sum(tested), .groups = "drop") %>%
    mutate(variant = name)
}
optimum <- optima$all
shown <- c("floor free", "floor 0.10", "no floor", "window floors")
lines <- bind_rows(
  yearly(plug_in_p(point, a_pyr), a_pyr, "wb_v5 (plug-in, no shift)"),
  yearly(plug_in_p(point, a_pyr, optimum[1], optimum[2]), a_pyr,
         sprintf("block optimum (Ds %.2f, Df %.2f)", optimum[1], optimum[2])),
  bind_rows(lapply(shown, function(v) {
    yearly(by_window_p[[v]]$p, a_pyr, paste("by window,", v))
  })))
line_colours <- setNames(c(grey(0.3), "#D55E00", "#0072B2", "#009E73",
                           "#CC79A7", "#E69F00"),
                         unique(lines$variant))
trend_panel <- ggplot(mapping = aes(year)) +
  annotate("rect", xmin = c(2009.5, 2018.5), xmax = c(2015.5, 2024.5),
           ymin = -Inf, ymax = Inf, fill = grey(0.95)) +
  geom_line(aes(y = predicted, colour = variant), data = lines,
            linewidth = 0.7) +
  geom_point(aes(y = observed, size = tested),
             data = filter(observed, fit == profile_fit), shape = 21,
             fill = grey(0.8), colour = grey(0.1), stroke = 0.4) +
  scale_colour_manual(values = line_colours, name = NULL) +
  scale_size_area(max_size = 4, name = "mosquitoes tested",
                  breaks = c(2000, 10000, 20000)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = NULL, y = "LLIN-pyrethroid mortality",
       title = "b. A's LLIN pyrethroids: observed and site-matched predictions") +
  theme_minimal(base_size = 10)
variant_colours <- c("floor free" = "#0072B2", "wb_v5 floor" = "#56B4E9",
                     "floor 0.15" = "#E69F00", "floor 0.10" = "#009E73",
                     "floor 0.05" = "#F0E442", "no floor" = "#CC79A7",
                     "window floors" = "#D55E00")
window_panel <- by_window %>%
  mutate(variant = factor(variant, names(variant_colours)),
         window = factor(window, windows),
         lower = exp(ds - 1.96 * se_ds), upper = exp(ds + 1.96 * se_ds)) %>%
  ggplot(aes(window, multiplier_shift, colour = variant)) +
  geom_hline(yintercept = 1, colour = grey(0.5)) +
  geom_pointrange(aes(ymin = lower, ymax = upper),
                  position = position_dodge(width = 0.6), size = 0.25) +
  scale_y_log10(limits = c(0.05, 20), oob = scales::squish) +
  scale_colour_manual(values = variant_colours, name = NULL) +
  labs(x = NULL, y = "selection multiplier exp(Ds_w)\nrelative to wb_v5's exp(u_s) C_t",
       title = "c. selection by window (95% from the Hessian; squished at 0.05 and 20)") +
  theme_minimal(base_size = 10)
gains <- bind_rows(
  tibble(variant = "block optimum", log_lik_gain = profiles$log_lik_gain[
    profiles$chains == "all"]),
  by_window %>% distinct(variant, log_lik_gain)) %>%
  mutate(variant = factor(variant, c("block optimum", names(variant_colours))))
gain_panel <- ggplot(gains, aes(log_lik_gain, variant)) +
  geom_col(fill = grey(0.6), width = 0.6) +
  geom_vline(xintercept = 0, colour = grey(0.3)) +
  geom_text(aes(label = sprintf("%+.0f", log_lik_gain),
                hjust = ifelse(log_lik_gain < 0, 1.1, -0.1)), size = 3) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  labs(x = "log likelihood gain over wb_v5 (A's pyrethroids)", y = NULL,
       title = "d. what each variant gains") +
  theme_minimal(base_size = 10)
figure <- (heat_panel | trend_panel) / (window_panel | gain_panel) +
  plot_annotation(
    title = "Region A: does a block shift of wb_v5's smooths, or selection by window, fit the continued decline?",
    caption = paste(
      "wb_v5's terms held at their posterior means at each bioassay (l0, exp(u_s) C_t, t kappa, the logit floor); weighted binomial",
      "log likelihood of A's pyrethroid bioassays at the replicate rho.\nDs shifts u_s (every class), Df the logit floor",
      "(pyrethroids and DDT). By window: one Ds per window (before 2010, 2010-15, 2016-18, 2019-24), with the floor free",
      "(one Df), at wb_v5's,\nfixed at 0.15, 0.10 or 0.05, none, or one Df per window too."),
    theme = theme(plot.caption = element_text(hjust = 0, size = 8.5)))
ggsave(file.path(figure_dir, "v5_profile_A.png"), figure, width = 14,
       height = 10, dpi = 130, bg = "white")


# 5. maps of the smooths over West Africa --------------------------------------------

# per fit and chain set, the posterior mean of u_s (as the multiplier) and of
# the pyrethroid floor over the West Africa grid (fit_terms())
sf_use_s2(FALSE)
borders <- readRDS("data/clean/country_borders.RDS")
outlines <- suppressMessages(
  borders %>%
    inner_join(transmute(regions, country_name = country, block = region),
               by = "country_name") %>%
    filter(block %in% blocks) %>%
    st_make_valid() %>%
    group_by(block) %>%
    summarise(geometry = st_union(geometry), .groups = "drop"))
missing <- setdiff(regions$country[regions$region %in% blocks],
                   borders$country_name)
if (length(missing) > 0) report("no borders for %s", toString(missing))
centres <- suppressWarnings(st_point_on_surface(outlines))
map_sets <- bind_rows(lapply(names(results), function(label) {
  tibble(fit = label, chains = names(results[[label]]$west))
}))
map_panel <- function(label, set, kind) {
  layer <- rast(west_grid)
  layer[] <- NA_real_
  layer[west_cells] <- results[[label]]$west[[set]][[kind]]
  scale <- if (kind == "u_s") {
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"),
                         limits = c(-2, 2), oob = scales::squish,
                         na.value = "transparent",
                         labels = function(v) sprintf("%.2g", exp(v)),
                         name = "selection\nmultiplier\nexp(u_s)")
  } else {
    scale_fill_viridis_c(option = "magma", limits = c(0, 0.6),
                         oob = scales::squish, na.value = "transparent",
                         labels = scales::percent,
                         name = "pyrethroid\nfloor f")
  }
  ggplot() +
    geom_sf(data = borders, fill = grey(0.92), colour = NA) +
    geom_spatraster(data = layer) +
    geom_sf(data = borders, fill = NA, colour = grey(0.5), linewidth = 0.15) +
    geom_sf(data = outlines, fill = NA, colour = "black", linewidth = 0.6) +
    geom_sf_label(aes(label = block), data = centres, size = 3.5,
                  fontface = "bold", label.size = 0, alpha = 0.7) +
    scale +
    coord_sf(xlim = west_xlim, ylim = west_ylim, expand = FALSE) +
    labs(title = sprintf("%s, %s: %s", label,
                         if (set == "all") "all usable chains" else set,
                         if (kind == "u_s") "selection multiplier exp(u_s)" else
                           "pyrethroid floor plogis(intercept + u_f)"),
         x = NULL, y = NULL) +
    theme_ir_maps() +
    theme(plot.title = element_markdown(size = 10))
}
map_panels <- list()
for (i in seq_len(nrow(map_sets))) {
  for (kind in c("u_s", "floor")) {
    map_panels[[length(map_panels) + 1]] <- map_panel(map_sets$fit[i],
                                                      map_sets$chains[i], kind)
  }
}
figure <- wrap_plots(map_panels, ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "The latent smooths over West Africa: posterior means, wb_v5 (all chains and each) and V5 (chains 1 and 2)",
    caption = paste(
      "Posterior means over the draws of each chain set (250 per chain). Left: u_s as the multiplier exp(u_s) of the cumulative",
      "log fitness, clipped at 0.14 and 7.4 (blue: stronger selection);\nright: the pyrethroid floor plogis(floor_intercept + u_f),",
      "which, unlike u_f, compares between chains. Regions A, B and C of outputs/species_runs/regions/regions.csv outlined."),
    theme = theme(plot.caption = element_text(hjust = 0, size = 8.5)))
ggsave(file.path(figure_dir, "smooth_west.png"), figure, width = 12,
       height = 15, dpi = 110, bg = "white")
report("written %s and figures in %s; peak memory %.1f GB",
       file.path(output_dir, "v5_region_a_*.csv"), figure_dir, peak_memory_gb())
