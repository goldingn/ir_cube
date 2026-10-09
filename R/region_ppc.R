# Posterior predictive checks of the regional LLIN-pyrethroid diagnostics of
# R/region_diagnostics.R (#47), and a check of the beta-binomial error model:
# where a model misses a region's observed pooled mortality, is the miss
# larger than its predictive spread, and is it a miss of the mean or of the
# shape of the bioassays' distribution? Does the overdispersion, fixed per
# insecticide type, pull the fitted mean away from the data where mortality
# is low (or high)?
#
#   USE_CHAINS="V5=1,2" Rscript R/region_ppc.R [<label>=<fitted_model.RData> ...]
#
# The fits are ref_f0, V3f, V4_class and V5 (`fits` below, as
# R/region_trends.R), or those given. The bioassays are the modelled
# LLIN-pyrethroid bioassays (alpha-cypermethrin, deltamethrin, permethrin) of
# every year, each country's region from
# outputs/species_runs/regions/regions.csv (R/diagnostic_regions.R).
#
# 1. Posterior predictive checks. For each fit, from about n_draws (200)
# posterior draws, the same number from each usable chain (even_draws(),
# R/species_fit_helpers.R; with USE_CHAINS, as R/species_misfit.R used for
# V5), the predicted mortality p_i at every bioassay (dynamical_logit(),
# R/dynamical_predictions.R, as R/species_misfit.R) and the draw's rho per
# type.
# Cached in outputs/species_runs/regions/ppc_draws_<label>.rds (recomputed
# when the fit is newer). Per draw, one replicate of every bioassay from the
# model's likelihood (betabinomial_p_rho(), R/functions.R, rho the intra-class
# correlation): q_i ~ Beta(p_i phi, (1 - p_i) phi) with phi = 1 / rho_type - 1,
# then died_i ~ Binomial(n_i, q_i), with p_i clamped to [1e-12, 1 - 1e-12] as
# R/species_misfit.R. Per region and year, and per region and window
# (2010-2015, 2016-2018, 2019-2024; also for all of Africa, and for region A
# with Ghana), the statistics
#   pooled   died over tested
#   mean     the unweighted mean of died / n
#   median   the median of died / n
#   le05     the share of bioassays with died / n <= 0.05
#   gt30     the share with died / n > 0.30
#   ge98     the share with died / n >= 0.98
# of the observed bioassays and of each replicate: the observed statistic, the
# predictive median and 50% and 90% intervals, and the two-sided tail
# probability 2 min(u, 1 - u), u = P(T_rep < T_obs) + P(T_rep = T_obs) / 2
# (mid-p, for the ties of the shares). For pooled and mean, `expected` is the
# posterior mean of the site-matched prediction (sum n_i p_i / sum n_i, or the
# mean p_i), as R/region_diagnostics.R. Per region and window, the same check
# mean-matched (check "mean-matched"; "as fitted" otherwise): each draw's
# predictions at the region-window's bioassays shifted by one constant on the
# logit scale so that their expected pooled mortality is the observed, which
# leaves the shape of the predictive distribution at the right mean.
#
# 2. The error model. Per region and window (and region A with Ghana), per
# insecticide type and for the three types together (a mean per type), a
# constant-mean beta-binomial fitted by maximum likelihood (i) with rho fixed
# at V4_class's posterior mean rho for the type (`fixed_rho_fit`), (ii) with
# rho free (one shift of the types' logit rho: one type's rho is free, and
# with the types together their ratios stay V4_class's, so that (i) is nested
# in (ii)): the fitted mean (types together: the mean of the fitted means over
# the bioassays, to set against the unweighted observed mean), the free rho
# with a 90% Wald interval on the logit scale, and the likelihood-ratio
# statistic of the free rho. A constant mean puts the spread of the true
# mortality between the sites into rho; so, with the types together, the
# same pair of fits with V4_class's posterior mean prediction as an offset
# instead of a constant mean, logit mu_i = logit p_i + c_type (offset_*),
# which leaves rho only what the model's own predictions leave. Then (iii) a
# mean-dependent rho, rho_k(mu) = plogis(a_k + b log(4 mu (1 - mu))), fitted
# across the region-window-type cells, each with its own constant mean
# (mean_dependent), and a weighted quadratic fit of the logit free rho on the
# observed mean.
#
# Writes
#   figures/species_runs/regions/region_ppc_pooled.png  observed pooled
#       mortality per region and year, and each model's 50% and 90%
#       predictive intervals
#   figures/species_runs/regions/region_ppc_shape.png  per window, the
#       pooled mortality and the shape statistics (median, le05, gt30),
#       observed against the predictive intervals, per region and model
#   figures/species_runs/regions/region_ppc_shape_matched.png  the same,
#       mean-matched (median, le05, gt30, ge98)
#   figures/species_runs/regions/region_error_model.png  the fitted means
#       under fixed, free and mean-dependent rho against the observed, and
#       the free rho against the observed mean
#   outputs/species_runs/regions/ppc_years.csv     per region, year, model
#                                                  and statistic
#   outputs/species_runs/regions/ppc_windows.csv   per region, window, model,
#                                                  check and statistic
#   outputs/species_runs/regions/error_model.csv          per region, window
#                                                         and type
#   outputs/species_runs/regions/error_model_pooled.csv   the types together
# Plain R with greta loaded (the draws are greta objects); one fit loaded at a
# time; about 5 GB at peak (V3f) and 2 minutes, then 35 s from the caches.

suppressMessages({
  library(greta)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(ggtext)
})
source("R/functions.R")
source("R/dynamical_predictions.R")
source("R/species_fit_helpers.R")

n_draws <- 200
seed <- 47
scratchpad <- paste0("/tmp/claude-1000/-home-nick-Dropbox-github-ir-cube/",
                     "be75c64a-3bb7-4b3e-a81c-c664fe72f5e2/scratchpad/species")
fits <- c(
  ref_f0 = paste0("../ir_cube_netscreen/outputs/pod_jobs/dh270_lin_f0_full/",
                  "temporary/fitted_model.RData"),
  V3f = "outputs/pod_jobs/sp_v3_floor/temporary/fitted_model.RData",
  V4_class = file.path(scratchpad,
                       "local_sp_v4_class/temporary/fitted_model.RData"),
  V5 = file.path(scratchpad, "local_sp_v5/temporary/fitted_model.RData"))
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) > 0) {
  fits <- character(0)
  for (argument in arguments) {
    parts <- strsplit(argument, "=", fixed = TRUE)[[1]]
    stopifnot(length(parts) == 2)
    fits[[parts[1]]] <- parts[2]
  }
}
stopifnot(all(file.exists(fits)))
labels <- names(fits)
# the fit whose posterior mean rho per type the error-model check fixes
fixed_rho_fit <- Sys.getenv("FIXED_RHO_FIT", "V4_class")
stopifnot(fixed_rho_fit %in% labels)
# Okabe-Ito, as R/region_diagnostics.R
fit_colours <- c(ref_f0 = grey(0.3), V3f = "#E69F00", V4_class = "#CC79A7",
                 V5 = "#009E73")
fit_colours <- c(fit_colours, V4 = "#0072B2")
llin_pyrethroids <- c("Alpha-cypermethrin", "Deltamethrin", "Permethrin")
windows <- list(`2010-2015` = 2010:2015, `2016-2018` = 2016:2018,
                `2019-2024` = 2019:2024)
statistic_names <- c(pooled = "pooled mortality",
                     mean = "unweighted mean of died / n",
                     median = "median of died / n",
                     le05 = "share with died / n <= 0.05",
                     gt30 = "share with died / n > 0.30",
                     ge98 = "share with died / n >= 0.98")
output_dir <- "outputs/species_runs/regions"
figure_dir <- "figures/species_runs/regions"
regions <- read.csv(file.path(output_dir, "regions.csv"))
region_counts <- read.csv(file.path(output_dir, "region_counts.csv"))
region_letters <- sort(unique(regions$region))
# the extra groups of the window summaries: region A with Ghana (where the
# low-mortality misfit was first seen), and all regions
extra_groups <- c(`A+Ghana` = "A + Ghana", Africa = "Africa")
dpi <- 120
options(width = 220)
clamp <- function(x) pmin(pmax(x, 1e-12), 1 - 1e-12)


# per fit: the posterior draws of p_i and rho at the bioassays --------------------

fit_draws <- function(label, file) {
  cache <- file.path(output_dir, sprintf("ppc_draws_%s.rds", label))
  if (file.exists(cache) && file.mtime(cache) > file.mtime(file)) {
    cached <- readRDS(cache)
    if (cached$n_draws_requested == n_draws &&
        identical(cached$chains_requested, named_chains(label))) {
      return(cached)
    }
  }
  fit <- load_fit(file)
  df <- fit$df
  llin <- which(df$insecticide_type %in% llin_pyrethroids)
  chosen <- even_draws(fit, n_draws, label)
  parameters <- fit_parameter_draws(fit, chosen$index)
  time <- system.time(
    logit <- dynamical_logit(parameters, df[llin, ], df, fit$x_cell_years,
                             fit$cell_years_index)
  )[["elapsed"]]
  report("%s: predictions at %d LLIN-pyrethroid bioassays x %d draws in %.0f s",
         label, ncol(logit), nrow(logit), time)
  rho <- parameters$rho_types
  colnames(rho) <- fit$types
  out <- list(bioassays = as.data.frame(df[llin, c("cell", "year_start",
                                                   "country_name",
                                                   "insecticide_type",
                                                   "type_id", "died",
                                                   "mosquito_number")]),
              p = plogis(logit), rho = rho,
              chains = sort(unique(chosen$chain)),
              n_draws_requested = n_draws,
              chains_requested = named_chains(label))
  saveRDS(out, cache)
  report("%s saved; peak memory %.1f GB", label, peak_memory_gb())
  rm(fit, df, parameters, logit)
  invisible(gc())
  out
}

draws <- lapply(setNames(nm = labels), function(label) {
  fit_draws(label, fits[[label]])
})

# the bioassays: the same in every fit (up to the type ids)
bioassays <- draws[[1]]$bioassays
for (label in labels) {
  b <- draws[[label]]$bioassays
  stopifnot(identical(b$cell, bioassays$cell),
            identical(b$year_start, bioassays$year_start),
            identical(b$insecticide_type, bioassays$insecticide_type),
            identical(b$died, bioassays$died),
            identical(b$mosquito_number, bioassays$mosquito_number))
}
bioassays <- bioassays %>%
  mutate(region = regions$region[match(country_name, regions$country)],
         window = NA_character_)
for (w in names(windows)) {
  bioassays$window[bioassays$year_start %in% windows[[w]]] <- w
}
stopifnot(!anyNA(bioassays$region))
died <- bioassays$died
n <- bioassays$mosquito_number

# posterior mean rho per type, per fit
rho_means <- bind_rows(lapply(labels, function(label) {
  rho_draws <- draws[[label]]$rho[, llin_pyrethroids]
  tibble(model = label, insecticide_type = llin_pyrethroids,
         rho = colMeans(rho_draws),
         lower = apply(rho_draws, 2, quantile, 0.05),
         upper = apply(rho_draws, 2, quantile, 0.95))
}))
cat("\nposterior mean rho (90% interval) of the LLIN pyrethroids\n")
print(as.data.frame(rho_means), digits = 3)


# the groups --------------------------------------------------------------------------

# each group's bioassays, as rows of `bioassays`
year_groups <- bioassays %>%
  mutate(index = row_number()) %>%
  group_by(region, year = year_start) %>%
  summarise(index = list(index), .groups = "drop")
window_groups <- bind_rows(
  bioassays %>%
    mutate(index = row_number()) %>%
    filter(!is.na(window)) %>%
    group_by(region, window) %>%
    summarise(index = list(index), .groups = "drop"),
  lapply(names(windows), function(w) {
    in_window <- bioassays$window %in% w
    tibble(region = names(extra_groups), window = w,
           index = list(which(in_window & (bioassays$region == "A" |
                                             bioassays$country_name ==
                                             "Ghana")),
                        which(in_window)))
  }))


# 1. the posterior predictive checks -----------------------------------------------------

# one replicate of every bioassay per draw, draws x bioassays
simulate_died <- function(d, seed) {
  set.seed(seed)
  p <- clamp(d$p)
  n_sims <- nrow(p)
  phi <- 1 / d$rho[, d$bioassays$type_id, drop = FALSE] - 1
  q <- rbeta(length(p), p * phi, (1 - p) * phi)
  simulated <- rbinom(length(p), rep(n, each = n_sims), q)
  stopifnot(!anyNA(simulated))
  matrix(simulated, n_sims)
}

# the statistics of died (replicates x bioassays) out of n_group each
tolerance <- 1e-9
statistics <- function(died_group, n_group) {
  f <- sweep(died_group, 2, n_group, "/")
  cbind(pooled = rowSums(died_group) / sum(n_group),
        mean = rowMeans(f),
        median = matrixStats::rowMedians(f),
        le05 = rowMeans(f <= 0.05 + tolerance),
        gt30 = rowMeans(f > 0.30 + tolerance),
        ge98 = rowMeans(f >= 0.98 - tolerance))
}

# per statistic of one group (its bioassays `index`): observed, expected (from
# the draws' predictions p, when given), predictive quantiles of the
# statistics of `replicated` (replicates x the group's bioassays) and the
# two-sided mid-p tail probability
check_group <- function(index, replicated, p = NULL) {
  observed <- statistics(matrix(died[index], 1), n[index])[1, ]
  replicated <- statistics(replicated, n[index])
  expected <- if (!is.null(p)) {
    c(pooled = mean(p[, index, drop = FALSE] %*% n[index]) / sum(n[index]),
      mean = mean(p[, index, drop = FALSE]))
  } else {
    c(pooled = NA_real_, mean = NA_real_)
  }
  below <- colMeans(sweep(replicated, 2, observed, "<"))
  equal <- colMeans(sweep(replicated, 2, observed, "=="))
  u <- below + equal / 2
  tibble(statistic = colnames(replicated),
         bioassays = length(index),
         tested = sum(n[index]),
         observed = observed,
         expected = expected[colnames(replicated)],
         median = apply(replicated, 2, median),
         lower90 = apply(replicated, 2, quantile, 0.05),
         upper90 = apply(replicated, 2, quantile, 0.95),
         lower50 = apply(replicated, 2, quantile, 0.25),
         upper50 = apply(replicated, 2, quantile, 0.75),
         tail = pmin(1, 2 * pmin(u, 1 - u)))
}

# replicates of the group's bioassays `index` (replicates x bioassays) with
# the draws' predictions shifted to the observed level: per draw, the shift c
# of their logits that makes the expected pooled mortality the observed one,
# sum n_i plogis(logit p_i + c) = sum died_i. What is left is the shape of the
# predictive distribution at the right mean
matched_replicates <- function(d, index) {
  logit_p <- qlogis(clamp(d$p[, index, drop = FALSE]))
  target <- min(max(sum(died[index]), 0.5), sum(n[index]) - 0.5)
  shift <- vapply(seq_len(nrow(logit_p)), function(s) {
    uniroot(function(c) sum(n[index] * plogis(logit_p[s, ] + c)) - target,
            c(-30, 30), tol = 1e-8)$root
  }, numeric(1))
  # shift recycles down the rows (draws)
  p <- clamp(plogis(logit_p + shift))
  phi <- 1 / d$rho[, d$bioassays$type_id[index], drop = FALSE] - 1
  q <- rbeta(length(p), p * phi, (1 - p) * phi)
  matrix(rbinom(length(p), rep(n[index], each = nrow(p)), q), nrow(p))
}

with_groups <- function(groups, columns, checks) {
  bind_cols(groups[rep(seq_len(nrow(groups)), each = length(statistic_names)),
                   columns],
            bind_rows(checks))
}
matched_groups <- filter(window_groups, region != "Africa")
ppc <- lapply(setNames(nm = labels), function(label) {
  d <- draws[[label]]
  simulated <- simulate_died(d, seed)
  as_fitted <- function(groups) {
    lapply(groups$index, function(index) {
      check_group(index, simulated[, index, drop = FALSE], d$p)
    })
  }
  set.seed(seed + 1)
  matched <- lapply(matched_groups$index, function(index) {
    check_group(index, matched_replicates(d, index))
  })
  list(years = with_groups(year_groups, c("region", "year"),
                           as_fitted(year_groups)) %>%
         mutate(model = label, .before = 1),
       windows = bind_rows(
         with_groups(window_groups, c("region", "window"),
                     as_fitted(window_groups)) %>%
           mutate(check = "as fitted"),
         with_groups(matched_groups, c("region", "window"), matched) %>%
           mutate(check = "mean-matched")) %>%
         mutate(model = label, .before = 1))
})
outside <- function(x) x$observed < x$lower90 | x$observed > x$upper90
ppc_years <- bind_rows(lapply(ppc, `[[`, "years")) %>%
  mutate(outside90 = outside(.))
ppc_windows <- bind_rows(lapply(ppc, `[[`, "windows")) %>%
  relocate(check, .after = model) %>%
  mutate(outside90 = outside(.))
write.csv(ppc_years, file.path(output_dir, "ppc_years.csv"), row.names = FALSE)
write.csv(ppc_windows, file.path(output_dir, "ppc_windows.csv"),
          row.names = FALSE)

# printed: per window and region, each statistic observed and each model's
# median [90%] (* outside)
printable <- function(rows) {
  rows %>%
    mutate(cell = sprintf("%.2f [%.2f,%.2f]%s", median, lower90, upper90,
                          ifelse(outside90, "*", " "))) %>%
    select(window, region, bioassays, statistic, observed, model, cell) %>%
    mutate(observed = sprintf("%.2f", observed)) %>%
    pivot_wider(names_from = model, values_from = cell) %>%
    arrange(window, region, match(statistic, names(statistic_names))) %>%
    as.data.frame()
}
cat("\nper window: observed, and each model's predictive median [90%]",
    "(* observed outside)\n")
print(printable(filter(ppc_windows, check == "as fitted")), right = FALSE)
cat("\nthe same, mean-matched: each draw's predictions shifted on the logit",
    "scale to the observed pooled mortality\n")
print(printable(filter(ppc_windows, check == "mean-matched",
                       statistic != "pooled")), right = FALSE)

cat("\nyears 2005-2024 with the observed pooled mortality outside the 90%",
    "predictive interval, per model and region\n")
year_failures <- ppc_years %>%
  filter(statistic == "pooled", year %in% 2005:2024) %>%
  group_by(model, region) %>%
  summarise(years = n(), outside = sum(outside90),
            which = paste(sprintf("%d (%+.0f)", year[outside90],
                                  100 * (observed - median)[outside90]),
                          collapse = ", "),
            .groups = "drop")
print(as.data.frame(year_failures), right = FALSE)


# the pooled figure ---------------------------------------------------------------------

years <- 2005:2024
offsets <- setNames(seq(-0.27, 0.27, length.out = length(labels)), labels)
pooled_years <- ppc_years %>%
  filter(statistic == "pooled", year %in% years) %>%
  mutate(model = factor(model, labels),
         x = year + offsets[as.character(model)])
observed_years <- pooled_years %>% distinct(region, year, observed, tested)
pooled_current <- ppc_windows %>%
  filter(check == "as fitted", statistic == "pooled", window == "2019-2024")

model_key <- function(text) {
  keys <- sprintf("<span style='color:%s'>**%s** %s</span>",
                  fit_colours[labels], labels, text)
  # three to a line, so that more than four models fit the panel
  lines <- split(keys, ceiling(seq_along(keys) / 3))
  paste(vapply(lines, paste, "", collapse = ", "), collapse = "<br>")
}
region_title <- function(r) {
  counts <- region_counts[region_counts$region == r, ]
  title <- sub("^([A-Z]): ", "**\\1**: ", counts$name)
  if (!is.na(counts$attached) && nzchar(counts$attached)) {
    title <- sprintf("%s (+ %s*)", title, counts$attached)
  }
  paste(strwrap(title, 58), collapse = "<br>")
}
pooled_panel <- function(r) {
  current <- pooled_current %>% filter(region == r) %>%
    arrange(match(model, labels))
  misses <- pooled_years %>% filter(region == r) %>%
    group_by(model) %>%
    summarise(text = sprintf("%d/%d", sum(outside90), n()),
              .groups = "drop") %>%
    arrange(match(model, labels))
  subtitle <- sprintf(
    "2019-2024: observed %.0f%%, predicted median [90%%]:<br>%s<br>years outside the 90%%: %s",
    100 * current$observed[1],
    model_key(sprintf("%.0f%% [%.0f-%.0f]", 100 * current$median,
                      100 * current$lower90, 100 * current$upper90)),
    model_key(misses$text))
  ggplot(mapping = aes(x = year)) +
    annotate("rect", xmin = c(2009.5, 2018.5), xmax = c(2015.5, 2024.5),
             ymin = -Inf, ymax = Inf, fill = grey(0.94)) +
    geom_linerange(aes(x, ymin = lower90, ymax = upper90, colour = model),
                   data = filter(pooled_years, region == r),
                   linewidth = 0.45) +
    geom_linerange(aes(x, ymin = lower50, ymax = upper50, colour = model),
                   data = filter(pooled_years, region == r),
                   linewidth = 1.6) +
    geom_point(aes(y = observed), data = filter(observed_years, region == r),
               shape = 21, fill = "black", colour = "white", size = 2.2,
               stroke = 0.4) +
    scale_colour_manual(values = fit_colours,
                        name = "predictive\n50% and 90%") +
    scale_x_continuous(limits = range(years) + c(-0.5, 0.5),
                       breaks = seq(2005, 2025, by = 5)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                       breaks = seq(0, 1, by = 0.25)) +
    labs(x = NULL, y = "pooled LLIN-pyrethroid mortality",
         title = region_title(r), subtitle = subtitle) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_markdown(size = 12, lineheight = 1.1),
          plot.subtitle = element_markdown(size = 9, lineheight = 1.2),
          panel.grid.minor = element_blank())
}
pooled_figure <- wrap_plots(lapply(region_letters, pooled_panel), ncol = 3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Posterior predictive check of the pooled LLIN-pyrethroid mortality per region and year",
    caption = paste0(
      "Black points: observed mortality pooled over the region's ",
      "alpha-cypermethrin, deltamethrin and permethrin bioassays in the ",
      "year (died / tested). Bars: each model's posterior predictive 50% ",
      "(thick) and 90% (thin) intervals of the same\nstatistic, from one ",
      "beta-binomial replicate of every bioassay (its own sample size, ",
      "predicted mortality and the type's rho) per posterior draw, about ",
      n_draws, " draws (V5: chains ",
      toString(draws[["V5"]]$chains), "). Subtitles: the 2019-2024 window, ",
      "predictive median [90%];\nthe years 2005-2024 with the observed ",
      "outside the 90% interval. Shaded: the windows 2010-2015 and ",
      "2019-2024. Regions: R/diagnostic_regions.R (*: no 2019-2024 ",
      "bioassays)."),
    theme = theme(plot.title = element_text(size = 15),
                  plot.caption = element_text(hjust = 0, size = 10)))
ggsave(file.path(figure_dir, "region_ppc_pooled.png"), pooled_figure,
       width = 17, height = 15, dpi = dpi, bg = "white")


# the shape figure ----------------------------------------------------------------------

shape_statistics <- c(pooled = "pooled mortality", median = "median",
                      le05 = "share <= 5%", gt30 = "share > 30%",
                      ge98 = "share >= 98%")
group_levels <- c(region_letters, "A+Ghana")
shape_figure <- function(which_check, statistics_shown, title, note) {
  shape <- ppc_windows %>%
    filter(check == which_check, statistic %in% statistics_shown,
           region %in% group_levels) %>%
    mutate(statistic = factor(shape_statistics[statistic],
                              shape_statistics[statistics_shown]),
           group = factor(region, rev(group_levels)),
           y = as.numeric(group) - offsets[model] * 1.2,
           model = factor(model, labels))
  ggplot(shape) +
    geom_linerange(aes(xmin = lower90, xmax = upper90, y = y, colour = model),
                   linewidth = 0.45) +
    geom_linerange(aes(xmin = lower50, xmax = upper50, y = y, colour = model),
                   linewidth = 1.6) +
    geom_point(aes(observed, as.numeric(group)),
               data = distinct(shape, statistic, window, group, observed),
               shape = 21, fill = "black", colour = "white", size = 2.4,
               stroke = 0.4) +
    facet_grid(window ~ statistic, scales = "free_x") +
    scale_y_continuous(breaks = seq_along(levels(shape$group)),
                       labels = levels(shape$group)) +
    scale_x_continuous(labels = scales::percent) +
    scale_colour_manual(values = fit_colours,
                        name = "predictive\n50% and 90%") +
    labs(x = NULL, y = NULL, title = title,
         caption = paste0(
           "Black points: the observed statistic of the bioassays of the ",
           "region and window (pooled: died / tested; the others of each ",
           "bioassay's died / n). Bars: each model's posterior predictive ",
           "50% (thick) and 90% (thin)\nintervals, from one beta-binomial ",
           "replicate of every bioassay per posterior draw (about ",
           n_draws, " draws; V5: chains ",
           toString(draws[["V5"]]$chains), ")", note, ". A+Ghana: region A ",
           "with Ghana (in region C). Regions: R/diagnostic_regions.R.")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(colour = grey(0.9)),
          panel.spacing.x = unit(1.2, "lines"),
          strip.text = element_text(size = 12, hjust = 0),
          plot.title = element_text(size = 15),
          plot.caption = element_text(hjust = 0, size = 10))
}
ggsave(file.path(figure_dir, "region_ppc_shape.png"),
       shape_figure("as fitted", c("pooled", "median", "le05", "gt30"),
                    "Posterior predictive check of the shape of the LLIN-pyrethroid bioassays per region and window",
                    ""),
       width = 17, height = 13, dpi = dpi, bg = "white")
ggsave(file.path(figure_dir, "region_ppc_shape_matched.png"),
       shape_figure("mean-matched", c("median", "le05", "gt30", "ge98"),
                    "The shape at the right mean: posterior predictive check with each draw's predictions shifted to the observed pooled mortality",
                    paste0(",\nwith each draw's predicted mortality at the ",
                           "region-window's bioassays shifted on the logit ",
                           "scale so that its expected pooled mortality is ",
                           "the observed")),
       width = 17, height = 13, dpi = dpi, bg = "white")


# 2. the error model ---------------------------------------------------------------------

# beta-binomial log likelihood with mean mu and intra-class correlation rho,
# as betabinomial_p_rho() (R/functions.R): phi = 1 / rho - 1
bb_loglik <- function(died, n, mu, rho) {
  mu <- clamp(mu)
  phi <- 1 / rho - 1
  sum(extraDistr::dbbinom(died, n, alpha = mu * phi, beta = (1 - mu) * phi,
                          log = TRUE))
}
fixed_rho <- rho_means %>% filter(model == fixed_rho_fit) %>%
  select(insecticide_type, rho) %>% tibble::deframe()
# V4_class's posterior mean prediction at each bioassay, for the offset fits
offset_logit <- qlogis(clamp(colMeans(draws[[fixed_rho_fit]]$p)))

# maximum likelihood fits of logit mu_i = theta_k(i) + offset_i, one theta per
# type k present, with rho fixed per type (rho_fixed, by type) or freed by one
# shift delta of the types' logit rho (so that one type's rho is free, and the
# fixed rho is nested in the free); returns the fitted mean at each bioassay,
# the free rho (the mean over the bioassays of the shifted rho) and its 90%
# Wald interval, and the log likelihoods
fit_constant <- function(died, n, type, offset, rho_fixed) {
  types <- sort(unique(type))
  k <- match(type, types)
  logit_rho <- qlogis(rho_fixed[type])
  # with rho fixed, the types are independent
  theta_fixed <- vapply(seq_along(types), function(j) {
    in_type <- k == j
    optimize(function(theta) {
      -bb_loglik(died[in_type], n[in_type], plogis(theta + offset[in_type]),
                 rho_fixed[[types[j]]])
    }, c(-15, 15))$minimum
  }, numeric(1))
  ll_fixed <- bb_loglik(died, n, plogis(theta_fixed[k] + offset),
                        plogis(logit_rho))
  negative <- function(par) {
    -bb_loglik(died, n, plogis(par[k] + offset),
               plogis(logit_rho + par[length(par)]))
  }
  free <- optim(c(theta_fixed, 0), negative, method = "L-BFGS-B",
                lower = c(rep(-15, length(types)), -9),
                upper = c(rep(15, length(types)), 9))
  hessian <- optimHess(free$par, negative)
  se <- tryCatch(sqrt(diag(solve(hessian)))[length(free$par)],
                 error = function(e) NA_real_)
  delta <- free$par[length(free$par)]
  shifted <- function(by) mean(plogis(logit_rho + by))
  list(mu_fixed = plogis(theta_fixed[k] + offset),
       mu_free = plogis(free$par[k] + offset),
       rho_free = shifted(delta),
       rho_lower = shifted(delta - qnorm(0.95) * se),
       rho_upper = shifted(delta + qnorm(0.95) * se),
       ll_fixed = ll_fixed, ll_free = -free$value,
       converged = free$convergence == 0)
}

error_groups <- window_groups %>% filter(region %in% group_levels)
min_bioassays <- 5
error_rows <- function(index, type_name = NULL, offset = 0) {
  d <- died[index]
  m <- n[index]
  type <- bioassays$insecticide_type[index]
  if (length(index) < min_bioassays) return(NULL)
  fitted <- fit_constant(d, m, type, rep_len(offset, length(index)),
                         fixed_rho)
  tibble(bioassays = length(index), tested = sum(m),
         observed_pooled = sum(d) / sum(m), observed_mean = mean(d / m),
         rho_fixed = mean(fixed_rho[type]),
         mean_fixed = mean(fitted$mu_fixed), mean_free = mean(fitted$mu_free),
         rho_free = fitted$rho_free, rho_free_lower = fitted$rho_lower,
         rho_free_upper = fitted$rho_upper,
         lr = 2 * (fitted$ll_free - fitted$ll_fixed),
         ll_fixed = fitted$ll_fixed, ll_free = fitted$ll_free,
         converged = fitted$converged)
}
error_model <- bind_rows(lapply(seq_len(nrow(error_groups)), function(g) {
  index <- error_groups$index[[g]]
  bind_rows(lapply(llin_pyrethroids, function(type_name) {
    rows <- error_rows(index[bioassays$insecticide_type[index] == type_name],
                       type_name)
    if (!is.null(rows)) {
      mutate(rows, region = error_groups$region[g],
             window = error_groups$window[g],
             insecticide_type = type_name, .before = 1)
    }
  }))
}))
error_model_pooled <- bind_rows(lapply(seq_len(nrow(error_groups)), function(g) {
  index <- error_groups$index[[g]]
  constant <- error_rows(index)
  offset <- error_rows(index, offset = offset_logit[index])
  if (is.null(constant)) return(NULL)
  bind_cols(tibble(region = error_groups$region[g],
                   window = error_groups$window[g]),
            constant %>% select(-ll_fixed, -ll_free),
            offset %>% select(offset_mean_fixed = mean_fixed,
                              offset_mean_free = mean_free,
                              offset_rho_free = rho_free,
                              offset_lr = lr),
            offset_model_mean = mean(plogis(offset_logit[index])))
}))

# A mean-dependent rho: rho_k(mu) = plogis(a_k + b log(4 mu (1 - mu))), a_k
# per type, b common, fitted to the region-window-type cells of the regions
# (each with its own constant mean, profiled out) with at least
# min_bioassays bioassays. b > 0 makes rho smaller towards 0 and 1; a
# constant overdispersion on the logit scale (a logit-normal mean) gives
# about rho = sigma^2 mu (1 - mu), i.e. b near 1 where rho is small
rho_dependent <- function(mu, type, par) {
  plogis(par[match(type, llin_pyrethroids)] +
           par[length(llin_pyrethroids) + 1] * log(4 * mu * (1 - mu)))
}
dependent_mean <- function(died, n, type, par) {
  fit <- optimize(function(theta) {
    mu <- clamp(plogis(theta))
    -bb_loglik(died, n, mu, rho_dependent(mu, type, par))
  }, c(-15, 15))
  list(mu = plogis(fit$minimum), ll = -fit$objective)
}
cell_index <- lapply(seq_len(nrow(error_model)), function(r) {
  index <- error_groups$index[[which(
    error_groups$region == error_model$region[r] &
      error_groups$window == error_model$window[r])]]
  index[bioassays$insecticide_type[index] == error_model$insecticide_type[r]]
})
in_regions <- which(error_model$region %in% region_letters)
dependent_fit <- optim(
  c(rep(qlogis(0.33), length(llin_pyrethroids)), 0),
  function(par) {
    -sum(vapply(cell_index[in_regions], function(index) {
      dependent_mean(died[index], n[index], bioassays$insecticide_type[index],
                     par)$ll
    }, numeric(1)))
  }, method = "Nelder-Mead", control = list(maxit = 2000, reltol = 1e-10))
dependent_par <- setNames(dependent_fit$par,
                          c(paste0("a_", llin_pyrethroids), "b"))
cells_dependent <- lapply(cell_index, function(index) {
  dependent_mean(died[index], n[index], bioassays$insecticide_type[index],
                 dependent_par)
})
error_model <- error_model %>%
  mutate(mean_dependent = vapply(cells_dependent, `[[`, numeric(1), "mu"),
         rho_dependent = rho_dependent(mean_dependent, insecticide_type,
                                       dependent_par),
         ll_dependent = vapply(cells_dependent, `[[`, numeric(1), "ll"),
         .after = rho_free_upper)
# with the types together: each type's mean under the mean-dependent rho,
# averaged over the bioassays
error_model_pooled$mean_dependent <- vapply(
  seq_len(nrow(error_model_pooled)), function(r) {
    index <- error_groups$index[[which(
      error_groups$region == error_model_pooled$region[r] &
        error_groups$window == error_model_pooled$window[r])]]
    type <- bioassays$insecticide_type[index]
    mu <- numeric(length(index))
    for (k in unique(type)) {
      at <- type == k
      mu[at] <- dependent_mean(died[index][at], n[index][at], type[at],
                               dependent_par)$mu
    }
    mean(mu)
  }, numeric(1))
error_model_pooled <- relocate(error_model_pooled, mean_dependent,
                               .after = mean_free)

write.csv(error_model, file.path(output_dir, "error_model.csv"),
          row.names = FALSE)
write.csv(error_model_pooled, file.path(output_dir, "error_model_pooled.csv"),
          row.names = FALSE)
cat(sprintf("\nerror model: rho fixed at %s's posterior mean (%s)\n",
            fixed_rho_fit,
            paste(sprintf("%s %.3f", names(fixed_rho), fixed_rho),
                  collapse = ", ")))
cat("\nthe types together (a mean per type; free: one shift of the types'",
    "logit rho)\n")
print(as.data.frame(error_model_pooled %>% select(-converged)), digits = 3)
cat("\nper type\n")
print(as.data.frame(error_model %>% select(-converged)), digits = 3)
if (!all(error_model$converged, error_model_pooled$converged)) {
  report("some free-rho fits did not converge")
}
cat("\nmean-dependent rho, rho_k(mu) = plogis(a_k + b log(4 mu (1 - mu))):\n")
print(round(dependent_par, 3))
grid <- expand.grid(mu = c(0.05, 0.15, 0.3, 0.5, 0.7, 0.9, 0.97),
                    type = llin_pyrethroids, stringsAsFactors = FALSE)
grid$rho <- rho_dependent(grid$mu, grid$type, dependent_par)
print(tidyr::pivot_wider(grid, names_from = mu, values_from = rho),
      digits = 3)
cells <- error_model[in_regions, ]
cat(sprintf(paste0("log likelihood over the %d region-window-type cells: ",
                   "fixed rho %.1f; mean-dependent rho (4 parameters) %.1f; ",
                   "free rho per cell (%d parameters) %.1f\n"),
            nrow(cells), sum(cells$ll_fixed), sum(cells$ll_dependent),
            nrow(cells), sum(cells$ll_free)))
cat("mean absolute gap of the fitted mean to the observed unweighted mean",
    "(percentage points), the cells with at least 20 bioassays:\n")
print(cells %>% filter(bioassays >= 20) %>%
        mutate(band = cut(observed_mean, c(0, 0.3, 0.8, 1),
                          include.lowest = TRUE)) %>%
        group_by(band) %>%
        summarise(cells = n(),
                  fixed = mean(100 * (mean_fixed - observed_mean)),
                  free = mean(100 * (mean_free - observed_mean)),
                  dependent = mean(100 * (mean_dependent - observed_mean)),
                  abs_fixed = mean(abs(100 * (mean_fixed - observed_mean))),
                  abs_free = mean(abs(100 * (mean_free - observed_mean))),
                  abs_dependent = mean(abs(100 * (mean_dependent -
                                                    observed_mean))),
                  .groups = "drop") %>%
        as.data.frame(), digits = 3)

# does the free rho follow the mean? A weighted fit of its logit on the
# observed mean and its square (weights from the Wald intervals), per type
free_rho_trend <- error_model %>%
  filter(bioassays >= 20, region %in% region_letters,
         is.finite(rho_free_lower)) %>%
  mutate(logit_rho = qlogis(rho_free),
         se = (qlogis(rho_free_upper) - qlogis(rho_free_lower)) /
           (2 * qnorm(0.95)))
trend_fit <- lm(logit_rho ~ observed_mean + I(observed_mean ^ 2),
                data = free_rho_trend, weights = 1 / se ^ 2)
cat("\nlogit free rho on the observed mean (types separately, region-windows",
    "with at least 20 bioassays; weighted by the inverse Wald variance)\n")
print(summary(trend_fit)$coefficients, digits = 3)


# the error-model figure ------------------------------------------------------------------

window_colours <- c(`2010-2015` = "#56B4E9", `2016-2018` = "#E69F00",
                    `2019-2024` = "#D55E00")
type_colours <- c(`Alpha-cypermethrin` = "#009E73", Deltamethrin = "#0072B2",
                  Permethrin = "#CC79A7")
gaps <- error_model_pooled %>%
  mutate(group = sprintf("%s %s", region, sub("20(..)-20(..)", "\\1-\\2",
                                                window))) %>%
  select(group, region, window, bioassays, observed_mean, mean_fixed,
         mean_free, mean_dependent, offset_mean_fixed, offset_mean_free) %>%
  pivot_longer(c(mean_fixed, mean_free, mean_dependent, offset_mean_fixed,
                 offset_mean_free),
               names_to = "fit", values_to = "fitted") %>%
  mutate(rho = case_when(grepl("fixed", fit) ~
                           sprintf("fixed (%s)", fixed_rho_fit),
                         grepl("free", fit) ~ "free",
                         TRUE ~ "mean-dependent"),
         mean_model = ifelse(grepl("offset", fit),
                             sprintf("offset: %s's prediction", fixed_rho_fit),
                             "constant mean per type"),
         gap = 100 * (fitted - observed_mean))
gap_panel <- ggplot(gaps, aes(observed_mean, gap)) +
  geom_hline(yintercept = 0, colour = grey(0.4)) +
  geom_line(aes(group = group), colour = grey(0.7), linewidth = 0.3) +
  geom_point(aes(colour = window, shape = rho, size = bioassays),
             stroke = 0.8) +
  ggrepel::geom_text_repel(aes(label = group),
                           data = filter(gaps, grepl("fixed", rho)),
                           size = 2.8, colour = grey(0.3), min.segment.length = 0.2,
                           max.overlaps = 30, seed = 1) +
  facet_wrap(~ mean_model, nrow = 1) +
  scale_colour_manual(values = window_colours) +
  scale_shape_manual(values = c(16, 1, 2)) +
  scale_size_area(max_size = 5, breaks = c(50, 200, 800)) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "observed unweighted mean of died / n",
       y = "fitted mean minus observed\n(percentage points)",
       shape = "rho", size = "bioassays",
       title = "a. Maximum likelihood mean of the region-window's LLIN-pyrethroid bioassays, with rho fixed, free and mean-dependent") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 12, hjust = 0))
rho_panel <- ggplot(filter(error_model, bioassays >= 20,
                           region %in% region_letters),
                    aes(observed_mean, rho_free)) +
  geom_hline(aes(yintercept = rho, colour = insecticide_type),
             data = filter(rho_means, model == fixed_rho_fit),
             linetype = "dashed") +
  geom_linerange(aes(ymin = rho_free_lower, ymax = rho_free_upper,
                     colour = insecticide_type), linewidth = 0.4) +
  geom_line(aes(mu, rho, colour = insecticide_type),
            data = expand.grid(mu = seq(0.01, 0.99, by = 0.01),
                               type = llin_pyrethroids,
                               stringsAsFactors = FALSE) %>%
              mutate(rho = rho_dependent(mu, type, dependent_par),
                     insecticide_type = type),
            linewidth = 0.8) +
  geom_point(aes(colour = insecticide_type, size = bioassays)) +
  geom_point(aes(observed_mean, rho_free, size = bioassays),
             data = filter(error_model_pooled, bioassays >= 20,
                           region %in% region_letters), shape = 4,
             colour = "black") +
  scale_colour_manual(values = type_colours, name = NULL) +
  scale_size_area(max_size = 5, breaks = c(50, 200, 800)) +
  scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "observed unweighted mean of died / n",
       y = "free rho (constant mean)", size = "bioassays",
       title = "b. Free rho against the mean mortality, per region, window and type (crosses: the types together)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())
error_figure <- (gap_panel / rho_panel) +
  plot_annotation(
    title = "Does the fixed overdispersion pull the fitted mean away from the data?",
    caption = paste0(
      "Beta-binomial maximum likelihood fits to the LLIN-pyrethroid ",
      "bioassays of each region and window (and region A with Ghana), a ",
      "mean per insecticide type, with rho fixed at ", fixed_rho_fit,
      "'s posterior mean per type\n(dashed in b), free (one shift of the ",
      "types' logit rho), or mean-dependent (rho_k(mu) = plogis(a_k + b ",
      "log(4 mu (1 - mu))), fitted across the cells; lines in b). a: the ",
      "mean of the fitted means over the\nbioassays minus the observed ",
      "unweighted mean; left, a constant mean per type; right, ",
      fixed_rho_fit, "'s posterior mean prediction as an offset (a shift ",
      "of its logit per type). b: the free rho of a constant\nmean per ",
      "type, with its 90% Wald interval, for the region-windows with at ",
      "least 20 bioassays of the type; a constant mean puts the spread of ",
      "the true mortality between sites and years into rho."),
    theme = theme(plot.title = element_text(size = 15),
                  plot.caption = element_text(hjust = 0, size = 10)))
ggsave(file.path(figure_dir, "region_error_model.png"), error_figure,
       width = 15, height = 13, dpi = dpi, bg = "white")
report("written %s, %s, %s, %s, and the tables in %s",
       file.path(figure_dir, "region_ppc_pooled.png"),
       file.path(figure_dir, "region_ppc_shape.png"),
       file.path(figure_dir, "region_ppc_shape_matched.png"),
       file.path(figure_dir, "region_error_model.png"), output_dir)
