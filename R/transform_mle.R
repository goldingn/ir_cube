# Screen of fixed covariate transforms by maximum likelihood (#37): for net
# use, the encounter transform 1 - exp(-n log(2) / n_half) against linear net
# use; for population density d, the encounter transform 1 - exp(-d log(2) /
# d_half) for several d_half; each on the within-pixel pairs of a fit with no
# mortality floor (R/change_residual_screen.R).
#
#   Rscript R/transform_mle.R <fit RData> <label> [gap]
#
# reads outputs/hinge_screen/<label>[_gap<g>]_change_pairs.csv (default gap
# 5). For a pair of sampled years y1 < y2 at a pixel, for one insecticide
# type, the model gives, with mortality the susceptible fraction q,
#
#   logit q(y2) = logit q(y1)
#                 - sum_{t = y1 + 1}^{y2} log(1 + o_t + exp(b) g_t T(x_t))
#                 - (y2 - y1) kappa
#
# T the candidate transform of covariate x; o_t the selection of the fit's
# other columns at the type's posterior mean exp(beta), held fixed; kappa the
# type's posterior mean reversion; g_t the column's time trend (g_dom(t) for
# population, 1 for net use). exp(b) is fitted per type; logit q(y1) is a
# free parameter of each pair, profiled out, so that only the change within
# the pair informs the fit. Deaths are binomial on effective counts (divided
# by 1 + (n - 1) rho, the fit's per-type rho). Every candidate has one
# parameter per type, so their log likelihoods compare directly.
#
# Writes outputs/hinge_screen/<label>_transform_loglik.csv (each candidate's
# log likelihood and coefficient by type) and
# figures/transform_mle_<label>_<covariate>.png (each type's fitted selection
# per year under the linear model and the best candidates, over the pairs'
# partial residuals).

arguments <- commandArgs(trailingOnly = TRUE)
fit_file <- arguments[1]
label <- arguments[2]
gap <- if (length(arguments) >= 3) as.integer(arguments[3]) else 5L
output_dir <- "outputs/hinge_screen"
n_cores <- min(8L, parallel::detectCores())
candidates <- list(
  nets = c(linear = Inf, n_half_0.01 = 0.01, n_half_0.02 = 0.02,
           n_half_0.05 = 0.05, n_half_0.1 = 0.1,
           n_half_0.2 = 0.2, n_half_0.3 = 0.3, n_half_0.5 = 0.5,
           n_half_1 = 1),
  pop = c(d_half_50 = 50, d_half_100 = 100, d_half_200 = 200,
          d_half_270 = 270, d_half_500 = 500, d_half_1000 = 1000,
          d_half_3000 = 3000))

source("R/packages.R")
source("R/functions.R")
source("R/model_covariates.R")
source("R/dynamical_model.R")
source("R/dynamical_predictions.R")

f <- new.env()
load(fit_file, envir = f)
stopifnot(isFALSE(f$model_options$mortality_floor))
fold <- list(draws = f$draws, options = f$model_options,
             x_cells_init = f$x_cells_init)
n_total <- sum(vapply(fold$draws, nrow, integer(1)))
par <- dynamical_parameter_draws(fold, f$classes_index, f$types,
                                 draw_index = round(seq(1, n_total,
                                                        length.out = 400)))
effect <- apply(par$effect_type, c(2, 3), mean)
dimnames(effect) <- list(colnames(f$x_cell_years), f$types)
kappa <- if (is.null(par$kappa_type)) rep(0, length(f$types)) else
  colMeans(par$kappa_type)
rho <- colMeans(par$rho_types)
names(kappa) <- names(rho) <- f$types
x_cell_years <- f$x_cell_years
n_years <- max(f$cell_years_index$year_id)
x_row <- matrix(NA_integer_, max(f$cell_years_index$cell_id), n_years)
x_row[cbind(f$cell_years_index$cell_id, f$cell_years_index$year_id)] <-
  seq_len(nrow(f$cell_years_index))
baseline_year <- f$baseline_year
pop_column <- grep("^pop_[a-z]+:g_dom$", colnames(x_cell_years), value = TRUE)
# population density (people per km2), cells x years, the input of the
# population transform
density <- pop_density_matrix(f$unique_cells, baseline_year,
                              baseline_year + n_years - 1)
rm(f, par, fold)
invisible(gc())
g_dom <- function(year) pmax((year - 1995) / 30, 0)

# the encounter transform of x at half-saturation x_half (Inf: x itself)
encounter <- function(x, x_half) {
  if (is.infinite(x_half)) x else -expm1(-x * log(2) / x_half)
}

stem <- if (gap > 1) sprintf("%s_gap%d", label, gap) else label
pairs <- read.csv(file.path(output_dir, sprintf("%s_change_pairs.csv", stem)))

# One type's pairs: effective counts of both years, and pairs x years
# matrices (0 beyond each span) of the raw covariate, its time trend and the
# other columns' selection
type_data <- function(type, covariate) {
  d <- pairs[pairs$insecticide_type == type, ]
  span <- d$year_id_2 - d$year_id_1
  width <- max(span)
  years <- outer(d$year_id_1, seq_len(width), "+")
  mask <- outer(span, seq_len(width), ">=")
  years[!mask] <- d$year_id_1[row(years)][!mask] + 1
  rows <- matrix(x_row[cbind(rep(d$cell_id, width), as.vector(years))],
                 nrow(d))
  calendar <- baseline_year + years - 1
  if (covariate == "nets") {
    column <- "nets"
    x <- matrix(x_cell_years[rows, "nets"], nrow(d))
    g <- 1 + 0 * calendar
  } else {
    column <- pop_column
    x <- matrix(density[cbind(rep(d$cell_id, width), as.vector(years))],
                nrow(d))
    g <- g_dom(calendar)
  }
  others <- setdiff(colnames(x_cell_years), column)
  o <- matrix(x_cell_years[rows, others, drop = FALSE] %*%
                effect[others, type], nrow(d))
  stopifnot(!anyNA(x[mask]), !anyNA(o[mask]))
  zero_outside <- function(v) ifelse(mask, v, 0)
  inflation <- function(n) 1 + (n - 1) * rho[[type]]
  list(x = zero_outside(x), g = zero_outside(g), o = zero_outside(o),
       mask = mask, span = span, pair = d$pair,
       y1 = d$died_1 / inflation(d$n_1), m1 = d$n_1 / inflation(d$n_1),
       y2 = d$died_2 / inflation(d$n_2), m2 = d$n_2 / inflation(d$n_2),
       kappa = kappa[[type]])
}

binomial_loglik <- function(y, m, p) {
  y * log(p) + (m - y) * log1p(-p)
}

# The log likelihood of a type's pairs given each year's selection term
# (pairs x years), with each pair's starting logit profiled out by Newton's
# method, vectorised over pairs (the log likelihood is concave in it)
profile_loglik <- function(d, selection) {
  s <- rowSums(log1p(d$o + d$g * selection) * d$mask) + d$span * d$kappa
  if (any(!is.finite(s))) {
    return(-1e12)
  }
  a <- qlogis(pmin(pmax((d$y1 + d$y2 + 0.5) / (d$m1 + d$m2 + 1), 1e-4),
                   1 - 1e-4))
  for (iteration in 1:30) {
    p1 <- plogis(a)
    p2 <- plogis(a - s)
    gradient <- (d$y1 - d$m1 * p1) + (d$y2 - d$m2 * p2)
    curvature <- d$m1 * p1 * (1 - p1) + d$m2 * p2 * (1 - p2)
    step <- gradient / pmax(curvature, 1e-8)
    a <- pmin(pmax(a + step, -20), 20)
    if (max(abs(step)) < 1e-8) break
  }
  p1 <- pmin(pmax(plogis(a), 1e-12), 1 - 1e-12)
  p2 <- pmin(pmax(plogis(a - s), 1e-12), 1 - 1e-12)
  sum(binomial_loglik(d$y1, d$m1, p1) + binomial_loglik(d$y2, d$m2, p2))
}

# the maximum likelihood over b of a type's pairs under transform half-point
# x_half
fit_candidate <- function(d, x_half) {
  tx <- encounter(d$x, x_half) * d$mask
  o <- optimize(function(b) -profile_loglik(d, exp(b) * tx), c(-15, 5))
  list(loglik = -o$objective, coefficient = exp(o$minimum))
}

types <- sort(unique(pairs$insecticide_type))
classes <- setNames(ifelse(types == "DDT", "DDT",
                           pairs$insecticide_class[match(types,
                                                         pairs$insecticide_type)]),
                    types)
results <- list()
for (covariate in names(candidates)) {
  data <- setNames(lapply(types, type_data, covariate = covariate), types)
  for (name in names(candidates[[covariate]])) {
    x_half <- candidates[[covariate]][[name]]
    fits <- parallel::mclapply(data, fit_candidate, x_half = x_half,
                               mc.cores = n_cores)
    results[[length(results) + 1]] <- data.frame(
      covariate = covariate, candidate = name, x_half = x_half,
      type = types, class = classes[types],
      loglik = vapply(fits, `[[`, numeric(1), "loglik"),
      coefficient = vapply(fits, `[[`, numeric(1), "coefficient"))
  }
  cat(covariate, "done\n")
}
results <- do.call(rbind, results)
rownames(results) <- NULL
write.csv(results, file.path(output_dir, sprintf("%s_transform_loglik.csv",
                                                 stem)), row.names = FALSE)

# log likelihood relative to the fit's own transform, by candidate: total,
# pyrethroids, and each type
reference <- c(nets = "linear", pop = "d_half_270")
summary_table <- results |>
  group_by(covariate, type) |>
  mutate(delta = loglik - loglik[candidate == reference[covariate[1]]]) |>
  ungroup()
wide <- summary_table |>
  group_by(covariate, candidate) |>
  summarise(total = sum(delta),
            pyrethroids = sum(delta[class == "Pyrethroids"]),
            .groups = "drop") |>
  left_join(summary_table |>
              select(covariate, candidate, type, delta) |>
              tidyr::pivot_wider(names_from = type, values_from = delta),
            by = c("covariate", "candidate"))
options(width = 250)
print(as.data.frame(mutate(wide, across(where(is.numeric), ~ round(.x, 1)))),
      row.names = FALSE)

# plots: each type's fitted selection per year (exp(b) T(x), per unit of the
# time trend for population) under the fit's own transform and the two best
# candidates by total log likelihood, over the pairs' partial residuals (the
# observed fall per year minus the fit's predicted fall plus the fit's own
# effect of the covariate, per unit of the span's mean trend)
span <- pairs$year_2 - pairs$year_1
fall <- (-pairs$change_obs + pairs$change_pred) / span
for (covariate in names(candidates)) {
  data <- setNames(lapply(types, type_data, covariate = covariate), types)
  points <- bind_rows(lapply(types, function(type) {
    d <- data[[type]]
    i <- match(d$pair, pairs$pair)
    width <- rowSums(d$mask)
    x_mean <- rowSums(d$x) / width
    g_mean <- rowSums(d$g) / width
    fit_column <- if (covariate == "nets") d$x else
      d$g * encounter(d$x, 270)
    tibble(type = type, x = x_mean,
           y = (fall[i] + effect[if (covariate == "nets") "nets" else
                                   pop_column, type] *
                  rowSums(fit_column) / width) / g_mean)
  })) |>
    filter(is.finite(y))
  best <- wide |>
    filter(covariate == !!covariate) |>
    arrange(desc(total)) |>
    pull(candidate)
  shown <- unique(c(reference[[covariate]], head(best, 2)))
  at <- seq(0, quantile(points$x, 0.99), length.out = 200)
  curves <- results |>
    filter(covariate == !!covariate, candidate %in% shown) |>
    rowwise() |>
    reframe(type = type, candidate = candidate, x = at,
            y = coefficient * encounter(at, x_half))
  labels <- summary_table |>
    filter(covariate == !!covariate, candidate %in% shown) |>
    group_by(type) |>
    summarise(text = paste(sprintf("%s %+.1f", candidate, delta),
                           collapse = "\n"), .groups = "drop")
  limits <- quantile(points$y, c(0.02, 0.98))
  p <- ggplot(points, aes(x, y)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_point(size = 0.4, alpha = 0.25, colour = "grey30") +
    geom_line(data = curves, aes(colour = candidate), linewidth = 0.8) +
    geom_text(data = labels, aes(x = -Inf, y = Inf, label = text),
              hjust = -0.05, vjust = 1.1, size = 2.5, lineheight = 0.9) +
    facet_wrap(~ type, ncol = 3) +
    coord_cartesian(ylim = limits) +
    labs(x = if (covariate == "nets") "net use (mean over the pair's span)" else
           "population density, people per km2 (mean over the pair's span)",
         y = "selection per year (logit mortality)", colour = NULL,
         title = sprintf("Transform screen, %s, %d-year pairs, %s", label,
                         gap, covariate),
         subtitle = paste("Lines: fitted selection per year; text: log",
                          "likelihood relative to the fit's transform;",
                          "points: pairs' partial residuals")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  if (covariate == "pop") {
    p <- p + scale_x_continuous(transform = scales::pseudo_log_trans(sigma = 10),
                                breaks = c(0, 10, 100, 1000, 10000))
  }
  ggsave(sprintf("figures/transform_mle_%s_%s.png", stem, covariate), p,
         width = 10, height = 11, dpi = 130, bg = "white")
}
