# Hinge selection by maximum likelihood (#37): the number and places of hinge
# knots for a selection covariate, on the within-pixel pairs of a fit with no
# mortality floor (R/change_residual_screen.R), with the dynamical model's own
# transform and constraint.
#
#   Rscript R/hinge_mle.R <fit RData> <label> [gap ...]
#
# reads outputs/hinge_screen/<label>[_gap<g>]_change_pairs.csv for each gap
# (default 5).
#
# For a pair of sampled years y1 < y2 at a pixel, for one insecticide type, the
# model (R/dynamical_model.R) gives, with mortality the susceptible fraction q,
#
#   logit q(y2) = logit q(y1)
#                 - sum_{t = y1 + 1}^{y2} log(1 + o_t + g_t f(x_t)) - (y2 - y1) kappa
#
#   f(x) = exp(b) x + sum_j exp(c_j) min(x, k_j)
#
# o_t: the selection of the fit's other columns, sum of x_kt exp(beta_k) at
# the type's posterior mean exp(beta), held fixed; kappa: the type's posterior
# mean reversion; g_t: the column's time trend (g_dom(t) for population, 1 for
# net use). f has positive coefficients, so it is increasing and concave, as
# in the model, and there is no intercept (no selection without the
# covariates). logit q(y1) is a free parameter of each pair, profiled out, so
# that only the change within the pair informs f. The two years' pooled deaths
# are binomial on effective counts (deaths and tested divided by
# 1 + (n - 1) rho, the fit's per-type beta-binomial rho), approximating the
# beta-binomial.
#
# Knots are shared across all insecticides; b and c_j are fitted per type. For
# 0 to 3 hinges the knots maximise the likelihood summed over types, over a
# grid (steps of 0.1 within the 1-99% range of x, distinct), with the
# types fitted in parallel. Each type's score is its AIC, -2 log
# likelihood + 2 (1 + hinges) (the pair parameters are common to every model
# and cancel); a class's is the sum over its types.
#
# Writes outputs/hinge_screen/<label>_mle_knots.csv, <label>_mle_by_group.csv
# (each type's and class's AIC change from the linear model, and the fitted
# coefficients) and <label>_mle_selected.csv.

arguments <- commandArgs(trailingOnly = TRUE)
fit_file <- arguments[1]
label <- arguments[2]
gaps <- if (length(arguments) >= 3) as.integer(arguments[-(1:2)]) else
  5L
output_dir <- "outputs/hinge_screen"
# IR_CUBE_PAIR_FRACTION < 1 fits a random share of each type's pairs, for
# timing; the outputs then carry "_sub<fraction>" in their names
pair_fraction <- as.numeric(Sys.getenv("IR_CUBE_PAIR_FRACTION", "1"))
output_label <- if (pair_fraction < 1) sprintf("%s_sub%g", label,
                                                pair_fraction) else label
# the largest number of hinges, and the knot grid's step, settable with
# IR_CUBE_MAX_HINGES and IR_CUBE_KNOT_STEP; other values put "_h<max>_step<step>"
# in the output names
max_hinges <- as.integer(Sys.getenv("IR_CUBE_MAX_HINGES", "3"))
knot_step <- as.numeric(Sys.getenv("IR_CUBE_KNOT_STEP", "0.1"))
if (max_hinges != 3 || knot_step != 0.1) {
  output_label <- sprintf("%s_h%d_step%g", output_label, max_hinges, knot_step)
}
# cores to fit the types on in parallel
n_cores <- min(8L, parallel::detectCores())

source("R/packages.R")
source("R/functions.R")
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
# posterior means: selection effects exp(beta), columns x types; reversion
# and rho, per type
effect <- apply(par$effect_type, c(2, 3), mean)
dimnames(effect) <- list(colnames(f$x_cell_years), f$types)
kappa <- if (is.null(par$kappa_type)) rep(0, length(f$types)) else
  colMeans(par$kappa_type)
rho <- colMeans(par$rho_types)
names(kappa) <- names(rho) <- f$types
x_cell_years <- f$x_cell_years
x_row <- matrix(NA_integer_, max(f$cell_years_index$cell_id),
                max(f$cell_years_index$year_id))
x_row[cbind(f$cell_years_index$cell_id, f$cell_years_index$year_id)] <-
  seq_len(nrow(f$cell_years_index))
baseline_year <- f$baseline_year
pop_column <- grep("^pop_[a-z]+:g_dom$", colnames(x_cell_years), value = TRUE)
rm(f, par, fold)
invisible(gc())
# g_dom(t) of the "linear_0_1" trend, 0 in 1995 and 1 in 2025
g_dom <- function(year) pmax((year - 1995) / 30, 0)

# The data of one type's pairs for covariate `covariate`: effective counts of
# both years, and pairs x years matrices (padded with 0 beyond each span) of
# the covariate on its 0-1 scale, its trend, the other columns' selection, and
# a mask of the years in the span
type_data <- function(pairs, type, covariate) {
  d <- pairs[pairs$insecticide_type == type, ]
  span <- d$year_id_2 - d$year_id_1
  width <- max(span)
  years <- outer(d$year_id_1, seq_len(width), "+")
  mask <- outer(span, seq_len(width), ">=")
  years[!mask] <- d$year_id_1[row(years)][!mask] + 1
  rows <- matrix(x_row[cbind(rep(d$cell_id, width), as.vector(years))],
                 nrow(d))
  column <- if (covariate == "nets") "nets" else pop_column
  calendar <- baseline_year + years - 1
  g <- if (covariate == "nets") 1 + 0 * calendar else g_dom(calendar)
  x <- matrix(x_cell_years[rows, column], nrow(d))
  if (covariate == "pop") x <- ifelse(g > 0, x / g, 0)
  others <- setdiff(colnames(x_cell_years), column)
  o <- matrix(x_cell_years[rows, others, drop = FALSE] %*%
                effect[others, type], nrow(d))
  inflation <- function(n) 1 + (n - 1) * rho[[type]]
  stopifnot(!anyNA(x[mask]), !anyNA(o[mask]))
  zero_outside <- function(v) ifelse(mask, v, 0)
  list(x = zero_outside(x), g = zero_outside(g), o = zero_outside(o),
       mask = mask, span = span,
       y1 = d$died_1 / inflation(d$n_1), m1 = d$n_1 / inflation(d$n_1),
       y2 = d$died_2 / inflation(d$n_2), m2 = d$n_2 / inflation(d$n_2),
       kappa = kappa[[type]], fit_slope = effect[column, type])
}

# The log likelihood of a type's pairs, with the starting logit of each pair
# profiled out by Newton's method (vectorised over pairs; the binomial log
# likelihood is concave in it)
binomial_loglik <- function(y, m, p) {
  y * log(p) + (m - y) * log1p(-p)
}
profile_loglik <- function(d, selection) {
  # the fall in logit mortality over each pair's span
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

# The maximum likelihood of a type's pairs with knots `knots`, over
# theta = (b, c_1, ..., c_m), the log coefficients
fit_type <- function(d, knots) {
  hinge <- lapply(knots, function(k) pmin(d$x, k))
  negative <- function(theta) {
    selection <- exp(theta[1]) * d$x
    for (j in seq_along(knots)) {
      selection <- selection + exp(theta[j + 1]) * hinge[[j]]
    }
    -profile_loglik(d, selection)
  }
  # Two starts, keeping the better: the type's linear fit with the hinge
  # coefficients near 0 (where the nested linear model ends, so a model with
  # hinges never scores below it), and the hinge coefficients at 0.01 (on the
  # log scale the gradient near 0 is too flat for the hinges to enter from
  # the first)
  slope <- if (is.null(d$linear_theta)) log(max(d$fit_slope, 1e-4)) else
    d$linear_theta
  starts <- if (length(knots) == 0) list(slope) else
    list(c(slope, rep(log(1e-6), length(knots))),
         c(slope, rep(log(1e-2), length(knots))))
  fits <- lapply(starts, function(start) {
    if (length(start) == 1) {
      o <- optim(start, negative, method = "Brent", lower = -15, upper = 5)
    } else {
      o <- optim(start, negative, method = "BFGS",
                 control = list(maxit = 500))
    }
    list(loglik = -o$value, theta = o$par)
  })
  fits[[which.max(vapply(fits, `[[`, numeric(1), "loglik"))]]
}

knot_sets <- function(grid, m) {
  if (m == 0) {
    return(list(numeric(0)))
  }
  sets <- combn(grid, m, simplify = FALSE)
  Filter(function(k) all(diff(k) >= knot_step - 1e-9), sets)
}

# fit every type, in parallel
fit_types <- function(data, knots) {
  parallel::mclapply(data, fit_type, knots = knots, mc.cores = n_cores)
}

# the knots of m hinges maximising the likelihood summed over types, on the
# grid
best_knots <- function(data, m, range_x) {
  if (m == 0) {
    return(numeric(0))
  }
  grid <- seq(ceiling(range_x[1] / knot_step) * knot_step,
              floor(range_x[2] / knot_step) * knot_step, by = knot_step)
  sets <- knot_sets(grid, m)
  loglik <- vapply(sets, function(k) {
    sum(vapply(fit_types(data, k), `[[`, numeric(1), "loglik"))
  }, numeric(1))
  sets[[which.max(loglik)]]
}

all_knots <- list()
all_groups <- list()
for (gap in gaps) {
  stem <- if (gap > 1) sprintf("%s_gap%d", label, gap) else label
  pairs <- read.csv(file.path(output_dir, sprintf("%s_change_pairs.csv",
                                                  stem)))
  if (pair_fraction < 1) {
    set.seed(37)
    pairs <- pairs[unlist(lapply(split(seq_len(nrow(pairs)),
                                       pairs$insecticide_type), function(i) {
      i[sample.int(length(i), ceiling(pair_fraction * length(i)))]
    })), ]
  }
  classes <- setNames(ifelse(pairs$insecticide_type == "DDT", "DDT",
                             pairs$insecticide_class),
                      pairs$insecticide_type)
  for (covariate in c("nets", "pop")) {
    types <- sort(unique(pairs$insecticide_type))
    data <- setNames(lapply(types, function(type) {
      type_data(pairs, type, covariate)
    }), types)
    x_all <- unlist(lapply(data, function(d) d$x[d$mask & d$g > 0]))
    range_x <- quantile(x_all, c(0.01, 0.99))
    linear <- fit_types(data, numeric(0))
    for (type in types) {
      data[[type]]$linear_theta <- linear[[type]]$theta
    }
    for (m in 0:max_hinges) {
      started <- Sys.time()
      knots <- best_knots(data, m, range_x)
      fits <- fit_types(data, knots)
      loglik <- vapply(fits, `[[`, numeric(1), "loglik")
      aic <- -2 * loglik + 2 * (1 + m)
      all_knots[[length(all_knots) + 1]] <- data.frame(
        gap = gap, covariate = covariate, hinges = m,
        knots = paste(sprintf("%.3f", knots), collapse = ","),
        total_aic = sum(aic) + 2 * m,
        seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
      cat(sprintf("gap %d, %s, %d hinges: %.0f s\n", gap, covariate, m,
                  all_knots[[length(all_knots)]]$seconds))
      type_class <- classes[types]
      all_groups[[length(all_groups) + 1]] <- rbind(
        data.frame(gap = gap, covariate = covariate, hinges = m,
                   knots = paste(sprintf("%.3f", knots), collapse = ","),
                   level = "type", group = types,
                   n = vapply(data, function(d) length(d$y1), numeric(1)),
                   aic = aic,
                   coefficients = vapply(fits, function(r) {
                     paste(sprintf("%.6g", exp(r$theta)), collapse = ",")
                   }, character(1))),
        data.frame(gap = gap, covariate = covariate, hinges = m,
                   knots = paste(sprintf("%.3f", knots), collapse = ","),
                   level = "class", group = names(tapply(aic, type_class, sum)),
                   n = as.vector(tapply(vapply(data, function(d) {
                     length(d$y1)
                   }, numeric(1)), type_class, sum)),
                   aic = as.vector(tapply(aic, type_class, sum)),
                   coefficients = ""))
    }
    cat(sprintf("gap %d, %s done\n", gap, covariate))
  }
}
knots <- do.call(rbind, all_knots)
groups <- do.call(rbind, all_groups)
key <- list(groups$gap, groups$covariate, groups$level, groups$group)
groups <- do.call(rbind, lapply(split(groups, key, drop = TRUE), function(r) {
  r$delta_aic <- r$aic - r$aic[r$hinges == 0]
  r
}))
rownames(groups) <- NULL
selected <- do.call(rbind, lapply(
  split(groups, list(groups$gap, groups$covariate, groups$level,
                     groups$group), drop = TRUE),
  function(r) r[which.min(r$aic), ]))
rownames(selected) <- NULL
write.csv(knots, file.path(output_dir, sprintf("%s_mle_knots.csv", output_label)),
          row.names = FALSE)
write.csv(groups, file.path(output_dir, sprintf("%s_mle_by_group.csv",
                                                output_label)), row.names = FALSE)
write.csv(selected, file.path(output_dir, sprintf("%s_mle_selected.csv",
                                                  output_label)), row.names = FALSE)
options(width = 200)
print(knots)
