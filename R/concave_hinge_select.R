# Monotone concave hinge selection (#37): the number and places of hinge
# knots for a selection covariate, fitted with the constraint the dynamical
# model imposes (every selection coefficient exp(beta) > 0, so the effect of a
# covariate, a x + sum_j b_j min(x, k_j), is increasing and concave), on the
# observed within-pixel changes of a fit's pairs (R/change_residual_screen.R).
#
#   Rscript R/concave_hinge_select.R <fit RData> <label> [gap ...]
#
# reads outputs/hinge_screen/<label>[_gap<g>]_change_pairs.csv and
# _change_pair_years.csv for each gap (default 1, 3, 5).
#
# The outcome of a pair is the partial residual of the covariate x: the
# observed fall in the empirical logit of mortality per year, minus the fit's
# predicted fall, plus the fit's own effect of x, e_x mean_t(g_t x_t), with
# e_x the type's posterior mean selection effect exp(beta) and g_t the
# column's time trend (g_dom(t) for population, 1 for net use). With no floor
# the model's logit mortality falls each year by log(1 + sum_k x_k exp(beta_k))
# plus reversion, about sum_k x_k exp(beta_k) for small effects, so the
# partial residual is the fall due to x, with the fit's other terms removed.
# It is regressed, weighted by its precision (the empirical-logit variances of
# the two years, inflated by the type's beta-binomial rho), on a free
# intercept (reversion and what the linearisation leaves) and the columns
# mean_t(g_t x_t) and mean_t(g_t min(x_t, k_j)), as the model builds them
# yearly, with non-negative coefficients (nnls_fit()).
#
# Knots are shared across all insecticides; the coefficients are fitted per
# type. For 0 to 3 hinges the knots minimise the weighted residual sum of
# squares summed over types: the best sets on a grid (steps of 0.025, within
# the 2.5-97.5% range of x, at least 0.025 apart) are refined by Nelder-Mead.
# Each type's score is its AIC, n log(WRSS / n) + 2 (columns), a class's the
# sum over its types, and the shared knots are counted once in the total.
#
# Writes outputs/hinge_screen/<label>_concave_hinges.csv (the knots of each
# model), <label>_concave_by_group.csv (each type's and class's AIC change
# from the linear model) and <label>_concave_selected.csv (the best number of
# hinges of each).

arguments <- commandArgs(trailingOnly = TRUE)
fit_file <- arguments[1]
label <- arguments[2]
gaps <- if (length(arguments) >= 3) as.integer(arguments[-(1:2)]) else
  c(1L, 3L, 5L)
output_dir <- "outputs/hinge_screen"
max_hinges <- 3
knot_step <- 0.025
n_refine <- 5

source("R/packages.R")
source("R/functions.R")
source("R/dynamical_model.R")
source("R/dynamical_predictions.R")

# each type's posterior mean selection effect exp(beta) of nets and of
# population
f <- new.env()
load(fit_file, envir = f)
fold <- list(draws = f$draws, options = f$model_options,
             x_cells_init = f$x_cells_init)
n_total <- sum(vapply(fold$draws, nrow, integer(1)))
par <- dynamical_parameter_draws(fold, f$classes_index, f$types,
                                 draw_index = round(seq(1, n_total,
                                                        length.out = 400)))
effect <- apply(par$effect_type, c(2, 3), mean)
dimnames(effect) <- list(colnames(f$x_cell_years), f$types)
pop_column <- grep("^pop_[a-z]+:g_dom$", rownames(effect), value = TRUE)
fit_effect <- list(nets = effect["nets", ], pop = effect[pop_column, ])
rm(f, par, fold)
invisible(gc())

# Non-negative least squares (Lawson and Hanson 1974, ch. 23): the b >= 0
# minimising |A b - y|^2, by an active set of positive coefficients grown one
# column at a time and shrunk whenever an unconstrained solution would go
# negative. Exact for the handful of columns here.
nnls_fit <- function(A, y, tolerance = 1e-10) {
  p <- ncol(A)
  b <- numeric(p)
  passive <- logical(p)
  w <- crossprod(A, y - A %*% b)
  for (iteration in seq_len(3 * p + 10)) {
    if (all(passive) || max(w[!passive]) <= tolerance) break
    j <- which(!passive)[which.max(w[!passive])]
    passive[j] <- TRUE
    repeat {
      s <- numeric(p)
      s[passive] <- qr.solve(A[, passive, drop = FALSE], y)
      if (all(s[passive] > tolerance)) break
      negative <- passive & s <= tolerance
      alpha <- min(b[negative] / (b[negative] - s[negative]))
      b <- b + alpha * (s - b)
      passive <- passive & b > tolerance
      b[!passive] <- 0
    }
    b <- s
    w <- crossprod(A, y - A %*% b)
  }
  list(coefficients = b, rss = sum((y - A %*% b) ^ 2))
}

# The design of a type's pairs: the intercept as two non-negative columns
# (+1 and -1, so it is free), then mean_t(g_t x_t) and mean_t(g_t min(x_t,
# k_j)) over each pair's years
pair_columns <- function(years, knots) {
  columns <- cbind(x = years$g * years$x,
                   vapply(knots, function(k) years$g * pmin(years$x, k),
                          numeric(nrow(years))))
  rowsum(columns, years$pair, reorder = FALSE) /
    as.vector(table(factor(years$pair, unique(years$pair))))
}

# weighted residual sum of squares and number of columns of each type, for
# the knots `knots`
type_fits <- function(data, knots) {
  t(vapply(data, function(d) {
    X <- cbind(1, -1, pair_columns(d$years, knots))
    fit <- nnls_fit(X * sqrt(d$w), d$y * sqrt(d$w))
    c(wrss = fit$rss, columns = ncol(X) - 1, n = length(d$y))
  }, numeric(3)))
}

knot_sets <- function(grid, m) {
  if (m == 0) {
    return(list(numeric(0)))
  }
  sets <- combn(grid, m, simplify = FALSE)
  Filter(function(k) all(diff(k) >= knot_step - 1e-9), sets)
}

# the knots of m hinges minimising the total weighted RSS: the best grid sets
# refined by Nelder-Mead, with the knots kept sorted and inside the range
best_knots <- function(data, m, range_x) {
  if (m == 0) {
    return(numeric(0))
  }
  grid <- seq(ceiling(range_x[1] / knot_step) * knot_step,
              floor(range_x[2] / knot_step) * knot_step, by = knot_step)
  sets <- knot_sets(grid, m)
  total <- function(k) sum(type_fits(data, k)[, "wrss"])
  rss <- vapply(sets, total, numeric(1))
  starts <- sets[order(rss)[seq_len(min(n_refine, length(sets)))]]
  objective <- function(k) {
    if (any(k < range_x[1] | k > range_x[2])) return(Inf)
    total(sort(k))
  }
  refined <- lapply(starts, function(k) {
    if (m == 1) {
      o <- optimize(objective, c(max(range_x[1], k - knot_step),
                                 min(range_x[2], k + knot_step)))
      list(par = o$minimum, value = o$objective)
    } else {
      optim(k, objective, method = "Nelder-Mead",
            control = list(reltol = 1e-10, maxit = 400))
    }
  })
  sort(refined[[which.min(vapply(refined, `[[`, numeric(1), "value"))]]$par)
}

all_knots <- list()
all_groups <- list()
for (gap in gaps) {
  stem <- if (gap > 1) sprintf("%s_gap%d", label, gap) else label
  pairs <- read.csv(file.path(output_dir, sprintf("%s_change_pairs.csv",
                                                  stem)))
  years <- read.csv(file.path(output_dir, sprintf("%s_change_pair_years.csv",
                                                  stem)))
  span <- pairs$year_2 - pairs$year_1
  # the empirical-logit variance of each year's pooled assays
  elogit_var <- function(died, n, rho) {
    (1 / (died + 0.5) + 1 / (n - died + 0.5)) * (1 + (n - 1) * rho)
  }
  pairs$w <- span ^ 2 / (elogit_var(pairs$died_1, pairs$n_1, pairs$rho) +
                           elogit_var(pairs$died_2, pairs$n_2, pairs$rho))
  pairs$fall_obs <- -pairs$change_obs / span
  pairs$fall_pred <- -pairs$change_pred / span
  pairs$class <- ifelse(pairs$insecticide_type == "DDT", "DDT",
                        pairs$insecticide_class)
  for (covariate in c("nets", "pop")) {
    years$x <- years[[covariate]]
    years$g <- if (covariate == "pop") years$g_dom else 1
    fit_column <- rowsum(years$g * years$x, years$pair, reorder = FALSE) /
      as.vector(table(factor(years$pair, unique(years$pair))))
    stopifnot(identical(as.integer(rownames(fit_column)), pairs$pair))
    pairs$y <- pairs$fall_obs - pairs$fall_pred +
      fit_effect[[covariate]][pairs$insecticide_type] * fit_column[, 1]
    data <- lapply(split(seq_len(nrow(pairs)), pairs$insecticide_type),
                   function(i) {
      list(y = pairs$y[i], w = pairs$w[i],
           years = years[years$pair %in% pairs$pair[i], c("pair", "x", "g")])
    })
    range_x <- quantile(years$x[years$g > 0], c(0.025, 0.975))
    for (m in 0:max_hinges) {
      knots <- best_knots(data, m, range_x)
      g <- type_fits(data, knots)
      aic <- g[, "n"] * log(g[, "wrss"] / g[, "n"]) + 2 * g[, "columns"]
      all_knots[[length(all_knots) + 1]] <- data.frame(
        gap = gap, covariate = covariate, hinges = m,
        knots = paste(sprintf("%.3f", knots), collapse = ","),
        total_aic = sum(aic) + 2 * m)
      types <- rownames(g)
      classes <- pairs$class[match(types, pairs$insecticide_type)]
      all_groups[[length(all_groups) + 1]] <- rbind(
        data.frame(gap = gap, covariate = covariate, hinges = m,
                   level = "type", group = types, n = g[, "n"], aic = aic),
        data.frame(gap = gap, covariate = covariate, hinges = m,
                   level = "class", group = names(tapply(aic, classes, sum)),
                   n = as.vector(tapply(g[, "n"], classes, sum)),
                   aic = as.vector(tapply(aic, classes, sum))))
    }
    cat(sprintf("gap %d, %s done\n", gap, covariate))
  }
}
knots <- do.call(rbind, all_knots)
groups <- do.call(rbind, all_groups)
groups <- do.call(rbind, lapply(
  split(groups, list(groups$gap, groups$covariate, groups$level,
                     groups$group), drop = TRUE),
  function(r) {
    r$delta_aic <- r$aic - r$aic[r$hinges == 0]
    r
  }))
rownames(groups) <- NULL
selected <- do.call(rbind, lapply(
  split(groups, list(groups$gap, groups$covariate, groups$level,
                     groups$group), drop = TRUE),
  function(r) r[which.min(r$aic), ]))
rownames(selected) <- NULL
write.csv(knots, file.path(output_dir, sprintf("%s_concave_hinges.csv",
                                               label)), row.names = FALSE)
write.csv(groups, file.path(output_dir, sprintf("%s_concave_by_group.csv",
                                                label)), row.names = FALSE)
write.csv(selected, file.path(output_dir, sprintf("%s_concave_selected.csv",
                                                  label)), row.names = FALSE)
options(width = 200)
print(knots)
