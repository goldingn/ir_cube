# Hinge selection (#37): the number and places of hinge knots for a selection
# covariate, chosen on the change residuals of a fit
# (R/change_residual_screen.R).
#
#   Rscript R/hinge_select.R <pairs csv> [covariate ...]
#
# For each covariate x (default nets and pop) separately, OLS regressions of
# the change residual z with 0 to 3 hinges:
#
#   z ~ type + type:x + sum_j type:min(x, k_j)
#
# hinge columns min(x, k) as in selection_design(hinges = ...), with an
# intercept, a slope and a coefficient on each hinge column per insecticide
# type, and the knots k_j shared across types. The knots of each model are
# those minimising the residual sum of squares over a grid (steps of 0.025,
# within the 2.5-97.5% range of x, at least 0.05 apart). Models are compared
# by AIC = n log(RSS / n) + 2 (coefficients + knots); BIC is also given. x is
# the covariate's mean over the years the change spans, so a hinge here is of
# that mean, not the mean of yearly hinges as in the model.
#
# Writes <pairs csv without .csv>_hinges.csv and prints the comparison.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) >= 1)
pairs_file <- arguments[1]
covariates <- if (length(arguments) >= 2) arguments[-1] else c("nets", "pop")
max_hinges <- 3
knot_step <- 0.025
min_spacing <- 0.05

pairs <- read.csv(pairs_file)
pairs <- pairs[is.finite(pairs$z), ]
type <- factor(pairs$insecticide_type)
n <- nrow(pairs)

# coefficients per insecticide type (default), or with IR_CUBE_HINGE_COEFS =
# "shared", one intercept, slope and hinge coefficient for all types
shared <- identical(Sys.getenv("IR_CUBE_HINGE_COEFS"), "shared")
type_indicators <- if (shared) matrix(1, n, 1) else model.matrix(~ type - 1)

# the design of the model with knots `knots` (possibly none) on x: each of
# the intercept, x and the hinge columns, per type
hinge_design <- function(x, knots) {
  columns <- cbind(1, x, vapply(knots, function(k) pmin(x, k),
                                numeric(length(x))))
  do.call(cbind, lapply(seq_len(ncol(columns)), function(j) {
    type_indicators * columns[, j]
  }))
}

# residual sum of squares of the model with knots `knots`, and its number of
# estimable coefficients
hinge_rss <- function(x, z, knots) {
  fit <- lm.fit(hinge_design(x, knots), z)
  c(rss = sum(fit$residuals ^ 2), coefficients = fit$rank)
}

# the knot sets on the grid with m knots, at least min_spacing apart
knot_sets <- function(grid, m) {
  if (m == 0) {
    return(list(numeric(0)))
  }
  sets <- combn(grid, m, simplify = FALSE)
  Filter(function(k) all(diff(k) >= min_spacing - 1e-9), sets)
}

results <- list()
for (covariate in covariates) {
  x <- pairs[[covariate]]
  range_x <- quantile(x, c(0.025, 0.975))
  grid <- seq(ceiling(range_x[1] / knot_step) * knot_step,
              floor(range_x[2] / knot_step) * knot_step, by = knot_step)
  for (m in 0:max_hinges) {
    sets <- knot_sets(grid, m)
    rss <- vapply(sets, function(k) hinge_rss(x, pairs$z, k)["rss"],
                  numeric(1))
    best <- sets[[which.min(rss)]]
    fit <- hinge_rss(x, pairs$z, best)
    parameters <- fit[["coefficients"]] + m
    results[[length(results) + 1]] <- data.frame(
      covariate = covariate, hinges = m,
      knots = if (m == 0) "" else paste(sprintf("%.3f", best), collapse = ","),
      rss = fit[["rss"]], parameters = parameters, n = n,
      aic = n * log(fit[["rss"]] / n) + 2 * parameters,
      bic = n * log(fit[["rss"]] / n) + log(n) * parameters,
      knot_sets_searched = length(sets))
  }
}
results <- do.call(rbind, results)
results <- do.call(rbind, lapply(split(results, results$covariate), function(r) {
  r$delta_aic <- r$aic - min(r$aic)
  r$delta_bic <- r$bic - min(r$bic)
  r
}))
rownames(results) <- NULL
write.csv(results, sub("\\.csv$", if (shared) "_hinges_shared.csv" else
                         "_hinges.csv", pairs_file), row.names = FALSE)
options(width = 200)
print(results[, c("covariate", "hinges", "knots", "rss", "parameters",
                  "delta_aic", "delta_bic")], digits = 4)
for (covariate in covariates) {
  r <- results[results$covariate == covariate, ]
  chosen <- r[which.min(r$aic), ]
  cat(sprintf("%s: AIC selects %d hinge(s)%s\n", covariate, chosen$hinges,
              if (chosen$hinges > 0) paste0(" at ", chosen$knots) else
                ", i.e. linear"))
}
