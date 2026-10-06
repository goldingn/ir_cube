# Hinge selection (#37): the number and places of hinge knots for a selection
# covariate, chosen on the change residuals of a fit
# (R/change_residual_screen.R).
#
#   Rscript R/hinge_select.R <pairs csv> [covariate ...]
#
# For each covariate x (default nets and pop) separately, OLS regressions of
# the change residual z with 0 to 3 hinges, fitted to each group g (an
# insecticide class, or type) separately:
#
#   z ~ 1 + x + sum_j min(x, k_j)
#
# hinge columns min(x, k) as in selection_design(hinges = ...), with the knots
# k_j shared across all groups. The knots of each model are those minimising
# the residual sum of squares, summed over groups, over a grid (steps of
# 0.025, within the 2.5-97.5% range of x, at least 0.05 apart). The score is
# the groups' AICs pooled, sum_g n_g log(RSS_g / n_g) + 2 (coefficients +
# knots), each group with its own residual variance; BIC likewise. Run with
# the groups as classes (DDT on its own) and as types. x is the covariate's
# mean over the years the change spans, so a hinge here is of that mean, not
# the mean of yearly hinges as in the model.
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
levels_of <- list(
  class = ifelse(pairs$insecticide_type == "DDT", "DDT",
                 pairs$insecticide_class),
  type = pairs$insecticide_type)

# the residual sum of squares and rank of the model with knots `knots`
# (possibly none) in each group
group_fits <- function(x, z, knots, rows) {
  t(vapply(rows, function(i) {
    design <- cbind(1, x[i], vapply(knots, function(k) pmin(x[i], k),
                                    numeric(length(i))))
    fit <- lm.fit(design, z[i])
    c(rss = sum(fit$residuals ^ 2), rank = fit$rank, n = length(i))
  }, numeric(3)))
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
for (level in names(levels_of)) {
  rows <- split(seq_len(nrow(pairs)), levels_of[[level]])
  for (covariate in covariates) {
    x <- pairs[[covariate]]
    range_x <- quantile(x, c(0.025, 0.975))
    grid <- seq(ceiling(range_x[1] / knot_step) * knot_step,
                floor(range_x[2] / knot_step) * knot_step, by = knot_step)
    for (m in 0:max_hinges) {
      sets <- knot_sets(grid, m)
      rss <- vapply(sets, function(k) {
        sum(group_fits(x, pairs$z, k, rows)[, "rss"])
      }, numeric(1))
      best <- sets[[which.min(rss)]]
      g <- group_fits(x, pairs$z, best, rows)
      parameters <- sum(g[, "rank"]) + m
      log_rss <- sum(g[, "n"] * log(g[, "rss"] / g[, "n"]))
      results[[length(results) + 1]] <- data.frame(
        level = level, covariate = covariate, hinges = m,
        knots = if (m == 0) "" else paste(sprintf("%.3f", best),
                                          collapse = ","),
        rss = sum(g[, "rss"]), parameters = parameters, n = sum(g[, "n"]),
        aic = log_rss + 2 * parameters,
        bic = log_rss + log(sum(g[, "n"])) * parameters,
        # each group's term of the pooled AIC, to read which groups drive it
        group_aic = paste(sprintf("%s %.1f", rownames(g),
                                  g[, "n"] * log(g[, "rss"] / g[, "n"]) +
                                    2 * g[, "rank"]), collapse = "; "))
    }
  }
}
results <- do.call(rbind, results)
results <- do.call(rbind, lapply(
  split(results, list(results$level, results$covariate), drop = TRUE),
  function(r) {
    r$delta_aic <- r$aic - min(r$aic)
    r$delta_bic <- r$bic - min(r$bic)
    r
  }))
rownames(results) <- NULL
write.csv(results, sub("\\.csv$", "_hinges.csv", pairs_file), row.names = FALSE)
options(width = 200)
print(results[, c("level", "covariate", "hinges", "knots", "parameters",
                  "delta_aic", "delta_bic")], digits = 4)
for (level in names(levels_of)) {
  for (covariate in covariates) {
    r <- results[results$level == level & results$covariate == covariate, ]
    chosen <- r[which.min(r$aic), ]
    cat(sprintf("%s level, %s: AIC selects %s\n", level, covariate,
                if (chosen$hinges == 0) "linear" else
                  sprintf("%d hinge(s) at %s", chosen$hinges, chosen$knots)))
  }
}
