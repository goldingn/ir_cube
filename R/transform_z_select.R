# Screen of fixed covariate transforms against the change residuals z of a
# fit (#37; R/change_residual_screen.R): which encounter transform
# T(x) = 1 - exp(-x log(2) / x_half) of net use, or of population density,
# best explains the misfit of the fit's own column x.
#
#   Rscript R/transform_z_select.R <label> [gap]
#
# reads outputs/hinge_screen/<label>[_gap<g>]_change_pairs.csv and
# _change_pair_years.csv (default gap 5). If the true effect is c T(x) and the
# fit has e x, the misfit in z is a mix of x and T(x), so each candidate is
# scored by OLS of z on, per insecticide type, an intercept, the fit's column
# and the candidate column (each a mean over the pair's years, the population
# columns times the time trend g_dom), as R/hinge_select.R does for hinges.
# Every candidate has the same number of parameters; the score is AIC pooled
# over the types, n log(RSS / n) summed, by class and type. Population
# density is recovered from the fit's d_half = 270 column.
#
# Writes outputs/hinge_screen/<label>[_gap<g>]_transform_z.csv.

arguments <- commandArgs(trailingOnly = TRUE)
label <- arguments[1]
gap <- if (length(arguments) >= 2) as.integer(arguments[2]) else 5L
fit_d_half <- 270
candidates <- list(
  nets = c(0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 1, 2),
  pop = c(25, 50, 100, 200, 270, 500, 1000, 3000, 10000))

stem <- if (gap > 1) sprintf("%s_gap%d", label, gap) else label
pairs <- read.csv(file.path("outputs/hinge_screen",
                            sprintf("%s_change_pairs.csv", stem)))
years <- read.csv(file.path("outputs/hinge_screen",
                            sprintf("%s_change_pair_years.csv", stem)))
encounter <- function(x, x_half) -expm1(-x * log(2) / x_half)
# density from the encounter column, capped where it is numerically 1
years$density <- -fit_d_half * log(pmax(1 - years$pop, 1e-12)) / log(2)
span_mean <- function(v) tapply(v, years$pair, mean)[as.character(pairs$pair)]

columns <- list(
  nets = list(fit = span_mean(years$nets),
              candidate = function(h) span_mean(encounter(years$nets, h))),
  pop = list(fit = span_mean(years$g_dom * years$pop),
             candidate = function(h) {
               span_mean(years$g_dom * encounter(years$density, h))
             }))
groups <- list(class = ifelse(pairs$insecticide_type == "DDT", "DDT",
                              pairs$insecticide_class),
               type = pairs$insecticide_type)

# n log(RSS / n) of the OLS of z on 1, the fit column and the candidate
# column, in each group
group_scores <- function(fit_column, candidate_column, group) {
  vapply(split(seq_along(pairs$z), group), function(i) {
    X <- cbind(1, fit_column[i], candidate_column[i])
    r <- lm.fit(X, pairs$z[i])$residuals
    length(i) * log(sum(r ^ 2) / length(i))
  }, numeric(1))
}

results <- list()
for (covariate in names(columns)) {
  for (h in candidates[[covariate]]) {
    tx <- columns[[covariate]]$candidate(h)
    for (level in names(groups)) {
      s <- group_scores(columns[[covariate]]$fit, tx, groups[[level]])
      results[[length(results) + 1]] <- data.frame(
        covariate = covariate, x_half = h, level = level, group = names(s),
        score = s)
    }
  }
}
results <- do.call(rbind, results)
# the score relative to the best candidate in each group (0 = best; a
# difference of 2 is the AIC worth of one parameter)
results$delta <- ave(results$score, results$covariate, results$level,
                     results$group, FUN = function(v) v - min(v))
write.csv(results, file.path("outputs/hinge_screen",
                             sprintf("%s_transform_z.csv", stem)),
          row.names = FALSE)
for (covariate in names(columns)) {
  for (level in names(groups)) {
    r <- results[results$covariate == covariate & results$level == level, ]
    w <- reshape(r[, c("x_half", "group", "delta")], idvar = "x_half",
                 timevar = "group", direction = "wide")
    names(w) <- sub("delta.", "", names(w))
    total <- tapply(r$score, r$x_half, sum)
    w$total <- round(total[as.character(w$x_half)] - min(total), 1)
    w[, -1] <- round(w[, -1], 1)
    cat(sprintf("\n%s, %s level (0 = best)\n", covariate, level))
    print(w, row.names = FALSE)
  }
}
