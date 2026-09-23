# Does predictive skill depend on distance to the nearest training observation?
#
# The answer, on the sub-national block folds, is that it cannot be determined:
# no model shows a distance effect that survives a joint test once the clustering
# of bioassays by pixel is accounted for. This script is the record of how that
# was established, and produces outputs/cv_distance.csv.
#
# The question matters because the raw distance-stratified skill table looks
# structured - a dip at short range, a peak near 200 km - but the strata are
# confounded: records within 50 km are 48% Kenya, records beyond 400 km are 60%
# Tanzania, and how far a held-out record sits from training data is largely a
# function of which country it is in and how large that country's blocks are.
#
# Three design decisions, each forced by something the data do:
#
# 1. The response is the per-record difference in squared error between a model
#    and the intercept null. Writing y = p + e for bioassay noise e, the e^2
#    terms are identical in the two squared errors and cancel, so no noise floor
#    has to be estimated and assay size stops being a confounder. The same
#    cancellation holds in aggregate: excess_model - excess_null is just
#    MSE_model - MSE_null.
#
# 2. Records are aggregated to (cell, insecticide, year), which is the level at
#    which both models' predictions are constant. Because the difference is
#    affine in the observed proportion, the group mean is an exact function of
#    (n, mean(y)): the reduction from 8,690 records to 6,987 strata loses
#    nothing. It also removes the replicate structure that a smoother would
#    otherwise absorb as if it were a distance trend.
#
# 3. There is no pixel term in the mean model, and this is deliberate. Distance
#    is very nearly a pixel-level covariate - its intraclass correlation by
#    pixel is 0.92 - so a random intercept per pixel absorbs the between-pixel
#    contrast that identifies the distance effect. Fitting one collapses the
#    effect to nothing, as does an explicitly between-pixel term in a Mundlak
#    decomposition. The clustering is therefore handled in the inference, by
#    resampling pixels, rather than in the mean.
#
# Distance enters as a five-level factor rather than a spline. An earlier
# version used a penalised thin-plate spline, but its effective degrees of
# freedom tracked the basis dimension (5.6, 9.6, 13.0 at k = 8, 20, 30) because
# distance is near-collinear with pixel identity, so a flexible smooth partly
# smooths pixel. Bins reach the same conclusion with no penalisation, no
# smoothing-parameter selection, and nothing to tune.

source("R/validation_functions.R")
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
})

set.seed(2026 - 9 - 24)
n_bootstrap <- 2000
distance_breaks <- c(0, 50, 100, 200, 400, Inf)

draws_dir <- "outputs/cv_draws"
experiment <- "spatial_blocks"

models <- c(nearest_neighbour = "1-NN",
            nearest_neighbour_oracle = "best-k NN",
            dynamical = "dynamical")

# per-record predictions, observations and geometry ------------------------

scores <- readRDS("outputs/predictable_variance_scores.RDS")

intercept <- bind_rows(lapply(
  unique(scores$fold[scores$experiment == experiment]),
  function(fold) {
    x <- readRDS(file.path(draws_dir, sprintf("intercept__%s__%s.rds",
                                              experiment, fold)))
    data.frame(experiment = experiment, fold = fold,
               row = seq_len(nrow(x$test_df)),
               p_intercept = colMeans(x$p_draws))
  }))

records <- scores %>%
  select(experiment, fold, row, cell, observed, died, mosquito_number, rho,
         country_name, insecticide_class, insecticide_type, year_start,
         distance_1, model, predicted) %>%
  pivot_wider(id_cols = c(experiment, fold, row, cell, observed, died,
                          mosquito_number, rho, country_name, insecticide_class,
                          insecticide_type, year_start, distance_1),
              names_from = model, values_from = predicted,
              names_prefix = "p_") %>%
  left_join(intercept, by = c("experiment", "fold", "row")) %>%
  filter(experiment == !!experiment,
         # countries with too few held-out records to support a fixed effect
         country_name %in% names(which(table(country_name) >= 50)))
stopifnot(!any(is.na(records$p_intercept)))

for (model in names(models)) {
  records[[paste0("d_", model)]] <-
    (records$observed - records[[paste0("p_", model)]]) ^ 2 -
    (records$observed - records$p_intercept) ^ 2
}

# the denominator: the intercept null's excess mean squared error over the whole
# experiment. One fixed number, so every distance bin is on a common scale
floor_all <- noise_floor_mse(records$died, records$mosquito_number, records$rho)
excess_null <- mean((records$observed - records$p_intercept) ^ 2) - floor_all
to_percent <- -100 / excess_null

# exact aggregation --------------------------------------------------------

strata <- records %>%
  group_by(cell, insecticide_type, year_start) %>%
  summarise(across(starts_with("d_"), mean),
            n = n(),
            distance = first(distance_1),
            country = first(country_name),
            insecticide = first(insecticide_class),
            fold = first(fold),
            .groups = "drop") %>%
  mutate(bin = cut(distance, distance_breaks, include.lowest = TRUE,
                   dig.lab = 4),
         year = factor(year_start),
         across(c(country, insecticide, fold, cell), factor))

cat(sprintf("%d records -> %d strata in %d pixels\n",
            nrow(records), nrow(strata), n_distinct(strata$cell)))

cat("\nstrata and pixels per distance bin:\n")
print(as.data.frame(strata %>% group_by(bin) %>%
  summarise(strata = n(), pixels = n_distinct(cell), assays = sum(n),
            median_km = round(median(distance)), .groups = "drop")),
  row.names = FALSE)

# the confounding that motivates the adjustment
icc <- function(y, g) {
  g <- factor(g); k <- nlevels(g); n <- length(y)
  ni <- as.numeric(table(g)); mi <- as.numeric(tapply(y, g, mean))
  msb <- sum(ni * (mi - mean(y)) ^ 2) / (k - 1)
  msw <- sum((y - mi[as.integer(g)]) ^ 2) / (n - k)
  m0 <- (n - sum(ni ^ 2) / n) / (k - 1)
  between <- max((msb - msw) / m0, 0)
  between / (between + msw)
}
cat(sprintf("\nintraclass correlation of log distance by pixel: %.3f\n",
            icc(log1p(strata$distance), strata$cell)))
cat("(this is why there is no pixel term in the mean model)\n")

# binned weighted least squares, with a pixel-cluster bootstrap -------------

cells <- levels(strata$cell)
rows_by_cell <- split(seq_len(nrow(strata)), strata$cell)

# standardised marginal mean per bin: every stratum is set to that bin, keeping
# its own country, insecticide, fold and year, and predictions are averaged with
# weights. Handles the empty country-by-distance cells that reweighting cannot
marginal <- function(fit, data) {
  vapply(levels(strata$bin), function(b) {
    newdata <- data
    newdata$bin <- factor(b, levels = levels(strata$bin))
    weighted.mean(predict(fit, newdata = newdata), newdata$n)
  }, numeric(1))
}

results <- list()
for (model in names(models)) {

  formula <- as.formula(sprintf(
    "d_%s ~ bin + country + insecticide + fold + year", model))
  fit <- lm(formula, data = strata, weights = n)
  point <- marginal(fit, strata)

  replicates <- t(replicate(n_bootstrap, {
    picked <- sample(cells, length(cells), replace = TRUE)
    resample <- strata[unlist(rows_by_cell[picked]), ]
    refit <- try(lm(formula, data = resample, weights = n), silent = TRUE)
    if (inherits(refit, "try-error")) {
      return(rep(NA_real_, nlevels(strata$bin)))
    }
    marginal(refit, resample)
  }))

  # joint test that the bins differ at all, against the last bin, using the
  # bootstrap covariance of the contrasts
  reference <- nlevels(strata$bin)
  contrasts <- to_percent * (point[-reference] - point[reference])
  contrast_replicates <- to_percent *
    (replicates[, -reference, drop = FALSE] - replicates[, reference])
  contrast_replicates <-
    contrast_replicates[complete.cases(contrast_replicates), , drop = FALSE]
  wald <- as.numeric(t(contrasts) %*% solve(cov(contrast_replicates)) %*%
                       contrasts)
  p_value <- pchisq(wald, reference - 1, lower.tail = FALSE)

  results[[model]] <- data.frame(
    model = models[model],
    bin = levels(strata$bin),
    estimate = to_percent * point,
    lower = to_percent * apply(replicates, 2, quantile, 0.975, na.rm = TRUE),
    upper = to_percent * apply(replicates, 2, quantile, 0.025, na.rm = TRUE),
    joint_chisq = wald,
    joint_df = reference - 1,
    joint_p = p_value,
    row.names = NULL)
}

distance_table <- bind_rows(results)
write.csv(distance_table, "outputs/cv_distance.csv", row.names = FALSE)

cat("\npercentage of the intercept null's excess mean squared error removed,",
    "\nwith 95% pixel-cluster bootstrap intervals:\n")
print(as.data.frame(distance_table %>%
  mutate(interval = sprintf("%6.1f [%6.1f, %6.1f]", estimate, lower, upper)) %>%
  select(bin, model, interval) %>%
  pivot_wider(names_from = model, values_from = interval)), row.names = FALSE)

cat("\njoint test that skill differs across distance bins:\n")
print(as.data.frame(distance_table %>%
  distinct(model, joint_chisq, joint_df, joint_p) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

# what a pixel term would do, shown rather than asserted ---------------------

cat("\nbin coefficients relative to the first bin, under three treatments of",
    "\nthe pixel, for the dynamical model (% of the null's excess MSE):\n")
formula <- d_dynamical ~ bin + country + insecticide + fold + year
show_bins <- function(coefficients, label) {
  keep <- grep("^bin", names(coefficients))
  cat(sprintf("  %-30s %s\n", label,
              paste(sprintf("%s %6.1f", sub("^bin", "", names(coefficients)[keep]),
                            to_percent * coefficients[keep]), collapse = "  ")))
}
show_bins(coef(lm(formula, data = strata, weights = n)),
          "pixel in the bootstrap only")
show_bins(fixef(lmer(update(formula, . ~ . + (1 | cell)), data = strata,
                     weights = n, REML = TRUE)),
          "random intercept per pixel")
show_bins(coef(lm(update(formula, . ~ . - country + cell), data = strata,
                  weights = n)),
          "pixel as a fixed factor")
cat("\nThe random intercept shrinks every contrast toward zero and the fixed\n",
    "factor scrambles them: both are conditioning on a factor the covariate is\n",
    "nearly nested within. This is the same effect seen with a spline, without\n",
    "any smoothing machinery involved.\n", sep = "")
