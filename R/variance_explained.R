# Out-of-sample variance explained, and the share of observed variance that
# bioassay sampling makes unexplainable, for every fitted fold and model.
#
# Two quantities, deliberately kept apart because they carry different kinds of
# uncertainty and the first must not depend on the second:
#
#   explained_i = 100 (1 - MSE_model / Var(y))
#       the fraction of observed variance in held-out mortality the model
#       accounts for. No noise floor enters it, so it does not inherit the
#       uncertainty in rho. Interval by resampling pixels, since bioassays
#       cluster hard by pixel and the assay count badly overstates the
#       information a fold holds.
#
#   noise = 100 floor / Var(y)
#       the share of that variance attributable to beta-binomial sampling in
#       the assay itself, which no model can explain. Interval from the
#       posterior of rho, estimated per insecticide type from replicate
#       bioassays (fig_illustrate_bioassay_variability.R).
#
#   ceiling = 100 - noise
#       the most any model could explain.
#
# The floor is the sampling variance of the observed proportion, p(1-p) k_i,
# with k_i = (1 + (m_i - 1) rho_t) / m_i the beta-binomial design effect over
# the assay size. What p(1-p) should be is the awkward part. Plugging in the
# observed proportion gives exactly zero wherever an assay reads 0% or 100% -
# which, with this much overdispersion, is a large share of the records for the
# insecticides that sit against the top of the scale - and it is noisy
# everywhere else. So p is instead regularised: each assay carries the
# information of m_eff = 1/k_i independent draws, so under a Beta(a_t, b_t)
# prior on the true proportion the posterior is
#   Beta(a_t + y_i m_eff, b_t + (1 - y_i) m_eff)
# and its mean of p(1-p) is AB / ((A+B)(A+B+1)), which is never zero. The prior
# is fitted by maximum likelihood per insecticide type over all of that type's
# held-out records, so the floor stays a property of the record rather than of
# whichever subset is being summarised.
#
# This is the honest estimator, not a fix for the saturated insecticides: for
# Fenitrothion in the interpolation fold the floor exceeds the observed
# variance however p is estimated, because a constant rho overstates assay
# noise as mortality approaches 100%. Those cells are excluded from the
# per-insecticide figure by the rule below rather than papered over.
#
# Writes outputs/cv_variance_explained.csv (pooled per experiment) and
# outputs/cv_variance_explained_by_insecticide.csv.

source("R/validation_functions.R")
suppressMessages({
  library(dplyr)
  library(tidyr)
})

set.seed(2026 - 9 - 24)
n_bootstrap <- 2000
n_posterior <- 4000

draws_dir <- "outputs/cv_draws"

# the three out-of-sample experiments worth reporting. Leave-one-country-out is
# excluded: it measures the difficulty of an entirely unsampled country, which
# the deployed model never faces, and its held-out bias tracks the country's own
# fitted effect at r = -0.94 (#12 review).
experiments <- list(
  list(label = "spatial interpolation", experiment = "spatial_interpolation",
       folds = "all"),
  list(label = "spatial extrapolation", experiment = "spatial_blocks",
       folds = c("1", "2")),
  # the two five-year rolling origins, pooled into one bar as the two spatial
  # blocks are. The 2020 three-year fold they replace both leaked training
  # records into the holdout and landed on the one pause in twenty years of
  # decline, which was also the thinnest window.
  list(label = "temporal change", experiment = "temporal_forecasting",
       folds = c("2014", "2018"))
)

# each contrast is the first model less the second, in percentage points of the
# observed variance
contrasts <- list(
  c("dynamical", "nearest_neighbour"),
  c("dynamical", "nearest_neighbour_oracle"),
  c("dynamical", "intercept"),
  c("nearest_neighbour", "intercept"))

models <- c(dynamical = "dynamical model",
            nearest_neighbour = "nearest recent survey",
            nearest_neighbour_oracle = "nearest surveys, best k",
            intercept = "insecticide mean")

# the same per-type overdispersion the rest of the pipeline scores against.
# rho_table, not rho_lookup: the latter is the function that supplies it.
rho_spec <- rho_lookup()
rho_key <- rho_spec$key
rho_table <- rho_spec$table %>%
  mutate(se_logit = (qlogis(rho_upper) - qlogis(rho_lower)) / (2 * 1.96))
cat("overdispersion used for the noise share:", rho_spec$source, "\n")


# assemble the held-out records, with one prediction column per model ---------

read_fold <- function(model, experiment, fold) {
  file <- file.path(draws_dir, sprintf("%s__%s__%s.rds", model, experiment, fold))
  if (!file.exists(file)) return(NULL)
  x <- readRDS(file)
  stopifnot(ncol(x$p_draws) == nrow(x$test_df))
  data.frame(row = seq_len(nrow(x$test_df)), model = model,
             predicted = colMeans(x$p_draws))
}

records <- bind_rows(lapply(experiments, function(spec) {
  bind_rows(lapply(spec$folds, function(fold) {

    reference <- readRDS(file.path(
      draws_dir, sprintf("dynamical__%s__%s.rds", spec$experiment, fold)))$test_df

    predictions <- bind_rows(lapply(names(models), read_fold,
                                    experiment = spec$experiment, fold = fold))
    wide <- predictions %>%
      pivot_wider(names_from = model, values_from = predicted,
                  names_prefix = "p_")
    stopifnot(nrow(wide) == nrow(reference))

    reference %>%
      transmute(experiment = spec$label, fold = fold,
                cell, insecticide_type, insecticide_class,
                died, mosquito_number,
                observed = died / mosquito_number) %>%
      bind_cols(wide %>% select(starts_with("p_")))
  }))
}))

# every model must have a prediction for every record, or a bar would be drawn
# from a different denominator than its neighbours
prediction_columns <- paste0("p_", names(models))
stopifnot(all(prediction_columns %in% names(records)),
          !any(is.na(records[, prediction_columns])))

records <- records %>%
  left_join(rho_table %>% select(key, rho_type = rho),
            by = setNames("key", rho_key))
stopifnot(!any(is.na(records$rho_type)))

cat(sprintf("%i held-out records across %i experiments\n",
            nrow(records), n_distinct(records$experiment)))


# the two quantities ---------------------------------------------------------

# beta-binomial design effect over the assay size: the variance of the observed
# proportion is p(1-p) k
design_effect <- function(m, rho) (1 + (m - 1) * rho) / m

# marginal likelihood of a Beta(a, b) prior on the true proportion, at the
# effective sample size the overdispersion leaves each assay
beta_prior_nll <- function(par, observed, m_eff) {
  a <- exp(par[1])
  b <- exp(par[2])
  successes <- observed * m_eff
  -sum(lbeta(a + successes, b + m_eff - successes) - lbeta(a, b))
}

fit_beta_priors <- function(data) {
  bind_rows(lapply(split(data, data[[rho_key]]), function(d) {
    m_eff <- 1 / design_effect(d$mosquito_number, d$rho_type)
    fit <- optim(c(0, 0), beta_prior_nll, observed = d$observed, m_eff = m_eff,
                 method = "Nelder-Mead", control = list(reltol = 1e-10))
    data.frame(key = d[[rho_key]][1], a = exp(fit$par[1]), b = exp(fit$par[2]),
               prior_mean = exp(fit$par[1]) / sum(exp(fit$par)),
               converged = fit$convergence == 0, records = nrow(d))
  }))
}

# per-record floor contributions at a given vector of rho values. The Beta
# prior is held at its fit under the point estimate of rho while rho varies for
# the interval: refitting it inside every posterior draw would cost thousands of
# optimisations for a second-order effect on the prior.
floor_contributions <- function(data, rho_vector) {
  k <- design_effect(data$mosquito_number, rho_vector)
  usable <- data$mosquito_number > 1 & is.finite(k) & k > 0 & k < 1
  out <- rep(NA_real_, nrow(data))
  m_eff <- 1 / k[usable]
  shape_a <- data$a[usable] + data$observed[usable] * m_eff
  shape_b <- data$b[usable] + (1 - data$observed[usable]) * m_eff
  out[usable] <- shape_a * shape_b /
    ((shape_a + shape_b) * (shape_a + shape_b + 1)) * k[usable]
  out
}

# fitted once, per insecticide type, over all of that type's held-out records
beta_priors <- fit_beta_priors(records)
stopifnot(all(beta_priors$converged))
records <- records %>%
  left_join(beta_priors %>% select(key, a, b), by = setNames("key", rho_key))
stopifnot(!any(is.na(records$a)), !any(is.na(records$b)))

cat("\nBeta prior on the true proportion, fitted per insecticide type:\n")
print(as.data.frame(beta_priors %>%
  transmute(type = key, records, a = round(a, 2), b = round(b, 2),
            prior_mean = round(prior_mean, 3))), row.names = FALSE)

explained_for <- function(data) {
  variance <- mean((data$observed - mean(data$observed)) ^ 2)
  vapply(prediction_columns, function(column) {
    100 * (1 - mean((data$observed - data[[column]]) ^ 2) / variance)
  }, numeric(1))
}

# bootstrap the explained fraction by resampling pixels
bootstrap_explained <- function(data) {
  cells <- unique(data$cell)
  rows_by_cell <- split(seq_len(nrow(data)), data$cell)
  replicates <- t(replicate(n_bootstrap, {
    picked <- sample(cells, length(cells), replace = TRUE)
    explained_for(data[unlist(rows_by_cell[as.character(picked)],
                              use.names = FALSE), ])
  }))
  colnames(replicates) <- prediction_columns
  replicates
}

# posterior draws of the noise share, propagating rho only. The logit-scale
# standard errors come from the fit that produced rho, so this carries the
# uncertainty in the overdispersion but treats Var(y) as known - that is where
# essentially all of the floor's uncertainty sits.
noise_posterior <- function(data) {
  variance <- mean((data$observed - mean(data$observed)) ^ 2)
  index <- match(data[[rho_key]], rho_table$key)
  centre <- qlogis(rho_table$rho)
  spread <- rho_table$se_logit
  vapply(seq_len(n_posterior), function(i) {
    drawn <- plogis(rnorm(nrow(rho_table), centre, spread))
    100 * mean(floor_contributions(data, drawn[index]), na.rm = TRUE) / variance
  }, numeric(1))
}

summarise_subset <- function(data, label, stratum = NA_character_) {

  variance <- mean((data$observed - mean(data$observed)) ^ 2)
  point <- explained_for(data)
  replicates <- bootstrap_explained(data)
  noise <- noise_posterior(data)
  noise_point <- 100 * mean(floor_contributions(data, data$rho_type),
                            na.rm = TRUE) / variance

  bind_rows(
    bind_rows(lapply(prediction_columns, function(column) {
      data.frame(
        quantity = models[sub("^p_", "", column)],
        kind = "model",
        estimate = point[[column]],
        lower = quantile(replicates[, column], 0.025, na.rm = TRUE),
        upper = quantile(replicates[, column], 0.975, na.rm = TRUE))
    })),
    # Paired contrasts, formed within each bootstrap replicate. Every model is
    # scored on the same resample, so the per-model intervals are strongly
    # correlated and whether two of them overlap says nothing about the
    # difference; this is the only thing that can settle a comparison. The
    # denominator cancels within a replicate, so the contrast is just
    # 100 (MSE_reference - MSE_model) / Var(y).
    bind_rows(lapply(contrasts, function(pair) {
      a <- paste0("p_", pair[1])
      b <- paste0("p_", pair[2])
      drawn <- replicates[, a] - replicates[, b]
      data.frame(quantity = paste(models[[pair[1]]], "-", models[[pair[2]]]),
                 kind = "contrast",
                 estimate = point[[a]] - point[[b]],
                 lower = quantile(drawn, 0.025, na.rm = TRUE),
                 upper = quantile(drawn, 0.975, na.rm = TRUE))
    })),
    data.frame(quantity = "bioassay variability", kind = "noise",
               estimate = noise_point,
               lower = quantile(noise, 0.025), upper = quantile(noise, 0.975))
  ) %>%
    mutate(experiment = label, stratum = stratum,
           assays = nrow(data), pixels = n_distinct(data$cell),
           variance = variance, .before = everything())
}

pooled <- bind_rows(lapply(split(records, records$experiment), function(data) {
  summarise_subset(data, data$experiment[1])
}))
row.names(pooled) <- NULL
write.csv(pooled, "outputs/cv_variance_explained.csv", row.names = FALSE)

# and per fold, for experiments built from more than one, so that a pooled bar
# can be checked against the folds it pools
by_fold <- bind_rows(lapply(
  split(records, paste(records$experiment, records$fold)), function(data) {
    if (n_distinct(records$fold[records$experiment == data$experiment[1]]) < 2) {
      return(NULL)
    }
    summarise_subset(data, data$experiment[1], data$fold[1])
  }))
row.names(by_fold) <- NULL
write.csv(by_fold, "outputs/cv_variance_explained_by_fold.csv", row.names = FALSE)

cat("\nper experiment, % of observed variance in held-out mortality:\n")
print(as.data.frame(pooled %>%
  mutate(value = sprintf("%5.1f [%5.1f, %5.1f]", estimate, lower, upper)) %>%
  select(experiment, assays, pixels, quantity, value) %>%
  pivot_wider(names_from = quantity, values_from = value)), row.names = FALSE)


cat("\npaired contrasts, within bootstrap replicates:\n")
print(as.data.frame(pooled %>%
  filter(kind == "contrast") %>%
  mutate(value = sprintf("%6.1f [%6.1f, %6.1f]", estimate, lower, upper)) %>%
  select(experiment, quantity, value) %>%
  pivot_wider(names_from = quantity, values_from = value)), row.names = FALSE)

cat("\nper fold within experiment:\n")
print(as.data.frame(by_fold %>%
  filter(kind != "contrast") %>%
  mutate(value = sprintf("%5.1f [%5.1f, %5.1f]", estimate, lower, upper)) %>%
  select(experiment, fold = stratum, assays, pixels, quantity, value) %>%
  pivot_wider(names_from = quantity, values_from = value)), row.names = FALSE)

cat("\npaired contrasts per fold:\n")
print(as.data.frame(by_fold %>%
  filter(kind == "contrast") %>%
  mutate(value = sprintf("%6.1f [%6.1f, %6.1f]", estimate, lower, upper)) %>%
  select(experiment, fold = stratum, quantity, value) %>%
  pivot_wider(names_from = quantity, values_from = value)), row.names = FALSE)


# and by insecticide ---------------------------------------------------------

# Every insecticide with held-out records is scored, and nothing is dropped
# from the table. What the figure shows is decided afterwards, by whether the
# cell can carry a reading at all, on two counts that are not substitutes:
#
#   ceiling >= min_ceiling
#       more than half the observed spread in held-out mortality has to be real
#       variation in resistance rather than assay noise. Where it is not, the
#       ratio is two small numbers divided by each other: Bendiocarb in the
#       interpolation fold has 154 assays in 67 pixels and a predictable
#       standard deviation of about 6 points of mortality, so a model would
#       have to land inside 6 points to score above zero. A sample-size rule
#       cannot see this, which is why the 20-assay threshold this replaces let
#       that cell through at -158%.
#
#   interval width <= max_ci_width
#       the estimate has to be resolvable, for both plotted models, or the
#       comparison the panel exists to make cannot be read. This is where the
#       pixel count enters, but only through its effect on precision:
#       Lambda-cyhalothrin in the forecast fold has 25 points of predictable
#       variation and 18 pixels to estimate it from.
min_ceiling <- 50
max_ci_width <- 100
shown_models <- c("dynamical model", "nearest recent survey")

by_insecticide <- bind_rows(lapply(
  split(records, paste(records$experiment, records$insecticide_type)),
  function(data) {
    # a single pixel leaves no between-cluster variation to bootstrap, and a
    # constant holdout leaves no denominator
    if (n_distinct(data$cell) < 2) return(NULL)
    if (mean((data$observed - mean(data$observed)) ^ 2) == 0) return(NULL)
    summarise_subset(data, data$experiment[1], data$insecticide_type[1])
  }))

by_insecticide <- by_insecticide %>%
  group_by(experiment, stratum) %>%
  mutate(
    ceiling = 100 - estimate[kind == "noise"],
    ci_width = max((upper - lower)[quantity %in% shown_models]),
    drop_reason = case_when(
      ceiling < min_ceiling & ci_width > max_ci_width ~ "no signal, imprecise",
      ceiling < min_ceiling                           ~ "no signal",
      ci_width > max_ci_width                         ~ "imprecise",
      TRUE                                            ~ NA_character_),
    shown = is.na(drop_reason)) %>%
  ungroup()
row.names(by_insecticide) <- NULL
write.csv(by_insecticide,
          "outputs/cv_variance_explained_by_insecticide.csv", row.names = FALSE)

cells <- by_insecticide %>% distinct(experiment, stratum, assays, pixels,
                                     variance, ceiling, ci_width, shown,
                                     drop_reason)
cat(sprintf("\nby insecticide: %i experiment-insecticide combinations, %i shown\n",
            nrow(cells), sum(cells$shown)))
cat(sprintf("thresholds: ceiling >= %i%%, 95%% interval width <= %i points\n",
            min_ceiling, max_ci_width))
print(as.data.frame(cells %>%
  transmute(experiment, insecticide = stratum, assays, pixels,
            sd_pp = round(100 * sqrt(variance), 1),
            ceiling = round(ceiling, 1),
            predictable_sd_pp = round(100 * sqrt(variance) *
                                        sqrt(pmax(ceiling, 0) / 100), 1),
            ci_width = round(ci_width, 1),
            shown = ifelse(shown, "yes", drop_reason)) %>%
  arrange(desc(shown == "yes"), experiment, insecticide)), row.names = FALSE)
