# Effective number of parameters of the stage-one dynamical model, and a
# grouped PSIS check of how much its training-assay predictions have already
# absorbed the data they are predicting (doc/two_stage_plan.md, step 6).
#
#   Rscript R/stage_one_effective_parameters.R <experiment> <fold>
#
# e.g. Rscript R/stage_one_effective_parameters.R spatial_interpolation all
#      Rscript R/stage_one_effective_parameters.R spatial_blocks 1
#      Rscript R/stage_one_effective_parameters.R temporal_forecasting 2014
#
# Why this is needed. The two-stage correction fits a residual model to the
# empirical logits of the training assays, offset by the dynamical model's
# posterior mean logit prediction there (m_ref). But the dynamical model was fit
# to those same assays, so m_ref has partly been pulled towards them, and the
# residuals z - m_ref understate the residual variance a new assay would show.
# How much depends on how flexible stage one is, locally: it has 689 nominal
# parameters, most of them shrunk hierarchically, so the nominal count says
# little. Two complementary measures are computed here from the model's own
# paired posterior draws of p and rho at the training assays:
#
#   1. Global effective parameter counts: pD (Spiegelhalter et al. 2002), pV
#      (Gelman et al.), p_WAIC and p_loo. These say how many degrees of freedom
#      the fit uses, overall and per insecticide type.
#   2. Grouped PSIS leave-out: the leave-one-group-out posterior mean of logit p
#      at each training assay, where the group is the whole pixel-year (most
#      assays share a cell-year with others of the same type, so leaving out a
#      single assay would still leave its twins in). The shift from m_ref,
#      relative to the residual SD, says directly how much the stage-one fit has
#      absorbed locally. A harsher variant leaves out the whole pixel (all
#      years).
#
# The per-assay leave-out means are saved so the stage-A fit can use them as an
# alternative m_ref, a cheap approximation to stacking stage one (refitting it
# K times is unaffordable at 60-95 h per MCMC run).
#
# Memory: loading a fold takes 4-8 GB and the draws x assays matrices ~0.5 GB
# each, so run one fold per process.

arguments <- commandArgs(trailingOnly = TRUE)
experiment_name <- arguments[1]
fold_name <- arguments[2]
stopifnot(!is.na(experiment_name), !is.na(fold_name))

start_time <- Sys.time()
report <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")
  flush(stdout())
}

suppressMessages({
  sink("/dev/null")
  source("R/validation_folds.R")
  source("R/validation_covariates.R")
  sink()
})
source("R/validation_functions.R")
source("R/dynamical_predictions.R")
suppressMessages(library(loo))

n_cores <- 4
output_dir <- "outputs/two_stage"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# the training set of the requested fold, found the same way run_one_fold.R
# found it when the fold was fitted
if (experiment_name == "spatial_interpolation") {
  training <- spatial_interpolation$training
  test <- spatial_interpolation$test
} else if (experiment_name == "spatial_blocks") {
  suppressMessages({
    sink("/dev/null")
    source("R/validation_blocks.R")
    sink()
  })
  index <- as.integer(fold_name)
  stopifnot(!is.na(index), index >= 1, index <= length(spatial_blocks))
  training <- spatial_blocks[[index]]$training
  test <- spatial_blocks[[index]]$test
} else if (experiment_name == "temporal_forecasting") {
  stopifnot(fold_name %in% names(temporal_forecasting_folds))
  training <- temporal_forecasting_folds[[fold_name]]$training
  test <- temporal_forecasting_folds[[fold_name]]$test
} else if (experiment_name == "spatial_extrapolation") {
  index <- match(fold_name, countries_to_validate)
  stopifnot(!is.na(index))
  training <- spatial_extrapolation$training[[index]]
  test <- spatial_extrapolation$test[[index]]
} else {
  stop("unknown experiment: ", experiment_name)
}

fold_file <- sprintf("outputs/cv_draws/dynamical__%s__%s.rds",
                     experiment_name, fold_name)
report("loading %s", fold_file)
fold <- readRDS(fold_file)

# The fold does not store its training set, but it does store its test set; if
# that matches the one rebuilt here, the fold definitions have not drifted since
# fitting and the training set rebuilt alongside it is the one the model saw.
stopifnot(identical(test$cell_id, fold$test_df$cell_id),
          identical(test$type_id, fold$test_df$type_id),
          identical(test$year_id, fold$test_df$year_id),
          identical(classes_index[training$type_id], training$class_id))

# don't need the stored held-out predictions, and they are large
fold$p_draws <- NULL
fold$p_draws_before <- NULL
invisible(gc())

draw_index <- paired_draw_index(fold)
n_draws <- length(draw_index)

# Chain of each paired draw. as.matrix() stacks the chains in order and the
# thinning keeps the order within a chain, so the thinned draws are still
# autocorrelated series per chain and loo's relative efficiencies can be
# estimated properly rather than assuming r_eff = 1.
chain_sizes <- vapply(fold$draws, nrow, integer(1))
chain_id <- findInterval(draw_index - 1, cumsum(chain_sizes)) + 1

draws_matrix <- as.matrix(fold$draws)[draw_index, , drop = FALSE]
rho_draws <- extract_parameter(draws_matrix, "rho_classes")
dim(rho_draws) <- c(n_draws, length(classes))
# these are the rho draws the held-out scoring paired with p (older folds
# stored them unthinned, so they are thinned by the same rule first)
if (!is.null(fold$rho_class_draws)) {
  stored_rho <- fold$rho_class_draws
  if (nrow(stored_rho) > n_draws) {
    stored_rho <- stored_rho[draw_index, , drop = FALSE]
  }
  stopifnot(isTRUE(all.equal(unname(rho_draws), unname(stored_rho))))
  rm(stored_rho)
}

# nominal parameter count: every named sampled quantity (689), plus the
# n_types elements of logit_init_mean, which were sampled but not named (see
# logit_init_mean_draws())
n_nominal <- ncol(draws_matrix) + length(types)

# Repeated draws. Some chains stick at a state for longer than the thinning
# interval (HMC rejections in a row), so a few percent of the paired draws are
# exact copies of the previous one. That is a valid, if inefficient, sample, but
# it breaks PSIS: the Pareto fit to the largest importance ratios sees tied
# values at its cutoff and fails, and loo reports k = Inf (or a spuriously large
# k once the ties are jittered apart), for groups whose weights are in fact
# well behaved. So the importance-sampling steps run on the distinct draws only
# and are expressed as shifts from the distinct-draw mean (see grouped_psis()),
# while the moment-based counts (pD, pV, p_WAIC) use every draw.
distinct_draw <- !duplicated(draws_matrix)

# data at the training assays
y <- training$died
n <- training$mosquito_number
class_id <- training$class_id
type_id <- training$type_id
n_train <- nrow(training)
# empirical logit, the stage-A response
z <- log((y + 0.5) / (n - y + 0.5))


# stage-one draws at the training assays -----------------------------------

report("predicting at %d training assays, %d draws", n_train, n_draws)
p_train <- dynamical_predictions(fold, training, df, x_cell_years,
                                 cell_years_index, classes_index, types,
                                 draw_index = draw_index)

# Logit p. p is computed as plogis() of a logit, which rounds to exactly 1 when
# the logit exceeds ~37, so it is clamped before qlogis(); this is the same
# clamp dbetabinom() applies in the likelihood, so the two agree.
clamp <- 1e-12
n_clamped <- sum(p_train > 1 - clamp | p_train < clamp)
logit_p <- qlogis(pmin(pmax(p_train, clamp), 1 - clamp))
m_ref <- colMeans(logit_p)
post_sd_logit <- matrixStats::colSds(logit_p)
p_bar <- colMeans(p_train)


# pointwise log-likelihood -------------------------------------------------

# beta-binomial log-likelihood of every training assay under every draw, in
# column chunks so the recycled data vectors stay small
pointwise_loglik <- function(p, rho, chunk = 2000) {
  out <- matrix(NA_real_, nrow(p), ncol(p))
  for (columns in split(seq_len(ncol(p)), ceiling(seq_len(ncol(p)) / chunk))) {
    out[, columns] <- dbetabinom(
      y = rep(y[columns], each = nrow(p)),
      size = rep(n[columns], each = nrow(p)),
      p = p[, columns, drop = FALSE],
      rho = rho[, class_id[columns], drop = FALSE],
      log = TRUE)
  }
  out
}

report("pointwise log-likelihood")
loglik <- pointwise_loglik(p_train, rho_draws)
stopifnot(all(is.finite(loglik)))
rm(p_train)
invisible(gc())


# plug-in deviances for pD ---------------------------------------------------

rho_bar <- matrix(colMeans(rho_draws), nrow = 1)

# (a) at the posterior mean of the parameters. A one-draw pseudo fold holding
# the parameter means is pushed through the same recursion. The means are of
# the quantities that were sampled (the named constrained values, and the raw
# free state for logit_init_mean, which is untransformed so its raw mean is its
# mean); the pseudo fold carries a model_info whose raw draws are the means of
# the raw draws, so logit_init_mean_draws() finds the block the usual way. Its
# identity check on beta_overall holds too, since raw and named agree drawwise.
model_info <- attr(fold$draws, "model_info")
raw <- do.call(rbind, lapply(model_info$raw_draws, as.matrix))[draw_index, ,
                                                               drop = FALSE]
mean_info <- model_info
mean_info$raw_draws <- list(coda::mcmc(matrix(colMeans(raw), nrow = 1,
                                              dimnames = list(NULL,
                                                              colnames(raw)))))
rm(raw)
mean_draws <- coda::mcmc.list(coda::mcmc(
  matrix(colMeans(draws_matrix), nrow = 1,
         dimnames = list(NULL, colnames(draws_matrix)))))
attr(mean_draws, "model_info") <- mean_info
mean_fold <- list(draws = mean_draws)
p_at_mean_theta <- dynamical_predictions(mean_fold, training, df, x_cell_years,
                                         cell_years_index, classes_index,
                                         types, draw_index = 1)

# the fold's draws are no longer needed; only the matrices derived from them
rm(fold, mean_fold, mean_draws, mean_info, model_info, draws_matrix)
invisible(gc())

loglik_mean_theta <- pointwise_loglik(p_at_mean_theta, rho_bar)
# (b) at the posterior mean of p, the cheaper, conventional-in-practice choice
loglik_mean_p <- pointwise_loglik(matrix(p_bar, nrow = 1), rho_bar)
# (c) at the posterior mean of logit p, i.e. at m_ref itself
loglik_mean_logit <- pointwise_loglik(matrix(plogis(m_ref), nrow = 1), rho_bar)


# PSIS-LOO, pointwise ----------------------------------------------------------

# relative efficiency of exp(log-lik); the column max is subtracted first so
# large groups cannot underflow, which rescales and so leaves ESS unchanged
relative_efficiency <- function(log_lik) {
  scaled <- exp(sweep(log_lik, 2, matrixStats::colMaxs(log_lik)))
  loo::relative_eff(scaled, chain_id = chain_id, cores = n_cores)
}

report("PSIS-LOO over single assays")
r_eff <- relative_efficiency(loglik)
# on the distinct draws (see above); r_eff comes from the full chains, since
# relative_eff() needs equal-length chains and the ties are part of what makes
# the chains inefficient
loo_fit <- loo::loo(loglik[distinct_draw, ], r_eff = r_eff, cores = n_cores)
pointwise_p_loo <- loo_fit$pointwise[, "p_loo"]
pointwise_k_single <- loo_fit$diagnostics$pareto_k
rm(loo_fit)


# effective parameters, overall and per type ----------------------------------

pointwise_var <- matrixStats::colVars(loglik)

effective_parameters <- function(columns, type_label) {
  deviance <- -2 * rowSums(loglik[, columns, drop = FALSE])
  mean_deviance <- mean(deviance)
  tibble(
    experiment = experiment_name,
    fold = fold_name,
    insecticide_type = type_label,
    n_assays = length(columns),
    n_pixel_years = n_distinct(training$cell_id[columns],
                               training$year_id[columns],
                               training$type_id[columns]),
    n_pixels = n_distinct(training$cell_id[columns],
                          training$type_id[columns]),
    n_nominal = n_nominal,
    mean_deviance = mean_deviance,
    pD_mean_theta = mean_deviance - -2 * sum(loglik_mean_theta[, columns]),
    pD_mean_p = mean_deviance - -2 * sum(loglik_mean_p[, columns]),
    pD_mean_logit = mean_deviance - -2 * sum(loglik_mean_logit[, columns]),
    pV = var(deviance) / 2,
    p_waic = sum(pointwise_var[columns]),
    p_loo = sum(pointwise_p_loo[columns]),
    elpd_loo_share_k_gt_0.7 = mean(pointwise_k_single[columns] > 0.7),
    r_eff_median = median(r_eff[columns])
  )
}

report("effective parameter counts")
parameter_table <- bind_rows(
  effective_parameters(seq_len(n_train), "all"),
  lapply(sort(unique(type_id)), function(k) {
    effective_parameters(which(type_id == k), types[k])
  })
)
print(parameter_table, width = Inf)


# grouped PSIS leave-out -------------------------------------------------------

# Leave out a whole group of assays: the log importance ratio for draw s is
# minus the group's summed log-likelihood, smoothed by PSIS. The leave-out mean
# of logit p at each assay is then the weighted mean over draws with its own
# group's weights. Where k > 0.7 the smoothed weights are unreliable and the
# plan falls back to m_ref, so `m_loo_fallback` does that.
grouped_psis <- function(group) {
  group_index <- match(group, unique(group))
  n_groups <- max(group_index)
  # draws x groups summed log-likelihood
  group_loglik <- t(rowsum(t(loglik), group_index, reorder = TRUE))
  group_r_eff <- relative_efficiency(group_loglik)
  psis_fit <- loo::psis(-group_loglik[distinct_draw, , drop = FALSE],
                        r_eff = group_r_eff, cores = n_cores)
  k <- loo::pareto_k_values(psis_fit)
  weights <- weights(psis_fit, log = FALSE, normalize = TRUE)
  rm(psis_fit, group_loglik)

  # the leave-out shift is taken against the unweighted mean over the same
  # distinct draws, then added to m_ref, so dropping the repeats does not itself
  # register as a shift
  m_loo <- numeric(n_train)
  for (columns in split(seq_len(n_train),
                        ceiling(seq_len(n_train) / 2000))) {
    logit_distinct <- logit_p[distinct_draw, columns, drop = FALSE]
    m_loo[columns] <- m_ref[columns] +
      colSums(weights[, group_index[columns], drop = FALSE] *
                logit_distinct) -
      colMeans(logit_distinct)
  }
  tibble(group_id = group_index,
         group_size = tabulate(group_index)[group_index],
         k = k[group_index],
         m_loo = m_loo,
         m_loo_fallback = ifelse(k > 0.7, m_ref, m_loo))
}

report("grouped PSIS, pixel-year groups")
by_pixel_year <- grouped_psis(paste(training$cell_id, training$year_id,
                                    training$type_id))
report("grouped PSIS, pixel groups")
by_pixel <- grouped_psis(paste(training$cell_id, training$type_id))

# per-assay table, for use as an alternative m_ref in stage A
assay_table <- training %>%
  select(any_of(c("cell_id", "country_id", "type_id", "insecticide_type",
                  "class_id", "year_id", "year_start", "died",
                  "mosquito_number"))) %>%
  mutate(training_row = row_number(),
         z = z,
         m_ref = m_ref,
         post_sd_logit = post_sd_logit,
         k_single = pointwise_k_single,
         group_pixel_year = by_pixel_year$group_id,
         group_size_pixel_year = by_pixel_year$group_size,
         k_pixel_year = by_pixel_year$k,
         m_loo_pixel_year = by_pixel_year$m_loo,
         m_loo_pixel_year_fallback = by_pixel_year$m_loo_fallback,
         group_pixel = by_pixel$group_id,
         k_pixel = by_pixel$k,
         m_loo_pixel = by_pixel$m_loo,
         m_loo_pixel_fallback = by_pixel$m_loo_fallback)
saveRDS(assay_table, file.path(output_dir,
                               sprintf("loo_m__%s__%s.rds",
                                       experiment_name, fold_name)))

# Summaries. The shift is judged against the residual SD of z - m_ref within
# type, the scale on which stage A works. The slope of the shift on the
# residual is a local "hat value": leaving a group out moves its prediction
# away from its data, so the slope is negative, and -slope is the share of the
# residual the stage-one fit had absorbed. The ratio of the leave-out residual
# SD to the in-sample one says how much stage A's residual variance is
# understated by using m_ref.
psis_summary <- function(result, columns, type_label, grouping) {
  residual <- z[columns] - m_ref[columns]
  residual_sd <- sd(residual)
  shift <- result$m_loo[columns] - m_ref[columns]
  good <- result$k[columns] <= 0.7
  groups <- !duplicated(result$group_id[columns])
  k_groups <- result$k[columns][groups]
  tibble(
    experiment = experiment_name,
    fold = fold_name,
    insecticide_type = type_label,
    grouping = grouping,
    n_assays = length(columns),
    n_groups = sum(groups),
    median_group_size = median(result$group_size[columns][groups]),
    k_median = median(k_groups),
    k_max = max(k_groups),
    share_groups_k_gt_0.5 = mean(k_groups > 0.5),
    share_groups_k_gt_0.7 = mean(k_groups > 0.7),
    share_assays_k_gt_0.7 = mean(!good),
    residual_sd = residual_sd,
    post_sd_rel = sqrt(mean(post_sd_logit[columns]^2)) / residual_sd,
    shift_mean_rel = mean(shift) / residual_sd,
    shift_rms_rel = sqrt(mean(shift^2)) / residual_sd,
    shift_rms_rel_k_ok = sqrt(mean(shift[good]^2)) / residual_sd,
    shift_p95_abs_rel = quantile(abs(shift), 0.95, names = FALSE) /
      residual_sd,
    shift_residual_slope = unname(coef(lm(shift ~ residual))[2]),
    loo_residual_sd_ratio =
      sd(z[columns] - result$m_loo_fallback[columns]) / residual_sd
  )
}

type_columns <- c(list(all = seq_len(n_train)),
                  setNames(lapply(sort(unique(type_id)),
                                  function(k) which(type_id == k)),
                           types[sort(unique(type_id))]))
summary_table <- bind_rows(
  lapply(names(type_columns), function(label) {
    bind_rows(psis_summary(by_pixel_year, type_columns[[label]], label,
                           "pixel_year"),
              psis_summary(by_pixel, type_columns[[label]], label, "pixel"))
  })
)
print(summary_table, width = Inf)


# write, replacing this fold's rows in the running tables ---------------------

runtime_minutes <- as.numeric(difftime(Sys.time(), start_time,
                                       units = "mins"))
parameter_table <- parameter_table %>%
  mutate(n_draws = n_draws,
         n_distinct_draws = sum(distinct_draw),
         n_logit_clamped = n_clamped,
         runtime_minutes = runtime_minutes)

upsert_csv <- function(new_rows, file) {
  path <- file.path(output_dir, file)
  if (file.exists(path)) {
    old <- read_csv(path, show_col_types = FALSE,
                    col_types = cols(fold = col_character()))
    new_rows <- bind_rows(
      old %>% filter(!(experiment == experiment_name & fold == fold_name)),
      new_rows)
  }
  write_csv(new_rows, path)
}
upsert_csv(parameter_table, "stage_one_effective_parameters.csv")
upsert_csv(summary_table, "stage_one_psis_summary.csv")

report("done in %.1f minutes", runtime_minutes)
