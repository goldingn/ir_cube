# Score the two-stage correction models against the dynamical model and the
# null models on the five cross-validation folds of #12 (issue #21).
#
#   Rscript R/two_stage_metrics.R
#
# Reads the #12 draws (outputs/cv_draws) and the two-stage draws
# (outputs/cv_draws_two_stage), both in the #12 fold format, and writes to
# outputs/two_stage/ and figures/two_stage/. Safe to run while the two-stage
# folds are still being produced: whatever exists is scored, and anything
# missing is skipped with a message.
#
# The scoring mirrors R/validation_metrics.R and R/variance_explained.R, with
# two deliberate differences:
#
#   1. Every model, the dynamical model and the nulls included, is scored at the
#      per-insecticide-type overdispersion from the hierarchical replicate fit
#      (outputs/bioassay_rho_hierarchical.csv), not the per-class estimate #12
#      uses. The two-stage models were fitted with the per-type values, so
#      scoring the others at anything else would make the comparison partly a
#      comparison of observation models. The dynamical and null numbers here
#      therefore differ slightly from outputs/cv_summary.csv.
#   2. Every difference is reported against the dynamical model, since the
#      question is whether the correction improves on it, and every difference
#      carries a paired pixel-cluster bootstrap interval. Skill is still anchored
#      on the intercept null, as in #12.
#
# validation_metrics.R is not sourced, because sourcing it scores every fold
# and writes #12's outputs. The pieces needed from it (thinning, the per-fold
# scoring, the pixel bootstrap) and from variance_explained.R (the regularised
# noise floor) are copied and adapted below; the primitives come from
# R/validation_functions.R unchanged.
#
# Only the five #12 folds are scored: spatial_interpolation/all,
# spatial_blocks/1 and 2, temporal_forecasting/2014 and 2018. Leave-one-country
# -out and the legacy 2020 forecasting fold are defunct (doc/two_stage_plan.md).

source("R/validation_functions.R")

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

start_time <- Sys.time()
report <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")
  flush(stdout())
}

draws_dir <- "outputs/cv_draws"
# overridable so the script can be exercised on test draws
two_stage_dir <- Sys.getenv("TWO_STAGE_DRAWS_DIR", "outputs/cv_draws_two_stage")
output_dir <- "outputs/two_stage"
figure_dir <- "figures/two_stage"
cache_dir <- file.path(output_dir, "score_cache")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

n_pit_reps <- 100
n_bootstrap <- 2000
coverage_levels <- seq(0.05, 0.95, by = 0.05)

# same thinning as #12: at roughly 50 draws per effective sample, every tenth of
# the 20,000 MCMC draws keeps almost all the information, and the scoring cost
# is linear in the draws
max_draws <- 2000
thin_draws <- function(x, maximum = max_draws) {
  if (nrow(x) <= maximum) return(x)
  keep <- round(seq(1, nrow(x), length.out = maximum))
  x[keep, , drop = FALSE]
}

# Bumped whenever the per-fold scoring changes, so stale cached scores are not
# reused. The cache exists because a dynamical fold is 1.2-3.5 GB on disk and
# 4-8 GB in memory; rescoring it every time a new two-stage fold lands would be
# minutes of I/O and a lot of shared RAM for an identical answer.
scoring_version <- "2026-09-25.1"

# The folds of #12, and the label each is pooled under. The two block folds are
# one experiment; the two forecasting origins are kept apart because they
# forecast different windows from different amounts of training data.
fold_specs <- tribble(
  ~experiment,             ~fold,  ~label,
  "spatial_interpolation", "all",  "spatial_interpolation",
  "spatial_blocks",        "1",    "spatial_blocks",
  "spatial_blocks",        "2",    "spatial_blocks",
  "temporal_forecasting",  "2014", "temporal_forecasting_2014",
  "temporal_forecasting",  "2018", "temporal_forecasting_2018"
)

# A fold is only scored when the dynamical model has been fitted to it, since
# every comparison is paired against it and the two-stage models cannot exist
# without it. Set TWO_STAGE_ALLOW_NO_DYNAMICAL=1 to score the nulls on such
# folds anyway (used to exercise the forecasting code before those fits land).
allow_no_dynamical <- Sys.getenv("TWO_STAGE_ALLOW_NO_DYNAMICAL") == "1"

null_models <- c("intercept", "nearest_neighbour", "nearest_neighbour_oracle")

# display names and colours, fixed per model so a colour never changes meaning
# between figures. The #12 models keep the colours of fig_variance_bars.R; the
# two-stage variants take warm hues, lighter for the leave-out-m_ref versions
model_labels <- c(
  dynamical = "dynamical model",
  two_stage_omega_u = "two-stage: ω + u",
  two_stage_omega_xi_u = "two-stage: ω + ξ + u",
  two_stage_omega_u_loo = "two-stage: ω + u (leave-out m)",
  two_stage_omega_xi_u_loo = "two-stage: ω + ξ + u (leave-out m)",
  nearest_neighbour_oracle = "nearest surveys, best k",
  nearest_neighbour = "nearest recent survey",
  intercept = "insecticide mean"
)
model_colours <- c(
  dynamical = "#2166AC",
  two_stage_omega_u = "#E08214",
  two_stage_omega_xi_u = "#B2182B",
  two_stage_omega_u_loo = "#FDB863",
  two_stage_omega_xi_u_loo = "#F4A582",
  nearest_neighbour_oracle = "#8073AC",
  nearest_neighbour = "#1B7837",
  intercept = grey(0.45)
)
label_for <- function(model) {
  ifelse(model %in% names(model_labels), model_labels[model], model)
}

theme_two_stage <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", hjust = 0),
        legend.position = "bottom")


# overdispersion ---------------------------------------------------------------

# per insecticide type, from the hierarchical MCMC fit to replicate bioassays
# (R/estimate_bioassay_rho.R); there is no class fallback on purpose, so a
# type missing from the table stops the script rather than being scored at a
# different rho from the rest
rho_table <- read.csv("outputs/bioassay_rho_hierarchical.csv")
if (any(rho_table$worst_rhat > 1.05)) {
  warning("the hierarchical rho fit has rhat above 1.05")
}
rho_for_type <- function(insecticide_type) {
  rho <- rho_table$rho[match(insecticide_type, rho_table$insecticide_type)]
  if (anyNA(rho)) {
    stop("no per-type rho for: ",
         paste(unique(insecticide_type[is.na(rho)]), collapse = ", "))
  }
  rho
}


# find the folds ---------------------------------------------------------------

# every model file for one fold, from both draws directories. The model name is
# taken from the file name, which for the two-stage files carries the variant
# and the leave-out flag (two_stage_<variant>[_loo])
find_model_files <- function(experiment, fold) {
  suffix <- sprintf("__%s__%s.rds", experiment, fold)
  base <- file.path(draws_dir, paste0(c("dynamical", null_models), suffix))
  two_stage <- if (dir.exists(two_stage_dir)) {
    list.files(two_stage_dir, full.names = TRUE,
               pattern = paste0("^two_stage_.*", gsub("\\.", "\\\\.", suffix),
                                "$"))
  } else {
    character(0)
  }
  files <- c(base, two_stage)
  files <- files[file.exists(files)]
  tibble(file = files,
         model = sub("__.*$", "", basename(files)))
}

inventory <- fold_specs %>%
  rowwise() %>%
  mutate(files = list(find_model_files(experiment, fold))) %>%
  ungroup() %>%
  unnest(files)

has_dynamical <- inventory %>%
  group_by(experiment, fold) %>%
  summarise(dynamical = any(model == "dynamical"), .groups = "drop")

skipped <- has_dynamical %>% filter(!dynamical)
if (nrow(skipped) > 0) {
  report("no dynamical draws yet for %s%s",
         paste(skipped$experiment, skipped$fold, sep = "/", collapse = ", "),
         if (allow_no_dynamical) "; scoring the nulls anyway" else
           "; skipping those folds")
}
if (!allow_no_dynamical) {
  inventory <- inventory %>%
    semi_join(has_dynamical %>% filter(dynamical), by = c("experiment", "fold"))
}
if (nrow(inventory) == 0) stop("nothing to score")

report("scoring %i model folds:", nrow(inventory))
print(as.data.frame(inventory %>%
  group_by(label, fold) %>%
  summarise(models = paste(model, collapse = ", "), .groups = "drop")),
  row.names = FALSE)


# memory -----------------------------------------------------------------------

# RAM is shared with long MCMC runs, and a dynamical fold takes 4-8 GB once
# loaded, so wait until there is room before reading a large file
available_gb <- function() {
  meminfo <- readLines("/proc/meminfo")
  line <- grep("^MemAvailable:", meminfo, value = TRUE)
  as.numeric(gsub("[^0-9]", "", line)) / 1024 ^ 2
}

wait_for_memory <- function(file, minimum_gb = 12, large_bytes = 5e8,
                            poll_seconds = 60, give_up_hours = 3) {
  if (file.size(file) < large_bytes) return(invisible())
  waited <- 0
  while (available_gb() < minimum_gb) {
    if (waited == 0) {
      report("waiting for %.0f GB free before reading %s (%.1f GB now)",
             minimum_gb, basename(file), available_gb())
    }
    if (waited > give_up_hours * 3600) {
      stop("gave up waiting for memory to read ", file)
    }
    Sys.sleep(poll_seconds)
    waited <- waited + poll_seconds
  }
  invisible()
}


# per record -------------------------------------------------------------------

# Read one fold and reduce it to what scoring needs, freeing the rest at once:
# the dynamical folds carry the greta draws and arrays, most of their size
extract_fold <- function(file) {
  wait_for_memory(file)
  fold <- readRDS(file)
  p_draws <- thin_draws(fold$p_draws)
  stopifnot(ncol(p_draws) == nrow(fold$test_df))

  # the model's own fitted overdispersion per record, a diagnostic only; which
  # form it takes depends on the model (see validation_metrics.R)
  rho_fitted <- if (!is.null(fold$rho_class_draws)) {
    colMeans(thin_draws(fold$rho_class_draws))[fold$class_id]
  } else if (!is.null(fold$rho_draws)) {
    colMeans(thin_draws(fold$rho_draws))
  } else if (!is.null(fold$rho_implied)) {
    rep(fold$rho_implied, nrow(fold$test_df))
  } else {
    rep(NA_real_, nrow(fold$test_df))
  }
  if (length(rho_fitted) == 1) rho_fitted <- rep(rho_fitted, ncol(p_draws))

  out <- list(test_df = as.data.frame(fold$test_df),
              p_draws = p_draws,
              rho_fitted = rho_fitted)
  rm(fold)
  gc(verbose = FALSE)
  out
}

score_fold <- function(file, model, experiment, fold, label) {

  x <- extract_fold(file)
  test <- x$test_df

  # a fixed seed per fold, so a cached fold and a fresh one give the same PIT
  # randomisation and simulated replicates
  set.seed(2026 + sum(utf8ToInt(paste(model, experiment, fold))))

  rho <- rho_for_type(test$insecticide_type)
  # byrow: ppd_summary would recycle a bare vector down the draws, not across
  # the records
  rho_draws <- matrix(rho, nrow = nrow(x$p_draws), ncol = nrow(test),
                      byrow = TRUE)

  summary <- ppd_summary(test$died, test$mosquito_number, x$p_draws, rho_draws)
  pit <- ppd_pit(summary, n_rep = n_pit_reps)
  sims <- ppd_simulate(test$mosquito_number, x$p_draws, rho_draws)
  crps <- ppd_crps(test$died, test$mosquito_number, sims)

  # years ahead of the last training year, for the forecasting folds: the
  # training data are year_start < cut, so the cut year itself is one year ahead
  horizon <- if (experiment == "temporal_forecasting") {
    test$year_start - as.integer(fold) + 1
  } else {
    NA_real_
  }

  scores <- summary %>%
    mutate(
      model = model,
      experiment = label,
      fold = fold,
      row = seq_len(nrow(test)),
      cell = test$cell,
      region = test$region,
      country_name = test$country_name,
      year_start = test$year_start,
      horizon = horizon,
      insecticide_type = test$insecticide_type,
      insecticide_class = test$insecticide_class,
      # one randomisation replicate, not the mean: the mean converges to the
      # mid-P value, which is not uniform for discrete data (#12 review)
      pit = pit[, 1],
      crps = crps,
      predicted_sd = apply(x$p_draws, 2, sd),
      rho = rho,
      rho_fitted = x$rho_fitted,
      .before = everything()
    )

  list(scores = scores, pit = pit)
}

# cached per fold, keyed on the source file's size and time and on the scoring
# version, so a regenerated draws file is always rescored
score_or_load <- function(file, model, experiment, fold, label) {
  info <- file.info(file)
  key <- list(version = scoring_version, size = info$size,
              mtime = as.numeric(info$mtime))
  cache_file <- file.path(cache_dir, sprintf("%s__%s__%s.rds", model,
                                             experiment, fold))
  if (file.exists(cache_file)) {
    cached <- readRDS(cache_file)
    if (identical(cached$key, key)) {
      report("cached   %s / %s / %s", model, experiment, fold)
      return(cached$value)
    }
  }
  report("scoring  %s / %s / %s (%.2f GB)", model, experiment, fold,
         info$size / 1e9)
  value <- score_fold(file, model, experiment, fold, label)
  saveRDS(list(key = key, value = value), cache_file)
  value
}

# smallest files first, so everything that fits is scored while waiting for
# room to read a dynamical fold
scored <- vector("list", nrow(inventory))
for (i in order(file.size(inventory$file))) {
  entry <- inventory[i, ]
  scored[[i]] <- score_or_load(entry$file, entry$model, entry$experiment,
                               entry$fold, entry$label)
}
names(scored) <- with(inventory, paste(model, label, fold, sep = "|"))


# align records across models ----------------------------------------------------

# Every comparison is paired, so each model's records must be the dynamical
# model's records in the same order. The #12 nulls were built from the same
# test_df; the two-stage files should be too, but that is checked rather than
# assumed, and a file whose records cannot be matched is dropped with a message.
record_key <- function(scores) {
  key <- paste(scores$cell, scores$year_start, scores$insecticide_type,
               scores$died, scores$mosquito_number)
  # an occurrence count separates genuine duplicate assays
  paste(key, ave(seq_along(key), key, FUN = seq_along))
}

aligned <- list()
for (fold_label in unique(paste(inventory$label, inventory$fold))) {
  in_fold <- which(paste(inventory$label, inventory$fold) == fold_label)
  models_here <- inventory$model[in_fold]
  reference_index <- in_fold[match(if ("dynamical" %in% models_here)
    "dynamical" else models_here[1], models_here)]
  reference_key <- record_key(scored[[reference_index]]$scores)

  for (i in in_fold) {
    entry <- scored[[i]]
    key <- record_key(entry$scores)
    if (identical(key, reference_key)) {
      aligned[[length(aligned) + 1]] <- entry
      next
    }
    order <- match(reference_key, key)
    if (anyNA(order) || length(key) != length(reference_key)) {
      message(sprintf(
        "dropping %s on %s: %i of %i held-out records do not match the dynamical model's",
        inventory$model[i], fold_label, sum(is.na(order)), length(reference_key)))
      next
    }
    report("reordered %s on %s to match the dynamical model",
           inventory$model[i], fold_label)
    entry$scores <- entry$scores[order, ] %>% mutate(row = seq_len(n()))
    entry$pit <- entry$pit[order, , drop = FALSE]
    aligned[[length(aligned) + 1]] <- entry
  }
}

all_scores <- bind_rows(lapply(aligned, `[[`, "scores"))
pit_matrices <- lapply(aligned, `[[`, "pit")
names(pit_matrices) <- vapply(aligned, function(entry) {
  with(entry$scores[1, ], paste(model, experiment, fold, sep = "|"))
}, character(1))


# the noise floor ----------------------------------------------------------------

# Regularised as in variance_explained.R: the per-record sampling variance of
# the observed proportion is p(1-p) k, with k = (1 + (m - 1) rho) / m, and p(1-p)
# is taken as its posterior mean under a Beta prior fitted per insecticide type
# to the held-out records, with each assay carrying 1/k effective draws. The raw
# plug-in estimator of #12 is exactly zero at 0% and 100% mortality, which is a
# large share of the records for the saturated insecticides. The floor is a
# property of the record, so it is computed once per record and any subset's
# floor is the mean over its records.
design_effect <- function(m, rho) (1 + (m - 1) * rho) / m

beta_prior_nll <- function(par, observed, m_eff) {
  a <- exp(par[1])
  b <- exp(par[2])
  successes <- observed * m_eff
  -sum(lbeta(a + successes, b + m_eff - successes) - lbeta(a, b))
}

# each held-out record once, whichever model it came from
records <- all_scores %>%
  distinct(experiment, fold, row, .keep_all = TRUE) %>%
  select(experiment, fold, row, insecticide_type, observed, mosquito_number,
         rho)

beta_priors <- bind_rows(lapply(split(records, records$insecticide_type),
                                function(d) {
  m_eff <- 1 / design_effect(d$mosquito_number, d$rho)
  fit <- optim(c(0, 0), beta_prior_nll, observed = d$observed, m_eff = m_eff,
               method = "Nelder-Mead", control = list(reltol = 1e-10))
  data.frame(insecticide_type = d$insecticide_type[1],
             a = exp(fit$par[1]), b = exp(fit$par[2]),
             converged = fit$convergence == 0)
}))
stopifnot(all(beta_priors$converged))

floor_contribution <- function(observed, m, rho, a, b) {
  k <- design_effect(m, rho)
  # a single mosquito carries no information about p(1-p): dropped, as in #12
  usable <- m > 1 & k > 0 & k < 1
  m_eff <- 1 / k
  shape_a <- a + observed * m_eff
  shape_b <- b + (1 - observed) * m_eff
  ifelse(usable,
         shape_a * shape_b / ((shape_a + shape_b) * (shape_a + shape_b + 1)) * k,
         NA_real_)
}

all_scores <- all_scores %>%
  left_join(beta_priors %>% select(insecticide_type, a, b),
            by = "insecticide_type") %>%
  mutate(floor = floor_contribution(observed, mosquito_number, rho, a, b),
         squared_error = (observed - predicted) ^ 2) %>%
  select(-a, -b)

write.csv(all_scores %>% mutate(model_label = label_for(model), .after = model),
          file.path(output_dir, "cv_scores_two_stage.csv"), row.names = FALSE)
report("wrote cv_scores_two_stage.csv (%i rows)", nrow(all_scores))


# summaries with a paired pixel bootstrap ------------------------------------------

# Bioassays cluster hard by pixel, so pixels are the resampling unit, and all
# models' records for a pixel move together: the models are scored on the same
# records, so the shared difficulty cancels in the differences, which are much
# better determined than either model's score alone.
#
# Everything reported is a ratio of sums over records, so the bootstrap works on
# per-pixel sums and a multinomial weight matrix rather than on resampled data
# frames: one matrix product per group instead of 2000 re-summaries.
#
# Per group and model:
#   elpd        mean beta-binomial log predictive density (higher is better)
#   crps        mean CRPS on the mortality scale (lower is better)
#   mse         mean squared error of the predictive mean
#   explained   100 (1 - mse / Var(y)), the out-of-sample variance explained;
#               no noise floor enters it (variance_explained.R)
#   ceiling     100 (1 - floor / Var(y)), the most any model could explain
#   skill       (mse_intercept - mse) / (mse_intercept - floor): 0 is the
#               insecticide mean, 1 is the noise floor (#12)
# and the difference of each from the dynamical model's, with 95% intervals and
# the bootstrap probability that the model beats the dynamical model.
summarise_group <- function(data, pit_list) {

  models <- sort(unique(data$model))
  reference_rows <- data %>% filter(model == models[1]) %>% arrange(row_id)
  cells <- sort(unique(reference_rows$cell))
  cell_index <- match(reference_rows$cell, cells)

  per_cell <- function(values) {
    as.vector(rowsum(values, cell_index, reorder = TRUE))
  }

  # record-level sums shared by every model
  common <- cbind(
    n = per_cell(rep(1, nrow(reference_rows))),
    y = per_cell(reference_rows$observed),
    y2 = per_cell(reference_rows$observed ^ 2),
    floor = per_cell(ifelse(is.na(reference_rows$floor), 0,
                            reference_rows$floor)),
    floor_n = per_cell(as.numeric(!is.na(reference_rows$floor)))
  )

  per_model <- lapply(models, function(m) {
    d <- data %>% filter(model == m) %>% arrange(row_id)
    stopifnot(identical(d$row_id, reference_rows$row_id))
    cbind(log = per_cell(d$log_score),
          crps = per_cell(d$crps),
          se = per_cell(d$squared_error))
  })
  names(per_model) <- models

  metrics_from <- function(weights) {
    # weights: n_cells x n_replicates; returns a list of n_replicates x models
    totals <- crossprod(weights, common)
    n <- totals[, "n"]
    variance <- totals[, "y2"] / n - (totals[, "y"] / n) ^ 2
    floor <- totals[, "floor"] / totals[, "floor_n"]
    out <- lapply(c("log", "crps", "se"), function(stat) {
      sapply(models, function(m) {
        crossprod(weights, per_model[[m]][, stat])[, 1] / n
      }, simplify = "matrix")
    })
    names(out) <- c("elpd", "crps", "mse")
    for (name in names(out)) {
      out[[name]] <- matrix(out[[name]], ncol = length(models),
                            dimnames = list(NULL, models))
    }
    out$explained <- 100 * (1 - out$mse / variance)
    out$skill <- if ("intercept" %in% models) {
      (out$mse[, "intercept"] - out$mse) / (out$mse[, "intercept"] - floor)
    } else {
      out$mse * NA
    }
    out$variance <- variance
    out$ceiling <- 100 * (1 - floor / variance)
    out
  }

  point <- metrics_from(matrix(1, nrow = length(cells), ncol = 1))
  weights <- rmultinom(n_bootstrap, length(cells), rep(1, length(cells)))
  replicates <- metrics_from(weights)

  interval <- function(x, probability) {
    apply(x, 2, quantile, probability, na.rm = TRUE)
  }

  result <- tibble(
    model = models,
    n = nrow(reference_rows),
    n_pixels = length(cells),
    variance = point$variance[1],
    ceiling = point$ceiling[1]
  )

  for (metric in c("elpd", "crps", "mse", "explained", "skill")) {
    result[[metric]] <- point[[metric]][1, ]
    result[[paste0(metric, "_lower")]] <- interval(replicates[[metric]], 0.025)
    result[[paste0(metric, "_upper")]] <- interval(replicates[[metric]], 0.975)
    if ("dynamical" %in% models) {
      difference <- replicates[[metric]] - replicates[[metric]][, "dynamical"]
      # higher is better for all but crps and mse
      better <- if (metric %in% c("crps", "mse")) difference < 0 else
        difference > 0
      result[[paste0("diff_", metric)]] <-
        point[[metric]][1, ] - point[[metric]][1, "dynamical"]
      result[[paste0("diff_", metric, "_lower")]] <- interval(difference, 0.025)
      result[[paste0("diff_", metric, "_upper")]] <- interval(difference, 0.975)
      result[[paste0("prob_better_", metric)]] <- colMeans(better)
    }
  }

  # calibration from the full PIT randomisation matrices
  calibration <- bind_rows(lapply(models, function(m) {
    d <- data %>% filter(model == m)
    pit <- do.call(rbind, lapply(split(d, d$source), function(part) {
      pit_list[[part$source[1]]][part$row, , drop = FALSE]
    }))
    coverage <- coverage_curve(pit, levels = c(0.5, 0.95))
    tibble(model = m,
           mean_pit = mean(pit),
           coverage_50 = coverage$empirical[1],
           coverage_95 = coverage$empirical[2],
           cvm = pit_statistic(pit, cvm_stat),
           bias = mean(d$predicted - d$observed))
  }))

  result %>% left_join(calibration, by = "model")
}

all_scores <- all_scores %>%
  mutate(source = paste(model, experiment, fold, sep = "|"),
         row_id = paste(fold, row))

# a pooled experiment includes only the models present on every one of its
# folds, so that each model is summarised over the same records
models_complete <- function(data) {
  folds_needed <- n_distinct(data$fold)
  data %>%
    group_by(model) %>%
    filter(n_distinct(fold) == folds_needed) %>%
    ungroup()
}

summarise_levels <- function(data, stratum_column = NULL) {
  groups <- bind_rows(
    data %>% mutate(fold_level = "pooled"),
    # separate folds only where an experiment has more than one
    data %>%
      group_by(experiment) %>%
      filter(n_distinct(fold) > 1) %>%
      ungroup() %>%
      mutate(fold_level = fold)
  )
  if (is.null(stratum_column)) {
    groups$stratum <- "all"
  } else {
    groups$stratum <- groups[[stratum_column]]
  }
  groups %>%
    group_by(experiment, fold_level, stratum) %>%
    group_modify(function(d, key) {
      if (key$fold_level == "pooled") d <- models_complete(d)
      reference <- d %>% filter(model == first(model))
      # nothing to bootstrap between, or no variance to explain
      if (n_distinct(reference$cell) < 2 || var(reference$observed) == 0) {
        return(tibble())
      }
      summarise_group(d, pit_matrices)
    }) %>%
    ungroup()
}

set.seed(2026 - 9 - 25)
report("bootstrapping summaries")
summary_table <- bind_rows(
  summarise_levels(all_scores),
  summarise_levels(all_scores, "insecticide_type")
) %>%
  mutate(model_label = label_for(model), .after = model) %>%
  rename(fold = fold_level)

write.csv(summary_table, file.path(output_dir, "cv_summary_two_stage.csv"),
          row.names = FALSE)
report("wrote cv_summary_two_stage.csv (%i rows)", nrow(summary_table))

headline <- summary_table %>%
  filter(fold == "pooled", stratum == "all")

cat("\nheadline, pooled over folds (per-type rho):\n")
print(as.data.frame(headline %>%
  transmute(experiment, model, n, n_pixels,
            elpd = round(elpd, 3),
            crps = round(crps, 4),
            explained = sprintf("%.1f [%.1f, %.1f]", explained,
                                explained_lower, explained_upper),
            ceiling = round(ceiling, 1),
            skill = round(skill, 2),
            coverage_95 = round(coverage_95, 3),
            d_elpd = sprintf("%+.3f [%+.3f, %+.3f]", diff_elpd,
                             diff_elpd_lower, diff_elpd_upper))),
  row.names = FALSE)


# coverage curves --------------------------------------------------------------------

coverage_table <- bind_rows(lapply(
  split(all_scores %>% distinct(model, experiment, fold, source),
        ~ model + experiment, drop = TRUE),
  function(entry) {
    pit <- do.call(rbind, pit_matrices[entry$source])
    coverage_curve(pit, levels = coverage_levels) %>%
      mutate(model = entry$model[1], experiment = entry$experiment[1],
             folds = paste(sort(entry$fold), collapse = "+"),
             .before = everything())
  })) %>%
  # the same completeness rule as the pooled summaries
  group_by(experiment) %>%
  filter(folds == folds[model == "dynamical"][1] | !any(model == "dynamical")) %>%
  ungroup() %>%
  mutate(model_label = label_for(model), .after = model)

write.csv(coverage_table, file.path(output_dir, "cv_coverage_two_stage.csv"),
          row.names = FALSE)


# skill by forecast horizon ------------------------------------------------------------

horizon_table <- all_scores %>%
  filter(!is.na(horizon))
if (nrow(horizon_table) > 0) {
  report("bootstrapping by forecast horizon")
  horizon_table <- horizon_table %>%
    mutate(stratum = paste("horizon", horizon)) %>%
    group_by(experiment, horizon) %>%
    group_modify(function(d, key) {
      reference <- d %>% filter(model == first(model))
      if (n_distinct(reference$cell) < 2) return(tibble())
      summarise_group(d, pit_matrices)
    }) %>%
    ungroup() %>%
    mutate(model_label = label_for(model), .after = model)
  write.csv(horizon_table, file.path(output_dir, "cv_horizon_two_stage.csv"),
            row.names = FALSE)
}


# figures ------------------------------------------------------------------------------

present_models <- names(model_labels)[names(model_labels) %in%
                                        unique(all_scores$model)]
present_models <- c(present_models,
                    setdiff(unique(all_scores$model), present_models))
colour_scale <- function(aesthetic = "colour") {
  values <- model_colours[present_models]
  values[is.na(values)] <- "black"
  names(values) <- label_for(present_models)
  scale_discrete_manual(aesthetic, values = values,
                        breaks = label_for(present_models), name = NULL)
}
as_model_factor <- function(model) {
  factor(label_for(model), levels = label_for(present_models))
}
experiment_labels <- c(
  spatial_interpolation = "spatial interpolation",
  spatial_blocks = "spatial blocks",
  temporal_forecasting_2014 = "forecast from 2014",
  temporal_forecasting_2018 = "forecast from 2018"
)
as_experiment_factor <- function(experiment) {
  factor(experiment_labels[experiment],
         levels = experiment_labels[experiment_labels %in%
                                      experiment_labels[unique(experiment)]])
}

save_figure <- function(plot, name, width, height) {
  ggsave(file.path(figure_dir, name), plot, width = width, height = height,
         dpi = 200, bg = "white")
  report("wrote %s", file.path(figure_dir, name))
}

# PIT histograms, pooling all randomisation replicates for a smoother histogram;
# a calibrated model is flat at 1. A U shape is too narrow, a hump too wide, and
# a slope a bias
pit_long <- bind_rows(lapply(
  split(all_scores %>% distinct(model, experiment, fold, source),
        ~ model + experiment, drop = TRUE),
  function(entry) {
    pit <- as.vector(do.call(rbind, pit_matrices[entry$source]))
    tibble(model = entry$model[1], experiment = entry$experiment[1],
           bin = cut(pit, seq(0, 1, by = 0.05), include.lowest = TRUE,
                     labels = FALSE)) %>%
      count(model, experiment, bin) %>%
      mutate(density = n / sum(n) / 0.05, pit = (bin - 0.5) * 0.05)
  }))

pit_plot <- ggplot(pit_long %>%
                     mutate(model = as_model_factor(model),
                            experiment = as_experiment_factor(experiment)),
                   aes(pit, density, fill = model)) +
  geom_col(width = 0.05, colour = "white", linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = "22", colour = grey(0.3)) +
  facet_grid(model ~ experiment, labeller = label_wrap_gen(18)) +
  colour_scale("fill") +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "randomised PIT", y = "density",
       title = "Probability integral transform of held-out bioassays",
       subtitle = "Flat at 1 when calibrated; all models scored at the per-type overdispersion") +
  theme_two_stage +
  theme(legend.position = "none", strip.text.y = element_text(angle = 0))
save_figure(pit_plot, "pit_histograms.png",
            width = 2 + 2 * n_distinct(pit_long$experiment),
            height = 1.5 + 1.2 * n_distinct(pit_long$model))

# coverage, as the gap from nominal so that small miscalibration is visible
coverage_plot <- ggplot(coverage_table %>%
                          mutate(model = as_model_factor(model),
                                 experiment = as_experiment_factor(experiment)),
                        aes(nominal, empirical - nominal, colour = model)) +
  geom_hline(yintercept = 0, colour = grey(0.4)) +
  geom_line(linewidth = 0.7) +
  geom_point(data = function(d) d %>% filter(abs(nominal - 0.95) < 1e-9),
             size = 2) +
  facet_wrap(~ experiment, nrow = 1) +
  colour_scale() +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
  labs(x = "nominal coverage of central predictive interval",
       y = "empirical minus nominal",
       title = "Coverage of held-out bioassays",
       subtitle = "Above zero: intervals too wide; below: too narrow. Point marks the 95% interval") +
  guides(colour = guide_legend(nrow = 2)) +
  theme_two_stage
save_figure(coverage_plot, "coverage_curves.png",
            width = 2 + 2.6 * n_distinct(coverage_table$experiment), height = 4.2)

# variance explained and differences from the dynamical model, pooled per
# experiment. Intervals are the paired pixel bootstrap
headline_plot_data <- headline %>%
  mutate(model = as_model_factor(model),
         experiment = as_experiment_factor(experiment))

explained_plot <- ggplot(headline_plot_data,
                         aes(explained, model, colour = model)) +
  geom_vline(xintercept = 0, colour = grey(0.4)) +
  geom_vline(aes(xintercept = ceiling), linetype = "22", colour = grey(0.4)) +
  geom_errorbar(aes(xmin = explained_lower, xmax = explained_upper),
                 width = 0, orientation = "y", linewidth = 0.7) +
  geom_point(size = 2.2) +
  geom_text(aes(label = sprintf("%.0f", explained)), vjust = -0.9,
            size = 2.8, colour = grey(0.2), show.legend = FALSE) +
  facet_wrap(~ experiment, nrow = 1, scales = "free_y") +
  colour_scale() +
  scale_y_discrete(limits = rev) +
  labs(x = "held-out variance explained (%)", y = NULL,
       title = "Out-of-sample variance explained",
       subtitle = "95% pixel-bootstrap intervals; dashed line is the ceiling set by bioassay noise") +
  theme_two_stage +
  theme(legend.position = "none")

difference_data <- headline %>%
  filter(model != "dynamical") %>%
  select(experiment, model, starts_with("diff_")) %>%
  pivot_longer(starts_with("diff_"), names_to = "name") %>%
  mutate(metric = sub("^diff_([a-z]+).*$", "\\1", name),
         part = case_when(grepl("_lower$", name) ~ "lower",
                          grepl("_upper$", name) ~ "upper",
                          TRUE ~ "estimate")) %>%
  select(-name) %>%
  pivot_wider(names_from = part, values_from = value) %>%
  filter(metric %in% c("elpd", "crps", "explained")) %>%
  # sign so that right is always better than the dynamical model
  mutate(across(c(estimate, lower, upper),
                ~ ifelse(metric == "crps", -.x, .x)),
         flipped_lower = pmin(lower, upper), upper = pmax(lower, upper),
         lower = flipped_lower,
         metric = factor(metric, c("elpd", "crps", "explained"),
                         c("log score", "CRPS (sign flipped)",
                           "variance explained (points)")),
         model = as_model_factor(model),
         experiment = as_experiment_factor(experiment))

difference_plot <- ggplot(difference_data, aes(estimate, model, colour = model)) +
  geom_vline(xintercept = 0, colour = grey(0.4)) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, orientation = "y", linewidth = 0.7) +
  geom_point(size = 2.2) +
  facet_grid(experiment ~ metric, scales = "free", space = "free_y",
             labeller = label_wrap_gen(20)) +
  colour_scale() +
  scale_y_discrete(limits = rev) +
  labs(x = "difference from the dynamical model (right is better)", y = NULL,
       title = "Each model against the dynamical model, on the same held-out bioassays",
       subtitle = "Mean per-assay difference, 95% paired pixel-bootstrap interval") +
  theme_two_stage +
  theme(legend.position = "none", strip.text.y = element_text(angle = 0))

save_figure(explained_plot, "variance_explained.png",
            width = 2.5 + 2.6 * n_distinct(headline$experiment), height = 3.4)
save_figure(difference_plot, "difference_vs_dynamical.png",
            width = 10,
            height = 1.5 + 1.1 * n_distinct(headline$experiment) *
              max(1, n_distinct(difference_data$model)) / 2)

if (is.data.frame(horizon_table) && nrow(horizon_table) > 0) {
  horizon_plot <- ggplot(horizon_table %>%
                           mutate(model = as_model_factor(model),
                                  experiment = as_experiment_factor(experiment)),
                         aes(horizon, explained, colour = model)) +
    geom_hline(yintercept = 0, colour = grey(0.4)) +
    geom_ribbon(aes(ymin = explained_lower, ymax = explained_upper,
                    fill = model), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2) +
    facet_wrap(~ experiment, nrow = 1) +
    colour_scale() +
    colour_scale("fill") +
    labs(x = "years ahead of the last training year",
         y = "held-out variance explained (%)",
         title = "Forecast skill by horizon",
         subtitle = "95% pixel-bootstrap intervals within each year") +
    guides(colour = guide_legend(nrow = 2), fill = "none") +
    theme_two_stage
  save_figure(horizon_plot, "skill_by_horizon.png",
              width = 2 + 3.2 * n_distinct(horizon_table$experiment),
              height = 4.2)

  if ("dynamical" %in% horizon_table$model &&
      any(grepl("^two_stage", horizon_table$model))) {
    gain_plot <- ggplot(horizon_table %>%
                          filter(grepl("^two_stage", model)) %>%
                          mutate(model = as_model_factor(model),
                                 experiment = as_experiment_factor(experiment)),
                        aes(horizon, diff_explained, colour = model)) +
      geom_hline(yintercept = 0, colour = grey(0.4)) +
      geom_pointrange(aes(ymin = diff_explained_lower,
                          ymax = diff_explained_upper),
                      position = position_dodge(width = 0.3)) +
      facet_wrap(~ experiment, nrow = 1) +
      colour_scale() +
      labs(x = "years ahead of the last training year",
           y = "gain in variance explained over\nthe dynamical model (points)",
           title = "What the correction adds, by forecast horizon") +
      theme_two_stage
    save_figure(gain_plot, "skill_gain_by_horizon.png",
                width = 2 + 3.2 * n_distinct(horizon_table$experiment),
                height = 4)
  }
}


# diagnostics: training residuals at the extremes ----------------------------------------

# Stage A treats the empirical logit z as Gaussian around m_ref + latent with
# variance v. That approximation is worst where mortality is 0% or 100%: the
# +0.5 continuity correction pins z there, and v is at its most wrong. If the
# standardised residuals (z - m_ref - latent) / sqrt(v) at the extremes are
# centred and of unit scale, like the rest, stage A is adequate; a systematic
# offset or a scale far from 1 there, especially one that varies with m or by
# region, is the case for stage B (PQL on the beta-binomial likelihood).
#
# The residual files are written by the stage-A fits; the columns are found by
# name, with a few alternatives, so a change of naming there does not break
# this. Absent files are skipped.

region_lookup <- NULL
get_region_lookup <- function() {
  # cell to region, from the full bioassay data, for residual files without a
  # region column; built lazily as it needs the raster mask
  if (is.null(region_lookup)) {
    mask <- terra::rast("data/clean/raster_mask.tif")
    data <- readRDS("data/clean/all_gambiae_complex_data.RDS")
    region_lookup <<- tibble(
      cell = terra::cellFromXY(mask, as.matrix(data[, c("longitude",
                                                         "latitude")])),
      region = data$region
    ) %>%
      filter(!is.na(cell)) %>%
      count(cell, region) %>%
      group_by(cell) %>%
      slice_max(n, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(cell, region)
  }
  region_lookup
}

pick_column <- function(data, candidates) {
  found <- intersect(candidates, names(data))
  if (length(found) == 0) NA_character_ else found[1]
}

read_residuals <- function(file) {
  x <- readRDS(file)
  if (!is.data.frame(x)) {
    frames <- Filter(is.data.frame, x)
    if (length(frames) == 0) {
      message("no data frame in ", basename(file), "; skipped")
      return(NULL)
    }
    x <- frames[[1]]
  }
  x <- as.data.frame(x)
  columns <- list(
    z = pick_column(x, c("z", "empirical_logit")),
    v = pick_column(x, c("v", "variance")),
    m_ref = pick_column(x, c("m_ref", "m")),
    latent = pick_column(x, c("latent", "latent_mean", "fitted_latent",
                              "correction", "fitted", "latent_fitted")),
    died = pick_column(x, c("died")),
    n = pick_column(x, c("mosquito_number", "n")),
    type = pick_column(x, c("insecticide_type", "type"))
  )
  missing <- names(columns)[is.na(unlist(columns))]
  if (any(missing %in% c("z", "v", "m_ref", "latent"))) {
    message(basename(file), " lacks ", paste(missing, collapse = ", "),
            "; skipped")
    return(NULL)
  }

  # mortality class from the counts if present, otherwise from z against the
  # continuity-corrected extremes, which needs the assay size
  if (!is.na(columns$died) && !is.na(columns$n)) {
    died <- x[[columns$died]]
    n <- x[[columns$n]]
  } else {
    message(basename(file), " has no died / mosquito_number; skipped")
    return(NULL)
  }

  parts <- strsplit(sub("\\.rds$", "", basename(file)), "__")[[1]]
  out <- tibble(
    model = parts[2],
    experiment = parts[3],
    fold = parts[4],
    insecticide_type = if (is.na(columns$type)) NA_character_ else
      x[[columns$type]],
    cell = if ("cell" %in% names(x)) x$cell else NA_real_,
    region = if ("region" %in% names(x)) x$region else NA_character_,
    m_ref = x[[columns$m_ref]],
    residual = (x[[columns$z]] - x[[columns$m_ref]] - x[[columns$latent]]) /
      sqrt(x[[columns$v]]),
    mortality = case_when(died == 0 ~ "0%",
                          died == n ~ "100%",
                          TRUE ~ "interior")
  )
  if (all(is.na(out$region)) && !all(is.na(out$cell))) {
    out <- out %>% select(-region) %>%
      left_join(get_region_lookup(), by = "cell")
  }
  out$source <- "observed"

  # The reference to read the observed residuals against. Conditioning on the
  # observed outcome is itself a selection on the residual: an assay only reads
  # 0% if its noise pulled it down, so residuals at 0% are negative, and at 100%
  # positive, even when stage A is exactly right. What would indicate a problem
  # is a departure from what the fitted model implies, so the same summaries are
  # computed on assays simulated from it: beta-binomial at the fitted mean
  # plogis(m_ref + latent) and the per-type rho, pushed through the same
  # empirical logit, and classed by their own simulated mortality. The latent
  # field is treated as known, which it is not, so this is a reference for the
  # selection effect rather than a calibrated test.
  if (!is.na(columns$type)) {
    size <- rep(n, n_residual_sims)
    rho <- rep(rho_for_type(x[[columns$type]]), n_residual_sims)
    fitted <- rep(x[[columns$m_ref]] + x[[columns$latent]], n_residual_sims)
    simulated_died <- rbetabinom(length(size), size, plogis(fitted), rho)
    simulated <- empirical_logit_z(simulated_died, size, rho)
    out <- bind_rows(
      out,
      out[rep(seq_len(nrow(out)), n_residual_sims), ] %>%
        mutate(source = "simulated from the stage-A fit",
               residual = (simulated$z - fitted) / sqrt(simulated$v),
               mortality = case_when(simulated_died == 0 ~ "0%",
                                     simulated_died == size ~ "100%",
                                     TRUE ~ "interior"))
    )
  }
  out
}

# the stage-A response, as in empirical_logit() in R/two_stage_correction.R,
# which is not sourced here since it loads TMB and the mesh code
empirical_logit_z <- function(died, n, rho) {
  list(z = log((died + 0.5) / (n - died + 0.5)),
       v = (1 / (died + 0.5) + 1 / (n - died + 0.5)) * (1 + (n - 1) * rho))
}
n_residual_sims <- 20

summarise_residuals <- function(data) {
  data %>%
    summarise(n = n(),
              mean = mean(residual),
              se_mean = sd(residual) / sqrt(n()),
              sd = sd(residual),
              rms = sqrt(mean(residual ^ 2)),
              beyond_2 = mean(abs(residual) > 2),
              .groups = "drop")
}

residual_diagnostics <- function(residual_files) {
  residuals <- bind_rows(lapply(residual_files, read_residuals))
  if (nrow(residuals) == 0) return(NULL)

  residuals <- residuals %>%
    mutate(
      region = ifelse(is.na(region), "unknown", region),
      # bins of the dynamical model's predicted mortality at the assay
      m_bin = cut(plogis(m_ref), c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1),
                  include.lowest = TRUE),
      mortality = factor(mortality, c("0%", "interior", "100%"))
    )

  table <- bind_rows(
    residuals %>%
      group_by(model, experiment, fold, source, mortality) %>%
      summarise_residuals() %>%
      mutate(grouping = "overall", group = "all"),
    residuals %>%
      group_by(model, experiment, fold, source, mortality, group = region) %>%
      summarise_residuals() %>%
      mutate(grouping = "region"),
    residuals %>%
      group_by(model, experiment, fold, source, mortality,
               group = as.character(m_bin)) %>%
      summarise_residuals() %>%
      mutate(grouping = "m_bin")
  ) %>%
    # assays per replicate, so the observed and simulated counts in each
    # class compare directly: more observed 0% or 100% readings than the fit
    # implies is itself a sign the Gaussian response misses the extremes
    mutate(assays = ifelse(source == "observed", n, n / n_residual_sims),
           .after = n) %>%
    relocate(grouping, group, .after = fold) %>%
    arrange(model, experiment, fold, grouping, group, mortality, source)

  list(residuals = residuals, table = table)
}

plot_residual_diagnostics <- function(table) {
  mortality_colours <- c("0%" = "#B2182B", "interior" = grey(0.45),
                         "100%" = "#2166AC")
  # bins in order of mortality, then regions alphabetically
  bin_levels <- levels(cut(0.5, c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1),
                           include.lowest = TRUE))
  group_levels <- c(bin_levels,
                    sort(unique(table$group[table$grouping == "region"])))
  # the mean (target 0) and the SD (target 1) side by side: an offset at the
  # extremes is a bias stage A cannot see, and an SD far from 1 is v being wrong
  plot_data <- table %>%
    filter(grouping %in% c("m_bin", "region"), n >= 5) %>%
    mutate(fold_label = paste0(label_for(model), ", ", experiment, " / ",
                               fold),
           group = factor(group, rev(group_levels)),
           grouping = factor(grouping, c("m_bin", "region"),
                             c("by dynamical-model\nmortality", "by region")))
  plot_data <- bind_rows(
    plot_data %>% mutate(metric = "mean (target 0)", value = mean,
                         lower = mean - 2 * se_mean,
                         upper = mean + 2 * se_mean, target = 0),
    plot_data %>% mutate(metric = "SD (target 1)", value = sd,
                         lower = NA, upper = NA, target = 1)
  )

  # observed filled, the simulated reference hollow at the same position: the
  # gap between the two is what matters, not the distance from the target line
  ggplot(plot_data, aes(value, group, colour = mortality, group = mortality)) +
    geom_vline(aes(xintercept = target), colour = grey(0.4)) +
    geom_errorbar(data = function(d) d %>% filter(source == "observed"),
                  aes(xmin = lower, xmax = upper), width = 0,
                  linewidth = 0.6, orientation = "y", na.rm = TRUE,
                  position = position_dodge(width = 0.6)) +
    geom_point(aes(shape = source), size = 2.2, stroke = 0.8,
               position = position_dodge(width = 0.6)) +
    scale_shape_manual(values = c("observed" = 16,
                                  "simulated from the stage-A fit" = 1),
                       name = NULL) +
    facet_grid(grouping ~ fold_label + metric, scales = "free",
               space = "free_y", labeller = label_wrap_gen(28)) +
    scale_colour_manual(values = mortality_colours, name = "observed mortality",
                        breaks = names(mortality_colours)) +
    labs(x = "standardised training residual, (z - m_ref - latent) / sqrt(v)",
         y = NULL,
         title = "Stage-A training residuals at 0% and 100% mortality",
         subtitle = "Mean \u00b1 2 SE, and SD. Conditioning on 0% or 100% biases residuals even under a correct model;\nthe case for stage B is a gap between observed (filled) and simulated from the fit (hollow)") +
    theme_two_stage +
    theme(strip.text.y = element_text(angle = 0))
}

residual_files <- list.files(output_dir, pattern = "^train_residuals__.*\\.rds$",
                             full.names = TRUE)
residual_files <- residual_files[grepl(
  paste0("__(", paste(paste(fold_specs$experiment, fold_specs$fold, sep = "__"),
                      collapse = "|"), ")\\.rds$"), residual_files)]

set.seed(2026 - 9 - 26)
if (length(residual_files) == 0) {
  report("no training residual files yet; skipping the residual diagnostics")
} else {
  diagnostics <- residual_diagnostics(residual_files)
  if (!is.null(diagnostics)) {
    write.csv(diagnostics$table,
              file.path(output_dir, "train_residual_extremes.csv"),
              row.names = FALSE)
    report("wrote train_residual_extremes.csv")

    cat("\ntraining residuals, standardised, by mortality class:\n")
    print(as.data.frame(diagnostics$table %>%
      filter(grouping == "overall") %>%
      select(-grouping, -group) %>%
      mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)

    residual_plot <- plot_residual_diagnostics(diagnostics$table)
    save_figure(residual_plot, "train_residual_extremes.png",
                width = 2.5 + 5 * n_distinct(paste(diagnostics$table$model,
                                                   diagnostics$table$experiment,
                                                   diagnostics$table$fold)),
                height = 7)
  }
}

report("done in %.1f minutes",
       as.numeric(difftime(Sys.time(), start_time, units = "mins")))
