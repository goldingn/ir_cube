# Summarise the fits of the #37 grid (R/net_grid.R), fetched with
# `irpod fetch <name>` to outputs/pod_jobs/<name>/.
#
#   Rscript R/net_grid_report.R
#
# For every fit found: each chain's mortality-floor mode and log posterior
# (chain_modes, R/chain_floor_modes.R), and for the folds the held-out scores
# with all chains and with the chains of each mode alone. Scores are at the
# external replicate-based rho, as R/validation_metrics.R scores them: the
# summed log predictive density (ELPD), the mean of the predictive mean's
# squared error, and the mean PIT (mid-P). Writes
#   outputs/net_grid/chain_modes.csv   one row per chain
#   outputs/net_grid/mode_summary.csv  per fit: chains in each mode, mean
#                                      log posterior of each mode and the
#                                      difference (low - high)
#   outputs/net_grid/cv_scores.csv     per fold and chain set
#   outputs/net_grid/cv_summary.csv    per variant and chain set, summed
#                                      over the folds of each experiment
# Loading a fold takes 4-8 GB; they are read one at a time.

suppressMessages({
  library(dplyr)
  library(tidyr)
})
source("R/net_grid.R")
source("R/validation_functions.R")
source("R/validation_scoring.R")
jobs_dir <- "outputs/pod_jobs"
output_dir <- "outputs/net_grid"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
rho_source <- rho_lookup()

jobs <- net_grid_jobs(code_ref = NA_character_)
jobs$dir <- file.path(jobs_dir, jobs$name)
jobs$file <- ifelse(
  jobs$fit == "full",
  file.path(jobs$dir, "temporary/fitted_model.RData"),
  file.path(jobs$dir, "outputs/cv_draws",
            sprintf("dynamical__%s__%s.rds", jobs$experiment, jobs$fold)))
found <- file.exists(jobs$file)
cat(sprintf("%d of %d grid fits found\n", sum(found), nrow(jobs)))
if (!any(found)) {
  quit(save = "no")
}

# The chain of each stored prediction draw: fit_fold() kept `stored` draws
# evenly spaced over the chains stacked in order
stored_draw_chain <- function(draws, n_stored) {
  chain <- draw_chain(draws)
  chain[round(seq(1, length(chain), length.out = n_stored))]
}

# the scores of a fold's test set from prediction draws p_draws
score_set <- function(test, p_draws) {
  rho <- rho_for_record(test, rho_source)
  summary <- ppd_summary(test$died, test$mosquito_number, p_draws, rho)
  data.frame(n_draws = nrow(p_draws), n_assays = nrow(test),
             elpd = sum(summary$log_score),
             mse = mean((summary$observed - summary$predicted) ^ 2),
             mean_pit = mean(summary$cdf_below + 0.5 * summary$pmf_at))
}

mode_rows <- list()
score_rows <- list()
for (i in which(found)) {
  job <- jobs[i, ]
  cat(format(Sys.time(), "%H:%M:%S"), job$name, "\n")
  if (job$fit == "full") {
    e <- new.env()
    load(job$file, envir = e)
    modes <- e$chain_modes
    rm(e)
  } else {
    fold <- readRDS(job$file)
    modes <- fold$chain_modes
    chain <- stored_draw_chain(fold$draws, nrow(fold$p_draws))
    sets <- list(all = unique(chain))
    if (!is.null(modes)) {
      for (m in unique(modes$mode)) {
        sets[[m]] <- modes$chain[modes$mode == m]
      }
    }
    for (set in names(sets)) {
      rows <- chain %in% sets[[set]]
      score_rows[[length(score_rows) + 1]] <- data.frame(
        variant = job$variant, fit = job$fit, experiment = fold$experiment,
        chains = set, chain_ids = paste(sets[[set]], collapse = ","),
        score_set(fold$test_df, fold$p_draws[rows, , drop = FALSE]))
    }
    rm(fold)
  }
  invisible(gc())
  if (!is.null(modes)) {
    mode_rows[[length(mode_rows) + 1]] <- data.frame(
      variant = job$variant, fit = job$fit, modes)
  }
}

chain_modes <- bind_rows(mode_rows)
write.csv(chain_modes, file.path(output_dir, "chain_modes.csv"),
          row.names = FALSE)
if (nrow(chain_modes) > 0) {
  mode_summary <- chain_modes |>
    group_by(variant, fit, mode) |>
    summarise(chains = paste(chain, collapse = ","), floor = mean(floor),
              log_posterior = mean(log_posterior), .groups = "drop") |>
    pivot_wider(names_from = mode,
                values_from = c(chains, floor, log_posterior)) |>
    mutate(across(any_of(c("chains_low", "chains_high")),
                  ~ replace_na(.x, "")))
  if (all(c("log_posterior_low", "log_posterior_high") %in%
          names(mode_summary))) {
    mode_summary$log_posterior_low_minus_high <-
      mode_summary$log_posterior_low - mode_summary$log_posterior_high
  }
  write.csv(mode_summary, file.path(output_dir, "mode_summary.csv"),
            row.names = FALSE)
  print(as.data.frame(mode_summary), digits = 6)
}

cv_scores <- bind_rows(score_rows)
write.csv(cv_scores, file.path(output_dir, "cv_scores.csv"), row.names = FALSE)
if (nrow(cv_scores) > 0) {
  cv_summary <- cv_scores |>
    group_by(variant, experiment, chains) |>
    summarise(folds = n(), n_assays = sum(n_assays), elpd = sum(elpd),
              mse = weighted.mean(mse, n_assays),
              mean_pit = weighted.mean(mean_pit, n_assays),
              .groups = "drop") |>
    arrange(experiment, chains, -elpd)
  write.csv(cv_summary, file.path(output_dir, "cv_summary.csv"),
            row.names = FALSE)
  print(as.data.frame(cv_summary), digits = 5)
}
