# Score the #34 experiments (R/two_stage_hyper_folds.R) on their held-out
# assays, as R/validation_metrics.R scores a fold: every assay at the
# external, replicate-based overdispersion of its type (rho_lookup()).
#
#   Rscript R/two_stage_hyper_scores.R
#
# Reads outputs/two_stage/hyper/parts/*.rds and writes, to
# outputs/two_stage/hyper/:
#   scores.csv   one row per method, fold and type:
#     elpd        mean log predictive probability of the observed count (nats
#                 per assay)
#     crps        mean CRPS of the observed mortality fraction (fraction units)
#     cov50,
#     cov95       expected coverage of the central 50% and 95% predictive
#                 intervals of the count (share of assays)
#     mean_pit    mean randomised PIT (0.5 if unbiased; above 0.5, the
#                 predictions are too low)
#     width95     mean width of the central 95% interval of the predicted
#                 mortality fraction, before assay noise (fraction units)
#     the hyperparameter summaries and timings of each part
#   pooled.csv   the same means pooled over the types every method in the fold
#                has, one row per method and fold

source("R/validation_scoring.R")
suppressMessages(library(dplyr))

hyper_dir <- "outputs/two_stage/hyper"
files <- list.files(file.path(hyper_dir, "parts"), pattern = "\\.rds$",
                    full.names = TRUE)
stopifnot(length(files) > 0)
rho_source <- rho_lookup()

scores <- bind_rows(lapply(files, function(file) {
  part <- readRDS(file)
  test <- part$test
  p_draws <- plogis(part$logit_draws)
  rho <- rho_for_record(test, rho_source)
  set.seed(string_seed(paste(basename(file), "score")))
  summary <- ppd_summary(test$died, test$mosquito_number, p_draws,
                         matrix(rho, nrow(p_draws), ncol(p_draws),
                                byrow = TRUE))
  pit <- ppd_pit(summary, n_rep = 1)
  sims <- ppd_simulate(test$mosquito_number, p_draws,
                       matrix(rho, nrow(p_draws), ncol(p_draws), byrow = TRUE))
  q <- apply(p_draws, 2, quantile, c(0.025, 0.975))
  bind_cols(
    part$summary[c("method", "experiment", "fold", "insecticide_type",
                   "n_test")],
    tibble(elpd = mean(log(summary$pmf_at)),
           crps = mean(ppd_crps(test$died, test$mosquito_number, sims)),
           cov50 = mean(expected_coverage(summary$cdf_below, summary$pmf_at,
                                          0.5)),
           cov95 = mean(expected_coverage(summary$cdf_below, summary$pmf_at,
                                          0.95)),
           mean_pit = mean(pit[, 1]),
           width95 = mean(q[2, ] - q[1, ])),
    select(part$summary, -c(method, experiment, fold, insecticide_type,
                            n_test))
  )
}))
write.csv(scores, file.path(hyper_dir, "scores.csv"), row.names = FALSE)

# pooled over the types that every method of the fold has, weighted by assays
pooled <- scores %>%
  group_by(experiment, fold) %>%
  filter(insecticide_type %in% Reduce(intersect,
                                      split(insecticide_type, method))) %>%
  group_by(experiment, fold, method) %>%
  summarise(n_types = n(), n_test_total = sum(n_test),
            across(c(elpd, crps, cov50, cov95, mean_pit, width95),
                   ~ weighted.mean(.x, n_test)),
            .groups = "drop")
write.csv(pooled, file.path(hyper_dir, "pooled.csv"), row.names = FALSE)
print(as.data.frame(pooled), digits = 3)
