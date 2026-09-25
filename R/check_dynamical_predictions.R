# Check R/dynamical_predictions.R reproduces the saved held-out predictions of a
# fold, draw for draw, and that the recursion is exactly additive on the logit
# scale. Run from the repo root:
#   Rscript R/check_dynamical_predictions.R outputs/cv_draws/dynamical__spatial_interpolation__all.rds interp
#   Rscript R/check_dynamical_predictions.R outputs/cv_draws/dynamical__temporal_forecasting__2020.rds forecast

fold_file <- commandArgs(trailingOnly = TRUE)[1]
experiment <- commandArgs(trailingOnly = TRUE)[2]
suppressMessages({
  sink("/dev/null"); source("R/validation_folds.R"); source("R/validation_covariates.R"); sink()
})
source("R/dynamical_predictions.R")

# layout checks the greta reshape relies on
stopifnot(n_times == 30,
          identical(as.integer(cell_years_index$cell_id), rep(seq_len(n_unique_cells), each = n_times)),
          identical(as.integer(cell_years_index$year_id), rep(seq_len(n_times), n_unique_cells)))

fold <- readRDS(fold_file)
cat("p_draws dim", dim(fold$p_draws), "\n")
draw_index <- paired_draw_index(fold)
p_ref <- fold$p_draws
if (nrow(p_ref) > 2000) p_ref <- p_ref[round(seq(1, nrow(p_ref), length.out = 2000)), ]
fold$p_draws <- NULL; fold$rho_draws <- NULL; fold$p_draws_before <- NULL; fold$rho_class_draws <- NULL
invisible(gc())

if (experiment == "interp") {
  training <- spatial_interpolation$training; test <- spatial_interpolation$test
} else {
  f <- temporal_forecasting_folds[["2020"]]; training <- f$training; test <- f$test
}
stopifnot(identical(test$cell_id, fold$test_df$cell_id), identical(test$year_id, fold$test_df$year_id),
          identical(test$type_id, fold$test_df$type_id))

t_test <- system.time(p_test <- dynamical_predictions(fold, fold$test_df, df, x_cell_years, cell_years_index, classes_index, types))
diff <- abs(p_test - p_ref)
cat(sprintf("TEST: n=%d, max abs diff %.3g, max rel logit diff %.3g, time %.1fs\n", ncol(p_test), max(diff),
            max(abs(qlogis(p_test) - qlogis(p_ref))), t_test[["elapsed"]]))
rm(p_test, p_ref, diff); invisible(gc(reset = TRUE))

cat(sprintf("TRAIN: %d rows, %d unique cell-type-years\n", nrow(training),
            nrow(distinct(training, cell_id, type_id, year_id))))
t_train <- system.time(p_train <- dynamical_predictions(fold, training, df, x_cell_years, cell_years_index, classes_index, types))
g <- gc()
cat(sprintf("TRAIN: time %.1fs, result %s, R max used since reset %.0f MB\n", t_train[["elapsed"]],
            format(object.size(p_train), units = "MB"), sum(g[, ncol(g)])))
cat("range", range(p_train), "NA", anyNA(p_train), "\n")
rm(p_train); invisible(gc())

# logit additivity against the direct (q-scale) recursion, as haploid_next does it
pars <- dynamical_parameter_draws(fold, df, classes_index, types, draw_index)
set.seed(1)
cases <- training %>% distinct(cell_id, type_id, country_id) %>% slice_sample(n = 6)
delta <- 0.5
worst <- 0; worst_vs_fn <- 0
for (i in seq_len(nrow(cases))) {
  c_id <- cases$cell_id[i]; k <- cases$type_id[i]; cc <- cases$country_id[i]
  x <- x_cell_years[cell_years_index$cell_id == c_id, , drop = FALSE]
  w <- 1 + x %*% t(pars$effect_type[, , k])      # n_times x draws
  l0 <- pars$logit_init[, cc, k]
  run <- function(q) { out <- matrix(NA, n_times, length(q)); for (t in seq_len(n_times)) { q <- q / (q + (1 - q) * w[t, ]); out[t, ] <- q }; out }
  q_base <- run(plogis(l0)); q_pert <- run(plogis(l0 + delta))
  d <- qlogis(q_pert) - qlogis(q_base)
  worst <- max(worst, abs(d - delta))
  rows <- tibble(cell_id = c_id, type_id = k, year_id = seq_len(n_times), country_id = cc)
  p_fn <- dynamical_predictions(fold, rows, df, x_cell_years, cell_years_index, classes_index, types, draw_index)
  worst_vs_fn <- max(worst_vs_fn, abs(p_fn - t(q_base)))
  cat(sprintf("case %d cell %d type %s: min q %.3g, max|dlogit - delta| %.3g\n", i, c_id, types[k], min(q_base), max(abs(d - delta))))
}
cat(sprintf("ADDITIVITY: worst |dlogit - %.1f| over all years/draws = %.3g; direct recursion vs function max abs diff %.3g\n", delta, worst, worst_vs_fn))
