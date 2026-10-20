# Check of the final model's CCD integration over the hyperparameters (#34)
# against importance sampling (R/two_stage_is.R), for one insecticide type's
# full-data fit saved by R/two_stage_maps.R (fit step).
#
#   Rscript R/check_two_stage_is.R <type>
#
# Run with single-threaded BLAS and one forked process per core, e.g.
#   OPENBLAS_NUM_THREADS=1 IR_CUBE_IS_CORES=8 Rscript R/check_two_stage_is.R Deltamethrin
#
# Writes to outputs/two_stage/maps/<type>/:
#   is_check_hyper.csv   per hyperparameter, the CCD and importance sampling
#                        posterior median and 95% interval;
#   is_check_fields.csv  for omega and for xi at the last data year T, the
#                        difference between the two posteriors' means at the
#                        mesh nodes, in units of the CCD posterior SD there
#                        (from 500 draws): its mean and maximum absolute value;
#   is_check_diagnostics.csv  Pareto k (reliable below 0.7), the effective
#                        sample size, and the effective sample size after each
#                        adaptation round.

source("R/two_stage_helpers.R")
source("R/two_stage_correction.R")
source("R/two_stage_is.R")

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1)
type <- arguments[1]
type_dir <- file.path("outputs/two_stage/maps", type)
fit <- readRDS(file.path(type_dir, "fit.rds"))
stopifnot(!is.null(fit$ccd))

set.seed(string_seed(paste("is check", type)))
is <- theta_is(fit, n = 200)
fit_is <- fit
fit_is$ccd <- is

# the hyperparameters' posterior median and 95% interval under each
weighted_quantile <- function(x, w, q) {
  o <- order(x)
  x[o][findInterval(q, cumsum(w[o]) / sum(w)) + 1]
}
names_hyper <- c("sigma_omega", "range_omega", "sigma_eta", "range_eta",
                 "persistence", vapply(correction_iid_terms, `[[`, "", "sd"))
hyper <- bind_rows(lapply(names_hyper, function(name) {
  bind_rows(lapply(list(ccd = fit$ccd, is = is), function(p) {
    tibble(q2.5 = weighted_quantile(p$hyper[[name]], p$weight, 0.025),
           q50 = weighted_quantile(p$hyper[[name]], p$weight, 0.5),
           q97.5 = weighted_quantile(p$hyper[[name]], p$weight, 0.975))
  }), .id = "method") %>%
    mutate(parameter = name, .before = 1)
}))

# the posterior means of omega and xi(T) at the nodes under each, against the
# CCD posterior SD there
mean_ccd <- correction_node_draws(fit, fit$T, mean = TRUE)
mean_is <- correction_node_draws(fit_is, fit$T, mean = TRUE)
draws <- correction_node_draws(fit, fit$T, 500)
field_check <- function(field, a, b, d) {
  z <- (a - b) / apply(d, 1, sd)
  tibble(field = field, mean_abs_z = mean(abs(z)), max_abs_z = max(abs(z)))
}
key <- as.character(fit$T)
fields <- bind_rows(
  field_check("omega", mean_ccd$omega[, 1], mean_is$omega[, 1], draws$omega),
  field_check(paste("xi", fit$T), mean_ccd$xi[[key]][, 1],
              mean_is$xi[[key]][, 1], draws$xi[[key]]))

diagnostics <- tibble(type = type, pareto_k = is$diagnostics$pareto_k,
                      ess = is$diagnostics$ess,
                      round_ess = paste(round(is$diagnostics$round_ess),
                                        collapse = " "),
                      n_failed = is$diagnostics$n_failed)
write.csv(hyper, file.path(type_dir, "is_check_hyper.csv"), row.names = FALSE)
write.csv(fields, file.path(type_dir, "is_check_fields.csv"),
          row.names = FALSE)
write.csv(diagnostics, file.path(type_dir, "is_check_diagnostics.csv"),
          row.names = FALSE)
print(as.data.frame(diagnostics))
print(as.data.frame(hyper), digits = 3)
print(as.data.frame(fields), digits = 3)
