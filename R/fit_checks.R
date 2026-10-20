# The chain checks of every fit in one table and one figure (check_chains(),
# R/chain_checks.R, stored with each fit), so that chains that disagree, and
# fits whose posteriors differ, are read in one place, not in each fit's log.
# Run after the full fit and the folds, before scoring.
#
#   Rscript R/fit_checks.R [fit file ...]
#
# Default: the full fit (temporary/fitted_model.RData) and every dynamical fold
# in outputs/cv_draws/. Writes outputs/fit_checks.csv, one row per fit:
#   fit             "full", or the fold's experiment and fold
#   chains          chains run
#   stuck           the stuck chains, left out of the rest
#   lp_rhat         rank-normalised Rhat of the log density
#   rhat_worst      worst rank-normalised Rhat over every element
#   above_1.01      number of elements with Rhat above 1.01, of `elements`
#   elements        number of elements
#   flagged         the chains disagree (check_chains())
#   chain_lp_means  each chain's mean log density (nats), in chain order
#   worst_elements  the 5 elements with the highest Rhat, with it
# and figures/fit_checks.png: compared_variables by fit, the median and 90%
# interval over the kept chains, and each chain's median. Loading a fold takes
# 4 to 8 GB. Run with the greta library (R/packages.R), which the saved
# objects need.

suppressMessages({
  library(greta)
  library(dplyr)
  library(ggplot2)
})
source("R/chain_checks.R")

arguments <- commandArgs(trailingOnly = TRUE)
files <- if (length(arguments) > 0) arguments else
  c("temporary/fitted_model.RData",
    Sys.glob("outputs/cv_draws/dynamical__*.rds"))
files <- files[file.exists(files)]
stopifnot(length(files) > 0)

checks <- lapply(files, function(file) {
  if (grepl("\\.RData$", file)) {
    fit <- new.env()
    load(file, envir = fit)
    check <- fit$chain_check
    name <- "full"
  } else {
    fold <- readRDS(file)
    check <- fold$chain_check
    # the forecasting folds' experiment labels end with the fold
    name <- paste(sub(paste0("_", fold$fold, "$"), "", fold$experiment),
                  fold$fold)
  }
  if (is.null(check)) stop(file, " has no chain check (fitted before #39)")
  c(check, fit = name)
})

table <- bind_rows(lapply(checks, function(check) tibble(
  fit = check$fit,
  chains = nrow(check$chains),
  stuck = paste(which(check$chains$stuck), collapse = ", "),
  lp_rhat = check$lp_rhat,
  rhat_worst = max(check$rhat, na.rm = TRUE),
  above_1.01 = sum(check$rhat > 1.01, na.rm = TRUE),
  elements = length(check$rhat),
  flagged = check$flagged,
  chain_lp_means = paste(sprintf("%.0f", check$chains$lp_mean),
                         collapse = " / "),
  worst_elements = worst_elements(check))))
write.csv(table, "outputs/fit_checks.csv", row.names = FALSE)
print(as.data.frame(table), digits = 4)

variables <- bind_rows(lapply(checks, function(check) {
  mutate(check$variables, fit = check$fit)
}))
panels <- setNames(sprintf("%s (%s)", names(compared_variables),
                           compared_variables), names(compared_variables))
dodge <- position_dodge(width = 0.7)
figure <- ggplot(filter(variables, chain == "kept"),
                 aes(y = label, colour = fit, group = fit)) +
  geom_linerange(aes(xmin = q05, xmax = q95), position = dodge) +
  geom_point(aes(x = q50), position = dodge, size = 1.2) +
  geom_point(aes(x = q50), data = filter(variables, chain != "kept"),
             position = dodge, shape = 1, size = 0.8) +
  facet_wrap(~ variable, scales = "free", labeller = as_labeller(panels)) +
  labs(x = paste("posterior median and 90% interval over the kept chains",
                 "(filled), and each chain's median (open), in the panel's",
                 "units"),
       y = NULL) +
  theme_minimal(base_size = 8)
ggsave("figures/fit_checks.png", figure, width = 11, height = 9, dpi = 200,
       bg = "white")
