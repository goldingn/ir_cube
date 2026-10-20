# The chain checks of every fit in one table and one figure (check_chains(),
# R/chain_checks.R, stored with each fit), so that chains that disagree, and
# fits whose posteriors differ, are read in one place, not in each fit's log.
# Run after the full fit and the folds, before scoring.
#
#   Rscript R/fit_checks.R [fit file ...]
#
# Default: the full fit (temporary/fitted_model.RData) and every dynamical fold
# in outputs/cv_draws/. A fit made before #39 has no chain check, and gets a
# row without one. Writes:
#   outputs/fit_checks.csv  one row per fit (columns below)
#   figures/fit_checks.png  compared_variables (R/chain_checks.R) by fit: the
#                           median and 90% interval over the kept chains, and
#                           each chain's median
# Loading a fold takes 4 to 8 GB. Run with the greta library (R/packages.R),
# which the saved objects need.

suppressMessages({
  library(greta)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})
# drop_chains(), and check_chains() (R/chain_checks.R)
source("R/dynamical_predictions.R")

arguments <- commandArgs(trailingOnly = TRUE)
files <- if (length(arguments) > 0) arguments else
  c("temporary/fitted_model.RData",
    Sys.glob("outputs/cv_draws/dynamical__*.rds"))
files <- files[file.exists(files)]
stopifnot(length(files) > 0)

# one fit's row, and its compared variables
check_fit <- function(file) {
  if (grepl("\\.RData$", file)) {
    saved <- new.env()
    load(file, envir = saved)
    saved <- list(name = "full", draws = saved$draws,
                  check = saved$chain_check, stuck = saved$stuck_chains)
  } else {
    fold <- readRDS(file)
    # the forecasting folds' experiment labels end with the fold
    experiment <- sub(paste0("_", fold$fold, "$"), "", fold$experiment)
    saved <- list(name = paste(experiment, fold$fold),
                  draws = fold$draws, check = fold$chain_check,
                  stuck = fold$stuck_chains)
    rm(fold)
  }
  # the kept chains: those not stuck, or all if every one is (a full fit
  # processed by R/drop_stuck_chains.R has only its kept chains)
  draws <- saved$draws
  stuck <- stuck_chains(draws)
  if (length(stuck) < length(draws)) {
    draws <- drop_chains(draws, setdiff(seq_along(draws), stuck))
  }
  rhat <- if (length(draws) > 1) {
    coda::gelman.diag(draws, multivariate = FALSE,
                      autoburnin = FALSE)$psrf[, 1]
  } else {
    NA_real_
  }
  check <- saved$check
  row <- tibble(
    fit = saved$name,
    chains = if (is.null(check)) length(saved$draws) else nrow(check$chains),
    stuck = paste(saved$stuck, collapse = ", "),
    rhat_worst = max(rhat, na.rm = TRUE),
    parameters_above_1.01 = sum(rhat > 1.01, na.rm = TRUE),
    parameters = length(rhat),
    lp_rhat = if (is.null(check)) NA_real_ else check$lp_rhat,
    chain_lp_means = if (is.null(check)) NA_character_ else
      paste(sprintf("%.0f", check$chains$lp_mean), collapse = " / "),
    flagged = if (is.null(check)) NA else check$flagged,
    differing_most = if (is.null(check)) NA_character_ else
      paste(sprintf("%s %.1f", check$differing$element,
                    check$differing$spread), collapse = "; "))
  variables <- if (!is.null(check)) mutate(check$variables, fit = saved$name)
  list(row = row, variables = variables)
}

checks <- lapply(files, function(file) {
  cat("checking", file, "\n")
  out <- check_fit(file)
  invisible(gc())
  out
})

# Columns of fit_checks.csv:
#   fit                    "full", or the fold's experiment and fold
#   chains                 chains run
#   stuck                  the stuck chains dropped (stuck_chains())
#   rhat_worst             worst Rhat (coda::gelman.diag()) over the raw
#                          parameters, kept chains
#   parameters_above_1.01  number of raw parameters with Rhat above 1.01
#   parameters             number of raw parameters
#   lp_rhat                rank-normalised Rhat of the log density, kept chains
#   chain_lp_means         each chain's mean log density (nats), in chain order
#   flagged                lp_rhat above 1.01: the chains disagree
#   differing_most         the elements whose chain means differ most, each
#                          with (max - min of the chain means) / within-chain
#                          sd (sds)
table <- bind_rows(lapply(checks, `[[`, "row"))
write.csv(table, "outputs/fit_checks.csv", row.names = FALSE)
options(width = 200)
print(as.data.frame(select(table, -differing_most)), digits = 4)
flagged <- filter(table, flagged %in% TRUE)
cat(sprintf("\n%d of %d fits flagged (chains disagree)%s\n", nrow(flagged),
            nrow(table), if (nrow(flagged) > 0)
              paste0(": ", paste(flagged$fit, collapse = ", ")) else ""))
for (i in seq_len(nrow(flagged))) {
  cat(sprintf("  %s: chain means %s; differing most: %s\n", flagged$fit[i],
              flagged$chain_lp_means[i], flagged$differing_most[i]))
}

# the figure, with each panel's units
variables <- bind_rows(lapply(checks, `[[`, "variables"))
if (nrow(variables) > 0) {
  titles <- c(
    floor_flat = "floor_flat: mean floor hazard per class\n(floor = 1 - exp(-hazard))",
    floor_sd_type = "floor_sd_type: sd of the types'\nlog floor hazard within class",
    floor_type_raw = "floor_type_raw: type's deviation of\nlog floor hazard (in floor_sd_type)",
    smooth_sd_selection = "smooth_sd_selection: sd of the\nsmooth of log selection strength",
    smooth_sd_floor = "smooth_sd_floor: sd of the\nsmooth of log floor hazard",
    logit_init_mean = "logit_init_mean: logit of the type's mean\nrelative initial susceptible fraction",
    logit_rho_type = "logit_rho_type: logit of the\ntype's overdispersion rho",
    reversion_rate = "reversion_rate: per class,\nlogit per year",
    sigma_class = "sigma_class: sd of the types' log\nselection effects within class")
  # each element by its labels alone (the panel names the variable)
  variables <- mutate(variables, label = ifelse(
    grepl("\\[", element), sub("^[^[]*\\[(.*)\\]$", "\\1", element), ""))
  kept <- filter(variables, chain == "kept")
  per_chain <- filter(variables, chain != "kept")
  dodge <- position_dodge(width = 0.7)
  figure <- ggplot(kept, aes(y = label, colour = fit, group = fit)) +
    geom_linerange(aes(xmin = q05, xmax = q95), position = dodge) +
    geom_point(aes(x = q50), position = dodge, size = 1.2) +
    geom_point(aes(x = q50), data = per_chain, position = dodge, shape = 1,
               size = 0.8) +
    facet_wrap(~ variable, scales = "free",
               labeller = as_labeller(titles[unique(variables$variable)])) +
    labs(x = "posterior median and 90% interval over the kept chains (filled),
and each chain's median (open), in the units of the panel's title",
         y = NULL, colour = "fit") +
    theme_minimal(base_size = 8)
  ggsave("figures/fit_checks.png", figure, width = 11, height = 9, dpi = 200,
         bg = "white")
}
