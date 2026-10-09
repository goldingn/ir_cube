# Convergence and posterior summaries of one saved full fit (#47; the fits of
# R/species_runs.R, or any temporary/fitted_model.RData of R/fit_model.R).
#
#   Rscript R/species_fit_summary.R <fitted_model.RData> <label>
#
# For the key parameters: beta_overall, beta_class (overall plus the class
# deviation, per covariate and class, by the model's own transforms,
# dynamical_terms(), whichever levels are centred, #48), sigma_overall,
# sigma_class, rho per type, the reversion rates, and those of #47 (gamma_*, delta_*, the floors,
# the latent smooths' sd and range),
# the posterior mean, sd and quantiles, R-hat and bulk and tail ESS
# (posterior::summarise_draws()), over the usable chains (all but those
# stuck, stuck_chains()) and, for R-hat, over all chains; gamma_* and delta_*
# also on the multiplier scale, exp(); and those of #47 within each floor
# mode of the usable chains (parameters_by_mode). Per chain: its floor means
# and mode, the share of distinct draws, whether it is stuck, and with
# chain_modes (R/chain_floor_modes.R) its log posterior. Stuck chains are
# reported, and left out of the summaries, never dropped from the fit. Fits
# saved before #47 load with their options completed
# (complete_model_options(), R/species_fit_helpers.R; FLOOR_PRIOR for their
# floor prior).
#
# Writes outputs/species_runs/summary/<label>.rds (parameters,
# parameters_by_mode, chains, fit); R/species_fit_summaries.R combines them.
# Plain R; under 2 GB and 15 s for a full fit.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2)
file <- arguments[1]
label <- arguments[2]

suppressMessages({
  library(greta)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(posterior)
})
source("R/dynamical_predictions.R")
source("R/species_fit_helpers.R")

output_dir <- "outputs/species_runs/summary"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

fit <- load_fit(file)
n_chains <- length(fit$draws)
report("%s: %s, %d chains x %d draws, %d parameters", label, file, n_chains,
       nrow(fit$draws[[1]]), ncol(fit$draws[[1]]))
if (smooth_on(fit$options) && smooth_range_fixed(fit$options$smooth)) {
  report("%s: the smooths' range is fixed at %.0f km", label,
         1000 * fit$options$smooth[["range"]])
}


# the key parameters, per chain ---------------------------------------------------

covariates <- colnames(fit$x_cell_years)
key_parameters <- function(chain) {
  m <- as.matrix(chain)
  p <- function(name) extract_parameter(m, name)
  n <- nrow(m)
  out <- list()
  beta_overall <- matrix(p("beta_overall"), n)
  sigma_overall <- matrix(p("sigma_overall"), n)
  sigma_class <- matrix(p("sigma_class"), n)
  # the class effects and rho per type, by dynamical_terms(): from the
  # standard normal deviations of the non-centred levels and the effects
  # themselves of the centred ones (#48)
  variables <- lapply(setNames(nm = unique(sub("\\[.*$", "", colnames(m)))),
                      extract_parameter, draws_matrix = m)
  terms <- dynamical_terms_draws(variables, fit$classes_index, fit$types,
                                 terms = c("beta_class", "rho_types"),
                                 options = fit$options)
  for (j in seq_along(covariates)) {
    out[[sprintf("beta_overall[%s]", covariates[j])]] <- beta_overall[, j]
    out[[sprintf("sigma_overall[%s]", covariates[j])]] <- sigma_overall[, j]
    out[[sprintf("sigma_class[%s]", covariates[j])]] <- sigma_class[, j]
    for (c in seq_along(fit$classes)) {
      out[[sprintf("beta_class[%s,%s]", covariates[j], fit$classes[c])]] <-
        terms$beta_class[, j, c]
    }
  }
  for (k in seq_along(fit$types)) {
    out[[sprintf("rho[%s]", fit$types[k])]] <- terms$rho_types[, k, 1]
  }
  if (any(grepl("^reversion_rate", colnames(m)))) {
    reversion <- matrix(p("reversion_rate"), n)
    for (c in seq_along(fit$classes)) {
      out[[sprintf("reversion_rate[%s]", fit$classes[c])]] <- reversion[, c]
    }
  }
  new <- grep("^(gamma_|delta_|floor_)|floor$", colnames(m), value = TRUE)
  new <- new[!grepl("^smooth_", new)]
  smooth <- smooth_on(fit$options)
  for (name in new) {
    out[[name]] <- m[, name]
    if (grepl("^(gamma_|delta_)", name)) {
      out[[sprintf("exp(%s)", name)]] <- exp(m[, name])
    }
    # the kdr-dependent floor at the mean kdr, or the floor of the latent
    # smooths where u_f is 0, per class with an intercept per class: from its
    # logit, floor_intercept, or itself, floor_flat (V5f)
    if (grepl("^floor_(intercept|flat)", name)) {
      index <- if (name %in% c("floor_intercept", "floor_flat")) "" else
        sprintf("[%s]", fit$classes[as.integer(sub("^.*\\[(\\d+),.*$", "\\1",
                                                    name))])
      out[[paste0(if (smooth) "floor_at_u0" else "floor_at_k0", index)]] <-
        if (grepl("^floor_flat", name)) m[, name] else plogis(m[, name])
    }
  }
  # the latent smooths' sd and range in km (none when it is fixed, which is
  # no variable)
  for (name in grep("^smooth_sd_", colnames(m), value = TRUE)) {
    out[[name]] <- m[, name]
  }
  for (name in grep("^smooth_inv_range_", colnames(m), value = TRUE)) {
    out[[sub("^smooth_inv_range_", "smooth_range_km_", name)]] <-
      1000 / m[, name]
  }
  do.call(cbind, out)
}
per_chain <- lapply(fit$draws, key_parameters)
# iterations x chains x variables
key_array <- function(chains) {
  a <- simplify2array(per_chain[chains])
  posterior::as_draws_array(aperm(a, c(1, 3, 2)))
}


# chains ----------------------------------------------------------------------------

detected_stuck <- stuck_chains(fit$draws)
usable <- usable_chains(fit$draws, label)
floors <- fit_floor_names(fit$draws)
chains <- tibble(
  label = label,
  chain = seq_len(n_chains),
  mode = unname(chain_floor_mode(fit$draws)),
  distinct_share = vapply(fit$draws, function(chain) {
    chain <- as.matrix(chain)
    nrow(unique(chain)) / nrow(chain)
  }, numeric(1)),
  stuck = chain %in% detected_stuck,
  dropped_by_drop_stuck_chains = chain %in% fit$recorded_stuck)
for (name in floors) {
  chains[[paste0("mean_", name)]] <- vapply(fit$draws, function(chain) {
    mean(as.matrix(chain)[, name])
  }, numeric(1))
}
if (!is.null(fit$chain_modes)) {
  chains <- left_join(chains,
                      select(fit$chain_modes, chain, starts_with("log_posterior")),
                      by = "chain")
}
report("chains: %s", paste(sprintf("%d %s%s", chains$chain, chains$mode,
                                   ifelse(chains$stuck, " (stuck)", "")),
                           collapse = ", "))


# parameters ------------------------------------------------------------------------

quantiles <- function(x) {
  quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), names = FALSE)
}
summary_usable <- posterior::summarise_draws(
  key_array(usable), mean, sd,
  ~ setNames(quantiles(.x), c("q2.5", "q25", "q50", "q75", "q97.5")),
  rhat = posterior::rhat, ess_bulk = posterior::ess_bulk,
  ess_tail = posterior::ess_tail)
rhat_all <- posterior::summarise_draws(key_array(seq_len(n_chains)),
                                       rhat_all_chains = posterior::rhat)
parameters <- summary_usable %>%
  left_join(rhat_all, by = "variable") %>%
  mutate(label = label, .before = 1) %>%
  rename(parameter = variable) %>%
  mutate(group = case_when(
    str_detect(parameter, "^exp\\(") ~ "multiplier",
    str_detect(parameter, "^(gamma_|delta_)") ~ "species and kdr",
    str_detect(parameter, "^smooth_") ~ "smooth",
    str_detect(parameter, "floor$|^floor_") ~ "floor",
    str_detect(parameter, "^rho") ~ "rho",
    str_detect(parameter, "^reversion") ~ "reversion",
    .default = "selection"), .after = parameter)

# the new parameters within each floor mode of the usable chains, where the
# chains are in more than one
modes <- split(usable, chains$mode[usable])
parameters_by_mode <- bind_rows(lapply(names(modes), function(this_mode) {
  posterior::summarise_draws(
    key_array(modes[[this_mode]]), mean, sd,
    ~ setNames(quantiles(.x), c("q2.5", "q25", "q50", "q75", "q97.5")),
    rhat = posterior::rhat) %>%
    rename(parameter = variable) %>%
    filter(str_detect(parameter,
                      "^exp\\(|^(gamma_|delta_|floor_|smooth_)|floor$")) %>%
    mutate(label = label, mode = this_mode,
           chains = toString(modes[[this_mode]]), .before = 1)
}))

chain_table <- chains
fit_row <- tibble(
  label = label,
  file = file,
  n_chains = n_chains,
  draws_per_chain = nrow(fit$draws[[1]]),
  stuck_chains = toString(detected_stuck),
  usable_chains = length(usable),
  modes = paste(chain_table$mode, collapse = ","),
  max_rhat = max(parameters$rhat, na.rm = TRUE),
  max_rhat_all_chains = max(parameters$rhat_all_chains, na.rm = TRUE),
  min_ess_bulk = min(parameters$ess_bulk, na.rm = TRUE),
  min_ess_tail = min(parameters$ess_tail, na.rm = TRUE),
  options = paste(deparse(fit$saved_options[setdiff(
    names(fit$saved_options), c("selection_columns", dynamical_built_options))]),
    collapse = ""))

saveRDS(list(parameters = parameters, parameters_by_mode = parameters_by_mode,
             chains = chains, fit = fit_row),
        file.path(output_dir, sprintf("%s.rds", label)))

options(width = 160)
print(as.data.frame(chains), digits = 4)
print(as.data.frame(fit_row %>% select(-file, -options)), digits = 4)
print(as.data.frame(parameters %>%
                      filter(group %in% c("multiplier", "floor", "species and kdr",
                                          "smooth", "reversion", "rho")) %>%
                      select(parameter, mean, q2.5, q50, q97.5, rhat,
                             rhat_all_chains, ess_bulk, ess_tail)),
      digits = 3)
worst <- parameters %>% arrange(desc(rhat)) %>% head(5)
report("worst R-hat (usable chains): %s",
       paste(sprintf("%s %.3f", worst$parameter, worst$rhat), collapse = ", "))
report("saved %s; peak memory %.1f GB",
       file.path(output_dir, sprintf("%s.rds", label)), peak_memory_gb())
