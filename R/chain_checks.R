# Whether the chains of a fit agree (#39). Every chain starts from the same
# initial values (dynamical_inits_file), so a chain that reaches another mode,
# or has not converged, shows only as disagreement between the chains.
# check_chains() summarises that for the fit's log, and the result is stored
# with the fit, for R/fit_checks.R, which tabulates every fit.

# The chains of `draws` that stopped moving: those whose share of distinct
# draws is below `min_distinct`. With one step size for all chains
# (windowed_hmc()), a chain can sit where every trajectory at that step size
# fails, and then repeats one draw: in the October 2026 run, one chain of
# spatial blocks fold 1 kept a single draw for all 3,000 samples.
stuck_chains <- function(draws, min_distinct = 0.5) {
  which(vapply(draws, function(chain) {
    chain <- as.matrix(chain)
    nrow(unique(chain)) / nrow(chain) < min_distinct
  }, logical(1)))
}

# The variables whose posteriors R/fit_checks.R compares between fits: the
# floor's class means and type deviations, the latent smooths' sds, the types'
# mean initial state and overdispersion, reversion per class, and the sds of
# the selection effects' class level.
compared_variables <- c("floor_flat", "floor_sd_type", "floor_type_raw",
                        "smooth_sd_selection", "smooth_sd_floor",
                        "logit_init_mean", "logit_rho_type", "reversion_rate",
                        "sigma_class")

# The model's log density at `n` evenly spaced draws of each chain of `draws`,
# as an n x chains matrix: greta's, on the free state, with the log Jacobian
# of the transforms, which is the density the sampler targets (as Stan's
# lp__). By the compiled function the sampler used, on all the chains' states
# at once, as it was called while sampling.
chain_log_density <- function(model, draws, n = 200) {
  free <- lapply(attr(draws, "model_info")$raw_draws, as.matrix)
  rows <- round(seq(1, nrow(free[[1]]), length.out = n))
  t(vapply(rows, function(i) {
    state <- do.call(rbind, lapply(free, function(chain) chain[i, ]))
    log_density <- model$dag$tf_log_prob_function_adjusted(
      tensorflow::tf$constant(state, dtype = model$dag$tf_float))
    as.numeric(log_density)
  }, numeric(length(free))))
}

# The chains of `draws` (greta draws), with the log density `lp`
# (chain_log_density()) and the labels of the variables' dimensions `labels`
# (built$variable_labels):
#   chains     per chain: whether it is stuck (stuck_chains()), and the mean
#              and sd of its log density
#   lp_rhat    over the kept chains (those not stuck, or all if every one
#              is), the rank-normalised Rhat of the log density
#              (posterior::rhat())
#   flagged    lp_rhat above 1.01 (Vehtari et al. 2021)
#   differing  the `top` elements whose chain means differ most, over the
#              kept chains: (max - min of the chain means) / the root mean
#              within-chain variance
#   variables  the 5%, 50% and 95% quantiles of each element of
#              compared_variables, per chain and over the kept chains ("kept")
check_chains <- function(draws, lp, labels, top = 5) {
  stuck <- seq_along(draws) %in% stuck_chains(draws)
  # if every chain is stuck, all are kept, as in fit_fold()
  kept <- if (all(stuck)) seq_along(draws) else which(!stuck)
  chains <- data.frame(chain = seq_along(draws), stuck = stuck,
                       lp_mean = colMeans(lp), lp_sd = apply(lp, 2, sd))
  lp_rhat <- if (length(kept) > 1) {
    posterior::rhat(lp[, kept, drop = FALSE])
  } else {
    NA_real_
  }

  draws <- lapply(draws, as.matrix)
  means <- vapply(draws[kept], colMeans, numeric(ncol(draws[[1]])))
  variances <- vapply(draws[kept], function(x) apply(x, 2, var),
                      numeric(ncol(draws[[1]])))
  spread <- (apply(as.matrix(means), 1, max) -
               apply(as.matrix(means), 1, min)) /
    sqrt(rowMeans(as.matrix(variances)))
  spread <- sort(spread[is.finite(spread)], decreasing = TRUE)
  spread <- spread[seq_len(min(top, length(spread)))]

  element <- colnames(draws[[1]])
  compared <- sub("\\[.*$", "", element) %in% compared_variables
  quantiles <- function(x, chain) {
    q <- apply(x[, compared, drop = FALSE], 2, quantile,
               c(0.05, 0.5, 0.95))
    data.frame(variable = sub("\\[.*$", "", element[compared]),
               element = label_elements(element[compared], labels),
               chain = chain, q05 = q[1, ], q50 = q[2, ], q95 = q[3, ],
               row.names = NULL)
  }
  variables <- rbind(
    do.call(rbind, Map(quantiles, draws, as.character(seq_along(draws)))),
    quantiles(do.call(rbind, draws[kept]), "kept"))

  list(chains = chains,
       lp_rhat = lp_rhat,
       flagged = isTRUE(lp_rhat > 1.01),
       differing = data.frame(element = label_elements(names(spread), labels),
                              spread = unname(spread)),
       variables = variables)
}

# Greta's names of variable elements, e.g. "init_country_level[12,7]", with
# their indices replaced by the labels of the variable's dimensions
# (`labels`), leaving out dimensions of length one.
label_elements <- function(elements, labels) {
  vapply(elements, function(element) {
    name <- sub("\\[.*$", "", element)
    if (name == element || is.null(labels[[name]])) return(element)
    index <- as.integer(strsplit(sub("^.*\\[(.*)\\]$", "\\1", element),
                                 ",")[[1]])
    dimensions <- labels[[name]]
    long <- lengths(dimensions) > 1
    parts <- mapply(function(d, i) d[i], dimensions[long], index[long])
    sprintf("%s[%s]", name, paste(parts, collapse = ", "))
  }, character(1), USE.NAMES = FALSE)
}

# Lines for the fit's log, from check_chains(): the log density's Rhat and
# each chain's mean (sd), and when flagged a WARNING naming the elements whose
# chain means differ most.
chain_check_lines <- function(check) {
  chains <- check$chains
  lines <- sprintf("chains: log density Rhat %.3f | chain means %s",
                   check$lp_rhat,
                   paste(sprintf("%.0f (sd %.0f)%s", chains$lp_mean,
                                 chains$lp_sd,
                                 ifelse(chains$stuck, " stuck", "")),
                         collapse = ", "))
  if (check$flagged) {
    lines <- c(lines, sprintf(
      "WARNING: chains disagree (log density Rhat %.3f); differing most: %s",
      check$lp_rhat,
      paste(sprintf("%s %.1f", check$differing$element,
                    check$differing$spread), collapse = ", ")))
  }
  lines
}
