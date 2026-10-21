# Whether the chains of a fit agree (#39). Every chain starts from the same
# initial values (dynamical_inits_file), so a chain that reaches another mode,
# or has not converged, shows only as disagreement between the chains.
# check_chains() summarises that for the fit's log, and the result is stored
# with the fit, for R/fit_checks.R, which tabulates every fit.

# The chains of `draws` that stopped moving: those whose share of distinct
# draws is below `min_distinct`. With one step size for all chains (as in
# greta's hmc(), and the windowed sampler before #41), a chain can sit where
# every trajectory at that step size fails, and then repeats one draw: in the
# October 2026 run, with the windowed sampler, one chain of spatial blocks fold
# 1 kept a single draw for all 3,000 samples.
stuck_chains <- function(draws, min_distinct = 0.5) {
  which(vapply(draws, function(chain) {
    chain <- as.matrix(chain)
    nrow(unique(chain)) / nrow(chain) < min_distinct
  }, logical(1)))
}

# The variables whose posteriors are flagged and plotted for every fit, with
# their units: the floor's class means and type deviations, the latent
# smooths' sds, the types' mean initial state and overdispersion, reversion
# per class, and the sds of the selection effects' class level.
compared_variables <- c(floor_flat = "hazard",
                        floor_sd_type = "sd of log hazard",
                        floor_type_raw = "sds",
                        smooth_sd_selection = "sd of log selection",
                        smooth_sd_floor = "sd of log hazard",
                        logit_init_mean = "logit",
                        logit_rho_type = "logit",
                        reversion_rate = "logit per year",
                        sigma_class = "sd of log selection effect")

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
# (built$variable_labels). Rhat is rank-normalised (posterior::rhat()), over
# the kept chains: those not stuck, or all if every one is.
#   chains     per chain: whether it is stuck, and the mean and sd of its log
#              density
#   lp_rhat    Rhat of the log density
#   rhat       Rhat of every element, named by its labels
#   flagged    lp_rhat, or the Rhat of an element of compared_variables, above
#              1.01 (Vehtari et al. 2021)
#   variables  the 5%, 50% and 95% quantiles of each element of
#              compared_variables (by its label), per chain and over the kept
#              chains ("kept")
check_chains <- function(draws, lp, labels) {
  stuck <- seq_along(draws) %in% stuck_chains(draws)
  kept <- if (all(stuck)) seq_along(draws) else which(!stuck)
  draws <- lapply(draws, as.matrix)
  variable <- sub("\\[.*$", "", colnames(draws[[1]]))
  label <- label_elements(colnames(draws[[1]]), labels)
  element <- ifelse(label == "", variable, sprintf("%s[%s]", variable, label))
  rhat <- setNames(apply(simplify2array(draws[kept]), 2, posterior::rhat),
                   element)
  lp_rhat <- posterior::rhat(lp[, kept, drop = FALSE])
  compared <- variable %in% names(compared_variables)
  quantiles <- function(x, chain) {
    q <- apply(x[, compared, drop = FALSE], 2, quantile, c(0.05, 0.5, 0.95))
    data.frame(variable = variable[compared], label = label[compared],
               chain = as.character(chain), q05 = q[1, ], q50 = q[2, ],
               q95 = q[3, ], row.names = NULL)
  }
  list(chains = data.frame(chain = seq_along(draws), stuck = stuck,
                           lp_mean = colMeans(lp), lp_sd = apply(lp, 2, sd)),
       lp_rhat = lp_rhat,
       rhat = rhat,
       flagged = isTRUE(max(lp_rhat, rhat[compared], na.rm = TRUE) > 1.01),
       variables = rbind(
         do.call(rbind, Map(quantiles, draws, seq_along(draws))),
         quantiles(do.call(rbind, draws[kept]), "kept")))
}

# The labels of greta's variable elements, e.g. "init_country_level[12,7]":
# their indices replaced by the labels of the variable's dimensions
# (`labels`), leaving out dimensions of length one, and "" for a scalar.
label_elements <- function(elements, labels) {
  vapply(elements, function(element) {
    name <- sub("\\[.*$", "", element)
    if (name == element) return("")
    index <- sub("^.*\\[(.*)\\]$", "\\1", element)
    dimensions <- labels[[name]]
    if (is.null(dimensions)) return(index)
    index <- as.integer(strsplit(index, ",")[[1]])
    long <- lengths(dimensions) > 1
    paste(mapply(function(d, i) d[i], dimensions[long], index[long]),
          collapse = ", ")
  }, character(1), USE.NAMES = FALSE)
}

# The `n` elements with the highest Rhat in check_chains()'s `check`, as text.
worst_elements <- function(check, n = 5) {
  worst <- head(sort(check$rhat, decreasing = TRUE), n)
  paste(sprintf("%s %.3f", names(worst), worst), collapse = "; ")
}

# Lines for the fit's log, from check_chains(): the Rhat of the log density
# and the worst element's, and each chain's mean log density (sd); when
# flagged, a WARNING naming the worst elements.
chain_check_lines <- function(check) {
  chains <- check$chains
  lines <- sprintf(paste("chains: rank Rhat of the log density %.3f, worst",
                         "element %.3f, %d of %d above 1.01 | chain mean log",
                         "densities %s"),
                   check$lp_rhat, max(check$rhat, na.rm = TRUE),
                   sum(check$rhat > 1.01, na.rm = TRUE), length(check$rhat),
                   paste(sprintf("%.0f (sd %.0f)%s", chains$lp_mean,
                                 chains$lp_sd,
                                 ifelse(chains$stuck, " stuck", "")),
                         collapse = ", "))
  if (check$flagged) {
    lines <- c(lines, paste("WARNING: chains disagree; worst elements:",
                            worst_elements(check)))
  }
  lines
}
