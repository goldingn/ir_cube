# Checks on R/validation_functions.R, run against simulated data where the
# truth is known. These verify that the beta-binomial primitives are correct,
# and that the calibration metrics behave as intended: flat for a model whose
# predictive distribution is the data-generating one, and deviating in the
# expected direction for models that are biased, overconfident or
# underconfident.

source("R/validation_functions.R")

set.seed(2026 - 8 - 31)

report <- function(label, pass, detail = "") {
  cat(sprintf("%-52s %s %s\n", label, ifelse(pass, "PASS", "FAIL"), detail))
  invisible(pass)
}


# beta-binomial primitives ------------------------------------------------

size <- 40
p <- 0.7
rho <- 0.12
support <- 0:size
pmf <- dbetabinom(support, size, p, rho)

report("pmf sums to one", isTRUE(all.equal(sum(pmf), 1)))
report("pmf mean matches size * p",
       isTRUE(all.equal(sum(support * pmf), size * p)))

# variance of a beta-binomial is n p (1 - p) [1 + (n - 1) rho]
variance_expected <- size * p * (1 - p) * (1 + (size - 1) * rho)
variance_observed <- sum((support - size * p) ^ 2 * pmf)
report("pmf variance matches n p q [1 + (n-1) rho]",
       isTRUE(all.equal(variance_observed, variance_expected)))

# as rho tends to zero the beta-binomial tends to the binomial
report("binomial limit as rho -> 0",
       isTRUE(all.equal(dbetabinom(support, size, p, 1e-9),
                        dbinom(support, size, p),
                        tolerance = 1e-5)))

report("cdf reaches one", isTRUE(all.equal(pbetabinom(size, size, p, rho), 1)))
report("cdf matches cumulative pmf",
       isTRUE(all.equal(pbetabinom(0:size, size, p, rho), cumsum(pmf))))

# sampled moments match
draws <- rbetabinom(2e5, size, p, rho)
report("rbetabinom mean", abs(mean(draws) - size * p) < 0.05,
       sprintf("(%.3f vs %.3f)", mean(draws), size * p))
report("rbetabinom variance",
       abs(var(draws) / variance_expected - 1) < 0.05,
       sprintf("(%.1f vs %.1f)", var(draws), variance_expected))

# CRPS against a brute force evaluation of the double sum
sims_test <- rnorm(300)
brute <- mean(abs(sims_test - 0.4)) -
  0.5 * mean(abs(outer(sims_test, sims_test, "-")))
report("crps_sample matches brute force double sum",
       isTRUE(all.equal(crps_sample(0.4, sims_test), brute)))


# behaviour of the calibration metrics ------------------------------------

n_obs <- 500
n_draws <- 300
rho_true <- 0.12
mosquito_number <- rep(100, n_obs)

# true population fractions, spanning the range seen in the data
p_true <- rbeta(n_obs, 6, 2)
died <- rbetabinom(n_obs, mosquito_number, p_true, rho_true)

# four candidate models, expressed as posterior draws. The first has the
# data-generating distribution as its predictive distribution
make_draws <- function(p, rho) {
  list(p = matrix(p, nrow = n_draws, ncol = n_obs, byrow = TRUE),
       rho = matrix(rho, nrow = n_draws, ncol = n_obs))
}

candidates <- list(
  calibrated = make_draws(p_true, rho_true),
  biased = make_draws(pmin(p_true + 0.08, 0.999), rho_true),
  overconfident = make_draws(p_true, rho_true / 4),
  underconfident = make_draws(p_true, rho_true * 4)
)

results <- lapply(
  candidates,
  function(candidate) {
    summary <- ppd_summary(died, mosquito_number, candidate$p, candidate$rho)
    pit <- ppd_pit(summary, n_rep = 20)
    coverage <- coverage_curve(pit, levels = c(0.5, 0.95))
    data.frame(
      mean_pit = mean(pit),
      coverage_50 = coverage$empirical[1],
      coverage_95 = coverage$empirical[2],
      cvm = pit_statistic(pit),
      elpd = mean(summary$log_score),
      mse = mean((summary$observed - summary$predicted) ^ 2)
    )
  }
)
results <- do.call(rbind, results)

cat("\n")
print(round(results, 4))
cat("\n")

null_band <- pit_null_band(n_obs, n_sim = 500)

report("calibrated: mean PIT near 0.5",
       abs(results["calibrated", "mean_pit"] - 0.5) < 0.03)
report("calibrated: 95% coverage near nominal",
       abs(results["calibrated", "coverage_95"] - 0.95) < 0.03)
report("calibrated: 50% coverage near nominal",
       abs(results["calibrated", "coverage_50"] - 0.5) < 0.06)
report("calibrated: CvM within null band",
       results["calibrated", "cvm"] < null_band[3],
       sprintf("(%.3f vs upper %.3f)", results["calibrated", "cvm"], null_band[3]))
report("biased: mean PIT shifted",
       abs(results["biased", "mean_pit"] - 0.5) > 0.05)
report("biased: CvM outside null band",
       results["biased", "cvm"] > null_band[3])
report("overconfident: coverage below nominal",
       results["overconfident", "coverage_95"] < 0.9)
report("underconfident: coverage above nominal",
       results["underconfident", "coverage_95"] > 0.98)
report("calibrated model has the best elpd",
       which.max(results$elpd) == 1)

# dispersion errors are invisible to a point summary: the biased model is the
# only one a mean squared error comparison can distinguish, which is the reason
# for scoring the distribution rather than the point
report("MSE cannot separate the dispersion errors",
       isTRUE(all.equal(results["overconfident", "mse"],
                        results["underconfident", "mse"],
                        tolerance = 1e-8)))


# the noise floor ----------------------------------------------------------

# with predictions equal to the truth, the mean squared error is entirely
# irreducible, so the floor estimator should recover it
mse_at_truth <- mean((died / mosquito_number - p_true) ^ 2)
floor_estimate <- noise_floor_mse(died, mosquito_number, rho_true)
report("noise floor recovers irreducible MSE",
       abs(floor_estimate / mse_at_truth - 1) < 0.1,
       sprintf("(%.5f vs %.5f)", floor_estimate, mse_at_truth))

# and the oracle CRPS should match the CRPS achieved by the calibrated model
sims <- ppd_simulate(mosquito_number, candidates$calibrated$p,
                     candidates$calibrated$rho)
crps_calibrated <- mean(ppd_crps(died, mosquito_number, sims))
crps_floor <- noise_floor_crps(mosquito_number, p_true, rho_true)
report("oracle CRPS matches calibrated model CRPS",
       abs(crps_floor / crps_calibrated - 1) < 0.05,
       sprintf("(%.4f vs %.4f)", crps_floor, crps_calibrated))


# reliability and aggregation ----------------------------------------------

reliability <- reliability_bins(colMeans(candidates$calibrated$p),
                                died / mosquito_number)
report("calibrated model is reliable across bins",
       max(abs(reliability$predicted - reliability$observed)) < 0.05,
       sprintf("(max deviation %.3f)",
               max(abs(reliability$predicted - reliability$observed))))

# pooled groups of assays: the group-level PIT should also be uniform
group <- rep(seq_len(n_obs / 5), each = 5)
aggregated <- ppd_aggregate(died, mosquito_number, group, sims)
report("aggregated predictive is calibrated",
       abs(mean(aggregated$pit) - 0.5) < 0.06,
       sprintf("(mean group PIT %.3f)", mean(aggregated$pit)))

# and pooling should shrink the scatter around the truth
scatter_single <- sd(died / mosquito_number - p_true)
scatter_pooled <- sd(aggregated$observed - aggregated$predicted)
report("pooling reduces scatter",
       scatter_pooled < scatter_single,
       sprintf("(%.4f vs %.4f, ratio %.2f)",
               scatter_pooled, scatter_single, scatter_pooled / scatter_single))


# agreement with the reference implementations ------------------------------

# The primitives above are written in base R so that the validation code
# carries no extra dependency. Where the reference packages are available,
# check that they agree.

if (requireNamespace("extraDistr", quietly = TRUE)) {

  agreement <- vapply(
    seq_along(c(0.2, 0.5, 0.9, 0.97)),
    function(i) {
      p <- c(0.2, 0.5, 0.9, 0.97)[i]
      rho <- c(0.05, 0.155, 0.3, 0.12)[i]
      a <- p * (1 / rho - 1)
      b <- a * (1 - p) / p
      density_matches <- all.equal(dbetabinom(0:100, 100, p, rho),
                                   extraDistr::dbbinom(0:100, 100, a, b))
      cdf_matches <- all.equal(pbetabinom(0:100, 100, p, rho),
                               extraDistr::pbbinom(0:100, 100, a, b),
                               tolerance = 1e-7)
      isTRUE(density_matches) && isTRUE(cdf_matches)
    },
    logical(1)
  )

  report("beta-binomial primitives agree with extraDistr", all(agreement))

}

if (requireNamespace("scoringRules", quietly = TRUE)) {

  crps_sims <- matrix(rnorm(2000), nrow = 500)
  crps_y <- c(0.3, -0.2, 1.1, 0.0)
  report(
    "crps_sample agrees with scoringRules",
    isTRUE(all.equal(
      vapply(seq_along(crps_y),
             function(i) crps_sample(crps_y[i], crps_sims[, i]),
             numeric(1)),
      scoringRules::crps_sample(y = crps_y, dat = t(crps_sims))
    ))
  )

}
