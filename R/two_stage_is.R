# Importance sampling from the posterior of the second stage's hyperparameters
# theta (#34): a check on the CCD integration of the final model
# (fit_correction(), ccd_points()), not part of it. The proposal starts from
# the moments of the CCD points and is adapted as in AMIS; the estimates come
# from a fresh batch from the fixed adapted proposal, weighted by the Laplace
# marginal posterior and Pareto-smoothed (PSIS), with Pareto k as the
# diagnostic (reliable below 0.7).
#
# API (source after R/two_stage_correction.R):
#
#   is <- theta_is(fit, n = 200, ...)    importance sampling for a fit:
#     theta, weight (resampled to equally weighted points, repeats merged),
#     hyper (the hyperparameters at each point), the draws, and diagnostics.
#     A fit with fit$ccd <- is predicts by the importance sampling posterior.

# log density of the multivariate t with df degrees of freedom, centre mu and
# scale matrix sigma, at the rows of x
log_dmvt <- function(x, mu, sigma, df) {
  d <- length(mu)
  R <- chol(sigma)
  delta <- colSums(backsolve(R, t(x) - mu, transpose = TRUE) ^ 2)
  lgamma((df + d) / 2) - lgamma(df / 2) - d / 2 * log(df * pi) -
    sum(log(diag(R))) - (df + d) / 2 * log1p(delta / df)
}

# The log marginal posterior at the rows of theta, in `cores` forked
# processes, each taking its share in a nearest-neighbour order, so that each
# inner solve starts from the mode at a nearby theta; -Inf where it fails
log_post_parallel <- function(obj, theta, covariance, cores) {
  shares <- split(seq_len(nrow(theta)), rep_len(seq_len(cores), nrow(theta)))
  parts <- parallel::mclapply(shares, function(rows) {
    rows <- rows[point_order(theta[rows, , drop = FALSE], covariance)]
    start <- Sys.time()
    out <- numeric(length(rows))
    for (r in seq_along(rows)) {
      f <- tryCatch(obj$fn(theta[rows[r], ]), error = function(e) NA_real_)
      out[r] <- if (is.finite(f)) -f else -Inf
      if (r %% 50 == 0) {
        report("is process %i: %i of %i draws, %.1f s per draw", Sys.getpid(),
               r, length(rows),
               as.numeric(Sys.time() - start, units = "secs") / r)
      }
    }
    list(rows = rows, log_post = out)
  }, mc.cores = cores, mc.preschedule = FALSE)
  if (!all(vapply(parts, is.list, logical(1)))) {
    stop("evaluating the marginal posterior failed in a forked process")
  }
  log_post <- numeric(nrow(theta))
  for (part in parts) log_post[part$rows] <- part$log_post
  log_post
}

# Importance sampling from the marginal posterior of theta, with a proposal
# adapted as in AMIS (Cornuet, Marin, Mira & Robert 2012), then a fresh batch
# from the fixed adapted proposal, so that the reported estimates are plain
# importance sampling.
#   adapt  `rounds` rounds of n_adapt / rounds draws, each from a multivariate
#          t with df degrees of freedom: the first with the weighted mean and
#          covariance of the CCD points (start), each later one with the
#          weighted mean and covariance of all the draws so far, weighted by
#          the exact marginal posterior density over the equal mixture of the
#          rounds' proposals;
#   final  n_final fresh draws from the equal mixture q of all rounds + 1
#          proposals (the first, the Laplace t, keeps the weights bounded
#          where the adapted ones fall short), weighted by p / q, the log
#          weights Pareto-smoothed (loo::psis()).
# Returns n equally weighted points, resampled systematically from the final
# batch (repeats merged, weight = share), and the diagnostics: Pareto k of
# the final batch (reliable below 0.7), its effective sample size, the ESS
# after each adaptation round, and the posterior mean, SD and Monte Carlo SE
# of each hyperparameter from the final batch
hyper_is <- function(obj, mode, start, n = 200, n_adapt = 1000,
                     n_final = 1000, rounds = 4, df = 3,
                     cores = as.integer(Sys.getenv("IR_CUBE_IS_CORES", "1"))) {
  d <- length(mode$mode)
  # the first proposal: the weighted moments of the CCD points (start: their
  # theta and weight), which follow the posterior's skew along each axis
  w0 <- start$weight / sum(start$weight)
  centre_0 <- colSums(start$theta * w0)
  scale_0 <- crossprod(sweep(start$theta, 2, centre_0) * sqrt(w0))
  # draws from component k of the proposals, and the log density of their
  # equal mixture
  draw_t <- function(n, centre, scale) {
    z <- matrix(rnorm(n * d), n) / sqrt(rchisq(n, df) / df)
    sweep(z %*% chol(scale), 2, centre, "+")
  }
  log_mixture <- function(x, centres, scales) {
    l <- sapply(seq_along(centres), function(k) {
      log_dmvt(x, centres[[k]], scales[[k]], df)
    })
    l <- matrix(l, nrow(x))
    apply(l, 1, function(v) max(v) + log(mean(exp(v - max(v)))))
  }
  normalise <- function(log_weight) {
    w <- exp(log_weight - max(log_weight))
    w / sum(w)
  }

  # adaptation
  per_round <- ceiling(n_adapt / rounds)
  centres <- list(centre_0)
  scales <- list(scale_0)
  theta <- matrix(numeric(0), 0, d)
  log_target <- numeric(0)
  round_ess <- numeric(rounds)
  for (r in seq_len(rounds)) {
    new <- draw_t(per_round, centres[[r]], scales[[r]])
    log_target <- c(log_target,
                    log_post_parallel(obj, new, scales[[r]], cores))
    theta <- rbind(theta, new)
    weight <- normalise(log_target - log_mixture(theta, centres, scales))
    round_ess[r] <- 1 / sum(weight ^ 2)
    report("is adaptation round %i: %i draws, ESS %.0f", r, nrow(theta),
           round_ess[r])
    centre <- colSums(theta * weight)
    covariance <- crossprod(sweep(theta, 2, centre) * sqrt(weight)) +
      diag(1e-6 * diag(mode$covariance), d)
    centres[[r + 1]] <- centre
    scales[[r + 1]] <- covariance
  }

  # the final batch, from the fixed mixture of all the proposals
  component <- sample(seq_along(centres), n_final, replace = TRUE)
  theta <- matrix(NA_real_, n_final, d)
  for (k in unique(component)) {
    rows <- which(component == k)
    theta[rows, ] <- draw_t(length(rows), centres[[k]], scales[[k]])
  }
  log_target <- log_post_parallel(obj, theta, scales[[length(scales)]], cores)
  log_weight <- log_target - log_mixture(theta, centres, scales)
  smoothed <- suppressWarnings(loo::psis(log_weight, r_eff = NA))
  weight <- as.vector(weights(smoothed, log = FALSE, normalize = TRUE))
  report("is final batch: %i draws, ESS %.0f, Pareto k %.2f", n_final,
         1 / sum(weight ^ 2), loo::pareto_k_values(smoothed))

  # summaries of each hyperparameter, with the Monte Carlo SE of a
  # self-normalised importance sampling mean
  centre <- colSums(theta * weight)
  deviation <- sweep(theta, 2, centre)
  summary <- data.frame(parameter = names(mode$mode) %||% seq_len(d),
                        mean = centre,
                        sd = sqrt(colSums(weight * deviation ^ 2)),
                        mcse = sqrt(colSums(weight ^ 2 * deviation ^ 2)))

  index <- ccd_allocation(weight, n)
  counts <- table(index)
  keep <- as.integer(names(counts))
  list(theta = theta[keep, , drop = FALSE],
       weight = as.vector(counts) / n,
       draws = theta, log_weight = log_weight, is_weight = weight,
       proposals = list(centres = centres, scales = scales, df = df),
       diagnostics = list(pareto_k = loo::pareto_k_values(smoothed),
                          ess = 1 / sum(weight ^ 2),
                          round_ess = round_ess,
                          n_failed = sum(!is.finite(log_target)),
                          summary = summary))
}

# Importance sampling for a fit, started from its CCD points, as fit$ccd is
# laid out (theta, weight, covariance, hyper), so that a copy of the fit with
# its ccd replaced predicts by it. cores: forked single-threaded processes
# evaluating the marginal posterior
theta_is <- function(fit, n = 200, ...) {
  obj <- correction_obj(fit)
  theta_hat <- setNames(hyper_vector(fit), names(obj$par))
  out <- hyper_is(obj, list(mode = theta_hat,
                            covariance = fit$ccd$covariance),
                  start = fit$ccd, n = n, ...)
  colnames(out$theta) <- names(theta_hat)
  out$covariance <- fit$ccd$covariance
  out$hyper <- bind_rows(lapply(seq_len(nrow(out$theta)), function(j) {
    as_tibble(correction_hyper(theta_par_list(fit, out$theta[j, ]),
                               obj$env$data))
  }))
  out
}
