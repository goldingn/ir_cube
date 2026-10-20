# Integrating the two-stage correction over its hyperparameters (#34, point 4).
#
# fit_correction() estimates the hyperparameters theta (the SDs and ranges of
# omega and eta, phi, and the iid SDs) once, by penalised maximum marginal
# likelihood, and predicts from p(latent | theta_hat, y): empirical Bayes.
# Here predictions instead mix p(latent | theta_j, y) over points theta_j of
# p(theta | y), with weights w_j:
#
#   is       importance sampling: a mixture of multivariate t proposals,
#            adapted as in AMIS from the moments of the CCD points, then a fresh
#            batch from it weighted by the exact (Laplace) marginal
#            posterior, the weights Pareto-smoothed (PSIS, Vehtari et al.
#            2024), and resampled to equally weighted points: the
#            reference, with Monte Carlo standard errors;
#   laplace  draws from the Gaussian approximation to p(theta | y) at its mode,
#            on the optimiser's (log / logit) scale, equal weights;
#   ccd      INLA's central composite design (Rue, Martino & Chopin 2009, sec.
#            6.5) on the same scale, with INLA's asymmetric scaling of each axis,
#            weighted by the posterior density at each point.
#
# Each works on one of two marginal likelihoods of theta
# (correction_marginal()):
#
#   pql            that of the final PQL expansion of the fit: the working
#                  response z and variance v held at their values at
#                  theta_hat, so that given theta the latent posterior is
#                  exactly Gaussian, N(mode(theta), H(theta)^-1);
#   beta_binomial  TMB's Laplace approximation of the beta-binomial likelihood
#                  of the counts itself (#34, point 2): given theta, the latent
#                  mode by Newton iterations, started from the last mode (so
#                  from a nearby theta), and N(mode, H^-1) at it. The
#                  cut-posterior shift per dynamical draw is then linearised,
#                  H^-1 A' W (m_ref - m), with W the negative second derivative
#                  of the log likelihood in lambda at the mode.
#
# Both include the priors and the Jacobians of the log / logit transforms
# (tmb/two_stage_correction.cpp).
#
# API (source after R/two_stage_correction.R):
#
#   obj <- correction_marginal(fit)          the TMB object of fit's marginal
#   hyper <- hyper_points(fit, obj, method, ...)
#     list(method, theta (points x hyperparameters, optimiser scale),
#     weight (sums to 1), mode, covariance, diagnostics, seconds)
#   conditional_fit(fit, obj, theta)         the correction_fit at theta: its
#     mode, H, factor and hyperparameters replaced; at theta_hat it is fit
#   predict_correction_hyper(fit, hyper, new, m_draws_train, m_draws_new,
#                            n_draws)        as predict_correction(), mixing
#     over the points of hyper: draw d uses dynamical draw ((d - 1) mod K) + 1
#     and point j(d), the points allotted to the draws by their weights

# An order of the points of theta (rows) in which each is near the last, a
# greedy nearest-neighbour chain from the first in the metric of covariance,
# so that each inner solve starts from the mode at a nearby theta
point_order <- function(theta, covariance) {
  z <- theta %*% t(solve(chol(covariance)))
  n <- nrow(z)
  order <- 1
  left <- seq_len(n)[-1]
  while (length(left) > 0) {
    last <- z[order[length(order)], ]
    d <- colSums((t(z[left, , drop = FALSE]) - last) ^ 2)
    order <- c(order, left[which.min(d)])
    left <- left[-which.min(d)]
  }
  order
}

# The mode of the marginal posterior of theta and the inverse of its Hessian
# there, from fit's estimate (which is the mode of an earlier PQL expansion,
# so a few iterations)
hyper_mode <- function(fit, obj, control = list(eval.max = 200,
                                                iter.max = 100)) {
  opt <- nlminb(hyper_vector(fit), obj$fn, obj$gr, control = control)
  hessian <- optimHess(opt$par, obj$fn, obj$gr)
  hessian <- (hessian + t(hessian)) / 2
  list(mode = opt$par, log_post = -opt$objective, hessian = hessian,
       covariance = solve(hessian), convergence = opt$convergence,
       message = opt$message)
}

# points of p(theta | y) by `method`; see the header ("mode" is the plug-in
# at the mode of this marginal). n is the number of points for gaussian and
# is (is resamples its weighted draws to n); ccd's number is set by the
# dimension
hyper_points <- function(fit, obj, method = c("mode", "gaussian", "ccd",
                                              "is"),
                         n = 100, ...) {
  method <- match.arg(method)
  start <- Sys.time()
  mode <- hyper_mode(fit, obj)
  out <- switch(method,
                mode = list(theta = matrix(mode$mode, 1), weight = 1),
                gaussian = hyper_gaussian(mode, n),
                ccd = hyper_ccd(obj, mode, ...),
                is = hyper_is(obj, mode, n, ...))
  colnames(out$theta) <- names(mode$mode)
  c(list(method = method, mode = mode$mode, covariance = mode$covariance,
         mode_convergence = mode$convergence),
    out,
    list(seconds = as.numeric(Sys.time() - start, units = "secs")))
}

# n draws from N(mode, covariance), equal weights
hyper_gaussian <- function(mode, n) {
  L <- t(chol(mode$covariance))
  z <- matrix(rnorm(n * length(mode$mode)), length(mode$mode))
  theta <- t(mode$mode + L %*% z)
  list(theta = theta, weight = rep(1 / n, n))
}

# The CCD design in d dimensions, as points x d in the standardised z scale:
# the centre, the 2d axial points at distance sqrt(d), and a two-level
# factorial of resolution at least V (full for d <= 4; the half fraction with
# the last factor the product of the others for d = 5-7; for d = 8 the quarter
# fraction with G = ABCD and H = ABEF), all scaled by f0
ccd_design <- function(d, f0 = 1.1) {
  stopifnot(d >= 2, d <= 8)
  base <- if (d <= 4) d else if (d <= 7) d - 1 else d - 2
  full <- as.matrix(expand.grid(rep(list(c(-1, 1)), base)))
  factorial <- if (d <= 4) {
    full
  } else if (d <= 7) {
    cbind(full, apply(full, 1, prod))
  } else {
    cbind(full, full[, 1] * full[, 2] * full[, 3] * full[, 4],
          full[, 1] * full[, 2] * full[, 5] * full[, 6])
  }
  axial <- rbind(diag(sqrt(d), d), diag(-sqrt(d), d))
  unname(rbind(0, axial, factorial) * f0)
}

# INLA's CCD integration: the design in z, where theta = mode + V S^1/2 z for
# the eigen-decomposition V S V' of the covariance, with each axis stretched
# separately on its positive and negative sides so that the log density falls
# by 1 at z = +-sqrt(2) as it would for a Gaussian (INLA's skewness
# correction). Weights: the posterior density times the design weight, 1 at
# the centre and exp(d f0^2 / 2) / ((n_p - 1)(f0^2 - 1)) at the other n_p - 1
# points, which all lie at radius f0 sqrt(d): with these, if the posterior is
# Gaussian the weighted points have its mean and covariance
hyper_ccd <- function(obj, mode, f0 = 1.1) {
  d <- length(mode$mode)
  eigen_cov <- eigen(mode$covariance, symmetric = TRUE)
  to_theta <- function(z) {
    as.vector(mode$mode + eigen_cov$vectors %*% (sqrt(eigen_cov$values) * z))
  }
  log_post <- function(theta) -obj$fn(theta)
  stretch <- matrix(1, d, 2, dimnames = list(NULL, c("negative", "positive")))
  for (i in seq_len(d)) {
    for (side in 1:2) {
      z <- rep(0, d)
      z[i] <- c(-1, 1)[side] * sqrt(2)
      drop <- mode$log_post - log_post(to_theta(z))
      stretch[i, side] <- if (is.finite(drop) && drop > 0) sqrt(1 / drop) else 1
    }
  }
  design <- ccd_design(d, f0)
  z <- design * ifelse(design < 0, rep(stretch[, 1], each = nrow(design)),
                       rep(stretch[, 2], each = nrow(design)))
  theta <- t(apply(z, 1, to_theta))
  log_density <- numeric(nrow(theta))
  for (j in point_order(theta, mode$covariance)) {
    log_density[j] <- log_post(theta[j, ]) - mode$log_post
  }
  n_points <- nrow(design)
  design_weight <- c(1, rep(exp(d * f0 ^ 2 / 2) /
                              ((n_points - 1) * (f0 ^ 2 - 1)), n_points - 1))
  weight <- design_weight * exp(log_density)
  list(theta = theta, weight = weight / sum(weight),
       diagnostics = list(stretch = stretch, log_density = log_density,
                          design_weight = design_weight))
}

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
#          t with df degrees of freedom: the first centred on the mode with the
#          covariance of the Gaussian approximation, each later one on the
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
hyper_is <- function(obj, mode, n = 200, n_adapt = 1000, n_final = 1000,
                     rounds = 4, df = 3, start = NULL,
                     cores = as.integer(Sys.getenv("IR_CUBE_IS_CORES", "1"))) {
  d <- length(mode$mode)
  # the first proposal: the weighted moments of the CCD points (start: their
  # theta and weight) when given, which follow the posterior's skew along
  # each axis; otherwise the Gaussian approximation at the mode
  if (is.null(start)) {
    centre_0 <- mode$mode
    scale_0 <- mode$covariance
  } else {
    w0 <- start$weight / sum(start$weight)
    centre_0 <- colSums(start$theta * w0)
    scale_0 <- crossprod(sweep(start$theta, 2, centre_0) * sqrt(w0))
  }
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

  index <- hyper_allocation(weight, n)
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

# The point of p(theta | y) each of n_draws draws uses: the points allotted by
# systematic resampling of their weights, in point order, so that each point's
# draws are contiguous and its conditional fit is made once
hyper_allocation <- function(weight, n_draws) {
  u <- (seq_len(n_draws) - runif(1)) / n_draws
  findInterval(u, cumsum(weight) / sum(weight)) + 1
}

# The fields of correction_node_draws() mixed over the points of hyper: the
# draws of point j come from conditional_fit() at it, with the cut-posterior
# shift for their dynamical draws (m_draws_train, n_draws x n_obs, or NULL).
# Returns the fields as correction_node_draws() does, and the point of each
# draw
correction_node_draws_hyper <- function(fit, obj, hyper, years, n_draws,
                                        m_draws_train = NULL,
                                        point = hyper_allocation(hyper$weight,
                                                                 n_draws)) {
  out <- NULL
  for (j in unique(point)) {
    draws <- which(point == j)
    fit_j <- conditional_fit(fit, obj, hyper$theta[j, ])
    fields <- correction_node_draws(
      fit_j, years, length(draws),
      if (is.null(m_draws_train)) NULL else
        m_draws_train[draws, , drop = FALSE])
    fields$theta <- NULL
    if (is.null(out)) {
      out <- list(omega = matrix(NA_real_, nrow(fields$omega), n_draws),
                  xi = lapply(fields$xi, function(x) {
                    matrix(NA_real_, nrow(x), n_draws)
                  }))
    }
    out$omega[, draws] <- fields$omega
    for (year in names(fields$xi)) out$xi[[year]][, draws] <- fields$xi[[year]]
  }
  out$point <- point
  out
}

# Draws x points predictive draws of new assays, as predict_correction(), with
# the hyperparameters integrated over the points of hyper. Each point's fit is
# made once; its draws are generated in batches of batch_size. The fresh iid
# noise uses each draw's own SDs
predict_correction_hyper <- function(fit, hyper, new, m_draws_train = NULL,
                                     m_draws_new = NULL, n_draws = 1000,
                                     batch_size = 100,
                                     obj = correction_marginal(fit)) {
  stopifnot(is.null(m_draws_train) == is.null(m_draws_new))
  point <- hyper_allocation(hyper$weight, n_draws)
  dynamical <- if (is.null(m_draws_train)) NULL else
    ((seq_len(n_draws) - 1) %% nrow(m_draws_train)) + 1
  draws <- matrix(NA_real_, n_draws, nrow(new))
  used <- unique(point)
  if (length(used) > 2) {
    used <- used[point_order(hyper$theta[used, , drop = FALSE],
                             hyper$covariance)]
  }
  # the conditional posterior mean of omega + xi at `new` at each point used,
  # for the Monte Carlo error of the integration over theta
  point_means <- matrix(NA_real_, length(used), nrow(new))
  for (i in seq_along(used)) {
    j <- used[i]
    fit_j <- conditional_fit(fit, obj, hyper$theta[j, ])
    point_means[i, ] <- project_correction(
      fit_j, correction_node_draws(fit_j, new$year, mean = TRUE), new)[, 1]
    rows <- which(point == j)
    for (batch in split(rows, ceiling(seq_along(rows) / batch_size))) {
      if (is.null(dynamical)) {
        fields <- correction_node_draws(fit_j, new$year, length(batch))
        m_new <- new$m
      } else {
        k <- dynamical[batch]
        fields <- correction_node_draws(fit_j, new$year, length(batch),
                                        m_draws_train[k, , drop = FALSE])
        m_new <- t(m_draws_new[k, , drop = FALSE])
      }
      draws[batch, ] <- t(m_new + project_correction(fit_j, fields, new,
                                                     noise = TRUE))
    }
  }
  attr(draws, "point") <- point
  attr(draws, "point_means") <- point_means
  attr(draws, "point_share") <- as.vector(table(factor(point, used))) /
    n_draws
  draws
}

# The Monte Carlo SE, from the integration over theta, of the predictive mean
# of omega + xi at each point of `new`, from predict_correction_hyper()'s
# draws: the weighted SD of the points' conditional means over the square
# root of the points' effective number (1 / sum of squared shares), or of
# n_eff when given (e.g. the importance sampling ESS, if smaller)
hyper_mcse <- function(draws, n_eff = NULL) {
  g <- attr(draws, "point_means")
  w <- attr(draws, "point_share")
  if (is.null(g) || nrow(g) < 2) return(rep(NA_real_, ncol(draws)))
  g_bar <- colSums(g * w)
  sd_w <- sqrt(colSums(w * sweep(g, 2, g_bar) ^ 2))
  sd_w / sqrt(min(c(1 / sum(w ^ 2), n_eff)))
}
