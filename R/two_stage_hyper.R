# Integrating the two-stage correction over its hyperparameters (#34, point 4).
#
# fit_correction() estimates the hyperparameters theta (the SDs and ranges of
# omega and eta, phi, and the iid SDs) once, by penalised maximum marginal
# likelihood, and predicts from p(latent | theta_hat, y): empirical Bayes.
# Here predictions instead mix p(latent | theta_j, y) over points theta_j of
# p(theta | y), with weights w_j:
#
#   nuts     NUTS draws of theta from the marginal posterior (tmbstan's
#            model, in coordinates whitened by the Gaussian approximation),
#            equal weights: the reference;
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

# the hyperparameters, on the optimiser's scale, in the template's order
correction_hyper_names <- c("log_sigma_omega", "log_kappa_omega",
                            "log_sigma_eta", "log_kappa_eta", "logit_phi",
                            "log_sigma_iid")

# theta of a fit, on the optimiser's scale, named as obj$par
hyper_vector <- function(fit) {
  unlist(fit$par_list[correction_hyper_names]) |>
    setNames(rep(correction_hyper_names,
                 lengths(fit$par_list[correction_hyper_names])))
}

# The TMB object of the fit's marginal likelihood of theta, started at its
# latent mode and theta_hat: "pql" at its final working response
# (fit$tmb_data), "beta_binomial" for the counts
correction_marginal <- function(fit, likelihood = c("pql", "beta_binomial")) {
  likelihood <- match.arg(likelihood)
  data <- fit$tmb_data
  data$likelihood <- if (likelihood == "pql") 0L else 1L
  correction_adfun(data, fit$par_list)
}

# W, the negative second derivative in lambda of the beta-binomial log
# likelihood of `died` of n with mean ilogit(lambda) and intra-class
# correlation rho (a = c p, b = c (1 - p), c = (1 - rho) / rho)
beta_binomial_weight <- function(lambda, died, n, rho) {
  p <- plogis(lambda)
  s <- p * (1 - p)
  c <- (1 - rho) / rho
  a <- c * p
  b <- c * (1 - p)
  d1 <- c * (digamma(died + a) - digamma(n - died + b) - digamma(a) +
               digamma(b))
  d2 <- c ^ 2 * (trigamma(died + a) + trigamma(n - died + b) - trigamma(a) -
                   trigamma(b))
  -(d2 * s ^ 2 + d1 * s * (1 - 2 * p))
}

# The correction_fit at theta (optimiser scale): the latent mode solved at
# theta, its Hessian, and the factor of that, updated from fit's symbolic
# factorisation (the sparsity pattern does not depend on theta). For the
# beta-binomial marginal, precision_obs becomes W at the mode, for the
# linearised cut-posterior shift
conditional_fit <- function(fit, obj, theta) {
  obj$fn(theta)
  par_full <- obj$env$last.par
  random <- obj$env$random
  H <- obj$env$spHess(par_full, random = TRUE)
  H <- Matrix::forceSymmetric(as(H, "CsparseMatrix"), uplo = "L")
  par_list <- obj$env$parList(par = par_full)
  fit$mode <- setNames(par_full[random], names(fit$mode))
  fit$H <- H
  fit$H_chol <- Matrix::update(fit$H_chol, H)
  fit$par_list <- par_list
  fit$hyper <- correction_hyper(par_list)
  data <- obj$env$data
  if (data$likelihood == 1) {
    lambda <- fit$m_ref + as.vector(fit$A_latent %*% fit$mode)
    fit$precision_obs <- beta_binomial_weight(lambda, data$died,
                                              data$n_mosquito, data$rho)
  }
  fit
}

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
# at the mode of this marginal). n is the number of
# points for laplace and nuts (nuts thins its draws to n); ccd's number is
# set by the dimension
hyper_points <- function(fit, obj, method = c("mode", "laplace", "ccd",
                                              "nuts"),
                         n = 100, ...) {
  method <- match.arg(method)
  start <- Sys.time()
  mode <- hyper_mode(fit, obj)
  out <- switch(method,
                mode = list(theta = matrix(mode$mode, 1), weight = 1),
                laplace = hyper_laplace(mode, n),
                ccd = hyper_ccd(obj, mode, ...),
                nuts = hyper_nuts(obj, mode, n, ...))
  colnames(out$theta) <- names(mode$mode)
  c(list(method = method, mode = mode$mode, covariance = mode$covariance,
         mode_convergence = mode$convergence),
    out,
    list(seconds = as.numeric(Sys.time() - start, units = "secs")))
}

# n draws from N(mode, covariance), equal weights
hyper_laplace <- function(mode, n) {
  L <- t(chol(mode$covariance))
  z <- matrix(rnorm(n * length(mode$mode)), length(mode$mode))
  theta <- t(mode$mode + L %*% z)
  list(theta = theta, weight = rep(1 / n, n))
}

# The CCD design in d dimensions, as points x d in the standardised z scale:
# the centre, the 2d axial points at distance sqrt(d), and a two-level
# factorial of resolution at least V (full for d <= 4, the half fraction with
# the last factor the product of the others for d = 5-7), all scaled by f0
ccd_design <- function(d, f0 = 1.1) {
  stopifnot(d >= 2, d <= 7)
  full <- as.matrix(expand.grid(rep(list(c(-1, 1)), if (d <= 4) d else d - 1)))
  factorial <- if (d <= 4) full else cbind(full, apply(full, 1, prod))
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

# NUTS on the marginal posterior of theta, in whitened coordinates z, theta =
# mode + L z with L L' the covariance of the Gaussian approximation, so that
# the target is near a standard normal from the first iteration and warmup
# need not find the scales (tmbstan's model, given the marginal's fn and gr in
# z; the Jacobian of the linear map is constant). Chains start near the mode;
# the draws are thinned evenly to n points. Returns the full draws and the
# chains' diagnostics too
hyper_nuts <- function(obj, mode, n, chains = 4, warmup = 500, iter = 1000,
                       cores = chains, adapt_delta = 0.9, init_sd = 0.5,
                       max_treedepth = 8) {
  suppressMessages(library(tmbstan))
  d <- length(mode$mode)
  L <- t(chol(mode$covariance))
  to_theta <- function(z) as.vector(mode$mode + L %*% z)
  # tmbstan's sampling method reads the DLL's name from fn's environment
  DLL <- obj$env$DLL
  model <- tmbstan:::tmbstan_model(
    par = setNames(rep(0, d), rep("z", d)),
    fn = function(z) obj$fn(to_theta(z)),
    gr = function(z) as.vector(obj$gr(to_theta(z)) %*% L))
  inits <- lapply(seq_len(chains), function(i) {
    list(z = numeric(), y = rnorm(d, 0, init_sd))
  })
  sampled <- rstan::sampling(model, chains = chains, iter = iter,
                             warmup = warmup, init = inits, cores = cores,
                             seed = sample.int(1e6, 1),
                             control = list(adapt_delta = adapt_delta,
                                            max_treedepth = max_treedepth))
  z_draws <- as.matrix(sampled)
  z_draws <- z_draws[, setdiff(colnames(z_draws), "lp__"), drop = FALSE]
  stopifnot(ncol(z_draws) == d)
  draws <- t(apply(z_draws, 1, to_theta))
  colnames(draws) <- names(mode$mode)
  summary <- posterior::summarise_draws(
    posterior::as_draws_array(as.array(sampled)), "rhat", "ess_bulk",
    "ess_tail")
  sampler <- rstan::get_sampler_params(sampled, inc_warmup = FALSE)
  keep <- round(seq(1, nrow(draws), length.out = n))
  list(theta = unname(draws[keep, , drop = FALSE]), weight = rep(1 / n, n),
       draws = draws,
       diagnostics = list(
         summary = summary,
         divergent = sum(vapply(sampler, function(s) sum(s[, "divergent__"]),
                                numeric(1))),
         treedepth = unlist(lapply(sampler, function(s) s[, "treedepth__"])),
         n_leapfrog = unlist(lapply(sampler,
                                    function(s) s[, "n_leapfrog__"]))))
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
  for (j in used) {
    fit_j <- conditional_fit(fit, obj, hyper$theta[j, ])
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
  draws
}
