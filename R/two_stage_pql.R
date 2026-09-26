# Stage B of the two-stage correction (#21, "Fitting cascade", step 3):
# penalised quasi-likelihood (PQL) on the beta-binomial counts, started from
# the stage-A fit of R/two_stage_correction.R.
#
# Stage A treats the empirical logit as Gaussian with a fixed variance. That
# is worst at 0% and 100% mortality, where the +0.5 continuity correction pins
# z and the fitted latent is not extreme enough (doc/two_stage_plan.md, "Stage
# B diagnostic"). PQL replaces (z, v) by the working response and variance of
# the quasi-binomial likelihood, re-expanded around the current latent mode
# lambda = m + omega + xi + u, p = logit^-1(lambda):
#
#   z_i = lambda_i + (y_i / n_i - p_i) / (p_i (1 - p_i))
#   v_i = (1 + (n_i - 1) rho) / (n_i p_i (1 - p_i))
#
# With the hyperparameters fixed, one pass is one Newton step on the penalised
# quasi-log-likelihood
#
#   sum_i [y_i lambda_i - n_i log(1 + e^lambda_i)] / (1 + (n_i - 1) rho)
#     - theta' Q theta / 2,
#
# which is concave, so a step that fails to increase it is halved (damping,
# recorded). Only D = diag(1 / v) changes between passes, so the Hessian
# H = Q + A' D A keeps its sparsity pattern and the stage-A symbolic Cholesky
# factorisation is reused (Matrix::update on the CHMfactor). Q is never formed:
# H_new = H_A + A' diag(D_new - D_A) A.
#
# Passes stop when max |change in lambda| at the data < 0.01. If lambda has
# then moved materially from stage A (RMS change at the data > 0.1), the
# hyperparameters are re-estimated once, by the stage-A marginal likelihood
# with the final working (z, v) as data, and the PQL passes are re-run at the
# new hyperparameters.
#
# The result is a correction_fit with mode, H, H_chol, precision_obs (= the
# final D) and the hyperparameters replaced, so predict_correction() uses it
# unchanged: joint latent draws from H, the AR(1) forecast, and the
# cut-posterior shift H^-1 A' D (m_ref - m_k) at the final D.
#
# Functions only; source R/two_stage_correction.R first.

stopifnot(exists("fit_correction"), exists("predict_correction"))

# log(1 + e^x) without overflow
log1pexp <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

# PQL working response and variance at the latent lambda. lambda is clamped to
# +-clamp for the expansion only, so that p (1 - p) stays representable; in the
# folds lambda never gets near it (stage A: |lambda| < 7), and the count of
# clamped assays is recorded
pql_working <- function(lambda, died, n, rho, clamp = 15) {
  lambda_c <- pmin(pmax(lambda, -clamp), clamp)
  p <- plogis(lambda_c)
  pq <- p * (1 - p)
  phi <- 1 + (n - 1) * rho
  list(z = lambda_c + (died / n - p) / pq,
       v = phi / (n * pq),
       D = n * pq / phi,
       n_clamped = sum(lambda != lambda_c))
}

# penalised quasi-log-likelihood (up to a constant), with theta' Q theta
# computed as theta' H_base theta - sum D_base (A theta)^2
pql_objective <- function(theta, eta, H_base, D_base, died, n, rho) {
  lambda_part <- sum((died * eta$lambda - n * log1pexp(eta$lambda)) /
                       (1 + (n - 1) * rho))
  quad <- sum(theta * as.vector(H_base %*% theta)) -
    sum(D_base * eta$a_theta ^ 2)
  lambda_part - quad / 2
}

# PQL passes at the hyperparameters of fit (a correction_fit), from its mode.
# died, n and rho are per training assay, in the fit's order. Returns the final
# mode, the Hessian and its factor at the final D, the final working values,
# and one row per pass
pql_passes <- function(fit, died, n, rho,
                       tol = 0.01, max_passes = 30, max_halvings = 10,
                       clamp = 15, verbose = FALSE) {
  A <- fit$A_latent
  m <- fit$m_ref
  H_base <- fit$H
  D_base <- fit$precision_obs
  chol <- fit$H_chol
  theta <- as.vector(fit$mode)

  state <- function(theta) {
    a_theta <- as.vector(A %*% theta)
    list(a_theta = a_theta, lambda = m + a_theta)
  }
  objective <- function(theta, eta) {
    pql_objective(theta, eta, H_base, D_base, died, n, rho)
  }
  # the Hessian at working weights D, factorised with the stage-A symbolic
  # analysis, and the Newton target theta = H^-1 A' D (z - m)
  newton <- function(eta) {
    w <- pql_working(eta$lambda, died, n, rho, clamp)
    H <- Matrix::forceSymmetric(
      H_base + Matrix::crossprod(A, Matrix::Diagonal(x = w$D - D_base) %*% A),
      uplo = "L")
    chol <- Matrix::update(chol, H)
    rhs <- Matrix::crossprod(A, w$D * (w$z - m))
    list(w = w, H = H, chol = chol,
         target = as.vector(Matrix::solve(chol, rhs, system = "A")))
  }

  eta <- state(theta)
  f <- objective(theta, eta)
  history <- list()
  converged <- FALSE
  for (pass in seq_len(max_passes)) {
    time <- system.time({
      step_target <- newton(eta)$target
      step <- 1
      halvings <- 0
      repeat {
        theta_try <- theta + step * (step_target - theta)
        eta_try <- state(theta_try)
        f_try <- objective(theta_try, eta_try)
        if (f_try >= f - 1e-10 * abs(f) || halvings >= max_halvings) break
        step <- step / 2
        halvings <- halvings + 1
      }
    })[["elapsed"]]
    change <- eta_try$lambda - eta$lambda
    theta <- theta_try
    eta <- eta_try
    f <- f_try
    history[[pass]] <- data.frame(pass = pass, max_change = max(abs(change)),
                                  rms_change = sqrt(mean(change ^ 2)),
                                  step = step, objective = f, seconds = time)
    if (verbose) {
      message(sprintf("PQL pass %2i: max |d lambda| %.4f, step %.3g, %.1f s",
                      pass, max(abs(change)), step, time))
    }
    if (max(abs(change)) < tol) {
      converged <- TRUE
      break
    }
  }

  # the final expansion: z, v and H at the final mode, and the mode solving
  # H theta = A' D (z - m) for them, so that the cut-posterior shift and the
  # latent draws are consistent with the saved z and D
  final <- newton(eta)
  final_change <- max(abs(as.vector(A %*% (final$target - theta))))
  theta <- final$target
  names(theta) <- names(fit$mode)

  list(mode = theta,
       lambda = m + as.vector(A %*% theta),
       H = final$H,
       H_chol = final$chol,
       working = final$w,
       history = do.call(rbind, history),
       passes = pass,
       converged = converged,
       final_change = final_change,
       damped = any(vapply(history, function(h) h$step < 1, logical(1))))
}

# Stage B from a stage-A fit. train is the data frame the stage-A fit was
# given (with died, mosquito_number and rho, in the same order); t0 and T as
# for fit_correction(). Returns the stage-B correction_fit, with $stage_b
# holding the diagnostics
fit_correction_pql <- function(fit_a, train,
                               tol = 0.01,
                               max_passes = 30,
                               refit_threshold = 0.1,
                               clamp = 15,
                               verbose = FALSE) {
  start_time <- Sys.time()
  died <- train$died
  n <- train$mosquito_number
  rho <- train$rho
  stopifnot(length(died) == fit_a$n_obs)

  lambda_a <- fit_a$m_ref + as.vector(fit_a$A_latent %*% fit_a$mode)

  first <- pql_passes(fit_a, died, n, rho, tol = tol, max_passes = max_passes,
                      clamp = clamp, verbose = verbose)
  rms_move <- sqrt(mean((first$lambda - lambda_a) ^ 2))
  time_first <- difftime(Sys.time(), start_time, units = "secs")

  # re-estimate the hyperparameters once if lambda moved materially: the
  # stage-A marginal likelihood with the working (z, v) as the response,
  # started at the stage-A hyperparameters, then PQL again at the new ones
  refit <- rms_move > refit_threshold
  base <- fit_a
  final <- first
  time_refit <- 0
  time_second <- 0
  if (refit) {
    hyper_names <- c("log_sigma_omega", "log_kappa_omega", "log_sigma_eta",
                     "log_kappa_eta", "logit_phi", "log_tau")
    start <- lapply(fit_a$par_list[hyper_names], as.numeric)
    train_w <- train
    train_w$z <- first$working$z
    train_w$v <- first$working$v
    t_refit <- Sys.time()
    base <- fit_correction(train_w, variant = fit_a$variant, t0 = fit_a$t0,
                           T = fit_a$T, mesh = fit_a$mesh,
                           mesh_xi = fit_a$mesh_xi, start = start)
    time_refit <- difftime(Sys.time(), t_refit, units = "secs")
    t_second <- Sys.time()
    final <- pql_passes(base, died, n, rho, tol = tol,
                        max_passes = max_passes, clamp = clamp,
                        verbose = verbose)
    time_second <- difftime(Sys.time(), t_second, units = "secs")
  }

  out <- base
  out$mode <- final$mode
  out$H <- final$H
  out$H_chol <- final$H_chol
  out$precision_obs <- final$working$D
  out$tmb_data$z <- final$working$z
  out$tmb_data$v <- final$working$v
  out$stage_b <- list(
    lambda_a = lambda_a,
    lambda = final$lambda,
    passes_first = first$passes,
    converged_first = first$converged,
    damped_first = first$damped,
    history_first = first$history,
    rms_move = rms_move,
    max_move = max(abs(first$lambda - lambda_a)),
    refit = refit,
    refit_threshold = refit_threshold,
    refit_convergence = if (refit) base$opt$convergence else NA_integer_,
    hyper_a = fit_a$hyper,
    passes_second = if (refit) final$passes else NA_integer_,
    converged_second = if (refit) final$converged else NA,
    damped_second = if (refit) final$damped else NA,
    history_second = if (refit) final$history else NULL,
    rms_move_second = if (refit)
      sqrt(mean((final$lambda - first$lambda) ^ 2)) else NA_real_,
    final_change = final$final_change,
    n_clamped = final$working$n_clamped,
    timings = c(passes_first = as.numeric(time_first),
                refit = as.numeric(time_refit),
                passes_second = as.numeric(time_second),
                total = as.numeric(difftime(Sys.time(), start_time,
                                            units = "secs")))
  )
  out
}
