# HMC with Stan-style windowed warmup, as a subclass of greta's hmc sampler.
#
# greta's hmc() estimates the diagonal of the mass matrix (diag_sd) between
# 10% and 40% of warmup from every warmup draw so far, including the transient
# from the initial values, and fixes it from 40%. This estimates it in windows
# that each start afresh: after an initial buffer (step size only), the metric
# is estimated from the draws in each window alone (all chains pooled), and at
# the end of each window the metric is set from them and step-size adaptation
# restarts. The last part of warmup adapts the step size only.
#
#   metric         "diag", a diagonal mass matrix, as hmc(), or "dense", a
#                  full one, which follows linear correlations between
#                  parameters: the sampler runs in the parameters
#                  premultiplied by the inverse Cholesky factor of the
#                  estimated posterior covariance
#   windows        fractions of warmup at which the metric windows end; the
#                  first starts at `buffer`
#   accept_target  target mean acceptance probability of the step-size
#                  adaptation (greta's hmc() uses 0.5)
#   shrink         dense only: weight on the diagonal of the estimated
#                  covariance, (1 - shrink) S + shrink diag(S), since S is
#                  estimated from few effective draws per parameter

# The R6 classes, built on first use so that sourcing this file does not need
# greta: windowed_hmc_class() with the diagonal metric, and
# dense_hmc_class() with the dense one, overriding the kernel.
windowed_hmc_class <- function() {
  greta_hmc_class <- greta:::hmc_sampler
  R6::R6Class(
    "windowed_hmc_sampler",
    inherit = greta_hmc_class,
    public = list(
      da_updates = 0,
      da_mu = 0,
      n_windows_done = 0,
      window_draws = list(),
      metric_log = list(),

      dense = function() {
        identical(self$parameters$metric, "dense")
      },

      update_welford = function() {
        # draws are collected in tune(), which knows the iteration
        invisible(NULL)
      },

      restart_step_size = function() {
        self$da_updates <- 0
        self$hbar <- 0
        self$log_epsilon_bar <- 0
        self$da_mu <- log(10 * self$parameters$epsilon)
      },

      tune = function(iterations_completed, total_iterations) {
        p <- self$parameters
        if (self$da_updates == 0 && self$n_windows_done == 0) {
          self$restart_step_size()
        }
        ends <- round(p$windows * total_iterations)
        n_ended <- sum(iterations_completed >= ends)
        if (iterations_completed > p$buffer * total_iterations &&
            self$n_windows_done < length(ends)) {
          self$window_draws[[length(self$window_draws) + 1]] <-
            do.call(rbind, self$last_burst_free_states)
        }
        if (n_ended > self$n_windows_done) {
          self$n_windows_done <- n_ended
          draws <- do.call(rbind, self$window_draws)
          self$window_draws <- list()
          n <- nrow(draws)
          if (!is.null(n) && n > 10) {
            self$set_metric(draws)
            self$metric_log[[length(self$metric_log) + 1]] <-
              list(iteration = iterations_completed, n = n,
                   epsilon = self$parameters$epsilon)
          }
          self$restart_step_size()
        }
        self$tune_step_size(iterations_completed == total_iterations)
      },

      set_metric = function(draws) {
        n <- nrow(draws)
        if (!self$dense()) {
          var <- apply(draws, 2, stats::var)
          var <- (n / (n + 5)) * var + 1e-3 * (5 / (n + 5))
          # the step size of parameter i is epsilon * diag_sd_i /
          # sum(diag_sd): keep it the same multiple of diag_sd_i
          self$parameters$epsilon <- self$parameters$epsilon *
            sum(sqrt(var)) / sum(self$parameters$diag_sd)
          self$parameters$diag_sd <- sqrt(var)
          return(invisible())
        }
        s <- stats::cov(draws)
        shrink <- self$parameters$shrink
        s <- (1 - shrink) * s + shrink * diag(diag(s))
        s <- (n / (n + 5)) * s + 1e-3 * (5 / (n + 5)) * diag(ncol(s))
        # the step size is in the whitened parameters, so it carries over
        self$parameters$chol <- t(chol(s))
      },

      tune_step_size = function(final) {
        kappa <- 0.75
        gamma <- 0.05
        t0 <- 10
        self$da_updates <- self$da_updates + 1
        t <- self$da_updates
        w1 <- 1 / (t + t0)
        self$hbar <- (1 - w1) * self$hbar +
          w1 * (self$accept_target - self$mean_accept_stat)
        log_epsilon <- self$da_mu - self$hbar * sqrt(t) / gamma
        w2 <- t^-kappa
        self$log_epsilon_bar <- w2 * log_epsilon +
          (1 - w2) * self$log_epsilon_bar
        self$parameters$epsilon <- exp(log_epsilon)
        if (final) {
          self$parameters$epsilon <- exp(self$log_epsilon_bar)
        }
      }
    )
  )
}

dense_hmc_class <- function() {
  windowed_class <- windowed_hmc_class()
  R6::R6Class(
    "dense_hmc_sampler",
    inherit = windowed_class,
    public = list(
      sampler_parameter_values = function() {
        if (is.null(self$parameters$chol)) {
          self$parameters$chol <- diag(self$n_free)
        }
        list(hmc_l = sample(seq(self$parameters$Lmin, self$parameters$Lmax),
                            1),
             hmc_epsilon = self$parameters$epsilon,
             # row-major, as tf$reshape reads it
             hmc_chol = matrix(c(t(self$parameters$chol))))
      },

      define_tf_kernel = function(sampler_param_vec) {
        tf <- tensorflow::tf
        tfp <- greta:::tfp
        dag <- self$model$dag
        n <- as.integer(self$n_free)
        hmc_l <- sampler_param_vec[0]
        hmc_epsilon <- sampler_param_vec[1]
        chol <- tf$reshape(sampler_param_vec[2:(1 + n * n)],
                           shape = list(n, n))
        inner <- tfp$mcmc$HamiltonianMonteCarlo(
          target_log_prob_fn = dag$tf_log_prob_function_adjusted,
          step_size = hmc_epsilon,
          num_leapfrog_steps = hmc_l)
        tfp$mcmc$TransformedTransitionKernel(
          inner_kernel = inner,
          bijector = tfp$bijectors$ScaleMatvecTriL(scale_tril = chol))
      },

      define_tf_draws = function(free_state, sampler_burst_length,
                                 sampler_thin, sampler_param_vec,
                                 sampler_seed) {
        tf <- tensorflow::tf
        tfp <- greta:::tfp
        sampler_kernel <- self$define_tf_kernel(sampler_param_vec)
        tfp$mcmc$sample_chain(
          num_results = tf$math$floordiv(sampler_burst_length, sampler_thin),
          current_state = free_state,
          kernel = sampler_kernel,
          # the acceptance statistics are those of the inner kernel
          trace_fn = function(current_state, kernel_results) {
            kernel_results$inner_results
          },
          num_burnin_steps = tf$constant(0L, dtype = tf$int32),
          num_steps_between_results = sampler_thin,
          parallel_iterations = 1L,
          seed = sampler_seed)
      }
    )
  )
}

# The sampler object for mcmc(), as greta's hmc() builds it.
windowed_hmc <- function(Lmin = 15, Lmax = 30, epsilon = 0.1, diag_sd = 1,
                         metric = c("diag", "dense"),
                         buffer = 0.15, windows = c(0.25, 0.45, 0.9),
                         accept_target = 0.5, shrink = 0.1) {
  metric <- match.arg(metric)
  parent_class <- if (metric == "dense") dense_hmc_class() else
    windowed_hmc_class()
  obj <- list(parameters = list(Lmin = Lmin, Lmax = Lmax, epsilon = epsilon,
                                diag_sd = diag_sd, metric = metric,
                                buffer = buffer, windows = windows,
                                shrink = shrink),
              name = "hmc",
              class = R6::R6Class("windowed_hmc_sampler_target",
                                  inherit = parent_class,
                                  public = list(accept_target =
                                                  accept_target)))
  class(obj) <- c("hmc sampler", "sampler")
  obj
}
