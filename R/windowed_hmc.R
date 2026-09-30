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
# The step size is adapted as greta's hmc() does, one for all chains from
# their mean acceptance, leaving out proposals whose acceptance is not finite.
# It can instead be adapted per chain (per_chain = TRUE), from each chain's
# own acceptance with a non-finite acceptance counted as a rejection. That
# freed a chain that moved once in 1,000 draws among 8 on a screening subset,
# but on the 2014 forecasting fold mixed worse (worst rank Rhat 1.120 against
# 1.049, one chain rejecting 35% of its proposals while sampling); counting
# non-finite acceptances as rejections with the step shared was also worse
# there (1.066), though runs of one setting vary by about as much.
#
#   windows        fractions of warmup at which the metric windows end; the
#                  first starts at `buffer`
#   accept_target  target mean acceptance probability of the step-size
#                  adaptation (greta's hmc() uses 0.5)
#   per_chain      TRUE to adapt a step size per chain

# The R6 class, built on first use so that sourcing this file does not need
# greta.
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
      chain_accept = NULL,
      metric_log = list(),

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
            var <- apply(draws, 2, stats::var)
            var <- (n / (n + 5)) * var + 1e-3 * (5 / (n + 5))
            # the step size of parameter i is epsilon * diag_sd_i /
            # sum(diag_sd): keep it the same multiple of diag_sd_i
            self$parameters$epsilon <- self$parameters$epsilon *
              sum(sqrt(var)) / sum(self$parameters$diag_sd)
            self$parameters$diag_sd <- sqrt(var)
            self$metric_log[[length(self$metric_log) + 1]] <-
              list(iteration = iterations_completed, n = n,
                   epsilon = self$parameters$epsilon)
          }
          self$restart_step_size()
        }
        self$tune_step_size(iterations_completed == total_iterations)
      },

      # dual averaging, per chain or shared
      tune_step_size = function(final) {
        kappa <- 0.75
        gamma <- 0.05
        t0 <- 10
        self$da_updates <- self$da_updates + 1
        t <- self$da_updates
        w1 <- 1 / (t + t0)
        accept <- if (isTRUE(self$parameters$per_chain)) self$chain_accept else
          self$mean_accept_stat
        self$hbar <- (1 - w1) * self$hbar +
          w1 * (self$accept_target - accept)
        log_epsilon <- self$da_mu - self$hbar * sqrt(t) / gamma
        w2 <- t^-kappa
        self$log_epsilon_bar <- w2 * log_epsilon +
          (1 - w2) * self$log_epsilon_bar
        self$parameters$epsilon <- exp(log_epsilon)
        if (final) {
          self$parameters$epsilon <- exp(self$log_epsilon_bar)
        }
      },

      sampler_parameter_values = function() {
        p <- self$parameters
        self$parameters$epsilon <- rep_len(p$epsilon, self$n_chains)
        list(hmc_l = sample(seq(p$Lmin, p$Lmax), 1),
             hmc_epsilon = matrix(self$parameters$epsilon),
             hmc_diag_sd = matrix(p$diag_sd))
      },

      define_tf_kernel = function(sampler_param_vec) {
        tf <- tensorflow::tf
        tfp <- greta:::tfp
        dag <- self$model$dag
        n <- as.integer(self$n_free)
        n_chains <- as.integer(self$n_chains)
        hmc_l <- sampler_param_vec[0]
        # tensorflow's R indexing here is 0-based and includes the end
        epsilon <- tf$reshape(sampler_param_vec[1:n_chains],
                              shape = list(n_chains, 1L))
        diag_sd <- tf$reshape(
          sampler_param_vec[(n_chains + 1):(n_chains + n)],
          shape = list(1L, n))
        step_sizes <- epsilon * diag_sd / tf$reduce_sum(diag_sd)
        tfp$mcmc$HamiltonianMonteCarlo(
          target_log_prob_fn = dag$tf_log_prob_function_adjusted,
          step_size = step_sizes,
          num_leapfrog_steps = hmc_l)
      },

      # greta's run_burst(), also recording the mean acceptance of each chain
      # with a non-finite acceptance counted as 0
      run_burst = function(n_samples, thin = 1L) {
        param_vec <- unlist(self$sampler_parameter_values())
        self$n_bursts <- self$n_bursts + 1L
        burst_seed <- c(self$seed, self$n_bursts)
        batch_results <- self$sample_carefully(
          free_state = self$free_state,
          sampler_burst_length = as.integer(n_samples),
          sampler_thin = as.integer(thin),
          sampler_param_vec = param_vec,
          sampler_seed = burst_seed)
        free_state_draws <- as.array(batch_results$all_states)
        if (length(dim(free_state_draws)) != 3) {
          dim(free_state_draws) <- c(1, dim(free_state_draws))
        }
        self$last_burst_free_states <- greta:::split_chains(free_state_draws)
        n_draws <- nrow(free_state_draws)
        if (n_draws > 0) {
          free_state <- free_state_draws[n_draws, , , drop = FALSE]
          dim(free_state) <- dim(free_state)[-1]
          self$free_state <- free_state
        }
        log_accept <- as.array(batch_results$trace$log_accept_ratio)
        is_accepted <- as.array(batch_results$trace$is_accepted)
        self$accept_history <- rbind(self$accept_history, is_accepted)
        accept <- pmin(1, exp(log_accept))
        # greta's: non-finite proposals left out
        self$mean_accept_stat <- mean(accept, na.rm = TRUE)
        accept[!is.finite(log_accept)] <- 0
        self$chain_accept <- colMeans(matrix(accept, ncol = self$n_chains))
        self$numerical_rejections <- self$numerical_rejections +
          sum(!is.finite(log_accept))
      }
    )
  )
}

# The sampler object for mcmc(), as greta's hmc() builds it.
windowed_hmc <- function(Lmin = 15, Lmax = 30, epsilon = 0.1, diag_sd = 1,
                         buffer = 0.15, windows = c(0.25, 0.45, 0.9),
                         accept_target = 0.5, per_chain = FALSE) {
  obj <- list(parameters = list(Lmin = Lmin, Lmax = Lmax, epsilon = epsilon,
                                diag_sd = diag_sd, buffer = buffer,
                                windows = windows, per_chain = per_chain),
              name = "hmc",
              class = R6::R6Class("windowed_hmc_sampler_target",
                                  inherit = windowed_hmc_class(),
                                  public = list(accept_target =
                                                  accept_target)))
  class(obj) <- c("hmc sampler", "sampler")
  obj
}
