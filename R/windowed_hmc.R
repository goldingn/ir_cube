# HMC with Stan-style windowed warmup, as a subclass of greta's hmc sampler.
#
# greta's hmc() estimates the diagonal of the mass matrix (diag_sd) between
# 10% and 40% of warmup from every warmup draw so far, including the transient
# from the initial values, and fixes it from 40%. This estimates it in windows
# that each start afresh: after an initial buffer (step size only), the metric
# is estimated from the draws in each window alone, and at the end of each
# window diag_sd is set from them and step-size adaptation restarts. The last
# part of warmup adapts the step size only.
#
#   windows      fractions of warmup at which the metric windows end; the
#                first starts at `buffer`
#   accept_target  target mean acceptance probability of the step-size
#                  adaptation (greta's hmc() uses 0.5)

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
    metric_log = list(),

    update_welford = function() {
      # accumulated in tune(), which knows the iteration
      invisible(NULL)
    },

    accumulate = function() {
      trace_matrix <- do.call(rbind, self$last_burst_free_states)
      if (self$welford_state$count == 0) {
        self$welford_state <- list(count = 0,
                                   mean = rep(0, ncol(trace_matrix)),
                                   m2 = rep(0, ncol(trace_matrix)))
      }
      count <- self$welford_state$count
      mean <- self$welford_state$mean
      m2 <- self$welford_state$m2
      for (i in seq_len(nrow(trace_matrix))) {
        count <- count + 1
        delta <- trace_matrix[i, ] - mean
        mean <- mean + delta / count
        m2 <- m2 + delta * (trace_matrix[i, ] - mean)
      }
      self$welford_state <- list(count = count, mean = mean, m2 = m2)
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
        self$accumulate()
      }
      if (n_ended > self$n_windows_done) {
        self$n_windows_done <- n_ended
        n <- self$welford_state$count
        if (n > 10) {
          var <- self$welford_state$m2 / (n - 1)
          var <- (n / (n + 5)) * var + 1e-3 * (5 / (n + 5))
          # the step size of parameter i is epsilon * diag_sd_i / sum(diag_sd):
          # keep it the same multiple of diag_sd_i
          self$parameters$epsilon <- self$parameters$epsilon *
            sum(sqrt(var)) / sum(self$parameters$diag_sd)
          self$parameters$diag_sd <- sqrt(var)
          self$metric_log[[length(self$metric_log) + 1]] <-
            list(iteration = iterations_completed, n = n,
                 epsilon = self$parameters$epsilon, diag_sd = sqrt(var))
        }
        self$welford_state <- list(count = 0, mean = 0, m2 = 0)
        self$restart_step_size()
      }
      self$tune_step_size(iterations_completed == total_iterations)
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
    },

    n_windows_done = 0
  )
)
}

# The sampler object for mcmc(), as greta's hmc() builds it.
windowed_hmc <- function(Lmin = 15, Lmax = 30, epsilon = 0.1, diag_sd = 1,
                         buffer = 0.15, windows = c(0.25, 0.45, 0.9),
                         accept_target = 0.5) {
  obj <- list(parameters = list(Lmin = Lmin, Lmax = Lmax, epsilon = epsilon,
                                diag_sd = diag_sd, buffer = buffer,
                                windows = windows),
              name = "hmc",
              class = R6::R6Class("windowed_hmc_sampler_target",
                                  inherit = windowed_hmc_class(),
                                  public = list(accept_target =
                                                  accept_target)))
  class(obj) <- c("hmc sampler", "sampler")
  obj
}
