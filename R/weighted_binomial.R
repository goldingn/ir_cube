# The weighted binomial likelihood of the dynamical model (#47;
# dynamical_model_options(likelihood = "weighted_binomial")): the binomial log
# likelihood of each bioassay, weighted by its design effect,
#   log L_i = w_i [y_i log p_i + (n_i - y_i) log(1 - p_i)],
#   w_i = 1 / (1 + (n_i - 1) rho_type(i)),
# without the binomial coefficient, which is constant. rho_type is fixed at
# the intra-bioassay correlation estimated from replicate bioassays, so it
# cannot absorb misfit of the model. This is the quasi-binomial: each
# bioassay counts as n_i w_i independent mosquitoes. Unlike the beta-binomial
# with an estimated rho, its fitted mean is unbiased whatever spread between
# pixel-years the model does not explain (simulations of 8 October 2026). It
# is not a normalised density of the data, so its log likelihood is not
# comparable with a beta-binomial fit's, and data are simulated from the
# beta-binomial at the replicate rho instead.
# Functions and settings only. Sourced by R/dynamical_model.R.


# the replicate rho ---------------------------------------------------------

# rho per insecticide type, fixed. The file is a copy (8 October 2026) of
# outputs/bioassay_rho_hierarchical.csv, written by
# R/fig_illustrate_bioassay_variability.R: the posterior mean rho per type
# of a beta-binomial fitted by MCMC to groups of replicate bioassays (the same
# cell, year and type), with rho per type nested in class (columns
# insecticide_type, insecticide_class, rho, rho_lower, rho_upper, n_groups,
# n_assays, worst_rhat).
replicate_rho_file <- "data/clean/bioassay_rho_replicate.csv"

# The replicate rho of each of `types`, named, in their order
replicate_rho <- function(types, file = replicate_rho_file) {
  table <- utils::read.csv(file)
  missing <- setdiff(types, table$insecticide_type)
  if (length(missing) > 0) {
    stop("no replicate rho in ", file, " for ", toString(missing))
  }
  if (!all(table$worst_rhat < 1.05)) {
    stop("the replicate rho fit in ", file, " has not converged")
  }
  rho <- table$rho[match(types, table$insecticide_type)]
  stopifnot(all(rho > 0 & rho < 1))
  setNames(rho, types)
}

# The design-effect weight of bioassays of n mosquitoes with intra-bioassay
# correlation rho, 1 / (1 + (n - 1) rho): the variance of died / n is
# p (1 - p) (1 + (n - 1) rho) / n, that of n w independent mosquitoes
design_effect_weight <- function(n, rho) {
  1 / (1 + (n - 1) * rho)
}


# log p and log(1 - p), without losing precision ------------------------------

# log p and log(1 - p) for bioassay mortality p = f + (1 - f) q, q = ilogit(l),
# for a floor f (NULL for none), as list(log_p, log_not_p). Computed from the
# logit l, so that neither loses precision where p is near 0 or 1:
#   log q = -softplus(-l), log(1 - q) = -softplus(l)
# (softplus(x) = log(1 + e^x), stable in both tails). So
#   log(1 - p) = log(1 - f) - softplus(l)
# never takes 1 minus a number near 1, which rounds to 0 (and log(1 - p) to
# -Inf) for l above about 37 in double precision; and without a floor,
#   log p = -softplus(-l)
# never takes the log of a q that has underflowed. With a floor p >= f, so
# log(f + (1 - f) q) is safe as it is. For greta arrays or plain R
# (conformable l and f).
floored_log_probs <- function(l, floor = NULL) {
  softplus <- if (inherits(l, "greta_array")) greta::log1pe else
    function(x) -stats::plogis(-x, log.p = TRUE)
  log_q <- -softplus(-l)
  log_not_q <- -softplus(l)
  if (is.null(floor)) {
    return(list(log_p = log_q, log_not_p = log_not_q))
  }
  list(log_p = log(floor + (1 - floor) * exp(log_q)),
       log_not_p = log(1 - floor) + log_not_q)
}

# log p and log(1 - p) for the mixture p = s p_a + (1 - s) p_b of two
# trajectories (the species model: s the arabiensis share, a arabiensis, b
# the other members of the complex), from theirs (floored_log_probs()), as
# list(log_p, log_not_p). 1 - p is the mixture s (1 - p_a) + (1 - s)(1 - p_b)
# of the trajectories' own 1 - p: sums of non-negative terms, so neither
# cancels, and each 1 - p comes from its log, never from 1 minus p. Only
# where both trajectories' p (or 1 - p) are below 1e-308 does the log
# underflow. For greta arrays or plain R.
mixture_log_probs <- function(share, a, b) {
  list(log_p = log(share * exp(a$log_p) + (1 - share) * exp(b$log_p)),
       log_not_p = log(share * exp(a$log_not_p) +
                         (1 - share) * exp(b$log_not_p)))
}


# the likelihood --------------------------------------------------------------

# The weighted binomial log likelihood of `died` of `n` (in plain R), from log
# p and log(1 - p) (floored_log_probs()) and the weights w
# (design_effect_weight()): w [died log p + (n - died) log(1 - p)], with
# 0 log 0 = 0, as the limit of the binomial's
weighted_binomial_log_lik <- function(died, n, log_p, log_not_p, weight) {
  died_term <- ifelse(died == 0, 0, died * log_p)
  survived_term <- ifelse(n - died == 0, 0, (n - died) * log_not_p)
  weight * (died_term + survived_term)
}

# The weighted binomial as a greta distribution, for
#   distribution(died) <- weighted_binomial(n, log_p, log_not_p, weight)
# with n and weight data and log_p and log_not_p greta arrays
# (floored_log_probs()), each one per bioassay.
#
# greta has no weighted likelihood, and no way to add an arbitrary term to the
# log density: every term comes from a distribution, assigned to data with
# distribution(). So this is a distribution class of its own, built as greta
# builds its own: an R6 class inheriting greta's distribution_node, whose
# tf_distrib() returns the log density function, as
# greta:::beta_binomial_distribution's does. It uses greta internals
# (distribution_node, as.greta_array(), check_dims()), as the model's own ops
# use greta:::op() (closed_form_states()). The alternative, greta's binomial()
# with the non-integer counts w y of w n, would compute log p and log(1 - p)
# from p itself, losing the precision kept here, and add a binomial
# coefficient of non-integer counts. The class's methods name every function
# with its package (tensorflow::tf), as the parent of the methods'
# environment is the global environment, and a reloaded fit's model has
# none of this file's functions. Checked against plain R
# (weighted_binomial_log_lik()) in R/check_dynamical_model.R. The weighted
# binomial is not a model of the data, so its sample() simulates from the
# noise model the weights stand for, the beta-binomial with mean p and the
# replicate rho, recovered from each weight as rho = (1 / w - 1) / (n - 1)
# (0 for n = 1); greta needs a sample() for calculate() with nsim, as
# fit_model.R uses to cache the posterior means.
weighted_binomial_distribution <- R6::R6Class(
  "weighted_binomial_distribution",
  inherit = greta:::distribution_node,
  public = list(
    initialize = function(size, log_p, log_not_p, weight, dim = NULL) {
      size <- greta:::as.greta_array(size)
      log_p <- greta:::as.greta_array(log_p)
      log_not_p <- greta:::as.greta_array(log_not_p)
      weight <- greta:::as.greta_array(weight)
      dim <- greta:::check_dims(size, log_p, log_not_p, weight,
                                target_dim = dim)
      super$initialize(name = "weighted_binomial", dim = dim, discrete = TRUE)
      self$add_parameter(size, "size")
      self$add_parameter(log_p, "log_p")
      self$add_parameter(log_not_p, "log_not_p")
      self$add_parameter(weight, "weight")
    },
    tf_distrib = function(parameters, dag) {
      size <- parameters$size
      log_p <- parameters$log_p
      log_not_p <- parameters$log_not_p
      weight <- parameters$weight
      # x log p is 0 where x is 0, also where log p is -Inf (0 log 0 = 0), and
      # (size - x) log(1 - p) likewise
      log_prob <- function(x) {
        tf <- tensorflow::tf
        died <- tf$math$multiply_no_nan(log_p, x)
        survived <- tf$math$multiply_no_nan(log_not_p, size - x)
        weight * (died + survived)
      }
      sample <- function(seed) {
        tf <- tensorflow::tf
        tfp <- greta:::tfp
        dtype <- log_p$dtype
        p <- tf$exp(log_p)
        rho <- (1 / weight - 1) / tf$maximum(size - 1, tf$ones_like(size))
        # rho 0 (a binomial) as a tiny rho, for a finite Beta
        rho <- tf$maximum(rho, tf$constant(1e-12, dtype = dtype))
        phi <- 1 / rho - 1
        beta <- tfp$distributions$Beta(concentration1 = p * phi,
                                       concentration0 = (1 - p) * phi)
        q <- beta$sample(seed = seed)
        binomial <- tfp$distributions$Binomial(total_count = size, probs = q)
        binomial$sample(seed = seed)
      }
      list(log_prob = log_prob, sample = sample)
    }
  )
)

# The weighted binomial distribution (weighted_binomial_distribution), as a
# greta array, as greta's own distribution functions return theirs
weighted_binomial <- function(size, log_p, log_not_p, weight, dim = NULL) {
  distribution <- weighted_binomial_distribution$new(size, log_p, log_not_p,
                                                     weight, dim)
  greta:::as.greta_array(distribution$user_node)
}
