# The mortality-floor mode of each chain of a fit and its log posterior (#14,
# #37), computed while the model is still in the session: a saved draws object
# keeps the free states but not a model that can evaluate them. Sourced by
# fit_model.R and R/fit_validation_fold.R, after R/dynamical_model.R.

# One row per chain of `draws` (from run_dynamical_mcmc() on greta model
# `model`): the chain's mean mortality floor, its mode ("high" when the mean
# floor is above `high_floor`; the two modes sit near 0.002 and 0.26), and the
# mean and maximum over `n_per_chain` evenly spaced draws of the log posterior
# density, "unadjusted" (on the scale of the parameters, the density the
# modes are compared on) and "adjusted" (with the Jacobian of greta's
# transforms to the free state, the density HMC samples). Draws are evaluated
# `batch` at a time, which bounds the memory of the batched states. NULL for a
# model without a floor. With the species model's two floors (#47), floor is
# NA, other_floor and arabiensis_floor are the chain's mean floors, and mode
# is the mode of each, other members first, e.g. "low/high". With the
# kdr-dependent floor (#47), the floor is that at the mean kdr,
# plogis(floor_intercept) (floor_at_k0), one per class with floor = "class"
# (floor NA, and a mode per class in class order).
chain_floor_modes <- function(model, draws, n_per_chain = 60,
                              high_floor = 0.1, batch = 10) {
  floor_columns <- c(intersect(c("mortality_floor", "other_floor",
                                 "arabiensis_floor"),
                               colnames(draws[[1]])),
                     grep("^floor_intercept", colnames(draws[[1]]),
                          value = TRUE))
  if (length(floor_columns) == 0) {
    return(NULL)
  }
  raw <- attr(draws, "model_info")$raw_draws
  log_prob <- model$dag$generate_log_prob_function(which = "both")
  rows <- lapply(seq_along(draws), function(chain) {
    free <- as.matrix(raw[[chain]])
    index <- unique(round(seq(1, nrow(free), length.out = n_per_chain)))
    values <- lapply(split(index, ceiling(seq_along(index) / batch)),
                     function(i) {
      result <- log_prob(tensorflow::tf$constant(
        free[i, , drop = FALSE], dtype = tensorflow::tf$float64))
      cbind(unadjusted = as.numeric(result$unadjusted),
            adjusted = as.numeric(result$adjusted))
    })
    values <- do.call(rbind, values)
    floors <- colMeans(floor_values(as.matrix(draws[[chain]])[
      , floor_columns, drop = FALSE]))
    names(floors) <- sub("^floor_intercept", "floor_at_k0", names(floors))
    row <- data.frame(chain = chain,
                      floor = if (length(floors) == 1) floors[[1]] else NA,
                      mode = paste(ifelse(floors > high_floor, "high", "low"),
                                   collapse = "/"),
                      log_posterior = mean(values[, "unadjusted"]),
                      log_posterior_max = max(values[, "unadjusted"]),
                      log_posterior_adjusted = mean(values[, "adjusted"]),
                      n_evaluated = nrow(values))
    if (length(floors) > 1) {
      row <- cbind(row, as.list(floors))
    }
    row
  })
  do.call(rbind, rows)
}

# chain_floor_modes(), or NULL with a warning if it fails, so that a finished
# fit is never lost to its report
chain_floor_modes_safely <- function(model, draws, ...) {
  tryCatch({
    modes <- chain_floor_modes(model, draws, ...)
    if (!is.null(modes)) {
      print(modes, digits = 6)
    }
    modes
  }, error = function(e) {
    warning("chain_floor_modes() failed: ", conditionMessage(e))
    NULL
  })
}
