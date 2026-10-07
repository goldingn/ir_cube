# The mortality-floor mode of each chain of a fit and its log posterior (#14,
# #37), computed while the model is still in the session: a saved draws object
# keeps the free states but not a model that can evaluate them. Sourced by
# fit_model.R and R/fit_validation_fold.R, after R/dynamical_model.R.

# One row per chain of `draws` (from run_dynamical_mcmc() on greta model
# `model`): the chain's mean mortality floor, its mode ("high" when the mean
# floor is above `high_floor`; the two modes sit near 0.002 and 0.26), and the
# mean over `n_per_chain` evenly spaced draws of the log posterior density on
# the scale of the parameters (without the Jacobian of greta's transforms to
# the free state), the density the modes are compared on. Draws are evaluated
# `batch` at a time, which bounds the memory of the batched states. NULL for a
# model without a floor.
chain_floor_modes <- function(model, draws, n_per_chain = 60,
                              high_floor = 0.1, batch = 10) {
  floor_column <- grep("^mortality_floor", colnames(draws[[1]]), value = TRUE)
  if (length(floor_column) == 0) {
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
      as.numeric(result$unadjusted)
    })
    floor <- mean(as.matrix(draws[[chain]])[, floor_column])
    data.frame(chain = chain, floor = floor,
               mode = if (floor > high_floor) "high" else "low",
               log_posterior = mean(unlist(values)))
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
