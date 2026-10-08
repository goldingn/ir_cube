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
# model without a floor. With the species model's two floors (#47), floor is
# NA, other_floor and arabiensis_floor are the chain's mean floors, and mode
# is the mode of each, other members first, e.g. "low/high". With the
# kdr-dependent floor (#47), the floor is that at the mean kdr,
# plogis(floor_intercept) (floor_at_k0), one per class with floor = "class"
# (floor NA, and a mode per class in class order); with the latent smooths
# (V5), that where u_f is 0 (floor_at_u0), and the chain's mean sd and range
# (in km; none for a fixed range, which is no variable) of each smooth and
# shear loading, also for a model without a floor.
chain_floor_modes <- function(model, draws, n_per_chain = 60,
                              high_floor = 0.1, batch = 10) {
  floor_columns <- c(intersect(c("mortality_floor", "other_floor",
                                 "arabiensis_floor"),
                               colnames(draws[[1]])),
                     grep("^floor_intercept", colnames(draws[[1]]),
                          value = TRUE))
  smooth_columns <- grep("^smooth_(sd_|inv_range_|shear$)",
                         colnames(draws[[1]]), value = TRUE)
  if (length(floor_columns) == 0 && length(smooth_columns) == 0) {
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
    chain_draws <- as.matrix(draws[[chain]])
    floors <- colMeans(floor_values(chain_draws[, floor_columns,
                                                drop = FALSE]))
    names(floors) <- sub("^floor_intercept",
                         if (length(smooth_columns) > 0) "floor_at_u0" else
                           "floor_at_k0", names(floors))
    row <- data.frame(chain = chain,
                      floor = if (length(floors) == 1) floors[[1]] else NA,
                      mode = if (length(floors) == 0) "none" else
                        paste(ifelse(floors > high_floor, "high", "low"),
                              collapse = "/"),
                      log_posterior = mean(unlist(values)))
    if (length(floors) > 1) {
      row <- cbind(row, as.list(floors))
    }
    for (column in smooth_columns) {
      if (!grepl("^smooth_inv_range_", column)) {
        row[[column]] <- mean(chain_draws[, column])
      } else {
        row[[sub("^smooth_inv_range_", "smooth_range_km_", column)]] <-
          1000 * mean(1 / chain_draws[, column])
      }
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
