# Fit the dynamical model to one cross-validation training fold and return
# posterior predictive draws for the held-out data.
#
# The model is build_dynamical_model() in R/dynamical_model.R, which fit_model.R
# also uses: the likelihood here is restricted to the training fold. Rather than
# collapsing the posterior to a mean predicted fraction and a mean
# overdispersion, this returns the draws themselves, so that the held-out data
# can be scored against the full posterior predictive distribution.
#
# The predicted fraction and the overdispersion are drawn in a single
# calculate() call, so that each pair comes from the same posterior sample.
# Drawing them separately, as the earlier code did, breaks that pairing; it did
# not matter when only their means were used, but a predictive distribution
# needs them coupled.
#
# Because the greta arrays are local to this function they go out of scope when
# it returns, so the model does not need to be purged by hand between folds.

fit_fold <- function(train_df,
                     test_df,
                     before_df = NULL,
                     x_cell_years,
                     cell_years_index,
                     df,
                     classes_index,
                     types,
                     options = dynamical_model_options(),
                     x_cells_init = NULL,
                     settings = dynamical_mcmc_settings(),
                     stored_draws = 2000,
                     inits_file = dynamical_inits_file) {

  # the model, with the likelihood over the training fold (R/dynamical_model.R)
  built <- build_dynamical_model(train_df = train_df,
                                 df = df,
                                 x_cell_years = x_cell_years,
                                 cell_years_index = cell_years_index,
                                 classes_index = classes_index,
                                 types = types,
                                 options = options,
                                 x_cells_init = x_cells_init)

  # use cached posterior means as inits, and the sampler settings
  # (R/dynamical_model.R)
  inits_one <- dynamical_inits(readRDS(inits_file), built$variables,
                               columns = colnames(x_cell_years),
                               country_region_index =
                                 built$lookups$country_region_index,
                               init_covariate_centre =
                                 built$options$init_covariate_centre)
  draws <- run_dynamical_mcmc(built$model, built$variables, inits_one,
                              settings)

  # `n_samples` is taken in one call rather than accumulated in batches towards
  # an effective sample size target. The previous version topped up with
  # extra_samples() until `coda::effectiveSize(draws)` reached 1,000, but that
  # is the effective sample size of the ~690 raw hierarchical parameters, whose
  # minimum was 76-96 on every fold: the target was never reachable and the cap
  # always bound, so this was a fixed-length run with extra bookkeeping. Folds
  # compared with each other must share their sampling settings anyway, so
  # stopping when one fold happens to converge is not an option (#12 review).
  report <- function(...) {
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...), "\n")
    flush(stdout())
  }

  sampled <- settings$n_samples
  ess <- coda::effectiveSize(draws)
  report("sampled %d per chain | raw parameter ESS min %.0f median %.0f",
         sampled, min(ess, na.rm = TRUE), median(ess, na.rm = TRUE))

  # Predictions at the held-out data.
  #
  # These come from greta's calculate() applied to the draws object, which is
  # the supported way to predict from a fitted greta model. Note there is no
  # nsim argument: with nsim, calculate() returns an independent resample of the
  # posterior, which is a valid posterior sample but destroys the MCMC ordering,
  # so effective sample size cannot be recovered from it. Without nsim it
  # returns the draws in order, as an mcmc.list, and the predicted fractions can
  # be diagnosed like any other monitored quantity.
  population_mortality_vec_test <- built$mortality(test_df)
  # the overdispersion of each insecticide type (with rho = "class", its
  # class's)
  rho_types <- built$terms$rho_types

  # Optionally, predictions at a second set of cell-years. The forecasting
  # experiment is scored on the change in mortality between a window before the
  # cut and the holdout window, which differences the site level out and leaves
  # the local slope; that needs the model's prediction in the before window as
  # well as the holdout one. It has to be asked for here, because sampling
  # cannot be resumed in a later session, so a quantity not requested at fitting
  # time cannot be added to a finished fold without refitting (#12 review 5.1).
  report("computing predictions at %d held-out assays%s", nrow(test_df),
         if (is.null(before_df)) "" else
           sprintf(" and %d before-window records", nrow(before_df)))
  prediction_draws <- calculate(
    population_mortality_vec_test = population_mortality_vec_test,
    rho_types = rho_types,
    values = draws
  )

  # effective sample size of the quantities the validation metrics actually
  # consume, rather than of the raw model parameters
  ess_prediction <- coda::effectiveSize(prediction_draws)
  ess_p <- ess_prediction[grep("population_mortality_vec_test\\[",
                               names(ess_prediction))]
  ess_rho <- ess_prediction[grep("rho_types", names(ess_prediction))]

  report("prediction ESS: p median %.0f min %.0f | rho median %.0f min %.0f",
         median(ess_p, na.rm = TRUE), min(ess_p, na.rm = TRUE),
         median(ess_rho, na.rm = TRUE), min(ess_rho, na.rm = TRUE))

  # flatten the mcmc.list to a draws x quantity matrix, preserving order
  prediction_matrix <- as.matrix(prediction_draws)
  p_columns <- grep("population_mortality_vec_test\\[",
                    colnames(prediction_matrix))
  rho_columns <- grep("rho_types", colnames(prediction_matrix))
  p_draws <- prediction_matrix[, p_columns, drop = FALSE]
  rho_type_draws <- prediction_matrix[, rho_columns, drop = FALSE]
  rm(prediction_matrix, prediction_draws)
  invisible(gc())

  # The before-window predictions come from a second calculate() on the same
  # draws object. Splitting them off is purely a memory measure: the five-year
  # folds ask for four times as many predictions as the three-year one, and
  # holding the mcmc.list and its flattened matrix for both windows at once is
  # the largest allocation in the whole fit. It costs nothing in correctness,
  # because calculate(values = draws) with no `nsim` is a deterministic function
  # of the draws — every target is a deterministic function of the sampled
  # parameters — so draw i here is drawn from the same posterior sample as draw
  # i above, and the pairing the change score needs is preserved. (The pairing
  # that must not be broken is the one `nsim` breaks, by resampling.)
  p_draws_before <- NULL
  if (!is.null(before_df) && nrow(before_df) > 0) {
    population_mortality_vec_before <- built$mortality(before_df)
    before_draws <- calculate(
      population_mortality_vec_before = population_mortality_vec_before,
      values = draws
    )
    p_draws_before <- as.matrix(before_draws)
    rm(before_draws)
    invisible(gc())
  }

  # Thin the prediction draws before they are stored, on a single set of indices
  # so that p, rho and the before-window p stay drawn from the same posterior
  # samples. The scoring thins to 2,000 draws anyway, and so does the change
  # score, so nothing downstream sees a difference; what this avoids is holding
  # and writing twenty thousand draws of a five-year holdout, which is four
  # times the prediction volume of the three-year one and would put two
  # concurrent folds into swap. Effective sample size is measured above, on the
  # unthinned ordered draws, so the diagnostics are unaffected.
  keep_draws <- if (nrow(p_draws) > stored_draws) {
    round(seq(1, nrow(p_draws), length.out = stored_draws))
  } else {
    seq_len(nrow(p_draws))
  }
  p_draws <- p_draws[keep_draws, , drop = FALSE]
  rho_type_draws <- rho_type_draws[keep_draws, , drop = FALSE]
  if (!is.null(p_draws_before)) {
    p_draws_before <- p_draws_before[keep_draws, , drop = FALSE]
  }

  convergence <- coda::gelman.diag(draws,
                                   multivariate = FALSE,
                                   autoburnin = FALSE)$psrf
  report("Rhat worst %.3f, %d of %d parameters above 1.01",
         max(convergence[, 1], na.rm = TRUE),
         sum(convergence[, 1] > 1.01, na.rm = TRUE),
         nrow(convergence))

  list(# The draws object, and the greta arrays the predictions were computed
       # from. calculate(values = draws) works on a reloaded draws object,
       # recovering its targets from attr(draws, "model_info"), so prediction
       # survives a session; keeping the arrays explicitly makes that
       # independent of greta's internal layout.
       #
       # Sampling, by contrast, cannot be continued after a session ends:
       # extra_samples() on a reloaded draws object fails with "object is from
       # previous session and is now invalid", because the sampler state is
       # bound to the session that created it, and redefining the model gives
       # new nodes the draws cannot attach to. So a longer run has to be asked
       # for up front through `warmup` and `n_samples` — it cannot be added
       # to a fold afterwards. Note also that folds to be compared with each
       # other must share their sampling settings, so extending one fold means
       # refitting all of them.
       draws = draws,
       prediction_arrays = list(p = population_mortality_vec_test,
                                rho = rho_types),
       p_draws = p_draws,
       # the overdispersion is shared by every assay of an insecticide type, so
       # it is stored by type with the index needed to expand it, rather than
       # as one column per held-out assay. Folds fitted before #20 stored
       # rho_class_draws and class_id instead
       rho_type_draws = rho_type_draws,
       type_id = test_df$type_id,
       options = built$options,
       x_cells_init = x_cells_init,
       test_df = test_df,
       p_draws_before = p_draws_before,
       before_df = before_df,
       convergence = convergence,
       ess = ess,
       ess_p = ess_p,
       ess_rho = ess_rho,
       n_sampled = sampled,
       n_chains = settings$n_chains,
       settings = settings)

}

