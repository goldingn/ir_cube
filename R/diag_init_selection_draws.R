# The parts of the full fit (temporary/fitted_model.RData) that the
# initial-state and selection diagnostics need (R/diag_init_selection.R),
# saved small so the fit is loaded once:
#
#   Rscript R/diag_init_selection_draws.R
#
# Needs no python: every variable was passed to model(), so the draws name
# logit_init_mean, and the parameters are computed in plain R by
# dynamical_parameter_draws() (R/dynamical_predictions.R). Writes
# outputs/diag_init_selection_draws.RDS.

suppressMessages({
  library(tidyverse)
  library(coda)
  library(terra)
})
source("R/functions.R")
source("R/model_covariates.R")
source("R/dynamical_predictions.R")

fit_env <- new.env()
load("temporary/fitted_model.RData", envir = fit_env)
keep <- c("df", "types", "classes", "classes_index", "countries", "regions",
          "baseline_year", "unique_cells", "x_cells_init", "model_options",
          "x_cell_years")
fit <- mget(keep, envir = fit_env)
fold <- list(draws = fit_env$draws, options = fit_env$model_options,
             x_cells_init = fit_env$x_cells_init)
rm(fit_env)
invisible(gc())

options <- fold_options(fold)
draw_index <- paired_draw_index(fold)
parameters <- dynamical_parameter_draws(fold, df = fit$df,
                                        classes_index = fit$classes_index,
                                        types = fit$types,
                                        draw_index = draw_index,
                                        options = options)
parameters$covariate_names <- colnames(fit$x_cell_years)
parameters$draw_chain <- draw_chain(fold$draws)[draw_index]
fit$x_cell_years <- NULL

saveRDS(c(fit, list(options = options, parameters = parameters)),
        "outputs/diag_init_selection_draws.RDS")
cat("saved", parameters$n_draws, "draws\n")
