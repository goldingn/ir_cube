# Make the initial values for the fits (dynamical_inits_file, in
# R/dynamical_model.R) from a full fit: the posterior mean of every variable,
# over the chains R/drop_stuck_chains.R kept, with the labels of their
# dimensions as dimnames, by which dynamical_inits() matches them to a model.
#
#   Rscript R/make_dynamical_inits.R [fitted_model.RData]
#
# Run R/drop_stuck_chains.R on the fit first. Then commit the file, and name
# the fit in the comment above dynamical_inits_file. Run with the greta
# library (R/packages.R), which the saved objects need.

suppressMessages({
  library(greta)
  library(dplyr)
  library(stringr)
})
# extract_parameter(), and dynamical_inits_file
source("R/dynamical_predictions.R")

arguments <- commandArgs(trailingOnly = TRUE)
fit_file <- if (length(arguments) >= 1) arguments[1] else
  "temporary/fitted_model.RData"

fit <- new.env()
load(fit_file, envir = fit)
if (length(stuck_chains(fit$draws)) > 0) {
  stop(fit_file, " has stuck chains: run R/drop_stuck_chains.R on it first")
}

# every variable is a target of model(), so its draws are columns of `draws`
draws_matrix <- as.matrix(fit$draws)
labels <- fit$built$variable_labels
post_means <- lapply(setNames(nm = names(fit$built$variables)), function(name) {
  x <- extract_parameter(draws_matrix, name)
  means <- colMeans(matrix(x, nrow(x)))
  array(means, lengths(labels[[name]]), labels[[name]])
})

saveRDS(post_means, dynamical_inits_file)
cat("wrote", dynamical_inits_file, "from", fit_file, "(",
    length(fit$draws), "chains,", nrow(draws_matrix), "draws )\n")
