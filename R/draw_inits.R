# Initial values from single posterior draws of a full fit, one file per
# chain, in the format fit_model.R caches (dynamical_inits_file), so that
# dynamical_inits() matches them to any model by name: e.g. to start each
# chain of a cross-validation fold at its own draw of the full fit, in the
# posterior's typical set rather than at cached posterior means, and spread
# along the direction the chains mix slowest in (#47).
#
#   Rscript R/draw_inits.R <fit file> <prefix> [variable]
#
# e.g. Rscript R/draw_inits.R \
#        outputs/pod_jobs/v5h_bb_full/temporary/fitted_model.RData \
#        temporary/inits_v5h_draw
# writes temporary/inits_v5h_draw1.RDS, ..., one per chain of the fit (4),
# for
#   IR_CUBE_INITS=temporary/inits_v5h_draw1.RDS,...,temporary/inits_v5h_draw4.RDS
# which starts chain i of a 4-chain run at draw i (dynamical_chain_inits()).
#
# Draw i is from chain i of the fit, and the draws are spread along
# `variable` (default floor_intercept[1,1], the logit floor of the
# pyrethroids where the floor smooth is 0, along which V5h's 2018 forecasting
# fold mixed slowest): the chains are ranked by their mean of it, and the
# chain of rank r gives its draw nearest the (r - 1/2) / n quantile of the
# pooled draws. The draws are those of the fit's own chains (draws_all_chains
# if R/drop_stuck_chains.R dropped some), moved to the variables of the
# non-centred model (noncentred_draws()), the form dynamical_inits() reads.
# Plain R; about 2 GB.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) %in% 2:3)
file <- arguments[1]
prefix <- arguments[2]
variable <- if (length(arguments) == 3) arguments[3] else
  "floor_intercept[1,1]"

suppressMessages({
  library(greta)
  library(dplyr)
  library(stringr)
})
source("R/dynamical_predictions.R")

f <- new.env()
load(file, envir = f)
draws <- if (exists("draws_all_chains", envir = f, inherits = FALSE)) {
  f$draws_all_chains
} else {
  f$draws
}
n_chains <- length(draws)
stopifnot(variable %in% colnames(draws[[1]]))

# the draw of each chain: chains ranked by their mean of `variable`, each
# giving its draw nearest its quantile of the pooled draws
values <- lapply(draws, function(chain) as.matrix(chain)[, variable])
targets <- quantile(unlist(values), (seq_len(n_chains) - 0.5) / n_chains,
                    names = FALSE)
rank <- rank(vapply(values, mean, numeric(1)), ties.method = "first")
iteration <- vapply(seq_len(n_chains), function(chain) {
  which.min(abs(values[[chain]] - targets[rank[chain]]))
}, integer(1))

dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)
for (chain in seq_len(n_chains)) {
  draw_matrix <- as.matrix(draws[[chain]])[iteration[chain], , drop = FALSE]
  names <- unique(sub("\\[.*$", "", colnames(draw_matrix)))
  v <- lapply(setNames(nm = names), extract_parameter,
              draws_matrix = draw_matrix)
  # those of the non-centred model, if the fit centred some levels
  v <- noncentred_draws(v, f$classes_index, f$model_options)
  # with the dimensions of the greta variables (vectors as one-column
  # matrices), as fit_model.R saves them
  one <- lapply(v, function(x) {
    dims <- dim(x)[-1]
    if (length(dims) == 1) dims <- c(dims, 1)
    array(x, dims)
  })
  attr(one, "columns") <- colnames(f$x_cell_years)
  attr(one, "levels") <- dynamical_lookups(f$df)$levels
  out <- sprintf("%s%i.RDS", prefix, chain)
  saveRDS(one, out)
  cat(sprintf("chain %i, iteration %i: %s %.3f (chain mean %.3f, target %.3f) -> %s\n",
              chain, iteration[chain], variable,
              values[[chain]][iteration[chain]], mean(values[[chain]]),
              targets[rank[chain]], out))
}
