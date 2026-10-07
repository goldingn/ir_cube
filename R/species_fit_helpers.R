# Helpers for evaluating saved full fits (temporary/fitted_model.RData of
# R/fit_model.R) of the species runs (#47; R/species_runs.R) and earlier fits:
# loading a fit with its options completed for the current code, its chains
# and their floor modes, and rebuilding its greta model to evaluate its log
# density at free states. Functions only; source after R/dynamical_model.R
# (and R/greta_setup.R's start_greta() before rebuilding a model).

# report() and peak_memory_gb()
source("R/two_stage_helpers.R")

# The floors of a fit: mortality_floor, or the species model's other_floor and
# arabiensis_floor, or the kdr-dependent floor's intercepts (floor_intercept,
# one or one per class: the logit floor at the mean kdr; #47). The scripts
# treat an intercept as its floor, plogis(floor_intercept) (floor_values(),
# R/dynamical_model.R), whose free state is the same, qlogis(floor)
floor_names <- c("mortality_floor", "other_floor", "arabiensis_floor")

# A fit's saved options, completed for the current code. Options added since
# the fit take the values that reproduce it: floor_prior (#47) the prior of
# the floor before #47, Beta(1, 49), unless `floor_prior` gives another (the
# fits before d2dee17, 2 October 2026, used Beta(1, 9)); species and kdr off;
# with kdr, no kdr-dependent floor.
# Settings of older code that the current code no longer has are dropped if
# their values are what the current code does, and stop the script otherwise.
complete_model_options <- function(options, floor_prior = c(1, 49)) {
  defaults <- dynamical_model_options()
  equivalent <- list(rho = "type", init_centred = "country")
  for (name in setdiff(names(options),
                       c(names(defaults), "init_covariate_centre"))) {
    if (!identical(options[[name]], equivalent[[name]])) {
      stop("the fit's option ", name, " = ", deparse(options[[name]]),
           " is not in the current model")
    }
  }
  design <- options$selection_columns
  design_equivalent <- list(net_transform = "linear", net_scale = NULL,
                            trend_after = "continue")
  for (name in setdiff(names(design), names(formals(selection_design)))) {
    if (!identical(design[[name]], design_equivalent[[name]])) {
      stop("the fit's selection setting ", name, " = ",
           deparse(design[[name]]), " is not in the current model")
    }
  }
  out <- options[intersect(names(options),
                           c(names(defaults), "init_covariate_centre"))]
  if (is.null(out$floor_prior)) out$floor_prior <- floor_prior
  if (is.null(out$species)) out$species <- FALSE
  if (is.null(out$kdr)) out$kdr <- FALSE
  # kdr options saved before the kdr-dependent floor have none
  if (is.list(out$kdr) && is.null(out$kdr$floor)) out$kdr$floor <- FALSE
  out$selection_columns <- complete_selection_design(design)
  out
}

# The floor prior for fits saved without one, from FLOOR_PRIOR (e.g. "1,9"),
# else Beta(1, 49)
floor_prior_override <- function() {
  value <- Sys.getenv("FLOOR_PRIOR")
  if (!nzchar(value)) return(c(1, 49))
  as.numeric(strsplit(value, ",")[[1]])
}

# A saved fit, as a list: every chain (draws_all_chains if
# R/drop_stuck_chains.R has dropped some, else draws), the chains recorded as
# stuck, chain_modes (R/chain_floor_modes.R; NULL for fits before it), the
# completed options, and the data, covariates and lookups build_dynamical_model()
# takes. Run with the greta library loaded (the draws are greta objects).
load_fit <- function(file, floor_prior = floor_prior_override()) {
  e <- new.env()
  load(file, envir = e)
  has <- function(name) exists(name, envir = e, inherits = FALSE)
  list(file = file,
       draws = if (has("draws_all_chains")) e$draws_all_chains else e$draws,
       recorded_stuck = if (has("stuck_chains")) e$stuck_chains else integer(0),
       chain_modes = if (has("chain_modes")) e$chain_modes,
       saved_options = e$model_options,
       options = complete_model_options(e$model_options, floor_prior),
       df = e$df,
       x_cell_years = e$x_cell_years,
       cell_years_index = e$cell_years_index,
       x_cells_init = e$x_cells_init,
       classes_index = e$classes_index,
       types = e$types,
       classes = e$classes,
       countries = e$countries,
       regions = e$regions,
       baseline_year = e$baseline_year)
}

# The floors of `draws` (a greta_mcmc_list) it has
fit_floor_names <- function(draws) {
  c(intersect(floor_names, colnames(draws[[1]])),
    grep("^floor_intercept", colnames(draws[[1]]), value = TRUE))
}

# The floor mode of each chain of `draws`: as chain_floor_modes(), "low" or
# "high" for each floor by its chain mean against `high_floor`, other members
# first with the species model ("low/high"); "none" for a fit without a floor
chain_floor_mode <- function(draws, high_floor = 0.1) {
  floors <- fit_floor_names(draws)
  vapply(draws, function(chain) {
    if (length(floors) == 0) return("none")
    means <- colMeans(floor_values(as.matrix(chain)[, floors, drop = FALSE]))
    paste(ifelse(means > high_floor, "high", "low"), collapse = "/")
  }, character(1))
}

# The chains of `draws` to use: all but those stuck (stuck_chains(),
# R/dynamical_predictions.R), reported
usable_chains <- function(draws, label = "") {
  stuck <- stuck_chains(draws)
  if (length(stuck) > 0) {
    report("%s: chain(s) %s stuck (fewer than half their draws distinct), left out",
           label, toString(stuck))
  }
  setdiff(seq_along(draws), stuck)
}

# The fit's greta model, rebuilt with its options, data and covariates
# (build_dynamical_model()), with `free_order`, the column of the fit's free
# states (attr(draws, "model_info")$raw_draws) for each column of the rebuilt
# model's: greta orders the variables in the free state as it reaches them in
# its graph, which changed with the code (r2_main, 3 October 2026, has another
# order), so they are matched by name through the fit's own model
# (model_info$model). Checked against the draws: the rebuilt model's values
# at the fit's free states must be the saved draws, for `n_check` draws of
# each chain.
rebuild_fit_model <- function(fit, n_check = 5) {
  built <- build_dynamical_model(train_df = fit$df,
                                 df = fit$df,
                                 x_cell_years = fit$x_cell_years,
                                 cell_years_index = fit$cell_years_index,
                                 classes_index = fit$classes_index,
                                 types = fit$types,
                                 options = fit$options,
                                 x_cells_init = fit$x_cells_init)
  new <- free_state_columns(built$model)
  old <- free_state_columns(attr(fit$draws, "model_info")$model)
  stopifnot(setequal(names(attr(new, "targets")), names(attr(old, "targets"))))
  free_order <- integer(sum(lengths(new)))
  for (name in names(attr(new, "targets"))) {
    free_order[new[[attr(new, "targets")[[name]]]]] <-
      old[[attr(old, "targets")[[name]]]]
  }
  built$free_order <- free_order
  raw <- attr(fit$draws, "model_info")$raw_draws
  difference <- 0
  for (chain in seq_along(raw)) {
    free <- as.matrix(raw[[chain]])[, free_order, drop = FALSE]
    index <- unique(round(seq(1, nrow(free), length.out = n_check)))
    trace <- built$model$dag$trace_values(free[index, , drop = FALSE])
    saved <- as.matrix(fit$draws[[chain]])[index, , drop = FALSE]
    stopifnot(identical(colnames(trace), colnames(saved)))
    difference <- max(difference, abs(trace - saved))
  }
  report("rebuilt model: %d free parameters (%s order); values at the fit's free states vs its draws: max abs diff %.3g",
         length(free_order),
         if (identical(free_order, seq_along(free_order))) "the same" else
           "another", difference)
  stopifnot(difference < 1e-8)
  built
}

# About n_draws draws of `fit`, the same number evenly spaced in each usable
# chain (loo::relative_eff() needs equal chains), as rows of
# as.matrix(fit$draws) (the chains stacked in order): list(index, chain)
even_draws <- function(fit, n_draws, label = "") {
  usable <- usable_chains(fit$draws, label)
  chain_rows <- split(seq_len(sum(vapply(fit$draws, nrow, integer(1)))),
                      draw_chain(fit$draws))
  per_chain <- min(floor(n_draws / length(usable)),
                   min(lengths(chain_rows[usable])))
  index <- unlist(lapply(chain_rows[usable], function(rows) {
    rows[round(seq(1, length(rows), length.out = per_chain))]
  }), use.names = FALSE)
  list(index = index, chain = draw_chain(fit$draws)[index])
}

# The posterior parameters of `fit` (dynamical_parameter_draws(),
# R/dynamical_predictions.R) at draws `index`
fit_parameter_draws <- function(fit, index) {
  dynamical_parameter_draws(
    list(draws = fit$draws, options = fit$options,
         x_cells_init = fit$x_cells_init),
    fit$classes_index, fit$types, draw_index = index, options = fit$options)
}

# The regions of the misfit and trend summaries (#47): the UN geoscheme's,
# with the Horn of Africa (and Sudan) split from Eastern Africa and the
# southern countries of Eastern Africa moved to Southern Africa, as
# kdr_region() on branch latent-kdr. `region` is the UNSD region of each
# record (df$region)
analysis_region <- function(country, region) {
  dplyr::case_when(
    country %in% c("Djibouti", "Eritrea", "Ethiopia", "Somalia",
                   "Sudan") ~ "Horn",
    country %in% c("Comoros", "Madagascar", "Malawi", "Mauritius",
                   "Mozambique", "Zambia", "Zimbabwe") ~ "Southern",
    region == "Western Africa" ~ "West",
    region == "Middle Africa" ~ "Central",
    region == "Eastern Africa" ~ "East",
    region == "Southern Africa" ~ "Southern",
    region == "Northern Africa" ~ "North")
}

# The free states of `fit`'s draws, chains `chains` stacked, in the order of
# the rebuilt model `built` (rebuild_fit_model())
free_states <- function(fit, built, chains = seq_along(fit$draws)) {
  raw <- attr(fit$draws, "model_info")$raw_draws
  do.call(rbind, lapply(chains, function(chain) {
    as.matrix(raw[[chain]])[, built$free_order, drop = FALSE]
  }))
}

# The free-state column of variable `name` of greta model `model`
# (free_state_columns(), R/dynamical_predictions.R)
free_column <- function(model, name) {
  columns <- free_state_columns(model)
  base <- sub("\\[.*$", "", name)
  block <- columns[[attr(columns, "targets")[[base]]]]
  if (base == name) {
    return(block)
  }
  # an element of a vector variable, e.g. floor_intercept[2,1], whose free
  # state is in element order
  index <- as.integer(strsplit(sub("^.*\\[(.*)\\]$", "\\1", name), ",")[[1]])
  stopifnot(length(index) == 2, index[2] == 1)
  block[index[1]]
}

# The log density of greta model `model` at free states `free` (states x
# free parameters), batched: adjusted (with the Jacobians of greta's
# transforms to the free state, the density HMC samples) and unadjusted (on
# the parameters' own scales), as a states x 2 matrix
# (compiled as a TensorFlow function, about 15 times faster than eager; each
# batch is padded to `batch` states, so it is traced once)
log_density_function <- function(model, batch = 25) {
  tf <- tensorflow::tf
  f <- model$dag$generate_log_prob_function(which = "both")
  compiled <- tf$`function`(function(x) {
    result <- f(x)
    list(result$adjusted, result$unadjusted)
  })
  n_free <- length(unlist(model$dag$example_parameters(free = TRUE)))
  function(free) {
    free <- matrix(free, ncol = n_free)
    out <- lapply(split(seq_len(nrow(free)),
                        ceiling(seq_len(nrow(free)) / batch)), function(i) {
      x <- free[c(i, rep(i[1], batch - length(i))), , drop = FALSE]
      result <- compiled(tf$constant(x, dtype = tf$float64))
      cbind(adjusted = as.numeric(as.array(result[[1]])),
            unadjusted = as.numeric(as.array(result[[2]])))[
              seq_along(i), , drop = FALSE]
    })
    do.call(rbind, out)
  }
}

# The adjusted log density of greta model `model` and its gradient at free
# states `free` (states x free parameters), as list(value, gradient): one
# TensorFlow call per batch; the gradient of each state's density is that of
# the batch's sum
log_density_gradient_function <- function(model) {
  tf <- tensorflow::tf
  `%as%` <- reticulate::`%as%`
  f <- model$dag$generate_log_prob_function(which = "adjusted")
  value_and_gradient <- tf$`function`(function(x) {
    with(tf$GradientTape() %as% tape, {
      tape$watch(x)
      value <- f(x)
    })
    list(value, tape$gradient(value, x))
  })
  function(free) {
    free <- if (is.matrix(free)) free else matrix(free, nrow = 1)
    result <- value_and_gradient(tf$constant(free, dtype = tf$float64))
    list(value = as.numeric(as.array(result[[1]])),
         gradient = matrix(as.array(result[[2]]), nrow(free)))
  }
}

# The floor's free state for floor value f: greta's transform of a beta
# variable to (0, 1) is the logistic, checked by rebuild checks of the floor
# (floor_free_check())
floor_free <- function(f) qlogis(f)

# Checks that setting the floor's free state to floor_free(f) gives floor f
floor_free_check <- function(model, free, name, f = 0.123) {
  free[free_column(model, name)] <- floor_free(f)
  trace <- model$dag$trace_values(matrix(free, nrow = 1))
  stopifnot(abs(floor_values(trace[, name, drop = FALSE])[1, 1] - f) < 1e-12)
  invisible(TRUE)
}
