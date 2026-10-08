# Recompute the dynamical model's predictions, as logit draws, at arbitrary
# cell-years from its saved posterior draws, in plain R.
#
# A saved draws object cannot be resumed to ask greta for more predictions, so
# the two-stage model (paired draw for draw with a fold's stored test
# predictions), R/predict.R and the figure scripts recompute them from the
# sampled parameters, with the model's own transforms (dynamical_terms() in
# R/dynamical_model.R) and its options (fold$options). The recursion is
# solved in closed form on the logit scale, as in closed_form_states():
#
#   logit q_t = logit q_0 - sum_{s <= t} log w_s - t kappa
#
# the state for year index t having had the fitness of years 1..t applied.
# The logit draws are of predicted mortality: q_t, or with a mortality floor
# f, f + (1 - f) q_t. With the species model (#47), the one cumulative log
# fitness gives the states of the other members of the complex and of
# arabiensis (outer_mortality() in R/dynamical_model.R), and the
# predictions are of their mixture at a share of arabiensis: each bioassay's
# own (arabiensis_share(), R/species.R) in dynamical_logit(), and one given
# per cell in dynamical_logit_cells() (r(x) for the whole complex, 0 or 1 for
# one species). With the kdr covariate (#47), each trajectory's cumulative log
# fitness and reversion are multiplied by its factors at the cell's kdr, by
# the same outer form, and with the latent smooths (V5), the cumulative log
# fitness by exp(u_s(x)) and the floor shifted by u_f(x) at the cell, from the
# basis at the cell (prediction_basis(), R/latent_smooth.R).

source("R/dynamical_model.R")
# thin_draws() and max_draws
source("R/validation_scoring.R")

# Draw indices into as.matrix(fold$draws) that the stored predictions
# correspond to: newer folds thin p_draws to `maximum` at fitting time with
# this rule, older folds are thinned by the same rule at scoring time. The
# rule is thin_draws() (R/validation_scoring.R). Chains recorded as stuck
# (fold$stuck_chains, set by R/drop_stuck_chains.R, which drops their rows
# from the stored predictions) are left out.
paired_draw_index <- function(fold, maximum = max_draws) {
  n_total <- sum(vapply(fold$draws, nrow, integer(1)))
  index <- thin_draws(matrix(seq_len(n_total)), maximum)[, 1]
  if (length(fold$stuck_chains) > 0) {
    index <- index[!draw_chain(fold$draws)[index] %in% fold$stuck_chains]
  }
  index
}

# The chain of each row of as.matrix(draws).
draw_chain <- function(draws) {
  rep(seq_along(draws), vapply(draws, nrow, integer(1)))
}

# The chains of `draws` that stopped moving: those whose share of distinct
# draws is below `min_distinct`. With one step size for all chains
# (windowed_hmc()), a chain can sit where every trajectory at that step size
# fails, and then repeats one draw: in the October 2026 run, one chain of
# spatial blocks fold 1 kept a single draw for all 3,000 samples.
stuck_chains <- function(draws, min_distinct = 0.5) {
  which(vapply(draws, function(chain) {
    chain <- as.matrix(chain)
    nrow(unique(chain)) / nrow(chain) < min_distinct
  }, logical(1)))
}

# `draws` (a greta_mcmc_list) with only the chains `keep`, in the raw
# free-state draws that calculate() reads as well.
drop_chains <- function(draws, keep) {
  out <- draws[keep]
  model_info <- attr(draws, "model_info")
  model_info$raw_draws <- model_info$raw_draws[keep]
  attr(out, "model_info") <- model_info
  class(out) <- class(draws)
  out
}

# One named parameter of a draws x parameter matrix, as a draws x
# dim(parameter) array. greta names elements by their R (column-major) index,
# e.g. "beta_type_raw[13,9]", so the index is parsed from the names
extract_parameter <- function(draws_matrix, name) {
  # a scalar has one column, named without an index
  if (name %in% colnames(draws_matrix)) {
    return(draws_matrix[, name, drop = FALSE])
  }
  columns <- grep(paste0("^", name, "\\["), colnames(draws_matrix))
  stopifnot(length(columns) > 0)
  index <- colnames(draws_matrix)[columns] %>%
    str_remove(paste0("^", name, "\\[")) %>%
    str_remove("\\]$") %>%
    str_split(",", simplify = TRUE)
  index <- matrix(as.integer(index), nrow = length(columns))
  dims <- apply(index, 2, max)
  # vectors are stored as column vectors, [i,1]
  if (length(dims) == 2 && dims[2] == 1) dims <- dims[1]
  linear <- if (length(dims) == 1) index[, 1] else
    index[, 1] + (index[, 2] - 1) * dims[1]
  out <- array(NA_real_, dim = c(nrow(draws_matrix), prod(dims)))
  out[, linear] <- draws_matrix[, columns]
  dim(out) <- c(nrow(draws_matrix), dims)
  out
}

# The columns of each variable's free state in the raw draws of greta model
# `model`, which concatenate them in the dag's variable order, named by the
# variables' TensorFlow names; attribute "targets" maps the model's target
# names to those names.
free_state_columns <- function(model) {
  dag <- model$dag
  sizes <- vapply(dag$example_parameters(free = TRUE), length, integer(1))
  columns <- Map(function(end, size) (end - size + 1):end, cumsum(sizes),
                 sizes)
  node_names <- vapply(dag$node_list, function(node) node$unique_name,
                       character(1))
  targets <- vapply(model$target_greta_arrays,
                    function(x) greta:::get_node(x)$unique_name,
                    character(1))
  attr(columns, "targets") <- setNames(
    dag$get_tf_names()[match(targets, node_names)], names(targets))
  columns
}

# dynamical_terms() for every draw of the variables `v` (a named list of
# draws x dim arrays), stacked as draws x dim(term) arrays, for the terms
# `terms`.
dynamical_terms_draws <- function(v, classes_index, types, terms, options) {
  n_draws <- nrow(v[[1]])
  out <- NULL
  for (i in seq_len(n_draws)) {
    v_i <- variables_at_draw(v, i)
    terms_i <- dynamical_terms(v_i, classes_index = classes_index,
                               types = types, options = options)[terms]
    if (is.null(out)) {
      out <- lapply(terms_i, function(x) {
        array(NA_real_, c(n_draws, dim(as.matrix(x))))
      })
    }
    for (name in terms) {
      out[[name]][i, , ] <- terms_i[[name]]
    }
  }
  out
}

# The parameters the predictions need, for the draws in `draw_index`:
#   effect_type          draws x n_covs x n_types, exp(beta_type)
#   logit_init_relative  draws x n_countries x n_types, the logit relative
#                        initial state (above init_frac_min) of each fitted
#                        country, to which a cell's initial-state covariate
#                        effects are added (dynamical_logit_cells())
#   init_coef            draws x n_init_covs x n_types, the coefficients of
#                        the initial-state covariates (NULL for none)
#   rho_types            draws x n_types, the observation overdispersion;
#                        with the weighted binomial likelihood (#47), the
#                        fixed replicate rho in every draw (fixed_rho_types())
#   mortality_floor      draws (NULL for none)
#   kappa_type           draws x n_types, the reversion kappa (<= 0; NULL for
#                        none, see reversion_kappa())
#   gamma_selection, gamma_cost, other_floor, arabiensis_floor
#                        draws, with the species model (#47; NULL without it,
#                        and gamma_cost and the floors NULL when not in the
#                        model; see dynamical_terms())
#   kdr_slopes           the slopes of the kdr covariate (#47; kdr_slope_names),
#                        a list of draws each, empty without it
#   floor_intercept      draws x 1 or classes, and
#   floor_kdr            draws, the kdr-dependent floor (#47; NULL without
#                        it; kdr_floor_value()); floor_intercept is that of
#                        the floor of the latent smooths too
#   smooth_weights       the weights of the latent smooths' basis functions
#                        (V5; smooth_weight_terms()), a list of draws x basis
#                        functions, named by smooth, empty without them
#   init_min             n_types, init_frac_min
#   x_cells_init         the fit's initial-state covariates, one row per
#                        cell_id (NULL for none)
#   variables            every variable, as draws x dim arrays, from which
#                        map_logit_init() draws countries without data
# and the options, types and classes_index.
dynamical_parameter_draws <- function(fold,
                                      classes_index,
                                      types,
                                      draw_index = paired_draw_index(fold),
                                      options = fold$options) {

  # as.matrix() on an mcmc.list stacks chains in order, which is how fit_fold()
  # flattened the calculate() output that p_draws came from
  draws_matrix <- as.matrix(fold$draws)[draw_index, , drop = FALSE]
  n_draws <- nrow(draws_matrix)
  n_types <- length(types)

  variable_names <- unique(sub("\\[.*$", "", colnames(draws_matrix)))
  v <- lapply(setNames(nm = variable_names), extract_parameter,
              draws_matrix = draws_matrix)
  stopifnot(!is.null(options),
            identical(dim(v$logit_init_mean), c(n_draws, n_types)),
            species_on(options) == !is.null(v$gamma_selection),
            kdr_on(options) ==
              any(kdr_slope_names %in% names(v)),
            smooth_on(options) == any(grepl("^smooth_raw_", names(v))))

  reversion <- !isFALSE(options$reversion)
  terms <- dynamical_terms_draws(
    v, classes_index, types,
    terms = c("beta_type", "logit_init_relative", "rho_types",
              if (reversion) "kappa_type"),
    options = options)

  list(effect_type = exp(terms$beta_type),
       logit_init_relative = terms$logit_init_relative,
       init_coef = if (!is.null(options$init_covariates)) {
         array(v$init_coef,
               c(n_draws, length(options$init_covariates), n_types),
               dimnames = list(NULL, options$init_covariates, types))
       },
       rho_types = matrix(terms$rho_types, n_draws),
       mortality_floor = if (isTRUE(options$mortality_floor)) {
         c(v$mortality_floor)
       },
       kappa_type = if (reversion) matrix(terms$kappa_type, n_draws),
       gamma_selection = if (!is.null(v$gamma_selection)) {
         c(v$gamma_selection)
       },
       gamma_cost = if (!is.null(v$gamma_cost)) c(v$gamma_cost),
       other_floor = if (!is.null(v$other_floor)) c(v$other_floor),
       arabiensis_floor = if (!is.null(v$arabiensis_floor)) {
         c(v$arabiensis_floor)
       },
       kdr_slopes = lapply(v[intersect(kdr_slope_names, names(v))], c),
       floor_intercept = if (!is.null(v$floor_intercept)) {
         matrix(v$floor_intercept, n_draws)
       },
       floor_kdr = if (!is.null(v$floor_kdr)) c(v$floor_kdr),
       smooth_weights = smooth_weight_terms(v, options),
       init_min = init_frac_constants(types)$min,
       x_cells_init = select_init_covariates(fold$x_cells_init, options),
       variables = v,
       options = options,
       types = types,
       classes_index = classes_index,
       n_draws = n_draws)
}

# `parameters` (dynamical_parameter_draws()) for its draws `draws` only
subset_draws <- function(parameters, draws) {
  rows <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.null(dim(x))) return(x[draws])
    do.call(`[`, c(list(x, draws), rep(list(TRUE), length(dim(x)) - 1),
                   drop = FALSE))
  }
  for (name in c("effect_type", "logit_init_relative", "init_coef",
                 "rho_types", "mortality_floor", "kappa_type",
                 "gamma_selection", "gamma_cost", "other_floor",
                 "arabiensis_floor", "floor_intercept", "floor_kdr")) {
    parameters[name] <- list(rows(parameters[[name]]))
  }
  parameters$kdr_slopes <- lapply(parameters$kdr_slopes, rows)
  parameters$smooth_weights <- lapply(parameters$smooth_weights, rows)
  parameters$variables <- lapply(parameters$variables, rows)
  parameters$n_draws <- length(draws)
  parameters
}

# Logit q_0 at cells, draws x cells, for type k, from the logit relative
# initial state at each cell's country (`logit_init`, draws x cells) and the
# cells' initial-state covariates x_init (cells x covariates, named columns;
# for fits with them), as logit_init_relative_rows() and closed_form_states()
# in the model
cell_logit_init <- function(parameters, k, logit_init, x_init = NULL) {
  init_coef <- parameters$init_coef
  if (!is.null(init_coef)) {
    stopifnot(!is.null(x_init), nrow(x_init) == ncol(logit_init))
    logit_init <- logit_init +
      matrix(init_coef[, , k], nrow = parameters$n_draws) %*%
      t(x_init[, dimnames(init_coef)[[2]], drop = FALSE])
  }
  floored_logit(logit_init, parameters$init_min[k])
}

# The recursion for cells of insecticide type k:
#   parameters  dynamical_parameter_draws() (or subset_draws())
#   logit_init  draws x cells, the logit relative initial state at each cell's
#               country (parameters$logit_init_relative, or map_logit_init())
#   x           cells x years x n_covs covariates, year index 1 the baseline
#               year, covariates in the column order of x_cell_years
#   years_keep  year indices to return
#   x_init      cells x initial-state covariates (named columns), for fits
#               with them
#   share       with the species model (#47) only, and needed there: the
#               arabiensis share of the predictions at each cell (length 1
#               or cells): r(x) for the whole complex (prediction_share()),
#               0 for the other members or 1 for arabiensis
#   kdr         with the kdr covariate (#47) only, and needed there: the
#               standardised kdr at each cell, cells x kdr_bands()
#               (prediction_kdr())
#   basis       with the latent smooths (V5) only, and needed there: their
#               centred basis at each cell, cells x basis functions
#               (prediction_basis())
# Returns a list named by years_keep of draws x cells logit mortality.
dynamical_logit_cells <- function(parameters, k, logit_init, x, years_keep,
                                  x_init = NULL, share = NULL, kdr = NULL,
                                  basis = NULL) {
  mix_trajectories(dynamical_trajectories(parameters, k, logit_init, x,
                                          years_keep, x_init, kdr, basis),
                   share)
}

# The logit mortality, by year, of the output of dynamical_trajectories() at
# the arabiensis share `share` (as dynamical_logit_cells() takes it): the
# trajectories themselves without the species model
mix_trajectories <- function(trajectories, share) {
  if (!is.list(trajectories[[1]])) {
    stopifnot(is.null(share))
    return(trajectories)
  }
  if (is.null(share)) {
    stop("the species model (#47) predicts at an arabiensis share: give ",
         "dynamical_logit_cells() a share, e.g. prediction_share()")
  }
  lapply(trajectories, function(year) {
    mixture_logit(year$arabiensis, year$other, share)
  })
}

# The recursion of dynamical_logit_cells(), returning a list named by
# years_keep of draws x cells logit mortality, or with the species model, of
# lists of two of them, "other" (the other members of the complex) and
# "arabiensis", whose log fitness from selection is the other members' times
# exp(gamma_selection) and reversion kappa theirs times exp(gamma_cost); each
# species has its own mortality floor, if any. With the species model, the
# kdr covariate or the latent smooths, this is the outer form of
# outer_mortality() (R/dynamical_model.R): one cumulative log fitness serves
# every trajectory.
dynamical_trajectories <- function(parameters, k, logit_init, x, years_keep,
                                   x_init = NULL, kdr = NULL, basis = NULL) {
  n_cells <- dim(x)[1]
  n_draws <- parameters$n_draws
  stopifnot(ncol(logit_init) == n_cells, nrow(logit_init) == n_draws,
            max(years_keep) <= dim(x)[2])

  logit_init <- cell_logit_init(parameters, k, logit_init, x_init)

  effect <- matrix(parameters$effect_type[, , k], nrow = n_draws)
  kappa <- parameters$kappa_type[, k]
  species <- species_on(parameters$options)
  outer <- species || kdr_on(parameters$options) ||
    smooth_on(parameters$options)
  if (smooth_on(parameters$options)) {
    if (is.null(basis)) {
      stop("the latent smooths (V5) need their basis at the cells: give ",
           "dynamical_logit_cells() basis, e.g. prediction_basis()")
    }
    stopifnot(nrow(basis) == n_cells)
  } else {
    stopifnot(is.null(basis))
  }
  if (kdr_on(parameters$options)) {
    if (is.null(kdr)) {
      stop("the kdr covariate (#47) needs the standardised kdr at the ",
           "cells: give dynamical_logit_cells() kdr, e.g. prediction_kdr()")
    }
    stopifnot(nrow(kdr) == n_cells)
  } else {
    stopifnot(is.null(kdr))
  }
  # the floor of one trajectory: the constant mortality_floor, or the
  # kdr-dependent floor at the cells, draws x cells (kdr_floor_value())
  floor <- parameters$mortality_floor
  if (!isFALSE(kdr_floor(parameters$options))) {
    class <- 1L
    k_floor <- kdr[, "complex"]
    if (identical(kdr_floor(parameters$options), "class")) {
      class <- parameters$classes_index[k]
      k_floor <- k_floor * parameters$options$kdr$floor_classes[class]
    }
    floor <- kdr_floor_value(parameters$floor_intercept[, class],
                             parameters$floor_kdr, k_floor)
  }
  # the latent smooths at the cells for type k, draws x cells, each NULL
  # where the model or k's class has none; and the floor of the smooths, at
  # k's class's intercept
  smooth <- list()
  if (smooth_on(parameters$options)) {
    class <- parameters$classes_index[k]
    for (kind in names(parameters$smooth_weights)) {
      if (smooth_class_weight(parameters$options, kind, class) == 1) {
        smooth[[kind]] <- parameters$smooth_weights[[kind]] %*% t(basis)
      }
    }
    if (smooth_floor_on(parameters$options)) {
      floor <- smooth_floor_value(
        parameters$floor_intercept[
          , smooth_intercept_index(parameters$options, class)],
        smooth$floor)
    }
  }
  # logit mortality of one trajectory of the outer form, from `rows`
  # (outer_logit()), with species offsets gamma (NULL for none), the kdr of
  # `band` and the slopes named with `suffix` (kdr_slope_names), and `floor`
  slopes <- parameters$kdr_slopes
  trajectory <- function(rows, gamma_selection, gamma_cost, band, suffix,
                         floor) {
    k_band <- if (!is.null(kdr)) kdr[, band]
    floored_logit(
      outer_logit(rows,
                  log_selection = outer_log_multiplier(
                    gamma_selection,
                    slopes[[paste0("delta_selection", suffix)]], k_band,
                    smooth$selection),
                  log_cost = outer_log_multiplier(
                    gamma_cost, slopes[[paste0("delta_cost", suffix)]],
                    k_band)),
      floor)
  }
  cumulative <- 0
  out <- list()
  for (t in seq_len(max(years_keep))) {
    # log1p because the selection term can be tiny
    log_w <- log1p(effect %*% t(matrix(x[, t, ], nrow = n_cells)))
    cumulative <- cumulative + log_w
    if (outer) {
      # cumulative is the log fitness alone; the multipliers are one per draw
      # (gamma), or per draw and cell (with kdr or the latent smooths)
      if (t %in% years_keep) {
        rows <- list(logit_init = logit_init, cumulative = cumulative,
                     reversion = if (!is.null(kappa)) t * kappa)
        out[[as.character(t)]] <- if (!species) {
          trajectory(rows, NULL, NULL, "complex", "", floor)
        } else {
          list(other = trajectory(rows, NULL, NULL, "other", "_other",
                                  parameters$other_floor),
               arabiensis = trajectory(rows, parameters$gamma_selection,
                                       parameters$gamma_cost, "arabiensis",
                                       "_arabiensis",
                                       parameters$arabiensis_floor))
        }
      }
      next
    }
    # reversion: - t kappa in year t (reversion_kappa()), one per draw
    if (!is.null(kappa)) cumulative <- cumulative + kappa
    if (t %in% years_keep) {
      out[[as.character(t)]] <- floored_logit(
        logit_init - cumulative, parameters$mortality_floor)
    }
  }
  out
}

# Draws of the dynamical model's logit predicted mortality at (cell, type,
# year) rows (cell_id, type_id, year_id, as in `df`), with the covariates of
# x_cell_years (one row per (cell_id, year_id) of cell_years_index). The
# initial condition of a cell is that of the country of its first record in
# the full `df`, as in the model (dynamical_lookups()), whatever the rows'
# country_id. Returns a draws x nrow(rows) matrix, paired with
# thin_draws(fold$p_draws) when `parameters` are at paired_draw_index().
# With the species model (#47), the prediction at each row is the mixture at
# its arabiensis share, `share` (one per row, or one for all), by default the
# bioassay's own (arabiensis_share(): the rows then need species and cell).
# With the kdr covariate, each cell's kdr is that of its mask cell in `df`,
# and with the latent smooths, so is each cell's basis.
dynamical_logit <- function(parameters, rows, df, x_cell_years,
                            cell_years_index, max_block = 2.5e7,
                            share = NULL) {

  n_draws <- parameters$n_draws
  n_times <- max(cell_years_index$year_id)
  if (any(rows$year_id > n_times | rows$year_id < 1)) {
    stop("rows fall outside the covariate years 1..", n_times)
  }
  cell_country <- dynamical_lookups(df)$cell_country_lookup
  stopifnot(!anyNA(cell_country[rows$cell_id]))
  # the standardised kdr at each cell_id, NULL without the kdr covariate, and
  # the centred basis of the latent smooths, NULL without them
  mask_cells <- df$cell[match(seq_len(max(df$cell_id)), df$cell_id)]
  kdr_cells <- prediction_kdr(parameters$options, mask_cells)
  basis_cells <- prediction_basis(parameters$options, mask_cells)

  # row of x_cell_years for each (cell, year)
  x_row <- matrix(NA_integer_, max(cell_years_index$cell_id), n_times)
  x_row[cbind(cell_years_index$cell_id, cell_years_index$year_id)] <-
    seq_len(nrow(cell_years_index))

  # assays sharing a (cell, type, year), and with the species model a share,
  # share a prediction, computed once
  keys <- paste(rows$cell_id, rows$type_id, rows$year_id)
  species <- species_on(parameters$options)
  if (species) {
    rows$share <- if (is.null(share)) {
      arabiensis_share(rows, parameters$options)
    } else {
      rep_len(share, nrow(rows))
    }
    keys <- paste(keys, rows$share)
  } else {
    stopifnot(is.null(share))
  }
  unique_rows <- rows[!duplicated(keys),
                      c("cell_id", "type_id", "year_id",
                        if (species) "share")]
  result <- matrix(NA_real_, n_draws, nrow(unique_rows))
  chunk_size <- max(1, floor(max_block / (n_times * n_draws)))

  for (k in sort(unique(unique_rows$type_id))) {
    cells_k <- sort(unique(unique_rows$cell_id[unique_rows$type_id == k]))
    for (cells in split(cells_k, ceiling(seq_along(cells_k) / chunk_size))) {
      target <- which(unique_rows$type_id == k &
                        unique_rows$cell_id %in% cells)
      years_keep <- sort(unique(unique_rows$year_id[target]))
      x_index <- x_row[cells, seq_len(max(years_keep)), drop = FALSE]
      stopifnot(!anyNA(x_index))
      x <- array(x_cell_years[as.vector(x_index), , drop = FALSE],
                 c(length(cells), max(years_keep), ncol(x_cell_years)))
      logit <- dynamical_trajectories(
        parameters, k,
        matrix(parameters$logit_init_relative[, cell_country[cells], k],
               nrow = n_draws),
        x, years_keep,
        x_init = parameters$x_cells_init[cells, , drop = FALSE],
        kdr = if (!is.null(kdr_cells)) kdr_cells[cells, , drop = FALSE],
        basis = if (!is.null(basis_cells)) {
          basis_cells[cells, , drop = FALSE]
        })
      for (t in years_keep) {
        at_t <- target[unique_rows$year_id[target] == t]
        columns <- match(unique_rows$cell_id[at_t], cells)
        logit_t <- logit[[as.character(t)]]
        result[, at_t] <- if (!species) {
          logit_t[, columns, drop = FALSE]
        } else {
          mixture_logit(logit_t$arabiensis[, columns, drop = FALSE],
                        logit_t$other[, columns, drop = FALSE],
                        unique_rows$share[at_t])
        }
      }
    }
  }

  result[, match(keys, keys[!duplicated(keys)]), drop = FALSE]
}
