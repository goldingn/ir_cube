# Recompute the dynamical model's predicted mortality at arbitrary cell-years
# from a saved cross-validation fold (or the full fit's draws), in plain R.
#
# A saved draws object cannot be resumed to ask for more predictions.
# calculate(values = draws) on a reloaded fold would work, but it needs the
# whole greta/TensorFlow stack. A second-stage model needs the first stage's
# predictions at the training assays too, paired draw for draw with the stored
# test predictions, so this recomputes them directly from the sampled
# parameters.
#
# The transforms from the sampled parameters to the selection effects and the
# initial states are dynamical_terms() in R/dynamical_model.R, the function the
# greta model itself is built with, applied here to one draw at a time. The
# draws hold every variable passed to model(). Folds fitted before
# logit_init_mean was passed to model() hold it only in the raw free-state
# draws (see logit_init_mean_draws()).
#
# The recursion is solved in closed form on the logit scale, as in the model
# (see closed_form_states() in R/dynamical_model.R):
#
#   logit q_t = logit q_0 - sum_{s <= t} log w_s
#
# where the state recorded for year_id t has had the fitness of years 1..t
# applied; the cumulative sum starts at the baseline year, not after it.
#
# Working on the logit scale is also what makes this cheap: the cumulative
# selection depends on the cell and type but not the country, so it is one
# matrix product and one cumulative sum per type, and many assays share a
# cell-year.

source("R/dynamical_model.R")

# Draw indices into as.matrix(fold$draws) that the stored predictions
# correspond to. Newer folds thin p_draws to `stored_draws` at fitting time with
# round(seq(1, n, length.out = stored_draws)); older folds stored all draws and
# the scoring (thin_draws() in validation_metrics.R) applies the same rule, so
# one rule reproduces the pairing for both.
paired_draw_index <- function(fold, maximum = 2000) {
  n_total <- sum(vapply(fold$draws, nrow, integer(1)))
  if (n_total <= maximum) {
    return(seq_len(n_total))
  }
  round(seq(1, n_total, length.out = maximum))
}

# Pull one named parameter out of a draws x parameter matrix and return it as a
# draws x dim(parameter) array. greta names elements by their R (column-major)
# index, e.g. "beta_type_raw[13,9]", so the index is parsed from the names
# rather than assumed from the column order.
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

# Draws of logit_init_mean (draws x n_types). Fits since logit_init_mean was
# passed to model() have it as named columns. In older fits it was a free
# variable not passed to model(), so it has no named columns; greta still
# sampled it, though, so its values are in the raw free-state draws,
# attr(draws, "model_info")$raw_draws, which concatenate every variable's free
# state in the dag's variable order. The block is located from the saved dag:
# it is the one variable node that is not among the model's named targets. An
# untruncated normal has an identity free-state transform, so its raw value is
# its value; that is checked here on beta_overall, the same kind of variable.
logit_init_mean_draws <- function(fold, draw_index = paired_draw_index(fold)) {
  named <- as.matrix(fold$draws)[draw_index, , drop = FALSE]
  if (any(grepl("^logit_init_mean\\[", colnames(named)))) {
    return(extract_parameter(named, "logit_init_mean"))
  }

  model_info <- attr(fold$draws, "model_info")
  dag <- model_info$model$dag
  free <- dag$example_parameters(free = TRUE)
  sizes <- vapply(free, length, integer(1))
  ends <- cumsum(sizes)
  columns_of <- function(tf_name) {
    i <- match(tf_name, names(free))
    (ends[i] - sizes[i] + 1):ends[i]
  }

  node_names <- vapply(dag$node_list, function(node) node$unique_name,
                       character(1))
  tf_names <- dag$get_tf_names()
  target_nodes <- vapply(model_info$model$target_greta_arrays,
                         function(x) greta:::get_node(x)$unique_name,
                         character(1))
  target_tf <- tf_names[match(target_nodes, node_names)]
  names(target_tf) <- names(target_nodes)
  untargeted <- setdiff(names(free), target_tf)
  if (length(untargeted) != 1) {
    stop("expected logit_init_mean to be the one variable not passed to ",
         "model(), found ", length(untargeted))
  }

  raw <- do.call(rbind, lapply(model_info$raw_draws, as.matrix))
  raw <- raw[draw_index, , drop = FALSE]
  stopifnot(ncol(raw) == sum(sizes))

  # the identity-transform check
  beta_overall_named <- named[, grep("^beta_overall\\[", colnames(named)),
                              drop = FALSE]
  stopifnot(isTRUE(all.equal(unname(raw[, columns_of(target_tf["beta_overall"]),
                                        drop = FALSE]),
                             unname(beta_overall_named))))

  raw[, columns_of(untargeted), drop = FALSE]
}

# Every variable in a draws x parameter matrix, as a named list of draws x
# dim(variable) arrays, adding logit_init_mean if it is not among them.
variable_draws <- function(draws_matrix, logit_init_mean = NULL) {
  variable_names <- unique(sub("\\[.*$", "", colnames(draws_matrix)))
  variables <- lapply(setNames(nm = variable_names), extract_parameter,
                      draws_matrix = draws_matrix)
  if (!"logit_init_mean" %in% variable_names) {
    stopifnot(!is.null(logit_init_mean),
              nrow(logit_init_mean) == nrow(draws_matrix))
    variables$logit_init_mean <- logit_init_mean
  }
  variables
}

# The model options a fold was fitted with (dynamical_model_options()). Folds
# fitted before the options were saved used the defaults of the time, which are
# the settings dynamical_model_options() gives with every term off.
fold_options <- function(fold) {
  if (!is.null(fold$options)) {
    return(complete_dynamical_model_options(fold$options))
  }
  dynamical_model_options(rho = "class", mortality_floor = FALSE,
                          init_covariates = NULL,
                          selection_columns = selection_design_untrended(),
                          reversion = FALSE)
}

# One draw of each variable, as the arrays dynamical_terms() takes.
one_draw <- function(variables, i) {
  lapply(variables, function(a) {
    d <- dim(a)[-1]
    array(a[i + (seq_len(prod(d)) - 1) * nrow(a)], d)
  })
}

# dynamical_terms() for every draw, stacked as draws x dim(term) arrays.
dynamical_terms_draws <- function(variables, classes_index,
                                  country_region_index, types,
                                  terms = c("beta_type", "logit_init_country"),
                                  options = dynamical_model_options()) {
  n_draws <- nrow(variables[[1]])
  out <- NULL
  for (i in seq_len(n_draws)) {
    terms_i <- dynamical_terms(one_draw(variables, i),
                               classes_index = classes_index,
                               country_region_index = country_region_index,
                               types = types,
                               options = options)[terms]
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

# Draws of the fitness effects and the country-level logit initial state, for
# the draws in `draw_index`.
#
#   effect_type      draws x n_covs x n_types   exp(beta_type)
#   logit_init       draws x n_countries x n_types, logit of q_0
#   rho_types        draws x n_types, the observation overdispersion
#   mortality_floor  draws, the floor on mortality (NULL for none)
#   logit_init_relative  draws x n_countries x n_types, and
#   kappa_type       draws x n_types, the per-year change in logit resistance
#                    from reversion (<= 0; NULL for none, see
#                    reversion_kappa())
#   init_coef        draws x n_init_covs x n_types (NULL for none): the parts
#                    of the initial state when it has covariates (#19), the
#                    logit relative initial state of each country and the
#                    coefficients of the covariates on it
dynamical_parameter_draws <- function(fold,
                                      df,
                                      classes_index,
                                      types,
                                      draw_index = paired_draw_index(fold),
                                      logit_init_mean = NULL,
                                      options = fold_options(fold)) {

  # as.matrix() on an mcmc.list stacks chains in order, which is exactly how
  # fit_fold() flattened the calculate() output that p_draws came from, so row
  # i here is the posterior sample behind row i of the unthinned p_draws
  draws_matrix <- as.matrix(fold$draws)[draw_index, , drop = FALSE]
  n_draws <- nrow(draws_matrix)
  n_types <- length(types)

  if (is.null(logit_init_mean) &&
      !any(grepl("^logit_init_mean\\[", colnames(draws_matrix)))) {
    logit_init_mean <- logit_init_mean_draws(fold, draw_index)
  }
  variables <- variable_draws(draws_matrix, logit_init_mean)
  stopifnot(identical(dim(variables$logit_init_mean), c(n_draws, n_types)))

  # the same lookup the model builds, from the full data rather than the fold
  terms <- dynamical_terms_draws(
    variables,
    classes_index = classes_index,
    country_region_index = dynamical_lookups(df)$country_region_index,
    types = types,
    terms = c("beta_type", "logit_init_country", "rho_types",
              "logit_init_relative",
              if (!isFALSE(options$reversion)) "kappa_type"),
    options = options)

  effect_type <- exp(terms$beta_type)
  logit_init <- terms$logit_init_country
  n_covs <- dim(effect_type)[2]

  list(effect_type = effect_type,
       logit_init = logit_init,
       rho_types = matrix(terms$rho_types, n_draws),
       mortality_floor = if (isTRUE(options$mortality_floor)) {
         c(variables$mortality_floor)
       },
       logit_init_relative = terms$logit_init_relative,
       kappa_type = if (!isFALSE(options$reversion)) {
         matrix(terms$kappa_type, n_draws)
       },
       init_coef = if (!is.null(options$init_covariates)) {
         array(variables$init_coef,
               c(n_draws, length(options$init_covariates), n_types))
       },
       n_draws = n_draws,
       n_covs = n_covs)
}

# Draws of predicted mortality at arbitrary (cell, type, year) rows.
#
#   rows              data frame with cell_id, type_id, year_id (as in `df`),
#                     and optionally country_id; if absent, the country is
#                     looked up from `df` by cell_id
#   x_cell_years      covariate matrix, one row per (cell_id, year_id) in
#                     `cell_years_index`, as built by validation_covariates.R.
#                     To predict at cells outside the data, append their rows
#                     with new cell ids and supply country_id in `rows`.
#   draw_index        rows of as.matrix(fold$draws) to use; the default
#                     pairs with the stored or scored test predictions
#   x_cells_init      initial-state covariates, one row per cell_id
#                     (init_covariate_matrix()), for folds fitted with them;
#                     by default those saved with the fold
#
# Returns a draws x nrow(rows) matrix, paired row for row with
# thin_draws(fold$p_draws) under the default draw_index.
dynamical_predictions <- function(fold,
                                  rows,
                                  df,
                                  x_cell_years,
                                  cell_years_index,
                                  classes_index,
                                  types,
                                  draw_index = paired_draw_index(fold),
                                  max_block = 2.5e7,
                                  x_cells_init = fold$x_cells_init) {

  options <- fold_options(fold)
  parameters <- dynamical_parameter_draws(fold,
                                          df = df,
                                          classes_index = classes_index,
                                          types = types,
                                          draw_index = draw_index,
                                          options = options)
  x_init <- select_init_covariates(x_cells_init, options)
  init_min <- init_frac_constants(types)$min
  n_draws <- parameters$n_draws
  n_times <- max(cell_years_index$year_id)

  if (any(rows$year_id > n_times | rows$year_id < 1)) {
    stop("rows fall outside the covariate years 1..", n_times,
         "; pad the covariate cubes further to predict there")
  }

  if (!"country_id" %in% names(rows)) {
    cell_country_lookup <- df %>%
      distinct(cell_id, .keep_all = TRUE) %>%
      arrange(cell_id) %>%
      select(cell_id, country_id)
    rows$country_id <- cell_country_lookup$country_id[
      match(rows$cell_id, cell_country_lookup$cell_id)]
  }
  stopifnot(!anyNA(rows$country_id))

  # row of x_cell_years for each (cell, year), so the covariates can be read in
  # year order for any subset of cells without assuming the matrix's layout
  x_row <- matrix(NA_integer_, max(cell_years_index$cell_id), n_times)
  x_row[cbind(cell_years_index$cell_id, cell_years_index$year_id)] <-
    seq_len(nrow(cell_years_index))

  # the predictions depend on (cell, country, type, year) only; assays sharing
  # those share a column, computed once
  keys <- paste(rows$cell_id, rows$country_id, rows$type_id, rows$year_id)
  unique_index <- !duplicated(keys)
  unique_rows <- rows[unique_index, c("cell_id", "country_id", "type_id",
                                      "year_id")]
  unique_rows$column <- seq_len(nrow(unique_rows))
  unique_result <- matrix(NA_real_, n_draws, nrow(unique_rows))

  for (k in sort(unique(unique_rows$type_id))) {
    rows_k <- unique_rows[unique_rows$type_id == k, ]
    cells_k <- sort(unique(rows_k$cell_id))
    effect_k <- parameters$effect_type[, , k]   # draws x n_covs

    # chunk cells so the n_times x cells x draws block stays bounded
    chunk_size <- max(1, floor(max_block / (n_times * n_draws)))
    chunks <- split(cells_k, ceiling(seq_along(cells_k) / chunk_size))

    for (cells in chunks) {
      in_chunk <- rows_k$cell_id %in% cells
      target <- rows_k[in_chunk, ]
      # years after the latest one asked for in this chunk are not needed
      n_years <- max(target$year_id)
      x_index <- as.vector(t(x_row[cells, seq_len(n_years), drop = FALSE]))
      stopifnot(!anyNA(x_index))

      # draws x (years * cells), year fastest: log fitness, then cumulative
      # over years within each cell. Draws first keeps each year's slice
      # contiguous for the running sum and the final extraction. log1p because
      # the selection term can be tiny.
      log_w <- log1p(effect_k %*% t(x_cell_years[x_index, , drop = FALSE]))
      dim(log_w) <- c(n_draws, n_years, length(cells))
      for (t in seq_len(n_years)[-1]) {
        log_w[, t, ] <- log_w[, t, ] + log_w[, t - 1, ]
      }
      dim(log_w) <- c(n_draws, n_years * length(cells))

      cumulative <- log_w[, target$year_id +
                            (match(target$cell_id, cells) - 1) * n_years,
                          drop = FALSE]
      logit_init <- if (is.null(x_init)) {
        matrix(parameters$logit_init[, target$country_id, k], nrow = n_draws)
      } else {
        # the country's relative initial state plus the cell's covariate
        # effects, then the transform (as logit_init_relative_rows() and
        # logit_init_from_relative() in the model)
        relative <- matrix(parameters$logit_init_relative[, target$country_id,
                                                          k],
                           nrow = n_draws) +
          matrix(parameters$init_coef[, , k], nrow = n_draws) %*%
          t(x_init[target$cell_id, , drop = FALSE])
        logit_init_from_relative(relative, init_min[k])
      }
      # reversion: - t kappa in year t (see reversion_kappa())
      if (!is.null(parameters$kappa_type)) {
        cumulative <- cumulative +
          outer(parameters$kappa_type[, k], target$year_id)
      }
      # mortality is the fraction susceptible, with the floor applied (one
      # per draw, down the rows)
      unique_result[, target$column] <- floored_mortality(
        plogis(logit_init - cumulative), parameters$mortality_floor)
      rm(log_w, cumulative, logit_init)
    }
  }

  result <- unique_result[, match(keys, keys[unique_index]), drop = FALSE]
  rm(unique_result)
  result
}
