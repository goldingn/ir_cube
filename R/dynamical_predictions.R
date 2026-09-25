# Recompute the dynamical model's predicted mortality at arbitrary cell-years
# from a saved cross-validation fold, in plain R.
#
# fit_fold() (R/fit_validation_fold.R) only asked greta for predictions at the
# held-out assays, and a saved draws object cannot be resumed to ask for more.
# calculate(values = draws) on a reloaded fold would work, but it needs the
# whole greta/TensorFlow stack and rebuilds the full cells x types x years
# state for every draw. A second-stage model needs the first stage's predictions
# at the training assays too, paired draw for draw with the stored test
# predictions, so this recomputes them directly from the sampled parameters.
#
# Everything the recursion needs is either in the draws or is data:
#
#   in the draws (constrained values of the arrays passed to model()):
#     beta_overall, sigma_overall, sigma_class   n_covs
#     beta_class_raw                             n_covs x n_classes
#     beta_type_raw                              n_covs x n_types
#     init_region_sd, init_country_sd            n_types
#     init_region_raw                            n_regions x n_types
#     init_country_raw                           n_countries x n_types
#     rho_classes                                n_classes (not needed here)
#
#   data / fixed: x_cell_years and its (cell_id, year_id) index, the
#   type -> class index, the cell -> country and country -> region lookups
#   (both built from the full `df`, not the fold), and the prior-derived
#   constants init_frac_min.
#
#   in the raw free-state draws only: logit_init_mean (n_types), which was not
#   passed to model() and so has no named columns (see logit_init_mean_draws()).
#
# The recursion is solved in closed form on the logit scale. With
# q_{t+1} = q_t / (q_t + (1 - q_t) w_t), the odds of being susceptible are
# divided by w_t each year, so
#
#   logit q_t = logit q_0 - sum_{s <= t} log w_s
#
# exactly. greta.dynamics stores the state *after* each iteration and uses the
# fitness for iteration i at time slice i, so the state recorded for year_id t
# has had the fitness of years 1..t applied; the cumulative sum starts at the
# baseline year, not after it. The state limits are the defaults (-Inf, Inf)
# and tol = 0 can never be met, so there is no clamping and no early stop to
# reproduce.
#
# Working on the logit scale is also what makes this cheap: the cumulative
# selection depends on the cell and type but not the country, so it is one
# matrix product and one cumulative sum per type, and many assays share a
# cell-year.

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

# logit_init_mean is a free variable of the model that was not passed to
# model(), so it has no named columns in the draws. greta still samples it,
# though, so its values are in the raw free-state draws,
# attr(draws, "model_info")$raw_draws, which concatenate every variable's free
# state in the dag's variable order. The block is located from the saved dag:
# it is the one variable node that is not among the model's named targets. An
# untruncated normal has an identity free-state transform, so its raw value is
# its value; that is checked here on beta_overall, the same kind of variable.
logit_init_mean_draws <- function(fold, draw_index = paired_draw_index(fold)) {
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
  stopifnot(length(untargeted) == 1)

  raw <- do.call(rbind, lapply(model_info$raw_draws, as.matrix))
  raw <- raw[draw_index, , drop = FALSE]
  stopifnot(ncol(raw) == sum(sizes))

  # the identity-transform check
  named <- as.matrix(fold$draws)[draw_index, , drop = FALSE]
  beta_overall_named <- named[, grep("^beta_overall\\[", colnames(named)),
                              drop = FALSE]
  stopifnot(isTRUE(all.equal(unname(raw[, columns_of(target_tf["beta_overall"]),
                                        drop = FALSE]),
                             unname(beta_overall_named))))

  raw[, columns_of(untargeted), drop = FALSE]
}

# Draws of the fitness effects and the country-level logit initial state, for
# the draws in `draw_index`.
#
#   effect_type      draws x n_covs x n_types   exp(beta_type)
#   logit_init       draws x n_countries x n_types, logit of q_0
dynamical_parameter_draws <- function(fold,
                                      df,
                                      classes_index,
                                      types,
                                      draw_index = paired_draw_index(fold),
                                      logit_init_mean = NULL) {

  # as.matrix() on an mcmc.list stacks chains in order, which is exactly how
  # fit_fold() flattened the calculate() output that p_draws came from, so row
  # i here is the posterior sample behind row i of the unthinned p_draws
  draws_matrix <- as.matrix(fold$draws)[draw_index, , drop = FALSE]
  n_draws <- nrow(draws_matrix)
  n_types <- length(types)

  # regression coefficients, doubly hierarchical: overall -> class -> type
  beta_overall <- extract_parameter(draws_matrix, "beta_overall")
  sigma_overall <- extract_parameter(draws_matrix, "sigma_overall")
  sigma_class <- extract_parameter(draws_matrix, "sigma_class")
  beta_class_raw <- extract_parameter(draws_matrix, "beta_class_raw")
  beta_type_raw <- extract_parameter(draws_matrix, "beta_type_raw")

  n_covs <- dim(beta_type_raw)[2]
  beta_class <- beta_class_raw * as.vector(sigma_overall) +
    as.vector(beta_overall)
  # the sweeps above and below recycle a draws x n_covs matrix down the
  # draws x n_covs x k array, which is the same draw and covariate throughout
  effect_type <- exp(beta_class[, , classes_index, drop = FALSE] +
                       beta_type_raw * as.vector(sigma_class))

  # initial state: logit relative position above init_frac_min, with region and
  # country deviations
  init_frac_prior <- ifelse(types == "DDT", 0.9, 0.95)
  init_frac_min <- ifelse(types == "DDT", 0.75, 0.9)
  init_range <- 1 - init_frac_min

  init_region_sd <- extract_parameter(draws_matrix, "init_region_sd")
  init_country_sd <- extract_parameter(draws_matrix, "init_country_sd")
  init_region_raw <- extract_parameter(draws_matrix, "init_region_raw")
  init_country_raw <- extract_parameter(draws_matrix, "init_country_raw")

  if (is.null(logit_init_mean)) {
    logit_init_mean <- logit_init_mean_draws(fold, draw_index)
  }
  stopifnot(identical(dim(logit_init_mean), c(n_draws, n_types)))

  # the same lookups fit_fold() builds, from the full data rather than the fold
  country_region_index <- df %>%
    group_by(country_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(country_id, region_id) %>%
    arrange(country_id) %>%
    pull(region_id)

  n_countries <- dim(init_country_raw)[2]
  logit_init_relative <- array(NA_real_, c(n_draws, n_countries, n_types))
  for (k in seq_len(n_types)) {
    logit_init_relative[, , k] <-
      init_country_raw[, , k] * init_country_sd[, k] +
      init_region_raw[, country_region_index, k] * init_region_sd[, k] +
      logit_init_mean[, k]
  }

  # logit of q_0 = min + range * ilogit(l), computed without forming q_0, which
  # can round to 1 in double precision when l is large:
  #   1 - q_0 = range * ilogit(-l)
  logit_init <- logit_init_relative
  for (k in seq_len(n_types)) {
    l <- logit_init_relative[, , k]
    logit_init[, , k] <- log(init_frac_min[k] + init_range[k] * plogis(l)) -
      log(init_range[k]) - plogis(-l, log.p = TRUE)
  }

  list(effect_type = effect_type,
       logit_init = logit_init,
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
                                  max_block = 2.5e7) {

  parameters <- dynamical_parameter_draws(fold,
                                          df = df,
                                          classes_index = classes_index,
                                          types = types,
                                          draw_index = draw_index)
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
      logit_init <- matrix(parameters$logit_init[, target$country_id, k],
                           nrow = n_draws)
      unique_result[, target$column] <- plogis(logit_init - cumulative)
      rm(log_w, cumulative, logit_init)
    }
  }

  result <- unique_result[, match(keys, keys[unique_index]), drop = FALSE]
  rm(unique_result)
  result
}
