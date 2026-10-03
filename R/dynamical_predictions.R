# Recompute the dynamical model's predictions, as logit draws, at arbitrary
# cell-years from its saved posterior draws, in plain R.
#
# fit_fold() (R/fit_validation_fold.R) only asked greta for predictions at the
# held-out assays, and a saved draws object cannot be resumed to ask for more.
# The two-stage model needs them at the training assays too, paired draw for
# draw with the stored test predictions, and on the map grid, so this
# recomputes them from the sampled parameters:
#
#   in the draws (constrained values of the arrays passed to model()):
#     beta_overall, sigma_overall, sigma_class   n_covs
#     beta_class_raw                             n_covs x n_classes
#     beta_type_raw                              n_covs x n_types
#     init_region_sd, init_country_sd            n_types
#     init_region_raw                            n_regions x n_types
#     init_country_raw                           n_countries x n_types
#   in the raw free-state draws only: logit_init_mean (n_types), which was not
#   passed to model() (see logit_init_mean_draws()).
#
# The recursion is solved in closed form on the logit scale. With
# q_{t+1} = q_t / (q_t + (1 - q_t) w_t), the odds of being susceptible are
# divided by w_t each year, so
#
#   logit q_t = logit q_0 - sum_{s <= t} log w_s
#
# exactly. greta.dynamics stores the state after each iteration and uses the
# fitness for iteration i at time slice i, so the state recorded for year index
# t has had the fitness of years 1..t applied. The state limits are the
# defaults (-Inf, Inf) and tol = 0 can never be met, so there is no clamping
# and no early stop to reproduce.
#
# Functions:
#   paired_draw_index(fold)                     draws paired with fold$p_draws
#   dynamical_parameter_draws(fold, classes_index, types, draw_index)
#   dynamical_logit_init(parameters, country, region)
#                                               logit q_0 per country, fresh
#                                               deviations for countries and
#                                               regions without data
#   dynamical_logit_cells(effect, logit_init, x, years_keep)
#                                               the recursion, for one type
#   dynamical_logit(parameters, rows, df, x_cell_years, cell_years_index)
#                                               at (cell, type, year) rows

# Draw indices into as.matrix(fold$draws) that the stored predictions
# correspond to: newer folds thin p_draws to `maximum` at fitting time with
# this rule, older folds are thinned by the same rule at scoring time. The
# rule is thin_draws() (R/validation_scoring.R, sourced by
# R/two_stage_helpers.R)
paired_draw_index <- function(fold, maximum = max_draws) {
  n_total <- sum(vapply(fold$draws, nrow, integer(1)))
  thin_draws(matrix(seq_len(n_total)), maximum)[, 1]
}

# One named parameter of a draws x parameter matrix, as a draws x
# dim(parameter) array. greta names elements by their R (column-major) index,
# e.g. "beta_type_raw[13,9]", so the index is parsed from the names
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

# logit_init_mean has no named columns in the draws, but greta still samples
# it, so its values are in the raw free-state draws,
# attr(draws, "model_info")$raw_draws, which concatenate every variable's free
# state in the dag's variable order. Its block is the one variable node that is
# not among the model's named targets. An untruncated normal has an identity
# free-state transform, so its raw value is its value; that is checked here on
# beta_overall, the same kind of variable
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

  named <- as.matrix(fold$draws)[draw_index, , drop = FALSE]
  beta_overall_named <- named[, grep("^beta_overall\\[", colnames(named)),
                              drop = FALSE]
  stopifnot(isTRUE(all.equal(unname(raw[, columns_of(target_tf["beta_overall"]),
                                        drop = FALSE]),
                             unname(beta_overall_named))))

  raw[, columns_of(untargeted), drop = FALSE]
}

# The parameter draws the predictions need, for the draws in `draw_index`:
#   effect_type   draws x n_covs x n_types, exp(beta_type)
#   init_*        the initial-condition components, combined by
#                 dynamical_logit_init()
dynamical_parameter_draws <- function(fold,
                                      classes_index,
                                      types,
                                      draw_index = paired_draw_index(fold),
                                      logit_init_mean = NULL) {

  # as.matrix() on an mcmc.list stacks chains in order, which is how fit_fold()
  # flattened the calculate() output that p_draws came from
  draws_matrix <- as.matrix(fold$draws)[draw_index, , drop = FALSE]

  # regression coefficients, doubly hierarchical: overall -> class -> type.
  # The sweeps recycle a draws x n_covs matrix down a draws x n_covs x k array
  beta_overall <- extract_parameter(draws_matrix, "beta_overall")
  sigma_overall <- extract_parameter(draws_matrix, "sigma_overall")
  sigma_class <- extract_parameter(draws_matrix, "sigma_class")
  beta_class <- extract_parameter(draws_matrix, "beta_class_raw") *
    as.vector(sigma_overall) + as.vector(beta_overall)
  effect_type <- exp(beta_class[, , classes_index, drop = FALSE] +
                       extract_parameter(draws_matrix, "beta_type_raw") *
                         as.vector(sigma_class))

  if (is.null(logit_init_mean)) {
    logit_init_mean <- logit_init_mean_draws(fold, draw_index)
  }
  stopifnot(identical(dim(logit_init_mean),
                      c(length(draw_index), length(types))))

  list(effect_type = effect_type,
       logit_init_mean = logit_init_mean,
       init_country_sd = extract_parameter(draws_matrix, "init_country_sd"),
       init_region_sd = extract_parameter(draws_matrix, "init_region_sd"),
       init_country_raw = extract_parameter(draws_matrix, "init_country_raw"),
       init_region_raw = extract_parameter(draws_matrix, "init_region_raw"),
       init_frac_min = ifelse(types == "DDT", 0.75, 0.9),
       n_draws = length(draw_index))
}

# Draws of logit q_0, the initial fraction susceptible, per country:
#   country  for each output country, its index among the fitted countries, or
#            NA for a country without data
#   region   for each output country, its index among the fitted regions;
#            an index beyond them is a region without data
# Countries and regions without data get fresh N(0, 1) raw deviations (the
# hierarchical model's prediction for a new one), one per posterior draw and
# type, drawn with the caller's RNG and shared by every country with the same
# index. Returns a draws x countries x types array.
dynamical_logit_init <- function(parameters, country, region) {
  stopifnot(length(country) == length(region), !anyNA(region))
  n_draws <- parameters$n_draws
  n_types <- length(parameters$init_frac_min)
  n_regions_fitted <- dim(parameters$init_region_raw)[2]

  fresh <- function(n) array(rnorm(n_draws * n * n_types),
                             c(n_draws, n, n_types))
  new_countries <- which(is.na(country))
  country_raw <- parameters$init_country_raw[, pmax(country, 1, na.rm = TRUE),
                                             , drop = FALSE]
  country_raw[, new_countries, ] <- fresh(length(new_countries))
  region_raw <- array(NA_real_, c(n_draws, max(region, n_regions_fitted),
                                  n_types))
  region_raw[, seq_len(n_regions_fitted), ] <- parameters$init_region_raw
  new_regions <- setdiff(unique(region), seq_len(n_regions_fitted))
  region_raw[, new_regions, ] <- fresh(length(new_regions))

  # logit of q_0 = min + range * ilogit(l), computed without forming q_0,
  # which can round to 1 when l is large: 1 - q_0 = range * ilogit(-l)
  logit_init <- array(NA_real_, c(n_draws, length(country), n_types))
  for (k in seq_len(n_types)) {
    l <- country_raw[, , k] * parameters$init_country_sd[, k] +
      region_raw[, region, k] * parameters$init_region_sd[, k] +
      parameters$logit_init_mean[, k]
    min_k <- parameters$init_frac_min[k]
    logit_init[, , k] <- log(min_k + (1 - min_k) * plogis(l)) -
      log(1 - min_k) - plogis(-l, log.p = TRUE)
  }
  logit_init
}

# The recursion for cells of one insecticide type:
#   effect      draws x n_covs, exp(beta_type) for the type
#   logit_init  draws x cells, logit q_0 at each cell's country
#   x           cells x years x n_covs covariates, year index 1 the baseline
#               year, covariates in the column order of x_cell_years
#   years_keep  year indices to return
# Returns a list named by years_keep of draws x cells logit q
dynamical_logit_cells <- function(effect, logit_init, x, years_keep) {
  n_cells <- dim(x)[1]
  stopifnot(ncol(logit_init) == n_cells, max(years_keep) <= dim(x)[2])
  cumulative <- 0
  out <- list()
  for (t in seq_len(max(years_keep))) {
    # log1p because the selection term can be tiny
    cumulative <- cumulative +
      log1p(effect %*% t(matrix(x[, t, ], nrow = n_cells)))
    if (t %in% years_keep) out[[as.character(t)]] <- logit_init - cumulative
  }
  out
}

# Draws of the dynamical model's logit predicted mortality at (cell, type,
# year) rows (cell_id, type_id, year_id, as in `df`), with the covariates of
# x_cell_years (one row per (cell_id, year_id) of cell_years_index). The
# initial condition of a cell is that of the country of its first record in
# the full `df`, as in fit_fold(). Returns a draws x nrow(rows) matrix, paired
# with thin_draws(fold$p_draws) when `parameters` are at paired_draw_index().
dynamical_logit <- function(parameters, rows, df, x_cell_years,
                            cell_years_index, max_block = 2.5e7) {

  n_draws <- parameters$n_draws
  n_times <- max(cell_years_index$year_id)
  if (any(rows$year_id > n_times | rows$year_id < 1)) {
    stop("rows fall outside the covariate years 1..", n_times)
  }

  # the lookups fit_fold() builds: each country's region and each cell's
  # country, from their first record in df
  country_region <- df %>%
    distinct(country_id, .keep_all = TRUE) %>%
    arrange(country_id)
  stopifnot(identical(country_region$country_id,
                      seq_len(dim(parameters$init_country_raw)[2])))
  logit_init <- dynamical_logit_init(parameters, country_region$country_id,
                                     country_region$region_id)
  cell_country <- df %>% distinct(cell_id, .keep_all = TRUE)
  stopifnot(all(rows$cell_id %in% cell_country$cell_id))

  # row of x_cell_years for each (cell, year)
  x_row <- matrix(NA_integer_, max(cell_years_index$cell_id), n_times)
  x_row[cbind(cell_years_index$cell_id, cell_years_index$year_id)] <-
    seq_len(nrow(cell_years_index))

  # assays sharing a (cell, type, year) share a prediction, computed once
  keys <- paste(rows$cell_id, rows$type_id, rows$year_id)
  unique_rows <- rows[!duplicated(keys), c("cell_id", "type_id", "year_id")]
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
      cell_init <- logit_init[, cell_country$country_id[match(cells,
                                                              cell_country$cell_id)],
                              k]
      logit <- dynamical_logit_cells(parameters$effect_type[, , k],
                                     matrix(cell_init, nrow = n_draws), x,
                                     years_keep)
      for (t in years_keep) {
        at_t <- target[unique_rows$year_id[target] == t]
        result[, at_t] <- logit[[as.character(t)]][
          , match(unique_rows$cell_id[at_t], cells), drop = FALSE]
      }
    }
  }

  result[, match(keys, keys[!duplicated(keys)]), drop = FALSE]
}
