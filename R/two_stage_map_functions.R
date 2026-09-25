# Helpers for mapping the two-stage model (#21) on the full prediction grid:
# the dynamical model's posterior draws at every mask cell, and the stage-A
# correction's latent fields projected from the mesh nodes to those cells.
#
# The dynamical part mirrors predict_batch() in R/predict.R, but in plain R on
# the logit scale (see R/dynamical_predictions.R for why the recursion has a
# closed form there), so that it can be paired draw for draw with the
# correction's cut-posterior mode shift, which predict.R's greta calculate()
# cannot. Functions only: source from the repo root after
# R/dynamical_predictions.R and R/two_stage_correction.R.


# covariates at the prediction cells ------------------------------------------

# All covariates the dynamical model uses, at the cells `cells` of the mask and
# for the years baseline_year..end_year, built the way predict.R builds
# x_cell_years_predict: the cubes are padded back to the baseline year and
# forward to end_year by repeating their first and last layers.
#
# Returned as a cells x years x 3 array of the time-varying covariates (nets,
# irs, pop, the column order of x_cell_years) and a cells x 10 matrix of the
# static crop covariates, rather than predict.R's long (cell, year) matrix: the
# long form is 53M rows x 13 columns (5.5 GB) for 1.48M cells and 36 years,
# whereas this keeps the crops once per cell and lets a chunk of cells be
# assembled one year at a time.
map_covariates <- function(cells, baseline_year = 1995, end_year = 2030) {

  read_cube <- function(file) {
    cube <- rast(file)
    cube <- pre_pad_cube(cube, baseline_year)
    cube <- post_pad_cube(cube, end_year)
    years <- as.numeric(str_sub(names(cube), start = -4L))
    cube <- cube[[years >= baseline_year & years <= end_year]]
    stopifnot(identical(as.numeric(str_sub(names(cube), start = -4L)),
                        as.numeric(baseline_year:end_year)))
    as.matrix(terra::extract(cube, cells))
  }

  nets <- read_cube("data/clean/net_use_cube.tif")
  irs <- read_cube("data/clean/irs_coverage_scaled_cube.tif")
  pop <- read_cube("data/clean/pop_scaled_cube.tif")

  time_varying <- array(NA_real_, c(length(cells), ncol(nets), 3),
                        dimnames = list(NULL, baseline_year:end_year,
                                        c("nets", "irs", "pop")))
  time_varying[, , 1] <- nets
  time_varying[, , 2] <- irs
  time_varying[, , 3] <- pop
  rm(nets, irs, pop)

  crops_group <- rast("data/clean/crop_group_scaled.tif")
  crops_all <- rast("data/clean/crop_scaled.tif")
  covs_flat <- c(crops_group,
                 crops_all$cotton,
                 crops_all$vegetables,
                 crops_all$rice)
  flat <- as.matrix(terra::extract(covs_flat, cells))

  list(time_varying = time_varying, flat = flat)
}


# initial conditions for every country ----------------------------------------

# Draws of logit q_0, the initial fraction susceptible, for every African
# country in the UNSD lookup (not only those with data), as in predict_batch():
# countries and regions without data get fresh N(0, 1) raw deviations, i.e. the
# hierarchical model's prediction for a new country or region. Returns a
# draws x countries x types array with dimnames[[2]] the country names.
#
# The country -> region mapping for prediction is the UNSD one, as in
# predict.R, where the fit took each country's region from its first record;
# the two agree for every observed country (checked in two_stage_maps.R).
#
# emulate_predict_fill = TRUE reproduces what predict.R actually computed, for
# checking against its saved maps only. There, the observed countries' (and
# regions') n x types matrix of raw deviations is assigned into the rows of a
# greta zeros() array with `pred[index, ] <- raw`; greta fills the target
# elements in column-major order from the source's elements in row-major
# order, so every country x type deviation lands on the wrong country and type
# (a transpose-and-refill of the block). The saved outputs/ir_maps therefore
# carry scrambled initial conditions; their covariate effects are unaffected.
map_logit_init <- function(draws_matrix,
                           logit_init_mean,
                           types,
                           countries,
                           regions,
                           lookup = country_region_lookup(),
                           seed = 1,
                           emulate_predict_fill = FALSE) {

  n_draws <- nrow(draws_matrix)
  n_types <- length(types)
  init_frac_min <- ifelse(types == "DDT", 0.75, 0.9)
  init_range <- 1 - init_frac_min

  all_countries <- unique(lookup$country_name)
  all_regions <- unique(lookup$region)

  init_region_sd <- extract_parameter(draws_matrix, "init_region_sd")
  init_country_sd <- extract_parameter(draws_matrix, "init_country_sd")
  init_region_raw <- extract_parameter(draws_matrix, "init_region_raw")
  init_country_raw <- extract_parameter(draws_matrix, "init_country_raw")

  if (emulate_predict_fill) {
    # row-major values into column-major positions, draw by draw
    refill <- function(a) {
      for (d in seq_len(dim(a)[1])) {
        a[d, , ] <- matrix(as.vector(t(a[d, , ])), dim(a)[2], dim(a)[3])
      }
      a
    }
    init_country_raw <- refill(init_country_raw)
    init_region_raw <- refill(init_region_raw)
  }

  # observed countries and regions keep their sampled deviations; the rest are
  # drawn from the prior, once per posterior draw and shared by every cell
  set.seed(seed)
  country_raw <- array(rnorm(n_draws * length(all_countries) * n_types),
                       c(n_draws, length(all_countries), n_types))
  region_raw <- array(rnorm(n_draws * length(all_regions) * n_types),
                      c(n_draws, length(all_regions), n_types))
  country_raw[, match(countries, all_countries), ] <- init_country_raw
  region_raw[, match(regions, all_regions), ] <- init_region_raw

  country_region <- match(lookup$region[match(all_countries,
                                              lookup$country_name)],
                          all_regions)

  logit_init <- array(NA_real_, c(n_draws, length(all_countries), n_types),
                      dimnames = list(NULL, all_countries, types))
  for (k in seq_len(n_types)) {
    l <- country_raw[, , k] * init_country_sd[, k] +
      region_raw[, country_region, k] * init_region_sd[, k] +
      logit_init_mean[, k]
    # logit of min + range * ilogit(l), without rounding q_0 to 1 (as in
    # dynamical_parameter_draws())
    logit_init[, , k] <- log(init_frac_min[k] + init_range[k] * plogis(l)) -
      log(init_range[k]) - plogis(-l, log.p = TRUE)
  }
  logit_init
}


# dynamical draws on a chunk of cells ------------------------------------------

# Logit q (predicted bioassay mortality) at a chunk of cells, for draws of one
# insecticide type, returned as a list over `years_keep` of cells x draws
# matrices.
#
#   effect       draws x n_covs: exp(beta_type) for this type
#   logit_init   cells x draws: logit q_0 at each cell's country
#   time_varying cells x years x 3 and flat cells x n_flat, from
#                map_covariates(), for this chunk
#
# logit q_t = logit q_0 - sum_{s <= t} log w_s, with the fitness of year 1
# (the baseline year) already applied to the year-1 state, exactly as
# dynamical_predictions() and greta.dynamics do.
dynamical_logit_chunk <- function(effect, logit_init, time_varying, flat,
                                  years, years_keep) {
  stopifnot(all(years_keep %in% years))
  last <- max(match(years_keep, years))
  effect_t <- t(effect)
  cumulative <- 0
  out <- list()
  for (t in seq_len(last)) {
    x_t <- cbind(time_varying[, t, ], flat)
    cumulative <- cumulative + log1p(x_t %*% effect_t)
    if (years[t] %in% years_keep) {
      out[[as.character(years[t])]] <- logit_init - cumulative
    }
  }
  out
}


# correction fields at the mesh nodes ----------------------------------------

# Node values of the correction omega + xi for each year in `years`, as
#   omega      n_nodes_omega x n_draws (the same in every year)
#   xi[[y]]    n_nodes_xi x n_draws
# from a matrix of latent vectors `theta` (n_latent x n_draws, in the fit's
# latent order). Years at or before t0 have xi = 0, years in (t0, T] read the
# fitted x, and later years run the AR(1) for eta forward from
# eta_T = x_T - x_{T-1}, adding `innovations[[h]]` (n_nodes_xi x n_draws,
# N(0, Q_eta^-1)) at horizon h and accumulating into xi, as
# predict_correction() does. With innovations = NULL the recursion carries the
# mean forward instead: eta_{T+h} = phi^h eta_T, so
#   xi_{T+h} = xi_T + eta_T phi (1 - phi^h) / (1 - phi),
# the plateauing forecast of the issue.
correction_node_fields <- function(fit, theta, years, innovations = NULL) {
  theta <- as.matrix(theta)
  n_nodes_xi <- fit$mesh_xi$n
  omega <- theta[fit$blocks$w_omega, , drop = FALSE]
  out <- list(omega = omega, xi = list())
  if (fit$variant != "omega_xi_u") {
    for (y in years) {
      out$xi[[as.character(y)]] <- matrix(0, n_nodes_xi, ncol(theta))
    }
    return(out)
  }

  x <- theta[fit$blocks$x, , drop = FALSE]
  x_col <- function(j) x[(j - 1) * n_nodes_xi + seq_len(n_nodes_xi), ,
                         drop = FALSE]
  xi_T <- x_col(fit$n_years)
  eta_T <- if (fit$n_years > 1) xi_T - x_col(fit$n_years - 1) else xi_T

  phi <- fit$hyper$phi
  for (y in years) {
    if (y <= fit$t0) {
      xi_y <- matrix(0, n_nodes_xi, ncol(theta))
    } else if (y <= fit$T) {
      xi_y <- x_col(y - fit$t0)
    } else {
      xi_y <- xi_T
      eta <- eta_T
      for (h in seq_len(y - fit$T)) {
        eta <- phi * eta
        if (!is.null(innovations)) {
          eta <- eta + sqrt(1 - phi ^ 2) * innovations[[h]]
        }
        xi_y <- xi_y + eta
      }
    }
    out$xi[[as.character(y)]] <- xi_y
  }
  out
}

# row standard deviations of a matrix
row_sds <- function(x) {
  n <- ncol(x)
  mu <- rowMeans(x)
  sqrt(pmax(rowSums(x ^ 2) - n * mu ^ 2, 0) / (n - 1))
}
