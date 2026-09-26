# What the second-stage correction says about the dynamical model (#21).
#
#   Rscript R/two_stage_covariate_diagnostics.R
#
# Run with OpenBLAS, as R/two_stage_maps.R, e.g.
#   LD_PRELOAD=.../libopenblas.so.0 OPENBLAS_NUM_THREADS=4 nice -n 10 Rscript ...
#
# The stage-A correction fitted to all data per insecticide type
# (R/two_stage_maps.R, outputs/two_stage/maps/<type>/fit.rds) is
#
#   z = m + omega(s) + xi(s, t) + u + e
#
# with m the dynamical model's posterior mean logit mortality, omega a static
# field (a correction to the initial state), xi an accumulated AR(1) field of
# annual anomalies eta (xi = 0 in 1995), and u pixel-year noise. Wherever the
# fitted fields line up with something the dynamical model knows about, the
# dynamical model is misspecified in a way it could fix itself. Three
# hypotheses:
#
#   (a) the initial state varies within countries with pre-2000 conditions
#       (#19): omega vs static covariates, within countries as well as across;
#   (b) the selection coefficients act through a transform of the covariates
#       (saturation, thresholds, lags, cumulative exposure): xi and eta vs the
#       time-varying covariates, and the shape of eta against net use;
#   (c) mortality has a floor (#14): the correction vs m itself.
#
# Everything is evaluated at the observed data (cells for omega, pixel-years for
# xi and eta), where the fields are informed; away from data they shrink to 0,
# so map-wide statistics are diluted, and are given only as a secondary check.
#
# Steps:
#
#   1. reload the full dynamical fit and compute, at every data cell-year of each
#      type, the posterior mean annual selection log w, its cumulative sum S
#      (the dynamical model's own trend) and the country initial state; and the
#      same at a sample of grid cells in the sampled countries (cached);
#   2. per type, project the fields' posterior mode to the assays, and refit
#      the latent fields (hyperparameters fixed at the fit's) to data simulated
#      under the null that the dynamical model is right, as a reference for
#      artefacts of the empirical logit and shrinkage (cached);
#   3. Spearman correlations with spatial block bootstrap CIs, to
#      outputs/two_stage/covariate_diagnostics_spearman.csv;
#   4. figures in figures/two_stage/: heatmaps, binned means, and the mean xi
#      trajectory by region.

n_boot <- 1000
# 500 km blocks are about the size of a country and ten times omega's range,
# but smaller than eta's (700-5000 km), so a coarser block size is reported too
block_sizes_km <- c(500, 1500)
n_null_sims <- 2
n_map_cells <- 20000
n_cores <- 4
flag_threshold <- 0.3
excess_threshold <- 0.2

report <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", sprintf(...), "\n")
  flush(stdout())
}

suppressMessages({
  sink("/dev/null")
  source("R/validation_folds.R")
  source("R/validation_covariates.R")
  sink()
})
source("R/dynamical_predictions.R")
source("R/two_stage_correction.R")
source("R/two_stage_map_functions.R")

cache_dir <- "outputs/two_stage/covariate_diagnostics"
figure_dir <- "figures/two_stage"
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

map_dir <- "outputs/two_stage/maps"
clamp <- 1e-12
safe_logit <- function(p) qlogis(pmin(pmax(p, clamp), 1 - clamp))
years_fit <- baseline_year:final_data_year
n_years_fit <- length(years_fit)

# the cell's country, as the dynamical fit looked it up (first record), so the
# initial state here is the one behind m
cell_country <- df %>%
  distinct(cell_id, .keep_all = TRUE) %>%
  arrange(cell_id) %>%
  pull(country_id)

# row of x_cell_years for each (cell_id, year_id)
x_row <- matrix(NA_integer_, max(cell_years_index$cell_id), n_years_fit)
x_row[cbind(cell_years_index$cell_id, cell_years_index$year_id)] <-
  seq_len(nrow(cell_years_index))

# country per grid cell, for the map-wide sample
mask_cells <- terra::cells(mask)
grid_country <- match(
  as.character(terra::extract(rast("data/clean/country_raster.tif"),
                              mask_cells)$country_name),
  countries)


# 1. dynamical components -----------------------------------------------------

# Posterior means are linear in log w, so the mean cumulative selection is the
# cumulative sum of the mean annual selection, and the mean logit mortality is
# init - S, up to the clamp m_ref applies where p rounds to 0 or 1
dynamical_file <- file.path(cache_dir, "dynamical_components.rds")
if (!file.exists(dynamical_file)) {

  fit_env <- new.env()
  load("temporary/fitted_model.RData", envir = fit_env)
  stopifnot(
    isTRUE(all.equal(fit_env$df, df)),
    identical(fit_env$types, types),
    identical(fit_env$countries, countries),
    isTRUE(all.equal(fit_env$x_cell_years, x_cell_years)),
    isTRUE(all.equal(fit_env$cell_years_index, cell_years_index))
  )
  fold <- list(draws = fit_env$draws)
  rm(fit_env)
  invisible(gc())

  draw_index <- paired_draw_index(fold)
  logit_init_mean <- logit_init_mean_draws(fold, draw_index)
  parameters <- dynamical_parameter_draws(fold, df = df,
                                          classes_index = classes_index,
                                          types = types,
                                          draw_index = draw_index,
                                          logit_init_mean = logit_init_mean)

  # m_ref exactly as R/two_stage_maps.R built it for the stage-A fits
  p_train <- dynamical_predictions(fold, select(df, -country_id), df,
                                   x_cell_years, cell_years_index,
                                   classes_index, types,
                                   draw_index = draw_index)
  m_ref <- colMeans(safe_logit(p_train))
  rm(p_train, fold)
  invisible(gc())

  # country x type posterior mean initial state (logit q0)
  init_mean <- apply(parameters$logit_init, c(2, 3), mean)
  dimnames(init_mean) <- list(countries, types)

  # mean annual log w at every data cell-year of each type, for all years
  # 1995-2024, so lags, increments and trajectories can be read off
  mean_log_w <- function(effect, x) {
    rows <- split(seq_len(nrow(x)), ceiling(seq_len(nrow(x)) / 20000))
    unlist(lapply(rows, function(r) {
      colMeans(log1p(effect %*% t(x[r, , drop = FALSE])))
    }))
  }
  selection <- bind_rows(lapply(seq_along(types), function(k) {
    cells_k <- sort(unique(df$cell_id[df$type_id == k]))
    index <- as.vector(t(x_row[cells_k, , drop = FALSE]))
    log_w <- mean_log_w(parameters$effect_type[, , k],
                        x_cell_years[index, , drop = FALSE])
    tibble(insecticide_type = types[k],
           cell_id = rep(cells_k, each = n_years_fit),
           year = rep(years_fit, length(cells_k)),
           log_w = log_w) %>%
      group_by(cell_id) %>%
      mutate(S = cumsum(log_w)) %>%
      ungroup()
  }))

  # the map-wide sample: grid cells in the countries with data for the type
  set.seed(19)
  map_sample <- lapply(seq_along(types), function(k) {
    type_countries <- unique(df$country_id[df$type_id == k])
    candidates <- which(grid_country %in% type_countries)
    sort(sample(candidates, min(n_map_cells, length(candidates))))
  })
  names(map_sample) <- types
  sample_union <- sort(unique(unlist(map_sample)))
  covariates_map <- map_covariates(mask_cells[sample_union], baseline_year,
                                   final_data_year)
  map_selection <- lapply(seq_along(types), function(k) {
    rows <- match(map_sample[[k]], sample_union)
    effect_t <- t(parameters$effect_type[, , k])
    log_w <- sapply(seq_len(n_years_fit), function(t) {
      x_t <- cbind(covariates_map$time_varying[rows, t, ],
                   covariates_map$flat[rows, , drop = FALSE])
      rowMeans(log1p(x_t %*% effect_t))
    })
    list(log_w = log_w, S = t(apply(log_w, 1, cumsum)))
  })
  names(map_selection) <- types

  saveRDS(list(m_ref = m_ref, init_mean = init_mean, selection = selection,
               map_sample = map_sample, sample_union = sample_union,
               covariates_map = covariates_map,
               map_selection = map_selection),
          dynamical_file)
  rm(parameters)
  invisible(gc())
}
dynamical <- readRDS(dynamical_file)
m_ref <- dynamical$m_ref
init_mean <- dynamical$init_mean
report("dynamical components loaded")

# the decomposition reproduces m_ref wherever the clamp does not bind
check <- tibble(insecticide_type = df$insecticide_type, cell_id = df$cell_id,
                year = df$year_start, country_id = df$country_id,
                m_ref = m_ref) %>%
  left_join(dynamical$selection, by = c("insecticide_type", "cell_id", "year"))
check$init <- init_mean[cbind(cell_country[check$cell_id],
                              match(check$insecticide_type, types))]
check_error <- abs(check$init - check$S - check$m_ref)[abs(check$m_ref) < 20]
report("init - S vs m_ref: max |difference| %.2g (|m_ref| < 20)",
       max(check_error))
stopifnot(max(check_error) < 0.05)
rm(check)


# 2. fields at the data, and the null refits ----------------------------------

hyperparameters <- read.csv(file.path(map_dir, "hyperparameters.csv"))
rho_table <- read.csv("outputs/bioassay_rho_hierarchical.csv")
rho_for_type <- setNames(rho_table$rho, rho_table$insecticide_type)

# stage-A training data of a type, in the order R/two_stage_maps.R used
type_training <- function(k) {
  rows_k <- which(df$type_id == k)
  train <- tibble(row = rows_k,
                  lon = df$longitude[rows_k], lat = df$latitude[rows_k],
                  year = df$year_start[rows_k], cell = df$cell[rows_k],
                  cell_id = df$cell_id[rows_k],
                  country_id = df$country_id[rows_k],
                  region = df$region[rows_k],
                  died = df$died[rows_k],
                  mosquito_number = df$mosquito_number[rows_k],
                  m = m_ref[rows_k])
  stage_a <- empirical_logit(train$died, train$mosquito_number,
                             rho_for_type[[types[k]]])
  train$z <- stage_a$z
  train$v <- stage_a$v
  train
}

# omega at each assay and xi / eta at each assay's year, from a latent vector
assay_fields <- function(fit, theta, coords, year) {
  A_omega <- mesh_basis(fit$mesh, coords)
  A_xi <- mesh_basis(fit$mesh_xi, coords)
  n_nodes_xi <- fit$mesh_xi$n
  x <- matrix(theta[fit$blocks$x], n_nodes_xi, fit$n_years)
  # xi at every year 1995..T at the assay locations; column 1 is t0 (0)
  xi_all <- cbind(0, as.matrix(A_xi %*% x))
  j <- year - fit$t0 + 1
  tibble(omega = as.vector(A_omega %*% theta[fit$blocks$w_omega]),
         xi = xi_all[cbind(seq_along(year), j)],
         eta = ifelse(j > 1, xi_all[cbind(seq_along(year), j)] -
                        xi_all[cbind(seq_along(year), pmax(j - 1, 1))], 0))
}

# Refit the latent fields with the hyperparameters fixed at the fit's: the
# posterior is then Gaussian and the inner Newton solve gives its mode. With
# the observed z this must reproduce the saved mode
refit_mode <- function(fit, train, z, v) {
  coords <- coords_km(train)
  design <- correction_design(fit$mesh, fit$mesh_xi, coords, train$year,
                              fit$t0, fit$T)
  u_index <- fit$pixel_years$u_index[match(paste(train$cell, train$year),
                                           paste(fit$pixel_years$cell,
                                                 fit$pixel_years$year))]
  priors <- correction_priors()
  data <- list(z = z, v = v, m = train$m,
               A_omega = design$A_omega, A_xi = design$A_xi,
               u_index = as.integer(u_index - 1),
               spde = correction_fem(fit$mesh),
               spde_xi = correction_fem(fit$mesh_xi),
               include_xi = 1L,
               pc_omega = priors$pc_omega, pc_eta = priors$pc_eta,
               pc_tau = priors$pc_tau,
               persistence_prior = priors$persistence_prior)
  parameters <- c(list(w_omega = rep(0, fit$mesh$n),
                       x = matrix(0, fit$mesh_xi$n, fit$n_years),
                       u = rep(0, nrow(fit$pixel_years))),
                  as.list(fit$opt$par))
  obj <- correction_adfun(data, parameters, "omega_xi_u", fix_hyper = TRUE)
  obj$fn(obj$par)
  obj$env$last.par[obj$env$random]
}

fields_file <- file.path(cache_dir, "fields.rds")
if (!file.exists(fields_file)) {
  fields <- list()
  for (k in seq_along(types)) {
    type <- types[k]
    fit <- readRDS(file.path(map_dir, type, "fit.rds"))
    fit$variant <- "omega_xi_u"
    train <- type_training(k)
    stopifnot(identical(fit$pixel_years$cell,
                        distinct(train, cell, year)$cell))
    coords <- coords_km(train)

    observed <- assay_fields(fit, fit$mode, coords, train$year)

    # trajectories of xi at every data cell of the type, 1995..T, at the mean
    # location of the cell's assays
    cell_coords <- train %>%
      mutate(x_km = coords[, 1], y_km = coords[, 2]) %>%
      group_by(cell_id) %>%
      summarise(x_km = mean(x_km), y_km = mean(y_km),
                lon = mean(lon), lat = mean(lat),
                region = first(region), .groups = "drop")
    x <- matrix(fit$mode[fit$blocks$x], fit$mesh_xi$n, fit$n_years)
    xi_cells <- cbind(0, as.matrix(
      mesh_basis(fit$mesh_xi, cbind(cell_coords$x_km, cell_coords$y_km)) %*%
        x))
    colnames(xi_cells) <- fit$t0:fit$T

    # the null: data from the dynamical model's m, with u ~ N(0, tau) per
    # pixel-year and beta-binomial assay noise at the type's rho, refitted
    # with the same hyperparameters. Its fields are pure artefact (the
    # empirical logit's bounds, shrinkage), so any association they show with
    # m or the covariates is a reference level for the observed ones
    set.seed(4000 + k)
    time_check <- system.time(
      check_mode <- refit_mode(fit, train, train$z, train$v))
    report("%-18s refit at fixed hyperparameters reproduces the mode: max |diff| %.2g (%.0f s)",
           type, max(abs(check_mode - fit$mode)), time_check[["elapsed"]])
    stopifnot(max(abs(check_mode - fit$mode)) < 1e-4)
    rho <- rho_for_type[[type]]
    null <- lapply(seq_len(n_null_sims), function(i) {
      u_index <- match(paste(train$cell, train$year),
                       paste(fit$pixel_years$cell, fit$pixel_years$year))
      u <- rnorm(nrow(fit$pixel_years), 0, fit$hyper$tau)[u_index]
      p <- pmin(pmax(plogis(train$m + u), 1e-9), 1 - 1e-9)
      died <- extraDistr::rbbinom(nrow(train), train$mosquito_number,
                                  p * (1 - rho) / rho,
                                  (1 - p) * (1 - rho) / rho)
      stage_a <- empirical_logit(died, train$mosquito_number, rho)
      mode <- refit_mode(fit, train, stage_a$z, stage_a$v)
      assay_fields(fit, mode, coords, train$year) %>%
        mutate(sim = i, z = stage_a$z)
    })

    fields[[type]] <- list(
      assays = bind_cols(train, observed),
      null = null,
      cell_coords = cell_coords,
      xi_cells = xi_cells,
      fit_light = list(mesh = fit$mesh, mesh_xi = fit$mesh_xi,
                       mode = fit$mode, blocks = fit$blocks, t0 = fit$t0,
                       T = fit$T, n_years = fit$n_years))
    report("%-18s fields and %i null refits done", type, n_null_sims)
  }
  saveRDS(fields, fields_file)
}
fields <- readRDS(fields_file)


# 3. analysis units and covariates ------------------------------------------

classes_of <- df %>% distinct(insecticide_type, insecticide_class)
class_of <- setNames(classes_of$insecticide_class, classes_of$insecticide_type)

x_at <- function(cell_id, year, column) {
  x_cell_years[cbind(x_row[cbind(cell_id, year - baseline_year + 1)],
                     match(column, colnames(x_cell_years)))]
}
# net use and IRS summed over 2000..year (the cubes are padded back from
# 2000 / 1997 by repeating the first layer, so earlier years are not data)
cumulative_since_2000 <- function(cell_id, year, column) {
  vapply(seq_along(cell_id), function(i) {
    if (year[i] < 2000) return(0)
    sum(x_at(rep(cell_id[i], year[i] - 1999), 2000:year[i], column))
  }, numeric(1))
}

static_covariates <- c("all crops", "cereal crops", "root crops",
                       "pulse crops", "oil crops", "fibre crops",
                       "other crops", "cotton", "vegetables", "rice",
                       "pop_2000", "nets_2000_02", "irs_2000_02",
                       "latitude", "longitude", "init", "m_mean")
dynamic_covariates <- c("nets", "irs", "pop", "nets_cum", "irs_cum",
                        "nets_lag1", "nets_lag2", "nets_lag3", "nets_change",
                        "log_w", "S", "m", "year")

# one row per (type, pixel-year), the unit of xi and eta; the observed and
# null fields are averaged over the pixel-year's assays
pixel_years <- bind_rows(lapply(types, function(type) {
  f <- fields[[type]]
  null <- bind_rows(f$null) %>%
    mutate(assay = rep(seq_len(nrow(f$assays)), n_null_sims)) %>%
    group_by(assay) %>%
    summarise(across(c(omega, xi, eta), mean), .groups = "drop")
  f$assays %>%
    mutate(omega_null = null$omega, xi_null = null$xi, eta_null = null$eta,
           smooth = omega + xi, smooth_null = omega_null + xi_null,
           x_km = coords_km(.)[, 1], y_km = coords_km(.)[, 2]) %>%
    group_by(cell_id, year) %>%
    summarise(insecticide_type = type,
              country_id = first(country_id), region = first(region),
              across(c(omega, xi, eta, smooth, omega_null, xi_null,
                       eta_null, smooth_null, m, z, x_km, y_km), mean),
              n_assays = n(), .groups = "drop")
})) %>%
  left_join(dynamical$selection, by = c("insecticide_type", "cell_id",
                                        "year")) %>%
  mutate(
    nets = x_at(cell_id, year, "nets"),
    irs = x_at(cell_id, year, "irs"),
    pop = x_at(cell_id, year, "pop"),
    nets_lag1 = x_at(cell_id, pmax(year - 1, baseline_year), "nets"),
    nets_lag2 = x_at(cell_id, pmax(year - 2, baseline_year), "nets"),
    nets_lag3 = x_at(cell_id, pmax(year - 3, baseline_year), "nets"),
    nets_change = nets - nets_lag1,
    nets_cum = cumulative_since_2000(cell_id, year, "nets"),
    irs_cum = cumulative_since_2000(cell_id, year, "irs"),
    insecticide_class = class_of[insecticide_type]
  )

# one row per (type, cell), the unit of omega
cells <- pixel_years %>%
  group_by(insecticide_type, insecticide_class, cell_id) %>%
  summarise(country_id = first(country_id), region = first(region),
            omega = weighted.mean(omega, n_assays),
            omega_null = weighted.mean(omega_null, n_assays),
            m_mean = weighted.mean(m, n_assays),
            x_km = mean(x_km), y_km = mean(y_km), .groups = "drop") %>%
  mutate(
    pop_2000 = x_at(cell_id, rep(2000, n()), "pop"),
    nets_2000_02 = (x_at(cell_id, rep(2000, n()), "nets") +
                      x_at(cell_id, rep(2001, n()), "nets") +
                      x_at(cell_id, rep(2002, n()), "nets")) / 3,
    irs_2000_02 = (x_at(cell_id, rep(2000, n()), "irs") +
                     x_at(cell_id, rep(2001, n()), "irs") +
                     x_at(cell_id, rep(2002, n()), "irs")) / 3,
    init = init_mean[cbind(cell_country[cell_id],
                           match(insecticide_type, types))]
  )
for (column in static_covariates[1:10]) {
  cells[[column]] <- x_at(cells$cell_id, rep(2000, nrow(cells)), column)
}
cell_lonlat <- bind_rows(lapply(types, function(type) {
  fields[[type]]$cell_coords %>%
    transmute(insecticide_type = type, cell_id, longitude = lon,
              latitude = lat)
}))
cells <- left_join(cells, cell_lonlat, by = c("insecticide_type", "cell_id"))

# xi = 0 identically in 1995 (the dynamical model's start), so those
# pixel-years carry no information on xi or eta
pixel_years <- filter(pixel_years, year > baseline_year)
report("%i cells and %i pixel-years (after 1995) over the types",
       nrow(cells), nrow(pixel_years))


# 4. Spearman correlations with a spatial block bootstrap -------------------

# Groups: each type, and the types of a class pooled. Pooled, ranks are taken
# within type so the types' different scales do not create correlation, and
# the correlation is of the within-type ranks.
groups <- c(setNames(as.list(types), types),
            list("Pyrethroids (pooled)" =
                   types[class_of[types] == "Pyrethroids"],
                 "Organophosphates (pooled)" =
                   types[class_of[types] == "Organophosphates"]))

# Ranks within strata, then centred within `centre_by`: for the plain
# correlation centre_by = strata; for the within-group one (country for
# omega, year for xi and eta) it removes the between-group variation, which
# is what the dynamical model's country initial states and common trend
# already have a chance to explain. Pearson correlation of the centred ranks
# is then Spearman's rho (stratified / within-group)
centred_ranks <- function(M, strata, centre_by) {
  R <- M
  for (s in unique(strata)) {
    rows <- which(strata == s)
    R[rows, ] <- apply(M[rows, , drop = FALSE], 2, rank)
  }
  means <- rowsum(R, centre_by, reorder = FALSE) /
    as.vector(table(factor(centre_by, levels = unique(centre_by))))
  R - means[match(centre_by, unique(centre_by)), , drop = FALSE]
}

rank_correlations <- function(Y, X, strata, within) {
  plain <- suppressWarnings(cor(centred_ranks(Y, strata, strata),
                                centred_ranks(X, strata, strata)))
  inside <- suppressWarnings(cor(centred_ranks(Y, strata, within),
                                 centred_ranks(X, strata, within)))
  list(plain = plain, within = inside)
}

# block ids on the equal-area projection, in km
block_id <- function(x_km, y_km, size) {
  paste(floor(x_km / size), floor(y_km / size))
}

# Spearman rho (and within-group rho) of each field in Y with each covariate
# in X, with percentile CIs from resampling spatial blocks with replacement.
# Y_null holds the null refits' fields (same columns, same rows): the same
# statistics on them are the level of association the artefacts alone produce
# (the empirical logit's bounds wherever m is extreme, and shrinkage), and the
# excess rho - rho_null is bootstrapped paired, on the same resamples
block_bootstrap <- function(Y, Y_null, X, strata, within, x_km, y_km,
                            size, n_boot, seed) {
  blocks <- block_id(x_km, y_km, size)
  rows_by_block <- split(seq_along(blocks), blocks)
  fields <- colnames(Y)
  YY <- cbind(Y, Y_null)
  colnames(YY) <- c(fields, paste0(fields, "_null"))
  statistics <- function(idx) {
    r <- rank_correlations(YY[idx, , drop = FALSE], X[idx, , drop = FALSE],
                           strata[idx], within[idx])
    null_rows <- paste0(fields, "_null")
    list(plain = r$plain[fields, , drop = FALSE],
         within = r$within[fields, , drop = FALSE],
         plain_null = r$plain[null_rows, , drop = FALSE],
         within_null = r$within[null_rows, , drop = FALSE],
         excess = r$plain[fields, , drop = FALSE] -
           r$plain[null_rows, , drop = FALSE],
         excess_within = r$within[fields, , drop = FALSE] -
           r$within[null_rows, , drop = FALSE])
  }
  estimate <- statistics(seq_len(nrow(Y)))
  set.seed(seed)
  seeds <- sample.int(1e8, n_boot)
  boot <- parallel::mclapply(seeds, function(s) {
    set.seed(s)
    statistics(unlist(rows_by_block[sample(length(rows_by_block),
                                           replace = TRUE)],
                      use.names = FALSE))
  }, mc.cores = n_cores)
  quantiles <- function(version, p) {
    arr <- simplify2array(lapply(boot, `[[`, version))
    out <- apply(arr, c(1, 2), quantile, p, na.rm = TRUE)
    dim(out) <- dim(estimate$plain)
    dimnames(out) <- dimnames(estimate$plain)
    out
  }
  long <- function(M, name) {
    rownames(M) <- fields
    as_tibble(M, rownames = "field") %>%
      pivot_longer(-field, names_to = "covariate", values_to = name)
  }
  reduce(list(long(estimate$plain, "rho"),
              long(quantiles("plain", 0.025), "rho_lo"),
              long(quantiles("plain", 0.975), "rho_hi"),
              long(estimate$within, "rho_within"),
              long(quantiles("within", 0.025), "rho_within_lo"),
              long(quantiles("within", 0.975), "rho_within_hi"),
              long(estimate$plain_null, "rho_null"),
              long(estimate$within_null, "rho_within_null"),
              long(estimate$excess, "excess"),
              long(quantiles("excess", 0.025), "excess_lo"),
              long(quantiles("excess", 0.975), "excess_hi"),
              long(estimate$excess_within, "excess_within"),
              long(quantiles("excess_within", 0.025), "excess_within_lo"),
              long(quantiles("excess_within", 0.975), "excess_within_hi")),
         left_join, by = c("field", "covariate")) %>%
    mutate(n_units = nrow(Y), n_blocks = length(rows_by_block))
}

spearman_file <- "outputs/two_stage/covariate_diagnostics_spearman.csv"
results <- list()
for (group in names(groups)) {
  seed <- 100 + match(group, names(groups))
  for (unit in c("cell", "pixel_year")) {
    if (unit == "cell") {
      d <- filter(cells, insecticide_type %in% groups[[group]])
      Y <- as.matrix(d[, "omega"])
      Y_null <- as.matrix(d[, "omega_null"])
      X <- as.matrix(d[, static_covariates])
      within <- paste(d$insecticide_type, d$country_id)
    } else {
      d <- filter(pixel_years, insecticide_type %in% groups[[group]])
      Y <- as.matrix(d[, c("xi", "eta", "smooth")])
      Y_null <- as.matrix(d[, c("xi_null", "eta_null", "smooth_null")])
      X <- as.matrix(d[, dynamic_covariates])
      within <- paste(d$insecticide_type, d$year)
    }
    strata <- d$insecticide_type
    for (size in block_sizes_km) {
      results[[length(results) + 1]] <-
        block_bootstrap(Y, Y_null, X, strata, within, d$x_km, d$y_km, size,
                        n_boot, seed) %>%
        mutate(group = group, unit = unit, block_km = size,
               .before = everything())
    }
  }
  report("%-26s bootstrapped", group)
}
spearman <- bind_rows(results) %>%
  mutate(ci_excludes_0 = rho_lo > 0 | rho_hi < 0,
         within_ci_excludes_0 = rho_within_lo > 0 | rho_within_hi < 0,
         flag = abs(rho) >= flag_threshold & ci_excludes_0,
         flag_within = abs(rho_within) >= flag_threshold &
           within_ci_excludes_0,
         # beyond the artefact: the excess over the null, at a lower
         # threshold since it is a difference of correlations
         excess_ci_excludes_0 = excess_lo > 0 | excess_hi < 0,
         excess_within_ci_excludes_0 = excess_within_lo > 0 |
           excess_within_hi < 0,
         flag_excess = abs(excess) >= excess_threshold & excess_ci_excludes_0,
         flag_excess_within = abs(excess_within) >= excess_threshold &
           excess_within_ci_excludes_0)


# map-wide secondary check: the fields at a sample of grid cells in the
# type's sampled countries (and every year for xi), no CIs
map_rows <- list()
for (k in seq_along(types)) {
  type <- types[k]
  f <- fields[[type]]$fit_light
  sample_k <- dynamical$map_sample[[type]]
  rows <- match(sample_k, dynamical$sample_union)
  xy <- terra::xyFromCell(mask, mask_cells[sample_k])
  coords <- project_km(xy[, 1], xy[, 2])
  omega <- as.vector(mesh_basis(f$mesh, coords) %*% f$mode[f$blocks$w_omega])
  x <- matrix(f$mode[f$blocks$x], f$mesh_xi$n, f$n_years)
  xi <- cbind(0, as.matrix(mesh_basis(f$mesh_xi, coords) %*% x))
  eta <- xi[, -1] - xi[, -ncol(xi)]
  xi <- xi[, -1]
  years_k <- (f$t0 + 1):f$T
  j <- years_k - baseline_year + 1
  tv <- dynamical$covariates_map$time_varying[rows, , , drop = FALSE]
  flat <- dynamical$covariates_map$flat[rows, , drop = FALSE]
  init <- init_mean[grid_country[sample_k], k]
  S <- dynamical$map_selection[[type]]$S
  cumulative <- function(v) {
    out <- t(apply(v[, (2000 - baseline_year + 1):n_years_fit], 1, cumsum))
    cbind(matrix(0, nrow(v), 2000 - baseline_year), out)
  }
  nets_cum <- cumulative(tv[, , "nets"])
  irs_cum <- cumulative(tv[, , "irs"])
  lag <- function(v, l) v[, pmax(j - l, 1)]
  X_static <- cbind(flat, pop_2000 = tv[, 6, "pop"],
                    nets_2000_02 = rowMeans(tv[, 6:8, "nets"]),
                    irs_2000_02 = rowMeans(tv[, 6:8, "irs"]),
                    latitude = xy[, 2], longitude = xy[, 1], init = init,
                    m_mean = rowMeans(init - S[, j]))
  X_static <- X_static[, static_covariates]
  X_dynamic <- cbind(nets = as.vector(tv[, j, "nets"]),
                     irs = as.vector(tv[, j, "irs"]),
                     pop = as.vector(tv[, j, "pop"]),
                     nets_cum = as.vector(nets_cum[, j]),
                     irs_cum = as.vector(irs_cum[, j]),
                     nets_lag1 = as.vector(lag(tv[, , "nets"], 1)),
                     nets_lag2 = as.vector(lag(tv[, , "nets"], 2)),
                     nets_lag3 = as.vector(lag(tv[, , "nets"], 3)),
                     nets_change = as.vector(tv[, j, "nets"] -
                                               lag(tv[, , "nets"], 1)),
                     log_w = as.vector(dynamical$map_selection[[type]]$log_w[, j]),
                     S = as.vector(S[, j]),
                     m = as.vector(init - S[, j]),
                     year = rep(years_k, each = length(rows)))
  spearman_map <- function(y, X, field) {
    tibble(group = type, field = field, covariate = colnames(X),
           rho_mapwide = as.vector(suppressWarnings(
             cor(y, X, method = "spearman", use = "pairwise.complete.obs"))))
  }
  map_rows[[type]] <- bind_rows(
    spearman_map(omega, X_static, "omega"),
    spearman_map(as.vector(xi), X_dynamic, "xi"),
    spearman_map(as.vector(eta), X_dynamic, "eta"),
    spearman_map(as.vector(xi + omega), X_dynamic, "smooth"))
}
spearman <- spearman %>%
  left_join(bind_rows(map_rows), by = c("group", "field", "covariate"))
write.csv(spearman, spearman_file, row.names = FALSE)
report("Spearman table written to %s", spearman_file)


# 5. figures ------------------------------------------------------------------

group_levels <- c(types[order(class_of[types], types)],
                  "Pyrethroids (pooled)", "Organophosphates (pooled)")
covariate_labels <- c(
  "all crops" = "all crops", "cereal crops" = "cereal crops",
  "root crops" = "root crops", "pulse crops" = "pulse crops",
  "oil crops" = "oil crops", "fibre crops" = "fibre crops",
  "other crops" = "other crops", cotton = "cotton",
  vegetables = "vegetables", rice = "rice",
  pop_2000 = "population 2000", nets_2000_02 = "net use 2000-02",
  irs_2000_02 = "IRS 2000-02", latitude = "latitude",
  longitude = "longitude", init = "country initial state",
  m_mean = "m (cell mean)",
  nets = "net use", irs = "IRS", pop = "population",
  nets_cum = "cumulative net use", irs_cum = "cumulative IRS",
  nets_lag1 = "net use, lag 1", nets_lag2 = "net use, lag 2",
  nets_lag3 = "net use, lag 3", nets_change = "net use change",
  log_w = "annual selection log w", S = "cumulative selection S",
  m = "m", year = "year")

# Diverging blue (positive: the data are more susceptible than the model) to
# red, with a light grey midpoint, as the correction maps
heatmap <- function(data, title, subtitle, file, height) {
  data <- data %>%
    filter(block_km == 500) %>%
    select(group, field, covariate, rho, ci_excludes_0, flag,
           rho_within, within_ci_excludes_0, flag_within,
           excess, excess_ci_excludes_0, flag_excess) %>%
    pivot_longer(c(rho, rho_within, excess), names_to = "version",
                 values_to = "value") %>%
    mutate(
      excludes = case_when(version == "rho" ~ ci_excludes_0,
                           version == "rho_within" ~ within_ci_excludes_0,
                           TRUE ~ excess_ci_excludes_0),
      flagged = case_when(version == "rho" ~ flag,
                          version == "rho_within" ~ flag_within,
                          TRUE ~ flag_excess),
      version = factor(recode(version, rho = "Spearman rho",
                              rho_within = "within group (country / year)",
                              excess = "rho minus null rho"),
                       levels = c("Spearman rho",
                                  "within group (country / year)",
                                  "rho minus null rho")),
      group = factor(group, levels = group_levels),
      covariate = factor(covariate_labels[covariate],
                         levels = rev(unique(covariate_labels))),
      label = ifelse(excludes, sprintf("%.2f", value), ""))
  ggplot(data, aes(x = group, y = covariate, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_tile(data = filter(data, flagged %in% TRUE), fill = NA,
              colour = "black", linewidth = 0.6) +
    geom_text(aes(label = label), size = 2.2, colour = grey(0.15)) +
    scale_fill_gradient2(name = "rho", low = "#b2182b", mid = "#f2f2f2",
                         high = "#2166ac", midpoint = 0,
                         limits = c(-0.8, 0.8), oob = scales::squish) +
    facet_grid(field ~ version, scales = "free_y", space = "free_y") +
    labs(x = NULL, y = NULL, title = title,
         subtitle = paste(strwrap(subtitle, 150), collapse = "\n")) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          panel.grid = element_blank(),
          plot.subtitle = element_text(size = 7.5))
  ggsave(file, bg = "white", width = 15, height = height, dpi = 200)
}

heatmap(
  filter(spearman, field == "omega"),
  title = "Static correction omega vs static covariates, at the data cells",
  subtitle = paste("Spearman rho over the type's data cells; within group =",
                   "ranks centred within country. Values shown where the",
                   "500 km block-bootstrap 95% CI excludes 0; boxed where",
                   "also |rho| >= 0.3 (|difference| >= 0.2 for rho minus",
                   "null rho, the null being fields refitted to data",
                   "simulated from m). Pooled columns rank within type."),
  file = file.path(figure_dir, "covariate_diagnostics_omega_heatmap.png"),
  height = 5)
heatmap(
  filter(spearman, field %in% c("xi", "eta")) %>%
    mutate(field = factor(field, levels = c("xi", "eta"))),
  title = paste("Accumulated anomaly xi and annual anomaly eta vs",
                "time-varying covariates, at the data pixel-years"),
  subtitle = paste("Spearman rho over the type's pixel-years after 1995;",
                   "within group = ranks centred within year. Values shown",
                   "where the 500 km block-bootstrap 95% CI excludes 0;",
                   "boxed where also |rho| >= 0.3 (|difference| >= 0.2 for",
                   "rho minus null rho, the null being fields refitted to",
                   "data simulated from m). log w, S and m are the",
                   "dynamical model's posterior means."),
  file = file.path(figure_dir, "covariate_diagnostics_xi_eta_heatmap.png"),
  height = 7.5)


# binned means with block-bootstrap CIs
binned_bootstrap <- function(d, y_columns, bin, size = 500, seed = 1) {
  blocks <- block_id(d$x_km, d$y_km, size)
  rows_by_block <- split(seq_len(nrow(d)), blocks)
  Y <- as.matrix(d[, y_columns])
  bin <- factor(bin)
  bin_means <- function(idx) {
    sums <- rowsum(Y[idx, , drop = FALSE], bin[idx], reorder = TRUE)
    counts <- as.vector(table(bin[idx]))
    out <- matrix(NA_real_, nlevels(bin), ncol(Y),
                  dimnames = list(levels(bin), y_columns))
    out[rownames(sums), ] <- sums / counts[counts > 0]
    out
  }
  estimate <- bin_means(seq_len(nrow(d)))
  set.seed(seed)
  boot <- simplify2array(lapply(seq_len(n_boot), function(b) {
    bin_means(unlist(rows_by_block[sample(length(rows_by_block),
                                          replace = TRUE)],
                     use.names = FALSE))
  }))
  # (built before the tibble, whose own `bin` column would mask the factor)
  n_bins <- nlevels(bin)
  counts <- as.vector(table(bin))
  tibble(bin = rep(levels(bin), ncol(Y)),
         series = rep(y_columns, each = n_bins),
         mean = as.vector(estimate),
         lo = as.vector(apply(boot, c(1, 2), quantile, 0.025, na.rm = TRUE)),
         hi = as.vector(apply(boot, c(1, 2), quantile, 0.975, na.rm = TRUE)),
         n = rep(counts, ncol(Y)))
}

# net use in fixed bins, cumulative net use in bins of net-use-years
nets_breaks <- c(-Inf, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, Inf)
nets_labels <- c("<0.05", "0.05-0.15", "0.15-0.25", "0.25-0.35",
                 "0.35-0.45", "0.45-0.55", "0.55-0.65", ">0.65")
cum_breaks <- c(-Inf, 0.25, 0.5, 1, 2, 3, 4, 6, 8, Inf)
cum_labels <- c("<0.25", "0.25-0.5", "0.5-1", "1-2", "2-3", "3-4", "4-6",
                "6-8", ">8")
m_breaks <- c(-Inf, -3, -2, -1, 0, 1, 2, 3, 4, 6, Inf)
m_labels <- c("<-3", "-3 to -2", "-2 to -1", "-1 to 0", "0 to 1", "1 to 2",
              "2 to 3", "3 to 4", "4 to 6", ">6")

binned_panels <- function(covariate, breaks, labels, series, series_labels,
                          transform, types_plot) {
  bind_rows(lapply(names(types_plot), function(group) {
    d <- filter(pixel_years, insecticide_type %in% types_plot[[group]]) %>%
      transform()
    bins <- cut(d[[covariate]], breaks, labels = labels)
    binned_bootstrap(d, series, bins, seed = 7) %>%
      mutate(group = group)
  })) %>%
    mutate(bin = factor(bin, levels = labels),
           series = factor(series_labels[series], levels = series_labels),
           group = factor(group, levels = names(types_plot)))
}

types_plot <- c(groups[c("Pyrethroids (pooled)", "Organophosphates (pooled)")],
                setNames(as.list(types), types))
types_plot <- types_plot[c("Pyrethroids (pooled)", "Deltamethrin",
                           "Permethrin", "Alpha-cypermethrin",
                           "Lambda-cyhalothrin", "Organophosphates (pooled)",
                           "DDT", "Bendiocarb")]
series_colours <- c("#0072B2", "#D55E00", "#009E73")

binned_plot <- function(binned, x_label, y_label, title, subtitle, file,
                        hline = TRUE) {
  dodge <- position_dodge(width = 0.5)
  p <- ggplot(binned, aes(x = bin, y = mean, colour = series,
                          group = series)) +
    geom_line(position = dodge, linewidth = 0.5) +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = dodge,
                    size = 0.25, linewidth = 0.5) +
    facet_wrap(~group, ncol = 4, scales = "free_y") +
    scale_colour_manual(values = series_colours, name = NULL) +
    labs(x = x_label, y = y_label, title = title,
         subtitle = paste(strwrap(subtitle, 160), collapse = "\n")) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          legend.position = "top", plot.subtitle = element_text(size = 7.5))
  if (hline) p <- p + geom_hline(yintercept = 0, colour = grey(0.5),
                                 linewidth = 0.3)
  ggsave(file, p, bg = "white", width = 12, height = 6.5, dpi = 200)
}

# (b) the annual selection against net use: the model's -log w (logit change
# per year from selection), and with eta added, the data's
binned_plot(
  binned_panels("nets", nets_breaks, nets_labels,
                c("minus_log_w", "corrected", "corrected_null"),
                c(minus_log_w = "dynamical: -log w",
                  corrected = "corrected: -log w + eta",
                  corrected_null = "null: -log w + eta (data simulated from m)"),
                function(d) mutate(d, minus_log_w = -log_w,
                                   corrected = -log_w + eta,
                                   corrected_null = -log_w + eta_null),
                types_plot),
  x_label = "net use in the year", y_label = "annual change in logit mortality",
  title = "Annual selection against net use, at the data pixel-years",
  subtitle = paste("Mean over pixel-years in each bin, 95% CI from a 500 km",
                   "spatial block bootstrap. The dynamical model's annual",
                   "selection includes all covariates (w linear in them); the",
                   "corrected series adds the fitted annual anomaly eta; the",
                   "null adds eta refitted to data simulated from m, i.e.",
                   "what the correction does when the model is right."),
  file = file.path(figure_dir, "covariate_diagnostics_binned_eta_nets.png"))

# ... and the cumulative change against cumulative net use
binned_plot(
  binned_panels("nets_cum", cum_breaks, cum_labels,
                c("minus_S", "corrected", "corrected_null"),
                c(minus_S = "dynamical: -S",
                  corrected = "corrected: -S + xi",
                  corrected_null = "null: -S + xi (data simulated from m)"),
                function(d) mutate(d, minus_S = -S, corrected = -S + xi,
                                   corrected_null = -S + xi_null),
                types_plot),
  x_label = "cumulative net use since 2000 (net-use years)",
  y_label = "cumulative change in logit mortality since 1995",
  title = "Cumulative selection against cumulative net use, at the data pixel-years",
  subtitle = paste("Mean over pixel-years in each bin, 95% CI from a 500 km",
                   "spatial block bootstrap. S = the dynamical model's",
                   "cumulative selection since 1995 (it also includes padded",
                   "1995-99 net use and the other covariates)."),
  file = file.path(figure_dir, "covariate_diagnostics_binned_xi_cumnets.png"))

# (c) the correction against m, observed and under the null
binned_plot(
  binned_panels("m", m_breaks, m_labels,
                c("smooth", "smooth_null", "z_minus_m"),
                c(smooth = "omega + xi (observed)",
                  smooth_null = "omega + xi (null: dynamical model true)",
                  z_minus_m = "z - m (raw residual)"),
                function(d) mutate(d, z_minus_m = z - m),
                types_plot),
  x_label = "dynamical prediction m (logit mortality)",
  y_label = "correction (logit)",
  title = "Correction against the dynamical prediction, at the data pixel-years",
  subtitle = paste("Mean over pixel-years in each bin, 95% CI from a 500 km",
                   "spatial block bootstrap. The null refits the fields",
                   "(fixed hyperparameters) to data simulated from m with",
                   "pixel-year and beta-binomial noise, so its trend is the",
                   "artefact of the empirical logit's bounds and shrinkage."),
  file = file.path(figure_dir, "covariate_diagnostics_binned_m.png"))

# (a) omega against the static covariates with the strongest within-country
# association (at least one type flagged or the largest |rho_within|)
top_static <- spearman %>%
  filter(field == "omega", block_km == 500, !group %in% names(groups)[10:11],
         !covariate %in% c("init", "m_mean", "latitude", "longitude")) %>%
  arrange(desc(abs(rho_within))) %>%
  slice_head(n = 6)
omega_binned <- bind_rows(lapply(seq_len(nrow(top_static)), function(i) {
  group <- top_static$group[i]
  covariate <- top_static$covariate[i]
  d <- filter(cells, insecticide_type == group) %>%
    group_by(country_id) %>%
    mutate(omega_within = omega - mean(omega)) %>%
    ungroup()
  breaks <- unique(quantile(d[[covariate]], seq(0, 1, 0.1)))
  bins <- cut(d[[covariate]], breaks, include.lowest = TRUE)
  levels(bins) <- sprintf("D%i", seq_len(nlevels(bins)))
  binned_bootstrap(d, c("omega", "omega_within"), bins, seed = 11) %>%
    mutate(group = sprintf("%s: %s (rho %.2f, within %.2f)", group,
                           covariate_labels[covariate], top_static$rho[i],
                           top_static$rho_within[i]))
})) %>%
  mutate(series = recode(series, omega = "omega",
                         omega_within = "omega minus country mean"),
         bin = factor(bin, levels = sprintf("D%i", 1:10)))
binned_plot(
  omega_binned, x_label = "covariate decile", y_label = "omega (logit)",
  title = "Static correction against its most associated static covariates",
  subtitle = paste("Mean over the type's data cells per decile of the",
                   "covariate, 95% CI from a 500 km spatial block bootstrap.",
                   "Panels: the six type x covariate pairs with the largest",
                   "within-country |rho|."),
  file = file.path(figure_dir, "covariate_diagnostics_binned_omega.png"))


# 4. temporal pattern: mean xi by year and region, over each type's data cells
# (a fixed set of locations, so changes in where the data are do not move it)
region_colours <- c("Eastern Africa" = "#0072B2",
                    "Middle Africa" = "#E69F00",
                    "Northern Africa" = "#999999",
                    "Southern Africa" = "#CC79A7",
                    "Western Africa" = "#009E73")
trajectories <- bind_rows(lapply(types, function(type) {
  f <- fields[[type]]
  as_tibble(f$xi_cells) %>%
    mutate(region = f$cell_coords$region) %>%
    pivot_longer(-region, names_to = "year", values_to = "xi") %>%
    mutate(year = as.numeric(year), insecticide_type = type)
})) %>%
  group_by(insecticide_type, region, year) %>%
  summarise(xi = mean(xi), n_cells = n(), .groups = "drop")
# data-weighted: the year's share of the type's data in the region, to see
# where the trajectory is informed
data_share <- pixel_years %>%
  count(insecticide_type, region, year, name = "n_pixel_years")
trajectories <- left_join(trajectories, data_share,
                          by = c("insecticide_type", "region", "year")) %>%
  mutate(n_pixel_years = replace_na(n_pixel_years, 0),
         insecticide_type = factor(insecticide_type,
                                   levels = group_levels[1:9]))
write.csv(trajectories,
          file.path(cache_dir, "xi_trajectories_by_region.csv"),
          row.names = FALSE)
ggplot(trajectories, aes(x = year, y = xi, colour = region)) +
  geom_hline(yintercept = 0, colour = grey(0.5), linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  geom_point(data = filter(trajectories, n_pixel_years > 0),
             aes(size = n_pixel_years), alpha = 0.6) +
  scale_size_area(max_size = 2.5, name = "pixel-years\nwith data") +
  scale_colour_manual(values = region_colours, name = NULL) +
  facet_wrap(~insecticide_type, ncol = 3) +
  labs(x = NULL, y = "mean xi (logit)",
       title = "Accumulated anomaly xi by year and region",
       subtitle = paste("Posterior mode of xi averaged over each type's data",
                        "cells in the region (fixed locations), 1995-T;",
                        "negative = less susceptible than the dynamical",
                        "model. Points mark years with data in the region.")) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "right", plot.subtitle = element_text(size = 7.5))
ggsave(file.path(figure_dir, "covariate_diagnostics_xi_by_year_region.png"),
       bg = "white", width = 11, height = 8, dpi = 200)

report("done")
