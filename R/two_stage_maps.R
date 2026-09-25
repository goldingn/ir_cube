# Maps of the two-stage model (#21): the stage-A correction fitted per
# insecticide type to ALL the bioassay data, on top of the full dynamical fit
# (R/fit_model.R), projected onto the full prediction grid for the years the
# dynamical-model map figures show (R/fig_ir_maps.R).
#
#   Rscript R/two_stage_maps.R
#
# Run with OpenBLAS for CHOLMOD's supernodal factorisation (reference BLAS is
# ~10x slower), e.g.
#   LD_PRELOAD=.../libopenblas.so.0 OPENBLAS_NUM_THREADS=4 nice -n 10 Rscript ...
#
# Steps:
#
#   1. load the full dynamical fit (temporary/fitted_model.RData) and check its
#      data and covariates are exactly the ones the current scripts build;
#      recompute its 2000 paired posterior draws of logit p at every assay
#      (R/dynamical_predictions.R), whose mean is m_ref;
#   2. recompute the dynamical model's posterior mean mortality at a sample of
#      mask cells and check it against the maps predict.R saved in
#      outputs/ir_maps. Those maps have scrambled country initial conditions
#      (a greta subassignment in predict.R, see map_logit_init()), so the check
#      emulates that to confirm the draws and covariates are the ones behind
#      them, and the maps here use the correct initial conditions;
#   3. per type, fit stage A (omega_xi_u, default meshes, per-type rho,
#      t0 = 1995, T = the type's last data year) to all its assays;
#   4. per type, on every mask cell and map year, compute
#        - the correction's posterior mean (omega + xi, logit scale): the
#          latent mode projected to the cells. u is left out: its mean is 0 away
#          from sampled pixel-years, and a map of the smooth fields is the
#          point; so all the mortality maps here are of the smooth part,
#          m + omega + xi, not the pixel-year target m + omega + xi + u
#        - its posterior SD, from 200 joint latent draws, conditional on m_ref
#          (the second stage's own uncertainty)
#        - the two-stage posterior mean mortality, E[ilogit(m + omega + xi)]
#          over 200 dynamical draws, each paired with a latent draw shifted by
#          the cut-posterior formula (as predict_correction() does), so it is
#          the cut posterior's mean of mortality rather than the plug-in
#          ilogit(logit(E[mortality]) + E[correction]), which is also saved
#        - the dynamical posterior mean mortality from the same 200 draws, and
#          the difference two-stage - dynamical in percentage points;
#      and write them as rasters to outputs/two_stage/maps/<type>/;
#   5. figures in figures/two_stage/, in the layout of R/fig_ir_maps.R.
#
# Beyond T (2024 for most types, 2023 for Fenitrothion and Malathion) xi is an
# AR(1) forecast: its mean is xi_T + eta_T phi (1 - phi^h) / (1 - phi), which
# plateaus at a rate set by phi, and its SD grows with the horizon h.
#
# Steps 3-4 skip a type whose rasters already exist unless `overwrite` is TRUE,
# so the figures can be redrawn without refitting.

overwrite <- FALSE

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

output_dir <- "outputs/two_stage/maps"
figure_dir <- "figures/two_stage"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

# the panel years of R/fig_ir_maps.R
map_years <- c(2000, 2005, 2010, 2015, 2020, 2025, 2030)
end_year <- max(map_years)
years_all <- baseline_year:end_year

# xi(., t0) = 0 at the dynamical model's start year, as in the CV runner
t0 <- baseline_year

# posterior draws for the maps. 200 keeps a chunk of cells x draws x years
# small; the Monte Carlo SE of a mean mortality is then ~0.3 percentage points
# at a posterior SD of 5
n_map_draws <- 200
chunk_size <- 25000

# the clamp of R/run_two_stage_folds.R, so m_ref here is on its footing
clamp <- 1e-12
safe_logit <- function(p) qlogis(pmin(pmax(p, clamp), 1 - clamp))
logit_max <- qlogis(1 - clamp)

peak_memory_gb <- function() {
  status <- readLines("/proc/self/status")
  as.numeric(gsub("\\D", "", grep("^VmHWM", status, value = TRUE))) / 1024 ^ 2
}
time_start <- Sys.time()

quantities <- c("two_stage_mortality", "two_stage_mortality_plugin",
                "dynamical_mortality", "difference_pp", "correction_mean",
                "correction_sd")
raster_file <- function(type, quantity) {
  file.path(output_dir, type, sprintf("%s.tif", quantity))
}
types_done <- types[vapply(types, function(type) {
  all(file.exists(raster_file(type, quantities)))
}, logical(1))]
types_to_fit <- if (overwrite) types else setdiff(types, types_done)


if (length(types_to_fit) > 0) {

  # 1. the full dynamical fit ----------------------------------------------

  # save.image() of R/fit_model.R, November 2025. Its data must be exactly what
  # the current scripts build (the data preparation has changed since), or its
  # draws cannot be paired with the current assays
  fit_env <- new.env()
  load("temporary/fitted_model.RData", envir = fit_env)
  stopifnot(
    isTRUE(all.equal(fit_env$df, df)),
    identical(fit_env$types, types),
    identical(fit_env$classes, classes),
    identical(fit_env$countries, countries),
    identical(fit_env$regions, regions),
    identical(fit_env$unique_cells, unique_cells),
    identical(fit_env$classes_index, classes_index),
    isTRUE(all.equal(fit_env$x_cell_years, x_cell_years)),
    isTRUE(all.equal(fit_env$cell_years_index, cell_years_index))
  )
  report("full fit: data, indexing and covariates match the current build (%i assays)",
         nrow(df))

  # dynamical_predictions() only needs `draws` from a fold
  fold <- list(draws = fit_env$draws)
  rm(fit_env)
  invisible(gc())

  # 2000 of the 16000 draws (8 chains x 2000), by the CV folds' thinning rule
  draw_index <- paired_draw_index(fold)
  draws_matrix <- as.matrix(fold$draws)[draw_index, , drop = FALSE]
  logit_init_mean <- logit_init_mean_draws(fold, draw_index)
  parameters <- dynamical_parameter_draws(fold, df = df,
                                          classes_index = classes_index,
                                          types = types,
                                          draw_index = draw_index,
                                          logit_init_mean = logit_init_mean)

  # draws at every assay, with the country looked up from the cell as the fit
  # did (see R/run_two_stage_folds.R on border cells)
  time_train <- system.time(
    p_train <- dynamical_predictions(fold, select(df, -country_id), df,
                                     x_cell_years, cell_years_index,
                                     classes_index, types,
                                     draw_index = draw_index)
  )
  logit_train <- safe_logit(p_train)
  rm(p_train)
  m_ref <- colMeans(logit_train)
  report("dynamical draws at %i assays x %i draws in %.0f s",
         ncol(logit_train), nrow(logit_train), time_train[["elapsed"]])

  # the map draws are an even subset of the 2000
  map_draws <- round(seq(1, length(draw_index), length.out = n_map_draws))
  logit_train_map <- logit_train[map_draws, , drop = FALSE]
  rm(logit_train)
  invisible(gc())


  # the prediction grid -------------------------------------------------------

  cells <- terra::cells(mask)
  n_cells <- length(cells)
  covariates <- map_covariates(cells, baseline_year, end_year)

  # the grid's covariates must be the model's at the data cells, 1995-2024
  data_cells <- match(unique_cells, cells)
  stopifnot(!anyNA(data_cells))
  x_check <- x_cell_years[order(cell_years_index$cell_id,
                                cell_years_index$year_id), ]
  n_fit_years <- max(cell_years_index$year_id)
  x_grid <- do.call(rbind, lapply(seq_along(unique_cells), function(i) {
    cbind(covariates$time_varying[data_cells[i], seq_len(n_fit_years), ],
          covariates$flat[rep(data_cells[i], n_fit_years), , drop = FALSE])
  }))
  stopifnot(max(abs(x_grid - x_check)) < 1e-12)
  rm(x_grid, x_check)

  # initial conditions: the cell's country from the country raster and the
  # region from the UNSD lookup, as in predict.R
  lookup <- country_region_lookup()
  country_raster <- rast("data/clean/country_raster.tif")
  cell_country <- as.character(terra::extract(country_raster,
                                              cells)$country_name)
  logit_init_all <- map_logit_init(draws_matrix, logit_init_mean, types,
                                   countries, regions, lookup)
  cell_country_index <- match(cell_country, dimnames(logit_init_all)[[2]])
  report("%i of %i grid cells are in countries outside the UNSD lookup (NA)",
         sum(is.na(cell_country_index)), n_cells)

  # for observed countries, the same initial state as the fit's (so the UNSD
  # regions agree with the regions of the data)
  stopifnot(isTRUE(all.equal(
    logit_init_all[, countries, , drop = FALSE],
    parameters$logit_init, check.attributes = FALSE)))

  xy <- terra::xyFromCell(mask, cells)
  coords_cells <- project_km(xy[, 1], xy[, 2])


  # 2. check against the saved dynamical maps ---------------------------------

  # predict.R's maps are posterior means over 500 random draws of the same fit;
  # the recomputation uses the 2000 thinned draws, so the two can agree only up
  # to Monte Carlo error, judged against the posterior SD.
  #
  # They do not agree as predict.R intended: its greta subassignment scrambles
  # the country and region initial conditions (see map_logit_init()), so its
  # saved maps have each country's initial state drawn from another country
  # and type. Recomputed with the correct initial conditions, the maps agree
  # with the observed mortality far better than the saved ones do (correlation
  # at the assays about 0.4-0.67 per type, against -0.13-0.26 for the saved
  # maps). So the check reproduces the saved maps with that scrambling
  # emulated, which confirms that the draws, covariates and recursion here are
  # the ones behind them; the maps below use the correct initial conditions
  set.seed(21)
  check_cells <- sort(sample(which(!is.na(cell_country_index)), 5000))
  check_years <- c(2000, 2010, 2020, 2030)
  logit_init_predict <- map_logit_init(draws_matrix, logit_init_mean, types,
                                       countries, regions, lookup,
                                       emulate_predict_fill = TRUE)
  check_rows <- list()
  for (k in seq_along(types)) {
    for (initial in c("correct", "predict_fill")) {
      init <- if (initial == "correct") logit_init_all else logit_init_predict
      dyn <- dynamical_logit_chunk(
        effect = parameters$effect_type[, , k],
        logit_init = t(init[, cell_country_index[check_cells], k]),
        time_varying = covariates$time_varying[check_cells, , , drop = FALSE],
        flat = covariates$flat[check_cells, , drop = FALSE],
        years = years_all, years_keep = check_years)
      for (y in check_years) {
        p <- plogis(dyn[[as.character(y)]])
        saved <- rast(sprintf("outputs/ir_maps/%s/ir_%i_susceptibility.tif",
                              types[k], y))[cells[check_cells]][, 1]
        check_rows[[length(check_rows) + 1]] <- tibble(
          initial = initial, insecticide_type = types[k], year = y,
          cell = cells[check_cells], recomputed = rowMeans(p), saved = saved,
          posterior_sd = row_sds(p))
      }
    }
  }
  rm(logit_init_predict)
  map_check <- bind_rows(check_rows) %>%
    mutate(difference = recomputed - saved,
           # SE of the difference of two Monte Carlo means (2000 and 500
           # draws), floored so that cells where p is ~1 in every draw do not
           # dominate
           z = difference / (pmax(posterior_sd, 1e-3) *
                               sqrt(1 / 2000 + 1 / 500)))
  map_check_summary <- map_check %>%
    group_by(initial, insecticide_type, year) %>%
    summarise(n = n(),
              n_na_saved = sum(is.na(saved)),
              correlation = cor(recomputed, saved),
              max_abs_difference = max(abs(difference), na.rm = TRUE),
              mean_difference = mean(difference, na.rm = TRUE),
              mean_z = mean(z, na.rm = TRUE),
              sd_z = sd(z, na.rm = TRUE),
              share_abs_z_above_4 = mean(abs(z) > 4, na.rm = TRUE),
              .groups = "drop")
  write.csv(map_check_summary, file.path(output_dir, "dynamical_map_check.csv"),
            row.names = FALSE)
  for (initial in c("correct", "predict_fill")) {
    check_i <- map_check[map_check$initial == initial, ]
    report(paste("recomputed (%s initial conditions) vs saved dynamical maps:",
                 "max |diff| %.4f, mean z %.3f, sd z %.2f, |z| > 4 in %.2f%%"),
           initial, max(abs(check_i$difference)), mean(check_i$z),
           sd(check_i$z), 100 * mean(abs(check_i$z) > 4))
  }
  # a mismatch in data, covariates or parameters would show as a bias or a
  # spread of z far beyond Monte Carlo error, not as a few large values. The
  # saved maps share one set of 500 draws across all cells, so their Monte
  # Carlo errors are correlated between cells and the mean z need not be 0
  # (it is about 0.2)
  check_fill <- map_check[map_check$initial == "predict_fill", ]
  stopifnot(all(map_check_summary$n_na_saved == 0),
            abs(mean(check_fill$z)) < 0.5,
            sd(check_fill$z) < 1.5,
            mean(abs(check_fill$z) > 4) < 0.005)
  rm(map_check, check_rows, check_fill, dyn)

  # only the map draws are needed from here on
  effect_map <- parameters$effect_type[map_draws, , , drop = FALSE]
  logit_init_map <- logit_init_all[map_draws, , , drop = FALSE]
  rm(parameters, logit_init_all, draws_matrix)
  invisible(gc())


  # 3-4. fit and map each type ------------------------------------------------

  rho_table <- read.csv("outputs/bioassay_rho_hierarchical.csv")
  rho_for_type <- setNames(rho_table$rho, rho_table$insecticide_type)
  stopifnot(all(types %in% names(rho_for_type)))

  # a raster with one layer per map year, from a cells x years matrix
  to_raster <- function(values) {
    full <- matrix(NA_real_, ncell(mask), ncol(values))
    full[cells, ] <- values
    r <- rast(mask, nlyrs = ncol(values))
    values(r) <- full
    names(r) <- map_years
    r
  }

  chunks <- split(seq_len(n_cells), ceiling(seq_len(n_cells) / chunk_size))

  for (type in types_to_fit) {

    k <- match(type, types)
    rho <- rho_for_type[[type]]
    rows_k <- which(df$type_id == k)

    train_k <- tibble(
      lon = df$longitude[rows_k],
      lat = df$latitude[rows_k],
      year = df$year_start[rows_k],
      cell = df$cell[rows_k],
      died = df$died[rows_k],
      mosquito_number = df$mosquito_number[rows_k],
      m = m_ref[rows_k],
      rho = rho
    )
    stage_a <- empirical_logit(train_k$died, train_k$mosquito_number, rho)
    train_k$z <- stage_a$z
    train_k$v <- stage_a$v
    T_k <- max(train_k$year)

    # the default meshes of the CV runner
    coords <- coords_km(train_k)
    mesh <- suppressMessages(build_correction_mesh(coords, verbose = FALSE))
    mesh_xi <- suppressMessages(build_correction_mesh(coords, max_nodes = 600,
                                                      verbose = FALSE))

    set.seed(2026 + k)
    time_fit <- system.time(
      fit <- fit_correction(train_k, variant = "omega_xi_u", t0 = t0, T = T_k,
                            mesh = mesh, mesh_xi = mesh_xi)
    )
    max_gradient <- max(abs(fit$obj$gr(fit$opt$par)))
    report(paste("%-18s n=%5i T=%i nodes=%i/%i conv=%i |grad|=%.1e",
                 "range_omega=%.0f sigma_omega=%.2f range_eta=%.0f",
                 "sigma_eta=%.3f phi=%.2f tau=%.3f fit %.0f s"),
           type, nrow(train_k), T_k, mesh$n, mesh_xi$n, fit$opt$convergence,
           max_gradient, fit$hyper$range_omega, fit$hyper$sigma_omega,
           fit$hyper$range_eta, fit$hyper$sigma_eta, fit$hyper$phi,
           fit$hyper$tau, time_fit[["elapsed"]])

    # latent draws at the mesh nodes: joint deviations from the mode (so omega
    # and xi keep their posterior correlation), without and with the
    # cut-posterior shift for each map draw's dynamical offset
    set.seed(3026 + k)
    time_predict <- system.time({
      n_latent <- length(fit$mode)
      theta_cond <- sample_latent_deviation(fit$H_chol, n_latent,
                                            n_map_draws) + fit$mode
      theta_full <- theta_cond +
        correction_mode_shift(fit, t(logit_train_map[, rows_k, drop = FALSE]))

      # AR(1) innovations beyond T, shared by both sets of draws so they differ
      # only by the shift
      max_horizon <- max(0, end_year - T_k)
      innovations <- NULL
      if (max_horizon > 0) {
        Q_eta <- matern_precision_r(fit$fem_xi, fit$hyper$kappa_eta,
                                    fit$hyper$sigma_eta)
        Q_eta_chol <- Matrix::Cholesky(Matrix::forceSymmetric(Q_eta),
                                       perm = TRUE, LDL = FALSE, super = TRUE)
        innovations <- lapply(seq_len(max_horizon), function(h) {
          sample_latent_deviation(Q_eta_chol, fit$mesh_xi$n, n_map_draws)
        })
      }
      nodes_mean <- correction_node_fields(fit, fit$mode, map_years)
      nodes_cond <- correction_node_fields(fit, theta_cond, map_years,
                                           innovations)
      nodes_full <- correction_node_fields(fit, theta_full, map_years,
                                           innovations)
      rm(theta_cond, theta_full, innovations)

      results <- lapply(setNames(quantities, quantities), function(q) {
        matrix(NA_real_, n_cells, length(map_years))
      })

      for (chunk in chunks) {
        ok <- chunk[!is.na(cell_country_index[chunk])]
        if (length(ok) == 0) next
        dyn <- dynamical_logit_chunk(
          effect = effect_map[, , k],
          logit_init = t(logit_init_map[, cell_country_index[ok], k]),
          time_varying = covariates$time_varying[ok, , , drop = FALSE],
          flat = covariates$flat[ok, , drop = FALSE],
          years = years_all, years_keep = map_years)
        A_omega <- mesh_basis(fit$mesh, coords_cells[ok, , drop = FALSE])
        A_xi <- mesh_basis(fit$mesh_xi, coords_cells[ok, , drop = FALSE])
        project <- function(nodes, y) {
          as.matrix(A_omega %*% nodes$omega +
                      A_xi %*% nodes$xi[[as.character(y)]])
        }
        for (j in seq_along(map_years)) {
          y <- map_years[j]
          m <- pmin(pmax(dyn[[as.character(y)]], -logit_max), logit_max)
          p_dyn <- rowMeans(plogis(m))
          correction_mean <- as.vector(project(nodes_mean, y))
          p_two_stage <- rowMeans(plogis(m + project(nodes_full, y)))
          results$dynamical_mortality[ok, j] <- p_dyn
          results$two_stage_mortality[ok, j] <- p_two_stage
          # plug-in: ilogit(logit(E[p_dyn]) + E[correction]), kept for
          # comparison with the cut posterior's mean
          results$two_stage_mortality_plugin[ok, j] <-
            plogis(qlogis(p_dyn) + correction_mean)
          results$difference_pp[ok, j] <- 100 * (p_two_stage - p_dyn)
          results$correction_mean[ok, j] <- correction_mean
          results$correction_sd[ok, j] <- row_sds(project(nodes_cond, y))
        }
        rm(dyn, m)
      }
    })

    dir.create(file.path(output_dir, type), showWarnings = FALSE,
               recursive = TRUE)
    for (q in quantities) {
      writeRaster(to_raster(results[[q]]), raster_file(type, q),
                  overwrite = TRUE, datatype = "FLT4S",
                  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3"))
    }

    # correlation of the fitted smooth correction (omega + xi at the mode, no
    # u) with the covariates and with m at the assays: structure the dynamical
    # model's covariates could have explained points to misspecification there
    smooth_blocks <- c(fit$blocks$w_omega, fit$blocks$x)
    smooth_train <- as.vector(fit$A_latent[, smooth_blocks] %*%
                                fit$mode[smooth_blocks])
    x_rows <- match(paste(df$cell_id[rows_k], df$year_id[rows_k]),
                    paste(cell_years_index$cell_id, cell_years_index$year_id))
    x_train <- cbind(x_cell_years[x_rows, , drop = FALSE],
                     m_ref = m_ref[rows_k], year = train_k$year)
    covariate_correlation <- tibble(
      insecticide_type = type,
      covariate = colnames(x_train),
      correlation = as.vector(cor(x_train, smooth_train)),
      # pixel-years, not assays, are the unit of the fields
      n_assays = nrow(train_k)
    )

    # a light copy of the fit (no TMB object or Hessian) for later use
    saveRDS(list(hyper = fit$hyper, mode = fit$mode, blocks = fit$blocks,
                 mesh = fit$mesh, mesh_xi = fit$mesh_xi, t0 = fit$t0,
                 T = fit$T, n_years = fit$n_years,
                 pixel_years = fit$pixel_years, opt = fit$opt,
                 covariate_correlation = covariate_correlation),
            file.path(output_dir, type, "fit.rds"))

    hyper_row <- tibble(
      insecticide_type = type, rho = rho, n_assays = nrow(train_k),
      n_pixel_years = nrow(fit$pixel_years), T = T_k,
      mesh_nodes = mesh$n, mesh_xi_nodes = mesh_xi$n,
      range_omega_km = fit$hyper$range_omega,
      sigma_omega = fit$hyper$sigma_omega,
      range_eta_km = fit$hyper$range_eta,
      sigma_eta = fit$hyper$sigma_eta,
      phi = fit$hyper$phi, persistence_years = fit$hyper$persistence,
      tau = fit$hyper$tau, convergence = fit$opt$convergence,
      max_gradient = max_gradient,
      time_fit_s = time_fit[["elapsed"]],
      time_map_s = time_predict[["elapsed"]]
    )
    write.csv(hyper_row, file.path(output_dir, type, "hyperparameters.csv"),
              row.names = FALSE)
    write.csv(covariate_correlation,
              file.path(output_dir, type, "covariate_correlation.csv"),
              row.names = FALSE)
    report("%-18s mapped in %.0f s; peak memory %.1f GB", type,
           time_predict[["elapsed"]], peak_memory_gb())

    rm(fit, results, nodes_mean, nodes_cond, nodes_full)
    invisible(gc())
  }
}

# per-type tables, combined
hyperparameters <- bind_rows(lapply(types, function(type) {
  read.csv(file.path(output_dir, type, "hyperparameters.csv"))
}))
write.csv(hyperparameters, file.path(output_dir, "hyperparameters.csv"),
          row.names = FALSE)
covariate_correlation <- bind_rows(lapply(types, function(type) {
  read.csv(file.path(output_dir, type, "covariate_correlation.csv"))
}))
write.csv(covariate_correlation,
          file.path(output_dir, "covariate_correlation.csv"),
          row.names = FALSE)
report("tables written; %.0f min so far",
       as.numeric(difftime(Sys.time(), time_start, units = "mins")))


# 5. figures -------------------------------------------------------------------

# the look of R/fig_ir_maps.R: grey Africa background, thin grey borders, masked
# to the limits of Pf transmission and water bodies, one panel per year in two
# rows with the legend in the eighth slot

borders <- readRDS("data/clean/gadm_polys.RDS")
pf_water_mask <- rast("data/clean/pfpr_water_mask.tif")

africa_bg <- geom_sf(data = borders,
                     linewidth = 0,
                     fill = grey(0.75))
border_col <- grey(0.4)
country_borders <- geom_sf(data = borders,
                           col = border_col,
                           linewidth = 0.1,
                           fill = "transparent")

# the per-insecticide colours of R/fig_ir_maps.R
insecticides_plot <- c("Alpha-cypermethrin",
                       "Deltamethrin",
                       "Lambda-cyhalothrin",
                       "Permethrin",
                       "Fenitrothion",
                       "Malathion",
                       "Pirimiphos-methyl",
                       "DDT",
                       "Bendiocarb")
insecticides_col <- setNames(rev(scales::hue_pal()(length(insecticides_plot))),
                             insecticides_plot)

read_map <- function(type, quantity) {
  r <- rast(raster_file(type, quantity))
  names(r) <- map_years
  terra::mask(r, pf_water_mask)
}

# diverging scale centred at 0. Positive = more susceptible (higher mortality)
# under the two-stage model than under the dynamical model
diverging_scale <- function(name, limit, labels = waiver()) {
  scale_fill_gradient2(
    name = name,
    low = "#b2182b",
    mid = "white",
    high = "#2166ac",
    midpoint = 0,
    limits = c(-limit, limit),
    oob = scales::squish,
    labels = labels,
    na.value = "transparent",
    guide = guide_colorbar(frame.colour = border_col,
                           frame.linewidth = 0.1))
}

year_panels <- function(raster, fill_scale, title, subtitle, file) {
  years_list <- lapply(seq_len(nlyr(raster)), function(i) {
    ggplot() +
      africa_bg +
      geom_spatraster(data = raster[[i]]) +
      country_borders +
      fill_scale +
      facet_wrap(~lyr, nrow = 1, ncol = 1) +
      theme_ir_maps() +
      theme(
        plot.margin = unit(rep(0, 4), "cm"),
        legend.text.position = "left",
        legend.ticks = element_blank()
      )
  })
  patchwork::wrap_plots(c(years_list, list(patchwork::guide_area()))) +
    patchwork::plot_layout(guides = "collect", nrow = 2) +
    patchwork::plot_annotation(title = title, subtitle = subtitle)
  ggsave(file, bg = "white", width = 13, height = 8, scale = 0.8, dpi = 300)
}

# common limits across types, so the correction and difference maps can be
# compared between insecticides: the 99.5th percentile of the absolute value
# over all types, years and (masked) cells, rounded
pooled_quantile <- function(quantity, prob = 0.995) {
  values <- unlist(lapply(types, function(type) {
    v <- values(read_map(type, quantity), mat = FALSE)
    abs(v[!is.na(v)])
  }))
  quantile(values, prob, names = FALSE)
}
# the difference made by the correction's mean alone:
# ilogit(logit(E[p_dyn]) + E[correction]) - E[p_dyn], in percentage points
read_plugin_difference <- function(type) {
  100 * (read_map(type, "two_stage_mortality_plugin") -
           read_map(type, "dynamical_mortality"))
}
correction_limit <- ceiling(pooled_quantile("correction_mean") * 4) / 4
difference_limit <- ceiling(pooled_quantile("difference_pp") / 5) * 5
sd_limit <- ceiling(pooled_quantile("correction_sd", 1) * 10) / 10

forecast_note <- function(type) {
  T_k <- hyperparameters$T[hyperparameters$insecticide_type == type]
  sprintf("data to %i; later years are the AR(1) forecast of xi", T_k)
}

for (type in insecticides_plot) {

  year_panels(
    read_map(type, "two_stage_mortality"),
    scale_fill_gradient(
      labels = scales::percent,
      name = "Susceptibility",
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      high = insecticides_col[[type]],
      low = "white",
      na.value = "transparent",
      guide = guide_colorbar(frame.colour = border_col,
                             frame.linewidth = 0.1)),
    title = sprintf("%s: two-stage model", type),
    subtitle = paste("Susceptibility of An. gambiae (s.l./s.s.) in WHO",
                     "bioassays; posterior mean of the smooth part (no",
                     "pixel-year effect u);", forecast_note(type)),
    file = file.path(figure_dir, sprintf("%s_two_stage_ir_map.png", type)))

  # the plug-in counterpart, which shows only the correction's mean: the
  # posterior mean above also carries the correction's variance, which pulls
  # mortality towards 50% wherever that variance is large, i.e. far from data
  # and increasingly so in later years (xi is close to a random walk)
  year_panels(
    read_map(type, "two_stage_mortality_plugin"),
    scale_fill_gradient(
      labels = scales::percent,
      name = "Susceptibility",
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      high = insecticides_col[[type]],
      low = "white",
      na.value = "transparent",
      guide = guide_colorbar(frame.colour = border_col,
                             frame.linewidth = 0.1)),
    title = sprintf("%s: two-stage model (plug-in)", type),
    subtitle = paste("ilogit(logit(dynamical posterior mean) + posterior mean",
                     "correction), no pixel-year effect u;",
                     forecast_note(type)),
    file = file.path(figure_dir,
                     sprintf("%s_two_stage_plugin_ir_map.png", type)))

  year_panels(
    read_map(type, "correction_mean"),
    diverging_scale("Correction<br>(logit)", correction_limit),
    title = sprintf("%s: second-stage correction", type),
    subtitle = paste("Posterior mean of omega + xi on the logit scale",
                     "(+ = more susceptible than the dynamical model);",
                     forecast_note(type)),
    file = file.path(figure_dir, sprintf("%s_correction_map.png", type)))

  year_panels(
    read_map(type, "difference_pp"),
    diverging_scale("Difference<br>(% points)", difference_limit),
    title = sprintf("%s: two-stage minus dynamical", type),
    subtitle = paste("Difference in posterior mean bioassay mortality",
                     "(including the shrinkage towards 50% from the",
                     "correction's variance), percentage points;",
                     forecast_note(type)),
    file = file.path(figure_dir, sprintf("%s_difference_map.png", type)))

  year_panels(
    read_plugin_difference(type),
    diverging_scale("Difference<br>(% points)", difference_limit),
    title = sprintf("%s: two-stage (plug-in) minus dynamical", type),
    subtitle = paste("Difference made by the posterior mean correction alone,",
                     "percentage points;", forecast_note(type)),
    file = file.path(figure_dir,
                     sprintf("%s_difference_plugin_map.png", type)))

  year_panels(
    read_map(type, "correction_sd"),
    scale_fill_gradient(
      name = "Correction SD<br>(logit)",
      limits = c(0, sd_limit),
      low = "white",
      high = grey(0.1),
      na.value = "transparent",
      guide = guide_colorbar(frame.colour = border_col,
                             frame.linewidth = 0.1)),
    title = sprintf("%s: second-stage correction SD", type),
    subtitle = paste("Posterior SD of omega + xi (logit scale), given the",
                     "dynamical posterior mean;", forecast_note(type)),
    file = file.path(figure_dir, sprintf("%s_correction_sd_map.png", type)))
}

# all types side by side, for one year, to compare the corrections' structure
# between insecticides
compare_year <- 2020
all_types <- rast(lapply(insecticides_plot, function(type) {
  r <- read_map(type, "correction_mean")[[as.character(compare_year)]]
  names(r) <- type
  r
}))
ggplot() +
  africa_bg +
  geom_spatraster(data = all_types) +
  country_borders +
  facet_wrap(~lyr, ncol = 3) +
  diverging_scale("Correction<br>(logit)", correction_limit) +
  theme_ir_maps() +
  theme(plot.margin = unit(rep(0, 4), "cm"),
        legend.ticks = element_blank()) +
  labs(title = sprintf("Second-stage correction in %i", compare_year),
       subtitle = "Posterior mean of omega + xi, logit scale")
ggsave(file.path(figure_dir,
                 sprintf("correction_all_types_%i.png", compare_year)),
       bg = "white", width = 10, height = 10, scale = 0.8, dpi = 300)

# fitted hyperparameters per type
hyperparameters %>%
  select(insecticide_type,
         `omega range (km)` = range_omega_km,
         `omega SD` = sigma_omega,
         `eta range (km)` = range_eta_km,
         `eta SD` = sigma_eta,
         `phi` = phi,
         `tau` = tau) %>%
  pivot_longer(-insecticide_type) %>%
  mutate(name = factor(name, levels = unique(name)),
         insecticide_type = factor(insecticide_type,
                                   levels = rev(insecticides_plot))) %>%
  ggplot(aes(x = value, y = insecticide_type,
             colour = insecticide_type)) +
  geom_point(size = 2) +
  facet_wrap(~name, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = insecticides_col, guide = "none") +
  labs(x = NULL, y = NULL,
       title = "Second-stage hyperparameters, fitted to all data") +
  theme_minimal()
ggsave(file.path(figure_dir, "map_hyperparameters.png"),
       bg = "white", width = 12, height = 3.5, dpi = 200)

report("done in %.0f min", as.numeric(difftime(Sys.time(), time_start,
                                               units = "mins")))
