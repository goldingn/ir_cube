# Posterior predicted mortality to the LLIN pyrethroids over West Africa, on a
# coarse grid, for the diagnostic maps of the West's pyrethroid plateau (#47;
# R/west_figures.R draws them), after checking the prediction path at the
# bioassays.
#
#   USE_CHAINS="V5=1,2" Rscript R/west_predictions.R [<label>=<fitted_model.RData> ...]
#
# The fits are ref_f0, V3f, V4_class and V5 (`fits` below), or those given.
# For each, from the n_check (500) posterior draws that R/species_misfit.R
# used (even_draws(), R/species_fit_helpers.R; with USE_CHAINS, as there):
#
# 1. The check. Predictions at n_rows bioassays (half of them LLIN-pyrethroid
#    bioassays in the West), by the map path: the cells' covariates, the
#    arabiensis share r(x) (prediction_share()), the standardised kdr
#    (prediction_kdr()) and the centred basis of the latent smooths
#    (prediction_basis()) from two_stage_cells() and two_stage_chunk()
#    (R/two_stage_predictions.R), and dynamical_logit_cells(), at the fit's
#    country of the cell (data_cell_country()). Against the posterior mean
#    predicted mortality of the same bioassays in
#    outputs/species_runs/misfit/<label>_bioassays.csv, from the same draws:
#    they must agree to rounding. The bioassays are of the complex (species
#    not identified), whose share in the misfit is r(x) at the cell, as on
#    the map. Also reported: the difference at the n_draws map draws (a
#    subset; Monte Carlo error), and at the country of the cell in
#    data/clean/country_raster.tif, which the map uses.
#    Writes outputs/species_runs/west/check_<label>.csv.
#
# 2. The grid. Every 4th cell of the model mask in each direction (about
#    0.17 degrees), in the West (analysis_region(): UNSD Western Africa,
#    data/clean/region_raster.tif) and within data/clean/pfpr_water_mask.tif.
#    At each, the posterior mean and sd, over n_draws (160) of the draws
#    above, of mortality averaged over the years of each window (2010-2015,
#    2019-2025) and over the three LLIN pyrethroids, weighted as the regional
#    trend figure (R/species_compare.R) weights the insecticides within the
#    West: by their mosquitoes tested in the West's modelled bioassays, all
#    years; and each pyrethroid's own posterior mean. With the species model,
#    of the whole complex at r(x). Countries take their fitted initial states
#    (all West countries with bioassays; others from the hierarchical prior,
#    map_logit_init(), seed 1, as R/predict.R).
#    Writes outputs/species_runs/west/grid_<label>.rds.
#
# Plain R; one fit at a time; about 4 GB for V5.

suppressMessages({
  library(greta)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(readr)
  library(terra)
})
source("R/functions.R")
source("R/two_stage_predictions.R")
source("R/species_fit_helpers.R")

n_check <- 500
n_draws <- 160
n_rows <- 400
aggregation <- 4
windows <- list(`2010-2015` = 2010:2015, `2019-2025` = 2019:2025)
llin_pyrethroids <- c("Alpha-cypermethrin", "Deltamethrin", "Permethrin")

scratchpad <- paste0("/tmp/claude-1000/-home-nick-Dropbox-github-ir-cube/",
                     "be75c64a-3bb7-4b3e-a81c-c664fe72f5e2/scratchpad/species")
fits <- c(
  ref_f0 = paste0("../ir_cube_netscreen/outputs/pod_jobs/dh270_lin_f0_full/",
                  "temporary/fitted_model.RData"),
  V3f = "outputs/pod_jobs/sp_v3_floor/temporary/fitted_model.RData",
  V4_class = file.path(scratchpad,
                       "local_sp_v4_class/temporary/fitted_model.RData"),
  V5 = file.path(scratchpad, "local_sp_v5/temporary/fitted_model.RData"))
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) > 0) {
  fits <- character(0)
  for (argument in arguments) {
    parts <- strsplit(argument, "=", fixed = TRUE)[[1]]
    stopifnot(length(parts) == 2)
    fits[[parts[1]]] <- parts[2]
  }
}
stopifnot(all(file.exists(fits)))

output_dir <- "outputs/species_runs/west"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# the grid -------------------------------------------------------------------------

mask <- rast("data/clean/raster_mask.tif")
lattice <- expand.grid(col = seq(2, ncol(mask), by = aggregation),
                       row = seq(2, nrow(mask), by = aggregation))
grid_cells <- terra::cellFromRowCol(mask, lattice$row, lattice$col)
region <- as.character(terra::extract(rast("data/clean/region_raster.tif"),
                                      grid_cells)$region)
in_mask <- !is.na(terra::extract(mask, grid_cells)[, 1]) &
  !is.na(terra::extract(rast("data/clean/pfpr_water_mask.tif"),
                        grid_cells)[, 1])
grid_cells <- grid_cells[in_mask & region %in% "Western Africa"]
grid_country <- as.character(terra::extract(
  rast("data/clean/country_raster.tif"), grid_cells)$country_name)
stopifnot(!anyNA(grid_country))
report("grid: %d cells in the West within the mask, every %dth cell (%.2f degrees)",
       length(grid_cells), aggregation, aggregation * res(mask)[1])


# per fit ---------------------------------------------------------------------------

# A setup for two_stage_cells() and two_stage_chunk() without a two-stage
# fit: the dynamical model's initial states (draws x countries x types, the
# countries named), years and options
dynamical_setup <- function(fit, logit_init, years) {
  list(logit_init = logit_init,
       baseline_year = fit$baseline_year,
       years = years,
       design = fit$options$selection_columns,
       options = fit$options)
}

# Logit mortality draws of type k by the map path at `chunk`
# (two_stage_chunk()), a list by the setup's years of draws x cells
map_logit <- function(parameters, setup, chunk, k, logit_init) {
  out <- dynamical_logit_cells(
    parameters, k, matrix(logit_init[, chunk$country, k], parameters$n_draws),
    chunk$x, setup$years - setup$baseline_year + 1, x_init = chunk$x_init,
    share = chunk$share, kdr = chunk$kdr, basis = chunk$basis)
  setNames(out, setup$years)
}

for (label in names(fits)) {
  time_start <- Sys.time()
  fit <- load_fit(fits[[label]])
  df <- fit$df
  types <- fit$types

  # the misfit's bioassays are the fit's, row for row
  misfit <- read.csv(file.path("outputs/species_runs/misfit",
                               sprintf("%s_bioassays.csv", label)))
  stopifnot(nrow(misfit) == nrow(df),
            identical(as.numeric(misfit$cell), as.numeric(df$cell)),
            identical(as.integer(misfit$year_start),
                      as.integer(df$year_start)),
            identical(misfit$insecticide_type, df$insecticide_type),
            identical(as.numeric(misfit$died), as.numeric(df$died)))

  chosen <- even_draws(fit, n_check, label)
  parameters <- fit_parameter_draws(fit, chosen$index)
  set.seed(1)
  logit_init <- map_logit_init(parameters, fit$countries, fit$regions, df)
  report("%s: %d draws from chains %s; options: species %s, kdr %s, smooth %s, floor %s",
         label, parameters$n_draws, toString(unique(chosen$chain)),
         species_on(fit$options), kdr_on(fit$options),
         smooth_on(fit$options),
         if (!is.null(parameters$mortality_floor)) "constant" else
           if (!is.null(parameters$floor_intercept)) "intercepts" else
             if (!is.null(parameters$other_floor)) "per species" else "none")

  # 1. the check ------------------------------------------------------------------

  complex <- is.na(arabiensis_identified(df$species))
  west_pyrethroid <- which(complex & df$insecticide_type %in% llin_pyrethroids &
                             analysis_region(df$country_name, df$region) ==
                             "West")
  elsewhere <- setdiff(which(complex), west_pyrethroid)
  set.seed(47)
  rows <- sort(c(sample(west_pyrethroid, n_rows / 2),
                 sample(elsewhere, n_rows / 2)))
  check <- df[rows, ] %>%
    transmute(row = rows, cell, cell_id, year_start, insecticide_type,
              type_id, country_name, region = analysis_region(country_name,
                                                              region),
              misfit_predicted = misfit$predicted[rows],
              fit_country = data_cell_country(df)[cell_id],
              raster_country = as.character(terra::extract(
                rast("data/clean/country_raster.tif"), cell)$country_name))
  stopifnot(!anyNA(check$raster_country))
  check$map_fit_country <- check$map_raster_country <- check$map_draws <-
    check$mc_se <- NA_real_
  check_years <- sort(unique(check$year_start))
  setup <- dynamical_setup(fit, logit_init, check_years)
  map_subset <- round(seq(1, parameters$n_draws, length.out = n_draws))
  map_parameters <- subset_draws(parameters, map_subset)
  for (by in c("fit_country", "raster_country")) {
    cells <- two_stage_cells(setup, check$cell, check[[by]])
    for (k in sort(unique(check$type_id))) {
      at <- which(check$type_id == k)
      chunk <- two_stage_chunk(setup, cells, at)
      logit <- map_logit(parameters, setup, chunk, k, logit_init)
      logit_map <- map_logit(map_parameters, setup, chunk, k,
                             logit_init[map_subset, , , drop = FALSE])
      for (j in seq_along(at)) {
        year <- as.character(check$year_start[at[j]])
        p <- plogis(logit[[year]][, j])
        p_map <- plogis(logit_map[[year]][, j])
        check[[paste0("map_", by)]][at[j]] <- mean(p)
        if (by == "fit_country") {
          check$map_draws[at[j]] <- mean(p_map)
          check$mc_se[at[j]] <- sd(p_map) / sqrt(n_draws)
        }
      }
    }
  }
  check <- check %>%
    mutate(difference = map_fit_country - misfit_predicted,
           difference_map_draws = map_draws - misfit_predicted,
           z_map_draws = difference_map_draws / mc_se,
           difference_raster_country = map_raster_country - misfit_predicted)
  write.csv(check, file.path(output_dir, sprintf("check_%s.csv", label)),
            row.names = FALSE)
  report("%s check, %d bioassays (%d West LLIN pyrethroids): same %d draws, max abs diff in mortality %.2g; %d map draws, max abs diff %.3f, |z| median %.2f, max %.2f; raster country differs from the fit's at %d, max abs diff there %.3f",
         label, nrow(check), sum(check$row %in% west_pyrethroid),
         parameters$n_draws, max(abs(check$difference)), n_draws,
         max(abs(check$difference_map_draws)),
         median(abs(check$z_map_draws)), max(abs(check$z_map_draws)),
         sum(check$fit_country != check$raster_country),
         max(abs(check$difference_raster_country)))
  if (max(abs(check$difference)) > 1e-8) {
    stop(label, ": the map path does not reproduce the misfit's predictions")
  }

  # 2. the grid -------------------------------------------------------------------

  # the insecticides weighted by their mosquitoes tested in the West, all years
  type_weights <- df %>%
    filter(insecticide_type %in% llin_pyrethroids,
           analysis_region(country_name, region) == "West") %>%
    group_by(insecticide_type) %>%
    summarise(tested = sum(mosquito_number), .groups = "drop") %>%
    mutate(weight = tested / sum(tested))
  years <- unlist(windows, use.names = FALSE)
  setup <- dynamical_setup(fit, logit_init[map_subset, , , drop = FALSE],
                           years)
  cells <- two_stage_cells(setup, grid_cells, grid_country)
  window_of <- rep(names(windows), lengths(windows))
  out <- list()
  for (chunk_rows in split(seq_along(grid_cells),
                           ceiling(seq_along(grid_cells) / 2000))) {
    chunk <- two_stage_chunk(setup, cells, chunk_rows)
    combined <- setNames(vector("list", length(windows)), names(windows))
    per_type <- list()
    for (i in seq_len(nrow(type_weights))) {
      k <- match(type_weights$insecticide_type[i], types)
      logit <- map_logit(map_parameters, setup, chunk, k, setup$logit_init)
      for (w in names(windows)) {
        # draws x cells, the mean over the window's years
        p <- Reduce(`+`, lapply(logit[window_of == w], plogis)) /
          length(windows[[w]])
        per_type[[paste(type_weights$insecticide_type[i], w)]] <- colMeans(p)
        combined[[w]] <- if (is.null(combined[[w]])) {
          type_weights$weight[i] * p
        } else {
          combined[[w]] + type_weights$weight[i] * p
        }
      }
    }
    out[[length(out) + 1]] <- bind_cols(
      tibble(cell = grid_cells[chunk_rows]),
      as_tibble(lapply(combined, colMeans)) %>%
        rename_with(~ paste("mean", .x)),
      as_tibble(lapply(combined, col_sds)) %>%
        rename_with(~ paste("sd", .x)),
      as_tibble(per_type))
  }
  xy <- terra::xyFromCell(mask, grid_cells)
  grid <- bind_rows(out) %>%
    mutate(longitude = xy[, 1], latitude = xy[, 2], country = grid_country,
           .after = cell) %>%
    pivot_longer(-c(cell, longitude, latitude, country),
                 names_to = c("quantity", "window"),
                 names_pattern = "^(.*) (\\d{4}-\\d{4})$") %>%
    pivot_wider(names_from = quantity, values_from = value)
  saveRDS(list(label = label, grid = grid, type_weights = type_weights,
               n_draws = n_draws, chains = unique(chosen$chain),
               windows = windows, aggregation = aggregation),
          file.path(output_dir, sprintf("grid_%s.rds", label)))
  report("%s: grid of %d cells x %d windows; weights %s; West mean 2019-2025 %.2f; done in %.1f min, peak memory %.1f GB",
         label, length(grid_cells), length(windows),
         paste(sprintf("%s %.2f", type_weights$insecticide_type,
                       type_weights$weight), collapse = ", "),
         mean(grid$mean[grid$window == "2019-2025"]),
         as.numeric(difftime(Sys.time(), time_start, units = "mins")),
         peak_memory_gb())
  rm(fit, df, parameters, map_parameters, logit_init, setup, cells, out, grid)
  invisible(gc())
}
