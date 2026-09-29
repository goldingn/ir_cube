# Covariates of the dynamical model that are built the same way for the fits,
# the folds and the maps. Functions only; needs terra. Source from the repo
# root.


# initial-state covariates (#19) -----------------------------------------------

# Static predictors of the logit initial fraction susceptible, which let the
# initial state vary within a country. Each is standardised over every cell of
# the mask, so the values at a cell do not depend on which cells have data, and
# the fits, folds and maps share them:
#   log_pop_2000  log population per cell in 2000, the earliest layer of the
#                 population cube (there is none before 2000)
#   all_crops     agricultural productivity: the "all crops" yield layer of
#                 crop_group_scaled.tif, as log(x + 1e-4). The layer is 0-1
#                 scaled and very skewed (56% zeros, mean 0.003, sd 0.012), so
#                 standardised untransformed it reaches 47 sd at the data
#                 cells; 1e-4 is about the 10th percentile of its non-zero
#                 values, and on this scale the data cells span -0.7 to 4.1 sd
init_covariate_names <- c("log_pop_2000", "all_crops")

init_covariate_layers <- function() {
  mask <- rast("data/clean/raster_mask.tif")
  log_pop <- log(rast("data/clean/pop_cube.tif")[["pop_2000"]])
  crops <- log(rast("data/clean/crop_group_scaled.tif")[["all crops"]] + 1e-4)
  layers <- terra::mask(c(log_pop, crops), mask)
  names(layers) <- init_covariate_names
  moments <- terra::global(layers, c("mean", "sd"), na.rm = TRUE)
  (layers - moments$mean) / moments$sd
}

# The standardised initial-state covariates at mask cells `cells`, as a
# cells x covariates matrix with named columns.
init_covariate_matrix <- function(cells, layers = init_covariate_layers()) {
  x <- as.matrix(terra::extract(layers, cells))
  colnames(x) <- names(layers)
  x
}


# selection covariates (#23) ---------------------------------------------------

# How the selection design matrix is built: the population column, raw or log
# (both min-max scaled to 0-1), and hinge columns min(x, k) for any of the
# time-varying covariates (nets, irs, pop, the last after its transform), with
# knots k on the covariate's 0-1 scale. Each hinge column enters the fitness
# with a positive coefficient, like the linear ones, so the effect of the
# covariate can saturate but not decrease. The default is the design of the
# fits to date: raw population, no hinges.
#   pop     "raw": pop_scaled_cube.tif, population min-max scaled;
#           "log": pop_log_scaled_cube.tif, log population min-max scaled
#           (write_log_pop_cubes())
#   hinges  named list of knots, e.g. list(nets = c(0.18, 0.35))
selection_design <- function(pop = c("raw", "log"), hinges = list()) {
  pop <- match.arg(pop)
  stopifnot(is.list(hinges),
            all(names(hinges) %in% c("nets", "irs", "pop")),
            all(vapply(hinges, function(k) {
              is.numeric(k) && all(k > 0 & k < 1)
            }, logical(1))))
  list(pop = pop, hinges = hinges)
}

# the columns of the selection design matrix, in order: the time-varying
# covariates, their hinge columns, then the static crop covariates
selection_column_names <- function(design = selection_design()) {
  pop_name <- c(raw = "pop", log = "log_pop")[[design$pop]]
  hinge_names <- unlist(lapply(names(design$hinges), function(name) {
    base <- if (name == "pop") pop_name else name
    paste0(base, "_min_", design$hinges[[name]])
  }))
  c("nets", "irs", pop_name, hinge_names, colnames(selection_flat_names()))
}

selection_flat_names <- function() {
  matrix(nrow = 0, ncol = 10,
         dimnames = list(NULL, c("all crops", "cereal crops", "root crops",
                                 "pulse crops", "oil crops", "fibre crops",
                                 "other crops", "cotton", "vegetables",
                                 "rice")))
}

# The time-varying selection covariates at mask cells `cells`, years
# baseline_year..end_year, as a cells x years x columns array: the cubes are
# padded back to baseline_year and forward to end_year by repeating their
# first and last layers, then the hinge columns are added.
selection_time_varying <- function(cells, baseline_year, end_year,
                                   design = selection_design()) {
  files <- c(nets = "data/clean/net_use_cube.tif",
             irs = "data/clean/irs_coverage_scaled_cube.tif",
             pop = c(raw = "data/clean/pop_scaled_cube.tif",
                     log = "data/clean/pop_log_scaled_cube.tif")[[design$pop]])
  if (!file.exists(files[["pop"]])) {
    stop(files[["pop"]], " not found; run write_log_pop_cubes()")
  }
  read_cube <- function(file) {
    cube <- rast(file)
    cube <- suppressWarnings(pre_pad_cube(cube, baseline_year))
    cube <- suppressWarnings(post_pad_cube(cube, end_year))
    years <- as.numeric(str_sub(names(cube), start = -4L))
    cube <- cube[[years >= baseline_year & years <= end_year]]
    stopifnot(identical(as.numeric(str_sub(names(cube), start = -4L)),
                        as.numeric(baseline_year:end_year)))
    as.matrix(terra::extract(cube, cells))
  }
  columns <- selection_column_names(design)
  n_time_varying <- 3 + length(unlist(design$hinges))
  out <- array(NA_real_,
               c(length(cells), end_year - baseline_year + 1, n_time_varying),
               dimnames = list(NULL, baseline_year:end_year,
                               columns[seq_len(n_time_varying)]))
  for (i in 1:3) {
    out[, , i] <- read_cube(files[[i]])
  }
  j <- 3
  for (name in names(design$hinges)) {
    for (k in design$hinges[[name]]) {
      j <- j + 1
      out[, , j] <- pmin(out[, , match(name, c("nets", "irs", "pop"))], k)
    }
  }
  out
}

# The static selection covariates (crop yields) at mask cells `cells`, as a
# cells x 10 matrix: the crop group totals, then cotton, vegetables and rice
selection_flat <- function(cells) {
  crops_group <- rast("data/clean/crop_group_scaled.tif")
  crops_all <- rast("data/clean/crop_scaled.tif")
  covs_flat <- c(crops_group,
                 crops_all$cotton,
                 crops_all$vegetables,
                 crops_all$rice)
  flat <- as.matrix(terra::extract(covs_flat, cells))
  stopifnot(identical(colnames(flat), colnames(selection_flat_names())))
  flat
}

# The selection design matrix at mask cells `cells` (cell_id = position in
# `cells`) for years baseline_year..end_year: x_cell_years, one row per
# (cell_id, year_id), cell-major, and its cell_years_index, as
# build_dynamical_model() takes them.
selection_design_matrix <- function(cells, baseline_year, end_year,
                                    design = selection_design()) {
  time_varying <- selection_time_varying(cells, baseline_year, end_year,
                                         design)
  flat <- selection_flat(cells)
  n_years <- dim(time_varying)[2]
  # cells x years x columns to (years x cells) x columns, year fastest
  long <- matrix(aperm(time_varying, c(2, 1, 3)),
                 ncol = dim(time_varying)[3])
  x_cell_years <- cbind(long, flat[rep(seq_along(cells), each = n_years), ,
                                   drop = FALSE])
  colnames(x_cell_years) <- selection_column_names(design)
  list(x_cell_years = x_cell_years,
       cell_years_index = tibble(
         cell_id = rep(seq_along(cells), each = n_years),
         year_id = rep(seq_len(n_years), length(cells))))
}

# Write min-max scaled log population cubes, pop_log_scaled_cube.tif
# (2000-2022) and pop_log_scaled_cube_future.tif (2023-2030), layers
# log_pop_<year>, as prep_rasters.R does with its log transform switched on,
# without rerunning it. prep_rasters.R scales over its whole population series
# (2000-2049), whose range is not in the saved cubes; it is recovered exactly
# from the saved raw and scaled cubes, a linear map of each other, and log is
# monotone, so the log range is the log of that range.
write_log_pop_cubes <- function() {
  pop <- rast("data/clean/pop_cube.tif")
  pop_scaled <- rast("data/clean/pop_scaled_cube.tif")
  pop_future <- rast("data/clean/pop_cube_future.tif")
  values <- na.omit(cbind(raw = c(terra::values(pop[[1]])),
                          scaled = c(terra::values(pop_scaled[[1]]))))
  fit <- lm(raw ~ scaled, data = as.data.frame(values))
  pop_min <- unname(coef(fit)[1])
  pop_max <- pop_min + unname(coef(fit)[2])
  stopifnot(max(abs(residuals(fit))) < 1e-6 * pop_max, pop_min > 0)
  log_range <- log(c(pop_min, pop_max))
  scale_log <- function(cube) {
    out <- (log(cube) - log_range[1]) / diff(log_range)
    names(out) <- str_replace(names(cube), "^pop_", "log_pop_")
    out
  }
  terra::writeRaster(scale_log(pop), "data/clean/pop_log_scaled_cube.tif",
                     overwrite = TRUE)
  terra::writeRaster(scale_log(pop_future),
                     "data/clean/pop_log_scaled_cube_future.tif",
                     overwrite = TRUE)
  invisible(log_range)
}
