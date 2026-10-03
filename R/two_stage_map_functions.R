# Helpers for running the dynamical model on the full prediction grid, for the
# two-stage maps (R/two_stage_maps.R) and the supplement figures
# (R/fig_two_stage_components.R). The recursion and the initial conditions are
# those of R/dynamical_predictions.R; this adds only the grid's covariates and
# the country -> region lookup for countries without data. Functions only:
# source from the repo root after R/dynamical_predictions.R.

# All covariates the dynamical model uses, at the cells `cells` of the mask and
# for the years baseline_year..end_year, built the way predict.R builds
# x_cell_years_predict: the cubes are padded back to the baseline year and
# forward to end_year by repeating their first and last layers. Returned as a
# cells x years x 3 array of the time-varying covariates (nets, irs, pop) and a
# cells x 10 matrix of the static crop covariates, rather than predict.R's long
# (cell, year) matrix (5.5 GB for the 1.48M cells and 36 years); map_x()
# assembles a chunk of cells in the column order of x_cell_years
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
  time_varying <- array(NA_real_, c(length(cells), ncol(nets), 3),
                        dimnames = list(NULL, baseline_year:end_year,
                                        c("nets", "irs", "pop")))
  time_varying[, , 1] <- nets
  rm(nets)
  time_varying[, , 2] <- read_cube("data/clean/irs_coverage_scaled_cube.tif")
  time_varying[, , 3] <- read_cube("data/clean/pop_scaled_cube.tif")

  crops_group <- rast("data/clean/crop_group_scaled.tif")
  crops_all <- rast("data/clean/crop_scaled.tif")
  flat <- as.matrix(terra::extract(c(crops_group, crops_all$cotton,
                                     crops_all$vegetables, crops_all$rice),
                                   cells))

  list(time_varying = time_varying, flat = flat)
}

# The covariates of rows `rows` of map_covariates() for its first n_years
# years, as the cells x years x n_covs array dynamical_logit_cells() takes
map_x <- function(covariates, rows, n_years) {
  flat <- covariates$flat[rows, , drop = FALSE]
  x <- array(NA_real_, c(length(rows), n_years, 3 + ncol(flat)))
  for (t in seq_len(n_years)) {
    x[, t, ] <- cbind(matrix(covariates$time_varying[rows, t, ],
                             nrow = length(rows)), flat)
  }
  x
}

# Draws of logit q_0 for every African country in the UNSD lookup, not only
# those with data (dynamical_logit_init(): fresh deviations for countries and
# regions without data, with the caller's RNG), as a draws x countries x types
# array with the country names as dimnames[[2]]. The region of a country is
# the UNSD one, as in predict.R; the fit took each country's region from its
# first record, and the two must agree for every country with data
map_logit_init <- function(parameters, countries, regions, df,
                           lookup = country_region_lookup()) {
  all_countries <- unique(lookup$country_name)
  unsd_region <- lookup$region[match(all_countries, lookup$country_name)]

  fit_region <- df %>%
    distinct(country_id, .keep_all = TRUE) %>%
    arrange(country_id)
  stopifnot(identical(fit_region$country_id, seq_along(countries)),
            identical(regions[fit_region$region_id],
                      unsd_region[match(countries, all_countries)]))

  new_regions <- setdiff(unique(unsd_region), regions)
  region <- match(unsd_region, c(regions, new_regions))
  logit_init <- dynamical_logit_init(parameters,
                                     match(all_countries, countries), region)
  dimnames(logit_init) <- list(NULL, all_countries, NULL)
  logit_init
}
