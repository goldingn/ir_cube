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

# How the selection design matrix is built. Net use and IRS enter as they are,
# with any hinge columns. Population and the crop layers are proxies for
# insecticide use outside vector control; with a trend they enter only as
# products z(x, t) = g(t) s(x, t) with a temporal trend g of that pressure, so
# that no covariate acts as a constant selection rate and there is no
# selection from them where g is 0. There are two trends: g_dom for
# population (domestic and public-health use) and g_ag for the crops
# (agricultural use).
#   pop          transform of population s(x, t), each on 0-1, with d the
#                population density in people per km2 (pop_cube.tif over the
#                cell area, 0 in the pixels WorldPop has as empty;
#                pop_density_matrix()):
#                "encounter": 1 - exp(-d / lambda), lambda =
#                pop_d_half / log(2), the chance of encountering at least one
#                source of exposure when sources are Poisson in number with
#                mean proportional to d;
#                "saturating": d / (d + pop_d_half);
#                "raw": pop_scaled_cube.tif, population per cell min-max
#                scaled; "log": pop_log_scaled_cube.tif, log population
#                min-max scaled (write_log_pop_cubes())
#   pop_d_half   the density (people per km2) at which "encounter" and
#                "saturating" are 0.5; not used by the others
#   hinges       named list of knots, e.g. list(nets = c(0.18, 0.35)): hinge
#                columns min(x, k) for any of nets, irs, pop (the last after
#                its transform), knots k on the 0-1 scale. Each enters the
#                fitness with a positive coefficient, like the linear ones, so
#                the effect can saturate but not decrease
#   trend_pop    g_dom(t), multiplying population and its hinges, and
#   trend_crops  g_ag(t), multiplying the crops, each one of: "none", the
#                columns as they are (the design of the fits to date);
#                "linear_0_1", linear in the year, 0 in trend_years[1] and 1
#                in trend_years[2]; or a regions x years matrix of g, rownames
#                the region names of country_region_lookup() and colnames the
#                years, covering every year the design is built for, each
#                cell taking its region's row (cell_regions()). g_ag is to be
#                estimated by region from FAOSTAT
#                (data/clean/faostat_regional_trend.csv)
#   trend_years  the years where a linear trend is 0 and 1. The first should
#                be the model's baseline year (1995), so these columns are 0
#                there
#   trend_after  a linear trend after trend_years[2]: "continue", on the same
#                line (g > 1, e.g. 1.17 in 2030); "cap", held at 1
# Net use, IRS and their hinges have no trend.
selection_design <- function(pop = c("raw", "log", "encounter",
                                     "saturating"),
                             pop_d_half = 50,
                             hinges = list(),
                             trend_pop = "none",
                             trend_crops = "none",
                             trend_years = c(1995, 2025),
                             trend_after = c("continue", "cap")) {
  pop <- match.arg(pop)
  trend_after <- match.arg(trend_after)
  stopifnot(is.list(hinges),
            all(names(hinges) %in% c("nets", "irs", "pop")),
            all(vapply(hinges, function(k) {
              is.numeric(k) && all(k > 0 & k < 1)
            }, logical(1))),
            is.numeric(pop_d_half), length(pop_d_half) == 1,
            pop_d_half > 0,
            is.numeric(trend_years), length(trend_years) == 2,
            trend_years[2] > trend_years[1])
  for (trend in list(trend_pop, trend_crops)) {
    if (is.matrix(trend)) {
      stopifnot(is.numeric(trend), !is.null(rownames(trend)),
                !anyNA(suppressWarnings(as.integer(colnames(trend)))))
    } else {
      stopifnot(length(trend) == 1, trend %in% c("none", "linear_0_1"))
    }
  }
  list(pop = pop, pop_d_half = pop_d_half, hinges = hinges, trend_pop = trend_pop,
       trend_crops = trend_crops, trend_years = trend_years,
       trend_after = trend_after)
}

# The design of the fits before the trend and the saturating population (raw
# population, no hinges, no trend), e.g. for matching the cached inits
selection_design_untrended <- function() {
  selection_design(pop = "raw", trend_pop = "none", trend_crops = "none")
}

# A saved design, completed as selection_design() would build it. Designs
# saved before the trends were added have only pop and hinges, and meant none.
complete_selection_design <- function(design) {
  for (name in c("trend_pop", "trend_crops")) {
    if (is.null(design[[name]])) {
      design[[name]] <- "none"
    }
  }
  do.call(selection_design, design)
}

# whether a trend (design$trend_pop or design$trend_crops) is on
has_trend <- function(trend) {
  is.matrix(trend) || trend != "none"
}

# the columns of the selection design matrix, in order: nets, irs, pop, the
# hinge columns, then the crops; with their trend, the pop and population
# hinge columns are named "<name>:g_dom" (e.g. "pop_enc:g_dom",
# "pop_enc_min_0.5:g_dom") and the crops "<name>:g_ag" ("all crops:g_ag")
selection_column_names <- function(design = selection_design()) {
  design <- complete_selection_design(design)
  pop_name <- c(raw = "pop", log = "log_pop",
                encounter = "pop_enc", saturating = "pop_sat")[[design$pop]]
  trended <- function(name, trend, suffix) {
    if (has_trend(trend)) paste0(name, ":", suffix) else name
  }
  hinge_names <- unlist(lapply(names(design$hinges), function(name) {
    if (name == "pop") {
      trended(paste0(pop_name, "_min_", design$hinges[[name]]),
              design$trend_pop, "g_dom")
    } else {
      paste0(name, "_min_", design$hinges[[name]])
    }
  }))
  c("nets", "irs", trended(pop_name, design$trend_pop, "g_dom"), hinge_names,
    trended(colnames(selection_flat_names()), design$trend_crops, "g_ag"))
}

selection_flat_names <- function() {
  matrix(nrow = 0, ncol = 10,
         dimnames = list(NULL, c("all crops", "cereal crops", "root crops",
                                 "pulse crops", "oil crops", "fibre crops",
                                 "other crops", "cotton", "vegetables",
                                 "rice")))
}

# A cube as a cells x years matrix at mask cells `cells`, years
# baseline_year..end_year, padded back to baseline_year and forward to
# end_year by repeating its first and last layers
read_padded_cube <- function(cube, cells, baseline_year, end_year) {
  if (is.character(cube)) {
    cube <- rast(cube)
  }
  cube <- suppressWarnings(pre_pad_cube(cube, baseline_year))
  cube <- suppressWarnings(post_pad_cube(cube, end_year))
  years <- as.numeric(str_sub(names(cube), start = -4L))
  cube <- cube[[years >= baseline_year & years <= end_year]]
  stopifnot(identical(as.numeric(str_sub(names(cube), start = -4L)),
                      as.numeric(baseline_year:end_year)))
  as.matrix(terra::extract(cube, cells))
}

# Population density, people per km2, at mask cells `cells` for years
# baseline_year..end_year (padded as read_padded_cube()), cells x years.
# prep_rasters.R filled the pixels WorldPop has as empty (0 or NA on land)
# with 1 person per cell, relevelled by a factor close to 1 in each year;
# these are set back to 0. The fill is each layer's most common value near 1
# (51,408 pixels in every year, against at most 3 for any other value there).
pop_density_matrix <- function(cells, baseline_year, end_year) {
  pop <- rast("data/clean/pop_cube.tif")
  fill <- vapply(seq_len(terra::nlyr(pop)), function(i) {
    v <- terra::values(pop[[i]], mat = FALSE)
    v <- v[!is.na(v) & v > 0.99 & v < 1.01]
    u <- unique(v)
    n <- tabulate(match(v, u))
    stopifnot(max(n) > 10000)
    u[which.max(n)]
  }, numeric(1))
  pop <- pop * (pop != fill)
  names(pop) <- names(rast("data/clean/pop_cube.tif"))
  area <- terra::cellSize(rast("data/clean/raster_mask.tif"), unit = "km")
  area_cells <- terra::extract(area, cells)[, 1]
  read_padded_cube(pop, cells, baseline_year, end_year) / area_cells
}

# The population column s(x, t) of the design, cells x years
selection_pop_matrix <- function(cells, baseline_year, end_year, design) {
  switch(
    design$pop,
    raw = read_padded_cube("data/clean/pop_scaled_cube.tif", cells,
                           baseline_year, end_year),
    log = {
      file <- "data/clean/pop_log_scaled_cube.tif"
      if (!file.exists(file)) {
        stop(file, " not found; run write_log_pop_cubes()")
      }
      read_padded_cube(file, cells, baseline_year, end_year)
    },
    encounter = {
      d <- pop_density_matrix(cells, baseline_year, end_year)
      -expm1(-d * log(2) / design$pop_d_half)
    },
    saturating = {
      d <- pop_density_matrix(cells, baseline_year, end_year)
      d / (d + design$pop_d_half)
    }
  )
}

# The region of each mask cell in `cells`, from the country raster and the
# UNSD lookup as predict.R assigns them; NA outside any country
cell_regions <- function(cells) {
  country <- as.character(terra::extract(
    rast("data/clean/country_raster.tif"), cells)$country_name)
  lookup <- country_region_lookup()
  lookup$region[match(country, lookup$country_name)]
}

# A trend g(t) (design$trend_pop or design$trend_crops) at mask cells `cells`
# for years baseline_year..end_year, as a cells x years matrix (see
# selection_design()). With a supplied regional trend, cells outside any
# region are NA.
selection_trend_matrix <- function(cells, baseline_year, end_year, trend,
                                   design) {
  years <- baseline_year:end_year
  if (is.matrix(trend)) {
    missing_years <- setdiff(years, as.integer(colnames(trend)))
    if (length(missing_years) > 0) {
      stop("the supplied trend has no values for ",
           toString(missing_years))
    }
    region <- cell_regions(cells)
    missing_regions <- setdiff(na.omit(unique(region)), rownames(trend))
    if (length(missing_regions) > 0) {
      stop("the supplied trend has no row for ", toString(missing_regions))
    }
    return(trend[match(region, rownames(trend)),
                 match(years, as.integer(colnames(trend))),
                 drop = FALSE])
  }
  stopifnot(trend == "linear_0_1")
  g <- (years - design$trend_years[1]) / diff(design$trend_years)
  g <- pmax(g, 0)
  if (design$trend_after == "cap") {
    g <- pmin(g, 1)
  }
  matrix(g, length(cells), length(years), byrow = TRUE)
}

# The time-varying selection covariates at mask cells `cells`, years
# baseline_year..end_year, as a cells x years x columns array: nets, irs and
# pop from the cubes (padded as read_padded_cube()), then the hinge columns,
# with trend_pop, pop and its hinges multiplied by g_dom, and with
# trend_crops, the crop products g_ag x crop after them (selection_static()
# then has none).
selection_time_varying <- function(cells, baseline_year, end_year,
                                   design = selection_design()) {
  design <- complete_selection_design(design)
  columns <- selection_column_names(design)
  n_time_varying <- 3 + length(unlist(design$hinges)) +
    if (has_trend(design$trend_crops)) ncol(selection_flat_names()) else 0
  out <- array(NA_real_,
               c(length(cells), end_year - baseline_year + 1, n_time_varying),
               dimnames = list(NULL, baseline_year:end_year,
                               columns[seq_len(n_time_varying)]))
  out[, , 1] <- read_padded_cube("data/clean/net_use_cube.tif", cells,
                                 baseline_year, end_year)
  out[, , 2] <- read_padded_cube("data/clean/irs_coverage_scaled_cube.tif",
                                 cells, baseline_year, end_year)
  out[, , 3] <- selection_pop_matrix(cells, baseline_year, end_year, design)
  j <- 3
  for (name in names(design$hinges)) {
    for (k in design$hinges[[name]]) {
      j <- j + 1
      out[, , j] <- pmin(out[, , match(name, c("nets", "irs", "pop"))], k)
    }
  }
  if (has_trend(design$trend_pop)) {
    g <- selection_trend_matrix(cells, baseline_year, end_year,
                                design$trend_pop, design)
    pop_columns <- c(3, 3 + which(rep(names(design$hinges),
                                      lengths(design$hinges)) == "pop"))
    for (p in pop_columns) {
      out[, , p] <- out[, , p] * g
    }
  }
  if (has_trend(design$trend_crops)) {
    g <- selection_trend_matrix(cells, baseline_year, end_year,
                                design$trend_crops, design)
    flat <- selection_flat(cells)
    for (i in seq_len(ncol(flat))) {
      j <- j + 1
      out[, , j] <- flat[, i] * g
    }
  }
  stopifnot(j == n_time_varying)
  out
}

# The crop layers at mask cells `cells`, as a cells x 10 matrix: the crop
# group totals, then cotton, vegetables and rice
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

# The static columns of the design at mask cells `cells`: the crop layers
# without trend_crops, and none (a cells x 0 matrix) with it, when the crop
# products are among the time-varying columns
selection_static <- function(cells, design = selection_design()) {
  design <- complete_selection_design(design)
  if (has_trend(design$trend_crops)) {
    return(matrix(numeric(0), length(cells), 0))
  }
  selection_flat(cells)
}

# The selection design matrix at mask cells `cells` (cell_id = position in
# `cells`) for years baseline_year..end_year: x_cell_years, one row per
# (cell_id, year_id), cell-major, and its cell_years_index, as
# build_dynamical_model() takes them.
selection_design_matrix <- function(cells, baseline_year, end_year,
                                    design = selection_design()) {
  time_varying <- selection_time_varying(cells, baseline_year, end_year,
                                         design)
  flat <- selection_static(cells, design)
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
