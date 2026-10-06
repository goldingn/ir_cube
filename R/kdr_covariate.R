# The kdr covariate of the dynamical model (#47): total kdr (995F + 995S) in
# 2015 from the joint allele model (branch latent-kdr), a static map per
# species group, as a covariate of the strength of selection and of the
# fitness cost. Functions and settings only; needs terra. Sourced by
# R/dynamical_model.R, after R/species.R.

# the map: bands "complex" (the whole complex), "arabiensis" and "other" (the
# other members of the complex), every mask cell of the model grid
kdr_total_file <- "data/clean/kdr_total_2015.tif"

# Settings of the kdr covariate, for dynamical_model_options(kdr = ). With it,
# the cumulative log fitness and the reversion of each trajectory are
# multiplied by exp(delta_selection k(x)) and exp(delta_cost k(x)) for a cell's
# standardised kdr k(x) (outer_mortality(), R/dynamical_model.R): without the
# species model, one pair of slopes and the "complex" band; with it, a pair
# for each species (kdr_slope_names) and each species' own band.
#   file    the map
#   floor   FALSE (the default) for a constant mortality floor (or none);
#           TRUE for a floor that depends on kdr, plogis(floor_intercept +
#           floor_kdr k(x)), in place of the constant mortality_floor; "class"
#           for one floor_intercept per insecticide class, with the kdr term
#           only for the classes kdr gives resistance to (kdr_floor_classes:
#           pyrethroids and DDT). Both need mortality_floor = TRUE and no
#           species model (check_dynamical_model_options())
# build_dynamical_model() adds centre and scale, the mean and sd of
# logit(kdr) of the "complex" band over the modelled bioassay cells, which
# standardise every band, so that the species' k are on one scale and their
# slopes comparable; and with floor = "class", floor_classes, whether each
# class (in class_id order) has the kdr term
kdr_options <- function(file = kdr_total_file, floor = FALSE) {
  list(file = file, floor = floor)
}

# the insecticide classes whose floor depends on kdr with floor = "class":
# kdr gives resistance to the pyrethroids and DDT
kdr_floor_classes <- c("Pyrethroids", "Organochlorines")

# whether a model's options have the kdr-dependent floor, and its kind: FALSE,
# TRUE (one intercept) or "class". kdr options saved before it have no floor
# element, and have none
kdr_floor <- function(options) {
  if (!kdr_on(options) || is.null(options$kdr$floor)) FALSE else
    options$kdr$floor
}

# The slopes of the kdr covariate (dynamical_variables()): one pair without
# the species model, and a pair per species with it; the cost slopes only
# with reversion
kdr_slope_names <- c("delta_selection", "delta_cost",
                     "delta_selection_other", "delta_cost_other",
                     "delta_selection_arabiensis", "delta_cost_arabiensis")

# whether a model's options have the kdr covariate on. Options saved before
# it have no kdr element, and are off
kdr_on <- function(options) {
  !is.null(options$kdr) && !isFALSE(options$kdr)
}

check_kdr_options <- function(kdr) {
  if (isFALSE(kdr) || is.null(kdr)) {
    return(invisible(kdr))
  }
  stopifnot(
    is.list(kdr),
    all(names(kdr_options()) %in% names(kdr)),
    all(names(kdr) %in% c(names(kdr_options()), "centre", "scale",
                          "floor_classes")),
    is.character(kdr$file), length(kdr$file) == 1,
    isFALSE(kdr$floor) || isTRUE(kdr$floor) || identical(kdr$floor, "class"),
    is.null(kdr$centre) || (is.numeric(kdr$centre) && length(kdr$centre) == 1),
    is.null(kdr$scale) || (is.numeric(kdr$scale) && length(kdr$scale) == 1 &&
                             kdr$scale > 0))
  invisible(kdr)
}

# the bands a model with `options` uses: each species' with the species model,
# else the whole complex's
kdr_bands <- function(options) {
  if (species_on(options)) c("other", "arabiensis") else "complex"
}

# p clamped to [eps, 1 - eps], so that its logit is finite. The map lies in
# 0.003-0.98, so this only guards against a map with values of 0 or 1
clamp_probability <- function(p, eps = 1e-4) {
  pmin(pmax(p, eps), 1 - eps)
}

# logit(kdr) of band `band` at mask cells `cells` (cell numbers of
# data/clean/raster_mask.tif), any missing cell filled from the nearest
# (filled_layer_values(), R/species.R)
kdr_logit_at <- function(cells, band, file = kdr_total_file) {
  values <- filled_layer_values(file, band)[cells]
  stopifnot(!anyNA(values))
  qlogis(clamp_probability(values))
}

# The kdr options with the standardisation added: the mean and sd of logit
# kdr of the "complex" band over the mask cells `cells` (the modelled
# bioassay cells, each once)
standardise_kdr <- function(kdr, cells) {
  l <- kdr_logit_at(cells, "complex", kdr$file)
  kdr$centre <- mean(l)
  kdr$scale <- stats::sd(l)
  kdr
}

# The kdr-dependent floor (kdr_options(floor = )) at rows with floor
# intercepts `intercept` (the intercept of each row's class), slope `slope`
# and standardised kdr `k` (0 where the row's class has no kdr term):
# plogis(intercept + slope k). For greta arrays (one element per row) or
# plain R (intercept and slope one per draw, k one per cell, giving draws x
# cells)
kdr_floor_value <- function(intercept, slope, k) {
  if (inherits(slope, "greta_array")) {
    return(ilogit(intercept + slope * k))
  }
  plogis(intercept + outer(slope, k))
}

# The standardised kdr k(x) of a model with `options` at mask cells `cells`,
# as a cells x bands matrix (kdr_bands()), or NULL without the kdr covariate
prediction_kdr <- function(options, cells) {
  if (!kdr_on(options)) {
    return(NULL)
  }
  kdr <- options$kdr
  if (is.null(kdr$centre)) {
    stop("the kdr options have no standardisation; build_dynamical_model() ",
         "adds it")
  }
  bands <- kdr_bands(options)
  out <- vapply(bands, function(band) {
    (kdr_logit_at(cells, band, kdr$file) - kdr$centre) / kdr$scale
  }, numeric(length(cells)))
  matrix(out, length(cells), length(bands), dimnames = list(NULL, bands))
}
