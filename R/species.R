# The An. arabiensis share of the An. gambiae complex (#47): the map of the
# arabiensis fraction r(x), the arabiensis share of each bioassay, and the
# helper the species model uses to mix its two trajectories (arabiensis and
# the other members of the complex). Functions and settings only; needs
# terra. Sourced by R/dynamical_model.R.


# settings ---------------------------------------------------------------------

# r(x), the arabiensis fraction of the complex: static, on the model grid
# (R/prep_arabiensis_fraction.R)
arabiensis_fraction_file <- "data/clean/arabiensis_fraction.tif"

# Settings of the species model, for dynamical_model_options(species = ).
# Arabiensis shares every parameter of the other members, with a multiplier
# exp(gamma_selection) on its log fitness from selection, a multiplier
# exp(gamma_cost) on the fitness cost (with reversion), and its own mortality
# floor (species_mortality(), R/dynamical_model.R). The other members' floor
# is the model's mortality_floor option.
#   arabiensis_floor  TRUE to estimate a mortality floor for arabiensis,
#                     FALSE for none
#   floor_prior       the Beta shape parameters of its prior, by default
#                     those of mortality_floor (dynamical_variables())
#   fraction_file     the map of r(x)
species_options <- function(arabiensis_floor = TRUE,
                            floor_prior = c(1, 49),
                            fraction_file = arabiensis_fraction_file) {
  list(arabiensis_floor = arabiensis_floor,
       floor_prior = floor_prior,
       fraction_file = fraction_file)
}

# whether a model's options have the species model on. Options saved before
# #47 have no species element, and are off
species_on <- function(options) {
  !is.null(options$species) && !isFALSE(options$species)
}

check_species_options <- function(species) {
  if (isFALSE(species)) {
    return(invisible(species))
  }
  stopifnot(
    is.list(species),
    setequal(names(species), names(species_options())),
    isTRUE(species$arabiensis_floor) || isFALSE(species$arabiensis_floor),
    is.numeric(species$floor_prior), length(species$floor_prior) == 2,
    all(species$floor_prior > 0),
    is.character(species$fraction_file), length(species$fraction_file) == 1)
  invisible(species)
}


# the map of r(x) -------------------------------------------------------------

# `values` (one per cell of the grid `template`, NA where missing) with the
# missing cells among `cells` given the value of their nearest non-missing
# cell, by distance on the grid with longitude scaled by the cosine of the
# latitude. The nearest non-missing cell always has a missing (or off-grid)
# neighbour, as a neighbour one step closer would otherwise be non-missing
# and nearer, so only those cells are searched.
fill_from_nearest_cells <- function(values, template, cells) {
  missing <- cells[is.na(values[cells])]
  if (length(missing) == 0) {
    return(values)
  }
  have <- terra::rast(template, nlyrs = 1)
  terra::values(have) <- as.numeric(!is.na(values))
  neighbours <- terra::values(terra::focal(have, w = 3, fun = "sum",
                                           na.rm = TRUE), mat = FALSE)
  edge <- which(!is.na(values) & neighbours < 9)
  xy_edge <- terra::xyFromCell(template, edge)
  xy_missing <- terra::xyFromCell(template, missing)
  for (i in seq_along(missing)) {
    scale <- cos(xy_missing[i, 2] * pi / 180)
    distance <- ((xy_edge[, 1] - xy_missing[i, 1]) * scale) ^ 2 +
      (xy_edge[, 2] - xy_missing[i, 2]) ^ 2
    values[missing[i]] <- values[edge[which.min(distance)]]
  }
  values
}

# r(x) at every cell of the model grid, with every mask cell filled
# (fill_from_nearest_cells()), as a vector indexed by cell number; cached for
# the session, as a plain vector so that forked workers can use it. The map
# has no value at 10,176 of the 1,479,742 mask cells (lake edges and coasts,
# and St Helena, 1,300 km from the nearest value), which hold 253 of the
# 27,865 modelled bioassays (41 cells)
arabiensis_fraction_cache <- new.env()
arabiensis_fraction_values <- function(file = arabiensis_fraction_file) {
  if (is.null(arabiensis_fraction_cache[[file]])) {
    if (!file.exists(file)) {
      stop(file, " not found; run R/prep_arabiensis_fraction.R")
    }
    mask <- terra::rast("data/clean/raster_mask.tif")
    layer <- terra::rast(file)
    stopifnot(terra::compareGeom(layer, mask))
    arabiensis_fraction_cache[[file]] <- fill_from_nearest_cells(
      terra::values(layer, mat = FALSE), layer, terra::cells(mask))
  }
  arabiensis_fraction_cache[[file]]
}

# r(x) at mask cells `cells` (cell numbers of data/clean/raster_mask.tif)
arabiensis_fraction_at <- function(cells, file = arabiensis_fraction_file) {
  r <- arabiensis_fraction_values(file)[cells]
  stopifnot(!anyNA(r))
  r
}

# r(x) at mask cells `cells`, the share of the complex-wide predictions of a
# fit with `options`, or NULL without the species model, as
# dynamical_logit_cells() takes `share`
prediction_share <- function(options, cells) {
  if (!species_on(options)) {
    return(NULL)
  }
  arabiensis_fraction_at(cells, options$species$fraction_file)
}


# the share of a bioassay -----------------------------------------------------

# What each species record of the modelled bioassays says about the arabiensis
# share: 1 for arabiensis, 0 for another member of the complex, NA for a
# record of the complex. A species is recorded only for a pool of one species,
# by studies that targeted it. "Anopheles gambiae" (144 modelled records, all
# from IR Mapper's species table) is NA: several of its studies are of
# An. gambiae s.l. by their titles, and some are where s.s. is rare (Zambia,
# Malawi).
complex_species_share <- c("Anopheles arabiensis" = 1,
                           "Anopheles gambiae s.s." = 0,
                           "Anopheles gambiae ss" = 0,
                           "Anopheles coluzzii" = 0,
                           "Anopheles merus" = 0,
                           "Anopheles quadriannulatus" = 0,
                           "Anopheles gambiae" = NA,
                           "gambiae complex" = NA)

arabiensis_identified <- function(species) {
  unknown <- setdiff(unique(species), names(complex_species_share))
  if (length(unknown) > 0) {
    stop("species not in complex_species_share: ", toString(unknown))
  }
  unname(complex_species_share[species])
}

# The arabiensis share a of each of `rows` (bioassays with species and cell),
# for a model with `options`: 1 or 0 where the species is recorded, and r(x)
# at the cell for a record of the complex. This assumes, untested, that the
# mosquitoes collected for a bioassay of the complex are a sample of the
# complex at that place with no bias towards either species; the studies that
# record a species targeted it, so they carry no information on that bias.
arabiensis_share <- function(rows, options) {
  if (is.null(rows$species) || is.null(rows$cell)) {
    stop("the species model (#47) needs the species and cell of each row")
  }
  share <- arabiensis_identified(rows$species)
  complex <- is.na(share)
  share[complex] <- arabiensis_fraction_at(rows$cell[complex],
                                           options$species$fraction_file)
  share
}


# mixing the trajectories ------------------------------------------------------

# Bioassay mortality is the mixture a pA + (1 - a) pG of the two species'
# mortality, for share a: each mosquito of the pool is arabiensis with
# probability a, independently, so the number dying is binomial with that
# mean, and the model's beta-binomial takes it as its mean.

# log(w exp(x) + (1 - w) exp(y)), elementwise, for weights w in (0, 1):
# shifted by the larger of x and y, so that the larger term is exp(0) = 1 and
# neither overflows or both underflow
weighted_log_sum_exp <- function(x, y, w) {
  m <- pmax(x, y)
  m + log(w * exp(x - m) + (1 - w) * exp(y - m))
}

# The logit of the mixture share pA + (1 - share) pG of two draws x n
# matrices of logits, share one per column, computed from the log
# probabilities of both outcomes so that it keeps its precision where p
# rounds to 0 or 1 (the logit draws reach beyond +-37). Where share is 0 or 1
# it is exactly the logit of that component.
mixture_logit <- function(logit_a, logit_g, share) {
  stopifnot(identical(dim(logit_a), dim(logit_g)),
            length(share) %in% c(1, ncol(logit_a)))
  w <- matrix(share, nrow(logit_a), ncol(logit_a), byrow = TRUE)
  out <- weighted_log_sum_exp(plogis(logit_a, log.p = TRUE),
                              plogis(logit_g, log.p = TRUE), w) -
    weighted_log_sum_exp(plogis(-logit_a, log.p = TRUE),
                         plogis(-logit_g, log.p = TRUE), w)
  out[w == 0] <- logit_g[w == 0]
  out[w == 1] <- logit_a[w == 1]
  out
}
