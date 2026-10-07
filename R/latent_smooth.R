# The latent spatial smooths of the dynamical model (V5, #47): two static
# surfaces fitted inside the model in place of the kdr covariate, u_s(x)
# scaling the strength of selection and u_f(x) shifting the mortality floor.
# Each is a low-rank Hilbert-space approximate Gaussian process (HSGP; Solin
# and Sarkka 2020, Riutort-Mayol et al. 2022) on a box around the modelled
# bioassay cells, in the Albers equal-area projection of Africa of the
# two-stage correction, in units of 1,000 km:
#   u(x) = sum_j sqrt(S(omega_j; sd, rho)) beta_j phi_j(x),  beta_j ~ N(0, 1)
# with phi_j the eigenfunctions of the Laplacian on the box (Dirichlet
# boundary), omega_j the square roots of their eigenvalues, and S the
# kernel's spectral density. The basis at the cells is data; only the
# weights depend on sd and rho. Functions and settings only; needs sf and
# terra. Sourced by R/dynamical_model.R, after R/kdr_covariate.R.


# settings ---------------------------------------------------------------------

# Settings of the latent smooths, for dynamical_model_options(smooth = ).
#   selection         the smooth u_s(x) of the strength of selection: the
#                     cumulative log fitness is multiplied by exp(u_s(x))
#                     (outer_mortality(), R/dynamical_model.R). TRUE (the
#                     default) for every class, "class" for the pyrethroids
#                     and DDT only (smooth_classes), FALSE for none
#   floor             the smooth u_f(x) of the logit mortality floor, f =
#                     plogis(floor_intercept + u_f(x)): "class" (the
#                     default) for the pyrethroids and DDT only, TRUE for
#                     every class, FALSE for none. Needs mortality_floor =
#                     TRUE
#   floor_intercepts  "class" (the default) for one floor_intercept per
#                     insecticide class, "one" for one for all. With
#                     mortality_floor = TRUE the floor is plogis(
#                     floor_intercept + u_f(x)), in place of the constant
#                     mortality_floor, with or without u_f; floor_intercept is
#                     the logit floor where u_f is 0, which is u_f's mean
#                     over the modelled cells
#   kernel            "matern52" (the default; Matern, smoothness 5/2) or
#                     "se" (squared exponential)
#   c                 the box's half-widths as multiples of the half-extents
#                     of the modelled bioassay cells; widened where needed to
#                     cover every cell of the mask (smooth_box())
#   m                 the number of basis functions per dimension (x, y), or
#                     NULL (the default) for the rule of Riutort-Mayol et al.
#                     at the prior's lower range (smooth_box())
#   range_prior       the penalised-complexity prior of each smooth's range
#                     rho (Fuglstad et al. 2019): P(rho < range_prior[1]) =
#                     range_prior[2], rho in 1,000 km
#   sd_prior          and of its marginal sd: P(sd > sd_prior[1]) =
#                     sd_prior[2]
# build_dynamical_model() adds the basis (smooth_box()): the projection, the
# box's origin and half-widths, m, the frequency indices of the basis
# functions and their means over the modelled cells (which centre them); and
# term_classes, whether each class (in class_id order) is one a "class"
# smooth applies to
smooth_options <- function(selection = TRUE,
                           floor = "class",
                           floor_intercepts = "class",
                           kernel = "matern52",
                           c = 1.5,
                           m = NULL,
                           range_prior = c(1, 0.05),
                           sd_prior = c(1, 0.05)) {
  list(selection = selection,
       floor = floor,
       floor_intercepts = floor_intercepts,
       kernel = kernel,
       c = c,
       m = m,
       range_prior = range_prior,
       sd_prior = sd_prior)
}

# the insecticide classes with the smooth terms where a smooth is "class":
# kdr gives resistance to the pyrethroids and DDT
smooth_classes <- kdr_floor_classes

# the elements build_dynamical_model() adds to the options (smooth_box())
smooth_basis_elements <- c("crs", "origin", "half_width", "indices",
                           "column_means", "term_classes")

# the projection of the coordinates: the two-stage correction's Albers
# equal-area conic for Africa (africa_equal_area_crs,
# R/two_stage_correction.R), in km
smooth_crs <- paste(
  "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25",
  "+x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"
)

# whether a model's options have the latent smooths on. Options saved before
# them have no smooth element, and are off
smooth_on <- function(options) {
  !is.null(options$smooth) && !isFALSE(options$smooth)
}

# the smooths a model with `options` has, of "selection" and "floor"
smooth_kinds <- function(options) {
  if (!smooth_on(options)) {
    return(character(0))
  }
  c("selection", "floor")[c(!isFALSE(options$smooth$selection),
                            !isFALSE(options$smooth$floor))]
}

# whether a model's options have the floor of the smooth model,
# plogis(floor_intercept + u_f(x)), with or without u_f
smooth_floor_on <- function(options) {
  smooth_on(options) && isTRUE(options$mortality_floor)
}

check_smooth_options <- function(smooth) {
  if (isFALSE(smooth) || is.null(smooth)) {
    return(invisible(smooth))
  }
  is_kind <- function(x) isFALSE(x) || isTRUE(x) || identical(x, "class")
  stopifnot(
    is.list(smooth),
    all(setdiff(names(smooth_options()), "m") %in% names(smooth)),
    all(names(smooth) %in% c(names(smooth_options()), smooth_basis_elements)),
    is_kind(smooth$selection), is_kind(smooth$floor),
    !isFALSE(smooth$selection) || !isFALSE(smooth$floor),
    smooth$floor_intercepts %in% c("class", "one"),
    smooth$kernel %in% names(hsgp_m_factor),
    is.numeric(smooth$c), length(smooth$c) == 1, smooth$c > 1,
    is.null(smooth$m) || (is.numeric(smooth$m) && length(smooth$m) == 2 &&
                            all(smooth$m >= 1)),
    is.numeric(smooth$range_prior), length(smooth$range_prior) == 2,
    all(smooth$range_prior > 0), smooth$range_prior[2] < 1,
    is.numeric(smooth$sd_prior), length(smooth$sd_prior) == 2,
    all(smooth$sd_prior > 0), smooth$sd_prior[2] < 1)
  invisible(smooth)
}


# priors -----------------------------------------------------------------------

# The rates of the penalised-complexity priors (Fuglstad et al. 2019) of a
# smooth's range rho and marginal sd in two dimensions: 1 / rho ~
# Exponential(range), i.e. density (range) rho^-2 exp(-range / rho), and sd ~
# Exponential(sd), with P(rho < rho_0) = alpha_rho and P(sd > sd_0) =
# alpha_sd. Derived for Matern fields; used for the squared exponential too
smooth_prior_rates <- function(smooth) {
  list(range = -log(smooth$range_prior[2]) * smooth$range_prior[1],
       sd = -log(smooth$sd_prior[2]) / smooth$sd_prior[1])
}


# the basis --------------------------------------------------------------------

# The number of basis functions per dimension needed for an accurate
# approximation at lengthscale ell, as a multiple of L / ell for box
# half-width L: the rules of Riutort-Mayol et al. (2022), m = 2.65 c S / ell
# for the Matern 5/2 and 1.75 c S / ell for the squared exponential, for data
# in [-S, S] and L = c S
hsgp_m_factor <- c(matern52 = 2.65, se = 1.75)

# Projected coordinates, in 1,000 km, of longitudes and latitudes
smooth_coords <- function(lon, lat, crs = smooth_crs) {
  xy <- sf::sf_project(from = "+proj=longlat +datum=WGS84 +no_defs",
                       to = crs, pts = cbind(lon, lat))
  xy / 1000
}

# Projected coordinates, in 1,000 km, of mask cells `cells` (cell numbers of
# data/clean/raster_mask.tif)
smooth_cell_coords <- function(cells, crs = smooth_crs) {
  xy <- terra::xyFromCell(terra::rast("data/clean/raster_mask.tif"), cells)
  smooth_coords(xy[, 1], xy[, 2], crs)
}

# the range of the projected coordinates of every mask cell, cached for the
# session
smooth_mask_cache <- new.env()
smooth_mask_range <- function(crs = smooth_crs) {
  if (is.null(smooth_mask_cache[[crs]])) {
    mask <- terra::rast("data/clean/raster_mask.tif")
    smooth_mask_cache[[crs]] <- apply(
      smooth_cell_coords(terra::cells(mask), crs), 2, range)
  }
  smooth_mask_cache[[crs]]
}

# The smooth options with the basis added, for the modelled bioassay cells
# `cells` (mask cells, each once) and the classes `classes` (in class_id
# order). The box is centred on the cells' midrange, its half-width in each
# dimension c times their half-extent, or more where the mask reaches beyond
# that, so that every mask cell can be predicted to. Unless given, m is
# hsgp_m_factor[kernel] L / ell_min per dimension, rounded up, for ell_min the
# lengthscale at the prior's lower range (rho = 2 ell: the Matern's
# correlation is about 0.13 at rho, and the squared exponential's 0.135). Of
# the m[1] x m[2] products of the one-dimensional eigenfunctions, only those
# with frequency up to the lower of the two dimensions' highest, omega_max =
# min_d pi m_d / (2 L_d), are kept (by the rule, about pi hsgp_m_factor /
# (2 ell_min) in both): the spectral density is isotropic, so the grid's
# corners beyond omega_max have less power than the frequencies just beyond
# the grid along each axis, which the rule already leaves out. The basis
# functions are ordered by frequency. A saved basis is kept, so a rebuilt fit
# has its own. With the defaults (c = 1.5, 1.64 in y to cover the mask; m =
# (29, 25), 539 of 725 kept), the covariance the basis implies at the
# modelled cells, centred, is within 2.4% of sd^2 of the Matern's for ranges
# of 1,000-2,000 km, but 21% at 500 km, where the basis is too coarse, and
# 18% at 4,000 km and 39% at 8,000 km, where the box is too narrow; c = 2
# (857 kept) gives 3.3% at 4,000 km and 20% at 8,000 km
# (R/check_latent_smooth.R).
smooth_box <- function(smooth, cells, classes) {
  smooth$term_classes <- classes %in% smooth_classes
  if (!is.null(smooth$indices)) {
    return(smooth)
  }
  smooth$crs <- smooth_crs
  coords <- smooth_cell_coords(cells, smooth$crs)
  lower <- apply(coords, 2, min)
  upper <- apply(coords, 2, max)
  smooth$origin <- (lower + upper) / 2
  mask_range <- smooth_mask_range(smooth$crs)
  reach <- pmax(abs(mask_range[1, ] - smooth$origin),
                abs(mask_range[2, ] - smooth$origin))
  smooth$half_width <- pmax(smooth$c * (upper - lower) / 2, reach)
  if (is.null(smooth$m)) {
    ell_min <- smooth$range_prior[1] / 2
    smooth$m <- ceiling(hsgp_m_factor[[smooth$kernel]] * smooth$half_width /
                          ell_min)
  }
  indices <- as.matrix(expand.grid(x = seq_len(smooth$m[1]),
                                   y = seq_len(smooth$m[2])))
  omega <- hsgp_frequencies(indices, smooth$half_width)
  omega_max <- min(pi * smooth$m / (2 * smooth$half_width))
  keep <- omega <= omega_max * (1 + 1e-12)
  smooth$indices <- indices[keep, , drop = FALSE][order(omega[keep]), ,
                                                  drop = FALSE]
  smooth$column_means <- colMeans(hsgp_basis(coords, smooth))
  smooth
}

# The frequencies omega_j = sqrt(lambda_j) of the basis functions with
# indices `indices` (one row per function, one column per dimension) on a box
# of half-widths `half_width`: lambda_j = sum_d (pi j_d / (2 L_d))^2
hsgp_frequencies <- function(indices, half_width) {
  sqrt(rowSums(sweep(indices, 2, pi / (2 * half_width), FUN = "*") ^ 2))
}

# The uncentred basis at projected coordinates `coords` (points x 2, in 1,000
# km), as a points x basis functions matrix: phi_j(x) = prod_d sqrt(1 / L_d)
# sin(pi j_d (x_d + L_d) / (2 L_d)), x relative to the box's origin
hsgp_basis <- function(coords, smooth) {
  x <- sweep(matrix(coords, ncol = 2), 2, smooth$origin)
  L <- smooth$half_width
  if (any(abs(x) > rep(L, each = nrow(x)) * (1 + 1e-9))) {
    stop("points outside the box of the latent smooths")
  }
  phi <- 1
  for (d in 1:2) {
    phi <- phi * sin(outer(x[, d] + L[d],
                           pi * smooth$indices[, d] / (2 * L[d]))) /
      sqrt(L[d])
  }
  phi
}

# The centred basis at projected coordinates `coords`: each function less its
# mean over the modelled bioassay cells, so that each smooth has mean 0 there
# and does not trade off with the overall strength of selection or the floor
# intercepts
smooth_basis_at <- function(smooth, coords) {
  sweep(hsgp_basis(coords, smooth), 2, smooth$column_means)
}

# The centred basis of a model with `options` at mask cells `cells`, as a
# cells x basis functions matrix, or NULL without the latent smooths
prediction_basis <- function(options, cells) {
  if (!smooth_on(options)) {
    return(NULL)
  }
  if (is.null(options$smooth$indices)) {
    stop("the smooth options have no basis; build_dynamical_model() adds it")
  }
  smooth_basis_at(options$smooth, smooth_cell_coords(cells, options$smooth$crs))
}


# the weights ------------------------------------------------------------------

# The square root of the kernel's spectral density in two dimensions at
# angular frequencies `omega`, for marginal sd `sd` and range rho = 1 /
# inv_range, lengthscale ell = rho / 2 (Rasmussen and Williams 2006, section
# 4.2, in angular frequency, with k(r) = (2 pi)^-2 int S(w) exp(i w.r) dw;
# R/check_latent_smooth.R checks the covariance they imply):
#   Matern 5/2  S(w) = sd^2 10 pi 5^(5/2) ell^2 (5 + ell^2 w^2)^(-7/2)
#   SE          S(w) = sd^2 2 pi ell^2 exp(-ell^2 w^2 / 2)
# For greta arrays (sd and inv_range scalars, giving one per frequency) or
# plain R (sd and inv_range one per draw, giving draws x frequencies)
smooth_sqrt_spectral <- function(omega, sd, inv_range, kernel) {
  ell <- 1 / (2 * inv_range)
  scaled <- if (inherits(ell, "greta_array")) ell ^ 2 * omega ^ 2 else
    outer(ell ^ 2, omega ^ 2)
  shape <- switch(kernel,
                  matern52 = sqrt(10 * pi) * 5 ^ (5 / 4) * (5 + scaled) ^ (-7 / 4),
                  se = sqrt(2 * pi) * exp(-scaled / 4))
  sd * ell * shape
}

# The weights sqrt(S(omega_j)) beta_j of a smooth's basis functions, from its
# standard normal raw weights `raw`, sd and inverse range. For greta arrays
# (raw one per basis function) or plain R (raw draws x basis functions, sd
# and inv_range one per draw)
smooth_weights <- function(raw, sd, inv_range, smooth) {
  omega <- hsgp_frequencies(smooth$indices, smooth$half_width)
  raw * smooth_sqrt_spectral(omega, sd, inv_range, smooth$kernel)
}

# The names of a smooth's variables (dynamical_variables()): its raw
# weights, sd and inverse range
smooth_variable_names <- function(kind) {
  parts <- c("raw", "sd", "inv_range")
  setNames(paste0("smooth_", parts, "_", kind), parts)
}

# The weights of each smooth of a model with `options`, from its variables
# `v` (greta arrays, or draws x dim arrays in plain R), as a named list
smooth_weight_terms <- function(v, options) {
  out <- list()
  for (kind in smooth_kinds(options)) {
    names <- smooth_variable_names(kind)
    raw <- v[[names[["raw"]]]]
    sd <- v[[names[["sd"]]]]
    inv_range <- v[[names[["inv_range"]]]]
    if (!inherits(raw, "greta_array")) {
      raw <- matrix(raw, nrow = length(c(sd)))
      sd <- c(sd)
      inv_range <- c(inv_range)
    }
    out[[kind]] <- smooth_weights(raw, sd, inv_range, options$smooth)
  }
  out
}


# the smooths at the rows ------------------------------------------------------

# The floor intercept of each class `class` (class_ids): its own, or the one
# for all
smooth_intercept_index <- function(options, class) {
  if (identical(options$smooth$floor_intercepts, "class")) class else
    rep(1L, length(class))
}

# Whether each class `class` has the smooth `kind`: every class (TRUE), or
# the pyrethroids and DDT ("class"), as 1 or 0
smooth_class_weight <- function(options, kind, class) {
  if (isTRUE(options$smooth[[kind]])) {
    return(rep(1, length(class)))
  }
  as.numeric(options$smooth$term_classes[class])
}

# The floor of the smooth model at rows with floor intercepts `intercept`
# (that of each row's class) and floor smooth `u` (0 where the row's class has
# none; NULL for no floor smooth): plogis(intercept + u). For greta arrays
# (one element per row) or plain R (intercept one per draw, u draws x cells,
# giving draws x cells, or one per draw without u)
smooth_floor_value <- function(intercept, u = NULL) {
  logit_floor <- if (is.null(u)) intercept else intercept + u
  if (inherits(logit_floor, "greta_array")) ilogit(logit_floor) else
    plogis(logit_floor)
}

# The posterior mean and sd of each latent smooth of a fit at projected
# coordinates `coords` (points x 2, in 1,000 km; smooth_coords()), from its
# parameters (dynamical_parameter_draws(), R/dynamical_predictions.R), as a
# list named by smooth of points x 2 matrices, columns mean and sd; the basis
# is made `chunk` points at a time
smooth_posterior_at <- function(parameters, coords, chunk = 5000) {
  smooth <- parameters$options$smooth
  out <- lapply(parameters$smooth_weights, function(weights) {
    matrix(NA_real_, nrow(coords), 2, dimnames = list(NULL, c("mean", "sd")))
  })
  for (rows in split(seq_len(nrow(coords)),
                     ceiling(seq_len(nrow(coords)) / chunk))) {
    basis <- smooth_basis_at(smooth, coords[rows, , drop = FALSE])
    for (kind in names(out)) {
      u <- parameters$smooth_weights[[kind]] %*% t(basis)
      mean <- colMeans(u)
      out[[kind]][rows, "mean"] <- mean
      out[[kind]][rows, "sd"] <- sqrt(colSums(sweep(u, 2, mean) ^ 2) /
                                        (nrow(u) - 1))
    }
  }
  out
}
