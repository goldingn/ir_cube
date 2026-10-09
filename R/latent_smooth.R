# The latent spatial smooths of the dynamical model (V5, #47): two static
# surfaces fitted inside the model in place of the kdr covariate, u_s(x)
# scaling the strength of selection and u_f(x) shifting the mortality floor,
# and optionally a third, u_init(x), of the initial state, in place of its
# hierarchy of regions and countries (smooth_options(init = TRUE)).
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
# The defaults are V5h's (#47), the default model's: the range fixed at
# 1,500 km and a half-normal prior of scale 0.5 on each smooth's sd. V5 is
# smooth_options(range = NULL, sd_prior = c(1, 0.05)), and V5r
# smooth_options(sd_prior = c(1, 0.05)).
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
#   shear             FALSE (the default), or TRUE for a selection smooth that
#                     shares the floor's: u_s(x) = v_s(x) + b u_f(x), v_s and
#                     u_f independent smooths and b (smooth_shear) one
#                     estimated loading, N(1, 0.5) a priori, before the
#                     class weighting of u_s. Absorbs the positive trade-off
#                     between the log selection multiplier and the logit
#                     floor, so that v_s and u_f are nearly uncorrelated in
#                     the posterior. Needs both smooths
#   floor_intercepts  "class" (the default) for one floor_intercept per
#                     insecticide class, "one" for one for all. With
#                     mortality_floor = TRUE the floor is plogis(
#                     floor_intercept + u_f(x)), in place of the constant
#                     mortality_floor, with or without u_f; floor_intercept is
#                     the logit floor where u_f is 0, which is u_f's mean
#                     over the modelled cells
#   kernel            "se" (the default; squared exponential) or "matern52"
#                     (Matern, smoothness 5/2)
#   c                 the box's half-widths as multiples of the half-extents
#                     of the modelled bioassay cells; widened where needed to
#                     cover every cell of the mask (smooth_box())
#   m                 the number of basis functions per dimension (x, y), or
#                     NULL (the default) for the rule of Riutort-Mayol et al.
#                     at range basis_range (smooth_box())
#   range             the range rho of both smooths, in 1,000 km, fixed: 1.5
#                     by default (V5r, V5h); or NULL to estimate each
#                     smooth's with range_prior (V5). A fixed range leaves no
#                     range variable in the model: the spectral weights use
#                     it as a constant. 1,500 km because the spatial
#                     correlation range of recent pyrethroid mortality is
#                     1,230 km [840, 1,800], and of block plateaus 1,500 km
#                     [1,100, 2,000] (R/plateau_range.R)
#   basis_range       the shortest range, in 1,000 km, the rule sets m for:
#                     the fixed range if there is one (by default 1.5), or
#                     else 1 (1,000 km); NULL for the fixed range, or else
#                     the prior's lower range, range_prior[1]
#                     (smooth_basis_range())
#   range_prior       the penalised-complexity prior of each smooth's range
#                     rho (Fuglstad et al. 2019): P(rho < range_prior[1]) =
#                     range_prior[2], rho in 1,000 km; unused with a fixed
#                     range
#   sd_prior          the prior of its marginal sd: list(family =
#                     "half_normal", scale = s) for sd ~ N(0, s^2) truncated
#                     to sd > 0, by default with scale 0.5 (V5h); or c(sd_0,
#                     alpha) for the penalised-complexity prior, exponential
#                     with P(sd > sd_0) = alpha (V5 and V5r, c(1, 0.05);
#                     smooth_sd_prior_family())
#   init              FALSE (the default), or TRUE for a smooth u_init(x) of
#                     the logit relative initial state (the position between
#                     init_frac_min and 1; init_frac_constants(),
#                     R/dynamical_model.R) in place of the hierarchy of
#                     regions and countries, with a loading lambda >= 0 per
#                     insecticide type:
#                       logit_init_relative(x, type) = logit_init_mean[type] +
#                         lambda[type] u_init(x) + the type's initial-state
#                         covariate effects at x
#                     u_init is one field shared by every type, with the
#                     other smooths' kernel, basis and range, but its marginal
#                     sd fixed at 1: the loadings carry the scale, each with
#                     the prior of a smooth's sd (sd_prior), so each is its
#                     type's sd of the field, and a type whose initial states
#                     do not follow the shared pattern shrinks towards 0.
#                     Positive loadings leave no sign-flip modes. The
#                     covariates are centred at their mean over the modelled
#                     cells (build_dynamical_model()), so logit_init_mean is
#                     the mean logit relative initial state over them
# build_dynamical_model() adds the basis (smooth_box()): the projection, the
# box's origin and half-widths, m, the frequency indices of the basis
# functions and their means over the modelled cells (which centre them); and
# term_classes, whether each class (in class_id order) is one a "class"
# smooth applies to
smooth_options <- function(selection = TRUE,
                           floor = "class",
                           shear = FALSE,
                           floor_intercepts = "class",
                           kernel = "se",
                           c = 2,
                           m = NULL,
                           range = 1.5,
                           basis_range = if (is.null(range)) 1 else range,
                           range_prior = c(1.5, 0.05),
                           sd_prior = list(family = "half_normal",
                                           scale = 0.5),
                           init = FALSE) {
  list(selection = selection,
       floor = floor,
       shear = shear,
       floor_intercepts = floor_intercepts,
       kernel = kernel,
       c = c,
       m = m,
       range = range,
       basis_range = basis_range,
       range_prior = range_prior,
       sd_prior = sd_prior,
       init = init)
}

# whether the smooths' range is fixed (smooth_options(range = )). Options
# saved before the option have no range element, and estimate it. Read as
# smooth[["range"]]: smooth$range would match range_prior partially
smooth_range_fixed <- function(smooth) {
  is.list(smooth) && !is.null(smooth[["range"]])
}

# The range the basis rule sets m for, in 1,000 km (smooth_box()):
# basis_range, or if it is NULL the fixed range, or else the prior's lower
# range
smooth_basis_range <- function(smooth) {
  if (!is.null(smooth$basis_range)) {
    return(smooth$basis_range)
  }
  if (smooth_range_fixed(smooth)) smooth[["range"]] else
    smooth$range_prior[1]
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

# the smooths a model with `options` has, of "selection", "floor" and "init"
smooth_kinds <- function(options) {
  if (!smooth_on(options)) {
    return(character(0))
  }
  c("selection", "floor", "init")[c(!isFALSE(options$smooth$selection),
                                    !isFALSE(options$smooth$floor),
                                    smooth_init_on(options))]
}

# whether a model's options have the smooth of the initial state, u_init(x),
# in place of the hierarchy of regions and countries (smooth_options(init =
# TRUE)). Options saved before it have no init element, and do not
smooth_init_on <- function(options) {
  smooth_on(options) && isTRUE(options$smooth[["init"]])
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
    all(setdiff(names(smooth_options()),
                c("m", "range", "basis_range", "shear", "init")) %in%
          names(smooth)),
    all(names(smooth) %in% c(names(smooth_options()), smooth_basis_elements)),
    is_kind(smooth$selection), is_kind(smooth$floor),
    is.null(smooth[["init"]]) || isFALSE(smooth[["init"]]) ||
      isTRUE(smooth[["init"]]),
    !isFALSE(smooth$selection) || !isFALSE(smooth$floor) ||
      isTRUE(smooth[["init"]]),
    is.null(smooth$shear) || isFALSE(smooth$shear) || isTRUE(smooth$shear),
    !isTRUE(smooth$shear) ||
      (!isFALSE(smooth$selection) && !isFALSE(smooth$floor)),
    smooth$floor_intercepts %in% c("class", "one"),
    smooth$kernel %in% names(hsgp_m_factor),
    is.numeric(smooth$c), length(smooth$c) == 1, smooth$c > 1,
    is.null(smooth$m) || (is.numeric(smooth$m) && length(smooth$m) == 2 &&
                            all(smooth$m >= 1)),
    !smooth_range_fixed(smooth) || (is.numeric(smooth[["range"]]) &&
                                      length(smooth[["range"]]) == 1 &&
                                      is.finite(smooth[["range"]]) &&
                                      smooth[["range"]] > 0),
    is.null(smooth$basis_range) || (is.numeric(smooth$basis_range) &&
                                      length(smooth$basis_range) == 1 &&
                                      smooth$basis_range > 0),
    is.numeric(smooth$range_prior), length(smooth$range_prior) == 2,
    all(smooth$range_prior > 0), smooth$range_prior[2] < 1,
    smooth_sd_prior_valid(smooth$sd_prior))
  invisible(smooth)
}

# whether `prior` is a valid sd_prior (smooth_options()): c(sd_0, alpha) with
# sd_0 > 0 and 0 < alpha < 1, or list(family = "half_normal", scale = s) with
# s > 0
smooth_sd_prior_valid <- function(prior) {
  if (is.list(prior)) {
    return(setequal(names(prior), c("family", "scale")) &&
             identical(prior$family, "half_normal") &&
             is.numeric(prior$scale) && length(prior$scale) == 1 &&
             is.finite(prior$scale) && prior$scale > 0)
  }
  is.numeric(prior) && length(prior) == 2 && all(prior > 0) && prior[2] < 1
}


# priors -----------------------------------------------------------------------

# The rates of the penalised-complexity priors (Fuglstad et al. 2019) of a
# smooth's range rho and marginal sd in two dimensions: 1 / rho ~
# Exponential(range), i.e. density (range) rho^-2 exp(-range / rho), and sd ~
# Exponential(sd), with P(rho < rho_0) = alpha_rho and P(sd > sd_0) =
# alpha_sd. Derived for Matern fields; used for the squared exponential too.
# With the half-normal prior of the sd, the sd rate is NA; the
# smooth_sd_prior_*() functions below handle either prior of the sd
smooth_prior_rates <- function(smooth) {
  list(range = -log(smooth$range_prior[2]) * smooth$range_prior[1],
       sd = if (smooth_sd_prior_family(smooth) == "pc") {
         -log(smooth$sd_prior[2]) / smooth$sd_prior[1]
       } else {
         NA_real_
       })
}

# The family of the prior of each smooth's marginal sd (smooth_options(
# sd_prior = )): "pc", the penalised-complexity prior, sd ~ Exponential (the
# default), or "half_normal", sd ~ N(0, scale^2) truncated to sd > 0. V5h's,
# with scale 0.5, is a standard weakly informative penalising prior for an
# sd, half-normal with scale 0.5 (Nick's choice), much stricter in the tail
# than the PC prior, which the data overrode (V5r's sds ran to 4-13): P(sd >
# 1) = 0.046 under both, but P(sd > 2) = 6e-5 against 0.0025
smooth_sd_prior_family <- function(smooth) {
  if (is.list(smooth$sd_prior)) smooth$sd_prior$family else "pc"
}

# The prior of a smooth's sd as a greta distribution of dimension `dim` (the
# sd, or the loadings per type of the smooth of the initial state, which
# have the sd's prior). Either way greta's free state is log sd
smooth_sd_prior_distribution <- function(smooth, dim = 1) {
  if (smooth_sd_prior_family(smooth) == "half_normal") {
    return(normal(0, smooth$sd_prior$scale, dim = dim,
                  truncation = c(0, Inf)))
  }
  exponential(smooth_prior_rates(smooth)$sd, dim = dim)
}

# The log density of the prior of a smooth's sd at `sd`, in plain R, as greta
# evaluates it: greta divides a truncated density by the probability of the
# interval, 1/2 for the half-normal, so its density is 2 N(sd; 0, s^2)
smooth_sd_log_prior <- function(sd, smooth) {
  if (smooth_sd_prior_family(smooth) == "half_normal") {
    return(log(2) + dnorm(sd, 0, smooth$sd_prior$scale, log = TRUE))
  }
  dexp(sd, smooth_prior_rates(smooth)$sd, log = TRUE)
}

# P(sd > x) under the prior of a smooth's sd
smooth_sd_prior_tail <- function(x, smooth) {
  if (smooth_sd_prior_family(smooth) == "half_normal") {
    return(2 * pnorm(-x / smooth$sd_prior$scale))
  }
  exp(-smooth_prior_rates(smooth)$sd * x)
}

# The quantile at probability p of the prior of a smooth's sd
smooth_sd_prior_quantile <- function(p, smooth) {
  if (smooth_sd_prior_family(smooth) == "half_normal") {
    return(smooth$sd_prior$scale * qnorm((1 + p) / 2))
  }
  qexp(p, smooth_prior_rates(smooth)$sd)
}

# A short description of the prior of a smooth's sd, for printed summaries
smooth_sd_prior_label <- function(smooth) {
  if (smooth_sd_prior_family(smooth) == "half_normal") {
    return(sprintf("sd ~ half-normal, scale %g", smooth$sd_prior$scale))
  }
  sprintf("sd ~ Exponential(%.3f), P(sd > %g) = %g",
          smooth_prior_rates(smooth)$sd, smooth$sd_prior[1],
          smooth$sd_prior[2])
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
# lengthscale at range smooth_basis_range(): basis_range, by default 1,000 km,
# or the fixed range (rho = 2 ell: the Matern's correlation is about 0.13 at
# rho, and the squared exponential's 0.135). Of
# the m[1] x m[2] products of the one-dimensional eigenfunctions, only those
# with frequency up to the lower of the two dimensions' highest, omega_max =
# min_d pi m_d / (2 L_d), are kept (by the rule, about pi hsgp_m_factor /
# (2 ell_min) in both): the spectral density is isotropic, so the grid's
# corners beyond omega_max have less power than the frequencies just beyond
# the grid along each axis, which the rule already leaves out. The basis
# functions are ordered by frequency. A saved basis is kept, so a rebuilt fit
# has its own. With the squared exponential, c = 2 and m set at 1,000 km
# (basis_range 1, as for V5's estimated range: m = (25, 20), 363 of 500
# kept), the covariance the basis implies at the modelled cells, centred, is
# within 2.4% of sd^2 of the kernel's for ranges of 1,000-4,000 km (0.02% at
# 1,500 km), but 27% at 8,000 km and 22% at 16,000 km, where the box is too
# narrow for the range. m set at 1,500 km ((17, 14), 164 of 238 kept), as
# for the fixed range of the defaults (V5r, V5h), is within 2.3% at 1,500 km
# (0.8% with all 238), but 19% at 1,000 km (R/check_latent_smooth.R).
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
    ell_min <- smooth_basis_range(smooth) / 2
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
# For greta arrays (sd a scalar, inv_range a scalar greta array or, for a
# fixed range, a number, giving one per frequency) or plain R (sd and
# inv_range one per draw, giving draws x frequencies)
smooth_sqrt_spectral <- function(omega, sd, inv_range, kernel) {
  ell <- 1 / (2 * inv_range)
  greta <- inherits(ell, "greta_array") || inherits(sd, "greta_array")
  scaled <- if (greta) ell ^ 2 * omega ^ 2 else outer(ell ^ 2, omega ^ 2)
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
  spectral <- smooth_sqrt_spectral(omega, sd, inv_range, smooth$kernel)
  # with sd and range both fixed (u_init at a fixed range), the spectral
  # weights are constants, one per frequency
  if (inherits(raw, "greta_array") && !inherits(spectral, "greta_array")) {
    spectral <- c(spectral)
  }
  raw * spectral
}

# The names of a smooth's variables (dynamical_variables()): its raw
# weights, sd and inverse range, and its loadings per type. With a fixed range
# (smooth_range_fixed()), the model has no inverse range variable; the smooth
# of the initial state ("init") has loadings and no sd (its sd is fixed at 1,
# smooth_sd()), and the others an sd and no loadings
smooth_variable_names <- function(kind) {
  parts <- c("raw", "sd", "inv_range", "loading")
  setNames(paste0("smooth_", parts, "_", kind), parts)
}

# The marginal sd of the smooth `kind` of a model, from its variables `v`:
# the variable (a greta array, or one per draw in plain R), or for the smooth
# of the initial state, 1, its loadings per type carrying the scale
# (smooth_options(init = ))
smooth_sd <- function(v, kind) {
  if (kind == "init") {
    return(1)
  }
  v[[smooth_variable_names(kind)[["sd"]]]]
}

# The inverse range of the smooth `kind` of a model with `options`, from its
# variables `v`: the variable (a greta array, or one per draw in plain R), or
# with a fixed range, the number 1 / range
smooth_inv_range <- function(v, options, kind) {
  if (smooth_range_fixed(options$smooth)) {
    return(1 / options$smooth[["range"]])
  }
  v[[smooth_variable_names(kind)[["inv_range"]]]]
}

# The range in km of the smooth `kind` of a fit with `options`, per draw, from
# its variables `v` (a list of draws x dim arrays, or a draws matrix with a
# column per variable): estimated, or the fixed range repeated
smooth_range_km <- function(v, options, kind) {
  get <- function(name) if (is.matrix(v)) v[, name] else c(v[[name]])
  names <- smooth_variable_names(kind)
  if (smooth_range_fixed(options$smooth)) {
    n_draws <- if (is.matrix(v)) nrow(v) else
      length(v[[names[["raw"]]]]) / nrow(options$smooth$indices)
    return(rep(1000 * options$smooth[["range"]], n_draws))
  }
  1000 / get(names[["inv_range"]])
}

# The prior of the shear loading b (smooth_options(shear = TRUE)): N(1, 0.5),
# as list(mean, sd). The block screen of #47 found the log selection
# multiplier rising with the logit floor at about 0.8 per unit, across
# blocks and within them
smooth_shear_prior <- list(mean = 1, sd = 0.5)

# The weights of each smooth of a model with `options`, from its variables
# `v` (greta arrays, or draws x dim arrays in plain R), as a named list. With
# the shear, the selection smooth u_s = v_s + b u_f shares the basis, so its
# weights are those of v_s plus b times the floor's; with `own`, also
# selection_own, those of v_s alone (for the maps and checks). The weights of
# the smooth of the initial state are those of u_init, at sd 1, before the
# loadings of the types (smooth_init_loadings())
smooth_weight_terms <- function(v, options, own = FALSE) {
  out <- list()
  for (kind in smooth_kinds(options)) {
    names <- smooth_variable_names(kind)
    raw <- v[[names[["raw"]]]]
    sd <- smooth_sd(v, kind)
    inv_range <- smooth_inv_range(v, options, kind)
    if (!inherits(raw, "greta_array")) {
      n_draws <- length(raw) / nrow(options$smooth$indices)
      raw <- matrix(raw, nrow = n_draws)
      sd <- rep_len(c(sd), n_draws)
      inv_range <- rep_len(c(inv_range), n_draws)
    }
    out[[kind]] <- smooth_weights(raw, sd, inv_range, options$smooth)
  }
  if (smooth_on(options) && isTRUE(options$smooth$shear)) {
    shear <- v$smooth_shear
    if (!inherits(shear, "greta_array")) shear <- c(shear)
    if (own) out$selection_own <- out$selection
    out$selection <- out$selection + shear * out$floor
  }
  out
}


# The loadings lambda of the smooth of the initial state on each type
# (smooth_options(init = TRUE)), from the variables `v`: a greta array, one
# per type, or in plain R, draws x types (a vector for one draw); NULL
# without the smooth
smooth_init_loadings <- function(v, options) {
  if (!smooth_init_on(options)) {
    return(NULL)
  }
  v[[smooth_variable_names("init")[["loading"]]]]
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
# list named by smooth of points x 2 matrices, columns mean and sd: u_s
# ("selection"), u_f ("floor"), with the shear v_s ("selection_own"), and
# u_init ("init", at sd 1, before the types' loadings),
# from `weights` (draws x basis functions each); the basis is made `chunk`
# points at a time
smooth_posterior_at <- function(parameters, coords, chunk = 5000,
                                weights = smooth_weight_terms(
                                  parameters$variables, parameters$options,
                                  own = TRUE)) {
  smooth <- parameters$options$smooth
  out <- lapply(weights, function(w) {
    matrix(NA_real_, nrow(coords), 2, dimnames = list(NULL, c("mean", "sd")))
  })
  for (rows in split(seq_len(nrow(coords)),
                     ceiling(seq_len(nrow(coords)) / chunk))) {
    basis <- smooth_basis_at(smooth, coords[rows, , drop = FALSE])
    for (kind in names(out)) {
      u <- weights[[kind]] %*% t(basis)
      mean <- colMeans(u)
      out[[kind]][rows, "mean"] <- mean
      out[[kind]][rows, "sd"] <- sqrt(colSums(sweep(u, 2, mean) ^ 2) /
                                        (nrow(u) - 1))
    }
  }
  out
}
