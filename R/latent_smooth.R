# The latent spatial smooths and the mortality floor of the dynamical model
# (#46). Two static surfaces, fitted inside the model, for the classes
# smooth_classes (the pyrethroids and DDT, to both of which kdr gives
# resistance):
#   u_s(x)  the log multiplier of the strength of selection: the cumulative
#           log fitness of a cell is multiplied by exp(u_s(x))
#           (closed_form_states(), R/dynamical_model.R)
#   u_f(x)  the log hazard of the mortality floor, 1 - exp(-lambda(x))
# Each is a low-rank Hilbert-space approximate Gaussian process (HSGP; Solin
# and Sarkka 2020, Riutort-Mayol et al. 2022) with the squared exponential
# kernel on a box around Africa, in the Albers equal-area projection of the
# two-stage correction (R/projection.R), in units of 1,000 km:
#   u(x) = sum_j sqrt(S(omega_j; sd, rho)) beta_j phi_j(x),  beta_j ~ N(0, 1)
# with phi_j the eigenfunctions of the Laplacian on the box (Dirichlet
# boundary), omega_j the square roots of their eigenvalues, and S the
# kernel's spectral density. The basis at the cells is data; only the
# weights depend on sd. Functions and settings only; needs sf and terra.
# Sourced by R/dynamical_model.R.

source("R/projection.R")


# settings ---------------------------------------------------------------------

# The classes with the smooths. The other classes' selection is set by the
# covariates alone, and their floors are flat
smooth_classes <- c("Pyrethroids", "Organochlorines")

# The smooths' range rho, in 1,000 km, fixed: the spatial correlation range of
# recent pyrethroid mortality is 1,230 km [840, 1,800], and of block plateaus
# 1,500 km [1,100, 2,000] (#47). The lengthscale is rho / 2
smooth_range <- 1.5

# The box's half-widths as multiples of the mask's half-extents. At 1.5 the
# covariance the basis implies at the modelled bioassay cells is within 2.2%
# of the kernel's (of sd^2, centred or not), as at 1.75 and 2; at 1.3 within
# 4%, and at 1.2 within 18%
smooth_box_factor <- 1.5

# The number of basis functions per dimension needed for an accurate
# approximation of the squared exponential at lengthscale ell, as a multiple
# of L / ell for box half-width L (Riutort-Mayol et al. 2022)
smooth_m_factor <- 1.75

# The scales of the half-normal priors, tight as complexity penalties: looser
# ones let the floors rise to 0.4-0.6 and the smooths' sds to 1-3, and the
# chains split into modes (#47)
#   smooth_sd_scale        each smooth's marginal sd
#   floor_hazard_scale     each class's mean floor hazard, over the cells with
#                          bioassays of smooth_classes for the classes with
#                          the floor smooth, and its flat hazard for the
#                          others: P(> 0.01) = 0.046, so the floors stay near
#                          0 unless the data demand them, and resistance
#                          rises by selection
#   floor_type_sd_scale    the sd of a type's log mean hazard about its
#                          class's: P(two types of a class have floors more
#                          than 3 times apart) = 0.05
smooth_sd_scale <- 0.025
floor_hazard_scale <- 0.005
floor_type_sd_scale <- 0.356


# the basis --------------------------------------------------------------------

# Projected coordinates, in 1,000 km, of mask cells `cells` (cell numbers of
# data/clean/raster_mask.tif)
smooth_cell_coords <- function(cells) {
  xy <- terra::xyFromCell(terra::rast("data/clean/raster_mask.tif"), cells)
  project_km(xy[, 1], xy[, 2]) / 1000
}

# The basis of the smooths, for the mask cells with bioassays of
# smooth_classes, `class_cells`, and `classes`, the class names in class_id
# order. build_dynamical_model() records it in the model's options
# (options$smooth), for the plain-R predictions:
#   origin, half_width  the box: centred on the mask's midrange, its
#                       half-width in each dimension smooth_box_factor times
#                       the mask's half-extent (prediction_mask_coords()), so
#                       the same for every fold and data update
#   range               the smooths' range, smooth_range
#   indices             the frequency indices of the basis functions, ordered
#                       by frequency. m per dimension is smooth_m_factor L /
#                       ell, rounded up, at the range's lengthscale (13, 15);
#                       of the m[1] x m[2] products of the one-dimensional
#                       eigenfunctions, only those with frequency up to the
#                       lower of the two dimensions' highest are kept (134),
#                       as the spectral density is isotropic
#   class_column_means  the basis functions' means over class_cells, which
#                       hsgp_basis() then subtracts, centring both smooths
#                       there: the selection smooth has mean 0 there, so it
#                       does not trade off with the overall strength of
#                       selection, and the floor smooth's level cancels in
#                       floor_log_hazard() anyway
#   floor_mean_basis    the basis at class_cells, sorted, over which the floor
#                       smooth's mean hazard is taken (floor_log_hazard())
#   term_classes        whether each class has the smooths
smooth_basis <- function(class_cells, classes) {
  stopifnot(length(class_cells) > 0)
  mask <- prediction_mask_coords() / 1000
  lower <- apply(mask, 2, min)
  upper <- apply(mask, 2, max)
  smooth <- list(origin = (lower + upper) / 2,
                 half_width = smooth_box_factor * (upper - lower) / 2,
                 range = smooth_range)
  m <- ceiling(smooth_m_factor * smooth$half_width / (smooth$range / 2))
  indices <- as.matrix(expand.grid(x = seq_len(m[1]), y = seq_len(m[2])))
  omega <- hsgp_frequencies(indices, smooth$half_width)
  omega_max <- min(pi * m / (2 * smooth$half_width))
  keep <- omega <= omega_max * (1 + 1e-12)
  smooth$indices <- indices[keep, , drop = FALSE][order(omega[keep]), ,
                                                  drop = FALSE]
  smooth$class_column_means <- colMeans(prediction_basis(smooth, class_cells))
  smooth$floor_mean_basis <- prediction_basis(smooth, sort(class_cells))
  smooth$term_classes <- classes %in% smooth_classes
  smooth
}

# Labels of the basis functions of `smooth` (smooth_basis()), by which
# dynamical_inits() matches cached weights to them: their frequency indices,
# with the box (origin and half-widths) and the range to 4 significant
# figures, so that the weights of another basis match none of them
smooth_basis_labels <- function(smooth) {
  box <- paste(signif(c(smooth$origin, smooth$half_width), 4), collapse = ",")
  sprintf("(%d,%d) box %s range %s", smooth$indices[, 1], smooth$indices[, 2],
          box, signif(smooth$range, 4))
}

# The frequencies omega_j = sqrt(lambda_j) of the basis functions with
# indices `indices` (one row per function, one column per dimension) on a box
# of half-widths `half_width`: lambda_j = sum_d (pi j_d / (2 L_d))^2
hsgp_frequencies <- function(indices, half_width) {
  sqrt(rowSums(sweep(indices, 2, pi / (2 * half_width), FUN = "*") ^ 2))
}

# The basis functions at projected coordinates `coords` (points x 2, in 1,000
# km), as a points x basis functions matrix: phi_j(x) = prod_d sqrt(1 / L_d)
# sin(pi j_d (x_d + L_d) / (2 L_d)), x relative to the box's origin, less
# their means over the cells with bioassays of smooth_classes once
# smooth_basis() has them (class_column_means)
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
  if (!is.null(smooth$class_column_means)) {
    phi <- sweep(phi, 2, smooth$class_column_means)
  }
  phi
}

# The basis functions at mask cells `cells`, for the basis `smooth`
# (options$smooth), as a cells x basis functions matrix
prediction_basis <- function(smooth, cells) {
  hsgp_basis(smooth_cell_coords(cells), smooth)
}


# the smooths and the floor ----------------------------------------------------

# The square root of the squared exponential's spectral density in two
# dimensions at angular frequencies `omega`, at marginal sd 1 and range
# `range` (Rasmussen and Williams 2006, section 4.2, in angular frequency,
# with k(r) = (2 pi)^-2 int S(w) exp(i w.r) dw):
#   S(w) = sd^2 2 pi ell^2 exp(-ell^2 w^2 / 2),  ell = rho / 2
smooth_sqrt_spectral <- function(omega, range) {
  ell <- range / 2
  ell * sqrt(2 * pi) * exp(-ell ^ 2 * omega ^ 2 / 4)
}

# The terms of the latent smooths and the floor, from the variables `v`, as
# greta arrays or one draw in plain R (vectors as vectors), for the basis
# `smooth` (options$smooth), and the class of each type, classes_index:
#   smooth_weights_selection, smooth_weights_floor
#                       the weights sqrt(S(omega_j)) beta_j of the basis
#                       functions, so that a smooth at points with basis phi
#                       (hsgp_basis()) is phi %*% weights
#   floor_log_hazard    the log floor hazard of each type where the floor
#                       smooth is 0 (floor_log_hazard())
smooth_terms <- function(v, smooth, classes_index) {
  spectral <- smooth_sqrt_spectral(hsgp_frequencies(smooth$indices,
                                                    smooth$half_width),
                                   smooth$range)
  weights_selection <- v$smooth_raw_selection * (v$smooth_sd_selection *
                                                    spectral)
  weights_floor <- v$smooth_raw_floor * (v$smooth_sd_floor * spectral)
  list(smooth_weights_selection = weights_selection,
       smooth_weights_floor = weights_floor,
       floor_log_hazard = floor_log_hazard(v, smooth, classes_index,
                                           weights_floor))
}

# The log of each type's floor hazard where the floor smooth is 0, from the
# variables `v` (greta arrays, or one draw in plain R) and the floor smooth's
# weights. Each type has a floor
#   f(x) = 1 - exp(-lambda(x)),  log lambda(x) = log hazard + u_f(x)
# with u_f only for the types of smooth_classes, whose mean hazard over the
# cells of floor_mean_basis is lambda_bar, the type's:
#   lambda(x) = lambda_bar exp(u_f(x)) / mean_c exp(u_f(c))
# so that the floor smooth's level cancels and the prior of the class's mean
# hazard is on the realised mean over those cells, whatever the smooth does;
# the other types' floors are flat at lambda_bar. The types' mean hazards
# are pooled within their class by the arithmetic mean:
#   lambda_bar = floor_flat[class] exp(d) / mean over the class's types of
#     exp(d),  d = floor_sd_type floor_type_raw
# so floor_flat is the arithmetic mean of its types' mean hazards, no type
# exceeds (number of types) floor_flat, and lowering one raises the others; a
# class of one type has its floor_flat. Pooled on the log scale instead, one
# type with a floor near 0 could hold the class's mean low while its other
# types rose
floor_log_hazard <- function(v, smooth, classes_index, floor_weights) {
  same <- outer(classes_index, classes_index, "==")
  averaging <- same / rowSums(same)
  d <- v$floor_sd_type * v$floor_type_raw
  log_lambda_bar <- log(v$floor_flat[classes_index]) + d -
    log(averaging %*% exp(d))
  u <- basis_product(smooth$floor_mean_basis, floor_weights)
  log_mean <- log(sum(exp(u))) - log(nrow(smooth$floor_mean_basis))
  log_lambda_bar - log_mean * as.numeric(smooth$term_classes[classes_index])
}

# The floor 1 - exp(-exp(eta)) at the log hazard eta, for greta arrays or
# plain R
floor_value <- function(eta) {
  if (inherits(eta, "greta_array")) icloglog(eta) else -expm1(-exp(eta))
}

# The product basis %*% weights of a data matrix `basis` (n x m) and weights
# (m x 1). For a greta array, as an n x 1 greta array by one greta op: greta's
# %*% would copy the data matrix to every chain's slice of the batch at every
# evaluation of the density, where the einsum of tf_basis_product()
# broadcasts it
basis_product <- function(basis, weights) {
  if (!inherits(weights, "greta_array")) {
    return(basis %*% weights)
  }
  op_env <- new.env(parent = globalenv())
  op_env$tf_basis_product <- tf_basis_product
  greta:::op("basis_product", weights,
             operation_args = list(basis = basis),
             tf_operation = "tf_basis_product",
             tf_function_env = op_env,
             dim = c(nrow(basis), 1L))
}

# The TensorFlow side of basis_product(): weights (B, m, 1), basis an R
# matrix (n, m); returns (B, n, 1)
tf_basis_product <- function(weights, basis) {
  tf <- tensorflow::tf
  basis <- tf$constant(basis, dtype = weights$dtype)
  tf$einsum("nm,bmk->bnk", basis, weights)
}
