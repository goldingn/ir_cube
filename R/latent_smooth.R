# The latent spatial smooths and the mortality floor of the dynamical model
# (#47). Two static surfaces, fitted inside the model, for the classes
# smooth_classes (the pyrethroids and DDT, to both of which kdr gives
# resistance):
#   u_s(x)  the log multiplier of the strength of selection: the cumulative
#           log fitness of a cell is multiplied by exp(u_s(x))
#           (closed_form_states(), R/dynamical_model.R)
#   u_f(x)  the log hazard of the mortality floor, 1 - exp(-lambda(x))
# Each is a low-rank Hilbert-space approximate Gaussian process (HSGP; Solin
# and Sarkka 2020, Riutort-Mayol et al. 2022) with the squared exponential
# kernel on a box around the modelled bioassay cells, in the Albers
# equal-area projection of Africa of the two-stage correction, in units of
# 1,000 km:
#   u(x) = sum_j sqrt(S(omega_j; sd, rho)) beta_j phi_j(x),  beta_j ~ N(0, 1)
# with phi_j the eigenfunctions of the Laplacian on the box (Dirichlet
# boundary), omega_j the square roots of their eigenvalues, and S the
# kernel's spectral density. The basis at the cells is data; only the
# weights depend on sd. Functions and settings only; needs sf and terra.
# Sourced by R/dynamical_model.R.


# settings ---------------------------------------------------------------------

# The classes with the smooths. The other classes' selection is set by the
# covariates alone, and their floors are flat
smooth_classes <- c("Pyrethroids", "Organochlorines")

# The smooths' range rho, in 1,000 km, fixed: the spatial correlation range of
# recent pyrethroid mortality is 1,230 km [840, 1,800], and of block plateaus
# 1,500 km [1,100, 2,000] (#47). The lengthscale is rho / 2
smooth_range <- 1.5

# The box's half-widths as multiples of the half-extents of the modelled
# bioassay cells, widened where needed to cover every cell of the mask
smooth_box_factor <- 2

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

# the projection of the coordinates: the two-stage correction's Albers
# equal-area conic for Africa (africa_equal_area_crs,
# R/two_stage_correction.R), in km
smooth_crs <- paste(
  "+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25",
  "+x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs"
)


# the basis --------------------------------------------------------------------

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

# The basis of the smooths, for the modelled bioassay cells `cells` (mask
# cells, each once) and those of them with bioassays of smooth_classes,
# `class_cells`; `classes` are the class names in class_id order and
# `classes_index` the class of each type. build_dynamical_model() records it
# in the model's options (options$smooth), for the plain-R predictions:
#   crs, origin, half_width  the box: centred on the cells' midrange, its
#                       half-width in each dimension smooth_box_factor times
#                       their half-extent, or more where the mask reaches
#                       beyond that, so that every mask cell can be predicted
#   indices             the frequency indices of the basis functions, ordered
#                       by frequency. m per dimension is smooth_m_factor L /
#                       ell, rounded up, at the fixed range's lengthscale
#                       ((17, 14) for the October 2026 data); of the m[1] x
#                       m[2] products of the one-dimensional eigenfunctions,
#                       only those with frequency up to the lower of the two
#                       dimensions' highest are kept (164), as the spectral
#                       density is isotropic. The covariance the basis implies
#                       at the modelled cells is within 2.3% of the kernel's
#                       at the fixed range
#   class_column_means  the basis functions' means over class_cells, which
#                       centre the selection smooth: it has mean 0 there, so it
#                       does not trade off with the overall strength of
#                       selection
#   floor_mean_cells    class_cells, sorted, over which the floor smooth's
#                       mean hazard is taken (floor_log_hazard())
#   term_classes        whether each class has the smooths
#   type_classes        the class of each type
smooth_basis <- function(cells, class_cells, classes, classes_index) {
  stopifnot(length(class_cells) > 0, all(class_cells %in% cells))
  smooth <- list(crs = smooth_crs)
  coords <- smooth_cell_coords(cells, smooth$crs)
  lower <- apply(coords, 2, min)
  upper <- apply(coords, 2, max)
  smooth$origin <- (lower + upper) / 2
  mask_range <- smooth_mask_range(smooth$crs)
  reach <- pmax(abs(mask_range[1, ] - smooth$origin),
                abs(mask_range[2, ] - smooth$origin))
  smooth$half_width <- pmax(smooth_box_factor * (upper - lower) / 2, reach)
  m <- ceiling(smooth_m_factor * smooth$half_width / (smooth_range / 2))
  indices <- as.matrix(expand.grid(x = seq_len(m[1]), y = seq_len(m[2])))
  omega <- hsgp_frequencies(indices, smooth$half_width)
  omega_max <- min(pi * m / (2 * smooth$half_width))
  keep <- omega <= omega_max * (1 + 1e-12)
  smooth$indices <- indices[keep, , drop = FALSE][order(omega[keep]), ,
                                                  drop = FALSE]
  smooth$class_column_means <- colMeans(hsgp_basis(
    smooth_cell_coords(class_cells, smooth$crs), smooth))
  smooth$floor_mean_cells <- sort(class_cells)
  smooth$term_classes <- classes %in% smooth_classes
  smooth$type_classes <- as.integer(classes_index)
  smooth
}

# The frequencies omega_j = sqrt(lambda_j) of the basis functions with
# indices `indices` (one row per function, one column per dimension) on a box
# of half-widths `half_width`: lambda_j = sum_d (pi j_d / (2 L_d))^2
hsgp_frequencies <- function(indices, half_width) {
  sqrt(rowSums(sweep(indices, 2, pi / (2 * half_width), FUN = "*") ^ 2))
}

# The basis functions at projected coordinates `coords` (points x 2, in 1,000
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

# The basis that the smooths' weights (smooth_weights()) multiply, at mask
# cells `cells`, for a model with `options` (options$smooth, set by
# build_dynamical_model()): the basis functions and a column of ones, as a
# cells x (basis functions + 1) matrix
prediction_basis <- function(options, cells) {
  smooth <- options$smooth
  if (is.null(smooth$indices)) {
    stop("the options have no basis; build_dynamical_model() adds it")
  }
  cbind(hsgp_basis(smooth_cell_coords(cells, smooth$crs), smooth), 1)
}


# the smooths and the floor ----------------------------------------------------

# The square root of the squared exponential's spectral density in two
# dimensions at angular frequencies `omega`, at marginal sd 1 and the fixed
# range (Rasmussen and Williams 2006, section 4.2, in angular frequency, with
# k(r) = (2 pi)^-2 int S(w) exp(i w.r) dw):
#   S(w) = sd^2 2 pi ell^2 exp(-ell^2 w^2 / 2),  ell = rho / 2
smooth_sqrt_spectral <- function(omega) {
  ell <- smooth_range / 2
  ell * sqrt(2 * pi) * exp(-ell ^ 2 * omega ^ 2 / 4)
}

# The weights of the selection and floor smooths, from the model's variables
# `v` (greta arrays, or draws x dim arrays in plain R), as list(selection,
# floor), for the basis of prediction_basis(): the weights sqrt(S(omega_j))
# beta_j of the basis functions, then a centring term, minus their product
# with the smooth's centre, which multiplies the basis's column of ones. The
# selection smooth is centred over the cells with bioassays of smooth_classes
# (class_column_means); the floor smooth is not centred, its level cancelling
# in floor_log_hazard(). A greta array is (basis functions + 1) x 1, and in
# plain R the weights are draws x (basis functions + 1)
smooth_weights <- function(v, smooth) {
  spectral <- smooth_sqrt_spectral(hsgp_frequencies(smooth$indices,
                                                    smooth$half_width))
  centres <- list(selection = smooth$class_column_means,
                  floor = numeric(length(spectral)))
  out <- list()
  for (kind in names(centres)) {
    raw <- v[[paste0("smooth_raw_", kind)]]
    sd <- v[[paste0("smooth_sd_", kind)]]
    centre <- centres[[kind]]
    if (inherits(raw, "greta_array")) {
      weights <- raw * (sd * spectral)
      out[[kind]] <- rbind(weights, -(matrix(centre, nrow = 1) %*% weights))
    } else {
      n_draws <- length(sd)
      weights <- matrix(raw, nrow = n_draws) *
        outer(c(sd), spectral)
      out[[kind]] <- cbind(weights, -(weights %*% centre))
    }
  }
  out
}

# The log of each type's floor hazard where the floor smooth is 0, as a greta
# array (types x 1), or in plain R draws x types, from the variables `v` and
# the floor smooth's weights (smooth_weights()). Each type has a floor
#   f(x) = 1 - exp(-lambda(x)),  log lambda(x) = log hazard + u_f(x)
# with u_f only for the types of smooth_classes, whose mean hazard over the
# cells floor_mean_cells is lambda_bar, the type's:
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
floor_log_hazard <- function(v, smooth, floor_weights) {
  classes <- smooth$type_classes
  same <- outer(classes, classes, "==")
  averaging <- same / rowSums(same)
  has_smooth <- as.numeric(smooth$term_classes[classes])
  basis <- prediction_basis(list(smooth = smooth), smooth$floor_mean_cells)
  if (inherits(v$floor_flat, "greta_array")) {
    d <- v$floor_sd_type * v$floor_type_raw
    log_lambda_bar <- log(v$floor_flat[classes]) + d -
      log(averaging %*% exp(d))
    u <- basis_product(basis, floor_weights)
    log_mean <- log(sum(exp(u))) - log(nrow(basis))
    return(log_lambda_bar - log_mean * has_smooth)
  }
  n_draws <- length(v$floor_sd_type)
  flat <- matrix(v$floor_flat, n_draws)
  d <- c(v$floor_sd_type) * matrix(v$floor_type_raw, n_draws)
  log_lambda_bar <- log(flat[, classes, drop = FALSE]) + d -
    log(exp(d) %*% t(averaging))
  u <- floor_weights %*% t(basis)
  u_max <- apply(u, 1, max)
  log_mean <- u_max + log(rowMeans(exp(u - u_max)))
  log_lambda_bar - outer(log_mean, has_smooth)
}

# The floor 1 - exp(-exp(log_hazard + u)), for greta arrays (one element per
# row) or plain R (log_hazard one per draw, u draws x cells, giving draws x
# cells; or one per draw without u)
floor_value <- function(log_hazard, u = NULL) {
  eta <- if (is.null(u)) log_hazard else log_hazard + u
  if (inherits(eta, "greta_array")) icloglog(eta) else -expm1(-exp(eta))
}

# The product basis %*% weights of a data matrix `basis` (n x m) and a greta
# array `weights` (m x 1), as an n x 1 greta array, by one greta op. greta's
# %*% would copy the data matrix to every chain's slice of the batch at every
# evaluation of the density; the einsum of tf_basis_product() broadcasts it
# instead
basis_product <- function(basis, weights) {
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
