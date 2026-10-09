# Check the Hilbert-space approximation of the latent smooths (V5, #47;
# R/latent_smooth.R): the covariance implied by the basis at the modelled
# bioassay cells against the exact kernel, at ranges across the prior, for
# all pairs of cells. The error relative to sd^2 does not depend on sd, so
# one sd serves. Reported as the maximum and root mean square absolute error
# over pairs, relative to sd^2, for
#   uncentred  sum_j S(omega_j) phi_j(x) phi_j(x') against k(x, x')
#   centred    the same with each basis function centred over the cells (as
#              in the model) against the exactly centred kernel, k(x, x')
#              less its means over x and over x' plus its mean over both
# at ranges from 1,000 to 16,000 km, and the fixed range if there is one; and
# at the range the basis is chosen for (smooth_basis_range()), for the full
# m[1] x m[2] grid of basis functions as well as the kept ones (smooth_box()).
#
#   Rscript R/check_latent_smooth.R ['<smooth options>']
# e.g.
#   Rscript R/check_latent_smooth.R 'smooth_options(kernel = "matern52")'
#   Rscript R/check_latent_smooth.R 'smooth_options(range = 1.5)'
#
# Plain R; about 1.5 GB and 1 minute.

arguments <- commandArgs(trailingOnly = TRUE)
smooth_expression <- if (length(arguments) >= 1) arguments[1] else
  "smooth_options()"

suppressMessages({
  sink("/dev/null")
  source("R/validation_folds.R")
  sink()
})
source("R/dynamical_model.R")

smooth <- eval(str2lang(smooth_expression))
check_smooth_options(smooth)
cells <- df$cell[match(seq_len(max(df$cell_id)), df$cell_id)]
smooth <- smooth_box(smooth, cells, classes)
coords <- smooth_cell_coords(cells, smooth$crs)
half_extent <- (apply(coords, 2, max) - apply(coords, 2, min)) / 2
rates <- smooth_prior_rates(smooth)
cat(sprintf("%s: %d modelled cells\n", smooth_expression, length(cells)))
cat(sprintf("box: origin (%.0f, %.0f) km, half-widths (%.0f, %.0f) km; the cells' half-extents (%.0f, %.0f) km, so c = (%.2f, %.2f)\n",
            1000 * smooth$origin[1], 1000 * smooth$origin[2],
            1000 * smooth$half_width[1], 1000 * smooth$half_width[2],
            1000 * half_extent[1], 1000 * half_extent[2],
            smooth$half_width[1] / half_extent[1],
            smooth$half_width[2] / half_extent[2]))
cat(sprintf("kernel %s, m = (%d, %d): %d basis functions kept of %d\n",
            smooth$kernel, smooth$m[1], smooth$m[2], nrow(smooth$indices),
            prod(smooth$m)))
basis_range <- smooth_basis_range(smooth)
cat(sprintf("basis range %.0f km; range %s\n", 1000 * basis_range,
            if (smooth_range_fixed(smooth)) {
              sprintf("fixed at %.0f km", 1000 * smooth[["range"]])
            } else {
              "estimated"
            }))
cat(sprintf("priors: 1 / range ~ Exponential(%.3f) (range in 1,000 km), %s\n",
            rates$range, smooth_sd_prior_label(smooth)))

# the exact kernel at distance d, for sd 1 and range rho
kernel <- function(d, rho) {
  ell <- rho / 2
  switch(smooth$kernel,
         matern52 = {
           r <- sqrt(5) * d / ell
           (1 + r + r ^ 2 / 3) * exp(-r)
         },
         se = exp(-d ^ 2 / (2 * ell ^ 2)))
}
centre_both <- function(k) {
  k <- sweep(k, 1, rowMeans(k))
  sweep(k, 2, colMeans(k))
}
distance <- as.matrix(dist(coords))

# the approximation errors of basis `basis` (cells x functions) with
# frequencies `omega`, at range rho
errors <- function(basis, omega, rho, centred) {
  s <- smooth_sqrt_spectral(omega, 1, 1 / rho, smooth$kernel)
  if (centred) basis <- sweep(basis, 2, colMeans(basis))
  weighted <- sweep(basis, 2, c(s), FUN = "*")
  approximate <- tcrossprod(weighted)
  exact <- kernel(distance, rho)
  if (centred) exact <- centre_both(exact)
  difference <- approximate - exact
  c(max = max(abs(difference)), rms = sqrt(mean(difference ^ 2)),
    variance = max(abs(diag(difference))))
}

basis <- hsgp_basis(coords, smooth)
omega <- hsgp_frequencies(smooth$indices, smooth$half_width)
ranges <- sort(unique(c(1, 1.5, 2, 4, 8, 16, basis_range,
                        smooth$range_prior[1], smooth[["range"]])))
rows <- list()
for (rho in ranges) {
  for (centred in c(FALSE, TRUE)) {
    e <- errors(basis, omega, rho, centred)
    rows[[length(rows) + 1]] <- data.frame(
      range_km = 1000 * rho,
      prior_quantile = exp(-rates$range / rho),
      basis = "kept", centred = centred,
      max_error = e[["max"]], rms_error = e[["rms"]],
      max_variance_error = e[["variance"]])
  }
}
# the full grid at the lower range
full <- smooth
full$indices <- as.matrix(expand.grid(x = seq_len(smooth$m[1]),
                                      y = seq_len(smooth$m[2])))
basis_full <- hsgp_basis(coords, full)
omega_full <- hsgp_frequencies(full$indices, full$half_width)
for (centred in c(FALSE, TRUE)) {
  e <- errors(basis_full, omega_full, basis_range, centred)
  rows[[length(rows) + 1]] <- data.frame(
    range_km = 1000 * basis_range,
    prior_quantile = exp(-rates$range / basis_range),
    basis = "full grid", centred = centred,
    max_error = e[["max"]], rms_error = e[["rms"]],
    max_variance_error = e[["variance"]])
}
result <- do.call(rbind, rows)
cat("\nerrors of the implied covariance over all pairs of modelled cells,",
    "relative to sd^2:\n")
print(result, digits = 3, row.names = FALSE)

