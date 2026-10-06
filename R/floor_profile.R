# The profile and Laplace approximation of the log posterior density of a
# saved full fit's mortality floor(s) (#14, #37, #47), to weigh the floor
# modes against each other.
#
#   Rscript R/floor_profile.R <fitted_model.RData> <label> [<floor> [<points>]]
#   Rscript R/floor_profile.R <fitted_model.RData> <label> summarise
#
# For each value f of a grid of the floor (PROFILE_GRID, comma-separated, or
# profile_grid below: 15 values, denser near 0), the floor is fixed at f and
# the log posterior density is maximised over every other free parameter, on
# greta's unconstrained scale and with the Jacobians of its transforms (the
# density of the free state, which HMC samples), by L-BFGS-B with
# TensorFlow's gradients, from the mean free state of each floor mode's
# usable chains (all but those stuck), then from the better optimum by Newton
# steps with the Hessian below until the largest gradient is under 0.01. The
# country levels are non-centred for this (build_dynamical_model(countries =
# "noncentred")), the same posterior in other coordinates: centred, the
# density is unbounded where a type's country levels can all meet their means
# as its init_country_sd goes to 0, and the optimiser runs off along it. With
# the species model's two floors, each is profiled in turn, the other free. At
# the optimum, H, the Hessian of the negative log density in the other free
# parameters, by central differences of the gradient, gives
#   profile height      the maximum, log p(u, theta*), u = logit f
#   -1/2 log |H|        the curvature correction
#   Laplace             their sum: the log marginal density of u, up to a
#                       constant, log p(u | y) ~ log p(u, theta*) + d/2 log 2
#                       pi - 1/2 log |H|
# The Laplace approximation depends on the coordinates; these are the
# non-centred ones. A non-positive-definite H (an optimum that is not a
# maximum, or one that has not converged) is reported, and has no Laplace
# value. Each grid point is saved as it is done, to
# outputs/species_runs/floor_profile/<label>/, and points already there are
# skipped, so a run can be resumed, and points split between processes:
# <floor> (e.g. mortality_floor, other_floor) and <points> (grid indices, e.g.
# 1-5 or 3,7) choose them. A point takes about 2 minutes (two L-BFGS-B runs of
# 10-20 s, and 2-3 Hessians of 30 s) and a process about 6.5 GB.
#
# `summarise` (plain R) combines the points, plots the profile height and the
# Laplace curve against the floor, and estimates the relative posterior mass
# of each floor mode: the Laplace density of u, interpolated by a natural
# spline over the grid, integrated over each bump (split at the minimum between
# the two highest maxima); mass beyond the grid is not counted. Writes
# outputs/species_runs/floor_profile/<label>.csv and
# figures/species_runs/floor_profile_<label>.png.
#
# Run with the greta 0.6 environment (doc/cv_run_plan.md); PROFILE_THREADS
# (default 4) TensorFlow threads; PROFILE_MAXIT (default 5000) the most
# L-BFGS-B iterations; PROFILE_NEWTON (default 10) the most Newton steps;
# FLOOR_PRIOR for fits saved without a floor prior.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) >= 2)
file <- arguments[1]
label <- arguments[2]
summarise_only <- length(arguments) >= 3 && arguments[3] == "summarise"

profile_grid <- c(0.0005, 0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.075,
                  0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5)
if (nzchar(Sys.getenv("PROFILE_GRID"))) {
  profile_grid <- as.numeric(strsplit(Sys.getenv("PROFILE_GRID"), ",")[[1]])
}
stopifnot(all(profile_grid > 0 & profile_grid < 1))
points_dir <- file.path("outputs/species_runs/floor_profile", label)
figure_dir <- "figures/species_runs"
dir.create(points_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
point_file <- function(name, i) {
  file.path(points_dir, sprintf("%s_%02d.rds",
                                gsub("[^A-Za-z0-9_]+", "_", name), i))
}


# summarise ------------------------------------------------------------------------

# The relative posterior mass of the bumps of a log density `l` on grid `u`:
# a natural spline through the points, on a fine grid, split at the minimum
# between the two highest local maxima. NA where fewer than 4 points are
# finite
mode_masses <- function(u, l) {
  ok <- is.finite(l)
  if (sum(ok) < 4) {
    return(tibble::tibble(bump = NA_character_, from = NA_real_,
                          to = NA_real_, mass = NA_real_))
  }
  spline <- splinefun(u[ok], l[ok], method = "natural")
  fine <- seq(min(u[ok]), max(u[ok]), length.out = 2000)
  value <- spline(fine)
  density <- exp(value - max(value))
  maxima <- which(diff(sign(diff(value))) < 0) + 1
  if (value[1] > value[2]) maxima <- c(1, maxima)
  if (value[length(value)] > value[length(value) - 1]) {
    maxima <- c(maxima, length(value))
  }
  top <- sort(maxima[order(value[maxima], decreasing = TRUE)][
    seq_len(min(2, length(maxima)))])
  cuts <- c(1, if (length(top) == 2) {
    top[1] - 1 + which.min(value[top[1]:top[2]])
  }, length(fine))
  trapezoid <- function(i) {
    sum(diff(fine[i]) * (head(density[i], -1) + tail(density[i], -1)) / 2)
  }
  masses <- vapply(seq_len(length(cuts) - 1), function(b) {
    trapezoid(cuts[b]:cuts[b + 1])
  }, numeric(1))
  tibble::tibble(bump = as.character(seq_along(masses)),
                 from = plogis(fine[head(cuts, -1)]),
                 to = plogis(fine[tail(cuts, -1)]),
                 peak = plogis(fine[top]),
                 mass = masses / sum(masses))
}

if (summarise_only) {
  suppressMessages({
    library(dplyr)
    library(ggplot2)
  })
  files <- list.files(points_dir, pattern = "\\.rds$", full.names = TRUE)
  stopifnot(length(files) > 0)
  points <- bind_rows(lapply(files, function(f) {
    point <- readRDS(f)
    as_tibble(point[setdiff(names(point), c("free", "starts",
                                            "eigenvalues"))])
  })) %>%
    arrange(floor, value)
  write.csv(points, file.path(dirname(points_dir), sprintf("%s.csv", label)),
            row.names = FALSE)
  options(width = 160)
  print(as.data.frame(points %>%
                        select(floor, value, height, half_log_det, laplace,
                               n_nonpositive, best_start, max_gradient,
                               converged, seconds)), digits = 6)
  masses <- points %>%
    group_by(floor) %>%
    group_modify(~ mode_masses(.x$u, .x$laplace)) %>%
    ungroup()
  cat("\nrelative posterior mass of each bump of the Laplace curve (within the grid):\n")
  print(as.data.frame(masses), digits = 4)
  write.csv(masses, file.path(dirname(points_dir),
                              sprintf("%s_masses.csv", label)),
            row.names = FALSE)
  curves <- points %>%
    group_by(floor) %>%
    mutate(`profile height` = height - max(height),
           Laplace = laplace - max(laplace, na.rm = TRUE)) %>%
    ungroup() %>%
    tidyr::pivot_longer(c(`profile height`, Laplace), names_to = "curve",
                        values_to = "relative")
  p <- ggplot(curves, aes(value, relative, colour = curve)) +
    geom_line() +
    geom_point(aes(shape = n_nonpositive > 0), size = 1.8) +
    facet_wrap(~ floor) +
    scale_x_log10() +
    coord_cartesian(ylim = c(-60, 2)) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                       name = "H not positive definite") +
    labs(x = "floor (log scale)",
         y = "log density of logit(floor), relative to the maximum",
         title = sprintf("%s: profile and Laplace approximation of the floor%s",
                         label, if (n_distinct(points$floor) > 1) "s" else ""),
         caption = paste("each point: the other free parameters maximised",
                         "with the floor fixed; Laplace adds -1/2 log|H|;",
                         "cut at -60")) +
    theme_bw(base_size = 9)
  ggsave(file.path(figure_dir, sprintf("floor_profile_%s.png", label)), p,
         width = 8, height = 4, dpi = 150)
  quit(save = "no")
}


# the grid points ---------------------------------------------------------------------

threads <- as.integer(Sys.getenv("PROFILE_THREADS", "4"))
maxit <- as.integer(Sys.getenv("PROFILE_MAXIT", "5000"))
newton_steps <- as.integer(Sys.getenv("PROFILE_NEWTON", "10"))
source("R/greta_setup.R")
start_greta(threads = threads)
suppressMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
})
source("R/functions.R")
source("R/dynamical_predictions.R")
source("R/species_fit_helpers.R")

fit <- load_fit(file)
floors <- fit_floor_names(fit$draws)
if (length(floors) == 0) {
  report("%s has no floor; nothing to profile", label)
  quit(save = "no")
}
# the fit's model, to check and read its free states, and the same posterior
# with non-centred country levels (build_dynamical_model(countries =
# "noncentred")), in which the density is maximised: centred, it is unbounded
# where a type's country levels can all meet their means as its
# init_country_sd goes to 0 (in r2_main, Pirimiphos-methyl's went to 5e-5)
built <- rebuild_fit_model(fit)
noncentred <- build_dynamical_model(train_df = fit$df,
                                    df = fit$df,
                                    x_cell_years = fit$x_cell_years,
                                    cell_years_index = fit$cell_years_index,
                                    classes_index = fit$classes_index,
                                    types = fit$types,
                                    options = fit$options,
                                    x_cells_init = fit$x_cells_init,
                                    countries = "noncentred")
value_gradient <- log_density_gradient_function(noncentred$model)

# Free states of the fit's model (centred, rows of `free`) in the non-centred
# model: every variable but the country levels as it is, and
# init_country_raw = (level - mean) / sd (country_level_prior()). Checked:
# the non-centred density is the centred one times prod sd, over countries and
# types
to_noncentred <- function(free) {
  free <- matrix(free, ncol = length(built$free_order))
  centred_columns <- free_state_columns(built$model)
  columns <- free_state_columns(noncentred$model)
  out <- matrix(NA_real_, nrow(free), sum(lengths(columns)))
  for (name in names(attr(columns, "targets"))) {
    if (name == "init_country_raw") next
    out[, columns[[attr(columns, "targets")[[name]]]]] <-
      free[, centred_columns[[attr(centred_columns, "targets")[[name]]]]]
  }
  values <- built$model$dag$trace_values(free)
  raw_column <- columns[[attr(columns, "targets")[["init_country_raw"]]]]
  log_jacobian <- numeric(nrow(free))
  raw <- list()
  for (i in seq_len(nrow(free))) {
    v <- lapply(setNames(nm = c("init_region_raw", "init_region_sd",
                                "logit_init_mean", "init_country_sd",
                                "init_country_level",
                                if (!is.null(fit$options$init_covariates)) {
                                  "init_coef"
                                })),
                function(name) {
                  x <- extract_parameter(values[i, , drop = FALSE], name)
                  array(x, dim(x)[-1])
                })
    prior <- country_level_prior(v, built$lookups$country_region_index,
                                 built$options)
    raw[[i]] <- (v$init_country_level - prior$mean) / prior$sd
    # greta lays a matrix variable out in the free state by rows
    out[i, raw_column] <- c(t(raw[[i]]))
    log_jacobian[i] <- sum(log(prior$sd))
  }
  check <- noncentred$model$dag$trace_values(out)
  raw_saved <- check[, grep("^init_country_raw\\[", colnames(check)),
                     drop = FALSE]
  stopifnot(max(abs(raw_saved - t(vapply(raw, c, numeric(length(raw_column)))))) <
              1e-12)
  density <- log_density_function(built$model)(free)[, "adjusted"]
  density_noncentred <- log_density_function(noncentred$model)(out)[,
                                                                 "adjusted"]
  difference <- max(abs(density_noncentred - density - log_jacobian))
  report("centred to non-centred free states: log density difference less log prod sd, max abs %.3g",
         difference)
  stopifnot(difference < 1e-6)
  out
}
usable <- usable_chains(fit$draws, label)
mode_of_chain <- chain_floor_mode(fit$draws)
modes <- split(usable, mode_of_chain[usable])
# each mode's mean free state, in the non-centred model
starts <- lapply(modes, function(chains) {
  c(to_noncentred(colMeans(free_states(fit, built, chains))))
})
report("%s: floors %s; starts from the mean of modes %s", label,
       toString(floors), paste(sprintf("%s (chains %s)", names(modes),
                                       vapply(modes, toString, "")),
                               collapse = ", "))

profile_floors <- if (length(arguments) >= 3) arguments[3] else floors
stopifnot(all(profile_floors %in% floors))
indices <- seq_along(profile_grid)
if (length(arguments) >= 4) {
  indices <- unlist(lapply(strsplit(arguments[4], ",")[[1]], function(part) {
    range <- as.integer(strsplit(part, "-")[[1]])
    seq(range[1], range[length(range)])
  }))
}
stopifnot(all(indices %in% seq_along(profile_grid)))

# Maximise the adjusted log density over every free parameter but `fixed`,
# from `start`, by L-BFGS-B. A state where the density is not finite (p
# rounding to 0 or 1 where a bioassay has both outcomes) gets a very large
# negative log density, so the line search steps back from it. (Whitening by
# the posterior covariance of the draws made it slower.)
optimise_others <- function(start, fixed) {
  others <- setdiff(seq_along(start), fixed)
  full <- start
  cache <- new.env()
  evaluate <- function(theta) {
    if (!identical(theta, cache$theta)) {
      full[others] <- theta
      result <- value_gradient(full)
      cache$theta <- theta
      if (is.finite(result$value) && all(is.finite(result$gradient))) {
        cache$value <- -result$value
        cache$gradient <- -result$gradient[1, others]
      } else {
        cache$value <- 1e100
        if (is.null(cache$gradient)) cache$gradient <- rep(0, length(others))
      }
    }
  }
  seconds <- system.time(
    opt <- optim(start[others],
                 fn = function(theta) { evaluate(theta); cache$value },
                 gr = function(theta) { evaluate(theta); cache$gradient },
                 method = "L-BFGS-B",
                 control = list(maxit = maxit, factr = 1e5, lmm = 20))
  )[["elapsed"]]
  full[others] <- opt$par
  result <- value_gradient(full)
  list(free = full, value = result$value,
       max_gradient = max(abs(result$gradient[1, others])),
       evaluations = opt$counts[["function"]],
       converged = opt$convergence == 0, message = opt$message,
       seconds = seconds, start_value = value_gradient(start)$value)
}

# Newton steps from `free` (the best L-BFGS-B optimum) in every free
# parameter but `fixed`, with the Hessian of hessian_others(), each step
# solving with the Hessian's eigenvalues floored at `min_eigenvalue` (so it is
# an ascent direction where the Hessian is not positive definite) and halved
# until the density increases; until the largest gradient is below
# `tolerance` or after `max_steps`. Returns the last state, its density and
# gradient, and the Hessian there, which the Laplace approximation uses
newton_others <- function(free, fixed, max_steps = newton_steps,
                          tolerance = 0.01, min_eigenvalue = 1e-3) {
  others <- setdiff(seq_along(free), fixed)
  current <- value_gradient(free)
  steps <- 0
  repeat {
    hessian <- hessian_others(free, fixed)
    gradient <- current$gradient[1, others]
    if (max(abs(gradient)) < tolerance || steps >= max_steps) break
    e <- eigen(hessian$H, symmetric = TRUE)
    direction <- c(e$vectors %*% (crossprod(e$vectors, gradient) /
                                    pmax(abs(e$values), min_eigenvalue)))
    step <- 1
    repeat {
      trial <- free
      trial[others] <- free[others] + step * direction
      result <- value_gradient(trial)
      if (is.finite(result$value) && result$value > current$value) break
      step <- step / 2
      if (step < 1e-6) break
    }
    if (step < 1e-6) break
    free <- trial
    current <- result
    steps <- steps + 1
  }
  list(free = free, value = current$value,
       max_gradient = max(abs(current$gradient[1, others])),
       steps = steps, hessian = hessian)
}

# H, the Hessian of the negative log density in every free parameter but
# `fixed`, at `free`, by central differences of the gradient with step `h`,
# `batch` coordinates (2 states each) per TensorFlow call, symmetrised; and
# its largest asymmetry before that, a check on the differences
hessian_others <- function(free, fixed, h = 1e-4, batch = 25) {
  others <- setdiff(seq_along(free), fixed)
  n <- length(others)
  H <- matrix(NA_real_, n, n)
  for (block in split(seq_len(n), ceiling(seq_len(n) / batch))) {
    states <- matrix(free, 2 * length(block), length(free), byrow = TRUE)
    for (b in seq_along(block)) {
      j <- others[block[b]]
      states[2 * b - 1, j] <- free[j] + h
      states[2 * b, j] <- free[j] - h
    }
    gradient <- value_gradient(states)$gradient[, others, drop = FALSE]
    for (b in seq_along(block)) {
      H[, block[b]] <- -(gradient[2 * b - 1, ] - gradient[2 * b, ]) / (2 * h)
    }
  }
  list(H = (H + t(H)) / 2, asymmetry = max(abs(H - t(H))))
}

for (name in profile_floors) {
  # in the non-centred model, whose free state is laid out differently
  column <- free_column(noncentred$model, name)
  floor_free_check(noncentred$model, starts[[1]], name)
  for (i in indices) {
    out <- point_file(name, i)
    if (file.exists(out)) {
      report("%s point %d (%g) done already", name, i, profile_grid[i])
      next
    }
    f <- profile_grid[i]
    time_start <- Sys.time()
    runs <- lapply(names(starts), function(mode) {
      start <- starts[[mode]]
      start[column] <- floor_free(f)
      run <- optimise_others(start, column)
      report("%s = %g from the %s mode: log density %.3f -> %.3f, %d evaluations, %.0f s, max |gradient| %.2g, %s",
             name, f, mode, run$start_value, run$value, run$evaluations,
             run$seconds, run$max_gradient, run$message)
      run
    })
    names(runs) <- names(starts)
    best <- names(runs)[which.max(vapply(runs, `[[`, numeric(1), "value"))]
    hessian_time <- system.time(
      polished <- newton_others(runs[[best]]$free, column))[["elapsed"]]
    report("%s = %g: %d Newton steps from the %s optimum, log density %.3f -> %.3f, max |gradient| %.2g -> %.2g, %.0f s",
           name, f, polished$steps, best, runs[[best]]$value, polished$value,
           runs[[best]]$max_gradient, polished$max_gradient, hessian_time)
    hessian <- polished$hessian
    eigenvalues <- eigen(hessian$H, symmetric = TRUE, only.values = TRUE)$values
    positive <- all(eigenvalues > 0)
    half_log_det <- if (positive) -0.5 * sum(log(eigenvalues)) else NA_real_
    floors_at <- floor_values(noncentred$model$dag$trace_values(
      matrix(polished$free, nrow = 1))[, floors, drop = FALSE])[1, ]
    stopifnot(abs(floors_at[[name]] - f) < 1e-10)
    point <- list(
      floor = name, index = i, value = f, u = qlogis(f),
      height = polished$value,
      half_log_det = half_log_det,
      laplace = polished$value + half_log_det,
      # on the floor's own scale, less log |d logit f / d f|
      laplace_floor_scale = polished$value + half_log_det - log(f * (1 - f)),
      n_nonpositive = sum(eigenvalues <= 0),
      min_eigenvalue = min(eigenvalues),
      max_eigenvalue = max(eigenvalues),
      hessian_asymmetry = hessian$asymmetry,
      best_start = best,
      start_heights = paste(sprintf("%s %.3f", names(runs),
                                    vapply(runs, `[[`, numeric(1), "value")),
                            collapse = "; "),
      other_floors = paste(sprintf("%s %.4g", setdiff(floors, name),
                                   floors_at[setdiff(floors, name)]),
                           collapse = "; "),
      lbfgs_max_gradient = runs[[best]]$max_gradient,
      lbfgs_converged = runs[[best]]$converged,
      lbfgs_evaluations = runs[[best]]$evaluations,
      newton_steps = polished$steps,
      max_gradient = polished$max_gradient,
      converged = polished$max_gradient < 0.01,
      optimise_seconds = sum(vapply(runs, `[[`, numeric(1), "seconds")),
      newton_seconds = hessian_time,
      seconds = as.numeric(difftime(Sys.time(), time_start, units = "secs")),
      peak_memory_gb = peak_memory_gb(),
      free = polished$free,
      eigenvalues = eigenvalues)
    saveRDS(point, out)
    invisible(gc())
    report("%s = %g: height %.3f, -1/2 log|H| %s, Laplace %s, %d non-positive eigenvalues; %.0f s (Newton and Hessians %.0f s); peak %.1f GB",
           name, f, point$height,
           format(point$half_log_det, digits = 8),
           format(point$laplace, digits = 10), point$n_nonpositive,
           point$seconds, hessian_time, point$peak_memory_gb)
  }
}
report("done; run with `summarise` to combine the points")
