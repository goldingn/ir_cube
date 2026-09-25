# Simulation checks of the two-stage correction model (#21,
# R/two_stage_correction.R).
#
# Data are simulated from the correction model itself, on a realistic domain:
# the bioassay sites of one insecticide, projected to km, with ~2000 assays over
# 20 years and some pixel-years assayed more than once. Both variants are then
# fitted and checked for:
#
#   1. recovery of the hyperparameters;
#   2. coverage of prediction intervals for the true latent lambda (excluding
#      assay noise), at held-out pixel-years in the training period
#      (interpolation), at training pixel-years, and 1-5 years beyond the last
#      data year (forecast);
#   3. the cut-posterior shift of the latent mode, H^-1 A' D (m_ref - m_new),
#      against refitting the latent mode in TMB with the perturbed m and the
#      hyperparameters held fixed: these should agree to numerical precision,
#      because the latent posterior is exactly Gaussian given the
#      hyperparameters;
#   4. that a constant offset in m is absorbed by omega (the initial condition
#      correction) near the data, rather than by xi or u, and passed through to
#      the predictions far from data.
#
# Runtimes of the fits and of prediction are reported.
#
# The truth is simulated on the fitting mesh, so the fields are exactly
# representable; this checks the implementation, not mesh discretisation error.
# Observations are Gaussian on the logit scale with the stage-A variances, i.e.
# the model is correctly specified: the empirical-logit bias at extreme
# mortality is a separate diagnostic (issue #21, fitting cascade step 2).

source("R/two_stage_correction.R")

set.seed(2026)

# domain and design --------------------------------------------------------

ir_africa <- readRDS("data/clean/all_gambiae_complex_data.RDS")
sites <- ir_africa %>%
  filter(insecticide_type == "Lambda-cyhalothrin",
         !is.na(longitude), !is.na(latitude)) %>%
  distinct(longitude, latitude)
sites_km <- project_km(sites$longitude, sites$latitude)
n_sites <- nrow(sites_km)

t0 <- 2000
T_data <- 2020
max_horizon <- 5
years_all <- (t0 + 1):(T_data + max_horizon)

n_obs <- 2000
n_pixel_years <- 1500

# truth
truth <- list(
  sigma_omega = 0.5,
  range_omega = 500,
  sigma_eta = 0.1,
  range_eta = 300,
  phi = 0.7,
  tau = 0.3
)
rho <- 0.15

# training pixel-years: sites and years at random over [t0, T], then assays
# drawn from them with replacement, so ~25% of pixel-years have replicates
pixel_years <- tibble(
  site = sample.int(n_sites, n_pixel_years, replace = TRUE),
  year = sample(t0:T_data, n_pixel_years, replace = TRUE)
) %>%
  distinct()
train <- pixel_years[sample.int(nrow(pixel_years), n_obs, replace = TRUE), ]

# held-out interpolation points: sites jittered by up to ~50 km, in new pixels
# (so fresh u), in years within the training period
n_interp <- 400
interp <- tibble(
  site = sample.int(n_sites, n_interp, replace = TRUE),
  year = sample((t0 + 1):T_data, n_interp, replace = TRUE),
  jitter_x = runif(n_interp, -35, 35),
  jitter_y = runif(n_interp, -35, 35),
  type = "interpolation"
)
# training pixel-years, where u has been observed
n_insample <- 300
insample <- train %>%
  distinct(site, year) %>%
  slice_sample(n = n_insample) %>%
  mutate(jitter_x = 0, jitter_y = 0, type = "training pixel-year")
# forecasts 1-5 years ahead, at training sites
n_forecast <- 500
forecast <- tibble(
  site = sample(unique(train$site), n_forecast, replace = TRUE),
  year = sample((T_data + 1):(T_data + max_horizon), n_forecast,
                replace = TRUE),
  jitter_x = 0, jitter_y = 0,
  type = "forecast"
)

add_coords <- function(df) {
  df %>%
    mutate(
      x_km = sites_km[site, 1] + if ("jitter_x" %in% names(df)) jitter_x else 0,
      y_km = sites_km[site, 2] + if ("jitter_y" %in% names(df)) jitter_y else 0
    )
}
train <- add_coords(train) %>% mutate(cell = site)
new <- bind_rows(interp, insample, forecast) %>%
  add_coords() %>%
  # jittered interpolation points are new pixels
  mutate(cell = ifelse(type == "interpolation", -row_number(), site))

# the meshes fit_correction() would build by default: omega on the finer mesh,
# xi on a coarser one
mesh <- build_correction_mesh(cbind(train$x_km, train$y_km))
mesh_xi <- build_correction_mesh(cbind(train$x_km, train$y_km),
                                 max_nodes = 600)
fem <- correction_fem(mesh)
fem_xi <- correction_fem(mesh_xi)
cat(sprintf(
  "domain: %i sites, %i assays, %i pixel-years, meshes %i (omega) and %i (xi) nodes\n",
  n_distinct(train$site), nrow(train), nrow(distinct(train, cell, year)),
  mesh$n, mesh_xi$n
))


# simulate the truth --------------------------------------------------------------

sample_gmrf <- function(Q, n = 1) {
  chol_Q <- Matrix::Cholesky(Matrix::forceSymmetric(Q), perm = TRUE,
                             LDL = FALSE)
  sample_latent_deviation(chol_Q, nrow(Q), n)
}

kappa <- function(range) sqrt(8) / range
w_true <- sample_gmrf(matern_precision_r(fem, kappa(truth$range_omega),
                                         truth$sigma_omega))[, 1]

# eta: stationary AR(1) in time of Matern fields, x = cumulative sum from t0
Q_eta <- matern_precision_r(fem_xi, kappa(truth$range_eta), truth$sigma_eta)
eps <- sample_gmrf(Q_eta, length(years_all))
eta_true <- eps
for (t in 2:length(years_all)) {
  eta_true[, t] <- truth$phi * eta_true[, t - 1] +
    sqrt(1 - truth$phi ^ 2) * eps[, t]
}
x_true <- t(apply(eta_true, 1, cumsum))

# a smooth dynamical model prediction, with a trend and spatial structure
m_fun <- function(x_km, y_km, year) {
  -1 + 0.12 * (year - t0) + 0.8 * sin(x_km / 700) * cos(y_km / 900)
}

latent_truth <- function(df, u) {
  A <- mesh_basis(mesh, cbind(df$x_km, df$y_km))
  A_xi <- mesh_basis(mesh_xi, cbind(df$x_km, df$y_km))
  omega <- as.vector(A %*% w_true)
  xi <- numeric(nrow(df))
  after <- df$year > t0
  for (i in which(after)) {
    xi[i] <- sum(A_xi[i, ] * x_true[, df$year[i] - t0])
  }
  m_fun(df$x_km, df$y_km, df$year) + omega + xi + u
}

# pixel-year effects, shared by every assay (and prediction) in a pixel-year
all_keys <- unique(c(paste(train$cell, train$year), paste(new$cell, new$year)))
u_true <- setNames(rnorm(length(all_keys), 0, truth$tau), all_keys)

train <- train %>%
  mutate(
    m = m_fun(x_km, y_km, year),
    lambda = latent_truth(., u_true[paste(cell, year)]),
    mosquito_number = sample(20:100, n(), replace = TRUE),
    # stage-A variances come from beta-binomial counts, then the response is
    # drawn from the Gaussian model with those variances
    died = {
      p <- plogis(lambda)
      a <- p * (1 - rho) / rho
      b <- (1 - p) * (1 - rho) / rho
      rbinom(n(), mosquito_number, rbeta(n(), a, b))
    },
    rho = rho,
    v = empirical_logit(died, mosquito_number, rho)$v,
    z = lambda + rnorm(n(), 0, sqrt(v))
  )
new <- new %>%
  mutate(m = m_fun(x_km, y_km, year),
         lambda = latent_truth(., u_true[paste(cell, year)]))


# fit both variants ----------------------------------------------------------------

fits <- list()
for (variant in c("omega_u", "omega_xi_u")) {
  fits[[variant]] <- fit_correction(train, variant = variant, t0 = t0,
                                    T = T_data, mesh = mesh,
                                    mesh_xi = mesh_xi)
  cat(sprintf("\n%s: nlminb convergence %i (%s), %i latent, fit %.1f s ",
              variant, fits[[variant]]$opt$convergence,
              fits[[variant]]$opt$message,
              length(fits[[variant]]$mode),
              fits[[variant]]$timings["total"]))
  cat(sprintf("(optimise %.1f s, Hessian + Cholesky %.1f s)\n",
              fits[[variant]]$timings["optimise"],
              fits[[variant]]$timings["hessian"]))
}

cat("\n1. hyperparameters (truth vs estimates)\n")
hyper_table <- tibble(parameter = names(truth), truth = unlist(truth)) %>%
  mutate(
    omega_u = sapply(parameter, function(p) {
      if (is.null(fits$omega_u$hyper[[p]])) NA else fits$omega_u$hyper[[p]]
    }),
    omega_xi_u = sapply(parameter, function(p) fits$omega_xi_u$hyper[[p]])
  )
print(hyper_table, digits = 3)


# prediction coverage -----------------------------------------------------------------

coverage_table <- function(draws, df) {
  levels <- c(0.5, 0.8, 0.95)
  res <- lapply(levels, function(level) {
    lower <- apply(draws, 2, quantile, (1 - level) / 2)
    upper <- apply(draws, 2, quantile, 1 - (1 - level) / 2)
    tapply(df$lambda >= lower & df$lambda <= upper, df$type, mean)
  })
  out <- do.call(cbind, res)
  colnames(out) <- paste0(levels * 100, "%")
  out <- cbind(out,
               rmse = tapply((colMeans(draws) - df$lambda) ^ 2, df$type,
                             function(x) sqrt(mean(x))),
               mean_sd = tapply(apply(draws, 2, sd), df$type, mean))
  round(out, 3)
}

cat("\n2. interval coverage of the true latent lambda (1000 draws)\n")
n_draws <- 1000
predictions <- list()
for (variant in names(fits)) {
  time <- system.time(
    predictions[[variant]] <- predict_correction(fits[[variant]], new,
                                                 n_draws = n_draws)
  )["elapsed"]
  cat(sprintf("\n%s (prediction at %i points: %.1f s)\n", variant, nrow(new),
              time))
  print(coverage_table(predictions[[variant]], new))
}
cat("\nforecast coverage (95%) by horizon, omega_xi_u:\n")
fc <- new$type == "forecast"
draws_fc <- predictions$omega_xi_u[, fc]
lower <- apply(draws_fc, 2, quantile, 0.025)
upper <- apply(draws_fc, 2, quantile, 0.975)
print(round(tapply(new$lambda[fc] >= lower & new$lambda[fc] <= upper,
                   new$year[fc] - T_data, mean), 3))
cat("mechanistic-only RMSE (m vs lambda) by type:\n")
print(round(tapply((new$m - new$lambda) ^ 2, new$type,
                   function(x) sqrt(mean(x))), 3))


# cut-posterior shift vs refit ------------------------------------------------------------

cat("\n3. cut-posterior shift vs refitting the latent mode with perturbed m\n")
m_perturbed <- train$m + 0.3 * sin(train$x_km / 400) +
  rnorm(nrow(train), 0, 0.1)
for (variant in names(fits)) {
  fit <- fits[[variant]]
  shifted <- fit$mode + correction_mode_shift(fit, m_perturbed)[, 1]
  data <- fit$tmb_data
  data$m <- m_perturbed
  obj <- correction_adfun(data, fit$par_list, variant, fix_hyper = TRUE)
  obj$fn(obj$par)
  refit <- obj$env$last.par[obj$env$random]
  cat(sprintf(
    "%s: max |shift formula - refit| = %.2e (max |shift| = %.3f)\n",
    variant, max(abs(shifted - refit)), max(abs(shifted - fit$mode))
  ))
}


# constant offset ---------------------------------------------------------------------------

cat("\n4. constant offset delta = 0.5 added to m at the training data\n")
delta <- 0.5
fit <- fits$omega_xi_u
shift <- correction_mode_shift(fit, train$m + delta)[, 1]
A_train <- fit$A_latent
contribution <- sapply(fit$blocks, function(idx) {
  as.vector(A_train[, idx, drop = FALSE] %*% shift[idx])
})
cat("mean change in each latent term at the training assays (as a fraction ",
    "of -delta):\n", sep = "")
print(round(colMeans(contribution) / -delta, 3))
# change in predicted lambda at held-out points = delta (in m) + latent change
A_new <- fm_basis(mesh, loc = cbind(new$x_km, new$y_km))
omega_change <- as.vector(A_new %*% shift[fit$blocks$w_omega])
cat("mean fraction of delta passed through to predictions of omega + m:\n")
print(round(tapply(1 + omega_change / delta, new$type, mean), 3))
# far from data: the mesh node furthest from any training site
node_distance <- apply(mesh$loc[, 1:2], 1, function(p) {
  min(sqrt((train$x_km - p[1]) ^ 2 + (train$y_km - p[2]) ^ 2))
})
far <- which.max(node_distance)
cat(sprintf("at the mesh node furthest from data (%.0f km): %.3f\n",
            node_distance[far], 1 + shift[fit$blocks$w_omega][far] / delta))


# mechanistic uncertainty -----------------------------------------------------------------------

cat("\n5. propagation of dynamical-model draws (cut posterior)\n")
n_m_draws <- 100
# dynamical draws: a per-draw logit offset plus a smooth per-draw trend error
m_error <- function(df, a, b) {
  outer(a, rep(1, nrow(df))) + outer(b, (df$year - t0) / 20)
}
a <- rnorm(n_m_draws, 0, 0.3)
b <- rnorm(n_m_draws, 0, 0.3)
m_draws_train <- sweep(m_error(train, a, b), 2, train$m, "+")
m_draws_new <- sweep(m_error(new, a, b), 2, new$m, "+")
time <- system.time(
  draws_m <- predict_correction(fit, new, m_draws_train = m_draws_train,
                                m_draws_new = m_draws_new, n_draws = n_draws)
)["elapsed"]
cat(sprintf("prediction with %i dynamical draws: %.1f s\n", n_m_draws, time))
cat("mean predictive SD by type, without vs with dynamical draws:\n")
print(round(rbind(
  without = tapply(apply(predictions$omega_xi_u, 2, sd), new$type, mean),
  with = tapply(apply(draws_m, 2, sd), new$type, mean),
  sd_of_m_draws = tapply(apply(m_draws_new, 2, sd), new$type, mean)
), 3))
