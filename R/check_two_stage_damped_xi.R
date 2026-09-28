# Checks of the damped accumulation option of the two-stage correction (#21;
# doc/two_stage_plan.md, "Damped accumulation"): fit_correction(damped_xi =
# TRUE), xi_t = psi xi_{t-1} + eta_t, and its use in stage B, prediction and
# the forecast beyond T.
#
#   Rscript R/check_two_stage_damped_xi.R
#
# 1. Defaults unchanged: the template and R code at DAMPED_REF (default HEAD)
#    and the current code, fitted to the same simulated data with the defaults
#    (undamped), with p, give the same objective, hyperparameters, latent mode,
#    Hessian, stage-B fit and predictive draws (including forecast years).
# 2. Recovery: beta-binomial counts simulated with omega + xi + u + p and a
#    damped xi (psi = 0.7), and with an undamped xi (psi = 1). Stage A and B
#    with psi estimated. Also: the cut-posterior shift at the final D matches a
#    refit, the forecast recursion matches its closed-form mean, and the
#    coverage of lambda at new pixels and in forecast years.

source("R/two_stage_correction.R")
source("R/two_stage_pql.R")
source("R/two_stage_map_functions.R")

scratch <- Sys.getenv("DAMPED_SCRATCH", tempdir())
ref <- Sys.getenv("DAMPED_REF", "HEAD")
set.seed(2029)


# simulated data -----------------------------------------------------------------

ir_africa <- readRDS("data/clean/all_gambiae_complex_data.RDS")
sites <- ir_africa %>%
  filter(insecticide_type == "Lambda-cyhalothrin",
         !is.na(longitude), !is.na(latitude)) %>%
  distinct(longitude, latitude)
sites_km <- project_km(sites$longitude, sites$latitude)
n_sites <- nrow(sites_km)

t0 <- 2000
T_data <- 2020
T_forecast <- 2025
truth_base <- list(sigma_omega = 0.5, range_omega = 500, sigma_eta = 0.3,
                   range_eta = 600, phi = 0.3, tau = 0.3, sigma_p = 0.3,
                   rho = 0.1)

revisited <- sample.int(n_sites, 200)
pixel_years <- tibble(
  site = c(sample.int(n_sites, 1200, replace = TRUE),
           sample(revisited, 1400, replace = TRUE)),
  year = sample(t0:T_data, 2600, replace = TRUE)
) %>%
  distinct()
n_py <- nrow(pixel_years)
train0 <- bind_rows(pixel_years,
                    pixel_years[sample.int(n_py, round(0.25 * n_py)), ]) %>%
  mutate(x_km = sites_km[site, 1], y_km = sites_km[site, 2], cell = site,
         mosquito_number = sample(20:100, n(), replace = TRUE))

# new points: near training sites (new pixels) in data years, and at training
# sites in forecast years
new0 <- bind_rows(
  tibble(site = sample.int(n_sites, 400, replace = TRUE),
         year = sample((t0 + 1):T_data, 400, replace = TRUE),
         jitter = 35),
  tibble(site = sample(unique(train0$site), 400, replace = TRUE),
         year = sample((T_data + 1):T_forecast, 400, replace = TRUE),
         jitter = 0)
) %>%
  mutate(x_km = sites_km[site, 1] + runif(n(), -jitter, jitter),
         y_km = sites_km[site, 2] + runif(n(), -jitter, jitter),
         cell = ifelse(jitter > 0, -row_number(), site))

mesh <- build_correction_mesh(cbind(train0$x_km, train0$y_km), verbose = FALSE)
mesh_xi <- build_correction_mesh(cbind(train0$x_km, train0$y_km),
                                 max_nodes = 600, verbose = FALSE)
fem <- correction_fem(mesh)
fem_xi <- correction_fem(mesh_xi)

sample_gmrf <- function(Q, n = 1) {
  sample_latent_deviation(Matrix::Cholesky(Matrix::forceSymmetric(Q),
                                           perm = TRUE, LDL = FALSE),
                          nrow(Q), n)
}
kappa <- function(range) sqrt(8) / range
m_fun <- function(x_km, y_km, year) {
  0.5 + 0.1 * (year - t0) - 1.5 * sin(x_km / 700) * cos(y_km / 900)
}
rbb <- function(size, p, rho) {
  a <- p * (1 / rho - 1)
  b <- (1 - p) * (1 / rho - 1)
  rbinom(length(size), size, rbeta(length(size), a, b))
}

# one simulated data set with damping psi (1 = undamped), xi to T_forecast
simulate <- function(psi, seed) {
  set.seed(seed)
  truth <- c(truth_base, psi = psi)
  w_true <- sample_gmrf(matern_precision_r(fem, kappa(truth$range_omega),
                                           truth$sigma_omega))[, 1]
  years_all <- (t0 + 1):T_forecast
  eps <- sample_gmrf(matern_precision_r(fem_xi, kappa(truth$range_eta),
                                        truth$sigma_eta), length(years_all))
  eta <- eps
  for (t in 2:length(years_all)) {
    eta[, t] <- truth$phi * eta[, t - 1] + sqrt(1 - truth$phi ^ 2) * eps[, t]
  }
  x_true <- eta
  for (t in 2:length(years_all)) x_true[, t] <- psi * x_true[, t - 1] + eta[, t]
  all_cells <- unique(c(train0$cell, new0$cell))
  p_true <- setNames(rnorm(length(all_cells), 0, truth$sigma_p), all_cells)
  all_py <- unique(c(paste(train0$cell, train0$year),
                     paste(new0$cell, new0$year)))
  u_true <- setNames(rnorm(length(all_py), 0, truth$tau), all_py)
  latent <- function(d) {
    A <- mesh_basis(mesh, cbind(d$x_km, d$y_km))
    A_xi <- mesh_basis(mesh_xi, cbind(d$x_km, d$y_km))
    xi <- vapply(seq_len(nrow(d)), function(i) {
      if (d$year[i] <= t0) 0 else sum(A_xi[i, ] * x_true[, d$year[i] - t0])
    }, numeric(1))
    m_fun(d$x_km, d$y_km, d$year) + as.vector(A %*% w_true) + xi +
      u_true[paste(d$cell, d$year)] + p_true[as.character(d$cell)]
  }
  train <- train0 %>%
    mutate(m = m_fun(x_km, y_km, year),
           lambda = latent(.),
           died = rbb(mosquito_number, plogis(lambda), truth$rho),
           rho = truth$rho)
  new <- new0 %>% mutate(m = m_fun(x_km, y_km, year), lambda = latent(.))
  list(train = train, new = new, truth = truth)
}

sim_damped <- simulate(0.7, 1)
sim_undamped <- simulate(1, 2)
train <- sim_damped$train
new <- sim_damped$new
cat(sprintf("simulated: %i assays in %i pixel-years, %i pixels\n",
            nrow(train), n_py, n_distinct(train$cell)))


# 1. defaults unchanged ---------------------------------------------------------------

cat("\n1. defaults against the code at", ref, "\n")
old_dir <- file.path(scratch, "damped_head")
dir.create(old_dir, showWarnings = FALSE, recursive = TRUE)
git_show <- function(file, out) {
  status <- system2("git", c("show", paste0(ref, ":", file)), stdout = out)
  stopifnot(status == 0)
}
old_cpp <- file.path(old_dir, "two_stage_correction_head_psi.cpp")
git_show("tmb/two_stage_correction.cpp", old_cpp)
git_show("R/two_stage_correction.R", file.path(old_dir, "correction.R"))
git_show("R/two_stage_pql.R", file.path(old_dir, "pql.R"))
git_show("R/two_stage_map_functions.R", file.path(old_dir, "map_functions.R"))
old <- new.env()
sys.source(file.path(old_dir, "correction.R"), envir = old)
assign("correction_template", old_cpp, envir = old)
sys.source(file.path(old_dir, "pql.R"), envir = old)
sys.source(file.path(old_dir, "map_functions.R"), envir = old)
stopifnot(!grepl("psi", paste(readLines(old_cpp), collapse = "\n")))

fits_a <- list(
  old = old$fit_correction(train, "omega_xi_u", t0 = t0, T = T_data,
                           mesh = mesh, mesh_xi = mesh_xi, pixel_effect = TRUE),
  new = fit_correction(train, "omega_xi_u", t0 = t0, T = T_data,
                       mesh = mesh, mesh_xi = mesh_xi, pixel_effect = TRUE)
)
fits_b <- list(old = old$fit_correction_pql(fits_a$old, train),
               new = fit_correction_pql(fits_a$new, train))
draws <- function(fits) {
  set.seed(1)
  d_old <- old$predict_correction(fits$old, new, n_draws = 200)
  set.seed(1)
  d_new <- predict_correction(fits$new, new, n_draws = 200)
  max(abs(d_old - d_new))
}
node_fields <- function(fits) {
  years <- c(2010, T_data, T_data + 3)
  a <- old$correction_node_fields(fits$old, fits$old$mode, years)
  b <- correction_node_fields(fits$new, fits$new$mode, years)
  max(abs(unlist(a$xi) - unlist(b$xi)))
}
differences <- c(
  objective_A = abs(fits_a$old$opt$objective - fits_a$new$opt$objective),
  hyper_A = max(abs(unlist(fits_a$old$hyper) - unlist(fits_a$new$hyper))),
  mode_A = max(abs(fits_a$old$mode - fits_a$new$mode)),
  hessian_A = max(abs(fits_a$old$H - fits_a$new$H)),
  hyper_B = max(abs(unlist(fits_b$old$hyper) - unlist(fits_b$new$hyper))),
  mode_B = max(abs(fits_b$old$mode - fits_b$new$mode)),
  D_B = max(abs(fits_b$old$precision_obs - fits_b$new$precision_obs)),
  draws_A = draws(fits_a),
  draws_B = draws(fits_b),
  node_fields_B = node_fields(fits_b)
)
print(signif(differences, 3))
stopifnot(identical(names(fits_a$old$hyper), names(fits_a$new$hyper)),
          all(differences == 0))
cat("defaults reproduce the old code exactly\n")
fit_undamped_b <- fits_b$new


# 2. recovery ---------------------------------------------------------------------

fit_damped <- function(d) {
  a <- fit_correction(d, "omega_xi_u", t0 = t0, T = T_data, mesh = mesh,
                      mesh_xi = mesh_xi, pixel_effect = TRUE,
                      damped_xi = TRUE, hyper_hessian = TRUE)
  b <- fit_correction_pql(a, d)
  list(a = a, b = b)
}
fits <- list(`truth psi = 0.7` = fit_damped(sim_damped$train),
             `truth psi = 1` = fit_damped(sim_undamped$train))

row_for <- function(name, fit, stage) {
  psi <- if (is.null(fit$hyper$psi)) 1 else fit$hyper$psi
  se <- fit$psi_se_logit
  tibble(setting = name, stage = stage, psi = psi,
         psi_lower = if (is.na(se)) NA else plogis(qlogis(psi) - 1.96 * se),
         psi_upper = if (is.na(se)) NA else plogis(qlogis(psi) + 1.96 * se),
         reversion_years = 1 / (1 - psi),
         phi = fit$hyper$phi, sigma_eta = fit$hyper$sigma_eta,
         range_eta = fit$hyper$range_eta, tau = fit$hyper$tau,
         sigma_p = fit$hyper$sigma_p, convergence = fit$opt$convergence)
}
recovery <- bind_rows(
  tibble(setting = "truth", stage = "", psi = 0.7, phi = truth_base$phi,
         sigma_eta = truth_base$sigma_eta, range_eta = truth_base$range_eta,
         tau = truth_base$tau, sigma_p = truth_base$sigma_p),
  row_for("undamped model, truth psi = 0.7", fit_undamped_b, "B"),
  bind_rows(lapply(names(fits), function(n) {
    bind_rows(row_for(n, fits[[n]]$a, "A"), row_for(n, fits[[n]]$b, "B"))
  }))
)
cat("\n2. recovery\n")
print(as.data.frame(recovery %>% mutate(across(where(is.numeric),
                                               ~ round(.x, 3)))),
      row.names = FALSE)

# the cut-posterior shift at the final D against a refit of the stage-B
# working model with the hyperparameters (psi included) fixed
f <- fits[["truth psi = 0.7"]]$b
m_perturbed <- train$m + 0.3 * sin(train$x_km / 400)
shifted <- f$mode + correction_mode_shift(f, m_perturbed)[, 1]
data <- f$tmb_data
data$m <- m_perturbed
obj <- correction_adfun(data, f$par_list, f$variant, fix_hyper = TRUE)
invisible(obj$fn(obj$par))
shift_error <- max(abs(shifted - obj$env$last.par[obj$env$random]))
data$m <- f$m_ref
obj <- correction_adfun(data, f$par_list, f$variant, fix_hyper = TRUE)
invisible(obj$fn(obj$par))
mode_error <- max(abs(f$mode - obj$env$last.par[obj$env$random]))
cat(sprintf("stage B, psi estimated: shift vs refit %.2e, mode vs refit %.2e\n",
            shift_error, mode_error))
stopifnot(shift_error < 1e-8, mode_error < 1e-8)

# the forecast: the mean of many simulated forecast draws from the mode against
# the closed-form mean recursion (innovations = NULL), at T + 1 and T + 5
theta <- matrix(f$mode, length(f$mode), 4000)
Q_eta <- matern_precision_r(f$fem_xi, f$hyper$kappa_eta, f$hyper$sigma_eta)
Q_chol <- Matrix::Cholesky(Matrix::forceSymmetric(Q_eta), perm = TRUE,
                           LDL = FALSE, super = TRUE)
innovations <- lapply(1:5, function(h)
  sample_latent_deviation(Q_chol, f$mesh_xi$n, 4000))
sim_nodes <- correction_node_fields(f, theta, T_data + c(1, 5), innovations)
mean_nodes <- correction_node_fields(f, f$mode, T_data + c(1, 5))
forecast_error <- sapply(1:2, function(j)
  max(abs(rowMeans(sim_nodes$xi[[j]]) - mean_nodes$xi[[j]][, 1])))
# the stationary SD of xi under the fitted (psi, phi, sigma_eta), against the
# simulated forecast SD far ahead (h = 5 is not yet stationary, so this is a
# loose check that the damping bounds the spread)
psi <- f$hyper$psi
phi <- f$hyper$phi
stationary_sd <- f$hyper$sigma_eta *
  sqrt((1 + psi * phi) / ((1 - psi ^ 2) * (1 - psi * phi)))
cat(sprintf(paste("forecast mean: max |simulated - closed form| %.3f (h = 1),",
                  "%.3f (h = 5); node SD at h = 5 median %.2f, stationary",
                  "SD %.2f\n"),
            forecast_error[1], forecast_error[2],
            median(row_sds(sim_nodes$xi[[2]])), stationary_sd))
stopifnot(all(forecast_error < 0.05))

# coverage of lambda: new pixels in data years, training pixels in forecast
# years
coverage <- function(fit, sim) {
  set.seed(3)
  d <- predict_correction(fit, sim$new, n_draws = 1000)
  lower <- apply(d, 2, quantile, c(0.025, 0.25))
  upper <- apply(d, 2, quantile, c(0.975, 0.75))
  set_of <- ifelse(sim$new$year > T_data, "forecast", "interpolation")
  bind_rows(lapply(split(seq_len(nrow(sim$new)), set_of), function(i) {
    tibble(cover50 = mean(sim$new$lambda[i] >= lower[2, i] &
                            sim$new$lambda[i] <= upper[2, i]),
           cover95 = mean(sim$new$lambda[i] >= lower[1, i] &
                            sim$new$lambda[i] <= upper[1, i]),
           rmse = sqrt(mean((colMeans(d[, i]) - sim$new$lambda[i]) ^ 2)),
           width95 = mean(upper[1, i] - lower[1, i]))
  }), .id = "set")
}
cat("\ncoverage of lambda, stage B:\n")
print(as.data.frame(bind_rows(
  coverage(fit_undamped_b, sim_damped) %>% mutate(fit = "undamped, truth 0.7"),
  coverage(fits[["truth psi = 0.7"]]$b, sim_damped) %>%
    mutate(fit = "damped, truth 0.7"),
  coverage(fits[["truth psi = 1"]]$b, sim_undamped) %>%
    mutate(fit = "damped, truth 1")
) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))), row.names = FALSE)
