# Checks of the joint rho option of the two-stage correction (#21;
# doc/two_stage_plan.md, "Joint rho"): fit_correction(estimate_rho = TRUE) and
# its use in stage B (R/two_stage_pql.R).
#
#   Rscript R/check_two_stage_joint_rho.R
#
# 1. Defaults unchanged: the template and R code at JOINT_RHO_REF (default
#    HEAD) and the current code, fitted to the same simulated data with the
#    defaults (fixed rho), with and without p, give the same objective,
#    hyperparameters, latent mode, Hessian, stage-B fit and predictive draws.
# 2. Recovery: beta-binomial counts simulated with omega + xi + u + p and a
#    known rho, with replicate assays in some pixel-years at about the real
#    data's rate. Stage A and stage B (+p) with rho estimated, the prior centred
#    on the true rho and on twice it (as if the external estimate were too
#    large), and with a vague prior; rho fixed at the truth for comparison.
#    Also: the cut-posterior shift at the final D matches a refit.

source("R/two_stage_correction.R")
source("R/two_stage_pql.R")

scratch <- Sys.getenv("JOINT_RHO_SCRATCH", tempdir())
ref <- Sys.getenv("JOINT_RHO_REF", "HEAD")
set.seed(2028)


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
truth <- list(sigma_omega = 0.5, range_omega = 500, sigma_eta = 0.1,
              range_eta = 300, phi = 0.7, tau = 0.3, sigma_p = 0.3, rho = 0.1)

# about 1.25 assays per pixel-year, as in the folds (24,318 assays in 19,453
# pixel-years on interpolation); a revisited subset of sites identifies p
revisited <- sample.int(n_sites, 200)
pixel_years <- tibble(
  site = c(sample.int(n_sites, 1200, replace = TRUE),
           sample(revisited, 1400, replace = TRUE)),
  year = sample(t0:T_data, 2600, replace = TRUE)
) %>%
  distinct()
n_py <- nrow(pixel_years)
train <- bind_rows(pixel_years,
                   pixel_years[sample.int(n_py, round(0.25 * n_py)), ]) %>%
  mutate(x_km = sites_km[site, 1], y_km = sites_km[site, 2], cell = site)

new <- tibble(site = sample.int(n_sites, 400, replace = TRUE),
              year = sample((t0 + 1):T_data, 400, replace = TRUE)) %>%
  mutate(x_km = sites_km[site, 1] + runif(n(), -35, 35),
         y_km = sites_km[site, 2] + runif(n(), -35, 35),
         cell = -row_number())

mesh <- build_correction_mesh(cbind(train$x_km, train$y_km), verbose = FALSE)
mesh_xi <- build_correction_mesh(cbind(train$x_km, train$y_km),
                                 max_nodes = 600, verbose = FALSE)
fem <- correction_fem(mesh)
fem_xi <- correction_fem(mesh_xi)

sample_gmrf <- function(Q, n = 1) {
  sample_latent_deviation(Matrix::Cholesky(Matrix::forceSymmetric(Q),
                                           perm = TRUE, LDL = FALSE),
                          nrow(Q), n)
}
kappa <- function(range) sqrt(8) / range
w_true <- sample_gmrf(matern_precision_r(fem, kappa(truth$range_omega),
                                         truth$sigma_omega))[, 1]
years_all <- (t0 + 1):T_data
eps <- sample_gmrf(matern_precision_r(fem_xi, kappa(truth$range_eta),
                                      truth$sigma_eta), length(years_all))
eta_true <- eps
for (t in 2:length(years_all)) {
  eta_true[, t] <- truth$phi * eta_true[, t - 1] +
    sqrt(1 - truth$phi ^ 2) * eps[, t]
}
x_true <- t(apply(eta_true, 1, cumsum))
# a spread of mortality with some saturation, as in the data
m_fun <- function(x_km, y_km, year) {
  0.5 + 0.1 * (year - t0) - 1.5 * sin(x_km / 700) * cos(y_km / 900)
}
all_cells <- unique(c(train$cell, new$cell))
p_true <- setNames(rnorm(length(all_cells), 0, truth$sigma_p), all_cells)
all_py <- unique(c(paste(train$cell, train$year), paste(new$cell, new$year)))
u_true <- setNames(rnorm(length(all_py), 0, truth$tau), all_py)

latent_truth <- function(d) {
  A <- mesh_basis(mesh, cbind(d$x_km, d$y_km))
  A_xi <- mesh_basis(mesh_xi, cbind(d$x_km, d$y_km))
  xi <- vapply(seq_len(nrow(d)), function(i) {
    if (d$year[i] <= t0) 0 else sum(A_xi[i, ] * x_true[, d$year[i] - t0])
  }, numeric(1))
  m_fun(d$x_km, d$y_km, d$year) + as.vector(A %*% w_true) + xi +
    u_true[paste(d$cell, d$year)] + p_true[as.character(d$cell)]
}
rbb <- function(size, p, rho) {
  a <- p * (1 / rho - 1)
  b <- (1 - p) * (1 / rho - 1)
  rbinom(length(size), size, rbeta(length(size), a, b))
}
train <- train %>%
  mutate(m = m_fun(x_km, y_km, year),
         lambda = latent_truth(.),
         mosquito_number = sample(20:100, n(), replace = TRUE),
         died = rbb(mosquito_number, plogis(lambda), truth$rho),
         rho = truth$rho)
new <- new %>% mutate(m = m_fun(x_km, y_km, year), lambda = latent_truth(.))
cat(sprintf(paste("simulated: %i assays in %i pixel-years (%i with > 1 assay),",
                  "%i pixels; %.0f%% at 0%%, %.0f%% at 100%%\n"),
            nrow(train), n_py,
            sum(table(paste(train$cell, train$year)) > 1),
            n_distinct(train$cell), 100 * mean(train$died == 0),
            100 * mean(train$died == train$mosquito_number)))


# 1. defaults unchanged ---------------------------------------------------------------

cat("\n1. defaults against the code at", ref, "\n")
old_dir <- file.path(scratch, "joint_rho_head")
dir.create(old_dir, showWarnings = FALSE, recursive = TRUE)
git_show <- function(file, out) {
  status <- system2("git", c("show", paste0(ref, ":", file)), stdout = out)
  stopifnot(status == 0)
}
old_cpp <- file.path(old_dir, "two_stage_correction_head_rho.cpp")
git_show("tmb/two_stage_correction.cpp", old_cpp)
git_show("R/two_stage_correction.R", file.path(old_dir, "correction.R"))
git_show("R/two_stage_pql.R", file.path(old_dir, "pql.R"))
old <- new.env()
sys.source(file.path(old_dir, "correction.R"), envir = old)
assign("correction_template", old_cpp, envir = old)
sys.source(file.path(old_dir, "pql.R"), envir = old)
stopifnot(!grepl("estimate_rho", paste(readLines(old_cpp), collapse = "\n")))

for (pixel_effect in c(FALSE, TRUE)) {
  fits_a <- list(
    old = old$fit_correction(train, "omega_xi_u", t0 = t0, T = T_data,
                             mesh = mesh, mesh_xi = mesh_xi,
                             pixel_effect = pixel_effect),
    new = fit_correction(train, "omega_xi_u", t0 = t0, T = T_data,
                         mesh = mesh, mesh_xi = mesh_xi,
                         pixel_effect = pixel_effect)
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
  differences <- c(
    objective_A = abs(fits_a$old$opt$objective - fits_a$new$opt$objective),
    hyper_A = max(abs(unlist(fits_a$old$hyper) - unlist(fits_a$new$hyper))),
    mode_A = max(abs(fits_a$old$mode - fits_a$new$mode)),
    hessian_A = max(abs(fits_a$old$H - fits_a$new$H)),
    hyper_B = max(abs(unlist(fits_b$old$hyper) - unlist(fits_b$new$hyper))),
    mode_B = max(abs(fits_b$old$mode - fits_b$new$mode)),
    D_B = max(abs(fits_b$old$precision_obs - fits_b$new$precision_obs)),
    draws_A = draws(fits_a),
    draws_B = draws(fits_b)
  )
  cat("pixel_effect =", pixel_effect, "\n")
  print(signif(differences, 3))
  stopifnot(identical(names(fits_a$old$hyper), names(fits_a$new$hyper)),
            all(differences == 0))
}
cat("defaults reproduce the old code exactly\n")
fit_fixed_b <- fits_b$new


# 2. recovery --------------------------------------------------------------------

cat("\n2. rho estimated, truth rho =", truth$rho, "\n")
fit_rho <- function(centre, sd = 0.56) {
  d <- mutate(train, rho = centre)
  a <- fit_correction(d, "omega_xi_u", t0 = t0, T = T_data, mesh = mesh,
                      mesh_xi = mesh_xi, pixel_effect = TRUE,
                      estimate_rho = TRUE,
                      priors = correction_priors(rho_logit_sd = sd))
  b <- fit_correction_pql(a, d)
  list(a = a, b = b)
}
settings <- list(
  "centre = truth" = list(centre = truth$rho, sd = 0.56),
  "centre = 2 x truth" = list(centre = 2 * truth$rho, sd = 0.56),
  "centre = 2 x truth, sd 3" = list(centre = 2 * truth$rho, sd = 3)
)
fits <- lapply(settings, function(s) fit_rho(s$centre, s$sd))

row_for <- function(name, fit, stage) {
  se <- fit$rho_se_logit
  tibble(setting = name, stage = stage,
         rho = if (is.null(fit$hyper$rho)) truth$rho else fit$hyper$rho,
         rho_lower = if (is.na(se)) NA else plogis(qlogis(rho) - 1.96 * se),
         rho_upper = if (is.na(se)) NA else plogis(qlogis(rho) + 1.96 * se),
         tau = fit$hyper$tau, sigma_p = fit$hyper$sigma_p,
         sigma_omega = fit$hyper$sigma_omega,
         sigma_eta = fit$hyper$sigma_eta,
         convergence = fit$opt$convergence)
}
recovery <- bind_rows(
  tibble(setting = "truth", stage = "", rho = truth$rho, tau = truth$tau,
         sigma_p = truth$sigma_p, sigma_omega = truth$sigma_omega,
         sigma_eta = truth$sigma_eta),
  row_for("rho fixed at truth", fit_fixed_b, "B"),
  bind_rows(lapply(names(fits), function(n) {
    bind_rows(row_for(n, fits[[n]]$a, "A"), row_for(n, fits[[n]]$b, "B"))
  }))
)
print(as.data.frame(recovery %>% mutate(across(where(is.numeric),
                                               ~ round(.x, 3)))),
      row.names = FALSE)

# the cut-posterior shift at the final D against a refit of the stage-B
# working model with the hyperparameters (rho included) fixed
f <- fits[["centre = 2 x truth"]]$b
m_perturbed <- train$m + 0.3 * sin(train$x_km / 400)
shifted <- f$mode + correction_mode_shift(f, m_perturbed)[, 1]
data <- f$tmb_data
data$m <- m_perturbed
obj <- correction_adfun(data, f$par_list, f$variant, fix_hyper = TRUE)
invisible(obj$fn(obj$par))
shift_error <- max(abs(shifted - obj$env$last.par[obj$env$random]))
# and the mode itself: the refit at m_ref reproduces the stage-B mode
data$m <- f$m_ref
obj <- correction_adfun(data, f$par_list, f$variant, fix_hyper = TRUE)
invisible(obj$fn(obj$par))
mode_error <- max(abs(f$mode - obj$env$last.par[obj$env$random]))
cat(sprintf("stage B, rho estimated: shift vs refit %.2e, mode vs refit %.2e\n",
            shift_error, mode_error))
stopifnot(shift_error < 1e-8, mode_error < 1e-8)

# coverage of lambda at new pixels
coverage <- function(fit) {
  set.seed(3)
  d <- predict_correction(fit, new, n_draws = 1000)
  lower <- apply(d, 2, quantile, c(0.025, 0.25))
  upper <- apply(d, 2, quantile, c(0.975, 0.75))
  c(cover50 = mean(new$lambda >= lower[2, ] & new$lambda <= upper[2, ]),
    cover95 = mean(new$lambda >= lower[1, ] & new$lambda <= upper[1, ]),
    rmse = sqrt(mean((colMeans(d) - new$lambda) ^ 2)))
}
cat("\ncoverage of lambda at new pixels, stage B:\n")
print(round(rbind(`rho fixed at truth` = coverage(fit_fixed_b),
                  t(sapply(fits, function(x) coverage(x$b)))), 3))
