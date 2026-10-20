# Simulation check of the final two-stage model (#21, #34, #55;
# R/two_stage_correction.R).
#
#   Rscript R/check_two_stage.R
#
# Data are simulated from the final model, omega + xi + u + p + sy, on the
# Lambda-cyhalothrin bioassay sites with the final model's meshes, as
# beta-binomial counts. The dynamical prediction m puts much of the latent
# lambda near mortality 1 (and some near 0), as in the real data, where the
# empirical logit of stage A is pinned by its +0.5 correction. A subset of
# sites is revisited often, so that p is identified, and each pixel-year
# belongs to a study, most studies covering several pixels in a year, so that
# sy is identified. The truth is simulated on the fitting meshes, so this
# checks the implementation, not discretisation.
#
#   1. the stage-A mode is H^-1 A' D (z - m): the Hessian identity the
#      cut-posterior shift relies on;
#   2. the cut-posterior shift H^-1 A' D (m_ref - m_new) of stage A against
#      refitting its latent mode in TMB with the perturbed m and the
#      hyperparameters fixed (exact for the Gaussian working model), and the
#      error of the linearised shift H^-1 A' W (m_ref - m_new) of the final
#      (Laplace) fit, which is read, not asserted;
#   3. the final fit's latent mode solves the score equation of the
#      beta-binomial log posterior;
#   4. the CCD design weights reproduce a Gaussian's mean and covariance;
#   5. recovery of the hyperparameters: at theta_hat and their CCD posterior
#      means and 95% intervals;
#   6. coverage, at new pixels, new years at training pixels, training
#      pixel-years and 1-5 year forecasts: of a new assay's lambda by
#      predict_correction() (fresh noise), stage A vs the final fit mixed over
#      the CCD points, and of the target m + omega + xi by the final fit's
#      draws without noise;
#   7. the posterior-mean fields (correction_node_draws(mean = TRUE)) against
#      the mean of the draws, for the smooth terms the maps show.
#
# Checks 1-4 are identities, met to ~1e-10 when the code is right; the script
# stops if any misses by more than `tolerance`. The others are statistical and
# are read, not asserted.

source("R/two_stage_correction.R")

set.seed(2027)
tolerance <- 1e-8

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
truth <- list(sigma_omega = 0.7, range_omega = 500, sigma_eta = 0.05,
              range_eta = 300, phi = 0.7, persistence = 1 / 0.3, tau = 0.3,
              sigma_p = 0.5, sigma_sy = 0.5)
rho <- 0.1


# design ------------------------------------------------------------------------

revisited <- sample.int(n_sites, 150)
pixel_years <- tibble(
  site = c(sample.int(n_sites, 700, replace = TRUE),
           sample(revisited, 800, replace = TRUE)),
  year = sample(t0:T_data, 1500, replace = TRUE)
) %>%
  distinct()
# studies: each pixel-year's, from 150 studies, a study's pixels in a year
# taken from one region (sites ordered by longitude), as surveys are
site_rank <- rank(sites_km[, 1], ties.method = "first")
pixel_years <- pixel_years %>%
  mutate(study = paste0("study ", pmin(150, ceiling(
    (site_rank[site] + runif(n(), -60, 60)) / n_sites * 150) %>% pmax(1))))
train <- pixel_years[sample.int(nrow(pixel_years), 2000, replace = TRUE), ]
train_sites <- unique(train$site)

held_out <- function(site, year, type, jitter = 0) {
  tibble(site = site, year = year, type = type,
         study = paste("new study", seq_along(site)),
         jitter_x = runif(length(site), -jitter, jitter),
         jitter_y = runif(length(site), -jitter, jitter))
}
new <- bind_rows(
  # sites jittered by up to ~50 km: new pixels, so fresh u, p and sy
  held_out(sample.int(n_sites, 400, replace = TRUE),
           sample((t0 + 1):T_data, 400, replace = TRUE), "new pixel", 35),
  # new years at revisited training pixels (by new studies): the truth shares
  # p with training
  held_out(sample(intersect(revisited, train_sites), 300, replace = TRUE),
           sample((t0 + 1):T_data, 300, replace = TRUE), "new year"),
  # training pixel-years: the truth shares u, p and sy with training
  distinct(train, site, year, study) %>%
    slice_sample(n = 300) %>%
    mutate(type = "training pixel-year", jitter_x = 0, jitter_y = 0),
  held_out(sample(train_sites, 500, replace = TRUE),
           sample((T_data + 1):(T_data + max_horizon), 500, replace = TRUE),
           "forecast")
)
train <- train %>%
  mutate(x_km = sites_km[site, 1], y_km = sites_km[site, 2], cell = site)
new <- new %>%
  mutate(x_km = sites_km[site, 1] + jitter_x,
         y_km = sites_km[site, 2] + jitter_y,
         cell = ifelse(type == "new pixel", -row_number(), site)) %>%
  filter(type == "training pixel-year" |
           !paste(cell, year) %in% paste(train$cell, train$year))

meshes <- build_correction_meshes(cbind(train$x_km, train$y_km),
                                  prediction_mask_coords())


# truth ---------------------------------------------------------------------------

sample_gmrf <- function(Q, n = 1) {
  sample_latent_deviation(Matrix::Cholesky(Matrix::forceSymmetric(Q),
                                           perm = TRUE, LDL = FALSE),
                          nrow(Q), n)
}
w_true <- sample_gmrf(matern_precision_r(correction_fem(meshes$omega),
                                         sqrt(8) / truth$range_omega,
                                         truth$sigma_omega))[, 1]
years_all <- (t0 + 1):(T_data + max_horizon)
eta_true <- sample_gmrf(matern_precision_r(correction_fem(meshes$xi),
                                           sqrt(8) / truth$range_eta,
                                           truth$sigma_eta), length(years_all))
for (t in 2:length(years_all)) {
  eta_true[, t] <- truth$phi * eta_true[, t - 1] +
    sqrt(1 - truth$phi ^ 2) * eta_true[, t]
}
x_true <- t(apply(eta_true, 1, cumsum))

# mostly susceptible, trending towards resistance, with a strong spatial
# gradient: true mortality spans ~0.1 to > 0.999
m_fun <- function(x_km, y_km, year) {
  4 - 0.12 * (year - t0) + 2 * sin(x_km / 700) * cos(y_km / 900)
}
u_keys <- unique(c(paste(train$cell, train$year), paste(new$cell, new$year)))
u_true <- setNames(rnorm(length(u_keys), 0, truth$tau), u_keys)
p_keys <- as.character(unique(c(train$cell, new$cell)))
p_true <- setNames(rnorm(length(p_keys), 0, truth$sigma_p), p_keys)
sy_keys <- unique(c(paste(train$study, train$year),
                    paste(new$study, new$year)))
sy_true <- setNames(rnorm(length(sy_keys), 0, truth$sigma_sy), sy_keys)

# the target m + omega + xi, and lambda = target + u + p + sy
smooth_truth <- function(d) {
  coords <- cbind(d$x_km, d$y_km)
  A_xi <- mesh_basis(meshes$xi, coords)
  xi <- vapply(seq_len(nrow(d)), function(i) {
    if (d$year[i] <= t0) 0 else sum(A_xi[i, ] * x_true[, d$year[i] - t0])
  }, numeric(1))
  m_fun(d$x_km, d$y_km, d$year) +
    as.vector(mesh_basis(meshes$omega, coords) %*% w_true) + xi
}
noise_truth <- function(d) {
  u_true[paste(d$cell, d$year)] + p_true[as.character(d$cell)] +
    sy_true[paste(d$study, d$year)]
}
train <- train %>%
  mutate(m = m_fun(x_km, y_km, year),
         lambda = smooth_truth(.) + noise_truth(.),
         mosquito_number = sample(c(rep(100, 6), 20:150), n(), replace = TRUE),
         died = rbinom(n(), mosquito_number,
                       rbeta(n(), plogis(lambda) * (1 - rho) / rho,
                             (1 - plogis(lambda)) * (1 - rho) / rho)),
         rho = rho)
# a new assay's lambda draws fresh noise (at a training pixel-year, the
# realised ones), the same as the training assays there
new <- new %>% mutate(m = m_fun(x_km, y_km, year), target = smooth_truth(.),
                      lambda = target + noise_truth(.))
cat(sprintf(paste("simulated %i assays (%.0f%% at 100%%, %.0f%% at 0%%),",
                  "%i pixels; meshes %i / %i nodes; %i held-out points\n"),
            nrow(train), 100 * mean(train$died == train$mosquito_number),
            100 * mean(train$died == 0), n_distinct(train$cell),
            meshes$omega$n, meshes$xi$n, nrow(new)))


# fits ------------------------------------------------------------------------------

fit_a <- fit_stage_a(train, t0 = t0, T = T_data, meshes = meshes)
fit <- fit_correction(train, t0 = t0, T = T_data, meshes = meshes)
cat(sprintf(paste("stage A: convergence %i, %.0f s. Final: convergence %i,",
                  "max gradient %.1e, %i CCD points; stage A %.0f s, Laplace",
                  "%.0f s, CCD %.0f s\n"),
            fit_a$opt$convergence, fit_a$timings[["total"]],
            fit$opt$convergence, fit$max_gradient, nrow(fit$ccd$theta),
            fit$timings[["stage_a"]], fit$timings[["laplace"]],
            fit$timings[["ccd"]]))

cat("\n1. stage-A mode vs H^-1 A' D (z - m): max |diff| = ")
rhs <- Matrix::crossprod(fit_a$A_latent, fit_a$precision_obs *
                           (fit_a$tmb_data$z - fit_a$m_ref))
check_1 <- max(abs(as.vector(Matrix::solve(fit_a$H_chol, rhs)) - fit_a$mode))
cat(sprintf("%.2e\n", check_1))
stopifnot(check_1 < tolerance)

m_perturbed <- train$m + 0.3 * sin(train$x_km / 400) + rnorm(nrow(train), 0, 0.1)
refit_mode <- function(f, likelihood) {
  data <- f$tmb_data
  data$m <- m_perturbed
  data$likelihood <- likelihood
  obj <- correction_adfun(data, f$par_list, fix_hyper = TRUE)
  invisible(obj$fn(obj$par))
  obj$env$last.par[obj$env$random]
}
cat("2. stage-A cut-posterior shift vs refit with perturbed m: max |diff| = ")
shifted_a <- fit_a$mode + correction_mode_shift(fit_a, m_perturbed)[, 1]
check_2 <- max(abs(shifted_a - refit_mode(fit_a, 0L)))
cat(sprintf("%.2e (max |shift| %.3f)\n", check_2,
            max(abs(shifted_a - fit_a$mode))))
stopifnot(check_2 < tolerance)
shifted <- fit$mode + correction_mode_shift(fit, m_perturbed)[, 1]
cat(sprintf(paste("   final fit, linearised shift vs refit: max |diff| %.3f,",
                  "RMS %.4f (max |shift| %.3f)\n"),
            max(abs(shifted - refit_mode(fit, 1L))),
            sqrt(mean((shifted - refit_mode(fit, 1L)) ^ 2)),
            max(abs(shifted - fit$mode))))

cat("3. beta-binomial score at the final fit's mode: max |gradient| = ")
A <- fit$A_latent
lambda_hat <- fit$m_ref + as.vector(A %*% fit$mode)
p_hat <- plogis(lambda_hat)
c_rho <- (1 - rho) / rho
dloglik <- c_rho * p_hat * (1 - p_hat) *
  (digamma(train$died + c_rho * p_hat) -
     digamma(train$mosquito_number - train$died + c_rho * (1 - p_hat)) -
     digamma(c_rho * p_hat) + digamma(c_rho * (1 - p_hat)))
# Q theta = H theta - A' W A theta, with the W that H was built with
q_theta <- as.vector(fit$H %*% fit$mode) -
  as.vector(Matrix::crossprod(A, fit$precision_obs *
                                as.vector(A %*% fit$mode)))
score <- as.vector(Matrix::crossprod(A, dloglik)) - q_theta
cat(sprintf("%.2e (max |Q theta| %.2e)\n", max(abs(score)), max(abs(q_theta))))
stopifnot(max(abs(score)) < 1e-6 * max(1, abs(q_theta)))

cat("4. CCD weights on a Gaussian: max |mean|, max |cov - I| = ")
d <- ncol(fit$ccd$theta)
z <- ccd_design(d)
n_points <- nrow(z)
w <- c(1, rep(exp(d * 1.1 ^ 2 / 2) / ((n_points - 1) * (1.1 ^ 2 - 1)),
              n_points - 1)) * exp(-rowSums(z ^ 2) / 2)
w <- w / sum(w)
z_mean <- colSums(z * w)
check_4 <- max(abs(crossprod(z * sqrt(w)) - tcrossprod(z_mean) - diag(d)))
cat(sprintf("%.2e, %.2e\n", max(abs(z_mean)), check_4))
stopifnot(max(abs(z_mean)) < tolerance, check_4 < tolerance)

cat("\n5. hyperparameters\n")
ccd_quantile <- function(name, q) {
  x <- fit$ccd$hyper[[name]]
  o <- order(x)
  x[o][findInterval(q, cumsum(fit$ccd$weight[o])) + 1]
}
print(tibble(parameter = names(truth), truth = unlist(truth),
             stage_a = unlist(fit_a$hyper[names(truth)]),
             theta_hat = unlist(fit$hyper[names(truth)]),
             ccd_mean = vapply(names(truth), function(x) {
               sum(fit$ccd$hyper[[x]] * fit$ccd$weight)
             }, numeric(1)),
             ccd_lower = vapply(names(truth), ccd_quantile, numeric(1), 0.025),
             ccd_upper = vapply(names(truth), ccd_quantile, numeric(1),
                                0.975)), digits = 3)
p_fitted <- fit$mode[fit$blocks$p]
cat(sprintf("cor(p-hat, p) = %.2f over %i pixels\n",
            cor(p_fitted, p_true[fit$iid_levels$p]), length(p_fitted)))
sy_fitted <- fit$mode[fit$blocks$sy]
cat(sprintf("cor(sy-hat, sy) = %.2f over %i study-years\n",
            cor(sy_fitted, sy_true[fit$iid_levels$sy]), length(sy_fitted)))


# coverage ------------------------------------------------------------------------

coverage <- function(draws, stage, truth = new$lambda) {
  groups <- c(split(seq_len(nrow(new)), new$type),
              list(`all, true mortality > 0.95` = which(plogis(truth) > 0.95),
                   all = seq_len(nrow(new))))
  bind_rows(lapply(names(groups), function(g) {
    i <- groups[[g]]
    d <- draws[, i, drop = FALSE]
    covered <- function(level) {
      lower <- apply(d, 2, quantile, (1 - level) / 2)
      upper <- apply(d, 2, quantile, 1 - (1 - level) / 2)
      mean(truth[i] >= lower & truth[i] <= upper)
    }
    tibble(stage = stage, group = g, n = length(i),
           cover_50 = covered(0.5), cover_95 = covered(0.95),
           bias = mean(colMeans(d) - truth[i]),
           rmse = sqrt(mean((colMeans(d) - truth[i]) ^ 2)),
           mean_sd = mean(apply(d, 2, sd)))
  }))
}

cat("\n6. coverage of a new assay's lambda (1000 draws)\n")
set.seed(1)
draws_a <- predict_correction(fit_a, new, n_draws = 1000)
set.seed(1)
draws_b <- predict_correction(fit, new, n_draws = 1000)
print(as.data.frame(bind_rows(coverage(draws_a, "A"),
                              coverage(draws_b, "final")) %>%
                      arrange(group, stage)),
      digits = 3, row.names = FALSE)
forecast <- new$type == "forecast"
lower <- apply(draws_b[, forecast], 2, quantile, 0.025)
upper <- apply(draws_b[, forecast], 2, quantile, 0.975)
cat("final fit's 95% coverage of forecasts by horizon:\n")
print(round(tapply(new$lambda[forecast] >= lower &
                     new$lambda[forecast] <= upper,
                   new$year[forecast] - T_data, mean), 3))


cat("coverage of the target m + omega + xi by the final fit's draws without",
    "noise:\n")
years <- sort(unique(new$year))
smooth_draws <- new$m + project_correction(
  fit, correction_node_draws(fit, years, 2000), new)
print(as.data.frame(coverage(t(smooth_draws), "final", new$target)),
      digits = 3, row.names = FALSE)


# mean fields ------------------------------------------------------------------------

cat("\n7. posterior-mean omega + xi vs the mean of 2000 draws: ")
smooth_mean <- new$m + project_correction(
  fit, correction_node_draws(fit, years, mean = TRUE), new)[, 1]
z <- (rowMeans(smooth_draws) - smooth_mean) /
  (apply(smooth_draws, 1, sd) / sqrt(2000))
cat(sprintf("mean z %.2f, sd z %.2f, max |z| %.1f\n", mean(z), sd(z),
            max(abs(z))))
