# Simulation check of the final two-stage model (#21, R/two_stage_correction.R
# and R/two_stage_pql.R).
#
#   Rscript R/check_two_stage.R
#
# Data are simulated from the final model, omega + xi + u + p, on the
# Lambda-cyhalothrin bioassay sites with the final model's meshes, as
# beta-binomial counts. The dynamical prediction m puts much of the latent
# lambda near mortality 1 (and some near 0), as in the real data, where the
# empirical logit of stage A is pinned by its +0.5 correction. A subset of
# sites is revisited often, so that p is identified. The truth is simulated on
# the fitting meshes, so this checks the implementation, not discretisation.
#
#   1. the stage-A mode is H^-1 A' D (z - m): the Hessian identity that PQL
#      and the cut-posterior shift rely on;
#   2. the cut-posterior shift H^-1 A' D (m_ref - m_new) of the stage-B fit
#      against refitting its latent mode in TMB with the perturbed m and the
#      hyperparameters fixed (exact for a Gaussian latent model);
#   3. the stage-B mode solves the penalised quasi-score equation, and the
#      bias of the fitted lambda at the training assays, stage A vs B;
#   4. recovery of the hyperparameters, sigma_p in particular;
#   5. coverage, at new pixels, new years at training pixels, training
#      pixel-years and 1-5 year forecasts: of a new assay's lambda by
#      predict_correction() (fresh u and p), stage A vs B, and of the target
#      m + omega + xi by its stage-B draws without noise;
#   6. the posterior-mean fields (correction_node_draws(mean = TRUE)) against
#      the mean of the draws, for the smooth terms the maps show.
#
# Checks 1-3 are identities, met to ~1e-11 when the code is right; the script
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
truth <- list(sigma_omega = 0.7, range_omega = 500, sigma_eta = 0.15,
              range_eta = 300, phi = 0.7, tau = 0.3, sigma_p = 0.5)
rho <- 0.1


# design ------------------------------------------------------------------------

revisited <- sample.int(n_sites, 150)
pixel_years <- tibble(
  site = c(sample.int(n_sites, 700, replace = TRUE),
           sample(revisited, 800, replace = TRUE)),
  year = sample(t0:T_data, 1500, replace = TRUE)
) %>%
  distinct()
train <- pixel_years[sample.int(nrow(pixel_years), 2000, replace = TRUE), ]
train_sites <- unique(train$site)

held_out <- function(site, year, type, jitter = 0) {
  tibble(site = site, year = year, type = type,
         jitter_x = runif(length(site), -jitter, jitter),
         jitter_y = runif(length(site), -jitter, jitter))
}
new <- bind_rows(
  # sites jittered by up to ~50 km: new pixels, so fresh u and p
  held_out(sample.int(n_sites, 400, replace = TRUE),
           sample((t0 + 1):T_data, 400, replace = TRUE), "new pixel", 35),
  # new years at revisited training pixels: the truth shares p with training
  held_out(sample(intersect(revisited, train_sites), 300, replace = TRUE),
           sample((t0 + 1):T_data, 300, replace = TRUE), "new year"),
  # training pixel-years: the truth shares u and p with training
  distinct(train, site, year) %>%
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

# the target m + omega + xi, and lambda = target + u + p
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
  u_true[paste(d$cell, d$year)] + p_true[as.character(d$cell)]
}
train <- train %>%
  mutate(m = m_fun(x_km, y_km, year),
         lambda = smooth_truth(.) + noise_truth(.),
         mosquito_number = sample(c(rep(100, 6), 20:150), n(), replace = TRUE),
         died = rbinom(n(), mosquito_number,
                       rbeta(n(), plogis(lambda) * (1 - rho) / rho,
                             (1 - plogis(lambda)) * (1 - rho) / rho)),
         rho = rho)
# a new assay's lambda draws fresh u and p (at a training pixel-year, the
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
fit_b <- fit_correction_pql(fit_a, train)
sb <- fit_b$stage_b
cat(sprintf(paste("stage A: convergence %i, %.0f s. Stage B: %i passes, RMS",
                  "move %.3f, refit %s (convergence %i), %i more passes,",
                  "%.0f s\n"),
            fit_a$opt$convergence, fit_a$timings[["total"]], sb$passes_first,
            sb$rms_move, sb$refit, fit_b$opt$convergence, sb$passes_second,
            sb$time_pql))

cat("\n1. stage-A mode vs H^-1 A' D (z - m): max |diff| = ")
rhs <- Matrix::crossprod(fit_a$A_latent, fit_a$precision_obs *
                           (fit_a$tmb_data$z - fit_a$m_ref))
check_1 <- max(abs(as.vector(Matrix::solve(fit_a$H_chol, rhs)) - fit_a$mode))
cat(sprintf("%.2e\n", check_1))
stopifnot(check_1 < tolerance)

cat("2. stage-B cut-posterior shift vs refit with perturbed m: max |diff| = ")
m_perturbed <- train$m + 0.3 * sin(train$x_km / 400) + rnorm(nrow(train), 0, 0.1)
shifted <- fit_b$mode + correction_mode_shift(fit_b, m_perturbed)[, 1]
data <- fit_b$tmb_data
data$m <- m_perturbed
obj <- correction_adfun(data, fit_b$par_list, fix_hyper = TRUE)
invisible(obj$fn(obj$par))
check_2 <- max(abs(shifted - obj$env$last.par[obj$env$random]))
cat(sprintf("%.2e (max |shift| %.3f)\n", check_2,
            max(abs(shifted - fit_b$mode))))
stopifnot(check_2 < tolerance)

cat("3. penalised quasi-score at the stage-B mode: max |gradient| = ")
A <- fit_b$A_latent
lambda_b <- fit_b$m_ref + as.vector(A %*% fit_b$mode)
# Q theta = H theta - A' D A theta, with the D that H was built with
q_theta <- as.vector(fit_b$H %*% fit_b$mode) -
  as.vector(Matrix::crossprod(A, fit_b$precision_obs *
                                as.vector(A %*% fit_b$mode)))
score <- as.vector(Matrix::crossprod(
  A, (train$died - train$mosquito_number * plogis(lambda_b)) /
    (1 + (train$mosquito_number - 1) * rho))) - q_theta
cat(sprintf("%.2e (max |Q theta| %.2e)\n", max(abs(score)), max(abs(q_theta))))
# relative to the size of the penalty term it balances
stopifnot(max(abs(score)) < tolerance * max(1, abs(q_theta)))
cat("fitted lambda minus the truth at the training assays, by true mortality:\n")
errors <- tibble(bin = cut(plogis(train$lambda),
                           c(0, 0.05, 0.5, 0.9, 0.95, 0.99, 1)),
                 a = sb$lambda_a - train$lambda,
                 b = lambda_b - train$lambda)
print(bind_rows(errors, mutate(errors, bin = "all")) %>%
        group_by(bin) %>%
        summarise(n = n(), bias_a = mean(a), bias_b = mean(b),
                  rmse_a = sqrt(mean(a ^ 2)), rmse_b = sqrt(mean(b ^ 2))),
      digits = 3)

cat("\n4. hyperparameters\n")
print(tibble(parameter = names(truth), truth = unlist(truth),
             stage_a = unlist(fit_a$hyper[names(truth)]),
             stage_b = unlist(fit_b$hyper[names(truth)])), digits = 3)
p_hat <- fit_b$mode[fit_b$blocks$p]
cat(sprintf("cor(p-hat, p) = %.2f over %i pixels\n",
            cor(p_hat, p_true[fit_b$iid_levels$p]), length(p_hat)))


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

cat("\n5. coverage of a new assay's lambda (1000 draws)\n")
set.seed(1)
draws_a <- predict_correction(fit_a, new, n_draws = 1000)
set.seed(1)
draws_b <- predict_correction(fit_b, new, n_draws = 1000)
print(as.data.frame(bind_rows(coverage(draws_a, "A"), coverage(draws_b, "B")) %>%
                      arrange(group, stage)),
      digits = 3, row.names = FALSE)
forecast <- new$type == "forecast"
lower <- apply(draws_b[, forecast], 2, quantile, 0.025)
upper <- apply(draws_b[, forecast], 2, quantile, 0.975)
cat("stage B 95% coverage of forecasts by horizon:\n")
print(round(tapply(new$lambda[forecast] >= lower &
                     new$lambda[forecast] <= upper,
                   new$year[forecast] - T_data, mean), 3))


cat("coverage of the target m + omega + xi by stage-B draws without noise:\n")
years <- sort(unique(new$year))
smooth_draws <- new$m + project_correction(
  fit_b, correction_node_draws(fit_b, years, 2000), new)
print(as.data.frame(coverage(t(smooth_draws), "B", new$target)), digits = 3,
      row.names = FALSE)


# mean fields ------------------------------------------------------------------------

cat("\n6. posterior-mean omega + xi vs the mean of 2000 draws: ")
smooth_mean <- new$m + project_correction(
  fit_b, correction_node_draws(fit_b, years, mean = TRUE), new)[, 1]
z <- (rowMeans(smooth_draws) - smooth_mean) /
  (apply(smooth_draws, 1, sd) / sqrt(2000))
cat(sprintf("mean z %.2f, sd z %.2f, max |z| %.1f\n", mean(z), sd(z),
            max(abs(z))))
