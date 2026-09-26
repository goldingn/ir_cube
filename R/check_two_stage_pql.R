# Simulation check of stage B (PQL, R/two_stage_pql.R) against stage A.
#
#   Rscript R/check_two_stage_pql.R
#
# Data are simulated from the correction model on the Lambda-cyhalothrin sites
# (as in R/check_two_stage_correction.R), but the assays are beta-binomial
# counts rather than Gaussian empirical logits, and the dynamical prediction m
# puts much of the latent lambda near p = 1 (and some near 0), as in the real
# data, where 29% of training assays read 100%. There the empirical logit is
# pinned by the +0.5 correction, so stage A is biased towards 0. Checks:
#
#   1. the stage-A mode is H^-1 A' D (z - m), i.e. the Gaussian priors are zero
#      mean and the Hessian identity stage B relies on holds;
#   2. the stage-B fixed point solves the penalised quasi-score equation;
#   3. bias and RMSE of the fitted latent lambda at the training assays, by
#      bin of true p, stage A vs stage B;
#   4. coverage of the true lambda by prediction intervals at held-out
#      pixel-years (interpolation), at training pixel-years and at 1-5 year
#      forecasts, stage A vs stage B, overall and where true p > 0.95.

source("R/two_stage_correction.R")
source("R/two_stage_pql.R")

set.seed(2027)

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

truth <- list(sigma_omega = 0.7, range_omega = 500, sigma_eta = 0.15,
              range_eta = 300, phi = 0.7, tau = 0.4)
rho <- 0.1

pixel_years <- tibble(
  site = sample.int(n_sites, n_pixel_years, replace = TRUE),
  year = sample(t0:T_data, n_pixel_years, replace = TRUE)
) %>%
  distinct()
train <- pixel_years[sample.int(nrow(pixel_years), n_obs, replace = TRUE), ]

n_interp <- 400
interp <- tibble(site = sample.int(n_sites, n_interp, replace = TRUE),
                 year = sample((t0 + 1):T_data, n_interp, replace = TRUE),
                 jitter_x = runif(n_interp, -35, 35),
                 jitter_y = runif(n_interp, -35, 35),
                 type = "interpolation")
insample <- train %>%
  distinct(site, year) %>%
  slice_sample(n = 400) %>%
  mutate(jitter_x = 0, jitter_y = 0, type = "training pixel-year")
n_forecast <- 500
forecast <- tibble(site = sample(unique(train$site), n_forecast, replace = TRUE),
                   year = sample((T_data + 1):(T_data + max_horizon),
                                 n_forecast, replace = TRUE),
                   jitter_x = 0, jitter_y = 0, type = "forecast")

add_coords <- function(df) {
  df %>% mutate(
    x_km = sites_km[site, 1] + if ("jitter_x" %in% names(df)) jitter_x else 0,
    y_km = sites_km[site, 2] + if ("jitter_y" %in% names(df)) jitter_y else 0)
}
train <- add_coords(train) %>% mutate(cell = site)
new <- bind_rows(interp, insample, forecast) %>%
  add_coords() %>%
  mutate(cell = ifelse(type == "interpolation", -row_number(), site))

mesh <- build_correction_mesh(cbind(train$x_km, train$y_km))
mesh_xi <- build_correction_mesh(cbind(train$x_km, train$y_km),
                                 max_nodes = 600)
fem <- correction_fem(mesh)
fem_xi <- correction_fem(mesh_xi)

sample_gmrf <- function(Q, n = 1) {
  chol_Q <- Matrix::Cholesky(Matrix::forceSymmetric(Q), perm = TRUE,
                             LDL = FALSE)
  sample_latent_deviation(chol_Q, nrow(Q), n)
}
kappa <- function(range) sqrt(8) / range
w_true <- sample_gmrf(matern_precision_r(fem, kappa(truth$range_omega),
                                         truth$sigma_omega))[, 1]
Q_eta <- matern_precision_r(fem_xi, kappa(truth$range_eta), truth$sigma_eta)
eps <- sample_gmrf(Q_eta, length(years_all))
eta_true <- eps
for (t in 2:length(years_all)) {
  eta_true[, t] <- truth$phi * eta_true[, t - 1] +
    sqrt(1 - truth$phi ^ 2) * eps[, t]
}
x_true <- t(apply(eta_true, 1, cumsum))

# mostly susceptible, trending towards resistance, with a strong spatial
# gradient: true p spans ~0.1 to >0.999
m_fun <- function(x_km, y_km, year) {
  4 - 0.12 * (year - t0) + 2 * sin(x_km / 700) * cos(y_km / 900)
}
latent_truth <- function(df, u) {
  A <- mesh_basis(mesh, cbind(df$x_km, df$y_km))
  A_xi <- mesh_basis(mesh_xi, cbind(df$x_km, df$y_km))
  omega <- as.vector(A %*% w_true)
  xi <- numeric(nrow(df))
  for (i in which(df$year > t0)) {
    xi[i] <- sum(A_xi[i, ] * x_true[, df$year[i] - t0])
  }
  m_fun(df$x_km, df$y_km, df$year) + omega + xi + u
}
all_keys <- unique(c(paste(train$cell, train$year), paste(new$cell, new$year)))
u_true <- setNames(rnorm(length(all_keys), 0, truth$tau), all_keys)

rbetabinom_sim <- function(n, p, rho) {
  a <- p * (1 - rho) / rho
  b <- (1 - p) * (1 - rho) / rho
  rbinom(length(n), n, rbeta(length(n), a, b))
}
train <- train %>%
  mutate(m = m_fun(x_km, y_km, year),
         lambda = latent_truth(., u_true[paste(cell, year)]),
         mosquito_number = sample(c(rep(100, 6), 20:150), n(), replace = TRUE),
         died = rbetabinom_sim(mosquito_number, plogis(lambda), rho),
         rho = rho)
new <- new %>%
  mutate(m = m_fun(x_km, y_km, year),
         lambda = latent_truth(., u_true[paste(cell, year)]))

cat(sprintf(paste("simulated %i assays: %.1f%% at 100%%, %.1f%% at 0%%;",
                  "true p > 0.95 at %.1f%%; meshes %i / %i nodes\n"),
            nrow(train), 100 * mean(train$died == train$mosquito_number),
            100 * mean(train$died == 0), 100 * mean(plogis(train$lambda) > 0.95),
            mesh$n, mesh_xi$n))


# fits ------------------------------------------------------------------------

fit_a <- fit_correction(train, variant = "omega_xi_u", t0 = t0, T = T_data,
                        mesh = mesh, mesh_xi = mesh_xi)
cat(sprintf("stage A: convergence %i, %.1f s\n", fit_a$opt$convergence,
            fit_a$timings["total"]))

cat("\n1. stage-A mode vs H^-1 A' D (z - m): max |diff| = ")
rhs <- Matrix::crossprod(fit_a$A_latent,
                         fit_a$precision_obs * (fit_a$tmb_data$z - fit_a$m_ref))
cat(sprintf("%.2e\n", max(abs(as.vector(Matrix::solve(fit_a$H_chol, rhs)) -
                                fit_a$mode))))

time_b <- system.time(
  fit_b <- fit_correction_pql(fit_a, train, verbose = TRUE)
)[["elapsed"]]
sb <- fit_b$stage_b
cat(sprintf(paste("stage B: %i passes (converged %s, damped %s); RMS move",
                  "%.3f, max %.3f; refit %s; second PQL %s passes; %.1f s\n"),
            sb$passes_first, sb$converged_first, sb$damped_first, sb$rms_move,
            sb$max_move, sb$refit, sb$passes_second, time_b))
print(sb$history_first, digits = 3)
if (sb$refit) print(sb$history_second, digits = 3)

cat("\nhyperparameters: truth, stage A, stage B\n")
print(tibble(parameter = names(truth), truth = unlist(truth),
             stage_a = sapply(names(truth), function(p) fit_a$hyper[[p]]),
             stage_b = sapply(names(truth), function(p) fit_b$hyper[[p]])),
      digits = 3)

cat("\n2. penalised quasi-score at the stage-B mode: max |gradient| = ")
A <- fit_b$A_latent
lambda_b <- fit_b$m_ref + as.vector(A %*% fit_b$mode)
phi_obs <- 1 + (train$mosquito_number - 1) * rho
# Q theta = H theta - A' D A theta, with the D that H was built with
q_theta <- as.vector(fit_b$H %*% fit_b$mode) -
  as.vector(Matrix::crossprod(A, fit_b$precision_obs *
                                as.vector(A %*% fit_b$mode)))
score <- as.vector(Matrix::crossprod(
  A, (train$died - train$mosquito_number * plogis(lambda_b)) / phi_obs)) -
  q_theta
cat(sprintf("%.2e (max |Q theta| = %.2e)\n", max(abs(score)),
            max(abs(q_theta))))


# bias at the training assays ------------------------------------------------------

cat("\n3. fitted lambda at the training assays minus the truth, by true p\n")
lambda_a <- sb$lambda_a
bias_table <- tibble(p_bin = cut(plogis(train$lambda),
                                 c(0, 0.05, 0.5, 0.9, 0.95, 0.99, 1)),
                     err_a = lambda_a - train$lambda,
                     err_b = lambda_b - train$lambda,
                     full = train$died == train$mosquito_number) %>%
  group_by(p_bin) %>%
  summarise(n = n(), at_100 = mean(full),
            bias_a = mean(err_a), bias_b = mean(err_b),
            rmse_a = sqrt(mean(err_a ^ 2)), rmse_b = sqrt(mean(err_b ^ 2)),
            .groups = "drop")
print(bias_table, digits = 3)
extreme <- plogis(train$lambda) > 0.95
cat(sprintf("true p > 0.95 (n = %i): bias A %+.3f, B %+.3f; RMSE A %.3f, B %.3f\n",
            sum(extreme), mean((lambda_a - train$lambda)[extreme]),
            mean((lambda_b - train$lambda)[extreme]),
            sqrt(mean((lambda_a - train$lambda)[extreme] ^ 2)),
            sqrt(mean((lambda_b - train$lambda)[extreme] ^ 2))))
cat(sprintf("all assays: RMSE A %.3f, B %.3f\n",
            sqrt(mean((lambda_a - train$lambda) ^ 2)),
            sqrt(mean((lambda_b - train$lambda) ^ 2))))


# prediction coverage ----------------------------------------------------------------

coverage_rows <- function(draws, df, stage) {
  groups <- c(split(seq_len(nrow(df)), df$type),
              list(`all, true p > 0.95` = which(plogis(df$lambda) > 0.95),
                   all = seq_len(nrow(df))))
  bind_rows(lapply(names(groups), function(g) {
    i <- groups[[g]]
    d <- draws[, i, drop = FALSE]
    covered <- function(level) {
      lower <- apply(d, 2, quantile, (1 - level) / 2)
      upper <- apply(d, 2, quantile, 1 - (1 - level) / 2)
      mean(df$lambda[i] >= lower & df$lambda[i] <= upper)
    }
    tibble(stage = stage, group = g, n = length(i),
           cover_50 = covered(0.5), cover_80 = covered(0.8),
           cover_95 = covered(0.95),
           bias = mean(colMeans(d) - df$lambda[i]),
           rmse = sqrt(mean((colMeans(d) - df$lambda[i]) ^ 2)))
  }))
}

cat("\n4. coverage of the true lambda (1000 draws each)\n")
set.seed(1)
draws_a <- predict_correction(fit_a, new, n_draws = 1000)
set.seed(1)
time_pred <- system.time(
  draws_b <- predict_correction(fit_b, new, n_draws = 1000)
)[["elapsed"]]
coverage <- bind_rows(coverage_rows(draws_a, new, "A"),
                      coverage_rows(draws_b, new, "B")) %>%
  arrange(group, stage)
print(as.data.frame(coverage), digits = 3, row.names = FALSE)
cat(sprintf("stage-B prediction at %i points: %.1f s\n", nrow(new), time_pred))
