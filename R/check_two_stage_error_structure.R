# Checks of the optional error-structure terms of the two-stage correction
# (#21; doc/two_stage_plan.md, "Error structure"): the static pixel effect p
# and the survey effect s.
#
#   Rscript R/check_two_stage_error_structure.R
#
# 1. Defaults unchanged: the template and R code as committed before these
#    terms (ERROR_STRUCTURE_REF) and the current code, fitted to the
#    same simulated data with the defaults (no p, no s), give the same
#    objective, hyperparameters, latent mode, stage-B (PQL) fit and predictive
#    draws. The old code is taken from `git show <ref>:<file>`
#    (ERROR_STRUCTURE_REF, default HEAD) and compiled under another name.
# 2. Recovery: data simulated with omega + xi + u + p + s; the +p+s fit
#    recovers sigma_p, sigma_s and the other hyperparameters, the cut-posterior
#    shift still matches a refit, and prediction intervals cover:
#    - the target lambda = m + omega + xi + u + p (survey = "none") at held-out
#      points: new pixels, and new years at training pixels (where p is known);
#    - a new assay's latent lambda + s (survey = "fresh") in new surveys.
#    The model without p and s is fitted to the same data for comparison.

source("R/two_stage_correction.R")
source("R/two_stage_pql.R")

scratch <- Sys.getenv("ERROR_STRUCTURE_SCRATCH", tempdir())
ref <- Sys.getenv("ERROR_STRUCTURE_REF", "HEAD")
set.seed(2027)


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
              range_eta = 300, phi = 0.7, tau = 0.2, sigma_p = 0.5,
              sigma_s = 0.4)
rho <- 0.15

# surveys: a spatial cluster of sites (k-means, 80 clusters) in one year, so a
# survey spans several pixels and a pixel is visited by several surveys
site_cluster <- kmeans(sites_km, centers = 80, nstart = 5)$cluster
# a subset of sites is revisited often, so p is identified from repeat years
revisited <- sample.int(n_sites, 150)
pixel_years <- tibble(
  site = c(sample.int(n_sites, 700, replace = TRUE),
           sample(revisited, 800, replace = TRUE)),
  year = sample(t0:T_data, 1500, replace = TRUE)
) %>%
  distinct()
train <- pixel_years[sample.int(nrow(pixel_years), 2000, replace = TRUE), ] %>%
  mutate(x_km = sites_km[site, 1], y_km = sites_km[site, 2], cell = site,
         survey = paste(site_cluster[site], year))

# held-out: new pixels (sites jittered, fresh p), and new years at revisited
# training pixels (p known), both in surveys not seen in training (odd
# cluster-years are relabelled as new surveys)
n_new <- 600
new <- bind_rows(
  tibble(site = sample.int(n_sites, n_new / 2, replace = TRUE),
         year = sample((t0 + 1):T_data, n_new / 2, replace = TRUE),
         jitter_x = runif(n_new / 2, -35, 35),
         jitter_y = runif(n_new / 2, -35, 35),
         type = "new pixel"),
  tibble(site = sample(intersect(revisited, train$site), n_new / 2,
                       replace = TRUE),
         year = sample((t0 + 1):T_data, n_new / 2, replace = TRUE),
         jitter_x = 0, jitter_y = 0,
         type = "training pixel, new year")
) %>%
  mutate(x_km = sites_km[site, 1] + jitter_x,
         y_km = sites_km[site, 2] + jitter_y,
         cell = ifelse(type == "new pixel", -row_number(), site),
         survey = paste("new", site_cluster[site], year))
# new-year rows must not reuse a training pixel-year (u would be known)
new <- new %>% filter(!paste(cell, year) %in% paste(train$cell, train$year))

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
m_fun <- function(x_km, y_km, year) {
  -1 + 0.12 * (year - t0) + 0.8 * sin(x_km / 700) * cos(y_km / 900)
}
all_cells <- unique(c(train$cell, new$cell))
p_true <- setNames(rnorm(length(all_cells), 0, truth$sigma_p), all_cells)
all_py <- unique(c(paste(train$cell, train$year), paste(new$cell, new$year)))
u_true <- setNames(rnorm(length(all_py), 0, truth$tau), all_py)
all_surveys <- unique(c(train$survey, new$survey))
s_true <- setNames(rnorm(length(all_surveys), 0, truth$sigma_s), all_surveys)

# the target lambda (no survey effect) and the assay latent lambda + s
latent_truth <- function(d) {
  A <- mesh_basis(mesh, cbind(d$x_km, d$y_km))
  A_xi <- mesh_basis(mesh_xi, cbind(d$x_km, d$y_km))
  xi <- vapply(seq_len(nrow(d)), function(i) {
    if (d$year[i] <= t0) 0 else sum(A_xi[i, ] * x_true[, d$year[i] - t0])
  }, numeric(1))
  m_fun(d$x_km, d$y_km, d$year) + as.vector(A %*% w_true) + xi +
    u_true[paste(d$cell, d$year)] + p_true[as.character(d$cell)]
}
train <- train %>%
  mutate(m = m_fun(x_km, y_km, year),
         lambda = latent_truth(.),
         lambda_assay = lambda + s_true[survey],
         mosquito_number = sample(20:100, n(), replace = TRUE),
         died = rbinom(n(), mosquito_number, plogis(lambda_assay)),
         rho = rho,
         v = empirical_logit(died, mosquito_number, rho)$v,
         z = lambda_assay + rnorm(n(), 0, sqrt(v)))
new <- new %>%
  mutate(m = m_fun(x_km, y_km, year),
         lambda = latent_truth(.),
         lambda_assay = lambda + s_true[survey])
cat(sprintf(paste("simulated: %i assays, %i pixels (%i with >1 year),",
                  "%i surveys; %i held-out points\n"),
            nrow(train), n_distinct(train$cell),
            sum(table(unique(train[, c("cell", "year")])$cell) > 1),
            n_distinct(train$survey), nrow(new)))


# 1. defaults unchanged ---------------------------------------------------------------

cat("\n1. defaults against the code at", ref, "\n")
old_dir <- file.path(scratch, "error_structure_head")
dir.create(old_dir, showWarnings = FALSE, recursive = TRUE)
git_show <- function(file, out) {
  status <- system2("git", c("show", paste0(ref, ":", file)), stdout = out)
  stopifnot(status == 0)
}
old_cpp <- file.path(old_dir, "two_stage_correction_head.cpp")
git_show("tmb/two_stage_correction.cpp", old_cpp)
git_show("R/two_stage_correction.R", file.path(old_dir, "correction.R"))
git_show("R/two_stage_pql.R", file.path(old_dir, "pql.R"))
old <- new.env()
sys.source(file.path(old_dir, "correction.R"), envir = old)
assign("correction_template", old_cpp, envir = old)
sys.source(file.path(old_dir, "pql.R"), envir = old)
stopifnot(!isTRUE(grepl("include_p", paste(readLines(old_cpp),
                                             collapse = "\n"))))

train_default <- select(train, -survey)
fits_default <- list(
  old = old$fit_correction(train_default, "omega_xi_u", t0 = t0, T = T_data,
                           mesh = mesh, mesh_xi = mesh_xi),
  new = fit_correction(train_default, "omega_xi_u", t0 = t0, T = T_data,
                       mesh = mesh, mesh_xi = mesh_xi)
)
fits_b <- list(
  old = old$fit_correction_pql(fits_default$old, train_default),
  new = fit_correction_pql(fits_default$new, train_default)
)
predict_both <- function(fits) {
  set.seed(1)
  draws_old <- old$predict_correction(fits$old, select(new, -survey),
                                      n_draws = 200)
  set.seed(1)
  draws_new <- predict_correction(fits$new, select(new, -survey),
                                  n_draws = 200)
  list(old = draws_old, new = draws_new)
}
draws_a <- predict_both(fits_default)
draws_b <- predict_both(fits_b)
differences <- c(
  objective_A = abs(fits_default$old$opt$objective -
                      fits_default$new$opt$objective),
  hyper_A = max(abs(unlist(fits_default$old$hyper) -
                      unlist(fits_default$new$hyper))),
  mode_A = max(abs(fits_default$old$mode - fits_default$new$mode)),
  hessian_A = max(abs(fits_default$old$H - fits_default$new$H)),
  hyper_B = max(abs(unlist(fits_b$old$hyper) - unlist(fits_b$new$hyper))),
  mode_B = max(abs(fits_b$old$mode - fits_b$new$mode)),
  draws_A = max(abs(draws_a$old - draws_a$new)),
  draws_B = max(abs(draws_b$old - draws_b$new))
)
print(signif(differences, 3))
stopifnot(identical(names(fits_default$old$blocks),
                    names(fits_default$new$blocks)),
          all(differences == 0))
cat("defaults reproduce the old code exactly\n")


# 2. recovery --------------------------------------------------------------------

cat("\n2. fits to data simulated with p and s\n")
fits <- list(
  omega_xi_u = fits_default$new,
  omega_xi_u_p = fit_correction(train, "omega_xi_u", t0 = t0, T = T_data,
                                mesh = mesh, mesh_xi = mesh_xi,
                                pixel_effect = TRUE),
  omega_xi_u_p_s = fit_correction(train, "omega_xi_u", t0 = t0, T = T_data,
                                  mesh = mesh, mesh_xi = mesh_xi,
                                  pixel_effect = TRUE, survey_effect = TRUE)
)
hyper_table <- tibble(parameter = names(truth), truth = unlist(truth))
for (name in names(fits)) {
  hyper_table[[name]] <- vapply(hyper_table$parameter, function(p) {
    value <- fits[[name]]$hyper[[p]]
    if (is.null(value)) NA_real_ else value
  }, numeric(1))
}
print(hyper_table, digits = 3)
cat("convergence:", vapply(fits, function(f) f$opt$convergence, 0), "\n")
cat("objective:", round(vapply(fits, function(f) f$opt$objective, 0), 2), "\n")

# latent effects against the truth at the training data
f <- fits$omega_xi_u_p_s
p_hat <- f$mode[f$blocks$p]
s_hat <- f$mode[f$blocks$s]
cat(sprintf("cor(p-hat, p) = %.2f over %i pixels, cor(s-hat, s) = %.2f over %i surveys\n",
            cor(p_hat, p_true[as.character(f$pixels$cell)]), length(p_hat),
            cor(s_hat, s_true[f$surveys$survey]), length(s_hat)))

# the cut-posterior shift against a refit with the hyperparameters fixed
m_perturbed <- train$m + 0.3 * sin(train$x_km / 400)
shifted <- f$mode + correction_mode_shift(f, m_perturbed)[, 1]
data <- f$tmb_data
data$m <- m_perturbed
obj <- correction_adfun(data, f$par_list, f$variant, fix_hyper = TRUE)
invisible(obj$fn(obj$par))
shift_error <- max(abs(shifted - obj$env$last.par[obj$env$random]))
cat(sprintf("cut-posterior shift vs refit, +p+s: max |difference| = %.2e\n",
            shift_error))
stopifnot(shift_error < 1e-8)

coverage <- function(draws, target, group) {
  out <- sapply(c(0.5, 0.95), function(level) {
    lower <- apply(draws, 2, quantile, (1 - level) / 2)
    upper <- apply(draws, 2, quantile, 1 - (1 - level) / 2)
    tapply(target >= lower & target <= upper, group, mean)
  })
  colnames(out) <- c("cover50", "cover95")
  cbind(out,
        rmse = tapply((colMeans(draws) - target) ^ 2, group,
                      function(x) sqrt(mean(x))),
        mean_sd = tapply(apply(draws, 2, sd), group, mean))
}
cat("\ncoverage of the target lambda (no s), survey = \"none\":\n")
for (name in names(fits)) {
  set.seed(3)
  draws <- predict_correction(fits[[name]], new, n_draws = 1000)
  cat(name, "\n")
  print(round(coverage(draws, new$lambda, new$type), 3))
}
cat("\ncoverage of a new assay's latent lambda + s:\n")
set.seed(3)
draws <- predict_correction(f, new, n_draws = 1000,
                            survey = c("none", "fresh"))
for (mode in names(draws)) {
  cat("+p+s, survey =", mode, "\n")
  print(round(coverage(draws[[mode]], new$lambda_assay, new$type), 3))
}
