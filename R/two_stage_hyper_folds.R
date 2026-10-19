# Experiments for #34: the two-stage correction with its hyperparameters
# estimated or integrated over in different ways, on one cross-validation fold
# and one insecticide type, scored on the fold's held-out assays.
#
#   Rscript R/two_stage_hyper_folds.R <experiment> <fold> <type index> <methods>
#
# e.g. Rscript R/two_stage_hyper_folds.R temporal_forecasting 2014 4 pql_plugin,la_ccd
#
# Needs the fold's stage-one cache (outputs/two_stage/stage_one__<key>.rds,
# from R/run_two_stage_folds.R). <methods> is a comma-separated list of
# plugin_once, the current algorithm (stage A, PQL, one re-estimation of
# theta, plug-in theta_hat; fit_correction()), or <latent>_<theta>:
#
#   latent  pql  PQL with theta re-estimated until it stops moving
#                (fit_correction(), max_refits = 10; #34, point 1); theta's
#                marginal is that of the final PQL working response;
#           la   the standalone Laplace fit (fit_correction_laplace(): stage
#                A, then theta by the Laplace approximation to the marginal
#                likelihood of the beta-binomial counts; #34, point 2)
#   theta   plugin  theta_hat;
#           gauss   100 draws from the Gaussian approximation to p(theta | y);
#           ccd     INLA's central composite design;
#           is      importance sampling, PSIS-smoothed, with Monte Carlo SEs:
#                   the reference (R/two_stage_hyper.R; #34, point 4).
#
# Each latent fit is saved once to outputs/two_stage/hyper/fits/
# <latent>__<key>__<k>.rds and reused. Each method saves
# outputs/two_stage/hyper/parts/<method>__<key>__<k>.rds: the held-out
# assays (test), their predictive draws of logit mortality (draws x assays,
# as in R/run_two_stage_folds.R: paired dynamical draw + omega + xi + fresh u
# and p), the points of theta and their weights, and a one-row summary.
# R/two_stage_hyper_scores.R scores them.
#
# is evaluates its proposal draws in IR_CUBE_IS_CORES forked processes. Run
# with OpenBLAS, as R/run_two_stage_folds.R.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 4)
experiment_name <- arguments[1]
fold_name <- arguments[2]
k <- as.integer(arguments[3])
methods <- strsplit(arguments[4], ",")[[1]]
all_methods <- c("plugin_once",
                 outer(c("pql", "la"), c("plugin", "gauss", "ccd", "is"),
                       paste, sep = "_"))
stopifnot(all(methods %in% all_methods))

source("R/two_stage_helpers.R")
source("R/two_stage_correction.R")
source("R/two_stage_hyper.R")

fold_key <- sprintf("%s__%s", experiment_name, fold_name)
hyper_dir <- "outputs/two_stage/hyper"
dir.create(file.path(hyper_dir, "fits"), showWarnings = FALSE,
           recursive = TRUE)
dir.create(file.path(hyper_dir, "parts"), showWarnings = FALSE)
base_file <- function(latent) {
  file.path(hyper_dir, "fits", sprintf("%s__%s__%i.rds", latent, fold_key, k))
}
part_file <- function(method) {
  file.path(hyper_dir, "parts", sprintf("%s__%s__%i.rds", method, fold_key, k))
}

cache <- readRDS(file.path("outputs/two_stage",
                           sprintf("stage_one__%s.rds", fold_key)))
types <- cache$types
stopifnot(!is.na(k), k >= 1, k <= length(types))
type <- types[k]
rho <- rho_for_record(tibble(insecticide_type = type), rho_lookup())
train_rows <- which(cache$training$type_id == k)
test_rows <- which(cache$test_df$type_id == k)
stopifnot(length(train_rows) > 0, length(test_rows) > 0)
m_draws_train <- cache$logit_train[, train_rows, drop = FALSE]
m_draws_test <- cache$logit_test[, test_rows, drop = FALSE]
train_k <- mutate(cache$training[train_rows, ], m = colMeans(m_draws_train),
                  rho = rho)
test <- cache$test_df[test_rows, c("insecticide_type", "insecticide_class",
                                   "country_name", "year_start", "cell",
                                   "longitude", "latitude", "died",
                                   "mosquito_number")]
test_k <- with(test, tibble(lon = longitude, lat = latitude,
                            year = year_start, cell = cell))
t0 <- cache$t0
rm(cache)
invisible(gc())
n_draws <- nrow(m_draws_train)
report("%s %s: %i training, %i held-out assays; methods %s", fold_key, type,
       nrow(train_k), nrow(test_k), paste(methods, collapse = ", "))

meshes <- build_correction_meshes(coords_km(train_k), prediction_mask_coords())

# the latent fit, made once
base_fit <- function(latent) {
  if (file.exists(base_file(latent))) return(readRDS(base_file(latent)))
  set.seed(string_seed(paste(fold_key, type, latent, "fit")))
  fit <- switch(latent,
                pql = fit_correction(train_k, t0 = t0, meshes = meshes,
                                     pql_args = list(max_refits = 10)),
                la = fit_correction_laplace(train_k, t0 = t0, meshes = meshes))
  fit$obj <- NULL
  saveRDS(fit, base_file(latent))
  fit
}

# the hyperparameters, as fit$hyper, at each point of theta
point_hyper <- function(fit, theta) {
  do.call(rbind, lapply(seq_len(nrow(theta)), function(j) {
    par_list <- fit$par_list
    index <- 0
    for (name in correction_hyper_names) {
      size <- length(par_list[[name]])
      par_list[[name]] <- unname(theta[j, index + seq_len(size)])
      index <- index + size
    }
    as.data.frame(correction_hyper(par_list))
  }))
}

# weighted quantiles of each hyperparameter over the points
hyper_quantiles <- function(hyper_table, weight, probs = c(0.05, 0.5, 0.95)) {
  out <- list()
  for (name in setdiff(names(hyper_table), c("kappa_omega", "kappa_eta"))) {
    x <- hyper_table[[name]]
    o <- order(x)
    cw <- cumsum(weight[o]) / sum(weight)
    q <- x[o][pmin(length(x), findInterval(probs, cw) + 1)]
    out[paste0(name, "_q", probs * 100)] <- q
  }
  tibble::as_tibble(out)
}

for (method in methods) {
  set.seed(string_seed(paste(fold_key, type, method)))
  start <- Sys.time()
  hyper <- NULL
  latent <- sub("_.*$", "", method)
  theta_method <- sub("^[a-z]+_", "", method)
  if (method == "plugin_once") {
    fit <- fit_correction(train_k, t0 = t0, meshes = meshes)
  } else {
    fit <- base_fit(latent)
  }
  fit_seconds <- as.numeric(Sys.time() - start, units = "secs")

  if (method == "plugin_once" || theta_method == "plugin") {
    lambda <- predict_correction(fit, test_k, m_draws_train = m_draws_train,
                                 m_draws_new = m_draws_test,
                                 n_draws = n_draws)
    theta <- matrix(hyper_vector(fit), 1)
    weight <- 1
  } else {
    obj <- correction_marginal(fit, if (latent == "la") "beta_binomial" else
      "pql")
    # the points are saved as soon as they are made, and reused by a rerun
    points_file <- file.path(hyper_dir, "fits",
                             sprintf("%s__%s__%i__points.rds", method,
                                     fold_key, k))
    hyper <- if (file.exists(points_file)) readRDS(points_file) else switch(
      theta_method,
      gauss = hyper_points(fit, obj, "gaussian", n = 100),
      ccd = hyper_points(fit, obj, "ccd"),
      is = hyper_points(fit, obj, "is", n = 200))
    saveRDS(hyper, points_file)
    report("%s: %i points of theta in %.0f s", method, nrow(hyper$theta),
           hyper$seconds)
    lambda <- predict_correction_hyper(fit, hyper, test_k,
                                       m_draws_train = m_draws_train,
                                       m_draws_new = m_draws_test,
                                       n_draws = n_draws, obj = obj)
    theta <- hyper$theta
    weight <- hyper$weight
    rm(obj)
  }
  stopifnot(!anyNA(lambda))
  seconds <- as.numeric(Sys.time() - start, units = "secs")
  mcse <- if (theta_method == "ccd") rep(NA_real_, ncol(lambda)) else
    hyper_mcse(lambda, if (theta_method == "is") hyper$diagnostics$ess)
  predictive_sd <- apply(lambda, 2, sd)

  hyper_table <- point_hyper(fit, theta)
  summary <- bind_cols(
    tibble(method = method, experiment = experiment_name, fold = fold_name,
           insecticide_type = type, rho = rho, n_train = nrow(train_k),
           n_test = nrow(test_k), n_points = nrow(theta),
           n_points_used = length(unique(attr(lambda, "point") %||% 1)),
           fit_seconds = fit_seconds, seconds = seconds,
           hyper_seconds = if (is.null(hyper)) NA else hyper$seconds,
           is_pareto_k = if (theta_method == "is")
             hyper$diagnostics$pareto_k else NA,
           is_ess = if (theta_method == "is") hyper$diagnostics$ess else NA,
           is_n_failed = if (theta_method == "is")
             hyper$diagnostics$n_failed else NA,
           # the Monte Carlo SE of the predictive mean of omega + xi from
           # the integration over theta, relative to the predictive SD of
           # the logit (median and max over the held-out assays)
           mcse_ratio_median = median(mcse / predictive_sd),
           mcse_ratio_max = max(mcse / predictive_sd),
           peak_memory_gb = peak_memory_gb()),
    hyper_quantiles(hyper_table, weight),
    fit_summary(fit) |> select(pql_refits, pql_last_theta_change,
                               convergence, max_gradient)
  )
  attr(lambda, "point") <- NULL
  attr(lambda, "point_means") <- NULL
  attr(lambda, "point_share") <- NULL
  saveRDS(list(summary = summary, test = test, logit_draws = lambda,
               theta = theta, weight = weight, hyper_table = hyper_table,
               hyper = hyper[setdiff(names(hyper), "theta")]),
          part_file(method))
  report("%s %s %s: %.0f s; omega range %.0f [%.0f, %.0f] km, eta range %.0f km, phi %.2f",
         fold_key, type, method, seconds, summary$range_omega_q50,
         summary$range_omega_q5, summary$range_omega_q95,
         summary$range_eta_q50, summary$phi_q50)
}
