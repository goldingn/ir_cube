# Experiments for #34: the two-stage correction with its hyperparameters
# estimated or integrated over in different ways, on one cross-validation fold
# and one insecticide type, scored on the fold's held-out assays.
#
#   Rscript R/two_stage_hyper_folds.R <experiment> <fold> <type index> <methods>
#
# e.g. Rscript R/two_stage_hyper_folds.R temporal_forecasting 2014 4 plugin,laplace,ccd
#
# Needs the fold's stage-one cache (outputs/two_stage/stage_one__<key>.rds,
# from R/run_two_stage_folds.R). <methods> is a comma-separated list of
#
#   plugin_once  the current algorithm: stage A, PQL, one re-estimation of
#                theta, plug-in theta_hat (fit_correction());
#   plugin       theta re-estimated until it stops moving (max_refits = 10),
#                plug-in (#34, point 1);
#   laplace      the plugin fit, integrated over draws from the Gaussian
#                approximation to p(theta | y) of the PQL marginal likelihood
#                (R/two_stage_hyper.R);
#   ccd          the same, over INLA's central composite design;
#   plugin_bb    the plug-in at the mode of the beta-binomial Laplace marginal
#                (#34, point 2), with its Laplace latent posterior;
#   laplace_bb,
#   ccd_bb       laplace and ccd on the beta-binomial Laplace marginal;
#   is           importance sampling of theta from the beta-binomial Laplace
#                marginal, PSIS-smoothed: the reference.
#
# All but plugin_once share the plugin fit, which is saved once to
# outputs/two_stage/hyper/fits/<key>__<k>.rds and reused. Each method saves
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
all_methods <- c("plugin_once", "plugin", "laplace", "ccd", "plugin_bb",
                 "laplace_bb", "ccd_bb", "is")
stopifnot(all(methods %in% all_methods))

source("R/two_stage_helpers.R")
source("R/two_stage_correction.R")
source("R/two_stage_hyper.R")

fold_key <- sprintf("%s__%s", experiment_name, fold_name)
hyper_dir <- "outputs/two_stage/hyper"
dir.create(file.path(hyper_dir, "fits"), showWarnings = FALSE,
           recursive = TRUE)
dir.create(file.path(hyper_dir, "parts"), showWarnings = FALSE)
base_file <- file.path(hyper_dir, "fits", sprintf("%s__%i.rds", fold_key, k))
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

# the plugin fit, made once
base_fit <- function() {
  if (file.exists(base_file)) return(readRDS(base_file))
  set.seed(string_seed(paste(fold_key, type, "plugin fit")))
  fit <- fit_correction(train_k, t0 = t0, meshes = meshes,
                        pql_args = list(max_refits = 10))
  fit$obj <- NULL
  saveRDS(fit, base_file)
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
  if (method == "plugin_once") {
    fit <- fit_correction(train_k, t0 = t0, meshes = meshes)
  } else {
    fit <- base_fit()
  }
  fit_seconds <- as.numeric(Sys.time() - start, units = "secs")

  if (method %in% c("plugin_once", "plugin")) {
    lambda <- predict_correction(fit, test_k, m_draws_train = m_draws_train,
                                 m_draws_new = m_draws_test,
                                 n_draws = n_draws)
    theta <- matrix(hyper_vector(fit), 1)
    weight <- 1
  } else {
    bb <- method %in% c("plugin_bb", "laplace_bb", "ccd_bb", "is")
    obj <- correction_marginal(fit, if (bb) "beta_binomial" else "pql")
    # the points are saved as soon as they are made, and reused by a rerun
    points_file <- file.path(hyper_dir, "fits",
                             sprintf("%s__%s__%i__points.rds", method,
                                     fold_key, k))
    hyper <- if (file.exists(points_file)) readRDS(points_file) else switch(
      sub("_bb$", "", method),
      plugin = hyper_points(fit, obj, "mode"),
      laplace = hyper_points(fit, obj, "laplace", n = 100),
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
  mcse <- if (method %in% c("ccd", "ccd_bb")) rep(NA_real_, ncol(lambda)) else
    hyper_mcse(lambda, if (method == "is") hyper$diagnostics$ess)
  predictive_sd <- apply(lambda, 2, sd)

  hyper_table <- point_hyper(fit, theta)
  summary <- bind_cols(
    tibble(method = method, experiment = experiment_name, fold = fold_name,
           insecticide_type = type, rho = rho, n_train = nrow(train_k),
           n_test = nrow(test_k), n_points = nrow(theta),
           n_points_used = length(unique(attr(lambda, "point") %||% 1)),
           fit_seconds = fit_seconds, seconds = seconds,
           hyper_seconds = if (is.null(hyper)) NA else hyper$seconds,
           is_pareto_k = if (method == "is")
             hyper$diagnostics$pareto_k else NA,
           is_ess = if (method == "is") hyper$diagnostics$ess else NA,
           is_n_failed = if (method == "is")
             hyper$diagnostics$n_failed else NA,
           # the Monte Carlo SE of the predictive mean of omega + xi from
           # the integration over theta, relative to the predictive SD of
           # the logit (median and max over the held-out assays)
           mcse_ratio_median = median(mcse / predictive_sd),
           mcse_ratio_max = max(mcse / predictive_sd),
           peak_memory_gb = peak_memory_gb()),
    hyper_quantiles(hyper_table, weight),
    fit_summary(fit) |> select(pql_refits, pql_last_theta_change,
                               convergence)
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
