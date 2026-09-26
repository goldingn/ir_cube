# Fit the stage-A two-stage correction (#21, R/two_stage_correction.R) inside
# one cross-validation fold of #12, and save its held-out predictive draws in
# the same format as the #12 folds, so they are scored by the same code.
#
#   Rscript R/run_two_stage_folds.R <experiment> <fold> [m_ref=mean|loo]
#     [mesh=<tag>] [variants=omega_u,omega_xi_u] [stage=A|B] [terms=p|p,s]
#
# e.g. Rscript R/run_two_stage_folds.R spatial_interpolation all
#      Rscript R/run_two_stage_folds.R spatial_blocks 1
#      Rscript R/run_two_stage_folds.R spatial_blocks 1 loo
#      Rscript R/run_two_stage_folds.R temporal_forecasting 2014
#      Rscript R/run_two_stage_folds.R spatial_interpolation all mesh=base \
#        variants=omega_xi_u
#
# mesh picks one of the named mesh configurations in mesh_configs below. The
# default is "omega5000_xi2500", the reported configuration
# (doc/two_stage_plan.md, "Mesh resolution results"); mesh=base reproduces the
# earlier runs. Every tag but base is appended to the model name, e.g.
# two_stage_omega_xi_u_mesh-omega5000_xi2500, so its draws, residuals and
# fit-summary rows sit alongside the base ones, and the base files keep their
# unsuffixed names.
#
# stage=B fits stage B (PQL on the beta-binomial counts, R/two_stage_pql.R) on
# top of each type's stage-A fit, and saves it as two_stage_<variant>_pql
# (plus the same suffixes), e.g. two_stage_omega_xi_u_pql_mesh-omega5000_xi2500.
# Only the stage-B model is saved: the stage-A draws already exist. Stage A is
# refitted from the hyperparameters recorded in fit_summary.csv for the same
# fold, variant and mesh where available (warm start; cold if that fails).
#
# terms adds the optional error-structure terms (doc/two_stage_plan.md, "Error
# structure"): p, a static per-pixel effect, and s, a survey effect (survey =
# citation x country x year, survey_id()). They go into the model name after
# the variant, e.g. two_stage_omega_xi_u_p_s_pql_mesh-omega5000_xi2500. With s,
# three draws files are saved from the same latent draws, differing only in
# the survey term of the held-out draws:
#   two_stage_<variant>_p_s...        a fresh N(0, sigma_s^2) per held-out
#                                     survey (the predictive distribution of a
#                                     new assay; the primary one);
#   two_stage_<variant>_p_s_snone...  no survey term (the map prediction);
#   two_stage_<variant>_p_s_spost...  held-out assays in a training survey take
#                                     its posterior draw, others a fresh one.
#
# Run with OpenBLAS for CHOLMOD's supernodal factorisation (reference BLAS is
# ~10x slower), e.g.
#   LD_PRELOAD=.../libopenblas.so.0 OPENBLAS_NUM_THREADS=4 nice -n 10 Rscript ...
#
# Per fold (doc/two_stage_plan.md, steps 1-4):
#
#   1. rebuild the fold's training and test sets exactly as run_one_fold.R did,
#      and check the test set against the saved dynamical fold's test_df;
#   2. recompute the dynamical model's 2000 paired posterior draws of logit p at
#      the training and the test assays (R/dynamical_predictions.R), checking
#      the test draws against the saved p_draws;
#   3. for each insecticide type, fit both variants (omega_u, omega_xi_u) to the
#      training assays with the posterior mean logit (or the PSIS leave-out
#      mean) as the offset m_ref, and predict at that type's test assays with
#      the paired dynamical draws, so every correction draw sits on its own
#      dynamical draw (the cut posterior);
#   4. save the draws, a per-type fit summary, and the training residuals the
#      diagnostics need.
#
# Memory: a saved dynamical fold is 4-8 GB in memory, so it is loaded only once
# enough is free, and dropped as soon as the dynamical draws are extracted.

arguments <- commandArgs(trailingOnly = TRUE)
experiment_name <- arguments[1]
fold_name <- arguments[2]
# optional arguments are key=value; a bare third argument is m_ref, as before
options <- list(m_ref = "mean", mesh = "omega5000_xi2500",
                variants = "omega_u,omega_xi_u", stage = "A", terms = "")
for (argument in arguments[-(1:2)]) {
  if (!grepl("=", argument)) argument <- paste0("m_ref=", argument)
  key <- sub("=.*$", "", argument)
  if (!key %in% names(options)) stop("unknown argument: ", argument)
  options[[key]] <- sub("^[^=]*=", "", argument)
}
m_ref_type <- options$m_ref
mesh_tag <- options$mesh
variants <- strsplit(options$variants, ",")[[1]]
stage <- options$stage
terms <- strsplit(options$terms, ",")[[1]]
pixel_effect <- "p" %in% terms
survey_effect <- "s" %in% terms
stopifnot(all(terms %in% c("p", "s")))
# the error-structure terms, as a suffix of the variant: "", "_p", "_s", "_p_s"
terms_suffix <- paste0(if (pixel_effect) "_p" else "",
                       if (survey_effect) "_s" else "")
# the held-out survey draws saved (see above), and their model-name suffixes
survey_modes <- if (survey_effect) c("fresh", "none", "posterior") else "none"
# (only with s: without it the single "none" set is the model itself)
survey_suffix <- if (survey_effect) {
  c(fresh = "", none = "_snone", posterior = "_spost")
} else {
  c(none = "")
}
stopifnot(stage %in% c("A", "B"),
          !is.na(experiment_name), !is.na(fold_name),
          m_ref_type %in% c("mean", "loo"),
          length(variants) > 0,
          all(variants %in% c("omega_u", "omega_xi_u")))

report <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", sprintf(...), "\n")
  flush(stdout())
}

suppressMessages({
  sink("/dev/null")
  source("R/validation_folds.R")
  source("R/validation_covariates.R")
  sink()
})
source("R/dynamical_predictions.R")
source("R/two_stage_correction.R")
source("R/two_stage_pql.R")

draws_dir <- "outputs/cv_draws"
output_draws_dir <- "outputs/cv_draws_two_stage"
output_dir <- "outputs/two_stage"
dir.create(output_draws_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Mesh configurations (doc/two_stage_plan.md, "Mesh resolution results").
# omega and xi give the arguments to build_correction_mesh() for each field;
# xi = "omega" puts xi on the omega mesh. "base" is the configuration of all
# the runs before the mesh experiment
mesh_configs <- list(
  base = list(omega = list(), xi = list(max_nodes = 600)),
  # xi mesh of at most 1200 nodes: roughly halves its data-triangle edges
  xi1200 = list(omega = list(), xi = list(max_nodes = 1200)),
  # xi on the omega mesh
  xifull = list(omega = list(), xi = "omega"),
  # finer omega mesh (15 km cutoff, inner edges of at most 150 km, at most
  # 5000 nodes), with xi on the 1200-node mesh
  omega5000 = list(omega = list(cutoff = 15, max_edge_inner = 150,
                                max_nodes = 5000),
                   xi = list(max_nodes = 1200)),
  # the finer omega mesh, with xi on the base omega mesh (as in xifull)
  omega5000_xi2500 = list(omega = list(cutoff = 15, max_edge_inner = 150,
                                       max_nodes = 5000),
                          xi = list(max_nodes = 2500))
)
if (!mesh_tag %in% names(mesh_configs)) {
  stop("unknown mesh configuration: ", mesh_tag, "; one of ",
       paste(names(mesh_configs), collapse = ", "))
}
mesh_config <- mesh_configs[[mesh_tag]]

model_suffix <- paste0(if (stage == "B") "_pql" else "",
                       if (m_ref_type == "loo") "_loo" else "",
                       if (mesh_tag == "base") "" else
                         paste0("_mesh-", mesh_tag))

# xi(., t0) = 0 at the dynamical model's start year: the correction to the
# initial condition is omega, and xi accumulates selection anomalies after it
t0 <- baseline_year

# p rounds to exactly 0 or 1 in double precision once |logit p| exceeds ~37, so
# it is clamped before qlogis(). 1e-12 (logit +-27.6) rather than a looser clamp
# because it is the one R/stage_one_effective_parameters.R uses, so m_ref here
# and its leave-out means are on the same footing; it is also what dbetabinom()
# applies in the likelihood. The empirical logit of an assay is bounded by
# log(n + 0.5) / 0.5 ~ 6.5 at n = 300, so a clamped m is a residual of -20 or
# so either way: such assays carry little weight at v = 1 / 0.5 ~ 2 and are
# counted in the summary
clamp <- 1e-12
safe_logit <- function(p) qlogis(pmin(pmax(p, clamp), 1 - clamp))

# wait until enough memory is available before loading a saved fold; other
# jobs share the machine
wait_for_memory <- function(needed_gb = 12, poll_seconds = 180) {
  repeat {
    meminfo <- readLines("/proc/meminfo")
    available_kb <- as.numeric(gsub("\\D", "",
                                    grep("^MemAvailable", meminfo,
                                         value = TRUE)))
    available_gb <- available_kb / 1024 ^ 2
    if (available_gb >= needed_gb) return(invisible(available_gb))
    report("%.1f GB available, waiting for %.0f GB", available_gb, needed_gb)
    Sys.sleep(poll_seconds)
  }
}

# peak resident memory of this process so far, from the kernel
peak_memory_gb <- function() {
  status <- readLines("/proc/self/status")
  as.numeric(gsub("\\D", "", grep("^VmHWM", status, value = TRUE))) / 1024 ^ 2
}


# 1. the fold's training and test sets -------------------------------------

# found the same way run_one_fold.R found them when the fold was fitted
if (experiment_name == "spatial_interpolation") {
  training <- spatial_interpolation$training
  test <- spatial_interpolation$test
} else if (experiment_name == "spatial_blocks") {
  suppressMessages({
    sink("/dev/null")
    source("R/validation_blocks.R")
    sink()
  })
  index <- as.integer(fold_name)
  stopifnot(!is.na(index), index >= 1, index <= length(spatial_blocks))
  training <- spatial_blocks[[index]]$training
  test <- spatial_blocks[[index]]$test
} else if (experiment_name == "temporal_forecasting") {
  stopifnot(fold_name %in% names(temporal_forecasting_folds))
  training <- temporal_forecasting_folds[[fold_name]]$training
  test <- temporal_forecasting_folds[[fold_name]]$test
} else {
  stop("unknown experiment: ", experiment_name)
}

report("%s / %s, m_ref = %s, mesh = %s, variants = %s, terms = %s: %i training, %i held-out assays",
       experiment_name, fold_name, m_ref_type, mesh_tag,
       paste(variants, collapse = ","),
       if (length(terms)) paste(terms, collapse = ",") else "none",
       nrow(training), nrow(test))

# surveys (citation x country x year) of the training and held-out assays
stopifnot("citation" %in% names(training), "citation" %in% names(test))
survey_train <- survey_id(training$citation, training$country_name,
                          training$year_start,
                          project_km(training$longitude, training$latitude))
survey_test <- survey_id(test$citation, test$country_name, test$year_start,
                         project_km(test$longitude, test$latitude))


# 2. paired dynamical draws at the training and test assays ----------------

fold_file <- file.path(draws_dir, sprintf("dynamical__%s__%s.rds",
                                          experiment_name, fold_name))
stopifnot(file.exists(fold_file))
wait_for_memory()
report("loading %s", fold_file)
fold <- readRDS(fold_file)

# the rebuilt held-out set must be the one the fold was fitted with, row for
# row, or the draws cannot be paired
test_df <- fold$test_df
stopifnot(
  nrow(test) == nrow(test_df),
  identical(as.integer(test$cell_id), as.integer(test_df$cell_id)),
  identical(as.integer(test$type_id), as.integer(test_df$type_id)),
  identical(as.integer(test$year_id), as.integer(test_df$year_id)),
  identical(as.numeric(test$died), as.numeric(test_df$died)),
  identical(as.numeric(test$mosquito_number),
            as.numeric(test_df$mosquito_number))
)
experiment_label <- fold$experiment

# the stored test draws, thinned by the same rule the pairing uses
p_saved <- fold$p_draws
draw_index <- paired_draw_index(fold)
if (nrow(p_saved) > length(draw_index)) {
  p_saved <- p_saved[draw_index, , drop = FALSE]
}
fold$p_draws <- NULL
fold$p_draws_before <- NULL
fold$rho_class_draws <- NULL
invisible(gc())

# The rows are passed without their country_id, so that dynamical_predictions()
# looks the country up from the cell. That is what fit_fold() did: the state is
# per cell, and its initial condition is the country of the cell's first record
# in `df`. A few cells straddle a border and carry records attributed to two
# countries; passing the record's own country_id there gives a prediction the
# fitted model never made (in spatial_blocks 1, four held-out Uganda assays in
# cells whose first record is Kenyan differed from the saved draws by up to 2.4
# on the logit scale)
cell_country <- df %>%
  distinct(cell_id, .keep_all = TRUE) %>%
  select(cell_id, country_id)
model_country <- function(rows) {
  cell_country$country_id[match(rows$cell_id, cell_country$cell_id)]
}
border_train <- training$country_id != model_country(training)
report("%i training and %i held-out assays are in cells the model assigns to another country",
       sum(border_train), sum(test_df$country_id != model_country(test_df)))
time_dynamical <- system.time({
  p_test <- dynamical_predictions(fold, select(test_df, -country_id), df,
                                  x_cell_years, cell_years_index,
                                  classes_index, types)
  p_train <- dynamical_predictions(fold, select(training, -country_id), df,
                                   x_cell_years, cell_years_index,
                                   classes_index, types)
})
rm(fold)
invisible(gc())
n_dyn_draws <- nrow(p_train)
report("dynamical draws: %i draws at %i training and %i test assays in %.0f s",
       n_dyn_draws, ncol(p_train), ncol(p_test), time_dynamical[["elapsed"]])

# the recomputed test draws must reproduce the saved ones. Compared on the logit
# scale, where the model is additive, after the clamp (both sides can round p
# to 1 at very susceptible cell-years)
logit_test <- safe_logit(p_test)
max_logit_diff <- max(abs(logit_test - safe_logit(p_saved)))
report("recomputed vs saved test draws: max |logit diff| = %.3g",
       max_logit_diff)
# 1e-6: the recursion over long forecast windows accumulates rounding (4.9e-8
# on temporal_forecasting 2014), still far below any effect on the scores
stopifnot(max_logit_diff < 1e-6)
rm(p_test)

logit_train <- safe_logit(p_train)
n_clamped_train <- colSums(p_train > 1 - clamp | p_train < clamp)
rm(p_train)
invisible(gc())

# m_ref: the posterior mean logit, or its approximate leave-one-pixel-year-out
# counterpart from the grouped PSIS in stage_one_effective_parameters.R, which
# already falls back to the posterior mean where Pareto k > 0.7
m_mean <- colMeans(logit_train)
m_ref_all <- m_mean
loo_used <- FALSE
if (m_ref_type == "loo") {
  loo_file <- file.path(output_dir, sprintf("loo_m__%s__%s.rds",
                                            experiment_name, fold_name))
  if (file.exists(loo_file)) {
    loo <- readRDS(loo_file)
    loo <- loo[order(loo$training_row), ]
    stopifnot(
      nrow(loo) == nrow(training),
      identical(as.integer(loo$cell_id), as.integer(training$cell_id)),
      identical(as.integer(loo$type_id), as.integer(training$type_id)),
      identical(as.integer(loo$year_id), as.integer(training$year_id)),
      identical(as.numeric(loo$died), as.numeric(training$died)),
      # same draws and the same clamp, so the same posterior mean, except
      # possibly at border cells if the leave-out file was computed with the
      # records' own country_id (see above)
      max(abs(loo$m_ref - m_mean)[!border_train]) < 1e-8
    )
    # where the leave-out file's posterior mean is not the model's, its
    # leave-out mean is not either, so those assays fall back to the mean
    loo_mismatch <- abs(loo$m_ref - m_mean) > 1e-8
    m_ref_all <- ifelse(loo_mismatch, m_mean, loo$m_loo_pixel_year_fallback)
    loo_used <- TRUE
    report(paste("m_ref: leave-one-pixel-year-out means from %s (%i assays",
                 "fall back to the mean at k > 0.7, %i more at border cells)"),
           loo_file, sum(loo$k_pixel_year > 0.7 & !loo_mismatch),
           sum(loo_mismatch))
    rm(loo)
  } else {
    report("m_ref: %s not found, using the posterior mean", loo_file)
  }
}


# 3. stage A, per insecticide type -----------------------------------------

rho_table <- read.csv("outputs/bioassay_rho_hierarchical.csv")
rho_for_type <- setNames(rho_table$rho, rho_table$insecticide_type)
stopifnot(all(types %in% names(rho_for_type)))

# held-out draws start as the dynamical model's own, and each type's columns are
# replaced by the corrected draws when its fit succeeds; so a type whose fit
# fails keeps the dynamical draws and is flagged
# one set per variant and saved survey mode
output_keys <- as.vector(outer(variants, survey_modes, paste, sep = "|"))
p_out <- lapply(setNames(output_keys, output_keys), function(v) p_saved)

# median length of the mesh edges among triangles whose three nodes all carry
# data weight, i.e. the resolution the fields actually have where the data are
median_data_edge_km <- function(mesh, coords) {
  A <- mesh_basis(mesh, coords)
  data_nodes <- which(Matrix::colSums(A) > 0)
  tv <- mesh$graph$tv
  tv <- tv[rowSums(matrix(tv %in% data_nodes, ncol = 3)) == 3, ,
           drop = FALSE]
  if (nrow(tv) == 0) return(NA_real_)
  edges <- rbind(tv[, 1:2], tv[, 2:3], tv[, c(3, 1)])
  edges <- unique(t(apply(edges, 1, sort)))
  loc <- mesh$loc[, 1:2]
  median(sqrt(rowSums((loc[edges[, 1], ] - loc[edges[, 2], ]) ^ 2)))
}

# stage B warm-starts stage A at the hyperparameters recorded for this fold,
# variant, m_ref and mesh, with the recorded objective for comparison
recorded_stage_a_start <- function(type, variant) {
  summary_file <- file.path(output_dir, "fit_summary.csv")
  if (!file.exists(summary_file)) return(NULL)
  recorded <- read.csv(summary_file, stringsAsFactors = FALSE)
  if (!"mesh_config" %in% names(recorded)) return(NULL)
  row <- recorded %>%
    filter(experiment == experiment_name, as.character(fold) == fold_name,
           variant == !!variant, m_ref == m_ref_type,
           mesh_config == mesh_tag, insecticide_type == type,
           convergence == 0)
  if (nrow(row) != 1) return(NULL)
  start <- list(log_sigma_omega = log(row$sigma_omega),
                log_kappa_omega = log(sqrt(8) / row$range_omega),
                log_tau = log(row$tau))
  if (variant == "omega_xi_u") {
    start <- c(start, list(log_sigma_eta = log(row$sigma_eta),
                           log_kappa_eta = log(sqrt(8) / row$range_eta),
                           logit_phi = qlogis(row$phi)))
  }
  list(start = start, objective = row$objective)
}

summaries <- list()
residuals <- lapply(setNames(variants, variants), function(v) list())

for (k in seq_along(types)) {

  type <- types[k]
  rho <- rho_for_type[[type]]
  train_rows <- which(training$type_id == k)
  test_rows <- which(test_df$type_id == k)

  # a type with no training assays in this fold cannot be corrected; its
  # held-out assays (if any) keep the dynamical draws
  if (length(train_rows) == 0) {
    report("%s: no training assays, %i held-out assays keep the dynamical draws",
           type, length(test_rows))
    for (variant in variants) {
      summaries[[length(summaries) + 1]] <- tibble(
        experiment = experiment_name, fold = fold_name,
        model = paste0("two_stage_", variant, terms_suffix, model_suffix),
        variant = paste0(variant, terms_suffix, if (stage == "B") "_pql" else ""),
        m_ref = if (loo_used) "loo" else "mean",
        mesh_config = mesh_tag, insecticide_type = type, rho = rho, n_train = 0L,
        n_test = length(test_rows), error = "no training assays",
        fallback_dynamical = length(test_rows) > 0)
    }
    next
  }

  train_k <- tibble(
    lon = training$longitude[train_rows],
    lat = training$latitude[train_rows],
    year = training$year_start[train_rows],
    cell = training$cell[train_rows],
    died = training$died[train_rows],
    mosquito_number = training$mosquito_number[train_rows],
    m = m_ref_all[train_rows],
    rho = rho,
    survey = survey_train[train_rows]
  )
  stage_a <- empirical_logit(train_k$died, train_k$mosquito_number, rho)
  train_k$z <- stage_a$z
  train_k$v <- stage_a$v

  # T is this type's last training year: xi is represented up to it, and the
  # AR(1) forecast inside predict_correction() carries it to later test years
  # (all of them in the forecasting folds, a few types' in the others)
  T_k <- max(train_k$year)
  test_k <- tibble(
    lon = test_df$longitude[test_rows],
    lat = test_df$latitude[test_rows],
    year = test_df$year_start[test_rows],
    cell = test_df$cell[test_rows],
    m = colMeans(logit_test[, test_rows, drop = FALSE]),
    survey = survey_test[test_rows]
  )

  # the dynamical draws paired with the stored ones. With m_ref = loo the
  # draws are recentred on the leave-out mean: the mode shift in
  # predict_correction() is H^-1 A' D (m_ref - m_draw), so uncentred draws
  # would average the correction back to the in-sample posterior mean and undo
  # the leave-out offset; recentring keeps the offset at m_ref and the draws'
  # spread as the mechanistic uncertainty
  m_draws_train <- logit_train[, train_rows, drop = FALSE]
  if (loo_used) {
    m_draws_train <- sweep(m_draws_train, 2,
                           m_ref_all[train_rows] - m_mean[train_rows], "+")
  }
  m_draws_test <- logit_test[, test_rows, drop = FALSE]

  # one pair of meshes per type, shared by both variants so they differ only in
  # the xi term
  time_mesh <- system.time({
    coords <- coords_km(train_k)
    mesh <- suppressMessages(do.call(build_correction_mesh,
                                     c(list(coords), mesh_config$omega,
                                       verbose = FALSE)))
    mesh_xi <- if (identical(mesh_config$xi, "omega")) mesh else
      suppressMessages(do.call(build_correction_mesh,
                               c(list(coords), mesh_config$xi,
                                 verbose = FALSE)))
    edge_km <- median_data_edge_km(mesh, coords)
    edge_xi_km <- median_data_edge_km(mesh_xi, coords)
  })

  for (variant in variants) {

    set.seed(2026 + k)
    fit <- NULL
    error_message <- NA_character_
    summary_variant <- paste0(variant, terms_suffix,
                              if (stage == "B") "_pql" else "")
    fit_stage_a <- function(start = list()) {
      tryCatch(
        withCallingHandlers(
          fit_correction(train_k,
                         variant = variant,
                         t0 = t0,
                         T = T_k,
                         mesh = mesh,
                         mesh_xi = if (variant == "omega_xi_u") mesh_xi else
                           NULL,
                         start = start,
                         pixel_effect = pixel_effect,
                         survey_effect = survey_effect),
          # non-convergence is recorded from opt below, not as a warning
          warning = function(w) {
            if (grepl("nlminb did not converge", conditionMessage(w))) {
              invokeRestart("muffleWarning")
            }
          }
        ),
        error = function(e) {
          error_message <<- conditionMessage(e)
          NULL
        }
      )
    }
    stage_a_start <- if (stage == "B") recorded_stage_a_start(type, variant) else
      NULL
    time_fit <- system.time({
      fit <- fit_stage_a(if (is.null(stage_a_start)) list() else
        stage_a_start$start)
      warm_start <- !is.null(stage_a_start)
      if (warm_start && (is.null(fit) || fit$opt$convergence != 0)) {
        error_message <- NA_character_
        fit <- fit_stage_a()
        warm_start <- FALSE
      }
    })
    stage_a_objective <- if (is.null(fit)) NA_real_ else fit$opt$objective

    # stage B on top of a converged stage-A fit
    time_pql <- c(elapsed = 0)
    if (stage == "B" && !is.null(fit) && fit$opt$convergence == 0) {
      time_pql <- system.time(
        fit <- tryCatch(fit_correction_pql(fit, train_k),
                        error = function(e) {
                          error_message <<- paste("PQL:", conditionMessage(e))
                          NULL
                        })
      )
    }
    stage_b <- if (is.null(fit)) NULL else fit$stage_b
    get_b <- function(name, missing = NA) {
      if (is.null(stage_b) || is.null(stage_b[[name]])) missing else
        stage_b[[name]]
    }

    converged <- !is.null(fit) && fit$opt$convergence == 0
    max_gradient <- if (is.null(fit)) NA_real_ else
      max(abs(fit$obj$gr(fit$opt$par)))

    # predict whenever the fit produced an optimum: a non-converged nlminb
    # still has a usable mode if the gradient is small, but its draws are only
    # used if the fit converged
    time_predict <- c(elapsed = 0)
    use_fit <- converged && length(test_rows) > 0
    if (use_fit) {
      time_predict <- system.time({
        lambda <- tryCatch(
          predict_correction(fit, test_k,
                             m_draws_train = m_draws_train,
                             m_draws_new = m_draws_test,
                             n_draws = n_dyn_draws,
                             batch_size = 100,
                             survey = survey_modes),
          error = function(e) {
            error_message <<- paste("prediction:", conditionMessage(e))
            NULL
          }
        )
      })
      if (!is.null(lambda)) {
        if (!is.list(lambda)) lambda <- list(lambda)
        names(lambda) <- survey_modes
        for (mode in survey_modes) {
          stopifnot(!anyNA(lambda[[mode]]))
          p_out[[paste(variant, mode, sep = "|")]][, test_rows] <-
            plogis(lambda[[mode]])
        }
        rm(lambda)
      } else {
        use_fit <- FALSE
      }
    }

    hyper <- if (is.null(fit)) list() else fit$hyper
    get_hyper <- function(name) {
      if (is.null(hyper[[name]])) NA_real_ else hyper[[name]]
    }
    summary_row <- tibble(
      experiment = experiment_name,
      fold = fold_name,
      model = paste0("two_stage_", variant, terms_suffix, model_suffix),
      variant = summary_variant,
      m_ref = if (loo_used) "loo" else "mean",
      mesh_config = mesh_tag,
      insecticide_type = type,
      rho = rho,
      n_train = nrow(train_k),
      n_pixel_years = n_distinct(paste(train_k$cell, train_k$year)),
      n_test = length(test_rows),
      n_m_clamped = sum(n_clamped_train[train_rows] > 0),
      t0 = t0,
      T = T_k,
      max_test_horizon = if (length(test_rows) > 0)
        max(0, max(test_k$year) - T_k) else NA_real_,
      mesh_nodes = mesh$n,
      mesh_cutoff_km = attr(mesh, "cutoff"),
      mesh_edge_km = edge_km,
      mesh_xi_nodes = if (variant == "omega_xi_u") mesh_xi$n else NA,
      mesh_xi_edge_km = if (variant == "omega_xi_u") edge_xi_km else NA,
      sigma_omega = get_hyper("sigma_omega"),
      range_omega = get_hyper("range_omega"),
      sigma_eta = get_hyper("sigma_eta"),
      range_eta = get_hyper("range_eta"),
      phi = get_hyper("phi"),
      persistence = get_hyper("persistence"),
      tau = get_hyper("tau"),
      sigma_p = get_hyper("sigma_p"),
      sigma_s = get_hyper("sigma_s"),
      n_pixels = n_distinct(train_k$cell),
      n_surveys = n_distinct(train_k$survey),
      share_test_pixel_in_train = if (length(test_rows) > 0)
        mean(test_k$cell %in% train_k$cell) else NA_real_,
      share_test_survey_in_train = if (length(test_rows) > 0)
        mean(test_k$survey %in% train_k$survey) else NA_real_,
      objective = if (is.null(fit)) NA_real_ else fit$opt$objective,
      convergence = if (is.null(fit)) NA_integer_ else fit$opt$convergence,
      nlminb_message = if (is.null(fit)) NA_character_ else fit$opt$message,
      iterations = if (is.null(fit)) NA_integer_ else fit$opt$iterations,
      max_gradient = max_gradient,
      # estimates at the edge of what the meshes and the data can support,
      # flagged rather than acted on: a range under two median data-triangle
      # edges is barely resolved by the mesh, one over 5000 km is a
      # continent-wide constant, and phi near 0 makes eta white noise in time
      # (xi a random walk), against a prior centred on five-year persistence
      flag_range_omega = case_when(
        get_hyper("range_omega") < 2 * edge_km ~ "below 2 mesh edges",
        get_hyper("range_omega") > 5000 ~ "above 5000 km",
        .default = ""),
      flag_range_eta = case_when(
        variant != "omega_xi_u" ~ "",
        get_hyper("range_eta") < 2 * edge_xi_km ~ "below 2 xi-mesh edges",
        get_hyper("range_eta") > 5000 ~ "above 5000 km",
        .default = ""),
      flag_phi = case_when(
        variant != "omega_xi_u" ~ "",
        get_hyper("phi") < 0.2 ~ "phi < 0.2",
        get_hyper("phi") > 0.98 ~ "phi > 0.98",
        .default = ""),
      error = error_message,
      fallback_dynamical = length(test_rows) > 0 && !use_fit,
      time_mesh_s = time_mesh[["elapsed"]],
      time_fit_s = time_fit[["elapsed"]],
      time_optimise_s = if (is.null(fit)) NA_real_ else
        fit$timings[["optimise"]],
      time_predict_s = time_predict[["elapsed"]],
      stage = stage,
      stage_a_warm_start = if (stage == "B") warm_start else NA,
      stage_a_objective = stage_a_objective,
      stage_a_objective_recorded = if (is.null(stage_a_start)) NA_real_ else
        stage_a_start$objective,
      pql_passes = get_b("passes_first", NA_integer_),
      pql_converged = get_b("converged_first"),
      pql_damped = get_b("damped_first"),
      pql_rms_move = get_b("rms_move", NA_real_),
      pql_max_move = get_b("max_move", NA_real_),
      pql_refit = get_b("refit"),
      pql_refit_convergence = get_b("refit_convergence", NA_integer_),
      pql_passes_second = get_b("passes_second", NA_integer_),
      pql_converged_second = get_b("converged_second"),
      pql_rms_move_second = get_b("rms_move_second", NA_real_),
      pql_n_clamped = get_b("n_clamped", NA_integer_),
      time_pql_s = time_pql[["elapsed"]],
      time_pql_refit_s = if (is.null(stage_b)) NA_real_ else
        stage_b$timings[["refit"]],
      peak_memory_gb = peak_memory_gb()
    )
    summaries[[length(summaries) + 1]] <- summary_row

    # training residuals for the diagnostics: the empirical logit, its
    # variance, the offset, and the fitted latent mean omega + xi + u at the
    # mode, so z - m_ref - latent is the fitted residual
    if (!is.null(fit)) {
      latent <- as.vector(fit$A_latent %*% fit$mode)
      residuals[[variant]][[type]] <- tibble(
        training_row = train_rows,
        insecticide_type = type,
        cell = train_k$cell,
        lon = train_k$lon,
        lat = train_k$lat,
        year = train_k$year,
        region = training$region[train_rows],
        country_name = training$country_name[train_rows],
        died = train_k$died,
        mosquito_number = train_k$mosquito_number,
        z = train_k$z,
        v = train_k$v,
        m_ref = train_k$m,
        m_mean = m_mean[train_rows],
        latent = latent,
        residual = train_k$z - train_k$m - latent
      )
      # stage B: z and v stay the stage-A empirical logit, so the extreme
      # residual diagnostic in two_stage_metrics.R compares like with like;
      # the latent is stage B's, and stage A's is kept beside it
      if (!is.null(stage_b)) {
        residuals[[variant]][[type]]$latent_stage_a <-
          stage_b$lambda_a - train_k$m
        residuals[[variant]][[type]]$z_pql <- fit$tmb_data$z
        residuals[[variant]][[type]]$v_pql <- fit$tmb_data$v
      }
    }

    report(paste("%-18s %-10s n=%5i test=%4i T=%i nodes=%i/%s",
                 "conv=%s range_omega=%.0f sigma_omega=%.2f tau=%.2f%s%s",
                 "fit %.0fs predict %.0fs%s%s"),
           type, variant, nrow(train_k), length(test_rows), T_k, mesh$n,
           if (variant == "omega_xi_u") mesh_xi$n else "-",
           if (is.null(fit)) "ERROR" else fit$opt$convergence,
           get_hyper("range_omega"), get_hyper("sigma_omega"),
           get_hyper("tau"),
           paste0(if (pixel_effect) sprintf(" sigma_p=%.2f",
                                            get_hyper("sigma_p")) else "",
                  if (survey_effect) sprintf(" sigma_s=%.2f",
                                             get_hyper("sigma_s")) else ""),
           if (variant == "omega_xi_u")
             sprintf(" range_eta=%.0f sigma_eta=%.3f phi=%.2f",
                     get_hyper("range_eta"), get_hyper("sigma_eta"),
                     get_hyper("phi")) else "",
           time_fit[["elapsed"]], time_predict[["elapsed"]],
           if (!is.null(stage_b))
             sprintf(paste(" | PQL %i passes, RMS move %.3f, refit %s",
                           "(+%s passes), %.0fs"),
                     stage_b$passes_first, stage_b$rms_move, stage_b$refit,
                     stage_b$passes_second, time_pql[["elapsed"]]) else "",
           if (summary_row$fallback_dynamical)
             paste(" FALLBACK TO DYNAMICAL:",
                   if (is.na(error_message)) fit$opt$message else
                     error_message) else "")

    rm(fit)
    invisible(gc())
  }
}


# 4. save --------------------------------------------------------------------

summary_table <- bind_rows(summaries)

rho_type <- tibble(insecticide_type = types, rho = rho_for_type[types])
rho_test <- rho_for_type[types[test_df$type_id]]

for (output_key in output_keys) {
  variant <- sub("\\|.*$", "", output_key)
  survey_mode <- sub("^.*\\|", "", output_key)
  model <- paste0("two_stage_", variant, terms_suffix,
                  survey_suffix[[survey_mode]], model_suffix)
  # (stage is also a column of summary_table, so the name is built outside)
  saved_variant <- paste0(variant, terms_suffix,
                          if (stage == "B") "_pql" else "")
  variant_summary <- filter(summary_table, variant == !!saved_variant) %>%
    mutate(model = !!model)
  object <- list(
    model = model,
    experiment = experiment_label,
    fold = fold_name,
    test_df = test_df,
    # the survey term in the held-out draws (see the header)
    survey_draws = if (survey_effect) survey_mode else NA_character_,
    p_draws = p_out[[output_key]],
    # every model is scored at the per-type replicate-based rho, which is also
    # what inflated v in the fit; as draws x assays so the scoring code's
    # rho_draws path applies it column by column
    rho_draws = matrix(rho_test, nrow = n_dyn_draws, ncol = nrow(test_df),
                       byrow = TRUE),
    rho_type = rho_type,
    m_ref = if (loo_used) "loo" else "mean",
    fallback_types = variant_summary$insecticide_type[
      variant_summary$fallback_dynamical],
    fit_summary = variant_summary
  )
  file <- file.path(output_draws_dir,
                    sprintf("%s__%s__%s.rds", model, experiment_name,
                            fold_name))
  saveRDS(object, file)
  report("saved %s", file)
  # residuals do not depend on the survey mode: saved once
  if (survey_mode != survey_modes[1]) next

  residual_file <- file.path(output_dir,
                             sprintf("train_residuals__%s__%s__%s.rds",
                                     model, experiment_name, fold_name))
  saveRDS(bind_rows(residuals[[variant]]), residual_file)
  report("saved %s", residual_file)
}

# one row per experiment x fold x variant x type x m_ref x mesh; a rerun
# replaces its own rows rather than duplicating them. Rows written before the
# mesh_config column existed are the base meshes
summary_file <- file.path(output_dir, "fit_summary.csv")
if (file.exists(summary_file)) {
  previous <- read.csv(summary_file, stringsAsFactors = FALSE) %>%
    mutate(across(any_of(c("fold", "nlminb_message", "error",
                           "flag_range_omega", "flag_range_eta",
                           "flag_phi")),
                  ~ ifelse(is.na(.x), NA_character_, as.character(.x))))
  if (!"mesh_config" %in% names(previous)) {
    previous$mesh_config <- NA_character_
  }
  previous <- previous %>%
    mutate(mesh_config = ifelse(is.na(mesh_config), "base", mesh_config),
           .after = m_ref) %>%
    anti_join(summary_table %>% select(experiment, fold, variant, m_ref,
                                       mesh_config, insecticide_type),
              by = c("experiment", "fold", "variant", "m_ref", "mesh_config",
                     "insecticide_type"))
  summary_table <- bind_rows(previous, summary_table)
}
write.csv(summary_table, summary_file, row.names = FALSE)
report("appended to %s; peak memory %.1f GB", summary_file, peak_memory_gb())
