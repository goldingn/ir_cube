# Is the reversion rate (#24) identifiable from these data? Simulate bioassay
# counts at the observed design under a known rate, refit with the rate
# estimated, and see whether it is recovered.
#
#   Rscript R/reversion_identifiability.R base  <rate> [threads] [warmup] [n_samples]
#   Rscript R/reversion_identifiability.R refit <rate> <seed> [threads] [warmup] [n_samples]
#   Rscript R/reversion_identifiability.R summarise
#
# <rate> is the per-year reversion rate on the logit scale (-kappa in
# reversion_kappa(), so >= 0): one value for all classes, or one per class,
# comma-separated in the order of `classes` below. 0 is no reversion.
#
# base: fit the default model (dynamical_model_options()) to the real data
# with reversion fixed at <rate>, and save the posterior means of its
# variables. These are the true values for the simulation at that rate, so
# that the simulated data resemble the real data whatever the rate: a larger
# rate is offset by stronger selection, as it would be in a fit.
#
# refit: simulate the bioassay counts at every observed record (same cell,
# year, type and number of mosquitoes) from the betabinomial with the base
# fit's mortality and rho, then fit the default model with the rate estimated
# (reversion = "estimated", half-normal(0, 0.1) per class), starting from the
# true values except the rate, which starts at 0.01 as in dynamical_inits().
# With REVERSION_SIM_START=cached, every variable starts where a fit to the
# real data would, from dynamical_inits() and temporary/inits.RDS, to check
# that the recovery does not depend on starting at the truth.
#
# summarise: per class, the true rate against its posterior, and the
# selection the fit attributes to the covariates against the truth, to
# outputs/reversion_identifiability.csv.
#
# Outputs of the fits are in outputs/reversion_sim/. Run under nice, with the
# greta 0.6 environment (doc/cv_run_plan.md, section 1).

arguments <- commandArgs(trailingOnly = TRUE)
mode <- arguments[1]
stopifnot(mode %in% c("base", "refit", "summarise"))
# REVERSION_SIM_DIR lets a smoke test write somewhere harmless
sim_dir <- Sys.getenv("REVERSION_SIM_DIR", "outputs/reversion_sim")
dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)

# per-year logit rate of the half-normal prior, as in dynamical_variables()
prior_sd <- 0.1

rate_label <- function(rate) {
  paste(formatC(rate, format = "f", digits = 3), collapse = "_")
}


# the observed design ---------------------------------------------------------

if (mode != "summarise") {
  rate <- as.numeric(strsplit(arguments[2], ",")[[1]])
  stopifnot(!anyNA(rate), all(rate >= 0))
  if (mode == "base") {
    seed <- NA
    extra <- arguments[-(1:2)]
  } else {
    seed <- as.integer(arguments[3])
    extra <- arguments[-(1:3)]
  }
  threads <- if (length(extra) >= 1) as.integer(extra[1]) else 3L
  warmup <- if (length(extra) >= 2) as.integer(extra[2]) else 1000L
  n_samples <- if (length(extra) >= 3) as.integer(extra[3]) else 1000L
  n_chains <- 4

  # threads before python, python before terra and sf (doc/cv_run_plan.md)
  source("R/greta_setup.R")
  start_greta(threads = threads)
}
source("R/packages.R")
source("R/functions.R")
source("R/bioassay_subset.R")
source("R/dynamical_model.R")

# the data and design as in fit_model.R
baseline_year <- 1995
final_data_year <- 2024
mask <- rast("data/clean/raster_mask.tif")
ir_africa <- readRDS(file = "data/clean/all_gambiae_complex_data.RDS")
df <- subset_modelled_bioassays(ir_africa,
                                mask,
                                insecticides_keep = modelled_insecticides,
                                baseline_year = baseline_year,
                                final_data_year = final_data_year)
classes <- unique(df$insecticide_class)
types <- unique(df$insecticide_type)
regions <- unique(df$region)
countries <- unique(df$country_name)
unique_cells <- unique(df$cell)
df <- df %>%
  mutate(
    cell_id = match(cell, unique_cells),
    region_id = match(region, regions),
    country_id = match(country_name, countries),
    class_id = match(insecticide_class, classes),
    type_id = match(insecticide_type, types)
  )
classes_index <- df %>%
  distinct(type_id, class_id) %>%
  arrange(type_id) %>%
  pull(class_id)
n_classes <- length(classes)

# non-centred, as when this was run: the simulation reads and writes the
# non-centred deviations, and temporary/inits.RDS has no logit_init_mean
default_options <- dynamical_model_options(init_centred = "none")
selection <- selection_design_matrix(unique_cells, baseline_year,
                                     final_data_year,
                                     default_options$selection_columns)
cell_years_index <- selection$cell_years_index
x_cell_years <- selection$x_cell_years
rm(selection)
x_cells_init <- init_covariate_matrix(unique_cells,
                                      default_options$selection_columns)
n_times <- max(cell_years_index$year_id)

# Plain-R versions of the model's quantities, for one set of variable values.
# The terms (dynamical_terms()) from the values:
terms_of <- function(values) {
  dynamical_terms(values, classes_index,
                  dynamical_lookups(df)$country_region_index, types,
                  default_options)
}
# the per-year reversion rate of each type
rate_type <- function(rate) {
  rep_len(rate, n_classes)[classes_index]
}
# covariates as cells x years x columns
x_cells <- aperm(array(x_cell_years,
                       c(n_times, nrow(x_cell_years) / n_times,
                         ncol(x_cell_years))),
                 c(2, 1, 3))
# the cumulative selection sum_{s <= t} log(1 + x_s' exp(beta_type)) at the
# cell, type and year of each record in `data`
record_cumulative_selection <- function(beta_type, data) {
  out <- numeric(nrow(data))
  for (type in unique(data$type_id)) {
    in_type <- which(data$type_id == type)
    cells <- sort(unique(data$cell_id[in_type]))
    x <- matrix(x_cells[cells, , , drop = FALSE], ncol = dim(x_cells)[3])
    log_w <- matrix(log1p(x %*% exp(beta_type[, type])), length(cells))
    cumulative <- t(apply(log_w, 1, cumsum))
    out[in_type] <- cumulative[cbind(match(data$cell_id[in_type], cells),
                                     data$year_id[in_type])]
  }
  out
}
# logit q_t without reversion at each record: the logit initial state
# (closed_form_states()) less the cumulative selection
record_logit_state <- function(terms, data) {
  lookups <- dynamical_lookups(df)
  country <- lookups$cell_country_lookup[data$cell_id]
  x_init <- select_init_covariates(x_cells_init, default_options)
  l <- logit_init_relative_rows(terms, country, data$type_id,
                                x_init[data$cell_id, , drop = FALSE])
  logit_init <- logit_init_from_relative(
    l, init_frac_constants(types)$min[data$type_id])
  logit_init - record_cumulative_selection(terms$beta_type, data)
}

# the model with reversion fixed at `rate` (0 for none) or estimated, with the
# likelihood over `data`
build <- function(data, reversion) {
  options <- default_options
  options$reversion <- reversion
  build_dynamical_model(train_df = data,
                        df = df,
                        x_cell_years = x_cell_years,
                        cell_years_index = cell_years_index,
                        classes_index = classes_index,
                        types = types,
                        options = options,
                        x_cells_init = x_cells_init)
}
fixed_reversion <- function(rate) {
  if (all(rate == 0)) FALSE else -rate
}

# the posterior means of the variables, in their greta dimensions
posterior_means <- function(built, draws) {
  posts <- do.call(calculate, c(built$variables,
                                list(values = draws, nsim = 1000)))
  lapply(posts, function(x) apply(x, 2:3, mean))
}

# worst Rhat of the variables
rhat_summary <- function(draws) {
  psrf <- coda::gelman.diag(draws, autoburnin = FALSE,
                            multivariate = FALSE)$psrf[, 1]
  psrf
}

sample_model <- function(built, inits) {
  mcmc(built$model,
       chains = n_chains,
       initial_values = replicate(n_chains, inits, simplify = FALSE),
       warmup = warmup,
       sampler = hmc(Lmin = 15, Lmax = 30),
       n_samples = n_samples)
}

report <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...), "\n")
  flush(stdout())
}


# base fit ------------------------------------------------------------------

if (mode == "base") {
  destination <- file.path(sim_dir, sprintf("base_%s.rds", rate_label(rate)))
  built <- build(df, fixed_reversion(rate))
  list2env(built$variables, environment())
  inits <- dynamical_inits(readRDS("temporary/inits.RDS"), built$variables,
                           columns = colnames(x_cell_years))
  report("base fit, rate %s | %d chains, %d warmup, %d samples, %d threads",
         toString(rate), n_chains, warmup, n_samples, threads)
  timing <- system.time(draws <- sample_model(built, inits))
  rhat <- rhat_summary(draws)
  report("done in %.1f h | worst Rhat %.3f, %d of %d above 1.1",
         timing[["elapsed"]] / 3600, max(rhat), sum(rhat > 1.1),
         length(rhat))
  saveRDS(list(rate = rate,
               means = posterior_means(built, draws),
               rhat = rhat,
               seconds = timing[["elapsed"]]),
          destination)
  cat("BASE COMPLETE\n")
}


# simulate and refit ----------------------------------------------------------

if (mode == "refit") {
  start <- Sys.getenv("REVERSION_SIM_START", "truth")
  stopifnot(start %in% c("truth", "cached"))
  destination <- file.path(sim_dir, sprintf("refit_%s_seed%d%s.rds",
                                            rate_label(rate), seed,
                                            if (start == "cached") "_cached"
                                            else ""))
  base <- readRDS(file.path(sim_dir, sprintf("base_%s.rds",
                                             rate_label(rate))))
  truth <- base$means

  # simulate: the mortality and rho at each record under the true values
  truth_terms <- terms_of(truth)
  p <- floored_mortality(
    plogis(record_logit_state(truth_terms, df) - df$year_id * -rate_type(rate)[df$type_id]),
    truth_terms$mortality_floor)
  rho <- truth_terms$rho_types[df$type_id]
  # betabinomial as betabinomial_p_rho(): a = p (1 / rho - 1), b = a (1 - p) / p
  set.seed(seed)
  a <- p * (1 / rho - 1)
  b <- (1 - p) * (1 / rho - 1)
  sim_df <- df %>%
    mutate(died = rbinom(n(), mosquito_number, rbeta(n(), a, b)))
  report("simulated %d records at rate %s, seed %d | mean mortality %.3f (observed %.3f)",
         nrow(sim_df), toString(rate), seed,
         sum(sim_df$died) / sum(sim_df$mosquito_number),
         sum(df$died) / sum(df$mosquito_number))

  built <- build(sim_df, "estimated")
  list2env(built$variables, environment())
  if (start == "truth") {
    inits <- truth[intersect(names(truth), names(built$variables))]
    inits$reversion_rate <- array(0.01, dim(built$variables$reversion_rate))
    inits <- do.call(greta::initials, inits)
  } else {
    inits <- dynamical_inits(readRDS("temporary/inits.RDS"), built$variables,
                             columns = colnames(x_cell_years))
  }
  report("refit with the rate estimated | %d chains, %d warmup, %d samples, %d threads",
         n_chains, warmup, n_samples, threads)
  timing <- system.time(draws <- sample_model(built, inits))
  rhat <- rhat_summary(draws)
  report("done in %.1f h | worst Rhat %.3f, %d of %d above 1.1 | reversion Rhat %s",
         timing[["elapsed"]] / 3600, max(rhat), sum(rhat > 1.1),
         length(rhat),
         toString(round(rhat[grep("reversion_rate", names(rhat))], 3)))

  # the posterior of the rate and of the selection coefficients, beta_type
  terms_draws <- calculate(reversion_rate = reversion_rate,
                           beta_type = built$terms$beta_type,
                           mortality_floor = built$terms$mortality_floor,
                           rho_types = built$terms$rho_types,
                           values = draws)
  saveRDS(list(rate = rate,
               seed = seed,
               start = start,
               truth = truth,
               sim_died = sim_df$died,
               draws = as.matrix(terms_draws),
               rhat = rhat,
               seconds = timing[["elapsed"]]),
          destination)
  cat("REFIT COMPLETE\n")
}


# summary ---------------------------------------------------------------------

if (mode == "summarise") {

  # The selection the covariates apply: for each record, the cumulative
  # log fitness sum_{s <= t} log(1 + x_s' exp(beta_type)) at its cell and type
  # up to its year, averaged over records by class. This is on the same
  # logit scale as the cumulative reversion t * rate, so a compensating bias
  # shows as both being too high together.
  mean_cumulative_selection <- function(beta_type) {
    c(tapply(record_cumulative_selection(beta_type, df), df$class_id, mean))
  }
  # the mean year index of the records by class, so t * rate averaged
  mean_t <- tapply(df$year_id, df$class_id, mean)

  prior_var <- prior_sd^2 * (1 - 2 / pi)
  files <- list.files(sim_dir, "^refit_.*\\.rds$", full.names = TRUE)
  rows <- lapply(files, function(file) {
    fit <- readRDS(file)
    true_rate <- rep_len(fit$rate, n_classes)
    true_selection <- mean_cumulative_selection(
      terms_of(fit$truth)$beta_type)
    d <- fit$draws
    rate_draws <- d[, grep("^reversion_rate", colnames(d)), drop = FALSE]
    n_covs <- ncol(x_cell_years)
    beta_columns <- grep("^beta_type", colnames(d))
    thin <- round(seq(1, nrow(d), length.out = 200))
    selection_draws <- t(vapply(thin, function(i) {
      mean_cumulative_selection(matrix(d[i, beta_columns], n_covs))
    }, numeric(n_classes)))
    rhat <- fit$rhat
    tibble(
      true_rate = paste(fit$rate, collapse = ","),
      seed = fit$seed,
      start = fit$start %||% "truth",
      insecticide_class = classes,
      class_true_rate = true_rate,
      rate_mean = colMeans(rate_draws),
      rate_lower = apply(rate_draws, 2, quantile, 0.025),
      rate_upper = apply(rate_draws, 2, quantile, 0.975),
      rate_sd = apply(rate_draws, 2, sd),
      # 1 - posterior variance / prior variance
      contraction = 1 - apply(rate_draws, 2, var) / prior_var,
      p_below_0.01 = colMeans(rate_draws < 0.01),
      p_above_0.03 = colMeans(rate_draws > 0.03),
      rate_rhat = rhat[grep("reversion_rate", names(rhat))],
      mean_t = c(mean_t),
      true_cumulative_selection = true_selection,
      cumulative_selection_mean = colMeans(selection_draws),
      cumulative_selection_lower = apply(selection_draws, 2, quantile, 0.025),
      cumulative_selection_upper = apply(selection_draws, 2, quantile, 0.975),
      worst_rhat = max(rhat),
      n_rhat_above_1.1 = sum(rhat > 1.1),
      hours = fit$seconds / 3600
    )
  })
  results <- bind_rows(rows)
  write_csv(results, "outputs/reversion_identifiability.csv")
  print(results, width = Inf)
}
