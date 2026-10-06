# Screen of the conventional-net weight w (#37): the relative selection of a
# conventional ITN (cITN) to an LLIN in the pyrethroid-only net column,
#   nets_w(x, t) = use (w cITN + LLIN) / total nets,
# which is linear in w, so nets_w = A + w B from the w = 0.25 and 0.47 cubes
# (R/prep_net_use_pyrethroid.R). For each w on a grid, and each insecticide
# type, the net coefficient is refitted on all the fit's bioassays with
#
#   logit q_i = init_i + a_{country, type}
#               - sum_{s <= t_i} log(1 + o_s + exp(b) nets_w(x_i, s)) - t_i kappa
#
# init_i the fit's initial state of the bioassay's cell and type (recovered
# from its predicted logit at the posterior mean), a a free offset per
# country and type profiled out (so the initial state can readjust to w), o
# the fit's other columns' selection at their posterior mean exp(beta), kappa
# the fit's reversion. Deaths are binomial on effective counts, each assay
# divided by 1 + (n - 1) rho. Every w has one parameter per type, so their log
# likelihoods compare directly.
#
#   Rscript R/net_w_screen.R <fit RData> <label> <w 0.47 cube>
#
# Writes outputs/hinge_screen/<label>_net_w.csv (log likelihood and
# coefficient by w and type), <label>_net_w_periods.csv (the level residual
# by period, for each w) and <label>_net_w_logits.rds (each bioassay's fitted
# logit, for each w).

arguments <- commandArgs(trailingOnly = TRUE)
fit_file <- arguments[1]
label <- arguments[2]
cube_047 <- arguments[3]
w_grid <- c(0.25, 0.47, 0.75, 1)
n_cores <- min(8L, parallel::detectCores())

source("R/packages.R")
source("R/functions.R")
source("R/model_covariates.R")
source("R/dynamical_model.R")
source("R/dynamical_predictions.R")

f <- new.env()
load(fit_file, envir = f)
stopifnot(isFALSE(f$model_options$mortality_floor),
          f$model_options$selection_columns$net_w == 0.25)
df <- f$df
fold <- list(draws = f$draws, options = f$model_options,
             x_cells_init = f$x_cells_init)
n_total <- sum(vapply(fold$draws, nrow, integer(1)))
par <- dynamical_parameter_draws(fold, f$classes_index, f$types,
                                 draw_index = round(seq(1, n_total,
                                                        length.out = 200)))
effect <- apply(par$effect_type, c(2, 3), mean)
dimnames(effect) <- list(colnames(f$x_cell_years), f$types)
kappa <- if (is.null(par$kappa_type)) rep(0, length(f$types)) else
  colMeans(par$kappa_type)
rho <- colMeans(par$rho_types)
names(kappa) <- names(rho) <- f$types
# predicted logit at the posterior mean of each draw's logit
predicted <- colMeans(dynamical_logit(par, select(df, cell_id, country_id,
                                                  type_id, year_id),
                                      df, f$x_cell_years, f$cell_years_index))
rm(par)
invisible(gc())

n_cells <- max(f$cell_years_index$cell_id)
n_years <- max(f$cell_years_index$year_id)
x_row <- matrix(NA_integer_, n_cells, n_years)
x_row[cbind(f$cell_years_index$cell_id, f$cell_years_index$year_id)] <-
  seq_len(nrow(f$cell_years_index))
column_matrix <- function(values) matrix(values[x_row], n_cells, n_years)
nets_025 <- column_matrix(f$x_cell_years[, "nets"])
# the w = 0.47 column, read and padded as the design reads it
nets_047 <- read_padded_cube(cube_047, f$unique_cells, f$baseline_year,
                             f$baseline_year + n_years - 1)
nets_047 <- matrix(nets_047, n_cells, n_years)
b_cube <- (nets_047 - nets_025) / (0.47 - 0.25)
a_cube <- nets_025 - 0.25 * b_cube
cat(sprintf("B: range %.3f to %.3f, share negative %.4f\n",
            min(b_cube), max(b_cube), mean(b_cube < -1e-9)))
nets_w <- function(w) pmax(a_cube + w * b_cube, 0)

# cumulative selection from 1995 to each year, cells x years
cumulative_selection <- function(selection) {
  t(apply(log1p(selection), 1, cumsum))
}

types <- f$types
type_data <- function(type) {
  i <- which(df$insecticide_type == type)
  k <- match(type, types)
  others <- setdiff(colnames(f$x_cell_years), "nets")
  o <- column_matrix(as.vector(f$x_cell_years[, others, drop = FALSE] %*%
                                 effect[others, type]))
  fit_selection <- o + effect["nets", type] * nets_025
  cells_years <- cbind(df$cell_id[i], df$year_id[i])
  s_fit <- cumulative_selection(fit_selection)[cells_years]
  # the fit's initial state of each bioassay, from its predicted logit
  init <- predicted[i] + s_fit + df$year_id[i] * kappa[[type]]
  inflation <- 1 + (df$mosquito_number[i] - 1) * rho[[type]]
  list(type = type, rows = i, o = o, cells_years = cells_years,
       init = init, kappa = kappa[[type]], year_id = df$year_id[i],
       group = as.integer(factor(df$country_name[i])),
       y = df$died[i] / inflation, m = df$mosquito_number[i] / inflation,
       fit_slope = effect["nets", type])
}

binomial_loglik <- function(y, m, p) y * log(p) + (m - y) * log1p(-p)

# log likelihood of a type's bioassays given the logit before the country
# offsets, with the offsets profiled out by Newton's method (vectorised over
# groups); returns the log likelihood and the fitted logits
profile_offsets <- function(d, eta) {
  a <- numeric(max(d$group))
  for (iteration in 1:40) {
    p <- plogis(eta + a[d$group])
    gradient <- rowsum(d$y - d$m * p, d$group)[, 1]
    curvature <- rowsum(d$m * p * (1 - p), d$group)[, 1]
    step <- gradient / pmax(curvature, 1e-8)
    a <- pmin(pmax(a + step, -15), 15)
    if (max(abs(step)) < 1e-8) break
  }
  p <- pmin(pmax(plogis(eta + a[d$group]), 1e-12), 1 - 1e-12)
  list(loglik = sum(binomial_loglik(d$y, d$m, p)),
       logit = eta + a[d$group])
}
type_eta <- function(d, nets, b) {
  s <- cumulative_selection(d$o + exp(b) * nets)[d$cells_years]
  d$init - s - d$year_id * d$kappa
}
fit_w <- function(d, w) {
  nets <- nets_w(w)
  o <- optimize(function(b) -profile_offsets(d, type_eta(d, nets, b))$loglik,
                c(-12, 3))
  best <- profile_offsets(d, type_eta(d, nets, o$minimum))
  list(loglik = best$loglik, coefficient = exp(o$minimum),
       logit = best$logit)
}

data <- setNames(lapply(types, type_data), types)
results <- list()
period_rows <- list()
saved_logits <- list()
bioassay_period <- cut(df$year_start, c(1994, 2004, 2009, 2014, 2019, 2025),
              labels = c("-2004", "2005-09", "2010-14", "2015-19", "2020-"))
empirical <- qlogis((df$died + 0.5) / (df$mosquito_number + 1))
weight <- 1 / ((1 / (df$died + 0.5) +
                  1 / (df$mosquito_number - df$died + 0.5)) *
                 (1 + (df$mosquito_number - 1) * rho[df$insecticide_type]))
for (w in w_grid) {
  fits <- parallel::mclapply(data, fit_w, w = w, mc.cores = n_cores)
  results[[length(results) + 1]] <- tibble(
    w = w, type = types,
    class = df$insecticide_class[match(types, df$insecticide_type)],
    loglik = vapply(fits, `[[`, numeric(1), "loglik"),
    coefficient = vapply(fits, `[[`, numeric(1), "coefficient"))
  logit <- numeric(nrow(df))
  for (type in types) logit[data[[type]]$rows] <- fits[[type]]$logit
  saved_logits[[as.character(w)]] <- logit
  pyrethroid <- df$insecticide_class == "Pyrethroids"
  residual <- empirical - logit
  period_rows[[length(period_rows) + 1]] <- tibble(
    w = w, period = levels(bioassay_period),
    mean_observed = as.vector(tapply(df$died[pyrethroid] /
                                       df$mosquito_number[pyrethroid],
                                     bioassay_period[pyrethroid], mean)),
    mean_predicted = as.vector(tapply(plogis(logit[pyrethroid]),
                                      bioassay_period[pyrethroid], mean)),
    weighted_residual = vapply(split(which(pyrethroid), bioassay_period[pyrethroid]),
                               function(i) weighted.mean(residual[i],
                                                         weight[i]),
                               numeric(1)))
  cat(sprintf("w = %.2f done\n", w))
}
results <- bind_rows(results)
periods <- bind_rows(period_rows)
write.csv(results, file.path("outputs/hinge_screen",
                             sprintf("%s_net_w.csv", label)),
          row.names = FALSE)
# each bioassay's fitted logit under each w, for R/misfit_maps.R
saveRDS(saved_logits, file.path("outputs/hinge_screen",
                                sprintf("%s_net_w_logits.rds", label)))
write.csv(periods, file.path("outputs/hinge_screen",
                             sprintf("%s_net_w_periods.csv", label)),
          row.names = FALSE)
summary_table <- results |>
  group_by(type) |>
  mutate(delta = loglik - loglik[w == 0.25]) |>
  group_by(w) |>
  summarise(total = sum(delta), pyrethroids = sum(delta[class == "Pyrethroids"]),
            deltamethrin = delta[type == "Deltamethrin"],
            deltamethrin_coefficient = coefficient[type == "Deltamethrin"],
            ddt = delta[type == "DDT"], .groups = "drop")
options(width = 200)
print(as.data.frame(mutate(summary_table, across(where(is.numeric),
                                                 ~ round(.x, 2)))))
print(as.data.frame(periods |>
                      mutate(across(where(is.numeric), ~ round(.x, 3))) |>
                      select(w, period, mean_observed, mean_predicted,
                             weighted_residual)))
