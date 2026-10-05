# Pre-screen of the net-use transform (#37), step 1: predictions and
# residuals of a saved full fit at every bioassay, under each mortality-floor
# mode, and the fit's covariates and selection effects, so that the GAMs
# (R/net_screen_residuals.R) and the pair analysis (R/net_screen_pairs.R) need
# not load the fit again.
#
#   Rscript R/net_screen_extract.R [fit file] [high-floor chains]
#
# Defaults: outputs/refit/fits/r2_main.RData (round 2 main fit, chain 4, in the
# high-floor mode, dropped from draws and kept in draws_all_chains) and chain 4.
# The low-floor mode is the fit's kept draws. About 1.5 GB.
#
# Writes outputs/net_screen/extract.rds: per assay (the fit's df), the
# posterior mean predicted mortality and rho under each mode, the empirical and
# predicted logits, Pearson residuals, and net use at the assay's cell-year
# (the model's nets column) with its running sum from 1995; the fit's
# x_cell_years and cell_years_index; and draws of the selection effects and
# reversion under the low-floor mode.

arguments <- commandArgs(trailingOnly = TRUE)
file <- if (length(arguments) >= 1) arguments[1] else
  "outputs/refit/fits/r2_main.RData"
high_chains <- if (length(arguments) >= 2) {
  as.integer(strsplit(arguments[2], ",")[[1]])
} else {
  4L
}
n_draws_keep <- 400
output_dir <- "outputs/net_screen"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  library(greta)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(terra)
})
source("R/functions.R")
source("R/dynamical_predictions.R")

f <- new.env()
load(file, envir = f)
df <- f$df
n_times <- max(f$cell_years_index$year_id)

# the draws of each mode: the kept draws are the low-floor mode; the
# high-floor chains are taken from draws_all_chains
mode_draws <- list(low = f$draws)
if (exists("draws_all_chains", envir = f, inherits = FALSE)) {
  mode_draws$high <- f$draws_all_chains[high_chains]
}

# net use (the model's nets column) at each cell-year and its running sum
x_row <- matrix(NA_integer_, max(f$cell_years_index$cell_id), n_times)
x_row[cbind(f$cell_years_index$cell_id, f$cell_years_index$year_id)] <-
  seq_len(nrow(f$cell_years_index))
nets <- matrix(f$x_cell_years[x_row, "nets"], nrow(x_row))
cum_nets <- t(apply(nets, 1, cumsum))

assay <- df |>
  transmute(row = row_number(), cell_id, cell, country_name, region,
            latitude, longitude, insecticide_class, insecticide_type,
            type_id, year_id, year = year_start, died,
            n = mosquito_number,
            nets = nets[cbind(cell_id, year_id)],
            cum_nets = cum_nets[cbind(cell_id, year_id)],
            irs = f$x_cell_years[x_row[cbind(cell_id, year_id)], "irs"])

modes <- list()
for (mode in names(mode_draws)) {
  fold <- list(draws = mode_draws[[mode]], options = f$model_options,
               x_cells_init = f$x_cells_init)
  n_total <- sum(vapply(fold$draws, nrow, integer(1)))
  index <- round(seq(1, n_total, length.out = n_draws_keep))
  par <- dynamical_parameter_draws(fold, f$classes_index, f$types,
                                   draw_index = index)
  logit <- dynamical_logit(par, select(df, cell_id, country_id, type_id,
                                       year_id),
                           df, f$x_cell_years, f$cell_years_index)
  p <- colMeans(plogis(logit))
  rho <- colMeans(par$rho_types)[df$type_id]
  n <- df$mosquito_number
  modes[[mode]] <- tibble(
    row = seq_len(nrow(df)),
    pred = p,
    pred_logit = colMeans(logit),
    rho = rho,
    # beta-binomial Pearson residual at the posterior mean
    pearson = (df$died - n * p) / sqrt(n * p * (1 - p) * (1 + (n - 1) * rho)),
    floor = mean(par$mortality_floor))
  cat(mode, ": chains", length(fold$draws), ", draws", par$n_draws,
      ", mean floor", signif(mean(par$mortality_floor), 3), "\n")
  if (mode == "low") {
    effects <- list(effect_type = par$effect_type, kappa_type = par$kappa_type,
                    rho_types = par$rho_types,
                    mortality_floor = par$mortality_floor)
  }
  rm(logit, par)
  invisible(gc())
}

saveRDS(list(file = file, assay = assay, modes = modes,
             x_cell_years = f$x_cell_years,
             cell_years_index = f$cell_years_index,
             unique_cells = f$unique_cells, types = f$types,
             classes = f$classes, classes_index = f$classes_index,
             baseline_year = f$baseline_year,
             options = f$model_options, effects = effects,
             high_chains = high_chains),
        file.path(output_dir, "extract.rds"))
cat("done\n")
