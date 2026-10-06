# The density d_half of the population transform (#37) at which the
# population covariate is most spread across the bioassays, so that the data
# carry the most information about population.
#
#   Rscript R/pop_d_half_spread.R
#
# At each modelled bioassay (subset_modelled_bioassays() of the cleaned data,
# as in R/fit_model.R) the encounter transform 1 - exp(-d log(2) / d_half)
# of the population density d at its pixel and year (pop_density_matrix(),
# R/model_covariates.R), and its sd over the assays, for d_half on a grid; also
# the sd of the trended column g_dom(t) s, which enters the model. Writes
# outputs/net_screen/pop_d_half_spread.csv and prints the maximum.

suppressMessages({
  library(dplyr)
  library(stringr)
  library(terra)
})
source("R/functions.R")
source("R/model_covariates.R")
source("R/bioassay_subset.R")
baseline_year <- 1995
assays <- subset_modelled_bioassays(
  readRDS("data/clean/all_gambiae_complex_data.RDS"),
  rast("data/clean/raster_mask.tif"), baseline_year = baseline_year)
index <- index_bioassays(assays)
d <- pop_density_matrix(index$unique_cells, baseline_year,
                        max(assays$year_start))
density <- d[cbind(index$df$cell_id, index$df$year_id)]
year <- index$df$year_start
design <- complete_selection_design(selection_design())
dir.create("outputs/net_screen", showWarnings = FALSE, recursive = TRUE)
g <- pmax((year - design$trend_years[1]) / diff(design$trend_years), 0)

grid <- sort(unique(c(round(exp(seq(log(5), log(1000), length.out = 60))),
                      50, 200, 270)))
spread <- tibble(d_half = grid) |>
  rowwise() |>
  mutate(s = list(-expm1(-density * log(2) / d_half)),
         sd = sd(s), mean = mean(s),
         share_above_0.9 = mean(s > 0.9),
         sd_trended = sd(g * s)) |>
  select(-s) |>
  ungroup()
write.csv(spread, "outputs/net_screen/pop_d_half_spread.csv",
          row.names = FALSE)
best <- spread[which.max(spread$sd), ]
cat(sprintf("assays %d, density median %.0f per km2 (IQR %.0f-%.0f)\n",
            length(density), median(density), quantile(density, 0.25),
            quantile(density, 0.75)))
cat(sprintf("sd maximal at d_half %g: sd %.3f; within 0.01 of it for d_half %g-%g\n",
            best$d_half, best$sd, min(spread$d_half[spread$sd > best$sd - 0.01]),
            max(spread$d_half[spread$sd > best$sd - 0.01])))
best_trended <- spread[which.max(spread$sd_trended), ]
cat(sprintf("trended column sd maximal at d_half %g: sd %.3f\n",
            best_trended$d_half, best_trended$sd_trended))
print(spread |> filter(d_half %in% c(5, 50, 100, 150, 200, best$d_half, 270,
                                     300, 400)), n = Inf)
