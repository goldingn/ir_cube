# The full-data fit's maps at the pixel-year of every assay it was fitted to:
# the in-sample predictions of the main cross-validation figure (#36, panels A
# and D, fig_variance_explained.R).
#
#   Rscript R/full_fit_maps_at_assays.R
#
# Computed by the published map's own code, two_stage_weighted_draws()
# (R/two_stage_predictions.R), with identity weights at the cells with an
# assay and the setup of the map step of R/two_stage_maps.R: every year from
# 1995 to 2030, 1000 paired draws in batches of 100, the same seeds, and each
# cell's country from the country raster. So each value is the published map
# at that pixel-year:
#
#   dynamical  posterior mean of ilogit(m), the dynamical_mortality map
#   two_stage  posterior mean of ilogit(m + omega + xi), without u and p: the
#              two_stage_mortality map
#
# The map rasters hold the dynamical map only at the panel years (2000,
# 2005, ...), so the values are computed here rather than read off them.
#
# Needs outputs/two_stage/maps/<type>/fit.rds and dynamical.rds
# (R/two_stage_maps.R), built from the same data. Cells are taken chunk_size
# at a time, so memory is set by one type's fit and a chunk's draws. Writes
# outputs/full_fit_maps_at_assays.csv.

source("R/packages.R")
source("R/functions.R")
source("R/bioassay_subset.R")
source("R/two_stage_predictions.R")

# the modelled assays, as R/fit_model.R builds them
baseline_year <- 1995
invisible(list2env(modelled_bioassays(baseline_year, 2024), environment()))

# the map step's setup (R/two_stage_maps.R)
years_all <- baseline_year:2030
n_map_draws <- 1000
map_batch_size <- 100
chunk_size <- 500

country_raster <- rast("data/clean/country_raster.tif")

maps_for_type <- function(type) {
  keys <- df %>%
    filter(insecticide_type == type) %>%
    distinct(cell, year_start)
  cells <- unique(keys$cell)
  country <- as.character(terra::extract(country_raster,
                                         cells)$country_name)
  setup <- two_stage_setup(type, years_all, df, n_draws = n_map_draws,
                           batch_size = map_batch_size,
                           baseline_year = baseline_year)

  # posterior means, cells x years
  means <- list(dynamical = matrix(NA_real_, length(cells), length(years_all)),
                two_stage = matrix(NA_real_, length(cells), length(years_all)))
  for (rows in split(seq_along(cells),
                     ceiling(seq_along(cells) / chunk_size))) {
    draws <- two_stage_weighted_draws(setup, cells[rows], country[rows],
                                      diag(length(rows)))
    for (model in names(means)) {
      means[[model]][rows, ] <- colMeans(draws[[model]])
    }
  }

  report("%-18s %5i cells, %5i pixel-years; peak memory %.1f GB", type,
         length(cells), nrow(keys), peak_memory_gb())
  at <- cbind(match(keys$cell, cells), match(keys$year_start, years_all))
  tibble(cell = keys$cell, year_start = keys$year_start,
         insecticide_type = type,
         dynamical = means$dynamical[at],
         two_stage = means$two_stage[at])
}

maps <- bind_rows(lapply(types, maps_for_type))
write.csv(maps, "outputs/full_fit_maps_at_assays.csv", row.names = FALSE)
report("wrote outputs/full_fit_maps_at_assays.csv: %i pixel-years",
       nrow(maps))
