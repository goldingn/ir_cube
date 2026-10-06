# The fraction of the An. gambiae complex that is An. arabiensis (#47), from
# the Vector Atlas species abundance maps (continental multispecies
# abundance, one band per species, 5 arcmin, static; figshare dataset
# 34044441, downloaded 6 October 2026). The complex members mapped there are
# arabiensis, coluzzii, gambiae, melas, merus and quadriannulatus; the
# fraction is arabiensis abundance over their sum, disaggregated to the
# model's 2.5 arcmin grid and masked to it.
#
#   Rscript R/prep_arabiensis_fraction.R
#
# Writes data/clean/arabiensis_fraction.tif.

source("R/packages.R")

abundance <- rast(
  "data/raw/va_species_maps_20261001/multispecies_abundance.tif")
complex_members <- c("arabiensis", "coluzzii", "gambiae", "melas", "merus",
                     "quadriannulatus")
stopifnot(all(complex_members %in% names(abundance)))

complex <- abundance[[complex_members]]
fraction <- complex[["arabiensis"]] / sum(complex)

mask <- rast("data/clean/raster_mask.tif")
fraction <- terra::disagg(fraction, fact = 2) |>
  terra::crop(mask) |>
  terra::resample(mask, method = "near") |>
  terra::mask(mask)
names(fraction) <- "arabiensis_fraction"

missing <- global(is.na(fraction) & !is.na(mask), "sum")[[1]]
cat(sprintf("mask cells with no fraction: %d of %d\n", missing,
            global(!is.na(mask), "sum")[[1]]))
writeRaster(fraction, "data/clean/arabiensis_fraction.tif", overwrite = TRUE)
