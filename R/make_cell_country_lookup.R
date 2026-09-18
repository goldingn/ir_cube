# One-off: which country each land cell of the mask belongs to. Cached, because
# the block construction needs true land areas rather than the convex hull of
# the data-bearing cells.
suppressMessages({library(terra); library(sf); library(dplyr)})
mask <- rast("data/clean/raster_mask.tif")
# the dissolved GADM geometry, not gadm_polys.RDS: that one carries the
# nearest-country fill of prep_admin.R, which would put ocean cells hundreds of
# km from a country into that country's area. See R/prep_country_borders.R
countries <- readRDS("data/clean/country_borders.RDS")
country_raster <- rasterize(vect(countries), mask, field = "country_name")
out <- data.frame(cell = which(!is.na(values(mask)))) %>%
  mutate(country_name = terra::extract(country_raster, cell)[[1]] %>%
           as.character()) %>%
  filter(!is.na(country_name))
saveRDS(out, "temporary/cell_country_lookup.RDS")
cat("land cells with a country:", nrow(out), "of", sum(!is.na(values(mask))), "\n")
print(head(sort(table(out$country_name), decreasing = TRUE), 8))
