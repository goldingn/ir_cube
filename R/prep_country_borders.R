# True country borders, for plotting, and an audit of the nearest-country fill.
#
# data/clean/gadm_polys.RDS is not a set of country borders. prep_admin.R
# rasterises GADM onto the 5 km mask, fills every mask cell that no country
# polygon covers with its nearest country, and then polygonises the raster
# back. The fill is deliberate — the model needs a country for every cell it
# predicts at — but the round trip leaves each filled offshore cell as a
# detached one-cell "part" of whichever country was nearest. Kenya comes out of
# it with 40 polygon parts, 36 of them more than 20 km from the mainland and
# every one an exact whole number of mask cells; a chain of them sits 270 to
# 290 km out in the Indian Ocean. Plotted, they read as bits of coastline
# annexed from the neighbours.
#
# So this writes the dissolved GADM geometry itself for plotting, and reports
# how far the filled cells sit from the country they were given, since
# predict.R:140 takes each prediction cell's country from the same raster.

source("R/packages.R")
source("R/functions.R")

mask <- rast("data/clean/raster_mask.tif")
africa_countries_un <- country_region_lookup()

# same harmonisation as prep_admin.R, reading the cached download
gadm <- geodata::world(resolution = 3, path = "data/raw") %>%
  st_as_sf() %>%
  mutate(
    country_name = case_when(
      NAME_0 == "Central African Republic" ~ "CAR",
      NAME_0 == "Côte d'Ivoire" ~ "Côte d’Ivoire",
      NAME_0 == "Democratic Republic of the Congo" ~ "DR Congo",
      NAME_0 == "São Tomé and Príncipe" ~ "Sao Tome & Principe",
      .default = NAME_0
    )
  ) %>%
  inner_join(africa_countries_un, by = join_by(country_name))

country_borders <- gadm %>%
  group_by(country_name, region) %>%
  summarise(.groups = "drop") %>%
  st_make_valid()

saveRDS(country_borders, "data/clean/country_borders.RDS")

cat("country borders written for", nrow(country_borders), "countries\n")
cat("polygon parts per country, worst five:\n")
print(head(sort(vapply(seq_len(nrow(country_borders)), function(i) {
  length(st_cast(st_cast(country_borders$geometry[i], "MULTIPOLYGON"), "POLYGON"))
}, numeric(1)) %>% setNames(country_borders$country_name), decreasing = TRUE), 5))


# how far are the filled cells from the country they were given? -----------

country_raster <- rast("data/clean/country_raster.tif")
land <- data.frame(cell = cells(mask)) %>%
  mutate(assigned = terra::extract(country_raster, cell)[[1]] %>% as.character()) %>%
  filter(!is.na(assigned))

coordinates <- terra::xyFromCell(mask, land$cell)
points <- st_as_sf(data.frame(coordinates), coords = c("x", "y"),
                   crs = st_crs(country_borders))

# distance from each cell to the border polygon it was assigned to, computed
# per country so the distance matrix stays small
land$km_from_assigned <- NA_real_
for (this_country in unique(land$assigned)) {
  index <- which(land$assigned == this_country)
  polygon <- country_borders %>% filter(country_name == this_country)
  if (nrow(polygon) == 0) next
  land$km_from_assigned[index] <-
    as.numeric(st_distance(points[index, ], polygon)) / 1000
}

cat("\nmask cells with a country:", nrow(land), "\n")
cat("cells not inside the country they were assigned to:",
    sum(land$km_from_assigned > 0, na.rm = TRUE),
    sprintf("(%.2f%%)\n", 100 * mean(land$km_from_assigned > 0, na.rm = TRUE)))
cat("of those, distance from the assigned country (km):\n")
print(round(quantile(land$km_from_assigned[land$km_from_assigned > 0],
                     c(0.5, 0.9, 0.99, 1), na.rm = TRUE), 1))
cat("cells more than 25 km outside their assigned country:",
    sum(land$km_from_assigned > 25, na.rm = TRUE), "\n")

far <- land %>%
  filter(km_from_assigned > 25) %>%
  count(assigned, name = "cells") %>%
  arrange(desc(cells))
if (nrow(far) > 0) {
  cat("\nby assigned country:\n")
  print(as.data.frame(far))
}

dir.create("outputs/review", showWarnings = FALSE, recursive = TRUE)
write.csv(land %>% filter(km_from_assigned > 5),
          "outputs/review/country_raster_fill_audit.csv", row.names = FALSE)
