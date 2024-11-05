# prep admin unit/region rasters and shapefiles

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load raster extent
mask <- rast("data/clean/raster_mask.tif")

# load all countries and all regions in Africa (per UN)
africa_countries_un <- country_region_lookup()

# harmonise these with the GADM layers. Need resolution <= 3 to include e.g.
# British Indian Ocean Territory, which is in the UN country list
gadm <- geodata::world(resolution = 3, path = "data/raw") %>%
  st_as_sf() %>%
  mutate(
    country_name = case_when(
      NAME_0 == "Central African Republic" ~ "CAR",
      # different apostrophe!
      NAME_0 == "Côte d'Ivoire" ~ "Côte d’Ivoire",
      NAME_0 == "Democratic Republic of the Congo" ~ "DR Congo",
      NAME_0 == "São Tomé and Príncipe" ~ "Sao Tome & Principe",
      .default = NAME_0
    )
  )

# subset down to the GADM shapefiles in the UN Africa region and the MAP mask
gadm_in_mask <- gadm %>%
  right_join(
    africa_countries_un,
    by = join_by(country_name)
  ) %>%
  mutate(
    mask_val = terra::extract(mask,
                              .,
                              fun = "mean",
                              na.rm = TRUE)[, 2],
  ) %>%
  filter(
    !is.na(mask_val)
  ) %>%
  select(-mask_val)

# now rasterise, mask, and crop, this to make a raster without extraneous bits.
gadm_raster_country <- gadm_in_mask %>%
  terra::rasterize(mask, field = "country_name")
gadm_raster_region <- gadm_in_mask %>%
  terra::rasterize(mask, field = "region")

# there are some cells missing from the raster, where they do not overlap with
# the countries
nrow(extract(mask, cells(mask)))
nrow(extract(gadm_raster_country, cells(gadm_raster_country)))
nrow(extract(gadm_raster_region, cells(gadm_raster_region)))

# find these cells
missing_cell_coords <- gadm_raster_country %>%
  extract(cells(mask), xy = TRUE) %>%
  filter(is.na(country_name)) %>%
  select(-country_name)

# find the nearest country/region to each cell
nearest_link <- missing_cell_coords %>%
  vect(geom = c("x", "y"),
       crs = crs(gadm_in_mask)) %>%
  nearest(vect(gadm_in_mask))
nearest_country <- gadm_in_mask$country_name[nearest_link$to_id]
nearest_region <- gadm_in_mask$region[nearest_link$to_id]

# find the corresponding codes
gadm_country_levels <- levels(gadm_raster_country)[[1]]
nearest_country_id <- gadm_country_levels$ID[match(nearest_country, gadm_country_levels$country_name)]
gadm_region_levels <- levels(gadm_raster_region)[[1]]
nearest_region_id <- gadm_region_levels$ID[match(nearest_region, gadm_region_levels$region)]

# put the values back in
missing_cell_index <- terra::cellFromXY(mask,
                                        missing_cell_coords)
gadm_raster_country[missing_cell_index] <- nearest_country_id
gadm_raster_region[missing_cell_index] <- nearest_region_id

# remask, for whatever reason
gadm_raster_country <- mask(gadm_raster_country, mask)
gadm_raster_region <- mask(gadm_raster_region, mask)

# check they are all filled in now
nrow(extract(mask, cells(mask)))
nrow(extract(gadm_raster_country, cells(gadm_raster_country)))
nrow(extract(gadm_raster_region, cells(gadm_raster_region)))

writeRaster(gadm_raster_country,
            "data/clean/country_raster.tif",
            overwrite = TRUE)

writeRaster(gadm_raster_region,
            "data/clean/region_raster.tif",
            overwrite = TRUE)

# then polygonise it again, for plotting
gadm_polys <- terra::as.polygons(gadm_raster_country) %>%
  st_as_sf() %>%
  left_join(
    africa_countries_un,
    by = join_by(country_name)
  )

saveRDS(gadm_polys,
        file = "data/clean/gadm_polys.RDS")
