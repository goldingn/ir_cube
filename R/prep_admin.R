# prep admin unit/region rasters and shapefiles

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load raster extent
mask <- rast("data/clean/raster_mask.tif")

# load all countries and all regions in Africa (per UN)
africa_countries_un <- country_region_lookup()

# harmonise these with the GADM layers. Need resolution <= 3 to include e.g.
# British Indian Ocean Terrotory, which is in the UN country list
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
