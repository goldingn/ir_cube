# Sub-national spatial block folds for cross-validation.
#
# Leave-one-country-out confounds spatial prediction with the country initial
# condition: a held-out country's `init_country_raw` has no data, reverts to the
# region prior, and that error is then amplified through 15-29 years of
# deterministic selection. Empirically the held-out bias across the six national
# folds correlates with the fitted country effect at r = -0.94. So those folds
# measure the difficulty of an entirely unsampled country, which is not a
# situation the deployed model faces.
#
# Holding out blocks *within* each country instead leaves two thirds of every
# country's records in training, so every country intercept stays identified and
# the test isolates spatial prediction.
#
# Only the six countries already selected for the national folds are split.
# Those were chosen by a criterion that matters here too — more than ten
# bioassays of every one of the nine insecticide types since 2010 — and keeping
# to them buys three things. The block experiment becomes directly comparable
# with the national one, differing in exactly the intended respect and no other:
# same countries, same insecticide coverage, and only whether the country's
# initial condition is identified. The held-out set keeps a controlled mix of
# insecticides, rather than becoming a weighted average over whatever the thin
# countries happen to hold. And each fold trains on 88% of the data rather than
# 67%, which is much closer to the production fit, so the result transfers to
# the deployed model more directly. Splitting all 34 data-bearing countries
# would cost the same three fits and buy precision that is not the binding
# constraint: the differences are already determined to a standard error of
# 0.003 to 0.006 in excess mean squared error.
#
# Every other country stays wholly in the training set of every fold, as do the
# other two blocks of each split country, so every country intercept is
# identified in every fold. That is the whole purpose of blocking rather than
# holding out countries.
#
# Requirements on the blocks, in priority order: each is one contiguous region;
# the three are large and even in area; and each carries at least a floor number
# of bioassays. Area comes before record balance deliberately — see the note on
# the objective below. No buffer is applied: the blocks are large enough that a
# few cells near a boundary cannot carry the result, and a buffer would remove
# training records from exactly the countries whose intercepts this design
# exists to keep identified.
#
# Sourcing this file defines `spatial_blocks`, a list of three training and test
# pairs. Running it as a script also writes the diagnostics table and the maps.

source("R/validation_folds.R")

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(sf)
  library(ggplot2)
})

n_blocks <- 3
# the countries to split: the same six the national folds hold out, chosen in
# validation_folds.R as those with more than ten bioassays of every insecticide
# type since 2010
countries_to_block <- countries_to_validate
# a country still needs enough cells to make three regions that are meaningfully
# regions, checked rather than assumed
min_cells_per_country <- 12
min_records_per_country <- 40
# floors on what a single block may carry, so that favouring area over record
# balance cannot leave a fold with too little to score
min_records_per_block <- 150
min_record_share <- 0.12
min_cells_per_block <- 10
# area of one cell of the mask at the equator, in square km; the mask is a
# 0.0416667 degree grid, so roughly 4.6 km on a side
cell_area_km2 <- (0.04166665 * 111.32) ^ 2
# candidate cut directions, in degrees
angles <- seq(0, 175, by = 5)
# equal area projection for Africa, so block areas are comparable
equal_area <- "+proj=laea +lat_0=0 +lon_0=20 +datum=WGS84 +units=km"

# One row per data-bearing cell, with the number of records it carries. A
# handful of cells straddle a border and carry records attributed to two
# countries; each cell is assigned to whichever country holds most of its
# records, so that a cell belongs to exactly one block
cells <- df %>%
  count(country_name, cell, name = "records") %>%
  group_by(cell) %>%
  summarise(country_name = country_name[which.max(records)],
            records = sum(records),
            .groups = "drop") %>%
  bind_cols(as.data.frame(terra::xyFromCell(mask, .$cell))) %>%
  rename(longitude = x, latitude = y)

projected <- cells %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(equal_area)
cells$east <- st_coordinates(projected)[, 1]
cells$north <- st_coordinates(projected)[, 2]

# Blocks are chosen to be large and even in *area*, subject to a floor on the
# number of bioassays each one carries.
#
# Balancing the record counts instead does not give usable blocks. Bioassay
# effort is wildly uneven in space — Kenya holds most of its records in one
# cluster west of Lake Victoria — so cutting at the record terciles puts one
# boundary straight through the densest cluster, and the block that gets a third
# of the records occupies a tiny area. The result is short separations exactly
# where most of the data is. Cutting at the *area* terciles instead lets a dense
# cluster sit whole inside one block, which both keeps the blocks large and
# lengthens the separation for the two sparse blocks, at the cost of an uneven
# split of the records between folds. Every record is still held out exactly
# once, so nothing is lost overall; only the per-fold record counts differ.
#
# Area is the number of land cells of the mask that fall in the block, not the
# convex hull of its data-bearing cells: the cut is applied to every land cell
# of the country, so the three areas are directly comparable and sum to the
# country. The floor on records stops a block being so empty that its fold
# cannot be scored.
#
# Maximising the smallest block's area is what "large and even" means here,
# since the three areas sum to the country: a cut that makes one block small
# necessarily makes another large.

# Two families of contiguous partition, both parameterised by where along a
# one-dimensional ordering the two cuts fall.
#
# Slabs: project onto a direction and cut across it. Three bands, each spanning
# the country's full width in the perpendicular direction.
#
# Sectors: take the bearing from the country's area centroid and cut on that.
# Three wedges meeting at the centre. For a compact country these are fatter
# than any slab, at the cost of cells near the centre being close to two other
# blocks, so which family wins is left to the objective.
slab_ordering <- function(east, north, angle) {
  radians <- angle * pi / 180
  east * cos(radians) + north * sin(radians)
}

sector_ordering <- function(east, north, angle, centre_east, centre_north) {
  bearing <- atan2(north - centre_north, east - centre_east)
  (bearing - angle * pi / 180) %% (2 * pi)
}

# Both families reduce to the same problem: a one-dimensional ordering of the
# country's land cells and of its data-bearing cells, and two cut points on it.
# That makes the search cheap. Sort both orderings once; the three areas are
# then just the three index ranges of the sorted land values, and the three
# record counts come from two binary searches into the sorted data values
# against a cumulative sum, so no candidate has to touch a cell.
#
# Maximising the smallest area over the two cut points puts them at exactly the
# area thirds, so the search over cut positions only does any work when the
# record floor binds and the split has to be pulled away from even.
cut_fractions <- sort(unique(c(seq(0.10, 0.90, by = 0.01), 1 / 3, 2 / 3)))

assign_blocks <- function(country_cells, country_land) {

  centre_east <- mean(country_land$east)
  centre_north <- mean(country_land$north)

  orderings <- c(
    lapply(angles, function(angle) list(
      family = "slab",
      parameter = angle,
      land = slab_ordering(country_land$east, country_land$north, angle),
      data = slab_ordering(country_cells$east, country_cells$north, angle)
    )),
    lapply(angles, function(angle) list(
      family = "sector",
      parameter = angle,
      land = sector_ordering(country_land$east, country_land$north, angle,
                             centre_east, centre_north),
      data = sector_ordering(country_cells$east, country_cells$north, angle,
                             centre_east, centre_north)
    ))
  )

  minimum_records <- max(min_records_per_block,
                         round(min_record_share * sum(country_cells$records)))
  n_land <- nrow(country_land)

  best <- NULL
  for (ordering in orderings) {

    land_sorted <- sort(ordering$land)
    data_order <- order(ordering$data)
    data_sorted <- ordering$data[data_order]
    record_cumulative <- c(0, cumsum(country_cells$records[data_order]))
    n_data <- length(data_sorted)

    indices <- unique(pmax(1L, pmin(n_land - 1L,
                                    round(cut_fractions * n_land))))

    for (a in seq_along(indices)) {
      for (b in seq_along(indices)) {
        if (b <= a) next
        k1 <- indices[a]
        k2 <- indices[b]
        areas <- c(k1, k2 - k1, n_land - k2)
        if (min(areas) <= 0) next
        # nothing to gain unless this beats the incumbent on the smallest block
        if (!is.null(best) && min(areas) <= best$min_area) next

        thresholds <- land_sorted[c(k1, k2)]
        if (thresholds[1] >= thresholds[2]) next

        # where the two thresholds fall among the sorted data cells
        d1 <- findInterval(thresholds[1], data_sorted)
        d2 <- findInterval(thresholds[2], data_sorted)
        if (min(c(d1, d2 - d1, n_data - d2)) < min_cells_per_block) next
        records <- c(record_cumulative[d1 + 1],
                     record_cumulative[d2 + 1] - record_cumulative[d1 + 1],
                     record_cumulative[n_data + 1] - record_cumulative[d2 + 1])
        if (min(records) < minimum_records) next

        best <- list(min_area = min(areas),
                     thresholds = thresholds,
                     ordering = ordering,
                     areas = areas,
                     records = records)
      }
    }
  }

  if (is.null(best)) {
    stop("no cut of ", country_cells$country_name[1],
         " meets the record floor; lower min_record_share")
  }

  list(block = findInterval(best$ordering$data, best$thresholds) + 1L,
       family = best$ordering$family,
       angle = best$ordering$parameter,
       smallest_area_km2 = min(best$areas) * cell_area_km2,
       area_evenness = min(best$areas) / max(best$areas),
       smallest_records = min(best$records))

}

splittable <- cells %>%
  group_by(country_name) %>%
  summarise(cells = n(), records = sum(records), .groups = "drop") %>%
  mutate(split = country_name %in% countries_to_block &
           cells >= min_cells_per_country &
           records >= min_records_per_country)

stopifnot(all(countries_to_block %in% splittable$country_name[splittable$split]))

cat(sprintf("\n%i countries with data, of which %i are split into blocks:\n",
            nrow(splittable), sum(splittable$split)))
print(as.data.frame(splittable %>% filter(split) %>% select(-split)))
cat(sprintf("the other %i countries (%i records) stay wholly in training\n",
            sum(!splittable$split),
            sum(splittable$records[!splittable$split])))

# Every land cell of the mask with the country it belongs to, so that block
# areas are true land areas. Built once by a rasterisation of the country
# polygons; see the note in the script that writes it
country_lookup_file <- "temporary/cell_country_lookup.RDS"
if (!file.exists(country_lookup_file)) {
  stop("missing ", country_lookup_file,
       "; build it with terra::rasterize of gadm_polys onto the mask")
}
land <- readRDS(country_lookup_file) %>%
  filter(country_name %in% splittable$country_name[splittable$split])
land_xy <- st_as_sf(
  data.frame(as.data.frame(terra::xyFromCell(mask, land$cell))),
  coords = c("x", "y"), crs = 4326
) %>%
  st_transform(equal_area) %>%
  st_coordinates()
land$east <- land_xy[, 1]
land$north <- land_xy[, 2]

assignments <- lapply(
  splittable$country_name[splittable$split],
  function(this_country) {
    country_cells <- cells %>% filter(country_name == this_country)
    result <- assign_blocks(country_cells,
                            land %>% filter(country_name == this_country))
    cat(sprintf("  %-14s %-7s angle %3i | smallest block %7.0f km2, ",
                this_country, result$family, result$angle,
                result$smallest_area_km2),
        sprintf("area evenness %.2f, smallest fold %i bioassays\n",
                result$area_evenness, result$smallest_records))
    country_cells %>%
      mutate(block = result$block,
             family = result$family,
             angle = result$angle,
             smallest_area_km2 = result$smallest_area_km2,
             area_evenness = result$area_evenness)
  }
)

cell_blocks <- bind_rows(assignments) %>%
  bind_rows(
    cells %>%
      filter(country_name %in% splittable$country_name[!splittable$split]) %>%
      mutate(block = NA_integer_, family = NA_character_, angle = NA_real_,
             smallest_area_km2 = NA_real_, area_evenness = NA_real_)
  )

stopifnot(
  nrow(cell_blocks) == nrow(cells),
  !anyDuplicated(cell_blocks$cell)
)


# the folds ----------------------------------------------------------------

block_of_cell <- setNames(cell_blocks$block, cell_blocks$cell)
df_blocked <- df %>%
  mutate(block = block_of_cell[as.character(cell)])

# Test sets are restricted to `test_min_year` and later, as the national folds
# are, so the two experiments are scored on the same era of data. Records at a
# held-out cell from before that year are dropped from the experiment rather
# than returned to training, which is also what the national folds do with
# pre-2010 records of a held-out country: putting them back would place the same
# pixel in both sets.
spatial_blocks <- lapply(seq_len(n_blocks), function(this_block) {
  held_out <- !is.na(df_blocked$block) & df_blocked$block == this_block
  list(
    training = df_blocked[!held_out, ] %>% select(-block),
    test = df_blocked[held_out & df_blocked$year_start >= test_min_year, ] %>%
      select(-block)
  )
})

cat("\nfold sizes:\n")
print(data.frame(
  fold = seq_len(n_blocks),
  training = vapply(spatial_blocks, function(x) nrow(x$training), numeric(1)),
  test = vapply(spatial_blocks, function(x) nrow(x$test), numeric(1)),
  test_cells = vapply(spatial_blocks, function(x) n_distinct(x$test$cell),
                      numeric(1)),
  test_countries = vapply(spatial_blocks,
                          function(x) n_distinct(x$test$country_name),
                          numeric(1)),
  test_insecticides = vapply(spatial_blocks,
                             function(x) n_distinct(x$test$insecticide_type),
                             numeric(1))
))

cat("\nheld-out bioassays per insecticide type, by fold:\n")
print(bind_rows(lapply(seq_len(n_blocks), function(i) {
  spatial_blocks[[i]]$test %>%
    count(insecticide_type) %>%
    mutate(fold = i)
})) %>%
  tidyr::pivot_wider(names_from = fold, values_from = n, names_prefix = "fold ") %>%
  as.data.frame())

# every blocked record is held out exactly once
stopifnot(
  sum(vapply(spatial_blocks, function(x) nrow(x$test), numeric(1))) ==
    sum(!is.na(df_blocked$block) & df_blocked$year_start >= test_min_year),
  all(vapply(spatial_blocks, function(x) {
    length(intersect(x$training$cell, x$test$cell)) == 0
  }, logical(1)))
)


# diagnostics and maps -----------------------------------------------------

if (identical(environment(), globalenv())) {

  # how far each held-out cell sits from the nearest training cell of the same
  # fold: the quantity the experiment exists to probe, and the one the
  # interpolation fold could only push to about 90 km
  separation <- bind_rows(lapply(seq_len(n_blocks), function(this_block) {
    test_cells <- cell_blocks %>% filter(!is.na(block), block == this_block)
    train_cells <- cell_blocks %>%
      filter(is.na(block) | block != this_block)
    distance <- fields::rdist.earth(
      as.matrix(test_cells[, c("longitude", "latitude")]),
      as.matrix(train_cells[, c("longitude", "latitude")]),
      miles = FALSE
    )
    test_cells %>%
      mutate(fold = this_block,
             distance_to_training = apply(distance, 1, min))
  }))

  cat("\nkm from each held-out pixel to the nearest training pixel:\n")
  print(as.data.frame(separation %>%
    group_by(fold) %>%
    summarise(pixels = n(),
              min = min(distance_to_training),
              q25 = quantile(distance_to_training, 0.25),
              median = median(distance_to_training),
              q75 = quantile(distance_to_training, 0.75),
              max = max(distance_to_training),
              .groups = "drop") %>%
    mutate(across(where(is.numeric), ~ round(.x, 1)))))

  per_country <- cell_blocks %>%
    filter(!is.na(block)) %>%
    group_by(country_name, family, angle, smallest_area_km2, area_evenness,
             block) %>%
    summarise(cells = n(), records = sum(records),
              span_km = max(diff(range(east)), diff(range(north))),
              .groups = "drop") %>%
    group_by(country_name) %>%
    mutate(record_share = records / sum(records)) %>%
    ungroup()

  block_diagnostics <- per_country %>%
    group_by(country_name, family, angle, smallest_area_km2, area_evenness) %>%
    summarise(cells = sum(cells), records = sum(records),
              worst_share = max(abs(record_share - 1 / n_blocks)),
              smallest_block_span_km = round(min(span_km)),
              .groups = "drop") %>%
    arrange(smallest_area_km2)

  write.csv(block_diagnostics, "outputs/cv_block_diagnostics.csv",
            row.names = FALSE)
  write.csv(cell_blocks %>% select(country_name, cell, longitude, latitude,
                                   records, block, family, angle,
                                   smallest_area_km2, area_evenness),
            "outputs/cv_block_assignments.csv", row.names = FALSE)

  cat("\nper country, smallest block area first:\n")
  print(as.data.frame(block_diagnostics %>%
    transmute(country_name, cells, records, family, angle,
              smallest_block_km2 = round(smallest_area_km2),
              area_evenness = round(area_evenness, 2),
              worst_record_share = round(worst_share, 3))),
    max = 400)

  # maps
  block_colours <- c("1" = "#1B7837", "2" = "#762A83", "3" = "#E08214")

  map_data <- cell_blocks %>%
    filter(!is.na(block)) %>%
    mutate(fold = factor(block))
  untested <- cell_blocks %>% filter(is.na(block))

  africa_map <- ggplot() +
    geom_sf(data = africa, fill = grey(0.97), colour = grey(0.8),
            linewidth = 0.2) +
    geom_sf(data = st_geometry(gadm_polys), fill = NA, colour = grey(0.85),
            linewidth = 0.15) +
    geom_point(data = untested,
               aes(x = longitude, y = latitude),
               colour = grey(0.6), size = 0.7, shape = 4) +
    geom_point(data = map_data,
               aes(x = longitude, y = latitude, colour = fold,
                   size = records),
               alpha = 0.8) +
    scale_colour_manual(values = block_colours, name = "held out in fold") +
    scale_size_area(max_size = 3, name = "bioassays") +
    coord_sf(expand = FALSE) +
    labs(
      x = "", y = "",
      title = "Sub-national spatial block folds",
      subtitle = paste("each country's data-bearing pixels cut into three",
                       "contiguous regions of about a third of its records;",
                       "\ncrosses are countries with too few pixels to split,",
                       "kept in training throughout")
    ) +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_blank())

  ggsave("figures/CV_spatial_blocks.png", africa_map, bg = "white",
         width = 9, height = 9)

  facet_map <- ggplot() +
    geom_sf(data = africa, fill = grey(0.97), colour = grey(0.8),
            linewidth = 0.2) +
    geom_point(data = map_data %>%
                 tidyr::crossing(panel = factor(seq_len(n_blocks))) %>%
                 mutate(role = ifelse(block == as.integer(panel),
                                      "held out", "training")),
               aes(x = longitude, y = latitude, colour = role),
               size = 0.7) +
    facet_wrap(~ panel, nrow = 1,
               labeller = labeller(panel = function(x) paste("fold", x))) +
    scale_colour_manual(values = c("held out" = "#B2182B",
                                   training = grey(0.7)), name = "") +
    coord_sf(expand = FALSE) +
    labs(x = "", y = "", title = "What each block fold holds out") +
    theme_minimal() +
    theme(legend.position = "bottom", panel.grid = element_blank(),
          axis.text = element_blank())

  ggsave("figures/CV_spatial_blocks_by_fold.png", facet_map, bg = "white",
         width = 12, height = 5)

  cat("\nmaps written to figures/CV_spatial_blocks.png and",
      "figures/CV_spatial_blocks_by_fold.png\n")

}
