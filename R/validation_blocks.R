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
# Three requirements on the blocks, in this order:
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

# area of the convex hull of a set of projected points, in square km. Reported
# as a diagnostic; it is not what the cut is chosen on, for the reason below
hull_area <- function(east, north) {
  if (length(east) < 3) {
    return(max(diff(range(east)), diff(range(north))) * 5)
  }
  points <- st_multipoint(cbind(east, north))
  area <- as.numeric(st_area(st_convex_hull(points)))
  max(area, max(diff(range(east)), diff(range(north))) * 5)
}

# The quantity a cut is chosen on: how far each cell sits from the nearest cell
# in a different block.
#
# Maximising the smallest block's area does not work. Rotating a cut barely
# changes the areas — they are three thirds of the same country however it is
# sliced — so that objective is nearly flat across angles and ends up choosing
# on density noise. What it does not see is block *shape*, and the cuts it
# picked were thin slices, which put a large share of cells within a few km of
# the next block.
#
# Scoring the separation directly fixes that, because it is what the experiment
# is for. It also chooses the axis sensibly without being told to: cutting a
# long country across its length gives three roughly square blocks, while
# cutting along its length gives three thin ribbons, and the ribbons score far
# worse.
#
# The 25th percentile rather than the minimum, because one unlucky pair of cells
# either side of a boundary should not decide the cut, and weighted by records
# because that is what the scoring will see.
separation_score <- function(distance, outside_distance, records, block) {
  separation <- vapply(seq_along(block), function(i) {
    other <- block != block[i]
    if (!any(other)) return(NA_real_)
    min(distance[i, other], outside_distance[i])
  }, numeric(1))
  if (anyNA(separation)) return(-Inf)
  as.numeric(quantile(rep(separation, records), 0.25))
}

# Two families of contiguous partition, both cut at record-count terciles so
# the blocks are balanced by construction.
#
# Slabs: project onto a direction and cut across it. Three bands, each spanning
# the country's full width in the perpendicular direction.
#
# Sectors: take the bearing from the record-weighted centroid and cut on that.
# Three wedges meeting at the centre. For a compact country these are fatter
# than any slab, at the cost of cells near the centre being close to two other
# blocks, so which family wins is a real question and is left to the score.
tercile_cut <- function(ordering, records) {
  cumulative <- cumsum(records[ordering])
  total <- sum(records)
  block <- findInterval(cumulative - records[ordering] / 2,
                        total * seq_len(n_blocks - 1) / n_blocks) + 1
  out <- integer(length(ordering))
  out[ordering] <- block
  out
}

slab_cuts <- function(country_cells) {
  lapply(angles, function(angle) {
    radians <- angle * pi / 180
    projection <- country_cells$east * cos(radians) +
      country_cells$north * sin(radians)
    list(block = tercile_cut(order(projection), country_cells$records),
         family = "slab",
         parameter = angle)
  })
}

sector_cuts <- function(country_cells) {
  centre_east <- weighted.mean(country_cells$east, country_cells$records)
  centre_north <- weighted.mean(country_cells$north, country_cells$records)
  bearing <- atan2(country_cells$north - centre_north,
                   country_cells$east - centre_east)
  lapply(angles, function(angle) {
    rotated <- (bearing - angle * pi / 180) %% (2 * pi)
    list(block = tercile_cut(order(rotated), country_cells$records),
         family = "sector",
         parameter = angle)
  })
}

# how far each block is from taking its share of the records
imbalance_of <- function(records, block) {
  shares <- tapply(records, block, sum) / sum(records)
  if (length(shares) < n_blocks) return(Inf)
  max(abs(shares - 1 / n_blocks))
}

max_imbalance <- 0.12

assign_blocks <- function(country_cells) {

  coordinates <- as.matrix(country_cells[, c("longitude", "latitude")])
  distance <- fields::rdist.earth(coordinates, coordinates, miles = FALSE)
  # distance to the nearest cell that is in training whatever the cut, which
  # caps how much separation any cut near a border can achieve
  outside_distance <- apply(
    fields::rdist.earth(coordinates,
                        as.matrix(always_training[, c("longitude",
                                                      "latitude")]),
                        miles = FALSE),
    1, min)

  candidates <- c(slab_cuts(country_cells), sector_cuts(country_cells))
  imbalance <- vapply(candidates,
                      function(x) imbalance_of(country_cells$records, x$block),
                      numeric(1))

  allowed <- which(imbalance <= max_imbalance)
  if (length(allowed) == 0) {
    # no cut balances the records this well; take the best balanced one and let
    # the diagnostics table show it
    allowed <- which.min(imbalance)
  }

  score <- vapply(allowed,
                  function(i) separation_score(distance,
                                               outside_distance,
                                               country_cells$records,
                                               candidates[[i]]$block),
                  numeric(1))
  best <- candidates[[allowed[which.max(score)]]]

  list(block = best$block,
       family = best$family,
       angle = best$parameter,
       separation_q25 = max(score))

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

# cells that are in the training set of every fold, whatever the cuts: they
# belong to the separation objective, because a block on a national border is
# genuinely close to training data on the other side of it
always_training <- cells %>%
  filter(!country_name %in% splittable$country_name[splittable$split])

assignments <- lapply(
  splittable$country_name[splittable$split],
  function(this_country) {
    country_cells <- cells %>% filter(country_name == this_country)
    result <- assign_blocks(country_cells)
    country_cells %>%
      mutate(block = result$block,
             family = result$family,
             angle = result$angle,
             separation_q25 = result$separation_q25)
  }
)

cell_blocks <- bind_rows(assignments) %>%
  bind_rows(
    cells %>%
      filter(country_name %in% splittable$country_name[!splittable$split]) %>%
      mutate(block = NA_integer_, family = NA_character_,
             angle = NA_real_, separation_q25 = NA_real_)
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
    group_by(country_name, family, angle, separation_q25, block) %>%
    summarise(cells = n(), records = sum(records),
              area = hull_area(east, north),
              span_km = max(diff(range(east)), diff(range(north))),
              .groups = "drop") %>%
    group_by(country_name) %>%
    mutate(record_share = records / sum(records)) %>%
    ungroup()

  block_diagnostics <- per_country %>%
    group_by(country_name, family, angle, separation_q25) %>%
    summarise(cells = sum(cells), records = sum(records),
              worst_share = max(abs(record_share - 1 / n_blocks)),
              smallest_block_km2 = round(min(area)),
              smallest_block_span_km = round(min(span_km)),
              .groups = "drop") %>%
    arrange(separation_q25)

  write.csv(block_diagnostics, "outputs/cv_block_diagnostics.csv",
            row.names = FALSE)
  write.csv(cell_blocks %>% select(country_name, cell, longitude, latitude,
                                   records, block, family, angle),
            "outputs/cv_block_assignments.csv", row.names = FALSE)

  cat("\nper country, shortest within-country separation first",
      "(record share of a third would be 0.333):\n")
  print(as.data.frame(block_diagnostics %>%
    transmute(country_name, cells, records, family, angle,
              separation_q25 = round(separation_q25, 1),
              worst_share = round(worst_share, 3),
              smallest_block_km2)),
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
