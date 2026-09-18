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
# Three requirements on the blocks, in this order:
#
#   contiguous            a block must be one connected region, not a scatter
#   about a third each    of that country's records, so the three folds are
#                         comparable and every record is held out exactly once
#   large spatial extent  a block must not collapse onto a small dense cluster,
#                         because the point of the experiment is to test
#                         prediction over distance
#
# Construction: within each country, project its data-bearing cell centroids
# onto a direction, sort along it, and cut at the record-count terciles. That
# gives three contiguous slabs, each spanning the country's full width in the
# perpendicular direction, and balanced in records by construction. The
# direction is chosen by searching over angles and taking the one that maximises
# the smallest block's area, which is what stops a slab collapsing onto a dense
# cluster.
#
# No buffer is applied. The blocks are large enough that a few cells near a
# boundary cannot carry the result, and a buffer would remove training records
# from the countries whose intercepts the design exists to keep identified.
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
# a country needs enough cells to make three regions that are meaningfully
# regions; below this it stays wholly in the training set for every fold
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

# area of the convex hull of a set of projected points, in square km. A
# degenerate set (one or two cells, or collinear cells) has zero area, so fall
# back on the extent along the longer axis times the cell size
hull_area <- function(east, north) {
  if (length(east) < 3) {
    return(max(diff(range(east)), diff(range(north))) * 5)
  }
  points <- st_multipoint(cbind(east, north))
  area <- as.numeric(st_area(st_convex_hull(points)))
  max(area, max(diff(range(east)), diff(range(north))) * 5)
}

# cut one country's cells into `n_blocks` contiguous slabs of roughly equal
# record count, along `angle`
slab_cut <- function(country_cells, angle) {
  radians <- angle * pi / 180
  projection <- country_cells$east * cos(radians) +
    country_cells$north * sin(radians)
  order_index <- order(projection)
  cumulative <- cumsum(country_cells$records[order_index])
  total <- sum(country_cells$records)
  # assign each cell to the tercile of the record count it falls in
  block <- findInterval(cumulative - country_cells$records[order_index] / 2,
                        total * seq_len(n_blocks - 1) / n_blocks) + 1
  out <- integer(nrow(country_cells))
  out[order_index] <- block
  out
}

# score a candidate cut: the area of its smallest block, penalised if any block
# takes too far from its share of the records
cut_score <- function(country_cells, block) {
  shares <- tapply(country_cells$records, block, sum) / sum(country_cells$records)
  if (length(shares) < n_blocks) return(-Inf)
  areas <- tapply(seq_along(block), block, function(index) {
    hull_area(country_cells$east[index], country_cells$north[index])
  })
  imbalance <- max(abs(shares - 1 / n_blocks))
  if (imbalance > 0.12) return(-Inf)
  min(unlist(areas))
}

assign_blocks <- function(country_cells) {
  scores <- vapply(angles,
                   function(angle) cut_score(country_cells,
                                             slab_cut(country_cells, angle)),
                   numeric(1))
  if (all(!is.finite(scores))) {
    # no angle balances the records; take the best balanced cut regardless of
    # area, and let the diagnostics table show it
    imbalance <- vapply(angles, function(angle) {
      shares <- tapply(country_cells$records,
                       slab_cut(country_cells, angle),
                       sum) / sum(country_cells$records)
      if (length(shares) < n_blocks) return(Inf)
      max(abs(shares - 1 / n_blocks))
    }, numeric(1))
    best <- angles[which.min(imbalance)]
  } else {
    best <- angles[which.max(scores)]
  }
  list(block = slab_cut(country_cells, best), angle = best)
}

splittable <- cells %>%
  group_by(country_name) %>%
  summarise(cells = n(), records = sum(records), .groups = "drop") %>%
  mutate(split = cells >= min_cells_per_country &
           records >= min_records_per_country)

cat(sprintf("\n%i countries with data: %i split into blocks, %i kept wholly in training\n",
            nrow(splittable), sum(splittable$split), sum(!splittable$split)))
cat("kept wholly in training (too few data-bearing cells):\n")
print(as.data.frame(splittable %>% filter(!split) %>% select(-split)))

assignments <- lapply(
  splittable$country_name[splittable$split],
  function(this_country) {
    country_cells <- cells %>% filter(country_name == this_country)
    result <- assign_blocks(country_cells)
    country_cells %>%
      mutate(block = result$block, angle = result$angle)
  }
)

cell_blocks <- bind_rows(assignments) %>%
  bind_rows(
    cells %>%
      filter(country_name %in% splittable$country_name[!splittable$split]) %>%
      mutate(block = NA_integer_, angle = NA_real_)
  )

stopifnot(
  nrow(cell_blocks) == nrow(cells),
  !anyDuplicated(cell_blocks$cell)
)


# the folds ----------------------------------------------------------------

block_of_cell <- setNames(cell_blocks$block, cell_blocks$cell)
df_blocked <- df %>%
  mutate(block = block_of_cell[as.character(cell)])

spatial_blocks <- lapply(seq_len(n_blocks), function(this_block) {
  list(
    training = df_blocked %>%
      filter(is.na(block) | block != this_block) %>%
      select(-block),
    test = df_blocked %>%
      filter(!is.na(block), block == this_block) %>%
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
                          numeric(1))
))

# every blocked record is held out exactly once
stopifnot(
  sum(vapply(spatial_blocks, function(x) nrow(x$test), numeric(1))) ==
    sum(!is.na(df_blocked$block)),
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
    group_by(country_name, angle, block) %>%
    summarise(cells = n(), records = sum(records),
              area = hull_area(east, north),
              span_km = max(diff(range(east)), diff(range(north))),
              .groups = "drop") %>%
    group_by(country_name) %>%
    mutate(record_share = records / sum(records)) %>%
    ungroup()

  block_diagnostics <- per_country %>%
    group_by(country_name, angle) %>%
    summarise(cells = sum(cells), records = sum(records),
              worst_share = max(abs(record_share - 1 / n_blocks)),
              smallest_block_km2 = round(min(area)),
              smallest_block_span_km = round(min(span_km)),
              .groups = "drop") %>%
    arrange(desc(worst_share))

  write.csv(block_diagnostics, "outputs/cv_block_diagnostics.csv",
            row.names = FALSE)
  write.csv(cell_blocks %>% select(country_name, cell, longitude, latitude,
                                   records, block, angle),
            "outputs/cv_block_assignments.csv", row.names = FALSE)

  cat("\nper country, worst record share imbalance first",
      "(a third would be 0.333):\n")
  print(as.data.frame(block_diagnostics %>%
    transmute(country_name, cells, records, angle,
              worst_share = round(worst_share, 3),
              smallest_block_km2, smallest_block_span_km)),
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
