# The latent smooths over West Africa and at the LLIN-pyrethroid bioassays of
# diagnostic regions A, B and C (outputs/species_runs/regions/regions.csv,
# R/diagnostic_regions.R; #47), for R/v5_region_a.R and R/smooth_west.R: the
# selection multiplier exp(u_s) and the pyrethroid floor
# plogis(floor_intercept + u_f), which, unlike u_f, compares between chains
# that split the floor differently between the two. Needs the packages and
# sources of those scripts (greta, dplyr, tibble, ggplot2, patchwork, terra,
# tidyterra, sf, ggtext; R/functions.R, R/dynamical_predictions.R,
# R/species_fit_helpers.R).

llin_pyrethroids <- c("Alpha-cypermethrin", "Deltamethrin", "Permethrin")
blocks <- c("A", "B", "C")
west_xlim <- c(-17.9, 24.5)
west_ylim <- c(4, 21)

# The grid of the West Africa maps: the cells of the water mask aggregated by
# 3 (about 14 km, as R/smooth_maps.R) in the maps' window, with the cells
# inside the mask and their coordinates
west_map_grid <- function() {
  grid <- terra::crop(
    terra::aggregate(terra::rast("data/clean/pfpr_water_mask.tif"), 3,
                     fun = "max", na.rm = TRUE),
    terra::ext(west_xlim[1] - 1, west_xlim[2] + 1, west_ylim[1] - 1,
               west_ylim[2] + 1))
  cells <- terra::cells(grid)
  list(grid = grid, cells = cells, xy = terra::xyFromCell(grid, cells))
}

# The chain sets of draws from chains `chain` (one per draw): all of them,
# then each chain on its own if `per_chain`, as a named list of draw indices
draw_chain_sets <- function(chain, per_chain = FALSE) {
  sets <- list(all = seq_along(chain))
  if (per_chain) {
    for (k in sort(unique(chain))) {
      sets[[paste0("chain ", k)]] <- which(chain == k)
    }
  }
  sets
}

# The fit's pyrethroid floor intercept, per draw (parameters from
# fit_parameter_draws())
pyrethroid_intercept <- function(parameters, options, classes) {
  parameters$floor_intercept[
    , smooth_intercept_index(options, match("Pyrethroids", classes))]
}

# Per chain set (draw_chain_sets()), the posterior means over the West
# Africa grid (west_map_grid()) of u_s and of the pyrethroid floor, with, if
# `sd`, the posterior sds of u_s and of the logit floor floor_intercept + u_f
# (the scale of the smooth)
west_smooth_summaries <- function(parameters, options, classes, xy, sets,
                                  sd = FALSE) {
  basis <- smooth_basis_at(options$smooth,
                           smooth_coords(xy[, 1], xy[, 2], options$smooth$crs))
  u_s <- parameters$smooth_weights$selection %*% t(basis)
  logit_floor <- pyrethroid_intercept(parameters, options, classes) +
    parameters$smooth_weights$floor %*% t(basis)
  rm(basis)
  floor <- plogis(logit_floor)
  column_sd <- function(m) {
    sqrt(colSums(sweep(m, 2, colMeans(m)) ^ 2) / (nrow(m) - 1))
  }
  lapply(sets, function(i) {
    out <- list(u_s = colMeans(u_s[i, , drop = FALSE]),
                floor = colMeans(floor[i, , drop = FALSE]))
    if (sd) {
      out$u_s_sd <- column_sd(u_s[i, , drop = FALSE])
      out$logit_floor_sd <- column_sd(logit_floor[i, , drop = FALSE])
    }
    out
  })
}

# The bioassays (rows of the fit's df) of regions A, B and C (`regions`, the
# regions table), with their region (block), whether an LLIN pyrethroid, and
# whether a pyrethroid
region_rows <- function(df, regions) {
  df %>%
    mutate(row = row_number()) %>%
    left_join(transmute(regions, country, block = region),
              by = c("country_name" = "country")) %>%
    filter(block %in% blocks) %>%
    mutate(llin = insecticide_type %in% llin_pyrethroids,
           pyrethroid = insecticide_class == "Pyrethroids")
}

# The smooths at bioassays `rows` (rows of the fit's df), draw by draw, as
# draws x bioassays matrices: us (u_s, 0 where the class has no selection
# smooth), uf (u_f, 0 where it has no floor smooth) and lf (the logit floor,
# floor_intercept + u_f), as R/v5_region_a.R's fit_terms()
bioassay_smooths <- function(parameters, options, rows) {
  cells <- unique(rows$cell)
  column <- match(rows$cell, cells)
  basis <- prediction_basis(options, cells)
  at_rows <- function(kind) {
    u <- parameters$smooth_weights[[kind]] %*% t(basis)
    u[, column, drop = FALSE] *
      rep(smooth_class_weight(options, kind, rows$class_id),
          each = nrow(u))
  }
  us <- at_rows("selection")
  uf <- at_rows("floor")
  intercept <- parameters$floor_intercept[
    , smooth_intercept_index(options, rows$class_id), drop = FALSE]
  list(us = us, uf = uf, lf = intercept + uf)
}

# per draw, the means of the columns `columns` of draws x bioassays matrix
# `m` weighted by `weight`
weighted_by_draw <- function(m, columns, weight) {
  c(m[, columns, drop = FALSE] %*% weight[columns]) / sum(weight[columns])
}
posterior <- function(x) {
  c(mean = mean(x), sd = sd(x), lower = unname(quantile(x, 0.025)),
    upper = unname(quantile(x, 0.975)))
}

# The smooths at the LLIN-pyrethroid bioassays of A (and of its countries),
# B and C of fit `label`, per chain set `sets` (draw_chain_sets()): per draw,
# the means weighted by mosquitoes tested of u_s, u_f, exp(u_s) and the
# pyrethroid floor f; their posterior mean, sd and 95% interval; and the
# spread between the bioassays of the posterior mean. From the bioassays
# `rows` (region_rows()), the smooths there `terms` (bioassay_smooths()) and
# the pyrethroid floor intercept per draw `intercept`
region_smooth_table <- function(label, rows, terms, sets, intercept) {
  a_llin <- which(rows$block == "A" & rows$llin)
  groups <- c(setNames(lapply(blocks, function(b) {
    which(rows$block == b & rows$llin)
  }), blocks),
  lapply(split(a_llin, rows$country_name[a_llin]), identity))
  names(groups)[-(1:3)] <- paste0("A: ", names(groups)[-(1:3)])
  weight <- rows$mosquito_number
  table <- list()
  for (chain_set in names(sets)) {
    draws <- sets[[chain_set]]
    quantities <- list(
      u_s = terms$us[draws, , drop = FALSE],
      u_f = terms$uf[draws, , drop = FALSE],
      multiplier = exp(terms$us[draws, , drop = FALSE]),
      floor = plogis(terms$lf[draws, , drop = FALSE]))
    for (group in names(groups)) {
      columns <- groups[[group]]
      for (quantity in names(quantities)) {
        by_draw <- weighted_by_draw(quantities[[quantity]], columns, weight)
        cell_means <- colMeans(quantities[[quantity]][, columns, drop = FALSE])
        w <- weight[columns] / sum(weight[columns])
        table[[length(table) + 1]] <- tibble(
          fit = label, chains = chain_set, group = group, quantity = quantity,
          bioassays = length(columns), cells = n_distinct(rows$cell[columns]),
          !!!as.list(posterior(by_draw)),
          spread = sqrt(sum(w * (cell_means - sum(w * cell_means)) ^ 2)))
      }
    }
    table[[length(table) + 1]] <- tibble(
      fit = label, chains = chain_set, group = "all",
      quantity = "pyrethroid floor intercept, plogis", bioassays = NA,
      cells = NA, !!!as.list(posterior(plogis(intercept[draws]))),
      spread = NA)
  }
  bind_rows(table)
}

# The layers under and over the maps: the country borders, and regions A, B
# and C outlined, with their labels' positions
west_map_layers <- function(regions) {
  sf_use_s2(FALSE)
  borders <- readRDS("data/clean/country_borders.RDS")
  outlines <- suppressMessages(
    borders %>%
      inner_join(transmute(regions, country_name = country, block = region),
                 by = "country_name") %>%
      filter(block %in% blocks) %>%
      st_make_valid() %>%
      group_by(block) %>%
      summarise(geometry = st_union(geometry), .groups = "drop"))
  missing <- setdiff(regions$country[regions$region %in% blocks],
                     borders$country_name)
  if (length(missing) > 0) report("no borders for %s", toString(missing))
  list(borders = borders, outlines = outlines,
       centres = suppressWarnings(st_point_on_surface(outlines)))
}

# A map over West Africa of `values` (one per cell of `map_grid`,
# west_map_grid()) of `kind`: "u_s" (as the multiplier exp(u_s), RdBu,
# clipped at 0.14 and 7.4), "floor" (the pyrethroid floor, magma, 0-60%), or
# "u_s_sd" and "logit_floor_sd" (posterior sds, YlGnBu, 0 to `sd_limit`),
# with `layers` (west_map_layers()) and, if given, the bioassay sites `sites`
# (a data frame of x and y) as points
west_map_panel <- function(values, kind, title, map_grid, layers,
                           sites = NULL, sd_limit = NULL) {
  layer <- rast(map_grid$grid)
  layer[] <- NA_real_
  layer[map_grid$cells] <- values
  scale <- switch(
    kind,
    u_s = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdBu"),
                               limits = c(-2, 2), oob = scales::squish,
                               na.value = "transparent",
                               labels = function(v) sprintf("%.2g", exp(v)),
                               name = "selection\nmultiplier\nexp(u_s)"),
    floor = scale_fill_viridis_c(option = "magma", limits = c(0, 0.6),
                                 oob = scales::squish, na.value = "transparent",
                                 labels = scales::percent,
                                 name = "pyrethroid\nfloor f"),
    scale_fill_gradientn(colours = RColorBrewer::brewer.pal(9, "YlGnBu"),
                         limits = c(0, sd_limit), oob = scales::squish,
                         na.value = "transparent",
                         name = "posterior sd\n(log or logit\nscale)"))
  p <- ggplot() +
    geom_sf(data = layers$borders, fill = grey(0.92), colour = NA) +
    geom_spatraster(data = layer) +
    geom_sf(data = layers$borders, fill = NA, colour = grey(0.5),
            linewidth = 0.15)
  if (!is.null(sites)) {
    p <- p + geom_point(aes(x, y), data = sites, size = 0.25, colour = "black",
                        alpha = 0.6)
  }
  p +
    geom_sf(data = layers$outlines, fill = NA, colour = "black",
            linewidth = 0.6) +
    geom_sf_label(aes(label = block), data = layers$centres, size = 3.5,
                  fontface = "bold", label.size = 0, alpha = 0.7) +
    scale +
    coord_sf(xlim = west_xlim, ylim = west_ylim, expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_ir_maps() +
    theme(plot.title = element_markdown(size = 10))
}

# The figure of the West Africa maps: one row per element of `maps`, each a
# list of label (its name in the panel titles) and values (one chain set's
# west_smooth_summaries()), with exp(u_s) left and the pyrethroid floor right
west_map_figure <- function(maps, map_grid, layers, title, caption) {
  panels <- list()
  for (map in maps) {
    panels[[length(panels) + 1]] <- west_map_panel(
      map$values$u_s, "u_s",
      sprintf("%s: selection multiplier exp(u_s)", map$label), map_grid,
      layers)
    panels[[length(panels) + 1]] <- west_map_panel(
      map$values$floor, "floor",
      sprintf("%s: pyrethroid floor plogis(intercept + u_f)", map$label),
      map_grid, layers)
  }
  wrap_plots(panels, ncol = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(title = title, caption = caption,
                    theme = theme(plot.caption = element_text(hjust = 0,
                                                              size = 8.5)))
}

# The panel title of chain set `set` of fit `label`
chain_set_title <- function(label, set) {
  sprintf("%s, %s", label, if (set == "all") "all usable chains" else set)
}
