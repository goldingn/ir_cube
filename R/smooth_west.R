# Maps of the latent smooths over West Africa for any smooth fits (#47), as
# R/v5_region_a.R's smooth_west.png (R/smooth_west_helpers.R): the posterior
# means of the selection multiplier exp(u_s) and of the pyrethroid floor
# plogis(floor_intercept + u_f), which, unlike u_f, compares between chains,
# over all usable chains of each fit and, for the fits PER_CHAIN names, over
# each chain too; with diagnostic regions A, B and C
# (outputs/species_runs/regions/regions.csv) outlined. And the smooths at the
# LLIN-pyrethroid bioassays of A (and of its countries), B and C, per fit and
# chain set, as R/v5_region_a.R's v5_region_a_smooths.csv.
#
#   USE_CHAINS="V5=1,2" PER_CHAIN="v5r_bb_full" Rscript R/smooth_west.R <name> \
#     <label>=<fitted_model.RData> [<label>=<fitted_model.RData> ...]
#
# The fits' rows in the order given; 250 draws per usable chain (even_draws(),
# R/species_fit_helpers.R; USE_CHAINS as there). PER_CHAIN: labels, comma
# separated. Writes
#   figures/species_runs/west/smooth_west_<name>.png     the means
#   figures/species_runs/west/smooth_west_<name>_sd.png  the first fit, all
#                                                        usable chains: the
#                                                        posterior sds of u_s
#                                                        and of the logit
#                                                        floor, with the
#                                                        bioassay sites
#   outputs/species_runs/west/<name>_region_smooths.csv  the smooths at the
#                                                        bioassays
# and prints, per fit, how far its chains agree over the West Africa grid.
# Plain R with greta loaded (the draws are greta objects); one fit at a time.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) >= 2, all(grepl("=", arguments[-1])))
name <- arguments[1]
fits <- setNames(sub("^[^=]*=", "", arguments[-1]),
                 sub("=.*$", "", arguments[-1]))
per_chain_fits <- trimws(strsplit(Sys.getenv("PER_CHAIN"), ",")[[1]])
draws_per_chain <- 250

suppressMessages({
  library(greta)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(terra)
  library(tidyterra)
  library(sf)
  library(ggtext)
})
source("R/functions.R")
source("R/dynamical_predictions.R")
source("R/species_fit_helpers.R")
source("R/smooth_west_helpers.R")

output_dir <- "outputs/species_runs/west"
figure_dir <- "figures/species_runs/west"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
options(width = 200)

regions <- read.csv("outputs/species_runs/regions/regions.csv")
map_grid <- west_map_grid()

# chains `chains` in words: "chain 1", "chains 1 and 2", "chains 1, 2 and 3"
chain_text <- function(chains) {
  if (length(chains) == 1) return(sprintf("chain %d", chains))
  sprintf("chains %s and %d", paste(head(chains, -1), collapse = ", "),
          tail(chains, 1))
}

# a fit's description for the captions: likelihood, the smooths' range and
# their sds (posterior medians, per chain)
describe_fit <- function(label, fit, parameters, chain) {
  options <- fit$options
  range <- if (smooth_range_fixed(options$smooth)) {
    sprintf("range fixed at %s km",
            format(1000 * options$smooth[["range"]], big.mark = ","))
  } else {
    sprintf("range (posterior median) %s", paste(vapply(
      c("selection", "floor"), function(kind) {
        sprintf("%s %s km", kind, format(round(median(smooth_range_km(
          parameters$variables, options, kind))), big.mark = ","))
      }, ""), collapse = ", "))
  }
  sds <- paste(vapply(c("selection", "floor"), function(kind) {
    sd <- c(smooth_sd(parameters$variables, kind))
    sprintf("%s %s", kind, paste(sprintf("%.1f", tapply(sd, chain, median)),
                                 collapse = ", "))
  }, ""), collapse = "; ")
  sprintf("%s: %s likelihood, %s, %s, smooth sds (posterior median by chain) %s",
          label, gsub("_", "-", options$likelihood),
          chain_text(sort(unique(chain))), range, sds)
}

west <- list()
fit_chains <- list()
tables <- list()
descriptions <- character(0)
sites <- NULL
for (label in names(fits)) {
  fit <- load_fit(fits[[label]])
  stopifnot(smooth_on(fit$options), smooth_floor_on(fit$options))
  used <- even_draws(fit, draws_per_chain *
                       length(usable_chains(fit$draws, label)), label)
  parameters <- fit_parameter_draws(fit, used$index)
  sets <- draw_chain_sets(used$chain, per_chain = label %in% per_chain_fits)
  first <- label == names(fits)[1]
  west[[label]] <- west_smooth_summaries(parameters, fit$options, fit$classes,
                                         map_grid$xy, sets, sd = first)
  rows <- region_rows(fit$df, regions)
  tables[[label]] <- region_smooth_table(
    label, rows, bioassay_smooths(parameters, fit$options, rows), sets,
    pyrethroid_intercept(parameters, fit$options, fit$classes))
  descriptions[[label]] <- describe_fit(label, fit, parameters, used$chain)
  fit_chains[[label]] <- sort(unique(used$chain))
  if (first) {
    # the bioassay sites at the centres of their grid cells: all, for u_s,
    # and those of the classes with the floor smooth, for the floor
    site_xy <- function(df) {
      cells <- unique(terra::cellFromXY(map_grid$grid,
                                        cbind(df$longitude, df$latitude)))
      as_tibble(terra::xyFromCell(map_grid$grid, cells[!is.na(cells)]))
    }
    floor_classes <- smooth_class_weight(fit$options, "floor",
                                         fit$df$class_id) == 1
    sites <- list(selection = site_xy(fit$df),
                  floor = site_xy(fit$df[floor_classes, ]))
  }
  report("%s: %d draws (%d chain sets) mapped at %d cells", label,
         parameters$n_draws, length(sets), length(map_grid$cells))
  rm(fit, parameters, rows)
  gc(verbose = FALSE)
}


# the smooths at the LLIN-pyrethroid bioassays ---------------------------------------

table <- bind_rows(tables)
write.csv(table, file.path(output_dir, sprintf("%s_region_smooths.csv", name)),
          row.names = FALSE)
cat("\nthe smooths at the LLIN-pyrethroid bioassays, weighted by mosquitoes tested: posterior mean (sd)\n")
print(as.data.frame(
  table %>%
    filter(quantity %in% c("multiplier", "floor"),
           group %in% blocks) %>%
    mutate(value = sprintf("%.2f (%.2f)", mean, sd)) %>%
    select(fit, chains, group, quantity, value) %>%
    pivot_wider(names_from = c(quantity, group), values_from = value)),
  right = FALSE)


# how far the chains agree over the West Africa grid ---------------------------------

# per chain set, the quantiles over the grid of exp(u_s) and the floor, and
# their roughness: the mean absolute difference of u_s and of the logit floor
# between neighbouring cells (about 14 km) of the grid
roughness <- function(values) {
  layer <- rast(map_grid$grid)
  layer[] <- NA_real_
  layer[map_grid$cells] <- values
  m <- as.matrix(layer, wide = TRUE)
  mean(abs(c(diff(m), diff(t(m)))), na.rm = TRUE)
}
cat("\nover the West Africa grid: quantiles (5%, 50%, 95%) and roughness (mean |difference| between neighbouring cells, log or logit scale)\n")
print(as.data.frame(bind_rows(lapply(names(west), function(label) {
  bind_rows(lapply(names(west[[label]]), function(set) {
    v <- west[[label]][[set]]
    tibble(fit = label, chains = set,
           multiplier = paste(sprintf("%.2f", exp(quantile(v$u_s, c(0.05, 0.5, 0.95)))),
                              collapse = " "),
           floor = paste(sprintf("%.2f", quantile(v$floor, c(0.05, 0.5, 0.95))),
                         collapse = " "),
           rough_u_s = roughness(v$u_s),
           rough_logit_floor = roughness(qlogis(v$floor)))
  }))
}))), digits = 3, right = FALSE)
cat("\nbetween chains, over the West Africa grid: correlation and mean |difference| of the chain means of u_s and of the logit floor\n")
for (label in names(west)) {
  chains <- setdiff(names(west[[label]]), "all")
  if (length(chains) < 2) next
  pairs <- combn(chains, 2)
  print(as.data.frame(bind_rows(lapply(seq_len(ncol(pairs)), function(i) {
    a <- west[[label]][[pairs[1, i]]]
    b <- west[[label]][[pairs[2, i]]]
    tibble(fit = label, pair = paste(pairs[, i], collapse = " / "),
           cor_u_s = cor(a$u_s, b$u_s),
           diff_u_s = mean(abs(a$u_s - b$u_s)),
           cor_logit_floor = cor(qlogis(a$floor), qlogis(b$floor)),
           diff_logit_floor = mean(abs(qlogis(a$floor) - qlogis(b$floor))))
  }))), digits = 3)
}


# the maps ---------------------------------------------------------------------------

layers <- west_map_layers(regions)
maps <- list()
for (label in names(west)) {
  for (set in names(west[[label]])) {
    maps[[length(maps) + 1]] <- list(label = chain_set_title(label, set),
                                     values = west[[label]][[set]])
  }
}
fit_sets <- vapply(names(west), function(label) {
  sprintf("%s (%s)", label, if (label %in% per_chain_fits) {
    "all chains and each"
  } else {
    chain_text(fit_chains[[label]])
  })
}, "")
figure <- west_map_figure(
  maps, map_grid, layers,
  title = sprintf("The latent smooths over West Africa: posterior means, %s",
                  paste(fit_sets, collapse = " and ")),
  caption = paste0(
    sprintf("Posterior means over the draws of each chain set (%d per chain). ",
            draws_per_chain),
    "Left: u_s as the multiplier exp(u_s) of the cumulative log fitness, clipped at 0.14 and 7.4 (blue: stronger selection);\n",
    "right: the pyrethroid floor plogis(floor_intercept + u_f), which, unlike u_f, compares between chains. ",
    "Regions A, B and C of outputs/species_runs/regions/regions.csv outlined.\n",
    paste(descriptions, collapse = "\n")))
ggsave(file.path(figure_dir, sprintf("smooth_west_%s.png", name)), figure,
       width = 12, height = 2.5 * length(maps), dpi = 110, bg = "white")

# the posterior sds, the first fit, all usable chains
label <- names(west)[1]
values <- west[[label]]$all
# one colour scale, to the 99th percentile of both (a few cells at the edge
# of the mask, far from the data, run higher)
sd_limit <- ceiling(10 * quantile(c(values$u_s_sd, values$logit_floor_sd),
                                  0.99)) / 10
sd_figure <- wrap_plots(
  west_map_panel(values$u_s_sd, "u_s_sd",
                 sprintf("%s: posterior sd of u_s (log multiplier)",
                         chain_set_title(label, "all")),
                 map_grid, layers, sites = sites$selection,
                 sd_limit = sd_limit),
  west_map_panel(values$logit_floor_sd, "logit_floor_sd",
                 sprintf("%s: posterior sd of the logit pyrethroid floor",
                         chain_set_title(label, "all")),
                 map_grid, layers, sites = sites$floor, sd_limit = sd_limit),
  ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = sprintf("The latent smooths over West Africa: posterior sds, %s, all usable chains",
                    label),
    caption = paste0(
      sprintf("Posterior sds over %d draws per chain, between chains as well as within them. ",
              draws_per_chain),
      "Left: u_s, with every bioassay site (the centre of its 14 km cell); right: the logit pyrethroid floor\n",
      "floor_intercept + u_f, with the sites of the classes with the floor smooth (pyrethroids and DDT). One colour scale, squished at ",
      sprintf("%.1f (the 99th percentile of both; max %.1f).\n", sd_limit,
              max(c(values$u_s_sd, values$logit_floor_sd))),
      descriptions[[label]]),
    theme = theme(plot.caption = element_text(hjust = 0, size = 8.5)))
ggsave(file.path(figure_dir, sprintf("smooth_west_%s_sd.png", name)),
       sd_figure, width = 12, height = 3.6, dpi = 110, bg = "white")
report("written %s and %s; peak memory %.1f GB",
       file.path(output_dir, sprintf("%s_region_smooths.csv", name)),
       file.path(figure_dir, sprintf("smooth_west_%s*.png", name)),
       peak_memory_gb())
