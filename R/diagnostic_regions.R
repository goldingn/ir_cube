# Diagnostic regions across Africa from current susceptibility to the LLIN
# pyrethroids (#47): groups of neighbouring countries whose current bioassay
# mortality is alike, in which to check that each model gets the broad trend,
# and above all the current level, right (R/region_diagnostics.R). The regions
# come from the data alone, not from any model's misfit, so they are fair to
# every model. The continental regions of analysis_region() are too coarse
# for this: one "West" hid Côte d'Ivoire, Burkina Faso and Ghana falling to
# about 16% mortality by 2019-2024 while Togo, Benin, Nigeria and Niger
# plateaued near 45% (R/west_figures.R).
#
#   Rscript R/diagnostic_regions.R [<number of regions>]
#
# The bioassays are the modelled LLIN-pyrethroid bioassays (alpha-cypermethrin,
# deltamethrin and permethrin, as the trend figure, R/species_compare.R), from
# outputs/species_runs/misfit/ref_f0_bioassays.csv (R/species_misfit.R; the
# data columns are the same in every fit's file). A site is a model cell
# within a country. The windows are current, 2019-2024 (the data end in 2024),
# and early, 2010-2015.
#
# 1. Country levels. A binomial GLMM of the current bioassays (lme4),
#      logit p = country + insecticide type + (1 | site) + (1 | bioassay),
#    gives the type effects (centred on the window's mosquitoes tested per
#    type) and the site and bioassay variances: bioassays cluster by site and
#    vary far beyond binomial between repeats, so a plain binomial standard
#    error is far too small. With those fixed, each site's log likelihood of a
#    level mu (its bioassays' binomial likelihoods integrated over the
#    bioassay and site effects, by convolution on a grid; site_curves()) is
#    computed once. The level m of any set of sites (a country, a merged group,
#    a resampled country) maximises the sum of its sites' curves plus a weakly
#    informative normal prior (sd prior_sd, centred on Africa's level), which
#    matters only where a few bioassays killed every mosquito (Eswatini); its
#    precision w is the curvature there. The countries' levels agree with the
#    GLMM's own estimates (printed).
# 2. Adjacency: countries within 1 km of each other in
#    data/clean/country_borders.RDS; islands linked to their nearest mainland
#    country; the countries with current bioassays must form one connected
#    graph.
# 3. Spatially constrained agglomerative clustering of the countries with
#    current bioassays (cluster()). Only adjacent groups merge, at the
#    precision-weighted Ward cost w_g w_h / (w_g + w_h) (m_g - m_h)^2, which is
#    about the chi-squared statistic (1 df) for the two groups having the same
#    level; a merged group's level is re-estimated from its pooled sites, not
#    averaged over its countries. First, while any group is below the minima
#    (min_current bioassays and min_sites sites in 2019-2024, min_early
#    bioassays in 2010-2015), the cheapest merge of such a group with a
#    neighbour; then the cheapest merge of any two adjacent groups, down to
#    one. The cut (the number of regions) is the given one, or else
#    choose_cut(): within k_range, the cut before the largest jump in merge
#    cost relative to every cost before it.
# 4. Countries without current bioassays (with bioassays in other years, or
#    with at least min_transmission cells in the limits of transmission) join
#    the adjacent region with which they share the longest border; flagged.
# 5. Stability: n_boot resamples of the current sites within each country
#    (the GLMM's type effects and variances held fixed), each clustered again
#    and cut at the same number of regions (or fewer, if fewer groups are
#    left after the minimum-size merges); how often each pair of countries
#    shares a region, how often each region is recovered intact, and the
#    adjusted Rand index against the chosen regions.
#
# Writes, in outputs/species_runs/regions/:
#   regions.csv       per country: region, flags (attached: no current
#                     bioassays; sparse: fewer than 3 current sites), the
#                     current level (logit, se, and the GLMM's), observed
#                     pooled mortality and counts in both windows, and the
#                     stability (the mean share of resamples that keep it with
#                     the other countries of its region, and its strongest
#                     link to a country of another region)
#   region_counts.csv per region: countries, level (logit, se), the share of
#                     resamples that recover it intact, observed pooled
#                     mortality and counts in both windows
#   merges.csv        the merge sequence: the two groups, their levels, the
#                     cost, the merged level and counts
#   partitions.csv    each country's group at every number of groups (to cut
#                     elsewhere)
#   coassignment.csv  per pair of countries, the share of resamples in which
#                     they share a region
#   glmm.rds          the GLMM (a cache)
# and in figures/species_runs/regions/:
#   regions_map.png   the regions in the limits of transmission, with the
#                     sites; each country's current observed mortality
#   merge_costs.png   merge cost by number of groups, and the cut
#   stability.png     the co-assignment matrix
# Plain R; under 1 GB and about 3 minutes, most of it the GLMM (cached in
# glmm.rds, after which under a minute).

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(sf)
  library(terra)
  library(lme4)
  library(ggtext)
  library(ggrepel)
})
source("R/functions.R")
# report()
source("R/species_fit_helpers.R")

arguments <- commandArgs(trailingOnly = TRUE)
n_regions <- if (length(arguments) > 0) as.integer(arguments[1]) else NA

llin_pyrethroids <- c("Alpha-cypermethrin", "Deltamethrin", "Permethrin")
windows <- list(early = 2010:2015, current = 2019:2024)
window_labels <- c(early = "2010-2015", current = "2019-2024")
# the minima per region: "about 100" bioassays in each window is taken as 90,
# so that a group a few bioassays short is not forced into a neighbour at a
# high cost (with 100: DR Congo, Uganda, Zambia and Angola, 98 bioassays in
# 2019-2024, into the East at a cost of 25, and Madagascar, 97, into the
# south at 31, leaving one region from Senegal to Somalia at 6 regions)
min_current <- 90
min_sites <- 10
min_early <- 90
prior_sd <- 2.5
n_boot <- 100
k_range <- 4:10
min_transmission <- 50
sparse_sites <- 3
output_dir <- "outputs/species_runs/regions"
figure_dir <- "figures/species_runs/regions"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dpi <- 120
options(width = 220)


# the bioassays ---------------------------------------------------------------------

bioassays <- read.csv("outputs/species_runs/misfit/ref_f0_bioassays.csv") %>%
  filter(insecticide_type %in% llin_pyrethroids) %>%
  transmute(country = country_name, cell, year = year_start,
            type = insecticide_type, died, tested = mosquito_number,
            window = case_when(year %in% windows$early ~ "early",
                               year %in% windows$current ~ "current"),
            site = paste(country, cell))
stopifnot(all(bioassays$died <= bioassays$tested),
          all(bioassays$died == round(bioassays$died)))
mask <- rast("data/clean/raster_mask.tif")
cell_xy <- terra::xyFromCell(mask, bioassays$cell)
bioassays$x <- cell_xy[, 1]
bioassays$y <- cell_xy[, 2]

# counts and observed pooled mortality per country and window
count_columns <- function(data) {
  data %>%
    summarise(bioassays = n(), sites = n_distinct(site),
              tested = sum(tested), observed = sum(died) / sum(tested),
              .groups = "drop")
}
country_counts <- bioassays %>%
  filter(!is.na(window)) %>%
  group_by(country, window) %>%
  count_columns() %>%
  pivot_wider(names_from = window,
              values_from = c(bioassays, sites, tested, observed),
              names_glue = "{.value}_{window}") %>%
  mutate(across(matches("^(bioassays|sites|tested)_"),
                ~ coalesce(as.numeric(.x), 0)))
current <- bioassays %>%
  filter(window == "current") %>%
  mutate(bioassay = factor(seq_len(n())), survived = tested - died,
         country_f = factor(country), type_f = factor(type))
report("current window: %d bioassays at %d sites in %d countries",
       nrow(current), n_distinct(current$site), n_distinct(current$country))


# 1. country levels ------------------------------------------------------------------

# the GLMM, cached in glmm_cache (refitted when the bioassays are newer)
contrasts(current$type_f) <- contr.sum(3)
input <- "outputs/species_runs/misfit/ref_f0_bioassays.csv"
glmm_cache <- file.path(output_dir, "glmm.rds")
if (file.exists(glmm_cache) && file.mtime(glmm_cache) > file.mtime(input)) {
  glmm <- readRDS(glmm_cache)
  time <- 0
} else {
  time <- system.time(
    glmm <- glmer(cbind(died, survived) ~ 0 + country_f + type_f +
                    (1 | site) + (1 | bioassay),
                  family = binomial, data = current,
                  control = glmerControl(optimizer = "bobyqa",
                                         optCtrl = list(maxfun = 1e5)))
  )[["elapsed"]]
  saveRDS(glmm, glmm_cache)
}
variances <- as.data.frame(VarCorr(glmm))
sd_site <- variances$sdcor[variances$grp == "site"]
sd_bioassay <- variances$sdcor[variances$grp == "bioassay"]
beta <- fixef(glmm)
type_levels <- levels(current$type_f)
type_raw <- c(beta[["type_f1"]], beta[["type_f2"]],
              -beta[["type_f1"]] - beta[["type_f2"]])
type_share <- tapply(current$tested, current$type_f, sum)[type_levels]
type_share <- type_share / sum(type_share)
type_offset <- setNames(type_raw - sum(type_share * type_raw), type_levels)
report("GLMM in %.0f s: sd site %.2f, sd bioassay %.2f; type offsets %s",
       time, sd_site, sd_bioassay,
       paste(sprintf("%s %+.2f", type_levels, type_offset), collapse = ", "))
# the GLMM's country levels at the window's type mix, with standard errors
country_names <- levels(current$country_f)
contrast_matrix <- matrix(0, length(country_names), length(beta),
                          dimnames = list(country_names, names(beta)))
for (i in seq_along(country_names)) {
  contrast_matrix[i, paste0("country_f", country_names[i])] <- 1
}
contrast_matrix[, "type_f1"] <- type_share[1] - type_share[3]
contrast_matrix[, "type_f2"] <- type_share[2] - type_share[3]
glmm_levels <- tibble(
  country = country_names,
  glmm_level = c(contrast_matrix %*% beta),
  glmm_se = sqrt(diag(contrast_matrix %*% as.matrix(vcov(glmm)) %*%
                        t(contrast_matrix))))

# The log likelihood of each site's current bioassays as a function of the
# level mu (on mu_grid), integrating each bioassay's binomial likelihood over
# its own effect, then the site's product over the site effect:
#   f_i(eta) = int Bin(died_i | tested_i, plogis(eta + e)) N(e; 0, sd_bioassay) de
#   L_s(mu)  = int prod_i f_i(mu + type_i + u) N(u; 0, sd_site) du
# both by convolution on x_grid (sums times the grid step)
x_grid <- seq(-20, 22, by = 0.025)
mu_grid <- seq(-8, 10, by = 0.02)
gaussian_kernel <- function(from, to, sd) {
  outer(from, to, function(a, b) dnorm(a - b, 0, sd)) * diff(from[1:2])
}
# log(sum(exp(log_values) * kernel)) per row, without underflow
log_convolve <- function(log_values, kernel) {
  row_max <- apply(log_values, 1, max)
  log(pmax(exp(log_values - row_max) %*% kernel, 1e-300)) + row_max
}
# the columns of m (values on x_grid) at x_grid + shift, linearly
# interpolated, constant beyond the ends
shift_columns <- function(m, shift) {
  position <- (x_grid + shift - x_grid[1]) / diff(x_grid[1:2]) + 1
  position <- pmin(pmax(position, 1), length(x_grid))
  low <- floor(position)
  high <- pmin(low + 1, length(x_grid))
  fraction <- rep(position - low, each = nrow(m))
  m[, low, drop = FALSE] * (1 - fraction) + m[, high, drop = FALSE] * fraction
}
site_curves <- function(data) {
  log_binomial <- lchoose(data$tested, data$died) +
    outer(data$died, plogis(x_grid, log.p = TRUE)) +
    outer(data$tested - data$died, plogis(-x_grid, log.p = TRUE))
  log_f <- log_convolve(log_binomial, gaussian_kernel(x_grid, x_grid,
                                                      sd_bioassay))
  shifted <- log_f
  for (type in type_levels) {
    rows <- data$type == type
    shifted[rows, ] <- shift_columns(log_f[rows, , drop = FALSE],
                                     type_offset[[type]])
  }
  by_site <- rowsum(shifted, data$site)
  curves <- log_convolve(by_site, gaussian_kernel(x_grid, mu_grid, sd_site))
  rownames(curves) <- rownames(by_site)
  curves
}
time <- system.time(curves <- site_curves(current))[["elapsed"]]
sites <- current %>%
  group_by(site) %>%
  summarise(country = first(country), bioassays = n(), .groups = "drop") %>%
  slice(match(rownames(curves), site))
report("site curves for %d sites in %.0f s", nrow(curves), time)

# the level (posterior mode) and precision of a summed curve; the prior is
# normal, prior_sd, centred on prior_mean (none if prior_sd is infinite)
mu_step <- diff(mu_grid[1:2])
level <- function(curve, sd = prior_sd) {
  if (is.finite(sd)) curve <- curve + dnorm(mu_grid, prior_mean, sd, log = TRUE)
  j <- min(max(which.max(curve), 2), length(mu_grid) - 1)
  second <- curve[j - 1] - 2 * curve[j] + curve[j + 1]
  c(m = mu_grid[j] + mu_step * (curve[j - 1] - curve[j + 1]) / (2 * second),
    w = -second / mu_step ^ 2)
}
prior_mean <- NA
prior_mean <- level(colSums(curves), Inf)[["m"]]
report("Africa's current level %.2f (%.0f%% at the median site); prior on each group's level N(%.2f, %.1f^2)",
       prior_mean, 100 * plogis(prior_mean), prior_mean, prior_sd)

clustered <- country_names
country_curves <- rowsum(curves, sites$country)[clustered, ]
levels_by_country <- t(apply(country_curves, 1, level))
estimates <- tibble(country = clustered, level = levels_by_country[, "m"],
                    se = 1 / sqrt(levels_by_country[, "w"])) %>%
  left_join(glmm_levels, by = "country") %>%
  left_join(country_counts, by = "country")
print(as.data.frame(estimates %>%
                      select(country, level, se, glmm_level, glmm_se,
                             observed_current, bioassays_current,
                             sites_current, bioassays_early)), digits = 2)
agree <- estimates %>% filter(sites_current >= sparse_sites)
report("against the GLMM (countries with >= %d sites): max |level diff| %.2f, se ratio %.2f-%.2f",
       sparse_sites, max(abs(agree$level - agree$glmm_level)),
       min(agree$se / agree$glmm_se), max(agree$se / agree$glmm_se))


# 2. adjacency -----------------------------------------------------------------------

borders <- readRDS("data/clean/country_borders.RDS")
transmission <- terra::freq(terra::mask(
  rast("data/clean/country_raster.tif"),
  rast("data/clean/pfpr_water_mask.tif")))
mapped <- union(unique(bioassays$country),
                transmission$value[transmission$count >= min_transmission])
borders <- borders %>% filter(country_name %in% mapped)
stopifnot(all(clustered %in% borders$country_name))
touching <- sf::st_is_within_distance(borders, borders, dist = 1000)
adjacency <- matrix(FALSE, nrow(borders), nrow(borders),
                    dimnames = list(borders$country_name, borders$country_name))
for (i in seq_along(touching)) adjacency[i, touching[[i]]] <- TRUE
diag(adjacency) <- FALSE
# islands to their nearest mainland country
islands <- names(which(rowSums(adjacency) == 0))
mainland <- setdiff(borders$country_name, islands)
distances <- sf::st_distance(borders[match(islands, borders$country_name), ],
                             borders[match(mainland, borders$country_name), ])
island_links <- tibble(island = islands,
                       mainland = mainland[apply(distances, 1, which.min)],
                       km = apply(distances, 1, min) / 1000)
for (i in seq_len(nrow(island_links))) {
  adjacency[island_links$island[i], island_links$mainland[i]] <- TRUE
  adjacency[island_links$mainland[i], island_links$island[i]] <- TRUE
}
report("islands linked: %s", paste(sprintf("%s to %s (%.0f km)",
                                           island_links$island,
                                           island_links$mainland,
                                           island_links$km), collapse = "; "))
# connected components of an adjacency matrix
components <- function(a) {
  group <- rep(NA_integer_, nrow(a))
  k <- 0
  for (start in seq_len(nrow(a))) {
    if (!is.na(group[start])) next
    k <- k + 1
    frontier <- start
    while (length(frontier) > 0) {
      group[frontier] <- k
      frontier <- setdiff(which(colSums(a[frontier, , drop = FALSE]) > 0),
                          which(!is.na(group)))
    }
  }
  setNames(group, rownames(a))
}
clustered_adjacency <- adjacency[clustered, clustered]
stopifnot(max(components(clustered_adjacency)) == 1)
report("adjacency: %d countries mapped, %d with current bioassays, %d links among them, connected",
       nrow(borders), length(clustered), sum(clustered_adjacency) / 2)

# shared border length between adjacent countries (km of each one's boundary
# within 2 km of the other), for attaching the countries without current data
sf::sf_use_s2(TRUE)
border_lines <- sf::st_boundary(borders)
shared_border <- function(a, b) {
  i <- match(a, borders$country_name)
  j <- match(b, borders$country_name)
  part <- suppressWarnings(sf::st_intersection(
    sf::st_geometry(border_lines)[i], sf::st_buffer(sf::st_geometry(borders)[j], 2000)))
  if (length(part) == 0) return(0)
  as.numeric(sum(sf::st_length(part))) / 1000
}


# 3. clustering ----------------------------------------------------------------------

# Spatially constrained agglomeration of the rows of `curves` (one per
# country, on mu_grid) with counts `counts` (columns current, sites, early)
# and adjacency `a`: the merges, in order, and the partition at every number
# of groups (a matrix, one row per number of groups, from the start down to
# 1, of each country's group, numbered by first appearance)
cluster <- function(curves, counts, a) {
  n <- nrow(curves)
  group <- setNames(seq_len(n), rownames(curves))
  estimate <- t(apply(curves, 1, level))
  partitions <- matrix(NA_integer_, n, n,
                       dimnames = list(seq_len(n), rownames(curves)))
  partitions[n, ] <- group
  merges <- vector("list", n - 1)
  alive <- rep(TRUE, n)
  for (step in seq_len(n - 1)) {
    small <- alive & (counts[, "current"] < min_current |
                        counts[, "sites"] < min_sites |
                        counts[, "early"] < min_early)
    pairs <- which(a & upper.tri(a), arr.ind = TRUE)
    if (any(small)) {
      pairs <- pairs[small[pairs[, 1]] | small[pairs[, 2]], , drop = FALSE]
    }
    stopifnot(nrow(pairs) > 0)
    w_g <- estimate[pairs[, 1], "w"]
    w_h <- estimate[pairs[, 2], "w"]
    cost <- w_g * w_h / (w_g + w_h) *
      (estimate[pairs[, 1], "m"] - estimate[pairs[, 2], "m"]) ^ 2
    best <- which.min(cost)
    g <- pairs[best, 1]
    h <- pairs[best, 2]
    before <- estimate[c(g, h), ]
    curves[g, ] <- curves[g, ] + curves[h, ]
    counts[g, ] <- counts[g, ] + counts[h, ]
    a[g, ] <- a[g, ] | a[h, ]
    a[, g] <- a[, g] | a[, h]
    a[h, ] <- FALSE
    a[, h] <- FALSE
    a[g, g] <- FALSE
    alive[h] <- FALSE
    members_g <- names(group)[group == g]
    members_h <- names(group)[group == h]
    group[group == h] <- g
    estimate[g, ] <- level(curves[g, ])
    estimate[h, ] <- NA
    partitions[n - step, ] <- match(group, unique(group))
    merges[[step]] <- tibble(
      step = step, groups_after = n - step,
      phase = if (any(small)) "minimum size" else "Ward",
      group_1 = paste(members_g, collapse = ", "),
      group_2 = paste(members_h, collapse = ", "),
      level_1 = before[1, "m"], se_1 = 1 / sqrt(before[1, "w"]),
      level_2 = before[2, "m"], se_2 = 1 / sqrt(before[2, "w"]),
      cost = cost[best],
      merged_level = estimate[g, "m"], merged_se = 1 / sqrt(estimate[g, "w"]),
      merged_current = counts[g, "current"], merged_sites = counts[g, "sites"],
      merged_early = counts[g, "early"])
  }
  list(merges = bind_rows(merges), partitions = partitions)
}

# The cut: within k_range, the number of groups K before the merge whose cost
# most exceeds every Ward merge before it (the largest ratio of the cost of
# going from K to K - 1 groups to the largest cost of the merges above K)
choose_cut <- function(merges) {
  ward <- merges %>% filter(phase == "Ward") %>% arrange(desc(groups_after))
  ward$k <- ward$groups_after + 1
  ward$ratio <- ward$cost / c(NA, cummax(ward$cost)[-nrow(ward)])
  candidates <- ward %>% filter(k %in% k_range, !is.na(ratio))
  candidates$k[which.max(candidates$ratio)]
}

count_matrix <- function(current_bioassays, current_sites) {
  cbind(current = current_bioassays, sites = current_sites,
        early = estimates$bioassays_early)
}
counts <- count_matrix(estimates$bioassays_current, estimates$sites_current)
rownames(counts) <- clustered
result <- cluster(country_curves, counts, clustered_adjacency)
merges <- result$merges
first_ward <- min(merges$step[merges$phase == "Ward"])
k_start <- merges$groups_after[first_ward] + 1
if (is.na(n_regions)) n_regions <- choose_cut(merges)
stopifnot(n_regions <= k_start)
print(as.data.frame(merges %>% select(step, groups_after, phase, group_1,
                                      group_2, level_1, level_2, cost,
                                      merged_current, merged_sites,
                                      merged_early)), digits = 3,
      right = FALSE)
report("%d minimum-size merges leave %d groups; cut at %d regions",
       first_ward - 1, k_start, n_regions)
partition <- result$partitions[n_regions, ]


# 4. countries without current data, and the region names ---------------------------

# the shared border length of each country without current bioassays with
# each of its neighbours (0 for an island link)
open_countries <- setdiff(borders$country_name, clustered)
neighbour_lengths <- lapply(setNames(nm = open_countries), function(country) {
  neighbours <- names(which(adjacency[country, ]))
  setNames(vapply(neighbours, function(b) shared_border(country, b),
                  numeric(1)), neighbours)
})
# Every mapped country's group, from `partition` (the countries with current
# bioassays): each other country joins the group with which it shares the
# longest border among its assigned neighbours, in passes until all are
# assigned
attach_open <- function(partition) {
  assignment <- setNames(rep(NA_integer_, nrow(borders)),
                         borders$country_name)
  assignment[names(partition)] <- partition
  repeat {
    open <- names(assignment)[is.na(assignment)]
    if (length(open) == 0) break
    joined <- vapply(open, function(country) {
      lengths <- neighbour_lengths[[country]]
      known <- names(lengths)[!is.na(assignment[names(lengths)])]
      if (length(known) == 0) return(NA_integer_)
      by_group <- tapply(lengths[known], assignment[known], sum)
      as.integer(names(by_group)[which.max(by_group)])
    }, integer(1))
    stopifnot(any(!is.na(joined)))
    assignment[open] <- joined
  }
  assignment
}
assignment <- attach_open(partition)
attached <- open_countries
attach_notes <- tibble(
  country = open_countries,
  attached_by = vapply(neighbour_lengths, function(lengths) {
    paste(sprintf("%s %.0f km", names(lengths), lengths), collapse = ", ")
  }, character(1)))

# regions lettered west to east by the mean longitude of their current
# sites, weighted by mosquitoes tested
order_west_east <- current %>%
  mutate(group = assignment[country]) %>%
  group_by(group) %>%
  summarise(x = sum(x * tested) / sum(tested), .groups = "drop") %>%
  arrange(x)
region_letter <- setNames(LETTERS[seq_len(nrow(order_west_east))],
                          order_west_east$group)
regions <- tibble(country = names(assignment),
                  region = unname(region_letter[as.character(assignment)]),
                  attached = country %in% attached) %>%
  left_join(estimates %>% select(country, level, se, glmm_level, glmm_se),
            by = "country") %>%
  left_join(country_counts, by = "country") %>%
  mutate(across(matches("^(bioassays|sites|tested)_"),
                ~ coalesce(.x, 0)),
         sparse = !attached & sites_current < sparse_sites) %>%
  left_join(attach_notes, by = "country")
# a region's name: its letter and its countries with current bioassays, the
# most bioassays first
region_names <- regions %>%
  filter(!attached) %>%
  arrange(region, desc(bioassays_current)) %>%
  group_by(region) %>%
  summarise(name = paste0(first(region), ": ",
                          paste(country, collapse = ", ")),
            .groups = "drop")
regions$region_name <- region_names$name[match(regions$region,
                                               region_names$region)]


# 5. stability -----------------------------------------------------------------------

site_index <- split(seq_len(nrow(sites)), sites$country)[clustered]
coassignment <- matrix(0, nrow(borders), nrow(borders),
                       dimnames = list(borders$country_name,
                                       borders$country_name))
adjusted_rand <- function(a, b) {
  table_ab <- table(a, b)
  pairs <- function(n) sum(n * (n - 1) / 2)
  index <- pairs(table_ab)
  expected <- pairs(rowSums(table_ab)) * pairs(colSums(table_ab)) /
    pairs(length(a))
  maximum <- (pairs(rowSums(table_ab)) + pairs(colSums(table_ab))) / 2
  (index - expected) / (maximum - expected)
}
set.seed(47)
rand <- numeric(n_boot)
boot_k <- integer(n_boot)
# each chosen region's countries with current bioassays, to count the
# resamples that recover it intact
group_key <- function(p) {
  vapply(split(names(p), p), function(m) paste(sort(m), collapse = "|"),
         character(1))
}
chosen_keys <- group_key(partition)
recovered <- setNames(numeric(length(chosen_keys)),
                      region_letter[names(chosen_keys)])
time <- system.time(for (b in seq_len(n_boot)) {
  weight <- numeric(nrow(sites))
  for (index in site_index) {
    drawn <- index[sample.int(length(index), length(index), replace = TRUE)]
    weight <- weight + tabulate(drawn, nrow(sites))
  }
  boot_curves <- rowsum(curves * weight, sites$country)[clustered, ]
  boot_counts <- count_matrix(
    c(rowsum(weight * sites$bioassays, sites$country)[clustered, ]),
    estimates$sites_current)
  rownames(boot_counts) <- clustered
  boot <- cluster(boot_curves, boot_counts, clustered_adjacency)
  boot_start <- boot$merges$groups_after[
    min(boot$merges$step[boot$merges$phase == "Ward"])] + 1
  boot_k[b] <- min(n_regions, boot_start)
  boot_partition <- boot$partitions[boot_k[b], ]
  rand[b] <- adjusted_rand(boot_partition, partition)
  recovered <- recovered + (chosen_keys %in% group_key(boot_partition))
  full <- attach_open(boot_partition)
  coassignment <- coassignment + outer(full, full, "==")
})[["elapsed"]]
coassignment <- coassignment / n_boot
recovered <- recovered / n_boot
report("%d resamples in %.0f s; adjusted Rand index against the chosen regions: median %.2f (%.2f-%.2f, 10%%-90%%)%s",
       n_boot, time, median(rand), quantile(rand, 0.1), quantile(rand, 0.9),
       if (any(boot_k < n_regions)) sprintf(
         "; %d resamples with fewer than %d groups after the minimum-size merges",
         sum(boot_k < n_regions), n_regions) else "")

# per country: the mean co-assignment with the other countries of its region,
# and the country of another region it most often joins
stability <- bind_rows(lapply(regions$country, function(country) {
  same <- setdiff(regions$country[regions$region ==
                                    regions$region[regions$country == country]],
                  country)
  other <- setdiff(regions$country, c(same, country))
  links <- coassignment[country, other]
  tibble(country = country,
         with_region = if (length(same) > 0) mean(coassignment[country, same])
           else NA,
         strongest_other = other[which.max(links)],
         strongest_other_share = max(links))
}))
regions <- regions %>% left_join(stability, by = "country")


# outputs ----------------------------------------------------------------------------

regions <- regions %>%
  arrange(region, attached, desc(bioassays_current)) %>%
  select(country, region, region_name, attached, sparse, level, se,
         glmm_level, glmm_se, observed_current, bioassays_current,
         sites_current, tested_current, observed_early, bioassays_early,
         sites_early, tested_early, with_region, strongest_other,
         strongest_other_share, attached_by)
write.csv(regions, file.path(output_dir, "regions.csv"), row.names = FALSE)

region_levels <- t(vapply(sort(unique(regions$region)), function(r) {
  countries <- intersect(regions$country[regions$region == r], clustered)
  level(colSums(country_curves[countries, , drop = FALSE]))
}, numeric(2)))
region_counts <- bioassays %>%
  filter(!is.na(window)) %>%
  mutate(region = regions$region[match(country, regions$country)]) %>%
  group_by(region, window) %>%
  count_columns() %>%
  pivot_wider(names_from = window,
              values_from = c(bioassays, sites, tested, observed),
              names_glue = "{.value}_{window}") %>%
  left_join(regions %>%
              group_by(region) %>%
              summarise(name = first(region_name),
                        countries = paste(country, collapse = ", "),
                        attached = paste(country[attached], collapse = ", "),
                        .groups = "drop"),
            by = "region") %>%
  mutate(level = region_levels[region, "m"],
         se = 1 / sqrt(region_levels[region, "w"]),
         recovered = unname(recovered[region])) %>%
  select(region, name, countries, attached, level, se, recovered,
         observed_current,
         bioassays_current, sites_current, tested_current, observed_early,
         bioassays_early, sites_early, tested_early)
write.csv(region_counts, file.path(output_dir, "region_counts.csv"),
          row.names = FALSE)
write.csv(merges, file.path(output_dir, "merges.csv"), row.names = FALSE)
partitions <- as.data.frame(t(result$partitions[k_start:1, ]))
names(partitions) <- paste0("groups_", k_start:1)
write.csv(cbind(country = clustered,
                region = regions$region[match(clustered, regions$country)],
                partitions),
          file.path(output_dir, "partitions.csv"), row.names = FALSE)
coassignment_long <- as.data.frame(as.table(coassignment)) %>%
  setNames(c("country_1", "country_2", "share")) %>%
  filter(as.character(country_1) < as.character(country_2))
write.csv(coassignment_long, file.path(output_dir, "coassignment.csv"),
          row.names = FALSE)
print(as.data.frame(region_counts %>% select(-countries)), digits = 3,
      right = FALSE)
print(as.data.frame(regions %>%
                      select(country, region, attached, sparse, level, se,
                             observed_current, bioassays_current,
                             sites_current, bioassays_early, with_region,
                             strongest_other, strongest_other_share)),
      digits = 2)


# figures ----------------------------------------------------------------------------

# distinct colours for neighbouring regions (Okabe-Ito and Tol)
region_colours <- setNames(
  c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#882255", "#D55E00",
    "#CC79A7", "#117733", "#999933", "#332288")[seq_len(n_regions)],
  LETTERS[seq_len(n_regions)])
annotation_theme <- theme(plot.caption = element_text(hjust = 0, size = 10),
                          plot.title = element_text(size = 15))

# the limits of transmission without water bodies, coarsened by 4, each cell
# with its country
aggregation <- 4
water_mask <- terra::aggregate(rast("data/clean/pfpr_water_mask.tif"),
                               aggregation, fun = "max", na.rm = TRUE)
country_raster <- terra::aggregate(rast("data/clean/country_raster.tif"),
                                   aggregation, fun = "modal", na.rm = TRUE)
tiles <- as.data.frame(c(water_mask, country_raster), xy = TRUE,
                       na.rm = TRUE) %>%
  setNames(c("x", "y", "mask", "country"))
tiles$country <- as.character(tiles$country)
tiles <- tiles %>%
  filter(country %in% regions$country) %>%
  left_join(regions %>% select(country, region, attached, observed_current),
            by = "country")
tile <- res(water_mask)
sf::sf_use_s2(FALSE)
all_borders <- readRDS("data/clean/country_borders.RDS")
region_outlines <- borders %>%
  mutate(region = regions$region[match(country_name, regions$country)]) %>%
  group_by(region) %>%
  summarise(.groups = "drop")
on_surface <- function(polygons) {
  points <- suppressWarnings(sf::st_point_on_surface(polygons))
  xy <- sf::st_coordinates(points)
  sf::st_drop_geometry(points) %>% mutate(x = xy[, 1], y = xy[, 2])
}
region_centres <- on_surface(region_outlines)
country_labels <- on_surface(borders) %>%
  select(country = country_name, x, y) %>%
  left_join(regions, by = "country") %>%
  mutate(label = sprintf("%s %.0f%% (%d)", country, 100 * observed_current,
                         as.integer(bioassays_current)))
site_points <- bioassays %>%
  filter(!is.na(window)) %>%
  group_by(site) %>%
  summarise(x = first(x), y = first(y),
            current = any(window == "current"), .groups = "drop")
xlim <- c(-18, 52)
ylim <- c(-35, 23)
base_map <- function(layers, title) {
  ggplot() +
    geom_sf(data = all_borders, fill = grey(0.95), colour = NA) +
    layers +
    geom_sf(data = all_borders, fill = NA, colour = grey(0.35),
            linewidth = 0.2) +
    geom_sf(data = region_outlines, fill = NA, colour = grey(0.05),
            linewidth = 0.75) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_ir_maps() +
    theme(plot.title = element_markdown(size = 13),
          panel.background = element_rect(fill = "#F4F8FB", colour = NA))
}
region_panel <- base_map(
  list(geom_tile(aes(x, y, fill = region), data = tiles, width = tile[1],
                 height = tile[2], alpha = 0.85),
       geom_tile(aes(x, y), data = filter(tiles, attached), fill = "white",
                 alpha = 0.5, width = tile[1], height = tile[2]),
       geom_point(aes(x, y), data = filter(site_points, !current),
                  shape = 1, size = 0.8, colour = grey(0.3), stroke = 0.3),
       geom_point(aes(x, y), data = filter(site_points, current),
                  shape = 16, size = 1.1, colour = "black"),
       geom_label_repel(aes(x, y, label = region, fill = region),
                        data = region_centres, size = 7, fontface = "bold",
                        label.padding = unit(3, "pt"), seed = 1,
                        min.segment.length = 0.3, box.padding = 0.6,
                        show.legend = FALSE)),
  sprintf("%d regions from current (2019-2024) LLIN-pyrethroid susceptibility",
          n_regions)) +
  scale_fill_manual(values = region_colours, guide = "none")
mortality_panel <- base_map(
  list(geom_tile(aes(x, y, fill = observed_current), data = tiles,
                 width = tile[1], height = tile[2]),
       geom_label_repel(aes(x, y, label = label),
                        data = filter(country_labels, !attached),
                        size = 2.9, fill = alpha("white", 0.8),
                        label.size = 0, label.padding = unit(1.2, "pt"),
                        min.segment.length = 0.2, box.padding = 0.15,
                        segment.colour = grey(0.2), segment.size = 0.3,
                        max.overlaps = Inf, seed = 1)),
  "Observed pooled mortality per country, 2019-2024 (bioassays)") +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(11, "RdYlBu"),
                       limits = c(0, 1), labels = scales::percent,
                       na.value = grey(0.8),
                       name = "pooled<br>mortality<br>2019-2024") +
  theme(legend.position = "inside",
        legend.position.inside = c(0.1, 0.3),
        legend.title = element_markdown(size = 11),
        legend.text = element_text(size = 10),
        legend.key.height = unit(28, "pt"))
# the key: per region, its countries and counts
key <- region_counts %>%
  mutate(index = seq_len(n()) - 1,
         column = index %/% ceiling(n_regions / 3),
         row = index %% ceiling(n_regions / 3),
         countries = sub("^[A-Z]: ", "", name),
         extra = if_else(attached == "", "",
                         paste0(" (+ ", attached, "*)")),
         text = sprintf(
           "**%s**: %s<br>2019-2024: %.0f%% (%d bioassays, %d sites); 2010-2015: %.0f%% (%d)",
           region, vapply(paste0(countries, extra), function(s) {
             paste(strwrap(s, 62), collapse = "<br>")
           }, character(1)),
           100 * observed_current, as.integer(bioassays_current),
           as.integer(sites_current), 100 * observed_early,
           as.integer(bioassays_early)))
key_panel <- ggplot(key) +
  geom_tile(aes(column * 10, -row, fill = region), width = 0.35,
            height = 0.5, colour = grey(0.2)) +
  geom_richtext(aes(column * 10 + 0.3, -row, label = text), hjust = 0,
                size = 3.3, fill = NA, label.colour = NA, lineheight = 1.1) +
  scale_fill_manual(values = region_colours, guide = "none") +
  scale_x_continuous(limits = c(-0.3, 30)) +
  scale_y_continuous(limits = c(-max(key$row) - 0.5, 0.5)) +
  theme_void()
map_figure <- (region_panel + mortality_panel) / key_panel +
  plot_layout(heights = c(1, 0.25)) +
  plot_annotation(
    title = "Diagnostic regions: neighbouring countries with alike current LLIN-pyrethroid susceptibility",
    caption = paste0(
      "Regions: spatially constrained clustering of the countries with ",
      "2019-2024 LLIN-pyrethroid bioassays (alpha-cypermethrin, ",
      "deltamethrin, permethrin), at a precision-weighted Ward cost on each ",
      "group's current level\n(GLMM: country + type + site + bioassay ",
      "effects); each region has at least ", min_current, " bioassays and ",
      min_sites, " sites in 2019-2024 and ", min_early, " bioassays in ",
      "2010-2015. Paler (grey on the right): no 2019-2024 bioassays,\n",
      "joined to the neighbouring region with the longest border (* in the ",
      "key). Points: sites (5 km cells) with bioassays in 2019-2024 ",
      "(filled) or only in 2010-2015 (open). Coloured: the limits of ",
      "transmission\nwithout water bodies. Right: observed mortality pooled ",
      "over the country's bioassays (died / tested), with the number of ",
      "bioassays. Key: observed pooled mortality and counts per window."),
    theme = annotation_theme)
ggsave(file.path(figure_dir, "regions_map.png"), map_figure, width = 17,
       height = 10.8, dpi = dpi, bg = "white")

# merge costs by the number of groups before the merge; a merged group is
# named by its regions at the cut if it is a union of them, else by its
# countries
group_label <- function(group) {
  members <- strsplit(group, ", ")[[1]]
  letters_in <- unique(regions$region[match(members, regions$country)])
  whole <- all(vapply(letters_in, function(r) {
    all(intersect(regions$country[regions$region == r], clustered) %in%
          members)
  }, logical(1)))
  if (whole) return(paste(sort(letters_in), collapse = ""))
  if (length(members) > 3) {
    members <- c(members[1:2], sprintf("+%d", length(members) - 2))
  }
  paste(members, collapse = ", ")
}
cost_data <- merges %>%
  mutate(k = groups_after + 1,
         label = if_else(phase == "Ward",
                         paste(vapply(group_1, group_label, character(1)),
                               "+", vapply(group_2, group_label,
                                           character(1))),
                         NA_character_))
cost_figure <- ggplot(cost_data, aes(k, cost)) +
  geom_hline(yintercept = qchisq(0.95, 1), linetype = "dotted",
             colour = grey(0.4)) +
  geom_vline(xintercept = n_regions + 0.5, colour = "#D55E00",
             linewidth = 0.6) +
  geom_line(data = filter(cost_data, phase == "Ward"), colour = grey(0.5)) +
  geom_point(aes(colour = phase), size = 3) +
  geom_text_repel(aes(label = label), size = 3.6, na.rm = TRUE,
                  min.segment.length = 0, max.overlaps = Inf,
                  box.padding = 0.5, seed = 1, direction = "y",
                  nudge_x = -1.5, hjust = 0) +
  annotate("text", x = n_regions + 0.7, y = min(cost_data$cost),
           label = sprintf("cut: %d regions", n_regions), hjust = 1,
           vjust = 0, colour = "#D55E00", size = 4.2) +
  annotate("text", x = max(cost_data$k), y = qchisq(0.95, 1), vjust = -0.5,
           hjust = 0, label = "3.84: chi-squared 5% (1 df)", size = 3.6,
           colour = grey(0.35)) +
  scale_x_reverse(breaks = seq(1, max(cost_data$k), by = 1)) +
  scale_y_log10() +
  scale_colour_manual(values = c("minimum size" = "#56B4E9",
                                 Ward = "black"), name = "merge") +
  labs(x = "number of groups before the merge",
       y = "merge cost (log scale)",
       title = "Merge cost of the spatially constrained clustering",
       caption = paste(
         "Cost: w_g w_h / (w_g + w_h) (m_g - m_h)^2, the precision-weighted",
         "Ward cost of merging adjacent groups g and h with current levels m",
         "(logit) and precisions w;\nabout the chi-squared statistic for",
         "the two having the same level. Blue: merges of groups below the",
         "minimum size, first (the cheapest such merge each time). Labels:",
         "the groups merged, by\nregion letters at the cut where they are",
         "unions of regions. The cut keeps the groups before the largest",
         "jump in cost relative to every earlier merge.")) +
  theme_minimal(base_size = 13) +
  theme(plot.caption = element_text(hjust = 0, size = 10),
        panel.grid.minor = element_blank(), legend.position = "top")
ggsave(file.path(figure_dir, "merge_costs.png"), cost_figure, width = 13,
       height = 7.5, dpi = dpi, bg = "white")

# the co-assignment matrix, countries by region and level
country_order <- regions %>%
  arrange(region, attached, level) %>%
  mutate(label = sprintf("%s  %s%s", region, country,
                         if_else(attached, "*", ""))) %>%
  select(country, label)
stability_data <- as.data.frame(as.table(coassignment)) %>%
  setNames(c("country_1", "country_2", "share")) %>%
  mutate(country_1 = factor(country_order$label[match(country_1,
                                                      country_order$country)],
                            country_order$label),
         country_2 = factor(country_order$label[match(country_2,
                                                      country_order$country)],
                            rev(country_order$label)))
boundaries <- cumsum(table(regions$region[match(country_order$country,
                                                regions$country)]))
stability_figure <- ggplot(stability_data,
                           aes(country_1, country_2, fill = share)) +
  geom_tile() +
  geom_vline(xintercept = boundaries[-length(boundaries)] + 0.5,
             linewidth = 0.4) +
  geom_hline(yintercept = length(country_order$label) -
               boundaries[-length(boundaries)] + 0.5, linewidth = 0.4) +
  scale_fill_gradient(low = "white", high = "#08306B", limits = c(0, 1),
                      labels = scales::percent,
                      name = "share of\nresamples\nin the same\nregion") +
  labs(x = NULL, y = NULL,
       title = sprintf("Stability of the %d regions: %d resamples of the sites within each country",
                       n_regions, n_boot),
       caption = sprintf(paste0(
         "Countries by region (A-%s west to east; lines), then those without ",
         "2019-2024 bioassays (*), then by current level. Each resample is ",
         "clustered again and cut at %d regions.\nAdjusted Rand index of ",
         "the resampled regions against the chosen ones: median %.2f ",
         "(10%%-90%%: %.2f-%.2f)."),
         LETTERS[n_regions], n_regions, median(rand), quantile(rand, 0.1),
         quantile(rand, 0.9))) +
  coord_equal() +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = element_blank(),
        plot.caption = element_text(hjust = 0, size = 10))
ggsave(file.path(figure_dir, "stability.png"), stability_figure, width = 13,
       height = 12.5, dpi = dpi, bg = "white")
report("written %s and %s; peak memory %.1f GB",
       toString(file.path(output_dir, c("regions.csv", "region_counts.csv",
                                        "merges.csv", "partitions.csv",
                                        "coassignment.csv"))),
       toString(file.path(figure_dir, c("regions_map.png", "merge_costs.png",
                                        "stability.png"))),
       peak_memory_gb())
