# A visual check (#46): the plateau map of a fit (the change in its pyrethroid
# level residual between two windows, from R/plateau_map.R on the
# net-transform-screen branch) beside simple surfaces from the kdr records
# (data/clean/kdr_records.csv, R/prep_kdr.R): the allele frequencies of
# Vgsc-995F and 995S, their sum (total kdr), and the share of An. arabiensis
# among the samples genotyped to species.
#
# Each surface is a quasi-binomial GAM (overdispersion as a scale, not a
# random effect of citation: studies are local, so a study effect absorbs the
# spatial pattern) with a spatial smooth and a linear year trend, predicted
# at 2015, and blank more than 2.5 degrees from that surface's records.
#
#   Rscript R/kdr_plateau_compare.R <plateau grid csv> <label>
#
# Writes figures/kdr_plateau_compare_<label>.png.

arguments <- commandArgs(trailingOnly = TRUE)
plateau_file <- arguments[1]
label <- arguments[2]

source("R/packages.R")
suppressMessages(library(mgcv))

grid <- read.csv(plateau_file)
records <- read.csv("data/clean/kdr_records.csv") |>
  mutate(citation = factor(citation), year_c = year_start - 2015)

# blank grid cells more than max_distance degrees from any of xy
near_records <- function(xy, max_distance = 2.5) {
  vapply(seq_len(nrow(grid)), function(i) {
    min((xy$longitude - grid$lon[i]) ^ 2 + (xy$latitude - grid$lat[i]) ^ 2)
  }, numeric(1)) <= max_distance ^ 2
}

# a binomial surface of successes out of trials, at 2015
fit_surface <- function(data) {
  model <- bam(cbind(successes, trials - successes) ~
                 s(longitude, latitude, bs = "tp", k = 80) + year_c,
               family = stats::quasibinomial(), data = data, method = "fREML")
  new <- data.frame(longitude = grid$lon, latitude = grid$lat, year_c = 0)
  p <- predict(model, new, type = "response")
  p[!near_records(data)] <- NA
  as.numeric(p)
}

allele_data <- function(v) {
  records |>
    filter(variant == v, !is.na(n_tested), n_tested > 0) |>
    mutate(trials = 2 * n_tested, successes = round(frequency * trials))
}
grid$f <- fit_surface(allele_data("995F"))
grid$s <- fit_surface(allele_data("995S"))
grid$total <- pmin(grid$f + grid$s, 1)

species_data <- records |>
  filter(species_group %in% c("gambiae_ss", "coluzzii", "arabiensis")) |>
  distinct(sample_id, .keep_all = TRUE) |>
  mutate(trials = 1, successes = as.numeric(species_group == "arabiensis"))
grid$arabiensis <- fit_surface(species_data)

panel <- function(column, title, diverging = FALSE) {
  p <- ggplot(grid, aes(lon, lat, fill = .data[[column]])) +
    geom_raster() +
    coord_equal() +
    labs(title = title, x = NULL, y = NULL, fill = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.text = element_blank())
  if (diverging) {
    p + scale_fill_distiller(palette = "RdBu", limits = c(-2, 2),
                             oob = scales::squish, na.value = "grey90")
  } else {
    p + scale_fill_viridis_c(limits = c(0, 1), na.value = "grey90")
  }
}
panels <- list(
  panel("change", "Missing plateau (red), level residual\n2012-16 to 2019-25",
        diverging = TRUE),
  panel("f", "Vgsc-995F frequency, 2015"),
  panel("s", "Vgsc-995S frequency, 2015"),
  panel("total", "Total kdr (995F + 995S), 2015"),
  panel("arabiensis", "Share An. arabiensis\n(samples genotyped to species)"))
p <- patchwork::wrap_plots(panels, ncol = 3) +
  patchwork::plot_annotation(
    title = sprintf("Plateau misfit (%s) and kdr records", label),
    subtitle = paste("Quasi-binomial GAMs, at 2015;",
                     "grey: more than 2.5 degrees from data"))
ggsave(sprintf("figures/kdr_plateau_compare_%s.png", label), p, width = 13,
       height = 9, dpi = 120, bg = "white")
