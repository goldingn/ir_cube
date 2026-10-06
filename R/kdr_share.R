# The share of kdr resistance alleles that are Vgsc-995S (L1014S, the East
# African kdr) rather than Vgsc-995F (L1014F, West African), as a candidate
# covariate for regional differences in net selection (#37), and a check of
# whether it explains the fit's net-use misfit.
#
#   Rscript R/kdr_share.R <fit RData> <label> [gap]
#
# From the Vector Atlas genotypic records (data/raw VA_FULL_DATA_20250721.csv):
# a binomial GAM of each allele's frequency on longitude and latitude (a
# thin-plate smooth), from the records that typed it (resistant-allele counts:
# frequency x 2 x the number tested), and the share s(x) = S / (F + S) of the
# two surfaces at each pixel. Then,
# for pyrethroids on the fit:
#   1. level residuals (empirical minus predicted logit of mortality) on
#      cumulative net use and its interaction with s, with country fixed
#      effects;
#   2. the change misfit of within-pixel pairs (observed minus predicted fall
#      in logit mortality per year) on mean net use and its interaction with
#      s, weighted by precision.
# A negative interaction in 2 (positive in 1) means the fit over-predicts net
# selection more where 995S dominates, which a term exp(beta + gamma s) on
# net selection could absorb.
#
# Writes outputs/hinge_screen/kdr_share_records.csv, kdr_995s_share.tif,
# <label>_kdr_share.csv and figures/kdr_share_map.png.

arguments <- commandArgs(trailingOnly = TRUE)
fit_file <- arguments[1]
label <- arguments[2]
gap <- if (length(arguments) >= 3) as.integer(arguments[3]) else 5L
va_file <- "/home/nick/Dropbox/github/ir_cube/data/raw/VA_FULL_DATA_20250721.csv"

source("R/packages.R")
source("R/functions.R")
source("R/dynamical_model.R")
source("R/dynamical_predictions.R")
suppressMessages(library(mgcv))

# kdr records -----------------------------------------------------------------
# Studies often type one allele only (995F in the West, 995S in the East), so
# a missing allele is unknown, not 0: each allele's frequency is smoothed from
# the records that typed it, and the share is S / (F + S) of the two surfaces
va <- read.csv(va_file)
genotypic <- va |>
  filter(insecticide_resistance_data %in% c("genotypic", "both")) |>
  transmute(country, year = year_start,
            latitude = latitude_1, longitude = longitude_1,
            tested = vgsc_number_of_mosquitoes_tested,
            f = as.numeric(vgsc995f_percent) / 100,
            s = as.numeric(vgsc995s_percent) / 100) |>
  filter(!is.na(tested), tested > 0, !is.na(latitude), !is.na(longitude)) |>
  # coordinates within Africa's bounding box (one record has a longitude of
  # about -43,000)
  filter(between(longitude, -26, 64), between(latitude, -36, 38)) |>
  mutate(alleles = 2 * tested)
allele_fit <- function(frequency) {
  r <- genotypic[!is.na(genotypic[[frequency]]), ]
  r$resistant <- round(r[[frequency]] * r$alleles)
  r$susceptible <- r$alleles - r$resistant
  list(records = r,
       fit = gam(cbind(resistant, susceptible) ~
                   s(longitude, latitude, k = 40),
                 family = quasibinomial(), data = r, method = "REML"))
}
allele <- list(f = allele_fit("f"), s = allele_fit("s"))
records <- bind_rows(lapply(names(allele), function(a) {
  mutate(allele[[a]]$records, allele = a)
}))
write.csv(records, "outputs/hinge_screen/kdr_share_records.csv",
          row.names = FALSE)
for (a in names(allele)) {
  r <- allele[[a]]$records
  cat(sprintf("995%s: %d records, %d countries; ", toupper(a), nrow(r),
              n_distinct(r$country)))
}
cat("\n")
predict_share <- function(longitude, latitude) {
  new <- data.frame(longitude = longitude, latitude = latitude)
  pf <- as.vector(predict(allele$f$fit, new, type = "response"))
  ps <- as.vector(predict(allele$s$fit, new, type = "response"))
  ps / (pf + ps)
}

# the fit and its bioassays ---------------------------------------------------
f <- new.env()
load(fit_file, envir = f)
df <- f$df
mask <- rast("data/clean/raster_mask.tif")
cell_xy <- terra::xyFromCell(mask, f$unique_cells)
cell_share <- predict_share(cell_xy[, 1], cell_xy[, 2])

# the share at every mask cell, as a raster and a map with the records
share_raster <- mask
mask_cells <- which(!is.na(terra::values(mask, mat = FALSE)))
mask_xy <- terra::xyFromCell(mask, mask_cells)
share_values <- rep(NA_real_, terra::ncell(mask))
share_values[mask_cells] <- predict_share(mask_xy[, 1], mask_xy[, 2])
terra::values(share_raster) <- share_values
names(share_raster) <- "kdr_995s_share"
terra::writeRaster(share_raster, "outputs/hinge_screen/kdr_995s_share.tif",
                   overwrite = TRUE)
grid <- as.data.frame(share_raster, xy = TRUE, na.rm = TRUE)
p <- ggplot(grid, aes(x, y, fill = kdr_995s_share)) +
  geom_raster() +
  geom_point(data = records, aes(longitude, latitude, shape = allele),
             inherit.aes = FALSE, size = 0.7, colour = "black") +
  scale_shape_manual(values = c(f = 1, s = 4),
                     labels = c(f = "995F typed", s = "995S typed")) +
  scale_fill_viridis_c(limits = c(0, 1)) +
  coord_equal() +
  labs(x = NULL, y = NULL, fill = "995S share", shape = NULL,
       title = "Share of kdr alleles that are Vgsc-995S (L1014S)",
       subtitle = paste("S / (F + S) from binomial GAMs of each allele's",
                        "frequency (Vector Atlas)")) +
  theme_minimal()
ggsave("figures/kdr_share_map.png", p, width = 8, height = 8, dpi = 130,
       bg = "white")

fold <- list(draws = f$draws, options = f$model_options,
             x_cells_init = f$x_cells_init)
n_total <- sum(vapply(fold$draws, nrow, integer(1)))
par <- dynamical_parameter_draws(fold, f$classes_index, f$types,
                                 draw_index = round(seq(1, n_total,
                                                        length.out = 200)))
logit <- dynamical_logit(par, select(df, cell_id, country_id, type_id,
                                     year_id),
                         df, f$x_cell_years, f$cell_years_index)
rho <- colMeans(par$rho_types)[df$type_id]
n_years <- max(f$cell_years_index$year_id)
x_row <- matrix(NA_integer_, max(f$cell_years_index$cell_id), n_years)
x_row[cbind(f$cell_years_index$cell_id, f$cell_years_index$year_id)] <-
  seq_len(nrow(f$cell_years_index))
nets <- matrix(f$x_cell_years[x_row, "nets"], nrow(x_row))
cum_nets <- t(apply(nets, 1, cumsum))

n <- df$mosquito_number
d <- tibble(class = df$insecticide_class, type = df$insecticide_type,
            country = df$country_name, region = df$region,
            residual = qlogis((df$died + 0.5) / (n + 1)) - colMeans(logit),
            cum_nets = cum_nets[cbind(df$cell_id, df$year_id)],
            share = cell_share[df$cell_id],
            w = 1 / ((1 / (df$died + 0.5) + 1 / (n - df$died + 0.5)) *
                       (1 + (n - 1) * rho))) |>
  filter(class == "Pyrethroids")
rm(logit)

stem <- if (gap > 1) sprintf("%s_gap%d", label, gap) else label
pairs <- read.csv(file.path("outputs/hinge_screen",
                            sprintf("%s_change_pairs.csv", stem))) |>
  filter(insecticide_class == "Pyrethroids")
span <- pairs$year_2 - pairs$year_1
elogit_var <- function(died, n, rho) {
  (1 / (died + 0.5) + 1 / (n - died + 0.5)) * (1 + (n - 1) * rho)
}
pairs <- pairs |>
  mutate(w = span ^ 2 / (elogit_var(died_1, n_1, rho) +
                           elogit_var(died_2, n_2, rho)),
         misfit = (-change_obs + change_pred) / span,
         share = cell_share[cell_id])

coefs <- function(fit, terms) {
  s <- summary(fit)$coefficients[terms, , drop = FALSE]
  tibble(term = terms, estimate = s[, "Estimate"], se = s[, "Std. Error"])
}
level_fit <- lm(residual ~ cum_nets * share + factor(country) + factor(type),
                data = d, weights = w)
change_fit <- lm(misfit ~ nets * share + factor(insecticide_type),
                 data = pairs, weights = w)
results <- bind_rows(
  coefs(level_fit, c("cum_nets", "cum_nets:share")) |>
    mutate(test = "level residual ~ cumulative nets x share"),
  coefs(change_fit, c("nets", "nets:share")) |>
    mutate(test = sprintf("change misfit ~ nets x share (%d-year pairs)",
                          gap)))
# the net slope at low and high share, by tercile of the share
tercile_slopes <- bind_rows(lapply(1:3, function(k) {
  breaks <- quantile(d$share, c(0, 1 / 3, 2 / 3, 1))
  g <- d[findInterval(d$share, breaks, rightmost.closed = TRUE) == k, ]
  h <- pairs[findInterval(pairs$share, breaks, rightmost.closed = TRUE) == k, ]
  lf <- coefs(lm(residual ~ cum_nets + factor(country) + factor(type),
                 data = g, weights = w), "cum_nets")
  cf <- coefs(lm(misfit ~ nets + factor(insecticide_type), data = h,
                 weights = w), "nets")
  tibble(share_tercile = k, share_range = sprintf("%.2f-%.2f", breaks[k],
                                                  breaks[k + 1]),
         level_n = nrow(g), level_slope = lf$estimate, level_se = lf$se,
         change_n = nrow(h), change_slope = cf$estimate, change_se = cf$se)
}))
write.csv(bind_rows(results, tercile_slopes),
          file.path("outputs/hinge_screen", sprintf("%s_kdr_share.csv", label)),
          row.names = FALSE)
options(width = 200)
print(as.data.frame(mutate(results, across(where(is.numeric), ~ signif(.x, 3)))),
      row.names = FALSE)
print(as.data.frame(mutate(tercile_slopes,
                           across(where(is.numeric), ~ signif(.x, 3)))),
      row.names = FALSE)
cat("\nshare at bioassays by region:\n")
print(tapply(d$share, d$region, function(v) round(quantile(v, c(.1, .5, .9)), 2)))
