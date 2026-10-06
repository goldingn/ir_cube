# A single latent spatial factor behind the kdr variants Vgsc-995F (L1014F)
# and Vgsc-995S (L1014S) (#37): a smooth West-East axis of kdr type, as a
# candidate covariate for regional differences in net selection, and a check
# of whether it explains the fit's net-use misfit.
#
#   Rscript R/kdr_factor.R <fit RData> <label> [gap]
#
# From the Vector Atlas genotypic records (data/raw VA_FULL_DATA_20250721.csv;
# a study may type one variant only), the resistant-allele counts of each
# variant v (frequency x 2 x the number tested) are binomial with
#
#   logit p_v(x, t) = alpha_v + beta_v (t - 2012) + lambda_v u(x)
#
# u(x) = B(x) theta, a thin-plate regression basis of k functions on
# longitude and latitude (mgcv's smoothCon, unpenalised, so k sets the
# smoothness), shared by the variants; lambda_F = 1 and lambda_S is free (u is
# identified up to that scale); the year trends keep the rise of kdr over time
# out of u. Fitted by maximum likelihood (BFGS) for k = 5, 8 and 12, compared
# by AIC. u at every mask cell, scaled to mean 0 and sd 1 over the bioassay
# pixels, is the covariate. Then, for pyrethroids on the fit, as in
# R/kdr_share.R: the level residual on cumulative net use x u, and the change
# misfit of within-pixel pairs on net use x u, and the net slopes by tercile
# of u.
#
# Writes outputs/hinge_screen/kdr_factor_fits.csv, kdr_factor.tif,
# <label>_kdr_factor.csv and figures/kdr_factor_map_k<k>.png.

arguments <- commandArgs(trailingOnly = TRUE)
fit_file <- arguments[1]
label <- arguments[2]
gap <- if (length(arguments) >= 3) as.integer(arguments[3]) else 5L
va_file <- "/home/nick/Dropbox/github/ir_cube/data/raw/VA_FULL_DATA_20250721.csv"
basis_sizes <- c(5, 8, 12)

source("R/packages.R")
source("R/functions.R")
source("R/dynamical_model.R")
source("R/dynamical_predictions.R")
suppressMessages(library(mgcv))

# variant records, one row per record and variant typed -------------------------
va <- read.csv(va_file)
genotypic <- va |>
  filter(insecticide_resistance_data %in% c("genotypic", "both")) |>
  transmute(country, year = year_start,
            latitude = latitude_1, longitude = longitude_1,
            tested = vgsc_number_of_mosquitoes_tested,
            f = as.numeric(vgsc995f_percent) / 100,
            s = as.numeric(vgsc995s_percent) / 100) |>
  filter(!is.na(tested), tested > 0, !is.na(year),
         between(longitude, -26, 64), between(latitude, -36, 38))
variants <- bind_rows(
  genotypic |> filter(!is.na(f)) |> mutate(variant = "F", frequency = f),
  genotypic |> filter(!is.na(s)) |> mutate(variant = "S", frequency = s)) |>
  mutate(n = 2 * tested, resistant = round(frequency * n),
         is_s = as.numeric(variant == "S"),
         year_c = year - 2012)
cat(sprintf("records: 995F %d, 995S %d\n", sum(variants$variant == "F"),
            sum(variants$variant == "S")))

# the factor model --------------------------------------------------------------
# the basis of the spatial factor, from the record locations; predict_basis()
# evaluates it anywhere
make_basis <- function(k) {
  sm <- smoothCon(s(longitude, latitude, k = k, bs = "tp"),
                  data = variants, absorb.cons = TRUE)[[1]]
  list(smooth = sm, X = PredictMat(sm, variants))
}
predict_basis <- function(basis, longitude, latitude) {
  PredictMat(basis$smooth, data.frame(longitude = longitude,
                                      latitude = latitude))
}
# the logit of each record's allele frequency under parameters par:
# alpha_F, alpha_S, beta_F, beta_S, log(-lambda_S)... lambda_S is left free
# in sign, so the fit can say whether the variants load in opposite directions
unpack <- function(par, k) {
  list(alpha = par[1:2], beta = par[3:4], lambda_s = par[5],
       theta = par[5 + seq_len(k - 1)])
}
record_logit <- function(p, X) {
  u <- as.vector(X %*% p$theta)
  lambda <- ifelse(variants$is_s == 1, p$lambda_s, 1)
  p$alpha[variants$is_s + 1] + p$beta[variants$is_s + 1] * variants$year_c +
    lambda * u
}
negative_loglik <- function(par, X, k) {
  eta <- record_logit(unpack(par, k), X)
  -sum(dbinom(variants$resistant, variants$n, plogis(eta), log = TRUE))
}
fits <- lapply(basis_sizes, function(k) {
  basis <- make_basis(k)
  start <- c(qlogis(0.3), qlogis(0.3), 0, 0, -1, rep(0, k - 1))
  o <- optim(start, negative_loglik, X = basis$X, k = k, method = "BFGS",
             control = list(maxit = 2000))
  list(k = k, basis = basis, par = unpack(o$par, k), loglik = -o$value,
       parameters = length(o$par), converged = o$convergence == 0)
})
fit_table <- bind_rows(lapply(fits, function(f) {
  tibble(k = f$k, loglik = f$loglik, parameters = f$parameters,
         aic = -2 * f$loglik + 2 * f$parameters, converged = f$converged,
         lambda_s = f$par$lambda_s, beta_f = f$par$beta[1],
         beta_s = f$par$beta[2])
}))
write.csv(fit_table, "outputs/hinge_screen/kdr_factor_fits.csv",
          row.names = FALSE)
print(as.data.frame(mutate(fit_table, across(where(is.numeric),
                                             ~ signif(.x, 4)))))

# the factor at the bioassays and every mask cell, and maps -------------------
f <- new.env()
load(fit_file, envir = f)
df <- f$df
mask <- rast("data/clean/raster_mask.tif")
mask_cells <- which(!is.na(terra::values(mask, mat = FALSE)))
mask_xy <- terra::xyFromCell(mask, mask_cells)
cell_xy <- terra::xyFromCell(mask, f$unique_cells)
bioassay_cells <- unique(df$cell_id)
factor_at <- function(fit, xy) {
  as.vector(predict_basis(fit$basis, xy[, 1], xy[, 2]) %*% fit$par$theta)
}
for (fit in fits) {
  u_cell <- factor_at(fit, cell_xy)
  centre <- mean(u_cell[df$cell_id])
  spread <- sd(u_cell[df$cell_id])
  values <- rep(NA_real_, terra::ncell(mask))
  values[mask_cells] <- (factor_at(fit, mask_xy) - centre) / spread
  u_raster <- mask
  terra::values(u_raster) <- values
  names(u_raster) <- "kdr_factor"
  if (fit$k == basis_sizes[which.min(fit_table$aic)]) {
    terra::writeRaster(u_raster, "outputs/hinge_screen/kdr_factor.tif",
                       overwrite = TRUE)
  }
  grid <- as.data.frame(u_raster, xy = TRUE, na.rm = TRUE)
  p <- ggplot(grid, aes(x, y, fill = kdr_factor)) +
    geom_raster() +
    geom_point(data = distinct(variants, longitude, latitude, variant),
               aes(longitude, latitude, shape = variant), inherit.aes = FALSE,
               size = 0.7, colour = "black") +
    scale_shape_manual(values = c(F = 1, S = 4),
                       labels = c(F = "995F typed", S = "995S typed")) +
    scale_fill_distiller(palette = "RdBu", limits = c(-3, 3),
                         oob = scales::squish) +
    coord_equal() +
    labs(x = NULL, y = NULL, fill = "factor u\n(sd units)", shape = NULL,
         title = sprintf("Latent kdr factor, %d basis functions (AIC %.0f)",
                         fit$k, -2 * fit$loglik + 2 * fit$parameters),
         subtitle = sprintf(paste("logit p_v = alpha_v + beta_v t + lambda_v u(x);",
                                  "lambda_F = 1, lambda_S = %.2f"),
                            fit$par$lambda_s)) +
    theme_minimal()
  ggsave(sprintf("figures/kdr_factor_map_k%d.png", fit$k), p, width = 8,
         height = 8, dpi = 130, bg = "white")
}

# the misfit tests with the selected factor -----------------------------------
best <- fits[[which.min(fit_table$aic)]]
u_cell <- factor_at(best, cell_xy)
u_cell <- (u_cell - mean(u_cell[df$cell_id])) / sd(u_cell[df$cell_id])

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
            country = df$country_name,
            residual = qlogis((df$died + 0.5) / (n + 1)) - colMeans(logit),
            cum_nets = cum_nets[cbind(df$cell_id, df$year_id)],
            u = u_cell[df$cell_id],
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
         u = u_cell[cell_id])
coefs <- function(fit, terms) {
  s <- summary(fit)$coefficients[terms, , drop = FALSE]
  tibble(term = terms, estimate = s[, "Estimate"], se = s[, "Std. Error"])
}
results <- bind_rows(
  coefs(lm(residual ~ cum_nets * u + factor(country) + factor(type), data = d,
           weights = w), c("cum_nets", "cum_nets:u")) |>
    mutate(test = "level residual ~ cumulative nets x u"),
  coefs(lm(misfit ~ nets * u + factor(insecticide_type), data = pairs,
           weights = w), c("nets", "nets:u")) |>
    mutate(test = sprintf("change misfit ~ nets x u (%d-year pairs)", gap)))
breaks <- quantile(d$u, c(0, 1 / 3, 2 / 3, 1))
terciles <- bind_rows(lapply(1:3, function(k) {
  g <- d[findInterval(d$u, breaks, rightmost.closed = TRUE) == k, ]
  h <- pairs[findInterval(pairs$u, breaks, rightmost.closed = TRUE) == k, ]
  lf <- coefs(lm(residual ~ cum_nets + factor(country) + factor(type),
                 data = g, weights = w), "cum_nets")
  cf <- coefs(lm(misfit ~ nets + factor(insecticide_type), data = h,
                 weights = w), "nets")
  tibble(u_tercile = k, u_range = sprintf("%.2f to %.2f", breaks[k],
                                          breaks[k + 1]),
         level_slope = lf$estimate, level_se = lf$se, change_n = nrow(h),
         change_slope = cf$estimate, change_se = cf$se)
}))
write.csv(bind_rows(results, terciles),
          file.path("outputs/hinge_screen", sprintf("%s_kdr_factor.csv", label)),
          row.names = FALSE)
options(width = 200)
cat(sprintf("\nselected k = %d\n", best$k))
print(as.data.frame(mutate(results, across(where(is.numeric), ~ signif(.x, 3)))),
      row.names = FALSE)
print(as.data.frame(mutate(terciles, across(where(is.numeric), ~ signif(.x, 3)))),
      row.names = FALSE)
