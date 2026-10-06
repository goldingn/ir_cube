# a joint model of the Vgsc-995F and 995S allele frequencies, with species
# terms chosen by cross-validation, and the static total kdr covariate (#46)

# load packages and functions
source("R/packages.R")
source("R/functions.R")
source("R/kdr_functions.R")
source("R/kdr_cv_functions.R")
source("R/kdr_joint_functions.R")

# Separate surfaces of 995F and 995S sum to more than 1 over Central Africa.
# Here the alleles are modelled together (R/kdr_joint_functions.R), so total
# kdr, p_F + p_S, is at most 1. We then test by study- and block-fold
# cross-validation whether species terms improve prediction: An. arabiensis
# offsets, with samples identified only to the complex either as their own
# level or as a mixture by the Vector Atlas arabiensis fraction r(x), and an
# An. arabiensis deviation smooth at two heavy penalties.
#
#   Rscript R/kdr_joint.R
#
# Writes outputs/kdr/kdr_joint_cv_scores.csv, outputs/kdr/kdr_joint_terms.csv,
# figures/kdr_joint_vs_separate.png, figures/kdr_joint_species.png and
# data/clean/kdr_total_2015.tif (bands complex, arabiensis, other: total kdr
# at 2015 for the complex as a mixture by r(x), for An. arabiensis, and for
# the rest of the complex).

# data -------------------------------------------------------------------------

records <- kdr_cv_records()
obs <- joint_observations(records) %>%
  mutate(
    r = arabiensis_fraction_at(longitude, latitude),
    block = paste(floor(longitude / block_degrees),
                  floor(latitude / block_degrees))
  )
print(count(obs, type, inconsistent_pair))
print(count(obs, species_type))

# the joint model without species terms, against the separate fits -------------

joint <- fit_joint(obs, joint_candidate())
cat(sprintf("joint model: %d iterations, edf %.1f, Pearson scale %.1f\n",
            joint$iterations, joint$edf, joint$scale))

separate <- list(
  F = fit_candidate(candidate(tp_smooth, sp = joint_sp$F, select = TRUE,
                              type = "thin plate"),
                    allele_data(records, "995F")),
  S = fit_candidate(candidate(tp_smooth, sp = joint_sp$S[1],
                              type = "thin plate"),
                    allele_data(records, "995S"))
)

# predict in chunks, as the thin-plate bases are large; f returns a table, or
# a list of tables
in_chunks <- function(new, f, size = 20000) {
  chunks <- split(seq_len(nrow(new)), ceiling(seq_len(nrow(new)) / size))
  results <- map(chunks, \(rows) f(new[rows, , drop = FALSE]))
  if (is.data.frame(results[[1]])) {
    bind_rows(results)
  } else {
    pmap(unname(results), bind_rows)
  }
}

# the map grid: the mask at a quarter of its resolution, at 2015, with the
# project's map mask and the country of each cell
grid <- kdr_prediction_grid()
grid_points <- as.matrix(grid[, c("longitude", "latitude")])
grid <- grid %>%
  mutate(
    in_mask = !is.na(terra::extract(
      terra::aggregate(rast("data/clean/pfpr_water_mask.tif"), 4,
                       fun = "max", na.rm = TRUE), grid_points)[[1]]),
    country = as.character(terra::extract(
      rast("data/clean/country_raster.tif"), grid_points)$country_name),
    region = kdr_region(country),
    r = arabiensis_fraction_at(longitude, latitude)
  )

joint_grid <- in_chunks(grid, \(chunk) predict_joint(joint, chunk))$g
grid <- grid %>%
  mutate(
    separate_F = as.numeric(predict(separate$F, grid, type = "response")),
    separate_S = as.numeric(predict(separate$S, grid, type = "response")),
    joint_F = joint_grid$p_F,
    joint_S = joint_grid$p_S,
    joint_total = joint_grid$total,
    joint_total_se = joint_grid$total_logit_se
  )

# where they differ: at the records typing one variant or both, and over the
# mask, where the separate fits sum to more than 1
joint_obs <- predict_observations(joint, obs)
record_comparison <- obs %>%
  mutate(
    joint_F = joint_obs$F,
    joint_S = joint_obs$S,
    separate_F = as.numeric(predict(separate$F, obs, type = "response")),
    separate_S = as.numeric(predict(separate$S, obs, type = "response"))
  ) %>%
  group_by(type) %>%
  summarise(
    median_difference_F = median(abs(joint_F - separate_F)),
    max_difference_F = max(abs(joint_F - separate_F)),
    median_difference_S = median(abs(joint_S - separate_S)),
    max_difference_S = max(abs(joint_S - separate_S)),
    .groups = "drop"
  )
print(as.data.frame(record_comparison))
cat(sprintf(paste("cells in the mask where the separate fits sum to more",
                  "than 1: %.1f%% (max sum %.2f)\n"),
            100 * mean((grid$separate_F + grid$separate_S)[grid$in_mask] > 1),
            max((grid$separate_F + grid$separate_S)[grid$in_mask])))

# species terms: cross-validation ---------------------------------------------

species_candidates <- list(
  none = joint_candidate("none"),
  own = joint_candidate("own"),
  mixture = joint_candidate("mixture"),
  mixture_deviation_10 = joint_candidate("mixture", deviation = 10),
  mixture_deviation_100 = joint_candidate("mixture", deviation = 100)
)

# the study and block folds, made as in R/kdr_cv_functions.R but of the joint
# observations (so a sample typing both variants is held out whole)
set.seed(fold_seed)
joint_folds <- list(study = group_folds(obs$study),
                    block = group_folds(obs$block))

# held-out deviance of every observation by every candidate, for one fold
joint_predict_fold <- function(fold, folds) {
  train <- obs[folds != fold, ]
  test <- obs[folds == fold, ]
  imap(species_candidates, \(candidate, name) {
    fit <- fit_joint(train, candidate)
    tibble(row = which(folds == fold), fold = fold, candidate = name,
           converged = fit$converged,
           deviance = joint_deviance(test, predict_observations(fit, test)))
  }) %>%
    bind_rows()
}

held_out <- imap(joint_folds, \(folds, fold_type) {
  mclapply(seq_len(n_folds), joint_predict_fold, folds = folds,
           mc.cores = cores) %>%
    bind_rows() %>%
    mutate(fold_type = fold_type)
}) %>%
  bind_rows() %>%
  mutate(m = obs$m[row])

cv_scores <- held_out %>%
  group_by(fold_type, candidate, fold) %>%
  summarise(total = sum(deviance), per_allele = sum(deviance / m),
            records = n(), converged = all(converged), .groups = "drop") %>%
  group_by(fold_type, fold) %>%
  mutate(
    total_fold_difference = total - total[candidate == "none"],
    per_allele_fold_difference = per_allele -
      per_allele[candidate == "none"]
  ) %>%
  group_by(fold_type, candidate) %>%
  summarise(
    total_deviance = sum(total),
    total_difference = sum(total_fold_difference),
    total_difference_se = sd(total_fold_difference) * sqrt(n()),
    mean_deviance_per_allele = sum(per_allele) / sum(records),
    per_allele_difference = sum(per_allele_fold_difference) / sum(records),
    per_allele_difference_se = sd(per_allele_fold_difference) * sqrt(n()) /
      sum(records),
    all_converged = all(converged),
    .groups = "drop"
  ) %>%
  mutate(candidate = factor(candidate, names(species_candidates))) %>%
  arrange(fold_type, candidate)

dir.create("outputs/kdr", showWarnings = FALSE, recursive = TRUE)
write.csv(cv_scores, "outputs/kdr/kdr_joint_cv_scores.csv", row.names = FALSE)
print(cv_scores %>%
        select(fold_type, candidate, total_difference, total_difference_se,
               per_allele_difference, per_allele_difference_se,
               all_converged) %>%
        mutate(across(where(is.numeric), \(x) signif(x, 3))) %>%
        as.data.frame())

# the candidate with the lowest mean deviance per allele, averaged over the
# fold types
chosen_name <- cv_scores %>%
  group_by(candidate) %>%
  summarise(score = mean(mean_deviance_per_allele)) %>%
  slice_min(score, n = 1) %>%
  pull(candidate) %>%
  as.character()
cat("chosen:", chosen_name, "\n")

# species terms: the fitted model ----------------------------------------------

# the offsets, from the mixture model and the chosen one, with standard errors
# scaled by the Pearson overdispersion
fits <- list(own = fit_joint(obs, species_candidates$own),
             mixture = fit_joint(obs, species_candidates$mixture))
fits[[chosen_name]] <- fits[[chosen_name]] %||%
  fit_joint(obs, species_candidates[[chosen_name]])
chosen <- fits[[chosen_name]]

terms <- imap(fits, \(fit, name) {
  keep <- grepl(":(arabiensis|complex|year)$", names(fit$theta))
  tibble(model = name,
         term = names(fit$theta)[keep],
         estimate = fit$theta[keep],
         se = sqrt(diag(fit$scaled_covariance))[keep],
         scale = fit$scale,
         edf = fit$edf)
}) %>%
  bind_rows()
write.csv(terms, "outputs/kdr/kdr_joint_terms.csv", row.names = FALSE)
print(terms %>%
        mutate(across(where(is.numeric), \(x) round(x, 3))) %>%
        as.data.frame())

# predictions of the chosen model on the map grid: the rest of the complex,
# An. arabiensis, and the complex as a mixture by r
chosen_grid <- in_chunks(grid, \(chunk) {
  predict_joint(chosen, chunk, r = chunk$r)
})
grid <- grid %>%
  mutate(
    total_other = chosen_grid$g$total,
    total_arabiensis = chosen_grid$a$total,
    total_complex = chosen_grid$complex$total,
    total_complex_se = chosen_grid$complex$total_logit_se,
    # the complex with the same species composition everywhere: the mean r
    # over the mask
    total_complex_fixed_r = (1 - mean(r[in_mask])) * total_other +
      mean(r[in_mask]) * total_arabiensis
  )

# How much of the West-East contrast in total kdr is species: the difference
# in mean total kdr over the mask between West and East Africa, for the
# complex, for each species group, and for the complex with the same
# composition everywhere.
contrast <- grid %>%
  filter(in_mask, region %in% c("West", "East")) %>%
  group_by(region) %>%
  summarise(across(c(total_complex, total_complex_fixed_r, total_other,
                     total_arabiensis, r), mean)) %>%
  pivot_longer(-region) %>%
  pivot_wider(names_from = region, values_from = value) %>%
  mutate(west_minus_east = West - East)
print(as.data.frame(contrast))

# maps -------------------------------------------------------------------------

borders <- readRDS("data/clean/country_borders.RDS")
map_panel <- function(column, title, limits = c(0, 1), option = "viridis",
                      name = NULL, oob = scales::squish, na.value = "grey50") {
  ggplot(filter(grid, in_mask), aes(longitude, latitude)) +
    geom_tile(aes(fill = .data[[column]]), width = 1 / 6, height = 1 / 6) +
    geom_sf(data = borders, inherit.aes = FALSE, fill = NA,
            colour = "white", linewidth = 0.1) +
    scale_fill_viridis_c(option = option, limits = limits, oob = oob,
                         na.value = na.value, name = name) +
    coord_sf(xlim = c(-18, 52), ylim = c(-35, 25)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text = element_blank())
}
difference_panel <- function(column, title) {
  ggplot(filter(grid, in_mask), aes(longitude, latitude)) +
    geom_tile(aes(fill = .data[[column]]), width = 1 / 6, height = 1 / 6) +
    geom_sf(data = borders, inherit.aes = FALSE, fill = NA,
            colour = "grey40", linewidth = 0.1) +
    scale_fill_distiller(palette = "RdBu", limits = c(-0.3, 0.3),
                         oob = scales::squish, direction = 1, name = NULL) +
    coord_sf(xlim = c(-18, 52), ylim = c(-35, 25)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(), axis.text = element_blank())
}

grid <- grid %>%
  mutate(
    difference_F = joint_F - separate_F,
    difference_S = joint_S - separate_S,
    separate_total = separate_F + separate_S
  )
comparison <- patchwork::wrap_plots(
  map_panel("separate_F", "995F, separate fit"),
  map_panel("joint_F", "995F, joint model"),
  difference_panel("difference_F", "995F, joint - separate"),
  map_panel("separate_S", "995S, separate fit"),
  map_panel("joint_S", "995S, joint model"),
  difference_panel("difference_S", "995S, joint - separate"),
  map_panel("separate_total", "995F + 995S, separate fits\n(over 1 in red)",
            oob = scales::censor, na.value = "#d62728"),
  map_panel("joint_total", "total kdr, joint model"),
  map_panel("joint_total_se", "total kdr, joint model:\nse (logit)",
            limits = c(0, 2), option = "magma"),
  ncol = 3
) +
  patchwork::plot_annotation(
    title = "Joint model of 995F and 995S against separate fits, 2015")
ggsave("figures/kdr_joint_vs_separate.png", comparison, width = 12,
       height = 10.5, dpi = 110, bg = "white")

species_maps <- patchwork::wrap_plots(
  map_panel("total_other",
            "total kdr, An. gambiae s.s.,\nAn. coluzzii and the rest"),
  map_panel("total_arabiensis", "total kdr, An. arabiensis"),
  map_panel("r", "An. arabiensis fraction r(x)\n(Vector Atlas)"),
  map_panel("total_complex", "total kdr, the complex\n(mixture by r)"),
  map_panel("total_complex_se", "total kdr, the complex:\nse (logit)",
            limits = c(0, 2), option = "magma"),
  ncol = 3
) +
  patchwork::plot_annotation(
    title = sprintf("Total kdr (995F + 995S) at 2015, joint model (%s)",
                    chosen_name))
ggsave("figures/kdr_joint_species.png", species_maps, width = 12, height = 7,
       dpi = 110, bg = "white")

# the covariate ----------------------------------------------------------------

# total kdr at 2015 on the model grid (2.5 arcmin), for every cell of the
# model mask
mask <- rast("data/clean/raster_mask.tif")
cells <- which(!is.na(values(mask)))
model_grid <- as.data.frame(xyFromCell(mask, cells)) %>%
  transmute(longitude = x, latitude = y, year_c = 0) %>%
  mutate(r = arabiensis_fraction_at(longitude, latitude))

# total kdr at some locations, for the complex as a mixture by r, for An.
# arabiensis and for the rest of the complex, without standard errors
predict_totals <- function(fit, new) {
  new$species_type <- "other"
  designs <- map(fit$smooths, \(smooth) {
    joint_design(new, smooth, fit$candidate)
  })
  theta <- fit$split_theta(fit$theta)
  total <- function(weight_a) {
    p <- joint_frequencies(theta, designs, weight_a, derivatives = FALSE)$p
    p$F + p$S
  }
  tibble(complex = total(new$r),
         arabiensis = total(rep(1, nrow(new))),
         other = total(rep(0, nrow(new))))
}

totals <- in_chunks(model_grid, \(chunk) predict_totals(chosen, chunk))
covariate <- rast(mask, nlyrs = 3)
names(covariate) <- c("complex", "arabiensis", "other")
band_values <- matrix(NA_real_, ncell(mask), 3)
band_values[cells, ] <- as.matrix(totals[, names(covariate)])
values(covariate) <- band_values
writeRaster(covariate, "data/clean/kdr_total_2015.tif", overwrite = TRUE)
print(global(covariate, c("min", "max", "mean"), na.rm = TRUE))
