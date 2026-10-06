# smooths for the kdr surfaces that extrapolate less: cross-validation of a
# thin plate with its null space penalised, and of Matern Gaussian processes
# with fixed smoothing (#46)

# load packages and functions
source("R/packages.R")
source("R/functions.R")
source("R/kdr_functions.R")
source("R/kdr_cv_functions.R")

# The thin plates R/kdr_cv.R picks (smoothing parameter 10^1 for 995F, 10^2
# for 995S) extrapolate along their linear null space: far from the records
# the surface follows a plane fitted across Africa. Here we test, on the same
# folds and scores:
#   (a) the thin plate with its null space penalised too (mgcv's `select`),
#       with the null-space smoothing parameter on a grid at the chosen
#       range-space one, so the plane can shrink away, and both on a coarse
#       2-D grid around them;
#   (b) Matern (kappa = 2.5) Gaussian processes without a linear trend, which
#       revert to the mean far from data, with ranges of 10, 20 and 40 degrees
#       and the smoothing parameter on a grid.
# Scores are differences from the chosen thin plate, with standard errors from
# their spread across folds.
#
#   Rscript R/kdr_cv_extrap.R
#
# Writes outputs/kdr/kdr_cv_extrap_scores.csv,
# outputs/kdr/kdr_cv_extrap_zones.csv and figures/kdr_cv_extrap_<variant>.png.

# the log10 smoothing parameters of the thin plates R/kdr_cv.R picks
chosen_log10_sp <- c(`995F` = 1, `995S` = 2)

# (a) log10 null-space smoothing parameters at the chosen range-space one, and
# the coarse 2-D grid around it
null_log10_sps <- seq(-2, 10, by = 2)
range_offsets <- c(-1, -0.5, 0.5, 1)
null_coarse <- c(0, 5, 10)

# (b) the GP ranges (degrees) and log10 smoothing parameters
gp_ranges <- c(10, 20, 40)
gp_log10_sps <- seq(-4, 3, by = 0.5)

# the candidates for a variant, named relative to its chosen thin plate
extrap_candidates <- function(v) {
  range_sp <- chosen_log10_sp[[v]]
  select_candidate <- function(range_offset, null_lsp) {
    candidate(tp_smooth, sp = 10 ^ c(range_sp + range_offset, null_lsp),
              select = TRUE, type = "thin plate, null space penalised")
  }
  coarse <- expand_grid(range_offset = range_offsets, null_lsp = null_coarse)
  gp <- expand_grid(range = gp_ranges, lsp = gp_log10_sps)
  c(
    list(tp_chosen = candidate(tp_smooth, sp = 10 ^ range_sp,
                               type = "thin plate, chosen")),
    map(null_log10_sps, \(null_lsp) select_candidate(0, null_lsp)) %>%
      setNames(sprintf("tp_null%+.0f", null_log10_sps)),
    pmap(coarse, select_candidate) %>%
      setNames(sprintf("tp_range%+.1f_null%+.0f", coarse$range_offset,
                       coarse$null_lsp)),
    pmap(gp, \(range, lsp) {
      candidate(gp_smooth(range), sp = 10 ^ lsp,
                type = sprintf("Matern GP, range %i", range))
    }) %>%
      setNames(sprintf("gp%i_sp%+.1f", gp$range, gp$lsp))
  )
}

# cross-validation -------------------------------------------------------------

records <- kdr_cv_records()
folds <- kdr_cv_folds(records)
held_out <- map(variants, \(v) {
  cv_variant(v, records, folds, extrap_candidates(v))
}) %>%
  bind_rows()

scores <- score_candidates(held_out, extrap_candidates("995F"),
                           reference = c(`995F` = "tp_chosen",
                                         `995S` = "tp_chosen")) %>%
  mutate(group = case_when(
    candidate == "tp_chosen" ~ "chosen",
    grepl("^tp", candidate) ~ "a",
    .default = "b"
  ))

dir.create("outputs/kdr", showWarnings = FALSE, recursive = TRUE)
write.csv(scores, "outputs/kdr/kdr_cv_extrap_scores.csv", row.names = FALSE)

# the best of (a) and of (b) for each variant: the lowest mean deviance per
# allele, averaged over the two fold types
best <- scores %>%
  group_by(variant, group, candidate, edf) %>%
  summarise(per_allele = mean(mean_deviance_per_allele), .groups = "drop") %>%
  group_by(variant, group) %>%
  slice_min(per_allele, n = 1) %>%
  ungroup()
print(as.data.frame(best))

# the scores of the best of each group, and of a few others, against the
# chosen thin plate
print(scores %>%
        group_by(variant, group) %>%
        filter(candidate %in% best$candidate |
                 candidate %in% c("tp_null+10", "tp_null+6") |
                 (group == "b" &
                    per_allele_difference == min(per_allele_difference))) %>%
        ungroup() %>%
        select(variant, fold_type, candidate, edf, total_difference,
               total_difference_se, per_allele_difference,
               per_allele_difference_se, failed_fits) %>%
        mutate(across(where(is.numeric), \(x) signif(x, 3))) %>%
        arrange(variant, candidate, fold_type) %>%
        as.data.frame())

# the best GP of each range
print(scores %>%
        filter(group == "b") %>%
        mutate(range = str_extract(candidate, "^gp[0-9]+")) %>%
        group_by(variant, fold_type, range) %>%
        slice_min(per_allele_difference, n = 1) %>%
        ungroup() %>%
        select(variant, fold_type, candidate, edf, per_allele_difference,
               per_allele_difference_se, total_difference,
               total_difference_se) %>%
        mutate(across(where(is.numeric), \(x) signif(x, 3))) %>%
        as.data.frame())

# maps -------------------------------------------------------------------------

# the 2015 predictions of the chosen thin plate and the best of (a) and (b),
# and their behaviour north of 18 N and south of 15 S
grid <- kdr_prediction_grid()
zones <- list()
for (v in variants) {
  candidates <- extrap_candidates(v)
  panels <- best %>%
    filter(variant == v) %>%
    mutate(
      group = factor(group, c("chosen", "a", "b")),
      panel = sprintf("%s: %s (edf %.0f)",
                      recode(as.character(group),
                             chosen = "chosen thin plate",
                             a = "(a) null space penalised",
                             b = "(b) Matern GP"),
                      candidate, edf)
    ) %>%
    arrange(group) %>%
    select(candidate, panel)
  data <- allele_data(records, v)
  maps <- predict_candidates(data, panels, candidates, grid)
  map_predictions(maps, data, v, sprintf("figures/kdr_cv_extrap_%s.png", v))
  zones[[v]] <- maps %>%
    mutate(
      variant = v,
      zone = case_when(latitude > 18 ~ "north of 18N",
                       latitude < -15 ~ "south of 15S",
                       .default = "between")
    ) %>%
    group_by(variant, candidate, zone) %>%
    summarise(
      mean_min = min(mean),
      mean_median = median(mean),
      mean_max = max(mean),
      se_median = median(se),
      se_max = max(se),
      .groups = "drop"
    )
}
zones <- bind_rows(zones)
write.csv(zones, "outputs/kdr/kdr_cv_extrap_zones.csv", row.names = FALSE)
print(zones %>%
        mutate(across(where(is.numeric), \(x) round(x, 2))) %>%
        as.data.frame())

# the records at the edges, and the mean frequency of the records overall
print(map(variants, \(v) {
  data <- allele_data(records, v)
  tibble(variant = v,
         north_of_18N = sum(data$latitude > 18),
         south_of_15S = sum(data$latitude < -15),
         pooled_frequency = sum(data$successes) / sum(data$trials))
}) %>%
  bind_rows())
