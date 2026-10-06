# choose the smoothness of the kdr allele frequency surfaces by
# cross-validation (#46)

# load packages and functions
source("R/packages.R")
source("R/functions.R")
source("R/kdr_functions.R")
source("R/kdr_cv_functions.R")

# The kdr surfaces (R/kdr_surface_check.R) are binomial GAMs of the 995F and
# 995S allele frequencies, with a spatial smooth and a linear year trend.
# Records cluster by study, and REML treats them as independent, so it fits
# differences between studies as spatial signal. Here we score smoothness by
# held-out binomial deviance, with two kinds of folds: studies held out
# together, and 5 x 5 degree blocks held out together. The candidates are a
# thin-plate spline with its smoothing parameter fixed on a grid, the same
# spline with REML smoothness selection, and Matern Gaussian process smooths
# (REML) of three ranges.
#
#   Rscript R/kdr_cv.R
#
# Writes outputs/kdr/kdr_cv_scores.csv, figures/kdr_cv_scores.png and
# figures/kdr_cv_<variant>.png.

# the smoothing parameters (log10) of the fixed thin-plate candidates, which
# span about 64 to 3.5 effective degrees of freedom on these data
log10_sps <- seq(-0.5, 4, by = 0.5)

candidates <- c(
  map(log10_sps,
      \(lsp) candidate(tp_smooth, sp = 10 ^ lsp,
                       type = "thin plate, fixed")) %>%
    setNames(sprintf("tp_sp%+.1f", log10_sps)),
  list(tp_reml = candidate(tp_smooth, type = "thin plate, REML")),
  map(c(5, 10, 20),
      \(range) candidate(gp_smooth(range), type = "Matern GP, REML")) %>%
    setNames(sprintf("gp%i", c(5, 10, 20)))
)

# cross-validation -------------------------------------------------------------

records <- kdr_cv_records()
folds <- kdr_cv_folds(records)
held_out <- map(variants, cv_variant, records = records, folds = folds,
                candidates = candidates) %>%
  bind_rows()

# scores -----------------------------------------------------------------------

scores <- score_candidates(held_out, candidates)

dir.create("outputs/kdr", showWarnings = FALSE, recursive = TRUE)
write.csv(scores, "outputs/kdr/kdr_cv_scores.csv", row.names = FALSE)

print(scores %>%
        select(variant, fold_type, candidate, edf, total_difference,
               total_difference_se, per_allele_difference,
               per_allele_difference_se, failed_fits) %>%
        mutate(across(where(is.numeric), \(x) signif(x, 3))) %>%
        as.data.frame())

# the candidate each fold type and score picks
picks <- bind_rows(
  scores %>%
    group_by(variant, fold_type) %>%
    slice_min(total_deviance, n = 1) %>%
    mutate(metric = "total deviance"),
  scores %>%
    group_by(variant, fold_type) %>%
    slice_min(mean_deviance_per_allele, n = 1) %>%
    mutate(metric = "mean deviance per allele")
) %>%
  ungroup() %>%
  select(variant, fold_type, metric, candidate, edf)
print(as.data.frame(picks))

# score against edf
score_plot_data <- bind_rows(
  scores %>%
    transmute(variant, fold_type, candidate, edf, type,
              metric = "total deviance",
              difference = total_difference, se = total_difference_se),
  scores %>%
    transmute(variant, fold_type, candidate, edf, type,
              metric = "mean deviance per allele",
              difference = per_allele_difference,
              se = per_allele_difference_se)
) %>%
  mutate(
    panel = sprintf("%s, %s folds", variant, fold_type),
    metric = factor(metric, c("total deviance", "mean deviance per allele")),
    marker = case_when(
      type == "thin plate, fixed" ~ "thin plate, fixed sp",
      candidate == "tp_reml" ~ "thin plate, REML",
      .default = sprintf("Matern GP, REML, range %s degrees",
                         str_remove(candidate, "gp"))
    ),
    # spread the GP candidates, whose edf are close, a little along the axis
    edf_plot = edf * case_when(candidate == "gp5" ~ 0.97,
                               candidate == "gp20" ~ 1.03,
                               .default = 1)
  )

score_plot <- ggplot(score_plot_data,
                     aes(edf_plot, difference, colour = marker)) +
  geom_line(data = filter(score_plot_data, type == "thin plate, fixed"),
            colour = grey(0.6)) +
  geom_linerange(aes(ymin = difference - se, ymax = difference + se)) +
  geom_point(size = 1.8) +
  facet_grid(metric ~ panel, scales = "free_y") +
  scale_x_log10() +
  scale_colour_manual(
    values = c("thin plate, fixed sp" = grey(0.3),
               "thin plate, REML" = "#c0392b",
               "Matern GP, REML, range 5 degrees" = "#85c1e9",
               "Matern GP, REML, range 10 degrees" = "#2e86c1",
               "Matern GP, REML, range 20 degrees" = "#1b4f72")
  ) +
  labs(x = "effective degrees of freedom of the spatial smooth (full data)",
       y = "held-out score, difference from the best (+/- 1 se)",
       colour = NULL,
       title = "Cross-validated smoothness of the kdr surfaces") +
  guides(colour = guide_legend(nrow = 2)) +
  theme_bw(base_size = 9) +
  theme(legend.position = "top", strip.background = element_blank())

ggsave("figures/kdr_cv_scores.png", score_plot, width = 11, height = 6,
       dpi = 150, bg = "white")

# maps -------------------------------------------------------------------------

# the predictions at 2015 of the candidates picked for each variant
grid <- kdr_prediction_grid()
for (v in variants) {
  panels <- picks %>%
    filter(variant == v) %>%
    group_by(candidate, edf) %>%
    summarise(picked_by = paste(sprintf("%s folds, %s", fold_type, metric),
                                collapse = "\n"),
              .groups = "drop") %>%
    arrange(desc(edf)) %>%
    mutate(panel = sprintf("%s (edf %.0f)\npicked by %s", candidate, edf,
                           picked_by))
  data <- allele_data(records, v)
  predict_candidates(data, panels, candidates, grid) %>%
    map_predictions(data, v, sprintf("figures/kdr_cv_%s.png", v))
}
