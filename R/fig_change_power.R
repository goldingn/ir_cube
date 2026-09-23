# How much data is there to validate a forecast of *change* in resistance, and
# when?
#
# The forecasting experiment asks whether the model predicts the change in
# mortality at a site, not its level. That can only be measured where the same
# pixel and insecticide were assayed in both windows, and at 2.2 and 1.8 assays
# per group the 2020-2022 test turned out to be far too thin to say anything
# about any single site. This finds where in the record there is enough paired
# data to make the test informative, so the forecast origin can be chosen rather
# than inherited from the end of the covariates (#12 review 5.1).
#
# Two views. The first is the literal question: between two years separated by a
# given gap, how much did mortality change at sites assayed in both, and how
# many such sites are there. The second slides a three-year before window and a
# three-year holdout window across the record, which is the design the
# experiment actually uses.

source("R/validation_functions.R")
source("R/validation_folds.R")

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(2026 - 9 - 23)
n_bootstrap <- 2000

rho_table <- read.csv("outputs/bioassay_rho.csv")
rho_for_class <- function(insecticide_class) {
  index <- match(insecticide_class, rho_table$insecticide_class)
  pooled <- rho_table$rho[rho_table$insecticide_class == "all"]
  ifelse(is.na(index), pooled, rho_table$rho[index])
}

# pooled mortality per (pixel, insecticide) in a set of years
window_totals <- function(years) {
  df %>%
    filter(year_start %in% years) %>%
    group_by(group = paste(cell, insecticide_type),
             insecticide_class) %>%
    summarise(assays = n(),
              died = sum(died),
              tested = sum(mosquito_number),
              sizes = list(mosquito_number),
              counts = list(died),
              .groups = "drop")
}

# the change between two windows, at pixels assayed in both
paired_change <- function(before_years, after_years) {

  b <- window_totals(before_years)
  a <- window_totals(after_years)
  shared <- intersect(b$group, a$group)
  if (length(shared) < 3) return(NULL)

  bi <- match(shared, b$group)
  ai <- match(shared, a$group)
  rho <- rho_for_class(a$insecticide_class[ai])

  delta <- a$died[ai] / a$tested[ai] - b$died[bi] / b$tested[bi]
  # irreducible variance of that difference: the two windows are independent
  # given the fractions, so their variances add
  floor_variance <- vapply(seq_along(shared), function(k) {
    noise_floor_var_pooled(a$counts[[ai[k]]], a$sizes[[ai[k]]], rho[k]) +
      noise_floor_var_pooled(b$counts[[bi[k]]], b$sizes[[bi[k]]], rho[k])
  }, numeric(1))

  usable <- is.finite(floor_variance)
  boot <- replicate(n_bootstrap,
                    mean(sample(delta, length(delta), replace = TRUE)))

  data.frame(
    groups = length(shared),
    assays = sum(a$assays[ai]) + sum(b$assays[bi]),
    assays_per_group = (sum(a$assays[ai]) + sum(b$assays[bi])) / length(shared),
    mean_change = mean(delta),
    lower = unname(quantile(boot, 0.025)),
    upper = unname(quantile(boot, 0.975)),
    # the smallest mean change this many groups could distinguish from zero,
    # given bioassay noise alone
    detectable = 1.96 * sqrt(mean(floor_variance[usable]) / sum(usable))
  )

}


# view one: single years, several gaps ------------------------------------

gaps <- c(1, 2, 3, 5)
later_years <- 2002:2024

single_year <- bind_rows(lapply(gaps, function(gap) {
  bind_rows(lapply(later_years, function(year) {
    out <- paired_change(year - gap, year)
    if (is.null(out)) return(NULL)
    out %>% mutate(gap = gap, year = year, .before = everything())
  }))
}))

write.csv(single_year, "outputs/change_power_years.csv", row.names = FALSE)


# view two: the main figure — five-year windows, sliding the cut -----------

# Five years rather than three. The gap between window midpoints is the window
# length, so a five-year window carries a signal roughly 5/3 larger, and it
# spans the 2018-2020 pause together with the decline either side of it rather
# than sitting inside the pause.
window <- 5
cut_years <- 2005:2019

sliding <- bind_rows(lapply(cut_years, function(cut) {
  out <- paired_change(seq(cut - window, cut - 1), seq(cut, cut + window - 1))
  if (is.null(out)) return(NULL)
  out %>% mutate(cut = cut, .before = everything())
}))

# the latest origin whose holdout window is covered by the covariate layers
latest_feasible <- final_year - window + 1
chosen_cuts <- c(2014, 2018)

write.csv(sliding %>% mutate(window = window,
                             covariates_cover_holdout = cut <= latest_feasible),
          "outputs/change_power_windows.csv", row.names = FALSE)

# figures ------------------------------------------------------------------

# the gap is ordinal, so a single hue light to dark, not four unrelated colours
gap_colours <- c("1" = "#9ECAE1", "2" = "#6BAED6", "3" = "#3182BD",
                 "5" = "#08519C")
base <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom")

# every figure is written twice: a 300 dpi raster to look at, and a vector PDF
# for the manuscript. The ggplot objects are also saved, so a panel can be
# restyled for the paper without recomputing anything
save_figure <- function(plot, name, width, height) {
  ggsave(file.path("figures", paste0(name, ".png")), plot, bg = "white",
         width = width, height = height, dpi = 300)
  ggsave(file.path("figures", paste0(name, ".pdf")), plot, bg = "white",
         width = width, height = height, device = cairo_pdf)
  invisible(plot)
}

single_year <- single_year %>% mutate(gap = factor(gap, levels = gaps))

# faceted by gap rather than overlaid: four ribbons on one set of axes is
# unreadable, and the comparison of interest is within a gap over time
gap_labels <- function(x) paste(x, "years apart")

change_panel <- single_year %>%
  ggplot(aes(x = year, y = mean_change, colour = gap, fill = gap)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = grey(0.45)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, colour = NA) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ gap, nrow = 1, labeller = labeller(gap = gap_labels)) +
  scale_colour_manual(values = gap_colours, guide = "none") +
  scale_fill_manual(values = gap_colours, guide = "none") +
  scale_x_continuous(breaks = seq(2004, 2024, by = 5)) +
  coord_cartesian(ylim = c(-0.3, 0.15)) +
  labs(x = NULL, y = "change in mortality",
       title = "Change in bioassay mortality at the same pixel and insecticide",
       subtitle = paste("mean over pixels assayed in both years, with a 95%",
                        "bootstrap interval; negative is rising resistance")) +
  base

count_panel <- single_year %>%
  ggplot(aes(x = year, y = groups, colour = gap)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  facet_wrap(~ gap, nrow = 1, labeller = labeller(gap = gap_labels)) +
  scale_colour_manual(values = gap_colours, guide = "none") +
  scale_x_continuous(breaks = seq(2004, 2024, by = 5)) +
  scale_y_continuous(trans = "log10") +
  labs(x = NULL, y = "pixel-insecticide pairs",
       subtitle = "how many pixels carry an assay in both years (log scale)") +
  base

slope_per_year <- 0.025
power_panel <- single_year %>%
  mutate(expected = slope_per_year * as.numeric(as.character(gap))) %>%
  ggplot(aes(x = year, colour = gap)) +
  geom_ribbon(aes(ymin = detectable, ymax = pmax(expected, detectable),
                  fill = gap), alpha = 0.2, colour = NA) +
  geom_line(aes(y = detectable), linewidth = 0.7) +
  geom_line(aes(y = expected), linewidth = 0.5, linetype = 2) +
  facet_wrap(~ gap, nrow = 1, labeller = labeller(gap = gap_labels)) +
  scale_colour_manual(values = gap_colours, guide = "none") +
  scale_fill_manual(values = gap_colours, guide = "none") +
  scale_x_continuous(breaks = seq(2004, 2024, by = 5)) +
  coord_cartesian(ylim = c(0, 0.16)) +
  labs(x = "later year of the pair", y = "change in mortality",
       subtitle = paste("solid: smallest mean change distinguishable from",
                        "zero given assay noise. dashed: the change expected",
                        "from\nthe historical decline of 2.5 points a year.",
                        "the test has power where the dashed line is above",
                        "the solid one")) +
  base

window_size_figure <- change_panel / count_panel / power_panel +
  plot_layout(heights = c(1, 0.8, 0.9))

save_figure(window_size_figure, "CV_change_window_size", width = 12, height = 10)
saveRDS(list(combined = window_size_figure, data = single_year),
        "outputs/figure_change_window_size.RDS")

# The main figure. Where in the record can a forecast of change be validated,
# and how much data is there to do it with.
chosen <- sliding %>% filter(cut %in% chosen_cuts)

sliding_panel <- sliding %>%
  ggplot(aes(x = cut)) +
  annotate("rect", xmin = latest_feasible + 0.5, xmax = max(cut_years) + 0.5,
           ymin = -Inf, ymax = Inf, fill = grey(0.92)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = grey(0.45)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.22, fill = "#3182BD") +
  geom_line(aes(y = detectable), linewidth = 0.5, colour = "#B2182B") +
  geom_line(aes(y = -detectable), linewidth = 0.5, colour = "#B2182B") +
  geom_line(aes(y = mean_change), linewidth = 0.9, colour = "#08519C") +
  geom_point(data = chosen, aes(y = mean_change), size = 3,
             colour = "#08519C") +
  geom_point(data = chosen, aes(y = mean_change), size = 1.4, colour = "white") +
  scale_x_continuous(breaks = cut_years[cut_years %% 2 == 1],
                     expand = expansion(add = 0.6)) +
  labs(x = NULL, y = "change in mortality",
       title = sprintf("Change in bioassay mortality over a %i-year forecast horizon", window),
       subtitle = paste0("blue: mean change at pixels assayed in both the ",
                         window, "-year holdout window and the ", window,
                         " years before it,\nwith a 95% bootstrap interval. ",
                         "red: what bioassay noise alone allows to be resolved.",
                         "\nringed points are the two origins chosen for ",
                         "validation; shading marks origins whose holdout ",
                         "window runs past the covariate layers")) +
  base

sliding_counts <- sliding %>%
  ggplot(aes(x = cut, y = groups)) +
  annotate("rect", xmin = latest_feasible + 0.5, xmax = max(cut_years) + 0.5,
           ymin = -Inf, ymax = Inf, fill = grey(0.92)) +
  geom_col(fill = "#6BAED6", width = 0.72) +
  geom_col(data = chosen, fill = "#08519C", width = 0.72) +
  geom_text(aes(label = groups), vjust = -0.45, size = 2.9,
            colour = grey(0.35)) +
  scale_x_continuous(breaks = cut_years[cut_years %% 2 == 1],
                     expand = expansion(add = 0.6)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "first year of the holdout window",
       y = "pixel-insecticide pairs",
       subtitle = "pixels carrying an assay in both windows") +
  base

change_power_figure <- sliding_panel / sliding_counts +
  plot_layout(heights = c(1.25, 0.75))

save_figure(change_power_figure, "CV_change_power", width = 9, height = 7.5)
saveRDS(list(combined = change_power_figure,
             change = sliding_panel,
             counts = sliding_counts,
             data = sliding),
        "outputs/figure_change_power.RDS")

cat(sprintf("\n%i-year windows, sliding the cut:\n", window))
print(as.data.frame(sliding %>%
        transmute(cut, groups, assays,
                  per_group = round(assays_per_group, 1),
                  mean_change = round(mean_change, 3),
                  ci = sprintf("[%+.3f, %+.3f]", lower, upper),
                  detectable = round(detectable, 3),
                  ratio = round(abs(mean_change) / detectable, 1),
                  covariates = cut <= latest_feasible)))
