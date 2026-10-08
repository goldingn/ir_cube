# Compare the average selection pressure of each driver on each insecticide
# type across fits (#47), from the draws R/fig_covariate_effects.R saves
# (outputs/species_runs/grids/selection_pressure_draws_<label>.rds): the
# posterior median and 90% interval of the pressure of every covariate
# together (overall), of the covariate groups (vector control, agriculture,
# human population) and of LLIN use and IRS coverage alone, by type, the fits
# side by side; for the fits with the latent smooths (V5) both with the
# smooth at 0 (beta alone) and post hoc (x exp(u_s(x)) per cell); and the
# share of LLIN use in each type's total pressure, which a multiplier common
# to every covariate at a cell leaves unchanged.
#
#   Rscript R/selection_effects_compare.R
#
# Writes figures/species_runs/grids/selection_effects_compare.png,
# figures/species_runs/grids/selection_effects_wb_v5_chains.png and
# outputs/species_runs/grids/selection_effects_compare.csv. Plain R; seconds.

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

input_dir <- "outputs/species_runs/grids"
figure_dir <- "figures/species_runs/grids"

# the fits: label, model and likelihood
fits <- tribble(
  ~label,        ~model,     ~likelihood,
  "ref_f0",      "ref_f0",   "beta-binomial",
  "V3f",         "V3f",      "beta-binomial",
  "V4_class",    "V4_class", "beta-binomial",
  "V5",          "V5",       "beta-binomial",
  "wb_ref",      "ref_f0",   "weighted binomial",
  "wb_v3f",      "V3f",      "weighted binomial",
  "wb_v4_class", "V4_class", "weighted binomial",
  "wb_v5",       "V5",       "weighted binomial",
  "wb_v5_ch1",   "V5 ch 1",  "weighted binomial",
  "wb_v5_ch2",   "V5 ch 2",  "weighted binomial",
  "wb_v5_ch3",   "V5 ch 3",  "weighted binomial",
  "wb_v5_ch4",   "V5 ch 4",  "weighted binomial")

# the insecticide types in class order, and their classes
type_order <- c("Alpha-cypermethrin", "Deltamethrin", "Lambda-cyhalothrin",
                "Permethrin", "DDT", "Bendiocarb", "Fenitrothion",
                "Malathion", "Pirimiphos-methyl")
pyrethroids <- type_order[1:4]

# the drivers: sums over covariates (by their column names) per draw
driver_columns <- function(covariates) {
  list(
    "Overall" = covariates,
    "Vector control" = c("nets", "irs"),
    "LLIN use" = "nets",
    "IRS coverage" = "irs",
    "Agriculture" = grep(":g_ag$", covariates, value = TRUE),
    "Human population" = grep("^pop", covariates, value = TRUE))
}

# per draw, the pressure of each driver on each type (draws x types), and
# the share of LLIN use in the total, for one fit's version
driver_draws <- function(x) {
  columns <- driver_columns(dimnames(x)[[2]])
  out <- lapply(columns, function(cols) {
    apply(x[, cols, , drop = FALSE], c(1, 3), sum)
  })
  out[["LLIN use share"]] <- out[["LLIN use"]] / out[["Overall"]]
  out
}

summaries <- list()
for (i in seq_len(nrow(fits))) {
  draws <- readRDS(file.path(input_dir, sprintf(
    "selection_pressure_draws_%s.rds", fits$label[i])))
  for (version in intersect(c("beta", "posthoc"), names(draws))) {
    drivers <- driver_draws(draws[[version]])
    for (driver in names(drivers)) {
      d <- drivers[[driver]]
      summaries[[length(summaries) + 1]] <- tibble(
        fits[i, ],
        version = version,
        driver = driver,
        insecticide = colnames(d),
        median = apply(d, 2, median),
        q05 = apply(d, 2, quantile, 0.05),
        q95 = apply(d, 2, quantile, 0.95))
    }
  }
}
summaries <- bind_rows(summaries) %>%
  mutate(entry = case_when(
    version == "posthoc" ~ paste0(model, ": post hoc"),
    str_starts(model, "V5") ~ paste0(model, ": beta alone"),
    .default = model))
write.csv(summaries, file.path(input_dir, "selection_effects_compare.csv"),
          row.names = FALSE)


# the comparison figure ---------------------------------------------------------

entries <- c("ref_f0", "V3f", "V4_class", "V5: beta alone", "V5: post hoc")
entry_colours <- c("ref_f0" = "#8c8c8c", "V3f" = "#2a78d6",
                   "V4_class" = "#1baf7a", "V5: beta alone" = "#eb6834",
                   "V5: post hoc" = "#eb6834")
entry_shapes <- c("ref_f0" = 16, "V3f" = 16, "V4_class" = 16,
                  "V5: beta alone" = 1, "V5: post hoc" = 16)
driver_order <- c("Overall", "Vector control", "LLIN use", "IRS coverage",
                  "Agriculture", "Human population")
floor_value <- 1e-4
dodge <- position_dodge(width = 0.8)

pooled <- summaries %>%
  filter(entry %in% entries, !str_detect(model, "ch")) %>%
  mutate(entry = factor(entry, levels = entries),
         insecticide = factor(insecticide, levels = type_order),
         likelihood = factor(likelihood,
                             levels = c("beta-binomial", "weighted binomial")))

# shaded behind the pyrethroids from `bottom`: 0 on a log scale, where it
# maps to -Inf, and -Inf on a linear one
dot_whisker <- function(data, bottom = -Inf) {
  ggplot(data, aes(x = insecticide, y = median, colour = entry,
                   shape = entry)) +
    annotate("rect", xmin = 0.5, xmax = 4.5, ymin = bottom, ymax = Inf,
             fill = grey(0.94)) +
    geom_linerange(aes(ymin = q05, ymax = q95), position = dodge,
                   linewidth = 0.5) +
    geom_point(position = dodge, size = 1.6, stroke = 0.7) +
    scale_colour_manual(values = entry_colours, name = NULL) +
    scale_shape_manual(values = entry_shapes, name = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "top",
          panel.spacing.y = unit(10, "pt"),
          strip.text.y = element_text(angle = 0, hjust = 0)) +
    xlab(NULL)
}

pressure_panel <- pooled %>%
  filter(driver %in% driver_order) %>%
  mutate(driver = factor(driver, levels = driver_order),
         q05 = pmax(q05, floor_value),
         median = pmax(median, floor_value)) %>%
  dot_whisker(bottom = 0) +
  facet_grid(driver ~ likelihood, scales = "free_y") +
  scale_y_log10(labels = function(x) format(x, scientific = FALSE,
                                            drop0trailing = TRUE)) +
  ylab("Average selection pressure (log scale; median, 90% interval)") +
  theme(axis.text.x = element_blank())

share_panel <- pooled %>%
  filter(driver == "LLIN use share") %>%
  mutate(driver = "LLIN use\nshare of\noverall") %>%
  dot_whisker() +
  facet_grid(driver ~ likelihood) +
  scale_y_continuous(limits = c(0, 1)) +
  ylab("Share") +
  theme(legend.position = "none", strip.text.x = element_blank())

compare_plot <- (pressure_panel / share_panel) +
  plot_layout(heights = c(6, 1.3)) +
  plot_annotation(
    title = "Average selection pressure by driver and insecticide type",
    subtitle = paste0(
      "Pyrethroids shaded. V5 hollow: latent smooth at 0 (beta alone); ",
      "V5 solid: x exp(u_s(x)) per cell, averaged post hoc.\n",
      "V5 beta-binomial: chains 1-2; weighted binomial V5: all 4 chains ",
      "(not converged). Values below 0.0001 drawn at 0.0001."))

ggsave(file.path(figure_dir, "selection_effects_compare.png"), compare_plot,
       bg = "white", width = 11, height = 14)


# wb_v5 by chain, against the fits without the smooths --------------------------

chain_entries <- c("V3f", "V4_class", paste("V5 ch", 1:4))
chains <- summaries %>%
  filter(likelihood == "weighted binomial",
         model %in% chain_entries,
         str_starts(model, "V5") | version == "beta",
         driver %in% c("Overall", "LLIN use", "LLIN use share")) %>%
  mutate(entry = factor(model, levels = chain_entries),
         version = factor(if_else(version == "beta", "beta alone",
                                  "post hoc"),
                          levels = c("beta alone", "post hoc")),
         driver = factor(driver, levels = c("Overall", "LLIN use",
                                            "LLIN use share")),
         insecticide = factor(insecticide, levels = type_order))
chain_colours <- c("V3f" = "#2a78d6", "V4_class" = "#1baf7a",
                   "V5 ch 1" = "#eb6834", "V5 ch 2" = "#e34948",
                   "V5 ch 3" = "#4a3aa7", "V5 ch 4" = "#e87ba4")
chains_plot <- ggplot(chains, aes(x = insecticide, y = median,
                                  colour = entry, shape = version)) +
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           fill = grey(0.94)) +
  geom_linerange(aes(ymin = q05, ymax = q95,
                     group = interaction(entry, version)),
                 position = position_dodge(width = 0.85), linewidth = 0.5) +
  geom_point(aes(group = interaction(entry, version)),
             position = position_dodge(width = 0.85), size = 1.5,
             stroke = 0.7) +
  scale_colour_manual(values = chain_colours, name = NULL) +
  scale_shape_manual(values = c("beta alone" = 1, "post hoc" = 16),
                     name = NULL) +
  facet_grid(driver ~ ., scales = "free_y") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "top",
        strip.text.y = element_text(angle = 0, hjust = 0)) +
  xlab(NULL) +
  ylab("Average selection pressure, or share (median, 90% interval)") +
  labs(title = "Weighted binomial V5 by chain, against V3f and V4_class",
       subtitle = "Pyrethroids shaded")
ggsave(file.path(figure_dir, "selection_effects_wb_v5_chains.png"),
       chains_plot, bg = "white", width = 10, height = 8)


# a few numbers -----------------------------------------------------------------

options(width = 160)
for (d in c("LLIN use", "LLIN use share", "Overall", "Human population",
            "IRS coverage", "Agriculture")) {
  cat(sprintf("\n%s, median (90%% interval):\n", d))
  summaries %>%
    filter(driver == d) %>%
    mutate(value = sprintf("%.3f (%.3f-%.3f)", median, q05, q95),
           insecticide = substr(insecticide, 1, 12)) %>%
    select(label, version, insecticide, value) %>%
    pivot_wider(names_from = insecticide, values_from = value) %>%
    as.data.frame() %>%
    print(row.names = FALSE)
}
