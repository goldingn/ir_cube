# Bar charts of out-of-sample variance explained, against the share that
# bioassay variability makes unexplainable.
#
# Each bar spans the full 100% of observed variance in held-out mortality. The
# bar is light grey, and the part bioassay sampling makes unexplainable is
# washed out toward white from the top, so what remains grey is the predictable
# part - the quantity the reader should be comparing against. The model's share
# is filled in colour from the bottom.
#
#   light grey, full height: observed variance in held-out mortality
#   washed to white from the top: the share attributable to bioassay sampling.
#     Fully washed over the smallest that share could be, half washed over the
#     95% interval, with a dashed rule at the resulting ceiling.
#   colour from the bottom: what the model explains. Solid to the lower bound
#     of its 95% interval, translucent to the upper, rule at the estimate.
#   grey left uncovered: real variation in resistance available to be explained
#     and not explained.
#
# Reading the interval as a translucent extension of the bar rather than as an
# error bar keeps the quantities on one additive scale, so the grey gap is
# always the shortfall. The wash is drawn over the model bars rather than under
# them, so a model reaching past its ceiling is washed out too - the visual
# signal that it is at the limit of what bioassay data can show.
source("R/packages.R")

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

source("R/fig_variance_bars.R")

pooled <- read.csv("outputs/cv_variance_explained.csv")
# every experiment-insecticide cell is scored in the csv; the "shown" flag
# marks the ones that can carry a reading, on the ceiling and interval-width
# rule set in variance_explained.R
by_insecticide_all <- read.csv("outputs/cv_variance_explained_by_insecticide.csv")
by_insecticide <- by_insecticide_all %>% filter(shown)

experiment_order <- c("spatial interpolation", "spatial extrapolation",
                      "temporal change")

model_axis_labels <- c(
  "insecticide mean"        = "insecticide\nmean",
  "nearest recent survey"   = "nearest\nrecent\nsurvey",
  "nearest surveys, best k" = "nearest\nsurveys,\nbest k",
  "dynamical model"         = "dynamical\nmodel")

# lay the bars out within an experiment, and give the noise block one position
# per experiment so it is drawn once per bar group
lay_out <- function(data, model_levels) {
  models <- data %>%
    filter(kind == "model", quantity %in% model_levels) %>%
    mutate(quantity = factor(quantity, levels = model_levels)) %>%
    arrange(experiment, quantity) %>%
    group_by(experiment) %>%
    mutate(position = row_number()) %>%
    ungroup()
  noise <- data %>%
    filter(kind == "noise") %>%
    select(experiment, stratum, estimate, lower, upper, kind, assays, pixels) %>%
    right_join(models %>% select(experiment, stratum, position),
               by = c("experiment", "stratum")) %>%
    mutate(kind = "noise")
  bind_rows(models, noise)
}

# The bars sit flush with the panel edges so that the strip label lines up
# with the first bar, which leaves the panel spacing as the only separation
# between experiments. It has to be a good multiple of the gap between bars
# within a panel or the grouping does not read, and since that gap is a fixed
# share of a panel whose width depends on how many bars it holds, the spacing
# has to be set per figure: 22 pt over four bars, 36 pt over two.
build_figure <- function(data, model_levels,
                         key_model = "dynamical model",
                         panel_spacing = 22,
                         wrap_key = FALSE) {
  laid <- lay_out(data, model_levels)
  laid$experiment <- factor(laid$experiment, levels = experiment_order)
  breaks <- laid %>% filter(kind == "model") %>%
    distinct(experiment, position, quantity)

  # the key sits in the last panel, against its dynamical model bar
  last_panel <- tail(levels(droplevels(laid$experiment)), 1)
  reference <- laid %>%
    filter(experiment == last_panel, kind == "model",
           quantity == key_model)
  reference_noise <- laid %>%
    filter(experiment == last_panel, kind == "noise") %>%
    slice(1)
  key_x <- max(laid$position[laid$experiment == last_panel]) + 0.6
  key <- region_key(reference$estimate[1], reference_noise$estimate[1], key_x,
                    wrap = wrap_key)
  key <- lapply(key, function(layer) {
    layer$data$experiment <- factor(last_panel, levels = experiment_order)
    layer
  })

  ggplot() +
    bar_layers(laid) +
    key +
    facet_wrap(~ experiment, nrow = 1) +
    scale_fill_manual(values = model_colours, breaks = model_levels) +
    scale_x_continuous(breaks = sort(unique(breaks$position)),
                       labels = function(x) {
                         unname(model_axis_labels[model_levels[x]])
                       },
                       expand = expansion(add = 0)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                       expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(xlim = c(1 - 0.4, length(model_levels) + 0.4),
                    clip = "off") +
    labs(x = NULL, y = "% of variance out of sample") +
    base_theme +
    theme(panel.spacing.x = unit(panel_spacing, "pt"))
}

main <- build_figure(
  pooled,
  c("nearest recent survey", "dynamical model"),
  panel_spacing = 36)

# For the caption. Out-of-sample predictive skill against the bioassay noise
# ceiling, for the dynamical model and the baseline a person would apply
# without one - the single nearest survey of the same insecticide in the most
# recent two years available at prediction time.
#
# Each bar is the whole variance in the held-out quantity across that
# experiment's bioassays - mortality for the two spatial experiments, and,
# once the rebuilt forecasting folds land, change in mortality for the
# temporal one, which is why the axis names neither.
# The washed-out band at the top is the share attributable
# to beta-binomial sampling within the assay, which no model can explain, so
# the dashed line is the ceiling on achievable skill; the coloured block from
# the bottom is the share the model explains, solid to the lower bound of its
# 95% interval and translucent to the upper, with a rule at the point estimate.
# The grey between them is real variation in resistance that is available to be
# explained and is not.
#
# Intervals: a pixel-cluster bootstrap over 2,000 replicates for the models,
# since bioassays cluster hard by pixel and the assay count overstates the
# information a fold carries; the posterior of the bioassay overdispersion for
# the ceiling. Sample sizes are in outputs/cv_variance_explained.csv - spatial
# interpolation 1,045 assays in 94 pixels, spatial extrapolation 8,694 in 905,
# temporal change 1,461 in 309.
#
# The three experiments hold out, respectively: records at pixels with no other
# data, contiguous sub-national regions covering about a third of each of six
# well-sampled countries, and a window of years after the training cut.
ggsave("figures/CV_variance_explained.png", main, width = 8.2, height = 4.2,
       dpi = 300, bg = "white")

all_models <- build_figure(
  pooled,
  c("insecticide mean", "nearest recent survey", "dynamical model",
    "nearest surveys, best k"),
  key_model = "nearest surveys, best k",
  # the key is aligned against the best-k bar, which leaves an unexplained
  # band too short to hold its label on one line
  wrap_key = c(FALSE, TRUE))

# For the caption. As the main figure, with two further reference models. The
# insecticide mean is the per-insecticide mean mortality fitted on the training
# data, which is the no-information baseline: it uses only which insecticide
# was tested. The nearest surveys at best k is the same nearest-neighbour rule
# averaged over whichever number of neighbours minimises its own error on the
# held-out records - hindsight the dynamical model is not given, so it is an
# upper bound on what any such rule could achieve rather than an achievable
# skill. The key is aligned against it for that reason.
ggsave("figures/CV_variance_explained_all_models.png", all_models,
       width = 9.6, height = 4.4, dpi = 300, bg = "white")

# per insecticide, one panel per experiment ----------------------------------

# every panel carries the same insecticides on its axis, so the panels stack in
# register; an insecticide with nothing to show in one experiment leaves a gap
# there. One with nothing to show in any of the three is dropped from the axis
# altogether rather than carried as three empty slots.
#
# grouped by insecticide class, as in fig_ir_maps.R and the other figures that
# arrange(desc(class), insecticide): the four pyrethroids, then the three
# organophosphates, then DDT and Bendiocarb
insecticide_class_order <- c(
  "Alpha-cypermethrin", "Deltamethrin", "Lambda-cyhalothrin", "Permethrin",
  "Fenitrothion", "Malathion", "Pirimiphos-methyl",
  "DDT", "Bendiocarb")
present <- unique(by_insecticide$stratum)
stopifnot(setequal(present, intersect(insecticide_class_order, present)))
insecticide_levels <- intersect(insecticide_class_order, present)

# the key is drawn past the right-hand end of the axis, so it has to describe
# the bar it sits beside. That means the rightmost insecticide needs something
# to show: the key goes in whichever experiment scores it highest, rather than
# always in the first panel.
key_insecticide <- tail(insecticide_levels, 1)
key_experiment <- by_insecticide %>%
  filter(kind == "model", quantity == "dynamical model",
         stratum == key_insecticide, is.finite(estimate)) %>%
  arrange(desc(estimate)) %>%
  slice(1) %>%
  pull(experiment)

panel_for <- function(experiment_label, show_key = FALSE) {

  data <- by_insecticide %>% filter(experiment == experiment_label)
  if (nrow(data) == 0) return(NULL)

  models <- data %>%
    filter(kind == "model",
           quantity %in% c("dynamical model", "nearest recent survey")) %>%
    mutate(quantity = factor(quantity,
             levels = c("nearest recent survey", "dynamical model"))) %>%
    arrange(stratum, quantity) %>%
    group_by(stratum) %>% mutate(offset = row_number()) %>% ungroup() %>%
    mutate(position = match(stratum, insecticide_levels) +
             (offset - 1.5) * 0.38)
  noise <- data %>% filter(kind == "noise") %>%
    select(stratum, estimate, lower, upper, assays, pixels) %>%
    right_join(models %>% select(stratum, position), by = "stratum") %>%
    mutate(kind = "noise")
  laid <- bind_rows(models %>% mutate(kind = "model"), noise)

  reference <- models %>%
    filter(quantity == "dynamical model", stratum == key_insecticide)
  reference_noise <- noise %>%
    filter(stratum == key_insecticide) %>% slice(1)
  # keyed off the full axis, not the last bar present, so the key sits in the
  # same place in every panel
  key_x <- length(insecticide_levels) + 0.45

  ggplot() +
    bar_layers(laid, width = 0.32) +
    (if (show_key) {
      region_key(reference$estimate[1], reference_noise$estimate[1], key_x,
                 label_size = 2.3, wrap = TRUE)
    } else NULL) +
    scale_fill_manual(values = model_colours,
                      breaks = c("nearest recent survey", "dynamical model")) +
    scale_x_continuous(
      breaks = seq_along(insecticide_levels),
      labels = gsub("-", "-\n", insecticide_levels),
      expand = expansion(add = 0)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                       expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(xlim = c(1 - 0.35, length(insecticide_levels) + 0.35),
                    clip = "off") +
    labs(x = NULL, y = "% of variance out of sample",
         subtitle = experiment_label) +
    base_theme +
    theme(axis.text.x = element_text(size = 7.5),
          # matched to the facet strips of the other two figures, and close
          # enough to the panel that it reads as its heading
          plot.subtitle = element_text(face = "bold", size = 11, hjust = 0,
                                       margin = margin(b = 5)),
          # room above for the heading, which would otherwise sit on the tick
          # labels of the panel above it
          plot.margin = margin(16, 40, 2, 6))
}

panels <- lapply(experiment_order, function(label) {
  panel_for(label, show_key = label == key_experiment)
})
panels <- panels[!vapply(panels, is.null, logical(1))]

by_type_figure <- wrap_plots(panels, ncol = 1)

# For the caption. The same comparison broken down by insecticide, one panel
# per experiment, grouped by insecticide class. Only cells where the comparison
# can be read are drawn: at least half the observed spread in held-out mortality
# has to be real variation rather than assay noise, and both models' 95%
# intervals have to be narrower than 100 points. That excludes 12 of the 27
# experiment-insecticide cells, and Fenitrothion and Malathion in every
# experiment - for Fenitrothion because its held-out mortality barely varies
# (standard deviation 2.3 points in the interpolation fold, 83% of assays
# reading exactly 100%), so there is nothing for any model to explain. The
# excluded cells and the reason for each are in
# outputs/cv_variance_explained_by_insecticide.csv.
#
# The ceiling is each insecticide's own, since bioassay
# overdispersion is estimated per type and ranges from 0.095 for
# Lambda-cyhalothrin to 0.250 for Alpha-cypermethrin. Every insecticide with
# held-out records is shown and none is excluded on precision, so the thinnest
# bars carry very wide intervals; assay and pixel counts per bar are in
# outputs/cv_variance_explained_by_insecticide.csv. Where a bar's ceiling is
# low the observed variation in mortality is mostly bioassay noise and there is
# correspondingly little for any model to explain.
ggsave("figures/CV_variance_explained_by_insecticide.png", by_type_figure,
       width = 10, height = 9, dpi = 300, bg = "white", limitsize = FALSE)

cat("\nper-insecticide panels: ", sum(by_insecticide_all$shown) /
      length(unique(by_insecticide_all$quantity)), " cells shown of ",
    nrow(distinct(by_insecticide_all, experiment, stratum)), "\n", sep = "")
print(as.data.frame(by_insecticide_all %>%
  filter(!shown) %>%
  distinct(experiment, stratum, assays, pixels, drop_reason)), row.names = FALSE)

saveRDS(list(main = main, all_models = all_models, by_insecticide = by_type_figure,
             pooled = pooled, by_insecticide_data = by_insecticide_all),
        "outputs/figure_variance_explained.RDS")

cat("written:\n",
    " figures/CV_variance_explained.png\n",
    " figures/CV_variance_explained_all_models.png\n",
    " figures/CV_variance_explained_by_insecticide.png\n", sep = "")
