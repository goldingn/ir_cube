# make a figure showing the fitted relationships between covariates and
# selection for each insecticide
#
#   Rscript R/fig_covariate_effects.R [<fitted_model.RData> <label>]
#
# Without arguments, the fit temporary/fitted_model.RData, writing
# figures/covariate_loadings.png and figures/covariate_loadings_effect_size.png.
# With a fit and a label (#47), writes
#   figures/species_runs/grids/covariate_loadings_<label>.png
#   figures/species_runs/grids/covariate_loadings_effect_size_<label>.png
#   outputs/species_runs/grids/selection_pressure_<label>.csv
#   outputs/species_runs/grids/selection_pressure_draws_<label>.rds
# the CSV the posterior mean, median and 90% interval of each average
# selection pressure (below), by covariate, covariate group and overall, and
# the RDS its draws, a list by version of draws x covariates x types arrays
# (with the chain of each draw).
# USE_CHAINS (named_chains(), R/species_fit_helpers.R) picks a fit's chains by
# label, e.g. USE_CHAINS="V5=1,2".
#
# The average selection pressure of covariate k on insecticide type t is
# exp(beta_kt) times the mean of covariate k over the modelled cell-years
# (x_cell_years: every modelled bioassay cell, every year from the baseline),
# per posterior draw; the grids show its posterior mean. With a spatial log
# multiplier of selection u(x), it is computed two ways:
#   beta     the multiplier at 0: exp(beta_kt) alone, as above
#   posthoc  per draw, each cell-year's covariate times the multiplier at its
#            cell, exp(u(x)), averaged over the same cell-years:
#            exp(beta_kt) mean(x_k exp(u(x)))
# The multiplier is the latent smooth of selection u_s(x) (V5; for the classes
# it applies to) or, for one trajectory with the kdr covariate (V2, V4), its
# kdr term delta_selection k(x) (for every class). The species model's
# multipliers differ by species and are not applied. For small per-year
# selection z = sum_k x_k exp(beta_k), exp(u) log(1 + z) is close to
# log(1 + exp(u) z), so the multiplier acts as a common shift of every beta at
# the cell. With a multiplier, the grid figure shows both, side by side.
#
# Plain R; about 4 GB and 1 minute for a fit with the latent smooths.

arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) %in% c(0, 2))
labelled <- length(arguments) == 2
file <- if (labelled) arguments[1] else "temporary/fitted_model.RData"
label <- if (labelled) arguments[2] else ""
n_draws <- 4000

suppressMessages({
  library(greta)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
source("R/dynamical_predictions.R")
source("R/species_fit_helpers.R")

if (labelled) {
  figure_dir <- "figures/species_runs/grids"
  output_dir <- "outputs/species_runs/grids"
  dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  loadings_file <- file.path(figure_dir,
                             sprintf("covariate_loadings_%s.png", label))
  effect_size_file <- file.path(
    figure_dir, sprintf("covariate_loadings_effect_size_%s.png", label))
  summary_file <- file.path(output_dir,
                            sprintf("selection_pressure_%s.csv", label))
  draws_file <- file.path(output_dir,
                          sprintf("selection_pressure_draws_%s.rds", label))
} else {
  loadings_file <- "figures/covariate_loadings.png"
  effect_size_file <- "figures/covariate_loadings_effect_size.png"
}

# the posterior draws of the selection effects, from the fit's own transforms
# (dynamical_parameter_draws(), R/dynamical_predictions.R), so that centred
# and non-centred fits are read alike: about n_draws, evenly spaced in each
# usable chain
fit <- load_fit(file)
draws_used <- even_draws(fit, n_draws, label)
parameters <- fit_parameter_draws(fit, draws_used$index)
effect <- parameters$effect_type
types <- fit$types
covariates <- colnames(fit$x_cell_years)
dimnames(effect) <- list(NULL, covariates, types)
report("%s: %d draws of chains %s", if (labelled) label else file,
       dim(effect)[1], toString(sort(unique(draws_used$chain))))

# the covariates' sums over the modelled cell-years at each cell (rows of
# x_cell_years are cell-major, years fastest; build_dynamical_model()), and
# their means over all of them
x_cell_years <- as.matrix(fit$x_cell_years)
cell_id <- fit$cell_years_index$cell_id
n_cells <- max(cell_id)
x_cell_sums <- rowsum(x_cell_years, cell_id, reorder = TRUE)
stopifnot(nrow(x_cell_sums) == n_cells)
data_mean <- colMeans(x_cell_years)

# the spatial log multiplier of selection at each modelled cell, draws x
# cells, and the types it applies to; NULL without one
cells <- fit$df$cell[match(seq_len(n_cells), fit$df$cell_id)]
multiplier <- NULL
if ("selection" %in% smooth_kinds(fit$options)) {
  basis <- prediction_basis(fit$options, cells)
  multiplier <- list(
    log = parameters$smooth_weights$selection %*% t(basis),
    types = if (isTRUE(fit$options$smooth$selection)) {
      rep(TRUE, length(types))
    } else {
      fit$options$smooth$term_classes[fit$classes_index]
    },
    name = "u_s(x)")
} else if (kdr_on(fit$options) && !species_on(fit$options) &&
           length(parameters$kdr_slopes$delta_selection) > 0) {
  k <- prediction_kdr(fit$options, cells)[, "complex"]
  multiplier <- list(log = outer(parameters$kdr_slopes$delta_selection, k),
                     types = rep(TRUE, length(types)),
                     name = "delta k(x)")
}

# the average selection pressure per draw, draws x covariates x types, with
# the multiplier at 0 ("beta") and, with one, applied post hoc ("posthoc")
pressure <- list(beta = sweep(effect, 2, data_mean, FUN = "*"))
if (!is.null(multiplier)) {
  weighted_mean <- exp(multiplier$log) %*% x_cell_sums / nrow(x_cell_years)
  posthoc <- pressure$beta
  for (t in which(multiplier$types)) {
    posthoc[, , t] <- effect[, , t] * weighted_mean
  }
  pressure$posthoc <- posthoc
  m_cells <- rowMeans(exp(multiplier$log))
  report("mean exp(%s) over the %d modelled cells: posterior median %.3g (90%%: %.3g to %.3g)",
         multiplier$name, n_cells, median(m_cells),
         quantile(m_cells, 0.05), quantile(m_cells, 0.95))
}

# plot labels of the covariates, and their groups: the population column
# carries its transform and trend, e.g. pop_enc:g_dom, and the crops their
# trend, e.g. rice:g_ag (selection_column_names())
covariate_label <- function(covariate) {
  case_when(
    str_detect(covariate, "^pop") ~ "Human population",
    covariate == "nets" ~ "LLIN use",
    covariate == "irs" ~ "IRS coverage",
    .default = str_to_sentence(str_remove(covariate, ":g_ag$"))
  )
}
covariate_group <- function(label) {
  case_when(
    label %in% c("LLIN use", "IRS coverage") ~ "Vector control",
    label == "Human population" ~ "Human population",
    .default = "Agriculture"
  )
}
labels <- covariate_label(covariates)
groups <- covariate_group(labels)

# posterior summaries of a draws x rows x types array, as a long tibble
summarise_draws <- function(x, level) {
  stat <- function(f, name) {
    as_tibble(apply(x, 2:3, f), rownames = "row") %>%
      pivot_longer(-row, names_to = "insecticide", values_to = name)
  }
  stat(mean, "mean") %>%
    left_join(stat(median, "median"), by = c("row", "insecticide")) %>%
    left_join(stat(function(v) quantile(v, 0.05, names = FALSE), "q05"),
              by = c("row", "insecticide")) %>%
    left_join(stat(function(v) quantile(v, 0.95, names = FALSE), "q95"),
              by = c("row", "insecticide")) %>%
    mutate(level = level, .before = everything())
}
# the sums of a draws x covariates x types array over the covariates in each
# level of `index`, per draw
sum_over <- function(x, index) {
  index <- factor(index)
  out <- vapply(levels(index), function(level) {
    apply(x[, index == level, , drop = FALSE], c(1, 3), sum)
  }, array(0, dim(x)[c(1, 3)]))
  out <- aperm(out, c(1, 3, 2))
  dimnames(out) <- list(NULL, levels(index), dimnames(x)[[3]])
  out
}
summaries <- bind_rows(lapply(names(pressure), function(version) {
  x <- pressure[[version]]
  dimnames(x) <- list(NULL, labels, types)
  bind_rows(
    summarise_draws(sum_over(x, rep("Overall", length(labels))), "overall"),
    summarise_draws(sum_over(x, groups), "group"),
    summarise_draws(x, "covariate")) %>%
    mutate(version = version, .before = everything())
})) %>%
  mutate(class = fit$classes[fit$classes_index[match(insecticide, types)]],
         .after = insecticide)
if (labelled) {
  write.csv(mutate(summaries, label = label, .before = everything()),
            summary_file, row.names = FALSE)
  saveRDS(c(pressure, list(chain = draws_used$chain)), draws_file)
}

# the posterior mean effect sizes, exp(beta), and the order of the
# covariates, by their mean pressure with the multiplier at 0
effect_size <- as_tibble(apply(effect, 2:3, mean)) %>%
  mutate(covariate = labels) %>%
  pivot_longer(-covariate, names_to = "insecticide", values_to = "mean")
covariate_order_scaled <- summaries %>%
  filter(version == "beta", level == "covariate") %>%
  group_by(row) %>%
  summarise(mean = mean(mean)) %>%
  arrange(mean) %>%
  pull(row)


# the grids ---------------------------------------------------------------------

# a grid of `stats` (insecticide, row and value), the rows in `order`
tile_grid <- function(stats, order, fill_name) {
  stats %>%
    mutate(row = factor(row, levels = order)) %>%
    ggplot(
      aes(
        x = insecticide,
        y = row,
        fill = value
      )
    ) +
    geom_tile(
      colour = grey(0.5)
    ) +
    scale_x_discrete(
      position = "top"
    ) +
    theme_minimal() +
    xlab("") +
    ylab("") +
    labs(fill = fill_name) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 0,
        vjust = 0)
    )
}

# The average contribution of each covariate to the selection coefficient for
# resistance to each insecticide. Computed as the posterior mean of the partial
# selection coefficient, multiplied by the mean over the dataset of the
# indicator variable for the covariate. This captures both the relative strength
# of the covariate in selecting for resistance, and its prevalence. Overall
# (A), by group of covariates (B) and by covariate (C), for one version, with
# the fill scale from 0 to its largest value
pressure_grids <- function(version, title = NULL) {
  stats <- summaries %>%
    filter(version == !!version) %>%
    transmute(level, row, insecticide, value = mean)
  fill_name <- "Average\nselection\npressure"
  # suppress some plot elements in each
  fill_scale <- scale_fill_gradient(low = "white",
                                    high = "red",
                                    limits = c(0, max(stats$value)))
  overall <- tile_grid(filter(stats, level == "overall"), "Overall",
                       fill_name)
  combined <- tile_grid(filter(stats, level == "group"),
                        c("Human population", "Agriculture",
                          "Vector control"),
                        fill_name)
  specific <- tile_grid(filter(stats, level == "covariate"),
                        covariate_order_scaled, fill_name)
  (overall +
      fill_scale +
      labs(title = title) +
      theme(legend.position = "none",
            plot.title.position = "panel")) +
    (combined +
       fill_scale +
       theme(axis.text.x = element_blank(),
             legend.position = "none",
             plot.title.position = "plot")) +
    (specific +
       fill_scale +
       theme(axis.text.x = element_blank(),
             plot.title.position = "plot")) +
    plot_layout(nrow = 3, heights = c(1, 3, 9))
}

# plot the combined and specific selection pressures: with a multiplier, the
# multiplier at 0 and post hoc side by side, each on its own scale
if (is.null(multiplier)) {
  loadings_plot <- pressure_grids("beta", if (labelled) label) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 12))
  width <- 8
} else {
  loadings_plot <- (pressure_grids(
    "beta", sprintf("%s: beta alone (%s = 0)", label, multiplier$name)) |
      pressure_grids(
        "posthoc", sprintf("%s: x exp(%s) per cell, post hoc", label,
                           multiplier$name))) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 12))
  width <- 16
}

ggsave(
  loadings_file,
  loadings_plot,
  bg = "white",
  scale = 0.8,
  width = width,
  height = 8
)


# do the same on the coefficients themselves
regression_grid <- effect_size %>%
  transmute(row = covariate, insecticide, value = mean) %>%
  tile_grid(covariate_order_scaled, "Effect\nsize") +
  scale_fill_gradient(
    low = "white",
    high = "red"
  ) +
  labs(title = if (labelled) label)

ggsave(
  effect_size_file,
  regression_grid,
  bg = "white",
  scale = 0.8,
  width = 8,
  height = 6
)
report("wrote %s and %s", loadings_file, effect_size_file)
