# Plots of the maximum-likelihood hinge fits (R/hinge_mle.R): for each
# insecticide type, the fitted yearly selection f(x) = exp(b) x +
# sum_j exp(c_j) min(x, k_j) of each model against the covariate, with the
# fit's own linear effect (dashed) and each pair's partial residual (points):
# its observed fall in the empirical logit of mortality per year, minus the
# fit's predicted fall, plus the fit's linear effect of the covariate, against
# the covariate's mean over the pair's span. For population both are per unit
# of the time trend g_dom (the partial residual divided by the span's mean
# g_dom), the scale of f.
#
#   Rscript R/hinge_mle_plot.R <fit RData> <label> <mle output label> [gap]
#
# Writes figures/hinge_mle_<mle output label>_<covariate>.png.

arguments <- commandArgs(trailingOnly = TRUE)
fit_file <- arguments[1]
label <- arguments[2]
mle_label <- arguments[3]
gap <- if (length(arguments) >= 4) as.integer(arguments[4]) else 5L
output_dir <- "outputs/hinge_screen"

source("R/packages.R")
source("R/functions.R")
source("R/dynamical_model.R")
source("R/dynamical_predictions.R")

f <- new.env()
load(fit_file, envir = f)
fold <- list(draws = f$draws, options = f$model_options,
             x_cells_init = f$x_cells_init)
n_total <- sum(vapply(fold$draws, nrow, integer(1)))
par <- dynamical_parameter_draws(fold, f$classes_index, f$types,
                                 draw_index = round(seq(1, n_total,
                                                        length.out = 400)))
effect <- apply(par$effect_type, c(2, 3), mean)
dimnames(effect) <- list(colnames(f$x_cell_years), f$types)
pop_column <- grep("^pop_[a-z]+:g_dom$", rownames(effect), value = TRUE)
fit_effect <- list(nets = effect["nets", ], pop = effect[pop_column, ])
rm(f, par, fold)

stem <- if (gap > 1) sprintf("%s_gap%d", label, gap) else label
pairs <- read.csv(file.path(output_dir, sprintf("%s_change_pairs.csv", stem)))
years <- read.csv(file.path(output_dir, sprintf("%s_change_pair_years.csv",
                                                stem)))
groups <- read.csv(file.path(output_dir, sprintf("%s_mle_by_group.csv",
                                                 mle_label)))
groups <- groups[groups$level == "type" & groups$gap == gap, ]
span <- pairs$year_2 - pairs$year_1
fall <- (-pairs$change_obs + pairs$change_pred) / span

for (covariate in c("nets", "pop")) {
  g <- if (covariate == "pop") years$g_dom else 1
  mean_by_pair <- function(v) {
    tapply(v, years$pair, mean)[as.character(pairs$pair)]
  }
  x_mean <- mean_by_pair(years[[covariate]])
  gx_mean <- mean_by_pair(g * years[[covariate]])
  g_mean <- mean_by_pair(g + 0 * years$pair)
  points <- tibble(
    type = pairs$insecticide_type, x = x_mean,
    y = (fall + fit_effect[[covariate]][pairs$insecticide_type] * gx_mean) /
      g_mean) |>
    filter(is.finite(y))
  rows <- groups[groups$covariate == covariate, ]
  at <- seq(0, quantile(points$x, 0.995), length.out = 200)
  curves <- bind_rows(lapply(seq_len(nrow(rows)), function(i) {
    coefficients <- as.numeric(strsplit(rows$coefficients[i], ",")[[1]])
    knots <- if (rows$knots[i] == "") numeric(0) else
      as.numeric(strsplit(rows$knots[i], ",")[[1]])
    y <- coefficients[1] * at
    for (j in seq_along(knots)) {
      y <- y + coefficients[j + 1] * pmin(at, knots[j])
    }
    tibble(type = rows$group[i], x = at, y = y,
           model = if (rows$hinges[i] == 0) "linear" else
             sprintf("%d hinge(s) at %s", rows$hinges[i], rows$knots[i]))
  }))
  # each type's AIC change from linear, by model
  labels <- rows |>
    group_by(type = group) |>
    summarise(text = paste(sprintf("%dh %+.1f", hinges[hinges > 0],
                                   delta_aic[hinges > 0]), collapse = "\n"),
              .groups = "drop") |>
    mutate(text = paste0("dAIC vs linear\n", text))
  reference <- tibble(type = names(fit_effect[[covariate]])) |>
    reframe(x = at, y = fit_effect[[covariate]][type] * at, .by = type)
  limits <- quantile(points$y, c(0.02, 0.98))
  p <- ggplot(points, aes(x, y)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_point(size = 0.4, alpha = 0.25, colour = "grey30") +
    geom_line(data = reference, linetype = "dashed", colour = "grey20") +
    geom_line(data = curves, aes(colour = model), linewidth = 0.8) +
    geom_text(data = labels, aes(x = -Inf, y = Inf, label = text),
              hjust = -0.05, vjust = 1.1, size = 2.6, lineheight = 0.9) +
    facet_wrap(~ type, ncol = 3) +
    coord_cartesian(ylim = limits) +
    labs(x = if (covariate == "nets") "net use (mean over the pair's span)" else
           "population, encounter transform (mean over the pair's span)",
         y = "selection per year (logit mortality)", colour = NULL,
         title = sprintf("Monotone concave hinge fits, %s, %d-year pairs, %s",
                         label, gap, covariate),
         subtitle = paste("Lines: fitted selection per year; dashed: the",
                          "fit's linear effect; points: pairs' partial",
                          "residuals")) +
    theme_minimal() +
    theme(legend.position = "bottom")
  ggsave(sprintf("figures/hinge_mle_%s_%s.png", mle_label, covariate), p,
         width = 10, height = 11, dpi = 130, bg = "white")
}
