# Plots of the hinge selection (R/hinge_select.R): in each insecticide class
# and type, every pair's change residual z against the covariate, with the
# fitted models with 0-3 hinges (the knots chosen at that level, shared across
# its groups), labelled with each group's change in AIC from the linear model.
#
#   Rscript R/hinge_select_plot.R <pairs csv> [hinge counts to draw]
#
# Writes figures/hinge_select_<pairs file label>_<covariate>.png and
# <pairs csv without .csv>_hinges_by_group.csv.

arguments <- commandArgs(trailingOnly = TRUE)
pairs_file <- arguments[1]
show_hinges <- if (length(arguments) >= 2) as.integer(arguments[-1]) else 0:3
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
pairs <- read.csv(pairs_file)
pairs <- pairs[is.finite(pairs$z), ]
hinges <- read.csv(sub("\\.csv$", "_hinges.csv", pairs_file))
label <- sub("_change_pairs\\.csv$", "", basename(pairs_file))
pairs$class <- ifelse(pairs$insecticide_type == "DDT", "DDT",
                      pairs$insecticide_class)

# the fitted line of one group with knots `knots`, and its AIC term
group_line <- function(d, covariate, knots) {
  x <- d[[covariate]]
  design <- function(x) cbind(1, x, vapply(knots, function(k) pmin(x, k),
                                           numeric(length(x))))
  fit <- lm.fit(design(x), d$z)
  at <- seq(min(x), max(x), length.out = 100)
  list(line = tibble(at = at, fit = c(design(at) %*%
                                        ifelse(is.na(fit$coefficients), 0,
                                               fit$coefficients))),
       aic = nrow(d) * log(sum(fit$residuals ^ 2) / nrow(d)) + 2 * fit$rank)
}

lines <- list()
scores <- list()
for (covariate in unique(hinges$covariate)) {
  for (level in c("class", "type")) {
    group_var <- if (level == "class") "class" else "insecticide_type"
    models <- hinges[hinges$level == level & hinges$covariate == covariate, ]
    for (group in sort(unique(pairs[[group_var]]))) {
      d <- pairs[pairs[[group_var]] == group, ]
      for (i in seq_len(nrow(models))) {
        knots <- if (models$knots[i] == "") numeric(0) else
          as.numeric(strsplit(models$knots[i], ",")[[1]])
        g <- group_line(d, covariate, knots)
        model <- sprintf("%d: %s", models$hinges[i],
                         if (length(knots) == 0) "linear" else
                           models$knots[i])
        lines[[length(lines) + 1]] <- mutate(
          g$line, covariate = covariate, level = level, group = group,
          hinges = models$hinges[i], model = model)
        scores[[length(scores) + 1]] <- tibble(
          covariate = covariate, level = level, group = group, n = nrow(d),
          hinges = models$hinges[i], knots = models$knots[i], aic = g$aic)
      }
    }
  }
}
lines <- bind_rows(lines)
scores <- bind_rows(scores) |>
  group_by(covariate, level, group) |>
  mutate(delta_aic_vs_linear = aic - aic[hinges == 0]) |>
  ungroup()
write.csv(scores, sub("\\.csv$", "_hinges_by_group.csv", pairs_file),
          row.names = FALSE)

for (covariate in unique(lines$covariate)) {
  panel <- function(level, group) sprintf("%s: %s", level, group)
  l <- lines[lines$covariate == covariate & lines$hinges %in% show_hinges, ] |>
    mutate(panel = panel(level, group))
  s <- scores[scores$covariate == covariate, ] |>
    group_by(level, group, n) |>
    summarise(text = paste(sprintf("%dh %+.1f", hinges[hinges %in% setdiff(show_hinges, 0)],
                                   delta_aic_vs_linear[hinges %in% setdiff(show_hinges, 0)]),
                           collapse = "\n"), .groups = "drop") |>
    mutate(panel = panel(level, group),
           text = sprintf("n = %d; dAIC vs linear\n%s", n, text))
  points <- bind_rows(
    pairs |> transmute(panel = panel("class", class), x = .data[[covariate]],
                       z),
    pairs |>
      transmute(panel = panel("type", insecticide_type),
                x = .data[[covariate]], z))
  panel_levels <- c(sort(unique(points$panel[startsWith(points$panel,
                                                        "class")])),
                    sort(unique(points$panel[startsWith(points$panel,
                                                        "type")])))
  l$panel <- factor(l$panel, panel_levels)
  s$panel <- factor(s$panel, panel_levels)
  points$panel <- factor(points$panel, panel_levels)
  models <- unique(l[, c("level", "hinges", "model")])
  file_suffix <- if (identical(sort(show_hinges), 0:3)) "" else
    paste0("_h", paste(show_hinges, collapse = ""))
  p <- ggplot(points, aes(x, z)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_point(size = 0.4, alpha = 0.25, colour = "grey30") +
    geom_line(data = l, aes(at, fit, colour = factor(hinges)),
              linewidth = 0.7) +
    geom_text(data = s, aes(x = -Inf, y = Inf, label = text), hjust = -0.05,
              vjust = 1.1, size = 2.4, lineheight = 0.9) +
    facet_wrap(~ panel, ncol = 3) +
    coord_cartesian(ylim = c(-3, 3)) +
    scale_colour_manual(values = c(`0` = "black", `1` = "#0072B2",
                                   `2` = "#D55E00", `3` = "#009E73")) +
    labs(x = covariate, y = "z = qnorm(PIT) of observed change",
         colour = "hinges",
         title = sprintf("Hinge models, %s, %s", label, covariate),
         subtitle = paste(
           "Knots shared across groups at each level:",
           paste(sprintf("%s %s", models$level, models$model),
                 collapse = "; "))) +
    theme_minimal() +
    theme(legend.position = "top",
          plot.subtitle = element_text(size = 7))
  ggsave(sprintf("figures/hinge_select_%s_%s%s.png", label, covariate,
                 file_suffix), p,
         width = 10, height = 14, dpi = 130, bg = "white")
}
