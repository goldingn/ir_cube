# The misfit of a fit's within-pixel changes on the selection scale (#37):
# for each pair (R/change_residual_screen.R), the observed fall in the
# empirical logit of mortality per year minus the fit's predicted fall
# (negative: the fit over-predicts selection), against net use and
# population, with a GAM weighted by the pair's precision (the empirical-logit
# variances of its two years, inflated by the type's beta-binomial rho). Unlike
# z (qnorm of the PIT), it is not divided by each pair's predictive spread,
# which depends on the mortality level and so on time.
#
#   Rscript R/change_misfit_plot.R <label> [gap]
#
# Writes figures/change_misfit_<label>[_gap<g>]_<covariate>.png.

arguments <- commandArgs(trailingOnly = TRUE)
label <- arguments[1]
gap <- if (length(arguments) >= 2) as.integer(arguments[2]) else 5L
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(mgcv)
})
stem <- if (gap > 1) sprintf("%s_gap%d", label, gap) else label
pairs <- read.csv(file.path("outputs/hinge_screen",
                            sprintf("%s_change_pairs.csv", stem)))
span <- pairs$year_2 - pairs$year_1
# the empirical-logit variance of a year's pooled assays
elogit_var <- function(died, n, rho) {
  (1 / (died + 0.5) + 1 / (n - died + 0.5)) * (1 + (n - 1) * rho)
}
pairs$w <- span ^ 2 / (elogit_var(pairs$died_1, pairs$n_1, pairs$rho) +
                         elogit_var(pairs$died_2, pairs$n_2, pairs$rho))
pairs$misfit <- (-pairs$change_obs + pairs$change_pred) / span

pyrethroid_types <- sort(unique(pairs$insecticide_type[
  pairs$insecticide_class == "Pyrethroids"]))
groups <- c("Pyrethroids", "DDT", "Carbamates", "Organophosphates",
            pyrethroid_types)
group_rows <- function(name) {
  if (name %in% pairs$group) pairs$group == name else
    pairs$insecticide_type == name
}
labels <- c(nets = "net use (mean over the pair's span)",
            pop = "population, encounter transform (mean over the pair's span)")
for (covariate in c("nets", "pop")) {
  points <- list()
  curves <- list()
  for (name in groups) {
    d <- pairs[group_rows(name), ]
    d$x <- d[[covariate]]
    fit <- gam(misfit ~ s(x, k = 8), data = d, weights = w / mean(w),
               method = "REML")
    at <- seq(quantile(d$x, 0.01), quantile(d$x, 0.99), length.out = 100)
    pr <- predict(fit, data.frame(x = at), se.fit = TRUE)
    curves[[name]] <- tibble(group = name, x = at, fit = pr$fit,
                             lower = pr$fit - 1.96 * pr$se.fit,
                             upper = pr$fit + 1.96 * pr$se.fit)
    points[[name]] <- tibble(group = name, x = d$x, misfit = d$misfit,
                             w = d$w)
  }
  points <- bind_rows(points) |> mutate(group = factor(group, groups))
  curves <- bind_rows(curves) |> mutate(group = factor(group, groups))
  p <- ggplot(points, aes(x, misfit)) +
    geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
    geom_point(aes(size = w), alpha = 0.2, colour = "grey30",
               show.legend = FALSE) +
    scale_size_area(max_size = 1.5) +
    geom_ribbon(data = curves, aes(x, ymin = lower, ymax = upper),
                inherit.aes = FALSE, fill = "#0072B2", alpha = 0.3) +
    geom_line(data = curves, aes(x, fit), colour = "#0072B2",
              linewidth = 0.8) +
    facet_wrap(~ group, ncol = 3) +
    coord_cartesian(ylim = c(-1.5, 1.5)) +
    labs(x = labels[[covariate]],
         y = "observed minus predicted fall in logit mortality per year",
         title = sprintf("Change misfit of %s, %d-year pairs, against %s",
                         label, gap, covariate),
         subtitle = paste("Negative: the fit over-predicts selection.",
                          "Points sized by precision; line: GAM weighted",
                          "by precision")) +
    theme_minimal()
  ggsave(sprintf("figures/change_misfit_%s_%s.png", stem, covariate), p,
         width = 10, height = 12, dpi = 130, bg = "white")
}
