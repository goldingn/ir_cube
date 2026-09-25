# The bar encoding shared by every variance-explained figure.
#
# Each bar spans the full 100% of observed variance in the held-out quantity.
# The bar is light grey, and the part bioassay sampling makes unexplainable is
# washed out toward white from the top, so what remains grey is the predictable
# part - the quantity the reader should be comparing against. The model's share
# is filled in colour from the bottom.
#
#   light grey, full height: observed variance out of sample
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
#
# Sourcing this file defines the colours, bar_layers(), region_key() and
# base_theme; it draws nothing on its own.
suppressMessages({
  library(dplyr)
  library(ggplot2)
})

# Blue for the mechanistic model and green for the baseline a person would
# actually apply, which are the two the reader is asked to compare; the other
# two models are greys so they read as reference rather than as competitors.
model_colours <- c(
  "dynamical model"         = "#2166AC",
  "nearest recent survey"   = "#1B7837",
  "nearest surveys, best k" = "#8073AC",
  "insecticide mean"        = grey(0.45))

bar_background <- grey(0.88)
wash_colour <- "white"
wash_solid <- 0.82
wash_interval <- 0.42

solid_alpha <- 1
interval_alpha <- 0.35

# one bar: the coloured base, its translucent interval, the rule at the point
# estimate, and the noise block from the top
bar_layers <- function(data, width = 0.8) {

  models <- data %>% filter(kind == "model")
  noise <- data %>% filter(kind == "noise")

  list(
    # the full bar: all the variance there is to explain
    geom_rect(data = noise,
              aes(xmin = position - width/2, xmax = position + width/2,
                  ymin = 0, ymax = 100),
              fill = bar_background, colour = NA),
    # from the bottom: what the model explains
    geom_rect(data = models,
              aes(xmin = position - width/2, xmax = position + width/2,
                  ymin = 0, ymax = pmax(lower, 0), fill = quantity),
              alpha = solid_alpha, colour = NA),
    geom_rect(data = models,
              aes(xmin = position - width/2, xmax = position + width/2,
                  ymin = pmax(lower, 0), ymax = upper, fill = quantity),
              alpha = interval_alpha, colour = NA),
    geom_segment(data = models,
                 aes(x = position - width/2, xend = position + width/2,
                     y = estimate, yend = estimate),
                 colour = "white", linewidth = 0.5),
    # washed toward white from the top, over everything below, so the
    # unpredictable share recedes and anything reaching into it is washed too
    geom_rect(data = noise,
              aes(xmin = position - width/2, xmax = position + width/2,
                  ymin = 100 - lower, ymax = 100),
              fill = wash_colour, alpha = wash_solid, colour = NA),
    geom_rect(data = noise,
              aes(xmin = position - width/2, xmax = position + width/2,
                  ymin = 100 - upper, ymax = 100 - lower),
              fill = wash_colour, alpha = wash_interval, colour = NA),
    # the ceiling: the most any model could explain
    geom_segment(data = noise,
                 aes(x = position - width/2, xend = position + width/2,
                     y = 100 - estimate, yend = 100 - estimate),
                 colour = grey(0.45), linewidth = 0.35, linetype = "22"),
    # the bar outline last, so the 100% extent stays legible
    geom_rect(data = noise,
              aes(xmin = position - width/2, xmax = position + width/2,
                  ymin = 0, ymax = 100),
              fill = NA, colour = grey(0.7), linewidth = 0.3)
  )
}

# The three regions labelled directly, as vertical bars in the space to the
# right of the last bar group, so that the meaning of the heights is read off
# the figure rather than from a caption. Drawn against a reference bar - the
# dynamical model in that panel - because the boundaries differ between bars.
region_key <- function(reference_estimate, reference_noise, x,
                       label_size = 2.9, wrap = FALSE) {
  ceiling_value <- 100 - reference_noise
  # a label that is longer than the band it names overruns its own key line,
  # so it is broken over two lines instead: that halves what it needs along
  # the axis, at the cost of a wider column in the margin. Which ones need it
  # depends on the panel height and on where the reference bar's boundaries
  # fall, so it is set per figure; the noise band is the shortest in every
  # figure, so that one is always broken.
  wrap <- rep_len(wrap, 3)
  wrap[3] <- TRUE
  labels <- ifelse(wrap,
                   sub(" ", "\n",
                       c("explained variance", "unexplained variance",
                         "bioassay noise")),
                   c("explained variance", "unexplained variance",
                     "bioassay noise"))
  regions <- data.frame(
    lower = c(0, reference_estimate, ceiling_value),
    upper = c(reference_estimate, ceiling_value, 100),
    label = labels,
    colour = c("#2166AC", grey(0.45), grey(0.62)))
  # each line runs the full height of its own region on the outer edge - 0%
  # for the explained share, 100% for the noise share - and is trimmed only
  # where it meets the next region, which is all the whitespace needed to read
  # them as three regions rather than one continuous rule
  gap <- 1.6
  regions$from <- c(regions$lower[1], regions$lower[2:3] + gap)
  regions$to <- c(regions$upper[1:2] - gap, regions$upper[3])
  # each label is justified to its region's outer edge, the way the region
  # itself is stacked: the explained share fills up from 0%, the noise share
  # down from 100%, and the unexplained share is what is left in the middle
  regions$anchor <- c(regions$from[1],
                      mean(c(regions$from[2], regions$to[2])),
                      regions$to[3])
  regions$hjust <- c(1, 0.5, 0)
  # every label sits the same distance off its line. vjust is a fraction of
  # the label's own height, so a two-line label has to be given half the
  # fraction to leave the same gap as a one-line one
  n_lines <- lengths(regmatches(regions$label, gregexpr("\n", regions$label))) + 1
  regions$vjust <- -0.55 / n_lines
  list(
    geom_segment(data = regions,
                 aes(x = x, xend = x, y = from, yend = to),
                 colour = regions$colour, linewidth = 0.9,
                 inherit.aes = FALSE),
    # set along the lines rather than beside them, so the key needs only a
    # line's width of margin. At angle 270 the text reads downward and its own
    # "up" direction points right, so vjust clears it of the line in units of
    # its own height - which, unlike an offset in data units, does not have to
    # be retuned for every panel width.
    geom_text(data = regions,
              aes(x = x, y = anchor, label = label, hjust = hjust,
                  vjust = vjust),
              angle = 270, size = label_size,
              lineheight = 0.9, colour = regions$colour, inherit.aes = FALSE)
  )
}

base_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 9),
        legend.position = "none",
        # flush with the panel edge, so the label starts at the first bar
        strip.text = element_text(face = "bold", hjust = 0,
                                  margin = margin(b = 5, l = 0)),
        panel.spacing.y = unit(14, "pt"),
        # the key is drawn past the right edge of the last panel, with
        # clipping off, so the margin has to carry it
        plot.margin = margin(6, 44, 6, 6))

