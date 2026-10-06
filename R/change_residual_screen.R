# Change-residual screen (#37): how well a saved full fit predicts the change
# in mortality between consecutive sampled years within a pixel, against net
# use and the transformed population covariate, to choose hinge knots
# (selection_design(hinges = ...), columns min(x, k)).
#
#   Rscript R/change_residual_screen.R [fit file] [label] [min gap, years]
#
# Defaults: outputs/refit/fits/r2_sens_dhalf200.RData (round 2 dynamical fit,
# population encounter transform with d_half = 200), label "dhalf200".
#
# Assays are pooled per pixel, insecticide type and year. Within each pixel
# and type, each pair of sampled years (y1, y2) at least min_gap years apart
# (default 1: consecutive sampled years; see gap_chain()) gives an observed
# change on the logit scale,
#
#   delta_obs = elogit(y2) - elogit(y1),   elogit = logit((died + 0.5) / (n + 1))
#
# and a predictive distribution of the same quantity: per posterior draw,
# every assay in both years is simulated from the beta-binomial with the
# draw's predicted mortality and rho, pooled, and differenced the same way.
# The residual is the randomised PIT of delta_obs in that distribution,
# gaussianised: z = qnorm(PIT). z > 0: mortality rose more (or fell less) than
# predicted, i.e. the fit over-predicts selection over the interval. A
# calibrated fit gives z ~ N(0, 1).
#
# The covariates of a pair are their means over the years the change spans,
# y1 + 1 to y2: the model's nets column, and population on its 0-1 transformed
# scale (the model's population column divided by its time trend g_dom(t)),
# the scale a hinge knot is set on. z is splined against both (gaussian GAM)
# by insecticide group.
#
# Writes outputs/hinge_screen/<label>_change_pairs.csv,
# <label>_change_curves.csv, <label>_change_gams.csv and
# figures/change_residuals_<label>.png.

arguments <- commandArgs(trailingOnly = TRUE)
file <- if (length(arguments) >= 1) arguments[1] else
  "outputs/refit/fits/r2_sens_dhalf200.RData"
label <- if (length(arguments) >= 2) arguments[2] else "dhalf200"
min_gap <- if (length(arguments) >= 3) as.integer(arguments[3]) else 1L
label <- if (min_gap > 1) sprintf("%s_gap%d", label, min_gap) else label
n_draws_keep <- 400
output_dir <- "outputs/hinge_screen"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(37)

source("R/packages.R")
source("R/functions.R")
source("R/dynamical_model.R")
source("R/dynamical_predictions.R")
suppressMessages(library(mgcv))

f <- new.env()
load(file, envir = f)
df <- f$df

# draws of predicted mortality and rho at each assay
fold <- list(draws = f$draws, options = f$model_options,
             x_cells_init = f$x_cells_init)
n_total <- sum(vapply(fold$draws, nrow, integer(1)))
par <- dynamical_parameter_draws(fold, f$classes_index, f$types,
                                 draw_index = round(seq(1, n_total,
                                                        length.out =
                                                          n_draws_keep)))
p <- plogis(dynamical_logit(par, select(df, cell_id, country_id, type_id,
                                        year_id),
                            df, f$x_cell_years, f$cell_years_index))
rho <- par$rho_types[, df$type_id, drop = FALSE]
floor_mean <- mean(par$mortality_floor)

# simulated deaths, draws x assays, from the beta-binomial with mean p and
# intra-assay correlation rho
size <- matrix(df$mosquito_number, nrow(p), ncol(p), byrow = TRUE)
shape <- 1 / rho - 1
died_sim <- matrix(extraDistr::rbbinom(length(p), size, p * shape,
                                       (1 - p) * shape),
                   nrow(p))
rm(size, shape, rho)
invisible(gc())

# pooled per pixel, type and year
key <- paste(df$cell_id, df$type_id, df$year_id)
pooled <- df |>
  mutate(key = key) |>
  group_by(key, cell_id, type_id, year_id, insecticide_class,
           insecticide_type) |>
  summarise(died = sum(died), n = sum(mosquito_number), .groups = "drop")
died_sim <- t(rowsum(t(died_sim), key))[, pooled$key, drop = FALSE]

# the empirical logit, kept finite at 0 and n deaths
elogit <- function(died, n) qlogis((died + 0.5) / (n + 1))

# Pairs of sampled years within each pixel and type, at least min_gap years
# apart: a chain from the earliest sampled year to the first one min_gap or
# more years later, and on from there, so that the intervals do not overlap
gap_chain <- function(years, min_gap) {
  keep <- 1L
  for (i in seq_along(years)[-1]) {
    if (years[i] - years[keep[length(keep)]] >= min_gap) keep <- c(keep, i)
  }
  keep
}
pairs <- pooled |>
  mutate(index = row_number()) |>
  arrange(cell_id, type_id, year_id) |>
  group_by(cell_id, type_id) |>
  filter(row_number() %in% gap_chain(year_id, min_gap)) |>
  mutate(index_1 = lag(index), year_id_1 = lag(year_id)) |>
  ungroup() |>
  filter(!is.na(index_1)) |>
  transmute(cell_id, type_id, insecticide_class, insecticide_type,
            year_id_1, year_id_2 = year_id, index_1, index_2 = index)

n_sims <- nrow(died_sim)
change_obs <- elogit(pooled$died[pairs$index_2], pooled$n[pairs$index_2]) -
  elogit(pooled$died[pairs$index_1], pooled$n[pairs$index_1])
change_sim <- elogit(died_sim[, pairs$index_2, drop = FALSE],
                     matrix(pooled$n[pairs$index_2], n_sims, nrow(pairs),
                            byrow = TRUE)) -
  elogit(died_sim[, pairs$index_1, drop = FALSE],
         matrix(pooled$n[pairs$index_1], n_sims, nrow(pairs), byrow = TRUE))
below <- colSums(change_sim < matrix(change_obs, n_sims, nrow(pairs),
                                     byrow = TRUE))
equal <- colSums(change_sim == matrix(change_obs, n_sims, nrow(pairs),
                                      byrow = TRUE))
# randomised PIT, so that ties (common at 0 and 100% mortality) spread evenly
pit <- (below + runif(nrow(pairs)) * (equal + 1)) / (n_sims + 1)

# covariates averaged over the years the change spans, y1 + 1 to y2
x_row <- matrix(NA_integer_, max(f$cell_years_index$cell_id),
                max(f$cell_years_index$year_id))
x_row[cbind(f$cell_years_index$cell_id, f$cell_years_index$year_id)] <-
  seq_len(nrow(f$cell_years_index))
pop_column <- grep("^pop_[a-z]+:g_dom$", colnames(f$x_cell_years),
                   value = TRUE)
stopifnot(length(pop_column) == 1)
# g_dom(t) of the "linear_0_1" trend, 0 in 1995 and 1 in 2025
g_dom <- function(year_id) pmax((f$baseline_year + year_id - 1 - 1995) / 30, 0)
span_mean <- function(column, cell, y1, y2, scale = function(y) 1) {
  years <- (y1 + 1):y2
  mean(f$x_cell_years[x_row[cell, years], column] / scale(years))
}

pairs <- pairs |>
  mutate(year_1 = f$baseline_year + year_id_1 - 1,
         year_2 = f$baseline_year + year_id_2 - 1,
         change_obs = change_obs,
         change_pred = colMeans(change_sim),
         pit = pit, z = qnorm(pit),
         nets = mapply(span_mean, "nets", cell_id, year_id_1, year_id_2),
         pop = mapply(span_mean, pop_column, cell_id, year_id_1, year_id_2,
                      MoreArgs = list(scale = g_dom)),
         group = case_when(insecticide_class == "Pyrethroids" ~ "Pyrethroids",
                           insecticide_type == "DDT" ~ "DDT",
                           TRUE ~ insecticide_class))
# the counts of each pair and its type's rho (posterior mean), and the yearly
# covariates over the span, for fits on the observed change itself
# (R/concave_hinge_select.R)
rho_type <- colMeans(par$rho_types)
pairs <- pairs |>
  mutate(pair = row_number(),
         died_1 = pooled$died[index_1], n_1 = pooled$n[index_1],
         died_2 = pooled$died[index_2], n_2 = pooled$n[index_2],
         rho = rho_type[type_id])
pair_years <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
  years <- (pairs$year_id_1[i] + 1):pairs$year_id_2[i]
  rows <- x_row[pairs$cell_id[i], years]
  tibble(pair = i, year = f$baseline_year + years - 1,
         nets = f$x_cell_years[rows, "nets"],
         # population on its 0-1 scale, and its time trend g_dom(t)
         pop = f$x_cell_years[rows, pop_column] / g_dom(years),
         g_dom = g_dom(years))
}))
rm(change_sim, died_sim, p)
invisible(gc())
write.csv(pairs, file.path(output_dir, sprintf("%s_change_pairs.csv", label)),
          row.names = FALSE)
write.csv(pair_years, file.path(output_dir,
                                sprintf("%s_change_pair_years.csv", label)),
          row.names = FALSE)

pyrethroid_types <- sort(unique(pairs$insecticide_type[
  pairs$insecticide_class == "Pyrethroids"]))
groups <- c("Pyrethroids", "DDT", "Carbamates", "Organophosphates",
            pyrethroid_types)
group_rows <- function(name) {
  if (name %in% pairs$group) pairs$group == name else
    pairs$insecticide_type == name
}

# the fitted z along one covariate, the other at its median
curve_of <- function(fit, d, term, other) {
  at <- seq(quantile(d[[term]], 0.005), quantile(d[[term]], 0.995),
            length.out = 60)
  new <- d[rep(1, length(at)), ]
  new[[term]] <- at
  new[[other]] <- median(d[[other]])
  pr <- predict(fit, new, se.fit = TRUE)
  tibble(term = term, at = at, estimate = pr$fit,
         lower = pr$fit - 1.96 * pr$se.fit, upper = pr$fit + 1.96 * pr$se.fit)
}
binned <- function(d, term) {
  d |>
    mutate(bin = ntile(.data[[term]], 10)) |>
    group_by(bin) |>
    summarise(at = mean(.data[[term]]), estimate = mean(z), pairs = n(),
              .groups = "drop") |>
    mutate(term = term)
}

curves <- list()
points <- list()
gams <- list()
for (name in groups) {
  d <- pairs[group_rows(name), ]
  fit <- gam(z ~ s(nets, k = 8) + s(pop, k = 8), data = d, method = "REML")
  s <- summary(fit)$s.table
  gams[[name]] <- tibble(
    group = name, pairs = nrow(d), mean_z = mean(d$z), sd_z = sd(d$z),
    edf_nets = s["s(nets)", "edf"], p_nets = s["s(nets)", "p-value"],
    edf_pop = s["s(pop)", "edf"], p_pop = s["s(pop)", "p-value"],
    dev_expl = summary(fit)$dev.expl)
  curves[[name]] <- bind_rows(curve_of(fit, d, "nets", "pop"),
                              curve_of(fit, d, "pop", "nets")) |>
    mutate(group = name)
  points[[name]] <- bind_rows(binned(d, "nets"), binned(d, "pop")) |>
    mutate(group = name)
}
curves <- bind_rows(curves)
points <- bind_rows(points)
gams <- bind_rows(gams)
write.csv(curves, file.path(output_dir, sprintf("%s_change_curves.csv",
                                                label)),
          row.names = FALSE)
write.csv(gams, file.path(output_dir, sprintf("%s_change_gams.csv", label)),
          row.names = FALSE)
options(width = 200)
print(gams, n = Inf)

# every pair, against each covariate
raw_points <- bind_rows(lapply(groups, function(name) {
  d <- pairs[group_rows(name), ]
  bind_rows(tibble(term = "nets", at = d$nets, z = d$z),
            tibble(term = "pop", at = d$pop, z = d$z)) |>
    mutate(group = factor(name, groups))
}))

term_labels <- c(nets = "net use (model's nets column), mean over interval",
                 pop = sprintf("population, encounter transform (%s)", label))
p <- ggplot(mutate(curves, group = factor(group, groups)), aes(at, estimate)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_point(data = raw_points, aes(at, z), size = 0.4, alpha = 0.25,
             colour = "grey30") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3,
              fill = "#0072B2") +
  geom_line(colour = "#0072B2", linewidth = 0.8) +
  facet_grid(group ~ term, scales = "free_x",
             labeller = labeller(term = term_labels)) +
  labs(x = NULL, y = "z = qnorm(PIT) of observed change",
       title = sprintf(paste("Change residuals of %s (mean floor %.3f),",
                             "sampled years >= %d apart within pixels"),
                       basename(file), floor_mean, min_gap),
       subtitle = paste("Positive: mortality rose more than predicted",
                        "(selection over-predicted). Points: every pair")) +
  theme_minimal()
ggsave(sprintf("figures/change_residuals_%s.png", label), p,
       width = 8, height = 2 + 1.6 * length(groups), dpi = 150, bg = "white")
