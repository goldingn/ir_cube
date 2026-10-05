# Hinge screen (#37): GAMs of the residuals of a saved full fit against net
# use and against the transformed population covariate, to choose hinge knots
# (selection_design(hinges = ...), columns min(x, k)).
#
#   Rscript R/hinge_screen.R [fit file] [label]
#
# Defaults: outputs/refit/fits/r2_sens_dhalf200.RData (round 2 dynamical fit,
# population encounter transform with d_half = 200), label "dhalf200".
#
# The fit's posterior mean prediction at each bioassay is taken as an offset
# in a quasibinomial GAM of observed mortality, so the smooths are the logit
# correction the fit needs. Assays at the same pixel, year and type share a
# prediction and are pooled into one record, weighted by its effective number
# of mosquitoes under the fit's rho, as in R/net_screen_residuals.R. Population
# is the model's population column divided by its time trend g_dom(t), i.e.
# the encounter-transformed density on the 0-1 scale that a hinge knot is set
# on. The curves are drawn on the absolute scale (intercept included, the
# other covariate at its median), not anchored at 0, with the pooled
# residuals binned by covariate decile as points.
#
# Writes outputs/hinge_screen/<label>_curves.csv, <label>_gams.csv and
# figures/hinge_screen_<label>.png.

arguments <- commandArgs(trailingOnly = TRUE)
file <- if (length(arguments) >= 1) arguments[1] else
  "outputs/refit/fits/r2_sens_dhalf200.RData"
label <- if (length(arguments) >= 2) arguments[2] else "dhalf200"
n_draws_keep <- 400
output_dir <- "outputs/hinge_screen"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

source("R/packages.R")
source("R/functions.R")
source("R/dynamical_model.R")
source("R/dynamical_predictions.R")
suppressMessages(library(mgcv))

f <- new.env()
load(file, envir = f)
df <- f$df

# posterior mean prediction at each assay, from the fit's kept draws
fold <- list(draws = f$draws, options = f$model_options,
             x_cells_init = f$x_cells_init)
n_total <- sum(vapply(fold$draws, nrow, integer(1)))
par <- dynamical_parameter_draws(fold, f$classes_index, f$types,
                                 draw_index = round(seq(1, n_total,
                                                        length.out =
                                                          n_draws_keep)))
logit <- dynamical_logit(par, select(df, cell_id, country_id, type_id,
                                     year_id),
                         df, f$x_cell_years, f$cell_years_index)
pred <- colMeans(plogis(logit))
rho <- colMeans(par$rho_types)[df$type_id]
floor_mean <- mean(par$mortality_floor)
rm(logit)
invisible(gc())

# covariates at each assay's cell-year
x_row <- match(paste(df$cell_id, df$year_id),
               paste(f$cell_years_index$cell_id, f$cell_years_index$year_id))
pop_column <- grep("^pop_[a-z]+:g_dom$", colnames(f$x_cell_years),
                   value = TRUE)
stopifnot(length(pop_column) == 1, !anyNA(x_row))
# g_dom(t) of the "linear_0_1" trend, 0 in 1995 and 1 in 2025
year <- f$baseline_year + df$year_id - 1
g_dom <- pmax((year - 1995) / 30, 0)

# The logit of a probability kept away from 0 and 1, for the offset
clamped_logit <- function(p, eps = 1e-6) {
  qlogis(pmin(pmax(p, eps), 1 - eps))
}

group_of <- function(class, type) {
  case_when(class == "Pyrethroids" ~ "Pyrethroids",
            type == "DDT" ~ "DDT",
            TRUE ~ class)
}

assay <- df |>
  transmute(cell_id, year_id, type_id, insecticide_class, insecticide_type,
            died, n = mosquito_number,
            nets = f$x_cell_years[x_row, "nets"],
            pop = f$x_cell_years[x_row, pop_column] / g_dom,
            pred = pred, rho = rho) |>
  filter(is.finite(pop)) |>
  mutate(group = group_of(insecticide_class, insecticide_type))

pyrethroid_types <- sort(unique(assay$insecticide_type[
  assay$insecticide_class == "Pyrethroids"]))
groups <- c("Pyrethroids", "DDT", "Carbamates", "Organophosphates",
            pyrethroid_types)
group_rows <- function(name) {
  if (name %in% assay$group) assay$group == name else
    assay$insecticide_type == name
}

group_data <- function(name) {
  assay[group_rows(name), ] |>
    group_by(cell_id, year_id, type_id, nets, pop, pred) |>
    summarise(n_eff_factor = sum(n * (1 + (n - 1) * rho)),
              died = sum(died), n = sum(n), .groups = "drop") |>
    mutate(cell_f = factor(cell_id),
           offset_logit = clamped_logit(pred),
           y = died / n,
           n_eff = n ^ 2 / n_eff_factor,
           # empirical logit residual, for the binned points
           resid = qlogis((died + 0.5) / (n + 1)) - offset_logit)
}

specs <- c(all_variation = "s(nets, k = 8) + s(pop, k = 8)",
           within_pixel = "s(nets, k = 8) + s(pop, k = 8) +
                           s(cell_f, bs = 're')")

# the fitted correction along one covariate, the other at its median, with
# the intercept and without the pixel effects
curve_of <- function(fit, d, term, other) {
  at <- seq(quantile(d[[term]], 0.005), quantile(d[[term]], 0.995),
            length.out = 60)
  new <- d[rep(1, length(at)), ]
  new[[term]] <- at
  new[[other]] <- median(d[[other]])
  new$offset_logit <- 0
  p <- predict(fit, new, se.fit = TRUE, exclude = "s(cell_f)")
  tibble(term = term, at = at, estimate = p$fit,
         lower = p$fit - 1.96 * p$se.fit, upper = p$fit + 1.96 * p$se.fit)
}

binned <- function(d, term) {
  d |>
    mutate(bin = ntile(.data[[term]], 10)) |>
    group_by(bin) |>
    summarise(at = weighted.mean(.data[[term]], n_eff),
              estimate = weighted.mean(resid, n_eff), records = n(),
              .groups = "drop") |>
    mutate(term = term)
}

# the GAMs take about 20 minutes; saved curves are reused unless
# IR_CUBE_REFIT_GAMS is set
curves_file <- file.path(output_dir, sprintf("%s_curves.csv", label))
refit_gams <- !file.exists(curves_file) ||
  nzchar(Sys.getenv("IR_CUBE_REFIT_GAMS"))
curves <- list()
points <- list()
gams <- list()
for (name in groups) {
  d <- group_data(name)
  for (model in if (refit_gams) names(specs) else character(0)) {
    fit <- bam(as.formula(paste("y ~ offset(offset_logit) +", specs[[model]])),
               data = d, weights = n_eff, family = quasibinomial(),
               discrete = TRUE)
    s <- summary(fit)$s.table
    gams[[length(gams) + 1]] <- tibble(
      group = name, model = model, records = nrow(d),
      edf_nets = s["s(nets)", "edf"], p_nets = s["s(nets)", "p-value"],
      edf_pop = s["s(pop)", "edf"], p_pop = s["s(pop)", "p-value"],
      dev_expl = summary(fit)$dev.expl)
    curves[[length(curves) + 1]] <- bind_rows(
      curve_of(fit, d, "nets", "pop"),
      curve_of(fit, d, "pop", "nets")) |>
      mutate(group = name, model = model)
    rm(fit)
  }
  points[[length(points) + 1]] <- bind_rows(binned(d, "nets"),
                                            binned(d, "pop")) |>
    mutate(group = name)
  cat(name, "done\n")
}
points <- bind_rows(points)
write.csv(points, file.path(output_dir, sprintf("%s_points.csv", label)),
          row.names = FALSE)
if (refit_gams) {
  curves <- bind_rows(curves)
  gams <- bind_rows(gams)
  write.csv(curves, curves_file, row.names = FALSE)
  write.csv(gams, file.path(output_dir, sprintf("%s_gams.csv", label)),
            row.names = FALSE)
} else {
  curves <- read.csv(curves_file)
  gams <- read.csv(file.path(output_dir, sprintf("%s_gams.csv", label)))
}
options(width = 200)
print(as_tibble(gams), n = Inf)

term_labels <- c(nets = "net use (model's nets column)",
                 pop = sprintf("population, encounter transform (%s)", label))
group_levels <- groups
p <- ggplot(mutate(curves, group = factor(group, group_levels)),
            aes(at, estimate)) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = 0.15) +
  geom_line(aes(colour = model)) +
  geom_point(data = mutate(points, group = factor(group, group_levels)),
             size = 0.9, colour = "grey20") +
  facet_grid(group ~ term, scales = "free",
             labeller = labeller(term = term_labels)) +
  scale_colour_manual(values = c(all_variation = "#0072B2",
                                 within_pixel = "#D55E00"),
                      aesthetics = c("colour", "fill")) +
  labs(x = NULL, y = "logit correction to predicted mortality",
       colour = NULL, fill = NULL,
       title = sprintf("Residuals of %s (mean floor %.3f)", basename(file),
                       floor_mean),
       subtitle = paste("Positive: observed mortality above the fit.",
                        "Points: pooled residuals by covariate decile")) +
  theme_minimal() +
  theme(legend.position = "top")
ggsave(sprintf("figures/hinge_screen_%s.png", label), p,
       width = 8, height = 2 + 1.6 * length(groups), dpi = 150, bg = "white")
