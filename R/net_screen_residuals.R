# Pre-screen of the net-use transform (#37), step 2: GAMs of the residuals of
# a full fit against net use and cumulative net use, per pixel-year, by
# insecticide group, from outputs/net_screen/extract.rds
# (R/net_screen_extract.R).
#
#   Rscript R/net_screen_residuals.R
#
# Residuals are taken on the logit scale against the fit's prediction: a
# quasibinomial GAM of observed mortality with the logit of the posterior mean
# prediction as an offset, so a smooth is the logit correction the fit would
# need. Bioassay-level noise enters through the records and the weights:
#   - replicate assays at the same pixel, year and type share noise beyond
#     rho (the two-stage model's u), so they are pooled into one record, the
#     unit treated as independent;
#   - each record is weighted by its effective number of mosquitoes,
#     n^2 / sum_i n_i (1 + (n_i - 1) rho), with the fit's per-type rho
#     (beta-binomial overdispersion within an assay), and the quasibinomial
#     dispersion takes up the rest.
# A second version adds an iid random effect per pixel, so that the smooths
# are informed by changes within pixels only, not by differences between
# places (country means mix places tested early and late). Beta-binomial
# Pearson residuals per record, as in R/floor_modes.R, are fitted the same way
# (gaussian, weighted by the number of assays pooled) as a check.
#
# Then, on the same offset model without random effects, two saturating
# families for the cumulative net selection, each with the linear one nested:
#   current use  the selection of a year is f(n) = n / (n + h); the fit's
#                cumulative selection e C (C = sum of n) is wrong by
#                a C + b sum f(n_t)
#   cumulative   the cumulative selection is F(C) = c (1 - exp(-C / c)),
#                linear for C << c and saturating at c; the fit is wrong by
#                a C + b F(C)
# profiled over h and c by the binomial deviance with the same weights.
#
# Writes outputs/net_screen/residual_gams.csv, residual_gam_curves.csv,
# residual_saturating_profile.csv and figures/net_screen_residual_gams.png.

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mgcv)
})
output_dir <- "outputs/net_screen"
x <- readRDS(file.path(output_dir, "extract.rds"))

# the groups: pyrethroids pooled and by type, DDT, the carbamate (bendiocarb)
# and the organophosphates pooled
group_of <- function(class, type) {
  case_when(class == "Pyrethroids" ~ "Pyrethroids",
            type == "DDT" ~ "DDT",
            TRUE ~ class)
}
assay <- x$assay |>
  mutate(group = group_of(insecticide_class, insecticide_type),
         cell_f = factor(cell_id))
groups <- c(list(Pyrethroids = "Pyrethroids", DDT = "DDT",
                 Carbamates = "Carbamates",
                 Organophosphates = "Organophosphates"),
            setNames(as.list(sort(unique(assay$insecticide_type[
              assay$insecticide_class == "Pyrethroids"]))),
              sort(unique(assay$insecticide_type[
                assay$insecticide_class == "Pyrethroids"]))))
group_rows <- function(name) {
  if (name %in% assay$group) assay$group == name else
    assay$insecticide_type == name
}

# The logit of a probability kept away from 0 and 1, for the offset: the
# posterior mean prediction can round to 1 in double precision at the most
# susceptible assays, and to the floor's tiny value at the most resistant
clamped_logit <- function(p, eps = 1e-6) {
  qlogis(pmin(pmax(p, eps), 1 - eps))
}

# the data of one group under one mode, pooled to one record per pixel, year
# and type (assays there share a prediction)
group_data <- function(name, mode) {
  m <- x$modes[[mode]]
  assay[group_rows(name), ] |>
    mutate(pred = m$pred[row], rho = m$rho[row]) |>
    group_by(cell_id, cell_f, year_id, type_id, nets, cum_nets, pred) |>
    summarise(n_eff_factor = sum(n * (1 + (n - 1) * rho)),
              rho = first(rho), assays = n(), died = sum(died), n = sum(n),
              .groups = "drop") |>
    mutate(offset_logit = clamped_logit(pred),
           y = died / n,
           n_eff = n ^ 2 / n_eff_factor,
           pearson = (died - n * pred) / sqrt(pred * (1 - pred) *
                                                n_eff_factor))
}

# the models: the formula's right-hand side and its family
gam_specs <- list(
  logit = list(rhs = "s(nets, k = 6) + s(cum_nets, k = 6)",
               family = quasibinomial()),
  logit_within_pixel = list(rhs = "s(nets, k = 6) + s(cum_nets, k = 6) +
                                   s(cell_f, bs = 're')",
                            family = quasibinomial()),
  pearson = list(rhs = "s(nets, k = 6) + s(cum_nets, k = 6)",
                 family = gaussian())
)

fit_gam <- function(d, spec) {
  if (identical(spec$family$family, "gaussian")) {
    bam(as.formula(paste("pearson ~", spec$rhs)), data = d,
        weights = assays, discrete = TRUE)
  } else {
    bam(as.formula(paste("y ~ offset(offset_logit) +", spec$rhs)), data = d,
        weights = n_eff, family = spec$family, discrete = TRUE)
  }
}

# a smooth's value at `at` less its value at `reference`, with a 95% interval
# from the lpmatrix
smooth_contrast <- function(fit, d, term, at, reference) {
  new <- d[rep(1, length(at) + 1), ]
  new[[term]] <- c(reference, at)
  other <- setdiff(c("nets", "cum_nets"), term)
  new[[other]] <- median(d[[other]])
  X <- predict(fit, new, type = "lpmatrix", exclude = "s(cell_f)")
  cols <- grep(sprintf("^s\\(%s\\)", term), colnames(X))
  D <- sweep(X[-1, cols, drop = FALSE], 2, X[1, cols], "-")
  est <- c(D %*% coef(fit)[cols])
  se <- sqrt(rowSums((D %*% vcov(fit)[cols, cols]) * D))
  tibble(term = term, at = at, reference = reference, estimate = est,
         lower = est - 1.96 * se, upper = est + 1.96 * se)
}

summary_rows <- list()
curve_rows <- list()
# every group and model in the low-floor mode (the main fit's); the pooled
# groups' logit model in the high-floor mode
runs <- bind_rows(
  expand.grid(mode = "low", group = names(groups), model = names(gam_specs),
              stringsAsFactors = FALSE),
  expand.grid(mode = "high", group = c("Pyrethroids", "DDT", "Carbamates",
                                       "Organophosphates"),
              model = "logit", stringsAsFactors = FALSE))
runs <- runs[runs$mode %in% names(x$modes), ]
for (r in seq_len(nrow(runs))) {
  mode <- runs$mode[r]
  name <- runs$group[r]
  model <- runs$model[r]
  d <- group_data(name, mode)
  fit <- fit_gam(d, gam_specs[[model]])
  fit_fixed <- if (model == "logit_within_pixel") {
    fit_gam(d, gam_specs$logit)
  } else {
    fit
  }
  s <- summary(fit)$s.table
  contrasts <- bind_rows(
    smooth_contrast(fit, d, "nets", c(0.2, 0.4, 0.6), 0),
    smooth_contrast(fit, d, "cum_nets", c(2, 4, 6, 8), 0))
  summary_rows[[length(summary_rows) + 1]] <- tibble(
    mode = mode, group = name, model = model, n_assays = sum(d$assays),
    n_records = nrow(d),
    n_pixels = n_distinct(d$cell_id),
    max_cum_nets = max(d$cum_nets),
    edf_nets = s["s(nets)", "edf"], p_nets = s["s(nets)", "p-value"],
    edf_cum_nets = s["s(cum_nets)", "edf"],
    p_cum_nets = s["s(cum_nets)", "p-value"],
    dev_expl_fixed = summary(fit_fixed)$dev.expl,
    dev_expl = summary(fit)$dev.expl) |>
    bind_cols(contrasts |>
                transmute(name = sprintf("%s_%g", term, at),
                          value = sprintf("%+.2f [%+.2f, %+.2f]",
                                          estimate, lower, upper)) |>
                pivot_wider())
  if (model != "pearson") {
    grid <- bind_rows(
      smooth_contrast(fit, d, "nets",
                      seq(0, quantile(d$nets, 0.99), length.out = 40)[-1],
                      0),
      smooth_contrast(fit, d, "cum_nets",
                      seq(0, quantile(d$cum_nets, 0.99),
                          length.out = 40)[-1], 0))
    curve_rows[[length(curve_rows) + 1]] <- grid |>
      mutate(mode = mode, group = name, model = model)
  }
  cat(mode, name, model, "done\n")
  rm(fit, fit_fixed)
  invisible(gc())
}
gams <- bind_rows(summary_rows)
curves <- bind_rows(curve_rows)
write.csv(gams, file.path(output_dir, "residual_gams.csv"), row.names = FALSE)
write.csv(curves, file.path(output_dir, "residual_gam_curves.csv"),
          row.names = FALSE)
options(width = 250)
print(gams |> filter(model == "logit") |>
        select(mode, group, n_records, edf_nets, edf_cum_nets, dev_expl_fixed,
               nets_0.6, cum_nets_4, cum_nets_8), n = Inf)
print(gams |> filter(model == "logit_within_pixel") |>
        select(mode, group, edf_nets, edf_cum_nets, nets_0.6, cum_nets_4,
               cum_nets_8), n = Inf)
print(gams |> filter(model == "pearson") |>
        select(mode, group, edf_nets, edf_cum_nets, dev_expl_fixed,
               nets_0.6, cum_nets_4, cum_nets_8), n = Inf)


# saturating families --------------------------------------------------------

# the running sum over years 1..t of f(n_s) at each assay's pixel, for the
# current-use family; nets_by_cell is cells x years
nets_by_cell <- matrix(NA_real_, max(x$cell_years_index$cell_id),
                       max(x$cell_years_index$year_id))
nets_by_cell[cbind(x$cell_years_index$cell_id, x$cell_years_index$year_id)] <-
  x$x_cell_years[, "nets"]
cumulative_of <- function(f, d) {
  cum <- t(apply(f(nets_by_cell), 1, cumsum))
  cum[cbind(d$cell_id, d$year_id)]
}
saturating_cumulative <- function(cum, c) c * -expm1(-cum / c)

# binomial deviance with the effective sample sizes, and its coefficients
offset_deviance <- function(d, z) {
  fit <- glm(y ~ offset(offset_logit) + cum_nets + z, data = d,
             weights = n_eff, family = quasibinomial())
  c(deviance = deviance(fit), a = coef(fit)[["cum_nets"]],
    b = coef(fit)[["z"]])
}
profile_rows <- list()
for (name in c("Pyrethroids", "DDT", "Carbamates", "Organophosphates")) {
  d <- group_data(name, "low")
  linear <- glm(y ~ offset(offset_logit) + cum_nets, data = d,
                weights = n_eff, family = quasibinomial())
  null <- glm(y ~ offset(offset_logit), data = d, weights = n_eff,
              family = quasibinomial())
  for (h in c(0.025, 0.05, 0.1, 0.2, 0.4, 0.8)) {
    z <- cumulative_of(function(n) n / (n + h), d)
    profile_rows[[length(profile_rows) + 1]] <- tibble(
      group = name, family = "current use n / (n + h)", parameter = h,
      !!!as.list(offset_deviance(d, z)))
  }
  for (c in c(1, 2, 3, 4, 6, 8, 12)) {
    z <- saturating_cumulative(d$cum_nets, c)
    profile_rows[[length(profile_rows) + 1]] <- tibble(
      group = name, family = "cumulative c (1 - exp(-C / c))", parameter = c,
      !!!as.list(offset_deviance(d, z)))
  }
  profile_rows[[length(profile_rows) + 1]] <- tibble(
    group = name, family = c("none", "linear in C"), parameter = NA,
    deviance = c(deviance(null), deviance(linear)),
    a = c(NA, coef(linear)[["cum_nets"]]), b = NA)
}
profile <- bind_rows(profile_rows) |>
  group_by(group) |>
  mutate(deviance_vs_none = deviance - deviance[family == "none"]) |>
  ungroup()
write.csv(profile, file.path(output_dir, "residual_saturating_profile.csv"),
          row.names = FALSE)
print(profile, n = Inf)


# figure ---------------------------------------------------------------------

model_labels <- c(logit = "all variation",
                  logit_within_pixel = "within pixels")
group_order <- c("Pyrethroids", "DDT", "Carbamates", "Organophosphates",
                 setdiff(names(groups), c("Pyrethroids", "DDT", "Carbamates",
                                          "Organophosphates")))
term_labels <- c(nets = "net use (pyrethroid-weighted)",
                 cum_nets = "cumulative net use since 1995")
plot_data <- curves |>
  filter(mode == "low") |>
  mutate(group = factor(group, group_order),
         model = factor(model_labels[model], model_labels),
         term = factor(term_labels[term], term_labels))
p <- ggplot(plot_data, aes(at, estimate, colour = model, fill = model)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.6) +
  facet_grid(group ~ term, scales = "free_x") +
  scale_colour_manual(values = c("#0072B2", "#D55E00"), name = NULL,
                      aesthetics = c("colour", "fill")) +
  coord_cartesian(ylim = c(-1.5, 2.5)) +
  labs(x = NULL,
       y = "logit correction to predicted mortality (vs no net use)",
       title = "Residuals of the round 2 main fit (low-floor mode)",
       subtitle = paste("Positive: observed mortality above the fit,",
                        "i.e. resistance over-predicted")) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "top", panel.grid.minor = element_blank())
ggsave("figures/net_screen_residual_gams.png", p, width = 7, height = 11,
       dpi = 150, bg = "white")
cat("done\n")
