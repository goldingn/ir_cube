# Between-place against within-place information on the net-use effect (#37):
# whether a fit's net coefficient is set by differences in level between
# places rather than by change within places.
#
#   Rscript R/between_within_nets.R <fit RData> <label>
#
# Per bioassay: the empirical logit of mortality, the fit's posterior mean
# predicted logit, and the pixel's cumulative net use (the model's nets column
# summed from the baseline year to the bioassay's year). Weighted (empirical-
# logit precision, inflated by the type's rho) regressions by insecticide
# group:
#   1. level residual (empirical minus predicted logit) on cumulative net use,
#      with country fixed effects: does the fit get levels right across
#      places with different net histories, within a country?
#   2. empirical logit on cumulative net use, with country x year fixed
#      effects: the between-site slope at a given time, to set against the
#      fit's coefficient (selection per unit of net use per year, so a slope
#      of -exp(beta) per unit of cumulative net use if levels were driven only
#      by nets) and the within-pixel change (R/change_misfit_plot.R).
#
# Writes outputs/hinge_screen/<label>_between_within.csv.

arguments <- commandArgs(trailingOnly = TRUE)
fit_file <- arguments[1]
label <- arguments[2]

source("R/packages.R")
source("R/functions.R")
source("R/dynamical_model.R")
source("R/dynamical_predictions.R")

f <- new.env()
load(fit_file, envir = f)
df <- f$df
fold <- list(draws = f$draws, options = f$model_options,
             x_cells_init = f$x_cells_init)
n_total <- sum(vapply(fold$draws, nrow, integer(1)))
par <- dynamical_parameter_draws(fold, f$classes_index, f$types,
                                 draw_index = round(seq(1, n_total,
                                                        length.out = 200)))
logit <- dynamical_logit(par, select(df, cell_id, country_id, type_id,
                                     year_id),
                         df, f$x_cell_years, f$cell_years_index)
effect <- apply(par$effect_type, c(2, 3), mean)
dimnames(effect) <- list(colnames(f$x_cell_years), f$types)
rho <- colMeans(par$rho_types)[df$type_id]

# cumulative net use at each bioassay's pixel and year
n_years <- max(f$cell_years_index$year_id)
x_row <- matrix(NA_integer_, max(f$cell_years_index$cell_id), n_years)
x_row[cbind(f$cell_years_index$cell_id, f$cell_years_index$year_id)] <-
  seq_len(nrow(f$cell_years_index))
nets <- matrix(f$x_cell_years[x_row, "nets"], nrow(x_row))
cum_nets <- t(apply(nets, 1, cumsum))

n <- df$mosquito_number
d <- tibble(
  type = df$insecticide_type, class = df$insecticide_class,
  country = df$country_name, year = df$year_start,
  empirical = qlogis((df$died + 0.5) / (n + 1)),
  predicted = colMeans(logit),
  cum_nets = cum_nets[cbind(df$cell_id, df$year_id)],
  nets = nets[cbind(df$cell_id, df$year_id)],
  w = 1 / ((1 / (df$died + 0.5) + 1 / (n - df$died + 0.5)) *
             (1 + (n - 1) * rho))) |>
  mutate(residual = empirical - predicted,
         group = case_when(class == "Pyrethroids" ~ "Pyrethroids",
                           type == "DDT" ~ "DDT",
                           TRUE ~ class))
rm(logit)

groups <- c("Pyrethroids", "DDT", "Carbamates", "Organophosphates",
            sort(unique(d$type[d$class == "Pyrethroids"])))
rows <- function(name) if (name %in% d$group) d$group == name else
  d$type == name
slope <- function(formula, data) {
  fit <- lm(formula, data = data, weights = w)
  s <- summary(fit)$coefficients["cum_nets", ]
  c(slope = s[["Estimate"]], se = s[["Std. Error"]])
}
results <- bind_rows(lapply(groups, function(name) {
  g <- d[rows(name), ]
  fit_effect <- if (name %in% f$types) effect["nets", name] else
    mean(effect["nets", unique(g$type)])
  # with insecticide type effects when the group has several types
  types <- if (n_distinct(g$type) > 1) " + factor(type)" else ""
  level <- slope(as.formula(paste("residual ~ cum_nets + factor(country)",
                                  types)), g)
  between <- slope(as.formula(paste(
    "empirical ~ cum_nets + factor(paste(country, year))", types)), g)
  tibble(group = name, assays = nrow(g),
         fit_net_effect = fit_effect,
         level_residual_slope = level[["slope"]],
         level_residual_se = level[["se"]],
         between_site_slope = between[["slope"]],
         between_site_se = between[["se"]])
}))
write.csv(results, file.path("outputs/hinge_screen",
                             sprintf("%s_between_within.csv", label)),
          row.names = FALSE)
options(width = 200)
print(as.data.frame(mutate(results, across(where(is.numeric), ~ signif(.x, 3)))),
      row.names = FALSE)
