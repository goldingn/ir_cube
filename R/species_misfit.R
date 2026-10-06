# Does the species tested explain the net-use misfit (#37)? On a saved fit,
# for the bioassays identified to species, An. arabiensis against An. gambiae
# s.s. and An. coluzzii:
#   1. level: the slope of the level residual (empirical minus the fit's
#      predicted logit of mortality) on cumulative net use, with country fixed
#      effects, by species group, overall and by region;
#   2. change: within-pixel pairs of sampled years (at least `min_gap` apart,
#      a chain as in R/change_residual_screen.R) for each species group
#      separately, and the slope of the misfit (observed minus predicted fall
#      in logit mortality per year; negative: the fit over-predicts selection)
#      on mean net use over the pair's span, weighted by precision.
#
#   Rscript R/species_misfit.R <fit RData> <label> [min gap]
#
# Writes outputs/hinge_screen/<label>_species_misfit.csv.

arguments <- commandArgs(trailingOnly = TRUE)
fit_file <- arguments[1]
label <- arguments[2]
min_gap <- if (length(arguments) >= 3) as.integer(arguments[3]) else 3L

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
rho <- colMeans(par$rho_types)[df$type_id]
n_years <- max(f$cell_years_index$year_id)
x_row <- matrix(NA_integer_, max(f$cell_years_index$cell_id), n_years)
x_row[cbind(f$cell_years_index$cell_id, f$cell_years_index$year_id)] <-
  seq_len(nrow(f$cell_years_index))
nets <- matrix(f$x_cell_years[x_row, "nets"], nrow(x_row))
cum_nets <- t(apply(nets, 1, cumsum))

n <- df$mosquito_number
d <- tibble(
  row = seq_len(nrow(df)), cell_id = df$cell_id, year_id = df$year_id,
  type = df$insecticide_type, class = df$insecticide_class,
  country = df$country_name, region = df$region,
  species = case_when(
    df$species == "Anopheles arabiensis" ~ "arabiensis",
    df$species %in% c("Anopheles gambiae", "Anopheles gambiae s.s.",
                      "Anopheles gambiae ss", "Anopheles coluzzii") ~
      "gambiae/coluzzii",
    TRUE ~ NA_character_),
  died = df$died, n = n,
  empirical = qlogis((df$died + 0.5) / (n + 1)),
  predicted = colMeans(logit),
  cum_nets = cum_nets[cbind(df$cell_id, df$year_id)],
  rho = rho) |>
  mutate(residual = empirical - predicted,
         w = 1 / ((1 / (died + 0.5) + 1 / (n - died + 0.5)) *
                    (1 + (n - 1) * rho))) |>
  filter(!is.na(species), class == "Pyrethroids")
rm(logit)

weighted_slope <- function(formula, data, term) {
  if (nrow(data) < 30) return(c(slope = NA, se = NA, n = nrow(data)))
  fit <- lm(formula, data = data, weights = w)
  s <- summary(fit)$coefficients
  if (!term %in% rownames(s)) return(c(slope = NA, se = NA, n = nrow(data)))
  c(slope = s[term, "Estimate"], se = s[term, "Std. Error"], n = nrow(data))
}
fixed <- function(data, base) {
  extra <- c(if (n_distinct(data$country) > 1) "factor(country)",
             if (n_distinct(data$type) > 1) "factor(type)")
  as.formula(paste(base, paste(c("", extra), collapse = " + ")))
}

# 1. level residuals on cumulative net use
level <- bind_rows(lapply(c("all", "Eastern Africa", "Western Africa",
                            "Middle Africa", "Northern Africa"), function(r) {
  bind_rows(lapply(c("arabiensis", "gambiae/coluzzii"), function(s) {
    g <- d[d$species == s & (r == "all" | d$region == r), ]
    v <- weighted_slope(fixed(g, "residual ~ cum_nets"), g, "cum_nets")
    tibble(test = "level residual ~ cumulative net use", region = r,
           species = s, slope = v[["slope"]], se = v[["se"]], n = v[["n"]])
  }))
}))

# 2. within-pixel changes by species group
pooled <- d |>
  group_by(species, cell_id, type, year_id, country, region) |>
  summarise(died = sum(died), n = sum(n),
            var_factor = sum(n * (1 + (n - 1) * rho)),
            predicted = mean(predicted), .groups = "drop") |>
  mutate(empirical = qlogis((died + 0.5) / (n + 1)),
         variance = (1 / (died + 0.5) + 1 / (n - died + 0.5)) *
           var_factor / n)
gap_chain <- function(years) {
  keep <- 1L
  for (i in seq_along(years)[-1]) {
    if (years[i] - years[keep[length(keep)]] >= min_gap) keep <- c(keep, i)
  }
  keep
}
pairs <- pooled |>
  arrange(species, cell_id, type, year_id) |>
  group_by(species, cell_id, type) |>
  filter(row_number() %in% gap_chain(year_id)) |>
  mutate(year_1 = lag(year_id), emp_1 = lag(empirical),
         pred_1 = lag(predicted), var_1 = lag(variance)) |>
  ungroup() |>
  filter(!is.na(year_1)) |>
  mutate(span = year_id - year_1,
         misfit = ((emp_1 - empirical) - (pred_1 - predicted)) / span,
         w = span ^ 2 / (var_1 + variance))
pairs$nets <- mapply(function(cell, y1, y2) mean(nets[cell, (y1 + 1):y2]),
                     pairs$cell_id, pairs$year_1, pairs$year_id)
change <- bind_rows(lapply(c("all", "Eastern Africa", "Western Africa",
                             "Northern Africa"), function(r) {
  bind_rows(lapply(c("arabiensis", "gambiae/coluzzii"), function(s) {
    g <- pairs[pairs$species == s & (r == "all" | pairs$region == r), ]
    v <- weighted_slope(misfit ~ nets, g, "nets")
    tibble(test = sprintf("change misfit ~ net use (pairs >= %d years)",
                          min_gap),
           region = r, species = s, slope = v[["slope"]], se = v[["se"]],
           n = v[["n"]],
           mean_misfit = if (nrow(g)) weighted.mean(g$misfit, g$w) else NA)
  }))
}))
results <- bind_rows(level, change)
write.csv(results, file.path("outputs/hinge_screen",
                             sprintf("%s_species_misfit.csv", label)),
          row.names = FALSE)
options(width = 200)
print(as.data.frame(mutate(results, across(where(is.numeric),
                                           ~ signif(.x, 3)))),
      row.names = FALSE)
