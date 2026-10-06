# Is West African selection attributed to nets that belongs to cotton (#37)?
# On a saved fit: each type's fitted selection effect of the cotton column
# (cotton area x the agricultural trend g_ag), and, for pyrethroids by
# region,
#   1. level residuals (empirical minus predicted logit of mortality) on
#      cumulative cotton exposure (the cotton column summed from the baseline
#      year), with country fixed effects, alone and with cumulative net use;
#   2. the change misfit of within-pixel pairs (observed minus predicted fall
#      in logit mortality per year; R/change_misfit_plot.R) on mean cotton
#      exposure over the pair's span, weighted by precision;
#   3. the correlation of cotton and net exposure at the bioassays.
#
#   Rscript R/cotton_misfit.R <fit RData> <label> [gap]
#
# Writes outputs/hinge_screen/<label>_cotton_misfit.csv.

arguments <- commandArgs(trailingOnly = TRUE)
fit_file <- arguments[1]
label <- arguments[2]
gap <- if (length(arguments) >= 3) as.integer(arguments[3]) else 5L

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
effect <- apply(par$effect_type, c(2, 3), mean)
dimnames(effect) <- list(colnames(f$x_cell_years), f$types)
cotton_column <- grep("^cotton", colnames(f$x_cell_years), value = TRUE)
cat("fitted selection effect exp(beta), by type:\n")
print(signif(effect[c("nets", cotton_column), ], 3))

logit <- dynamical_logit(par, select(df, cell_id, country_id, type_id,
                                     year_id),
                         df, f$x_cell_years, f$cell_years_index)
rho <- colMeans(par$rho_types)[df$type_id]
n_years <- max(f$cell_years_index$year_id)
x_row <- matrix(NA_integer_, max(f$cell_years_index$cell_id), n_years)
x_row[cbind(f$cell_years_index$cell_id, f$cell_years_index$year_id)] <-
  seq_len(nrow(f$cell_years_index))
column_matrix <- function(column) {
  matrix(f$x_cell_years[x_row, column], nrow(x_row))
}
nets <- column_matrix("nets")
cotton <- column_matrix(cotton_column)
cumulative <- function(m) t(apply(m, 1, cumsum))
cum_nets <- cumulative(nets)
cum_cotton <- cumulative(cotton)

n <- df$mosquito_number
d <- tibble(
  type = df$insecticide_type, class = df$insecticide_class,
  country = df$country_name, region = df$region,
  residual = qlogis((df$died + 0.5) / (n + 1)) - colMeans(logit),
  cum_nets = cum_nets[cbind(df$cell_id, df$year_id)],
  cum_cotton = cum_cotton[cbind(df$cell_id, df$year_id)],
  w = 1 / ((1 / (df$died + 0.5) + 1 / (n - df$died + 0.5)) *
             (1 + (n - 1) * rho))) |>
  filter(class == "Pyrethroids")
rm(logit)

slope <- function(formula, data, term) {
  fit <- lm(formula, data = data, weights = w)
  s <- summary(fit)$coefficients
  c(slope = s[term, "Estimate"], se = s[term, "Std. Error"])
}
regions <- c("all", "Western Africa", "Middle Africa", "Eastern Africa")
level <- bind_rows(lapply(regions, function(r) {
  g <- d[r == "all" | d$region == r, ]
  alone <- slope(residual ~ cum_cotton + factor(country) + factor(type), g,
                 "cum_cotton")
  both <- slope(residual ~ cum_cotton + cum_nets + factor(country) +
                  factor(type), g, "cum_cotton")
  tibble(test = "level residual ~ cumulative cotton", region = r,
         n = nrow(g), slope = alone[["slope"]], se = alone[["se"]],
         slope_with_nets = both[["slope"]], se_with_nets = both[["se"]],
         cor_cotton_nets = cor(g$cum_cotton, g$cum_nets),
         share_with_cotton = mean(g$cum_cotton > 0))
}))

stem <- if (gap > 1) sprintf("%s_gap%d", label, gap) else label
pairs <- read.csv(file.path("outputs/hinge_screen",
                            sprintf("%s_change_pairs.csv", stem)))
span <- pairs$year_2 - pairs$year_1
elogit_var <- function(died, n, rho) {
  (1 / (died + 0.5) + 1 / (n - died + 0.5)) * (1 + (n - 1) * rho)
}
pairs <- pairs |>
  mutate(w = span ^ 2 / (elogit_var(died_1, n_1, rho) +
                           elogit_var(died_2, n_2, rho)),
         misfit = (-change_obs + change_pred) / span,
         cotton = mapply(function(cell, y1, y2) mean(cotton[cell, (y1 + 1):y2]),
                         cell_id, year_id_1, year_id_2),
         region = df$region[match(cell_id, df$cell_id)]) |>
  filter(insecticide_class == "Pyrethroids")
change <- bind_rows(lapply(regions, function(r) {
  g <- pairs[r == "all" | pairs$region == r, ]
  s <- slope(misfit ~ cotton + nets, g, "cotton")
  tibble(test = sprintf("change misfit ~ cotton + nets (%d-year pairs)", gap),
         region = r, n = nrow(g), slope = s[["slope"]], se = s[["se"]],
         cor_cotton_nets = cor(g$cotton, g$nets),
         share_with_cotton = mean(g$cotton > 0))
}))
results <- bind_rows(level, change)
write.csv(results, file.path("outputs/hinge_screen",
                             sprintf("%s_cotton_misfit.csv", label)),
          row.names = FALSE)
options(width = 220)
print(as.data.frame(mutate(results, across(where(is.numeric),
                                           ~ signif(.x, 3)))),
      row.names = FALSE)
