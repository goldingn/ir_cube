# What explains a fit's pyrethroid level misfit (#37)? A national, a
# geographic, or a data effect? The response is each pyrethroid bioassay's
# level residual (empirical minus the fit's predicted logit of mortality),
# weighted by precision (empirical-logit variance inflated by the type's
# rho). Reported:
#   1. the weighted mean residual by country;
#   2. the deviance explained by country fixed effects, a smooth over
#      longitude and latitude, a random effect of citation (study), and the
#      candidate variables;
#   3. a multivariable regression on the candidates (data source, species
#      identification, WHO protocol version, sample size, a round sample size
#      of 100, year, net use, cumulative net use, MITN-filled pixel, IRS,
#      population, cotton, insecticide type), with and without country
#      effects, with a citation random effect; and the country effects with
#      and without the candidates.
#
#   Rscript R/misfit_regression.R <fit RData> <label> <filled raster> [from]
#
# With `from`, only bioassays from that year on.
#
# Writes outputs/hinge_screen/<label>_misfit_regression[_from<from>].txt.

arguments <- commandArgs(trailingOnly = TRUE)
fit_file <- arguments[1]
label <- arguments[2]
filled_file <- arguments[3]
from <- if (length(arguments) >= 4) as.integer(arguments[4]) else NA

source("R/packages.R")
source("R/functions.R")
source("R/dynamical_model.R")
source("R/dynamical_predictions.R")
suppressMessages(library(mgcv))

f <- new.env()
load(fit_file, envir = f)
df <- f$df
mask <- rast("data/clean/raster_mask.tif")
fold <- list(draws = f$draws, options = f$model_options,
             x_cells_init = f$x_cells_init)
n_total <- sum(vapply(fold$draws, nrow, integer(1)))
par <- dynamical_parameter_draws(fold, f$classes_index, f$types,
                                 draw_index = round(seq(1, n_total,
                                                        length.out = 200)))
predicted <- colMeans(dynamical_logit(par, select(df, cell_id, country_id,
                                                  type_id, year_id),
                                      df, f$x_cell_years, f$cell_years_index))
rho <- colMeans(par$rho_types)[df$type_id]

n_years <- max(f$cell_years_index$year_id)
x_row <- matrix(NA_integer_, max(f$cell_years_index$cell_id), n_years)
x_row[cbind(f$cell_years_index$cell_id, f$cell_years_index$year_id)] <-
  seq_len(nrow(f$cell_years_index))
column_at <- function(column) {
  m <- matrix(f$x_cell_years[x_row, column], nrow(x_row))
  list(now = m[cbind(df$cell_id, df$year_id)],
       cumulative = t(apply(m, 1, cumsum))[cbind(df$cell_id, df$year_id)])
}
pop_column <- grep("^pop_[a-z]+:g_dom$", colnames(f$x_cell_years),
                   value = TRUE)
cotton_column <- grep("^cotton", colnames(f$x_cell_years), value = TRUE)
net_columns <- column_at("nets")
irs_columns <- column_at("irs")
g_dom <- pmax((df$year_start - 1995) / 30, 0)
filled <- rast(filled_file)
filled_year <- pmin(pmax(df$year_start, 2000), terra::nlyr(filled) + 1999)
filled_value <- terra::extract(filled, f$unique_cells[df$cell_id])
filled_value <- as.matrix(filled_value)[cbind(seq_len(nrow(df)),
                                              filled_year - 1999)]

protocol <- function(test_type) {
  version <- suppressWarnings(as.integer(sub(".*bioassay (\\d{4}).*", "\\1",
                                             test_type)))
  case_when(grepl("bottle", test_type) ~ "bottle",
            !is.na(version) & version < 2013 ~ "WHO pre-2013",
            !is.na(version) & version >= 2013 & version < 2016 ~ "WHO 2013",
            !is.na(version) & version >= 2016 & version < 2030 ~ "WHO 2016+",
            TRUE ~ "WHO, version unknown")
}
n <- df$mosquito_number
d <- tibble(
  residual = qlogis((df$died + 0.5) / (n + 1)) - predicted,
  w = 1 / ((1 / (df$died + 0.5) + 1 / (n - df$died + 0.5)) *
             (1 + (n - 1) * rho)),
  class = df$insecticide_class, type = df$insecticide_type,
  country = df$country_name, citation = factor(df$citation),
  source = sub("_.*", "", df$source) |>
    recode(mapper = "IR Mapper", mtm = "MTM", va = "Vector Atlas"),
  species = case_when(
    df$species == "Anopheles arabiensis" ~ "arabiensis",
    df$species %in% c("Anopheles gambiae", "Anopheles gambiae s.s.",
                      "Anopheles gambiae ss", "Anopheles coluzzii") ~
      "gambiae/coluzzii",
    df$species == "gambiae complex" ~ "complex only",
    TRUE ~ "other"),
  protocol = protocol(df$test_type),
  log_n = log(n), n_100 = n == 100,
  year = df$year_start, nets = net_columns$now,
  cum_nets = net_columns$cumulative,
  irs = irs_columns$now, pop = f$x_cell_years[x_row[cbind(df$cell_id, df$year_id)],
                                      pop_column] / pmax(g_dom, 1e-6),
  cotton = f$x_cell_years[x_row[cbind(df$cell_id, df$year_id)],
                          cotton_column],
  net_filled = !is.na(filled_value) & filled_value > 0,
  lon = terra::xyFromCell(mask, f$unique_cells)[df$cell_id, 1],
  lat = terra::xyFromCell(mask, f$unique_cells)[df$cell_id, 2]) |>
  filter(class == "Pyrethroids", is.na(from) | year >= from) |>
  mutate(w = w / mean(w), country = factor(country),
         citation = droplevels(citation),
         source = relevel(factor(source), "MTM"),
         species = relevel(factor(species), "complex only"),
         protocol = relevel(factor(protocol), "WHO, version unknown"))

sink(file.path("outputs/hinge_screen",
               sprintf("%s_misfit_regression%s.txt", label,
                       if (is.na(from)) "" else paste0("_from", from))),
     split = TRUE)
cat(sprintf("%d pyrethroid bioassays, %d countries, %d citations\n\n",
            nrow(d), nlevels(d$country), nlevels(d$citation)))

# 1. by country
by_country <- d |>
  group_by(country) |>
  summarise(n = n(), mean_residual = weighted.mean(residual, w),
            .groups = "drop") |>
  arrange(mean_residual)
cat("weighted mean level residual by country (positive: the fit predicts",
    "too much resistance), countries with >= 50 bioassays:\n")
print(as.data.frame(filter(by_country, n >= 50) |>
                      mutate(mean_residual = round(mean_residual, 2))),
      row.names = FALSE)

# 2. deviance explained
candidates <- paste("source + species + protocol + log_n + n_100 +",
                    "s(year, k = 5) + nets + cum_nets + irs + pop + cotton +",
                    "net_filled + type")
models <- list(
  country = "residual ~ country",
  space = "residual ~ s(lon, lat, k = 80)",
  citation = "residual ~ s(citation, bs = 're')",
  candidates = paste("residual ~", candidates),
  "candidates + country" = paste("residual ~", candidates, "+ country"),
  "candidates + citation" = paste("residual ~", candidates,
                                  "+ s(citation, bs = 're')"),
  "country + citation" = "residual ~ country + s(citation, bs = 're')",
  "space + country" = "residual ~ country + s(lon, lat, k = 80)")
fits <- lapply(models, function(m) {
  bam(as.formula(m), data = d, weights = w, discrete = TRUE)
})
cat("\ndeviance explained (weighted):\n")
print(data.frame(model = names(fits),
                 dev_expl = round(vapply(fits, function(x) {
                   summary(x)$dev.expl
                 }, numeric(1)), 3),
                 aic = round(vapply(fits, AIC, numeric(1)))), row.names = FALSE)

# 3. the candidates' coefficients, without and with country effects
show <- function(fit, title) {
  s <- summary(fit)
  p <- s$p.table
  keep <- !grepl("^country|^\\(Intercept\\)|^type", rownames(p))
  cat(sprintf("\n%s\n", title))
  print(round(p[keep, c("Estimate", "Std. Error", "Pr(>|t|)"), drop = FALSE],
              3))
  if (nrow(s$s.table)) print(round(s$s.table, 3))
}
show(fits[["candidates + citation"]],
     "candidates, with a citation random effect")
show(fits[["candidates + country"]], "candidates, with country effects")

# the country effects of the most misfitted countries, without and with the
# candidates
country_effects <- function(fit) {
  p <- summary(fit)$p.table
  rows <- grep("^country", rownames(p))
  setNames(p[rows, "Estimate"], sub("^country", "", rownames(p)[rows]))
}
alone <- country_effects(fits$country)
adjusted <- country_effects(fits[["candidates + country"]])
baseline <- levels(d$country)[1]
shown <- by_country |> filter(n >= 50) |>
  filter(rank(mean_residual) <= 5 | rank(-mean_residual) <= 5)
cat(sprintf("\ncountry effects (relative to %s), alone and with the",
            baseline), "candidates:\n")
print(data.frame(country = shown$country, n = shown$n,
                 alone = round(alone[as.character(shown$country)], 2),
                 with_candidates = round(adjusted[as.character(shown$country)],
                                         2)),
      row.names = FALSE)
sink()
