# Pre-screen of the net-use transform (#37), step 3: the within-pixel pair
# analysis of #23 (R/selection_shape_analysis.R) redone on the trended design
# of the current fit, with saturating net-use transforms as variants.
#
#   Rscript R/net_screen_pairs.R
#
# Reads outputs/net_screen/extract.rds (R/net_screen_extract.R): the fit's
# bioassays and its design matrix, so the columns are exactly the model's
# (pyrethroid-weighted net use, IRS, population by the encounter transform
# times g_dom, the crops times g_ag), and its selection effects for the
# model-implied rates.
#
# As in #23: for consecutive records (pooled assays) at the same pixel and
# insecticide type, the per-year increase in logit resistance is the interval
# mean of log w_r, w_r = 1 + x_r' b, b >= 0, plus the reversion constant
# kappa <= 0 (#24; the model has it now). Fitted per class by weighted least
# squares with the variance model of #23, and compared by held-out log density
# over 1-degree blocks of pixels (10 folds, 3 repeats). The net-use column is
# replaced by each candidate transform; the others are unchanged:
#   linear          n, the model now
#   no_nets         no net-use column
#   linear_hinge_k  n and min(n, k), k at the quartiles of n over pair-years
#   hinge_k         min(n, k) alone
#   current_h       n / (n + h), saturating in the current year's use
#   cumulative_c    F(C_t) - F(C_{t-1}), F(C) = c (1 - exp(-C / c)), C the
#                   running sum of n from 1995: selection per unit of net use
#                   decays as exp(-C / c) with the cumulative use, and the
#                   cumulative selection saturates at c times the coefficient
#
# Writes outputs/net_screen/pairs_cv.csv, pairs_coefs.csv, pairs_net_slope.csv
# and figures/net_screen_pairs_cv.png.

suppressMessages({
  library(tidyverse)
})
output_dir <- "outputs/net_screen"
x <- readRDS(file.path(output_dir, "extract.rds"))
covariate_names <- colnames(x$x_cell_years)
stopifnot(covariate_names[1] == "nets")

# the replicate-assay rho per type (R/fig_illustrate_bioassay_variability.R),
# as in #23
rho_type <- read_csv("outputs/refit/bioassay_rho_hierarchical.csv",
                     show_col_types = FALSE) |>
  select(insecticide_type, rho)

covs <- bind_cols(x$cell_years_index, as_tibble(x$x_cell_years)) |>
  arrange(cell_id, year_id) |>
  group_by(cell_id) |>
  mutate(cum_nets = cumsum(nets)) |>
  ungroup()


# data: records and consecutive pairs, as #23 ---------------------------------

records <- x$assay |>
  left_join(rho_type, by = "insecticide_type") |>
  group_by(cell_id) |>
  mutate(latitude = first(latitude), longitude = first(longitude),
         country_name = first(country_name)) |>
  group_by(cell_id, latitude, longitude, country_name, insecticide_class,
           insecticide_type, type_id, year_id) |>
  summarise(n_eff_factor = sum(n * (1 + (n - 1) * rho)),
            died = sum(died), n = sum(n), .groups = "drop") |>
  mutate(p_hat = (died + 0.5) / (n + 1),
         logit_mort = qlogis(p_hat),
         var_logit = n_eff_factor / (n ^ 2 * p_hat * (1 - p_hat)))

pairs <- records |>
  arrange(cell_id, type_id, year_id) |>
  group_by(cell_id, type_id) |>
  mutate(year_1 = lag(year_id), logit_mort_1 = lag(logit_mort),
         var_logit_1 = lag(var_logit)) |>
  ungroup() |>
  filter(!is.na(year_1)) |>
  transmute(cell_id, latitude, longitude, country_name, insecticide_class,
            insecticide_type, type_id, year_1, year_2 = year_id,
            dt = year_2 - year_1,
            y = (logit_mort_1 - logit_mort) / dt,
            var_y = (var_logit_1 + var_logit) / dt ^ 2) |>
  mutate(pair_id = row_number(),
         block = paste(floor(latitude), floor(longitude)))

# the fitness of years year_1 + 1, ..., year_2 acts between the two records
pair_years <- pairs |>
  select(pair_id, cell_id, year_1, year_2) |>
  mutate(year_id = map2(year_1 + 1, year_2, seq)) |>
  unnest(year_id) |>
  left_join(covs, by = c("cell_id", "year_id")) |>
  select(-year_1, -year_2)
stopifnot(!anyNA(pair_years))
pairs <- pairs |>
  left_join(pair_years |> group_by(pair_id) |>
              summarise(nets = mean(nets), cum_nets_2 = last(cum_nets)),
            by = "pair_id")
pairs |> count(insecticide_class) |> print()


# the net-use transforms -------------------------------------------------------

# F(C) = c (1 - exp(-C / c)), the saturating cumulative net use; its
# increments are the yearly column. Written with expm1 so that it is accurate
# when C << c, where F(C) ~ C
saturating_cumulative <- function(cum, c) c * -expm1(-cum / c)
cumulative_column <- function(c) {
  cum <- pair_years$cum_nets
  saturating_cumulative(cum, c) - saturating_cumulative(cum - pair_years$nets,
                                                         c)
}

knots <- unname(quantile(pair_years$nets, c(0.25, 0.5, 0.75)))
other_columns <- as.matrix(pair_years[, covariate_names[-1]])
net_columns <- c(
  list(linear = cbind(nets = pair_years$nets),
       no_nets = NULL),
  setNames(lapply(seq_along(knots), function(i) {
    cbind(nets = pair_years$nets, hinge = pmin(pair_years$nets, knots[i]))
  }), sprintf("linear_hinge_%.2f", knots)),
  setNames(lapply(knots, function(k) cbind(nets = pmin(pair_years$nets, k))),
           sprintf("hinge_%.2f", knots)),
  setNames(lapply(c(0.05, 0.1, 0.2, 0.4), function(h) {
    cbind(nets = pair_years$nets / (pair_years$nets + h))
  }), sprintf("current_%g", c(0.05, 0.1, 0.2, 0.4))),
  setNames(lapply(c(1, 2, 3, 4, 6, 8), function(c) {
    cbind(nets = cumulative_column(c))
  }), sprintf("cumulative_%g", c(1, 2, 3, 4, 6, 8)))
)
bases <- lapply(net_columns, function(n) cbind(n, other_columns))


# fitting, as #23 ----------------------------------------------------------------

# weighted least squares for the per-year rate, the interval mean of
# log(1 + x' b) plus kappa, with b >= 0 and kappa <= 0, by L-BFGS-B with an
# analytic gradient. `rows` indexes pairs.
fit_pairs <- function(x, rows, weights) {
  keep <- pair_years$pair_id %in% rows
  x <- x[keep, , drop = FALSE]
  group <- match(pair_years$pair_id[keep], rows)
  dt <- pairs$dt[rows]
  y <- pairs$y[rows]
  w <- weights[rows]
  n_b <- ncol(x)
  rate <- function(par) {
    rowsum(log1p(x %*% par[seq_len(n_b)]), group)[, 1] / dt + par[n_b + 1]
  }
  objective <- function(par) sum(w * (y - rate(par)) ^ 2)
  gradient <- function(par) {
    resid <- y - rate(par)
    d_rate <- cbind(rowsum(x / c(1 + x %*% par[seq_len(n_b)]), group) / dt,
                    1)
    -2 * colSums(w * resid * d_rate)
  }
  fits <- lapply(c(0.05, 0.5), function(start) {
    optim(c(rep(start, n_b), 0), objective, gradient, method = "L-BFGS-B",
          lower = c(rep(0, n_b), -Inf), upper = c(rep(Inf, n_b), 0),
          control = list(maxit = 5000, factr = 1e5))
  })
  best <- fits[[which.min(vapply(fits, `[[`, 0, "value"))]]
  setNames(best$par, c(colnames(x), "kappa"))
}

predict_pairs <- function(x, par, rows) {
  keep <- pair_years$pair_id %in% rows
  group <- match(pair_years$pair_id[keep], rows)
  n_b <- ncol(x)
  rowsum(log1p(x[keep, , drop = FALSE] %*% par[seq_len(n_b)]),
         group)[, 1] / pairs$dt[rows] + par[n_b + 1]
}

class_rows <- split(pairs$pair_id, pairs$insecticide_class)

# the variance of each pair: phi times the approximate sampling variance plus
# 2 sigma2 / dt^2 for between-sample variation that rho does not cover,
# estimated per class under the linear variant and then held fixed (#23)
pair_var <- function(phi, sigma2) phi * pairs$var_y + 2 * sigma2 / pairs$dt ^ 2
dispersion <- map_dfr(class_rows, function(rows) {
  par <- c(1, 0.1)
  for (iter in 1:5) {
    fit <- fit_pairs(bases$linear, rows, 1 / pair_var(par[1], par[2]))
    resid <- pairs$y[rows] - predict_pairs(bases$linear, fit, rows)
    par <- optim(par, function(p) {
      -sum(dnorm(resid, 0, sqrt(pair_var(p[1], p[2])[rows]), log = TRUE))
    }, method = "L-BFGS-B", lower = c(1e-3, 0), upper = c(100, 20))$par
  }
  tibble(phi = par[1], sigma2 = par[2])
}, .id = "insecticide_class")
print(dispersion)
d <- dispersion[match(pairs$insecticide_class, dispersion$insecticide_class), ]
pairs$var_total <- pair_var(d$phi, d$sigma2)
weights <- 1 / pairs$var_total

coefs <- imap_dfr(bases, function(x, variant) {
  imap_dfr(class_rows, function(rows, class) {
    par <- fit_pairs(x, rows, weights)
    tibble(variant = variant, insecticide_class = class,
           term = names(par), estimate = unname(par))
  })
})
write_csv(coefs, file.path(output_dir, "pairs_coefs.csv"))
coefs |> filter(term %in% c("nets", "hinge", "kappa")) |>
  pivot_wider(names_from = term, values_from = estimate) |>
  print(n = Inf)


# cross-validation over 1-degree blocks -----------------------------------------

log_density <- function(rows, pred) {
  dnorm(pairs$y[rows], pred, sqrt(pairs$var_total[rows]), log = TRUE)
}
n_folds <- 10
n_repeats <- 3
blocks <- unique(pairs$block)
set.seed(2026)
pair_folds <- map(seq_len(n_repeats), function(r) {
  sample(rep_len(seq_len(n_folds), length(blocks)))[match(pairs$block,
                                                          blocks)]
})
cv_scores <- map_dfr(seq_len(n_repeats), function(r) {
  imap_dfr(bases, function(x, variant) {
    map_dfr(seq_len(n_folds), function(fold) {
      imap_dfr(class_rows, function(rows, class) {
        train <- rows[pair_folds[[r]][rows] != fold]
        test <- rows[pair_folds[[r]][rows] == fold]
        par <- fit_pairs(x, train, weights)
        tibble(repeat_id = r, variant = variant, insecticide_class = class,
               pair_id = test,
               log_density = log_density(test, predict_pairs(x, par, test)))
      })
    })
  })
})

# held-out log density summed over pairs (mean over repeats), its difference
# from the linear variant, and the standard error of that from the spread of
# per-block differences
cv_pair <- cv_scores |>
  group_by(variant, insecticide_class, pair_id) |>
  summarise(log_density = mean(log_density), .groups = "drop") |>
  left_join(pairs |> select(pair_id, block), by = "pair_id")
summarise_cv <- function(data) {
  base <- data |> filter(variant == "linear") |>
    select(pair_id, base = log_density)
  data |>
    left_join(base, by = "pair_id") |>
    group_by(variant, block) |>
    summarise(diff = sum(log_density - base),
              log_density = sum(log_density), .groups = "drop_last") |>
    summarise(log_density = sum(log_density), diff_vs_linear = sum(diff),
              se = sqrt(n()) * sd(diff), .groups = "drop")
}
cv_table <- bind_rows(
  cv_pair |> summarise_cv() |> mutate(insecticide_class = "all"),
  cv_pair |> split(~insecticide_class) |>
    imap_dfr(~ summarise_cv(.x) |> mutate(insecticide_class = .y))
) |>
  mutate(variant = factor(variant, names(bases))) |>
  arrange(insecticide_class, variant) |>
  left_join(pairs |> count(insecticide_class, name = "n_pairs") |>
              bind_rows(tibble(insecticide_class = "all",
                               n_pairs = nrow(pairs))),
            by = "insecticide_class")
write_csv(cv_table, file.path(output_dir, "pairs_cv.csv"))
print(cv_table, n = Inf)


# net-use slope: data against the fit ---------------------------------------------

# the per-year change against interval mean net use, for pyrethroid pairs, in
# the data and in the fit's implied rate (posterior mean over its draws of the
# interval mean of log w plus kappa), with a bootstrap over pixels
effect <- x$effects$effect_type
kappa <- x$effects$kappa_type
x_pairs <- as.matrix(pair_years[, covariate_names])
pair_type <- pairs$type_id[pair_years$pair_id]
model_rate <- matrix(NA_real_, dim(effect)[1], nrow(pairs))
for (s in seq_len(dim(effect)[1])) {
  log_w <- log1p(x_pairs %*% effect[s, , ])
  log_w <- log_w[cbind(seq_along(pair_type), pair_type)]
  model_rate[s, ] <- rowsum(log_w, pair_years$pair_id)[, 1] / pairs$dt +
    kappa[s, pairs$type_id]
}
pairs$model_rate <- colMeans(model_rate)
rm(model_rate)

cluster_boot_slope <- function(data, response, term, n_boot = 1000,
                               seed = 11) {
  fit_slope <- function(d) {
    unname(coef(lm(d[[response]] ~ d[[term]], weights = 1 / d$var_total))[2])
  }
  set.seed(seed)
  rows_by_cell <- split(seq_len(nrow(data)), data$cell_id)
  boots <- replicate(n_boot, {
    fit_slope(data[unlist(rows_by_cell[sample(names(rows_by_cell),
                                              replace = TRUE)]), ])
  })
  tibble(response = response, term = term, slope = fit_slope(data),
         lower = quantile(boots, 0.025), upper = quantile(boots, 0.975),
         n_pairs = nrow(data), n_cells = n_distinct(data$cell_id))
}
pyrethroid_pairs <- pairs |> filter(insecticide_class == "Pyrethroids")
net_slope <- bind_rows(
  cluster_boot_slope(pyrethroid_pairs, "y", "nets"),
  cluster_boot_slope(pyrethroid_pairs, "model_rate", "nets"),
  cluster_boot_slope(pyrethroid_pairs, "y", "cum_nets_2"),
  cluster_boot_slope(pyrethroid_pairs, "model_rate", "cum_nets_2"))
print(net_slope)
write_csv(net_slope, file.path(output_dir, "pairs_net_slope.csv"))

# binned observed and implied rates against cumulative net use, for the figure
binned <- pairs |>
  mutate(cum_bin = cut(cum_nets_2, c(-0.01, 1, 2, 3, 4, 6, 11))) |>
  group_by(insecticide_class, cum_bin) |>
  summarise(observed = weighted.mean(y, 1 / var_total),
            observed_se = sqrt(1 / sum(1 / var_total)),
            implied = weighted.mean(model_rate, 1 / var_total),
            n = n(), .groups = "drop")
write_csv(binned, file.path(output_dir, "pairs_binned.csv"))
print(binned, n = Inf)


# figure ----------------------------------------------------------------------

family_of <- function(variant) {
  case_when(variant %in% c("linear", "no_nets") ~ "reference",
            str_starts(variant, "linear_hinge") ~ "linear + hinge",
            str_starts(variant, "hinge") ~ "hinge",
            str_starts(variant, "current") ~ "current use n / (n + h)",
            str_starts(variant, "cumulative") ~ "cumulative c (1 - exp(-C/c))")
}
p <- cv_table |>
  filter(insecticide_class %in% c("all", "Pyrethroids", "Organochlorines")) |>
  mutate(family = family_of(as.character(variant)),
         variant = fct_rev(variant)) |>
  ggplot(aes(diff_vs_linear, variant, colour = family)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_pointrange(aes(xmin = diff_vs_linear - se,
                      xmax = diff_vs_linear + se), size = 0.25) +
  facet_wrap(~insecticide_class, scales = "free_x") +
  scale_colour_manual(values = c("#0072B2", "#D55E00", "#009E73",
                                 "#CC79A7", "grey40"), name = NULL) +
  labs(x = "held-out log density vs linear net use (± 1 se, block CV)",
       y = NULL,
       title = "Within-pixel pairs: net-use transforms, trended design") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave("figures/net_screen_pairs_cv.png", p, width = 9, height = 5.5,
       dpi = 150, bg = "white")
cat("done\n")
