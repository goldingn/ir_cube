# Paired-difference analysis of the shape of the covariate effects on selection
# (#23).
#
# The model's selection recursion is additive on the logit scale:
#   logit q_t2 - logit q_t1 = -sum_{r = t1 + 1}^{t2} log w_r,
#   w_r = 1 + x_r' b,  b >= 0,
# where q is the fraction susceptible (= expected bioassay mortality). So for
# repeat bioassays at the same pixel and insecticide, the per-year change in
# logit resistance is the mean of log w_r over the interval between them, with
# no dependence on the initial state. This fits that relationship directly to
# empirical logits, with the model's covariates, and asks whether adding
# saturating hinge bases min(x, k) (non-negative coefficients, so the effect
# stays non-decreasing and concave) improves cross-validated fit over the
# linear terms the model uses now.
#
# Outputs:
#   outputs/selection_shape_cv.csv      cross-validated comparison of variants
#   outputs/selection_shape_coefs.csv   coefficients of the variants, all data
#   outputs/selection_shape_net_slope.csv  univariate net-use slope check

library(tidyverse)
library(terra)
source("R/functions.R")

# build df and the cell-year covariates exactly as fit_model.R does, by
# evaluating its code up to the design matrix (skipping the package load, which
# pulls in greta)
fit_model_exprs <- parse("R/fit_model.R")
fit_model_text <- vapply(fit_model_exprs,
                         function(e) paste(deparse(e), collapse = " "),
                         "")
last_expr <- which(startsWith(fit_model_text, "x_cell_years <-"))
stopifnot(length(last_expr) == 1)
for (i in seq_len(last_expr)) {
  if (grepl("^source\\(\"R/(packages|functions).R\"\\)", fit_model_text[i])) {
    next
  }
  eval(fit_model_exprs[[i]], envir = globalenv())
}

# covariates for all cell-years with data, in the model's column order
covs <- bind_cols(cell_years_index, as_tibble(x_cell_years))
crop_names <- setdiff(colnames(x_cell_years), c("nets", "irs", "pop"))

# log population density, min-max scaled over the whole cube (2000-2030) as
# prep_rasters.R would with its log transform switched on. This is the column
# proposed to replace pop. The 2031-2049 projections that prep_rasters.R also
# scales over are not in data/clean; they move the scaling slightly, not the
# ranking.
pop_raw <- c(rast("data/clean/pop_cube.tif"),
             rast("data/clean/pop_cube_future.tif"))
log_pop_range <- range(log(global(pop_raw, "range", na.rm = TRUE)))
pop_raw <- pre_pad_cube(pop_raw[[paste0("pop_", 2000:final_data_year)]],
                        baseline_year)
pop_raw_extract <- terra::extract(pop_raw, unique_cells) %>%
  mutate(cell_id = seq_along(unique_cells)) %>%
  pivot_longer(-cell_id,
               names_prefix = "pop_",
               names_to = "year",
               values_to = "pop_raw") %>%
  mutate(year_id = as.numeric(year) - baseline_year + 1) %>%
  select(cell_id, year_id, pop_raw)
covs <- covs %>%
  left_join(pop_raw_extract, by = c("cell_id", "year_id")) %>%
  mutate(
    log_pop = (log(pop_raw) - log_pop_range[1]) / diff(log_pop_range),
    # a constant, for selection not tied to any covariate
    one = 1
  )
stopifnot(!anyNA(covs$log_pop))

# check the raw population reproduces the model's scaled column (a linear map,
# up to the per-year relevelling in prep_rasters.R)
stopifnot(cor(covs$pop_raw, covs$pop) > 0.999)


# data: one record per pixel, insecticide type and year

rho_type <- read_csv("outputs/bioassay_rho_hierarchical.csv",
                     show_col_types = FALSE) %>%
  select(insecticide_type, rho)

# empirical logit of mortality with a continuity correction, and its
# approximate variance: binomial variance on the logit scale, inflated for
# beta-binomial overdispersion within each assay, and combining the assays
# pooled into the record
records <- df %>%
  left_join(rho_type, by = "insecticide_type") %>%
  group_by(cell_id, cell, latitude, longitude,
           insecticide_class, insecticide_type, type_id, year_id) %>%
  summarise(
    died = sum(died),
    n = sum(mosquito_number),
    # sum of n_i (1 + (n_i - 1) rho) over pooled assays
    n_eff_factor = sum(mosquito_number * (1 + (mosquito_number - 1) * rho)),
    .groups = "drop"
  ) %>%
  # one location per pixel for the grouping below
  group_by(cell_id) %>%
  mutate(latitude = first(latitude),
         longitude = first(longitude)) %>%
  ungroup() %>%
  mutate(
    p_hat = (died + 0.5) / (n + 1),
    logit_mort = qlogis(p_hat),
    # var(logit p_hat) ~= var(p_hat) / (p (1 - p))^2, where
    # var(p_hat) = p (1 - p) sum n_i (1 + (n_i - 1) rho) / n^2
    var_logit = n_eff_factor / (n ^ 2 * p_hat * (1 - p_hat))
  ) %>%
  distinct(cell_id, type_id, year_id, .keep_all = TRUE)

# consecutive pairs of records at the same pixel and insecticide type. Using
# only consecutive pairs keeps each record in at most two differences.
pairs <- records %>%
  arrange(cell_id, type_id, year_id) %>%
  group_by(cell_id, type_id) %>%
  mutate(
    year_1 = lag(year_id),
    logit_mort_1 = lag(logit_mort),
    var_logit_1 = lag(var_logit),
    p_hat_1 = lag(p_hat)
  ) %>%
  ungroup() %>%
  filter(!is.na(year_1)) %>%
  transmute(
    cell_id, latitude, longitude,
    insecticide_class, insecticide_type, type_id,
    year_1,
    year_2 = year_id,
    dt = year_2 - year_1,
    p_hat_1,
    p_hat_2 = p_hat,
    # per-year increase in logit resistance = decrease in logit mortality
    y = (logit_mort_1 - logit_mort) / dt,
    var_y = (var_logit_1 + var_logit) / dt ^ 2
  ) %>%
  mutate(pair_id = row_number())

# cell-years in each interval: the fitness applied between the two records is
# that of years year_1 + 1, ..., year_2
pair_years <- pairs %>%
  select(pair_id, cell_id, year_1, year_2) %>%
  mutate(year_id = map2(year_1 + 1, year_2, seq)) %>%
  unnest(year_id) %>%
  left_join(covs, by = c("cell_id", "year_id")) %>%
  select(-year_1, -year_2)
stopifnot(!anyNA(pair_years))

# interval means of the covariates, for the univariate checks
pair_means <- pair_years %>%
  group_by(pair_id) %>%
  summarise(across(c(nets, irs, pop, log_pop, all_of(crop_names)), mean))
pairs <- pairs %>%
  left_join(pair_means, by = "pair_id")

pairs %>%
  count(insecticide_class, insecticide_type) %>%
  print()


# net-use slope checks. The figure quoted in #23 (-0.01 [-0.37, +0.34], model
# +0.61) used the 2014 forecasting fold's held-out pyrethroid assays paired with
# training assays at the same pixel, and that fold's posterior; it needs the
# fold's draws, so is not recomputed here. These are the same regression on all
# the data: first pyrethroid pairs 3-8 years apart (all pairs, not only
# consecutive ones), unweighted regression of the per-year change in logit
# resistance on net use averaged over the years y1, ..., y2, with a bootstrap
# over pixels

net_slope_pairs <- records %>%
  filter(insecticide_class == "Pyrethroids") %>%
  select(cell_id, type_id, year_id, logit_mort, var_logit) %>%
  inner_join(., ., by = c("cell_id", "type_id"),
             suffix = c("_1", "_2"),
             relationship = "many-to-many") %>%
  mutate(dt = year_id_2 - year_id_1) %>%
  filter(dt >= 3, dt <= 8) %>%
  mutate(
    y = (logit_mort_1 - logit_mort_2) / dt,
    var_y = (var_logit_1 + var_logit_2) / dt ^ 2,
    pair_id = row_number()
  )
net_slope_pairs <- net_slope_pairs %>%
  select(pair_id, cell_id, year_id_1, year_id_2) %>%
  mutate(year_id = map2(year_id_1, year_id_2, seq)) %>%
  unnest(year_id) %>%
  left_join(covs, by = c("cell_id", "year_id")) %>%
  group_by(pair_id) %>%
  summarise(nets = mean(nets)) %>%
  right_join(net_slope_pairs, by = "pair_id")

# the same regression on the consecutive pairs used below, weighted by the
# inverse approximate variance, with net use averaged over the interval the
# model applies (years y1 + 1, ..., y2)
net_slope_consecutive <- pairs %>%
  filter(insecticide_class == "Pyrethroids")

cluster_boot_slope <- function(data, weighted, n_boot = 1000, seed = 11) {
  fit_slope <- function(d) {
    w <- if (weighted) 1 / d$var_y else NULL
    unname(coef(lm(y ~ nets, data = d, weights = w))[2])
  }
  set.seed(seed)
  rows_by_cell <- split(seq_len(nrow(data)), data$cell_id)
  boots <- replicate(n_boot, {
    cells_boot <- sample(names(rows_by_cell), replace = TRUE)
    fit_slope(data[unlist(rows_by_cell[cells_boot]), ])
  })
  tibble(
    slope = fit_slope(data),
    lower = quantile(boots, 0.025),
    upper = quantile(boots, 0.975),
    n_pairs = nrow(data),
    n_cells = n_distinct(data$cell_id)
  )
}

net_slope <- bind_rows(
  cluster_boot_slope(net_slope_pairs, weighted = FALSE) %>%
    mutate(pairs = "all pairs 3-8 years apart, unweighted"),
  cluster_boot_slope(net_slope_pairs, weighted = TRUE) %>%
    mutate(pairs = "all pairs 3-8 years apart, inverse-variance weighted"),
  cluster_boot_slope(net_slope_consecutive, weighted = FALSE) %>%
    mutate(pairs = "consecutive pairs, unweighted"),
  cluster_boot_slope(net_slope_consecutive, weighted = TRUE) %>%
    mutate(pairs = "consecutive pairs, inverse-variance weighted")
) %>%
  relocate(pairs)
print(net_slope)


# model-implied per-year change for the same consecutive pairs, under the
# saved full-data fit (posterior mean over draws of the implied rate)

fitted <- new.env()
load("temporary/fitted_model.RData", envir = fitted)
stopifnot(identical(fitted$types, types),
          identical(fitted$classes_index, classes_index),
          identical(colnames(fitted$x_cell_years), colnames(x_cell_years)))
draws_mat <- as.matrix(fitted$draws)
draws_mat <- draws_mat[round(seq(1, nrow(draws_mat), length.out = 500)), ]
draw_array <- function(name, dim) {
  cols <- grep(paste0("^", name, "\\["), colnames(draws_mat))
  array(draws_mat[, cols], c(nrow(draws_mat), dim))
}
n_covs <- ncol(x_cell_years)
n_classes <- length(classes)
n_types <- length(types)
beta_overall <- draw_array("beta_overall", n_covs)
sigma_overall <- draw_array("sigma_overall", n_covs)
sigma_class <- draw_array("sigma_class", n_covs)
beta_class_raw <- draw_array("beta_class_raw", c(n_covs, n_classes))
beta_type_raw <- draw_array("beta_type_raw", c(n_covs, n_types))

x_pair_years <- as.matrix(pair_years[, colnames(x_cell_years)])
pair_type <- pairs$type_id[pair_years$pair_id]
model_rate <- matrix(NA, nrow(draws_mat), nrow(pairs))
for (s in seq_len(nrow(draws_mat))) {
  beta_class <- beta_overall[s, ] + sigma_overall[s, ] * beta_class_raw[s, , ]
  beta_type <- beta_class[, classes_index] +
    sigma_class[s, ] * beta_type_raw[s, , ]
  log_w <- log1p(x_pair_years %*% exp(beta_type))
  log_w <- log_w[cbind(seq_along(pair_type), pair_type)]
  model_rate[s, ] <- rowsum(log_w, pair_years$pair_id)[, 1] / pairs$dt
}
pairs$model_rate <- colMeans(model_rate)
rm(fitted, draws_mat, model_rate)

net_slope <- bind_rows(
  net_slope,
  cluster_boot_slope(pairs %>%
                       filter(insecticide_class == "Pyrethroids") %>%
                       mutate(y = model_rate),
                     weighted = TRUE) %>%
    mutate(pairs = "consecutive pairs, weighted, full-data model's implied rate")
)
print(net_slope)
write_csv(net_slope, "outputs/selection_shape_net_slope.csv")


# fit w = 1 + f(x)' b per insecticide class to the paired differences

# knots for the hinge bases min(x, k): quartiles of the covariate over the
# pair-years (for IRS, of its non-zero values; two thirds of pair-years have
# none)
knot_probs <- c(0.25, 0.5, 0.75)
knots <- list(
  nets = unname(quantile(pair_years$nets, knot_probs)),
  irs = unname(quantile(pair_years$irs[pair_years$irs > 0], knot_probs))
)
print(knots)

# variants: linear terms, hinge covariates and knot indices, and whether a
# constant (non-positive) change in logit resistance per year is allowed, as
# the reversion in #24 would add. log_pop spans 0.48-0.97 at the data, so it
# acts mostly as a constant; the _const variants add an explicit constant
# column to separate the two.
base_pop <- c("nets", "irs", "pop", crop_names)
base_log_pop <- c("nets", "irs", "log_pop", crop_names)
base_const <- c("one", "nets", "irs", "log_pop", crop_names)
variants <- list(
  linear_pop = list(linear = base_pop),
  linear_pop_const = list(linear = c("one", base_pop)),
  linear = list(linear = base_log_pop),
  linear_const = list(linear = base_const),
  nets_k1 = list(hinge = list(nets = 1)),
  nets_k2 = list(hinge = list(nets = 2)),
  nets_k3 = list(hinge = list(nets = 3)),
  nets_k123 = list(hinge = list(nets = 1:3)),
  irs_k1 = list(hinge = list(irs = 1)),
  irs_k2 = list(hinge = list(irs = 2)),
  irs_k3 = list(hinge = list(irs = 3)),
  irs_k123 = list(hinge = list(irs = 1:3)),
  both_k123 = list(hinge = list(nets = 1:3, irs = 1:3)),
  const_nets_k1 = list(linear = base_const, hinge = list(nets = 1)),
  const_nets_k123 = list(linear = base_const, hinge = list(nets = 1:3)),
  const_irs_k123 = list(linear = base_const, hinge = list(irs = 1:3)),
  linear_reversion = list(reversion = TRUE),
  const_reversion = list(linear = base_const, reversion = TRUE),
  nets_k123_reversion = list(hinge = list(nets = 1:3), reversion = TRUE)
)
variants <- map(variants, function(v) {
  list(linear = v$linear %||% base_log_pop,
       hinge = v$hinge %||% list(),
       reversion = v$reversion %||% FALSE)
})

# basis matrix over the pair-year rows
make_basis <- function(variant) {
  x <- as.matrix(pair_years[, variant$linear])
  for (cov in names(variant$hinge)) {
    for (k in knots[[cov]][variant$hinge[[cov]]]) {
      x <- cbind(x, pmin(pair_years[[cov]], k))
      colnames(x)[ncol(x)] <- sprintf("min(%s, %.3f)", cov, k)
    }
  }
  x
}

# weighted least squares for the interval mean of log w, with b >= 0 and the
# reversion constant <= 0, by L-BFGS-B with an analytic gradient. `rows`
# indexes pairs.
fit_pairs <- function(x, rows, weights, reversion) {
  keep <- pair_years$pair_id %in% rows
  x <- x[keep, , drop = FALSE]
  group <- match(pair_years$pair_id[keep], rows)
  dt <- pairs$dt[rows]
  y <- pairs$y[rows]
  w <- weights[rows]
  n_b <- ncol(x)
  predict_rate <- function(par) {
    rate <- rowsum(log1p(x %*% par[seq_len(n_b)]), group)[, 1] / dt
    if (reversion) rate <- rate + par[n_b + 1]
    rate
  }
  objective <- function(par) sum(w * (y - predict_rate(par)) ^ 2)
  gradient <- function(par) {
    resid <- y - predict_rate(par)
    d_rate <- rowsum(x / c(1 + x %*% par[seq_len(n_b)]), group) / dt
    if (reversion) d_rate <- cbind(d_rate, 1)
    -2 * colSums(w * resid * d_rate)
  }
  n_par <- n_b + reversion
  starts <- list(rep(0.05, n_par), rep(0.5, n_par))
  fits <- map(starts, function(start) {
    if (reversion) start[n_par] <- 0
    optim(start, objective, gradient,
          method = "L-BFGS-B",
          lower = c(rep(0, n_b), if (reversion) -Inf),
          upper = c(rep(Inf, n_b), if (reversion) 0),
          control = list(maxit = 5000, factr = 1e5))
  })
  best <- fits[[which.min(map_dbl(fits, "value"))]]
  par <- best$par
  names(par) <- c(colnames(x), if (reversion) "reversion")
  list(par = par, convergence = best$convergence)
}

# predict for any pairs from fitted parameters
predict_pairs <- function(x, par, rows, reversion) {
  keep <- pair_years$pair_id %in% rows
  group <- match(pair_years$pair_id[keep], rows)
  n_b <- ncol(x)
  rate <- rowsum(log1p(x[keep, , drop = FALSE] %*% par[seq_len(n_b)]),
                 group)[, 1] / pairs$dt[rows]
  if (reversion) rate <- rate + par[n_b + 1]
  rate
}

class_rows <- split(pairs$pair_id, pairs$insecticide_class)

# variance of each pair: the approximate sampling variance scaled by phi, plus
# an extra variance sigma2 per record for between-sample variation that the
# beta-binomial rho does not cover. phi < 1 is expected: the logit variance
# approximation overstates the scatter of records at or near 100% mortality.
# Both are estimated per class by maximum likelihood under the linear variant,
# alternating with the coefficients, then held fixed for every variant.
pair_var <- function(phi, sigma2) {
  phi * pairs$var_y + 2 * sigma2 / pairs$dt ^ 2
}
x_linear <- make_basis(variants$linear)
dispersion_class <- map_dfr(class_rows, function(rows) {
  par <- c(0, 0.1)
  for (iter in 1:5) {
    v <- pair_var(par[1], par[2])
    fit <- fit_pairs(x_linear, rows, 1 / v, reversion = FALSE)
    resid <- pairs$y[rows] - predict_pairs(x_linear, fit$par, rows, FALSE)
    par <- optim(par, function(p) {
      v <- pair_var(p[1], p[2])[rows]
      -sum(dnorm(resid, 0, sqrt(v), log = TRUE))
    }, method = "L-BFGS-B", lower = c(1e-3, 0), upper = c(100, 20))$par
  }
  tibble(phi = par[1], sigma2 = par[2])
}, .id = "insecticide_class")
print(dispersion_class)
pair_dispersion <- dispersion_class[match(pairs$insecticide_class,
                                          dispersion_class$insecticide_class), ]
pairs$var_total <- pair_var(pair_dispersion$phi, pair_dispersion$sigma2)
pair_weights <- 1 / pairs$var_total

# full-data fits of every variant, per class
bases <- map(variants, make_basis)
full_fits <- imap(variants, function(variant, name) {
  map(class_rows, function(rows) {
    fit_pairs(bases[[name]], rows, pair_weights, variant$reversion)
  })
})
stopifnot(all(unlist(map_depth(full_fits, 2, "convergence")) == 0))

# mean squared standardised residual under the linear variant (close to 1
# if the variance model is calibrated)
imap_dbl(class_rows, function(rows, class) {
  pred <- predict_pairs(bases$linear, full_fits$linear[[class]]$par, rows,
                        FALSE)
  mean((pairs$y[rows] - pred) ^ 2 / pairs$var_total[rows])
}) %>%
  print()

coefs <- imap_dfr(full_fits, function(fits, name) {
  imap_dfr(fits, function(fit, class) {
    tibble(variant = name, insecticide_class = class,
           term = names(fit$par), estimate = unname(fit$par))
  })
})
write_csv(coefs, "outputs/selection_shape_coefs.csv")


# cross-validation by blocks of pixels: 1-degree squares, all pairs in a
# square (every insecticide) held out together, 10 folds, 5 repeats

pairs$block <- paste(floor(pairs$latitude), floor(pairs$longitude))
blocks <- unique(pairs$block)
n_folds <- 10
n_repeats <- 5

log_density <- function(rows, pred) {
  dnorm(pairs$y[rows], pred, sqrt(pairs$var_total[rows]), log = TRUE)
}

set.seed(2026)
cv_scores <- map_dfr(seq_len(n_repeats), function(rep) {
  block_fold <- sample(rep_len(seq_len(n_folds), length(blocks)))
  pair_fold <- block_fold[match(pairs$block, blocks)]
  imap_dfr(variants, function(variant, name) {
    map_dfr(seq_len(n_folds), function(fold) {
      imap_dfr(class_rows, function(rows, class) {
        train <- rows[pair_fold[rows] != fold]
        test <- rows[pair_fold[rows] == fold]
        fit <- fit_pairs(bases[[name]], train, pair_weights, variant$reversion)
        pred <- predict_pairs(bases[[name]], fit$par, test, variant$reversion)
        tibble(repeat_id = rep, variant = name, insecticide_class = class,
               pair_id = test, pred = pred,
               log_density = log_density(test, pred))
      })
    })
  })
})

# in-sample log density of the saved full-data model's implied rates, as a
# reference (not cross-validated, so favoured)
model_reference <- pairs %>%
  mutate(log_density = log_density(pair_id, model_rate)) %>%
  group_by(insecticide_class) %>%
  summarise(log_density = sum(log_density))
print(model_reference)

# summarise: held-out log density summed over pairs (mean over repeats),
# difference from the linear variant, and its standard error from the spread
# of per-block differences
cv_pair <- cv_scores %>%
  group_by(variant, insecticide_class, pair_id) %>%
  summarise(log_density = mean(log_density), .groups = "drop") %>%
  left_join(pairs %>% select(pair_id, block), by = "pair_id")

summarise_cv <- function(data) {
  base <- data %>%
    filter(variant == "linear") %>%
    select(pair_id, base = log_density)
  data %>%
    left_join(base, by = "pair_id") %>%
    group_by(variant, block) %>%
    summarise(diff = sum(log_density - base),
              log_density = sum(log_density),
              .groups = "drop_last") %>%
    summarise(log_density = sum(log_density),
              diff_vs_linear = sum(diff),
              se = sqrt(n()) * sd(diff),
              .groups = "drop")
}
cv_table <- bind_rows(
  cv_pair %>%
    summarise_cv() %>%
    mutate(insecticide_class = "all"),
  cv_pair %>%
    split(.$insecticide_class) %>%
    imap_dfr(~ summarise_cv(.x) %>% mutate(insecticide_class = .y))
) %>%
  mutate(variant = factor(variant, names(variants))) %>%
  arrange(insecticide_class, variant) %>%
  left_join(
    pairs %>%
      count(insecticide_class, name = "n_pairs") %>%
      bind_rows(tibble(insecticide_class = "all", n_pairs = nrow(pairs))),
    by = "insecticide_class"
  )
print(cv_table, n = Inf)
write_csv(cv_table, "outputs/selection_shape_cv.csv")
