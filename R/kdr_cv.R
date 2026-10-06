# choose the smoothness of the kdr allele frequency surfaces by
# cross-validation (#46)

# load packages and functions
source("R/packages.R")
source("R/functions.R")
source("R/kdr_functions.R")
suppressMessages(library(mgcv))
library(parallel)

# The kdr surfaces (R/kdr_surface_check.R) are binomial GAMs of the 995F and
# 995S allele frequencies, with a spatial smooth and a linear year trend.
# Records cluster by study, and REML treats them as independent, so it fits
# differences between studies as spatial signal. Here we score smoothness by
# held-out binomial deviance, with two kinds of folds: studies held out
# together, and 5 x 5 degree blocks held out together. The candidates are a
# thin-plate spline with its smoothing parameter fixed on a grid, the same
# spline with REML smoothness selection, and Matern Gaussian process smooths
# (REML) of three ranges.
#
#   Rscript R/kdr_cv.R
#
# Writes outputs/kdr/kdr_cv_scores.csv, figures/kdr_cv_scores.png and
# figures/kdr_cv_<variant>.png.

n_folds <- 10
block_degrees <- 5
cores <- 4
set.seed(2026)

# the smoothing parameters (log10) of the fixed thin-plate candidates, which
# span about 64 to 3.5 effective degrees of freedom on these data
log10_sps <- seq(-0.5, 4, by = 0.5)

# the candidates: each a formula for the spatial smooth, and a fixed smoothing
# parameter (NULL for REML)
spatial_formula <- function(smooth) {
  as.formula(paste("cbind(successes, trials - successes) ~", smooth,
                   "+ year_c"))
}
tp_smooth <- "s(longitude, latitude, bs = 'tp', k = 80)"
gp_smooth <- function(range) {
  sprintf("s(longitude, latitude, bs = 'gp', k = 80, m = c(-4, %i))", range)
}
candidates <- c(
  map(log10_sps,
      \(lsp) list(formula = spatial_formula(tp_smooth), sp = 10 ^ lsp,
                  type = "thin plate, fixed")) %>%
    setNames(sprintf("tp_sp%+.1f", log10_sps)),
  list(tp_reml = list(formula = spatial_formula(tp_smooth), sp = NULL,
                      type = "thin plate, REML")),
  map(c(5, 10, 20),
      \(range) list(formula = spatial_formula(gp_smooth(range)), sp = NULL,
                    type = "Matern GP, REML")) %>%
    setNames(sprintf("gp%i", c(5, 10, 20)))
)

# fit a candidate to some records
fit_candidate <- function(candidate, data) {
  if (is.null(candidate$sp)) {
    gam(candidate$formula, family = quasibinomial(), data = data,
        method = "REML")
  } else {
    gam(candidate$formula, family = quasibinomial(), data = data,
        sp = candidate$sp)
  }
}

# the effective degrees of freedom of the spatial smooth (the total, less the
# intercept and the year trend)
spatial_edf <- function(model) {
  sum(model$edf) - 2
}

# binomial deviance of each record, with 0 log 0 = 0
binomial_deviance <- function(successes, trials, p) {
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  failures <- trials - successes
  2 * (if_else(successes > 0, successes * log(successes / (trials * p)), 0) +
         if_else(failures > 0, failures * log(failures / (trials * (1 - p))),
                 0))
}

# data and folds ---------------------------------------------------------------

records <- read_csv("data/clean/kdr_records.csv",
                    col_types = cols(variant = "c", .default = col_guess())) %>%
  mutate(year_c = year_start - 2015)

# The study of each record: its first author and publication year where we can
# find them, so that copies of a study held by two sources under different
# citations fall in the same fold, otherwise its citation.
studies <- kdr_study_lookup(records) %>%
  filter(!is.na(surname), surname != "", !is.na(publication_year)) %>%
  distinct(source, citation, .keep_all = TRUE) %>%
  transmute(source, citation, study = paste(surname, publication_year))
records <- records %>%
  left_join(studies, by = c("source", "citation")) %>%
  mutate(study = coalesce(study, citation))

# the records of a variant with a sample size, as allele counts
allele_data <- function(v) {
  records %>%
    filter(variant == v, !is.na(n_tested), n_tested > 0) %>%
    mutate(
      trials = 2 * n_tested,
      successes = round(frequency * trials),
      block = paste(floor(longitude / block_degrees),
                    floor(latitude / block_degrees))
    )
}

# assign groups to folds at random, in near-equal numbers of groups
group_folds <- function(groups) {
  unique_groups <- unique(groups)
  fold_of_group <- sample(rep_len(seq_len(n_folds), length(unique_groups)))
  fold_of_group[match(groups, unique_groups)]
}

# cross-validation -------------------------------------------------------------

# held-out predictions of every record by every candidate, for one fold
predict_fold <- function(fold, data, folds) {
  train <- data[folds != fold, ]
  test <- data[folds == fold, ]
  imap(candidates, \(candidate, name) {
    model <- tryCatch(fit_candidate(candidate, train),
                      error = function(e) NULL)
    p <- if (is.null(model)) {
      NA_real_
    } else {
      as.numeric(predict(model, test, type = "response"))
    }
    tibble(row = which(folds == fold), fold = fold, candidate = name, p = p)
  }) %>%
    bind_rows()
}

cv_variant <- function(v) {
  data <- allele_data(v)
  fold_sets <- list(study = group_folds(data$study),
                    block = group_folds(data$block))
  full_edf <- map_dbl(candidates, \(candidate) {
    spatial_edf(fit_candidate(candidate, data))
  })
  imap(fold_sets, \(folds, fold_type) {
    mclapply(seq_len(n_folds), predict_fold, data = data, folds = folds,
             mc.cores = cores) %>%
      bind_rows() %>%
      mutate(fold_type = fold_type)
  }) %>%
    bind_rows() %>%
    mutate(
      variant = v,
      successes = data$successes[row],
      trials = data$trials[row],
      deviance = binomial_deviance(successes, trials, p),
      edf = full_edf[candidate]
    )
}

held_out <- map(c("995F", "995S"), cv_variant) %>%
  bind_rows()

# scores -----------------------------------------------------------------------

# Per candidate: (i) the total held-out deviance, and (ii) the mean over
# records of the deviance per allele, so each record counts equally. For each,
# the difference from the best candidate, with a standard error from the
# spread of that difference across the folds.
fold_scores <- held_out %>%
  group_by(variant, fold_type, candidate, edf, fold) %>%
  summarise(
    total = sum(deviance),
    per_allele = sum(deviance / trials),
    records = n(),
    failed = any(is.na(p)),
    .groups = "drop"
  )

best <- fold_scores %>%
  group_by(variant, fold_type, candidate) %>%
  summarise(total = sum(total), per_allele = sum(per_allele),
            .groups = "drop_last") %>%
  summarise(best_total = candidate[which.min(total)],
            best_per_allele = candidate[which.min(per_allele)],
            .groups = "drop")

scores <- fold_scores %>%
  left_join(best, by = c("variant", "fold_type")) %>%
  group_by(variant, fold_type, fold) %>%
  mutate(
    total_fold_difference = total - total[candidate == best_total],
    per_allele_fold_difference = per_allele -
      per_allele[candidate == best_per_allele]
  ) %>%
  group_by(variant, fold_type, candidate, edf) %>%
  summarise(
    total_deviance = sum(total),
    total_difference = sum(total_fold_difference),
    total_difference_se = sd(total_fold_difference) * sqrt(n()),
    mean_deviance_per_allele = sum(per_allele) / sum(records),
    per_allele_difference = sum(per_allele_fold_difference) / sum(records),
    per_allele_difference_se = sd(per_allele_fold_difference) * sqrt(n()) /
      sum(records),
    failed_fits = sum(failed),
    .groups = "drop"
  ) %>%
  mutate(
    type = map_chr(candidate, \(name) candidates[[name]]$type)
  ) %>%
  arrange(variant, fold_type, desc(edf))

dir.create("outputs/kdr", showWarnings = FALSE, recursive = TRUE)
write.csv(scores, "outputs/kdr/kdr_cv_scores.csv", row.names = FALSE)

print(scores %>%
        select(variant, fold_type, candidate, edf, total_difference,
               total_difference_se, per_allele_difference,
               per_allele_difference_se, failed_fits) %>%
        mutate(across(where(is.numeric), \(x) signif(x, 3))) %>%
        as.data.frame())

# the candidate each fold type and score picks
picks <- best %>%
  pivot_longer(c(best_total, best_per_allele),
               names_to = "metric", values_to = "candidate") %>%
  mutate(
    metric = recode(metric,
                    best_total = "total deviance",
                    best_per_allele = "mean deviance per allele")
  ) %>%
  left_join(distinct(scores, variant, candidate, edf),
            by = c("variant", "candidate"))
print(as.data.frame(picks))

# score against edf
score_plot_data <- bind_rows(
  scores %>%
    transmute(variant, fold_type, candidate, edf, type,
              metric = "total deviance",
              difference = total_difference, se = total_difference_se),
  scores %>%
    transmute(variant, fold_type, candidate, edf, type,
              metric = "mean deviance per allele",
              difference = per_allele_difference,
              se = per_allele_difference_se)
) %>%
  mutate(
    panel = sprintf("%s, %s folds", variant, fold_type),
    metric = factor(metric, c("total deviance", "mean deviance per allele")),
    marker = case_when(
      type == "thin plate, fixed" ~ "thin plate, fixed sp",
      candidate == "tp_reml" ~ "thin plate, REML",
      .default = sprintf("Matern GP, REML, range %s degrees",
                         str_remove(candidate, "gp"))
    ),
    # spread the GP candidates, whose edf are close, a little along the axis
    edf_plot = edf * case_when(candidate == "gp5" ~ 0.97,
                               candidate == "gp20" ~ 1.03,
                               .default = 1)
  )

score_plot <- ggplot(score_plot_data,
                     aes(edf_plot, difference, colour = marker)) +
  geom_line(data = filter(score_plot_data, type == "thin plate, fixed"),
            colour = grey(0.6)) +
  geom_linerange(aes(ymin = difference - se, ymax = difference + se)) +
  geom_point(size = 1.8) +
  facet_grid(metric ~ panel, scales = "free_y") +
  scale_x_log10() +
  scale_colour_manual(
    values = c("thin plate, fixed sp" = grey(0.3),
               "thin plate, REML" = "#c0392b",
               "Matern GP, REML, range 5 degrees" = "#85c1e9",
               "Matern GP, REML, range 10 degrees" = "#2e86c1",
               "Matern GP, REML, range 20 degrees" = "#1b4f72")
  ) +
  labs(x = "effective degrees of freedom of the spatial smooth (full data)",
       y = "held-out score, difference from the best (+/- 1 se)",
       colour = NULL,
       title = "Cross-validated smoothness of the kdr surfaces") +
  guides(colour = guide_legend(nrow = 2)) +
  theme_bw(base_size = 9) +
  theme(legend.position = "top", strip.background = element_blank())

ggsave("figures/kdr_cv_scores.png", score_plot, width = 11, height = 6,
       dpi = 150, bg = "white")

# maps -------------------------------------------------------------------------

# The predictions at 2015 of the candidates picked for each variant, over the
# whole mask (no masking, so the extrapolation is visible), and their standard
# errors on the logit scale, with the record locations.
mask <- rast("data/clean/raster_mask.tif")
grid <- as.data.frame(terra::aggregate(mask, 4), xy = TRUE, na.rm = TRUE) %>%
  transmute(longitude = x, latitude = y, year_c = 0)

for (v in c("995F", "995S")) {
  data <- allele_data(v)
  variant_picks <- picks %>%
    filter(variant == v) %>%
    group_by(candidate, edf) %>%
    summarise(picked_by = paste(sprintf("%s folds, %s", fold_type, metric),
                                collapse = "\n"),
              .groups = "drop") %>%
    arrange(desc(edf)) %>%
    mutate(panel = sprintf("%s (edf %.0f)\npicked by %s", candidate, edf,
                           picked_by))
  maps <- pmap(variant_picks, \(candidate, panel, ...) {
    model <- fit_candidate(candidates[[candidate]], data)
    prediction <- predict(model, grid, se.fit = TRUE)
    grid %>%
      mutate(panel = panel,
             mean = plogis(as.numeric(prediction$fit)),
             se = as.numeric(prediction$se.fit))
  }) %>%
    bind_rows() %>%
    mutate(panel = factor(panel, variant_picks$panel))
  points <- distinct(data, longitude, latitude)
  base <- function(fill) {
    ggplot(maps, aes(longitude, latitude, fill = .data[[fill]])) +
      geom_raster() +
      geom_point(data = points, aes(longitude, latitude), inherit.aes = FALSE,
                 size = 0.15, colour = "grey50", alpha = 0.5) +
      facet_wrap(~ panel, nrow = 1) +
      coord_equal() +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_size = 9) +
      theme(panel.grid = element_blank(), axis.text = element_blank())
  }
  top <- base("mean") +
    scale_fill_viridis_c(limits = c(0, 1)) +
    labs(fill = "2015", title = sprintf("%s allele frequency, 2015", v))
  bottom <- base("se") +
    scale_fill_viridis_c(option = "magma", limits = c(0, 3),
                         oob = scales::squish) +
    labs(fill = "se (logit)", title = "standard error, logit scale")
  ggsave(sprintf("figures/kdr_cv_%s.png", v),
         patchwork::wrap_plots(top, bottom, ncol = 1),
         width = 3 + 3.5 * nrow(variant_picks), height = 8, dpi = 120,
         bg = "white")
}
