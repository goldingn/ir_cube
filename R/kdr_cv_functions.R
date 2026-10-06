# functions for cross-validating the kdr allele frequency surfaces (#46),
# used by R/kdr_cv.R and R/kdr_cv_extrap.R

suppressMessages(library(mgcv))
library(parallel)

# folds: studies, or 5 x 5 degree blocks, assigned at random to n_folds folds
n_folds <- 10
block_degrees <- 5
fold_seed <- 2026
cores <- 4
variants <- c("995F", "995S")

# Candidate surfaces are binomial GAMs of allele counts with a spatial smooth
# and a linear year trend. Each candidate is a list of the model formula, the
# smoothing parameters (NULL for REML), whether to penalise the smooth's null
# space (`select`), and a label for its type.
spatial_formula <- function(smooth) {
  as.formula(paste("cbind(successes, trials - successes) ~", smooth,
                   "+ year_c"))
}
tp_smooth <- "s(longitude, latitude, bs = 'tp', k = 80)"
gp_smooth <- function(range) {
  sprintf("s(longitude, latitude, bs = 'gp', k = 80, m = c(-4, %i))", range)
}
candidate <- function(smooth, sp = NULL, select = FALSE, type) {
  list(formula = spatial_formula(smooth), sp = sp, select = select,
       type = type)
}

# fit a candidate to some records
fit_candidate <- function(candidate, data) {
  if (is.null(candidate$sp)) {
    gam(candidate$formula, family = quasibinomial(), data = data,
        method = "REML", select = candidate$select)
  } else {
    gam(candidate$formula, family = quasibinomial(), data = data,
        sp = candidate$sp, select = candidate$select)
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

# The kdr records, with the study of each: its first author and publication
# year where we can find them, so that copies of a study held by two sources
# under different citations fall in the same fold, otherwise its citation.
kdr_cv_records <- function() {
  records <- read_csv("data/clean/kdr_records.csv",
                      col_types = cols(variant = "c",
                                       .default = col_guess())) %>%
    mutate(year_c = year_start - 2015)
  studies <- kdr_study_lookup(records) %>%
    filter(!is.na(surname), surname != "", !is.na(publication_year)) %>%
    distinct(source, citation, .keep_all = TRUE) %>%
    transmute(source, citation, study = paste(surname, publication_year))
  records %>%
    left_join(studies, by = c("source", "citation")) %>%
    mutate(study = coalesce(study, citation))
}

# the records of a variant with a sample size, as allele counts
allele_data <- function(records, v) {
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

# the study and block folds of each variant, the same in every script
kdr_cv_folds <- function(records) {
  set.seed(fold_seed)
  map(setNames(variants, variants), \(v) {
    data <- allele_data(records, v)
    list(study = group_folds(data$study),
         block = group_folds(data$block))
  })
}

# cross-validation -------------------------------------------------------------

# Held-out predictions of every record in one fold by every candidate. The
# fixed-smoothing candidates sharing a model formula share one model setup.
predict_fold <- function(fold, data, folds, candidates) {
  train <- data[folds != fold, ]
  test <- data[folds == fold, ]
  setups <- list()
  imap(candidates, \(candidate, name) {
    model <- tryCatch({
      if (is.null(candidate$sp)) {
        fit_candidate(candidate, train)
      } else {
        key <- paste(c(deparse(candidate$formula), candidate$select),
                     collapse = " ")
        if (is.null(setups[[key]])) {
          setups[[key]] <<- gam(candidate$formula, family = quasibinomial(),
                                data = train, select = candidate$select,
                                fit = FALSE)
        }
        gam(G = setups[[key]], sp = candidate$sp)
      }
    }, error = function(e) NULL)
    p <- if (is.null(model)) {
      NA_real_
    } else {
      as.numeric(predict(model, test, type = "response"))
    }
    tibble(row = which(folds == fold), fold = fold, candidate = name, p = p)
  }) %>%
    bind_rows()
}

# held-out predictions and deviances of every record of a variant, by every
# candidate and fold type, with each candidate's edf on the full data
cv_variant <- function(v, records, folds, candidates) {
  data <- allele_data(records, v)
  full_edf <- map_dbl(candidates, \(candidate) {
    spatial_edf(fit_candidate(candidate, data))
  })
  imap(folds[[v]], \(fold_set, fold_type) {
    mclapply(seq_len(n_folds), predict_fold, data = data, folds = fold_set,
             candidates = candidates, mc.cores = cores) %>%
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

# Per candidate: (i) the total held-out deviance, and (ii) the mean over
# records of the deviance per allele, so each record counts equally. For each,
# the difference from a comparison candidate - the best, or a reference
# candidate per variant (a named vector) - with a standard error from the
# spread of that difference across the folds.
score_candidates <- function(held_out, candidates, reference = NULL) {
  fold_scores <- held_out %>%
    group_by(variant, fold_type, candidate, edf, fold) %>%
    summarise(
      total = sum(deviance),
      per_allele = sum(deviance / trials),
      records = n(),
      failed = any(is.na(p)),
      .groups = "drop"
    )
  comparison <- fold_scores %>%
    group_by(variant, fold_type, candidate) %>%
    summarise(total = sum(total), per_allele = sum(per_allele),
              .groups = "drop_last") %>%
    summarise(compare_total = candidate[which.min(total)],
              compare_per_allele = candidate[which.min(per_allele)],
              .groups = "drop")
  if (!is.null(reference)) {
    comparison <- comparison %>%
      mutate(compare_total = reference[variant],
             compare_per_allele = reference[variant])
  }
  fold_scores %>%
    left_join(comparison, by = c("variant", "fold_type")) %>%
    group_by(variant, fold_type, fold) %>%
    mutate(
      total_fold_difference = total - total[candidate == compare_total],
      per_allele_fold_difference = per_allele -
        per_allele[candidate == compare_per_allele]
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
}

# maps -------------------------------------------------------------------------

# the prediction grid: the whole mask, at a quarter of its resolution, in 2015
kdr_prediction_grid <- function() {
  mask <- rast("data/clean/raster_mask.tif")
  as.data.frame(terra::aggregate(mask, 4), xy = TRUE, na.rm = TRUE) %>%
    transmute(longitude = x, latitude = y, year_c = 0)
}

# the 2015 prediction of each of `panels` (a table of candidate names and
# panel titles) fitted to all the records of a variant, and its standard error
# on the logit scale
predict_candidates <- function(data, panels, candidates, grid) {
  pmap(panels, \(candidate, panel, ...) {
    model <- fit_candidate(candidates[[candidate]], data)
    prediction <- predict(model, grid, se.fit = TRUE)
    grid %>%
      mutate(candidate = candidate,
             panel = panel,
             mean = plogis(as.numeric(prediction$fit)),
             se = as.numeric(prediction$se.fit))
  }) %>%
    bind_rows() %>%
    mutate(panel = factor(panel, panels$panel))
}

# map predictions over the whole mask (no masking, so the extrapolation is
# visible), and their standard errors, with the record locations
map_predictions <- function(maps, data, v, file) {
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
  ggsave(file, patchwork::wrap_plots(top, bottom, ncol = 1),
         width = 3 + 3.5 * n_distinct(maps$panel), height = 8, dpi = 120,
         bg = "white")
}
