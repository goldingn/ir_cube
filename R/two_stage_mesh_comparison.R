# Compare the two-stage correction across mesh configurations (#21,
# doc/two_stage_plan.md, "Mesh resolution results").
#
#   Rscript R/two_stage_mesh_comparison.R
#
# Reads the per-record scores written by R/two_stage_metrics.R
# (outputs/two_stage/cv_scores_two_stage.csv) and the fit summary, and reports
# each mesh configuration of omega_xi_u (models two_stage_omega_xi_u_mesh-<tag>)
# against the base meshes (two_stage_omega_xi_u), on the records both were
# scored on (and some finer configurations against each other): mean log score, CRPS and variance explained, with paired
# pixel-bootstrap intervals of the difference from base, as in
# two_stage_metrics.R. Coverage is read from cv_summary_two_stage.csv. Writes
# outputs/two_stage/mesh_comparison.csv and mesh_fit_summary.csv.

suppressMessages({
  library(dplyr)
  library(tidyr)
})

output_dir <- "outputs/two_stage"
base_model <- "two_stage_omega_xi_u"
n_bootstrap <- 2000

scores <- read.csv(file.path(output_dir, "cv_scores_two_stage.csv"),
                   stringsAsFactors = FALSE) %>%
  filter(model == base_model | startsWith(model, paste0(base_model, "_mesh-")) |
           model == "dynamical") %>%
  mutate(fold = as.character(fold))
summary_table <- read.csv(file.path(output_dir, "cv_summary_two_stage.csv"),
                          stringsAsFactors = FALSE)

compare <- function(data, candidate, reference_model = base_model) {
  # the folds on which both the candidate and the reference were run
  folds <- intersect(unique(data$fold[data$model == candidate]),
                     unique(data$fold[data$model == reference_model]))
  pair <- lapply(c(reference_model, candidate), function(m) {
    data %>% filter(model == m, fold %in% folds) %>% arrange(fold, row)
  })
  stopifnot(identical(pair[[1]]$row, pair[[2]]$row),
            identical(pair[[1]]$fold, pair[[2]]$fold))
  reference <- pair[[1]]
  cells <- sort(unique(reference$cell))
  index <- match(reference$cell, cells)
  per_cell <- function(x) as.vector(rowsum(x, index, reorder = TRUE))
  common <- cbind(n = per_cell(rep(1, nrow(reference))),
                  y = per_cell(reference$observed),
                  y2 = per_cell(reference$observed ^ 2))
  stats <- lapply(pair, function(d) cbind(log = per_cell(d$log_score),
                                          crps = per_cell(d$crps),
                                          se = per_cell(d$squared_error)))
  metrics <- function(w) {
    totals <- crossprod(w, common)
    n <- totals[, "n"]
    variance <- totals[, "y2"] / n - (totals[, "y"] / n) ^ 2
    lapply(stats, function(s) {
      t <- crossprod(w, s) / n
      cbind(elpd = t[, "log"], crps = t[, "crps"],
            explained = 100 * (1 - t[, "se"] / variance))
    })
  }
  point <- metrics(matrix(1, length(cells), 1))
  boot <- metrics(rmultinom(n_bootstrap, length(cells), rep(1, length(cells))))
  out <- tibble(model = candidate, reference_model = reference_model,
                folds = paste(folds, collapse = "+"),
                n = nrow(reference), n_pixels = length(cells))
  for (metric in c("elpd", "crps", "explained")) {
    out[[metric]] <- point[[2]][1, metric]
    out[[paste0(metric, "_ref")]] <- point[[1]][1, metric]
    difference <- boot[[2]][, metric] - boot[[1]][, metric]
    out[[paste0("diff_", metric)]] <- point[[2]][1, metric] -
      point[[1]][1, metric]
    out[[paste0("diff_", metric, "_lower")]] <- quantile(difference, 0.025)
    out[[paste0("diff_", metric, "_upper")]] <- quantile(difference, 0.975)
    out[[paste0("prob_better_", metric)]] <- mean(
      if (metric == "crps") difference < 0 else difference > 0)
  }
  out
}

set.seed(2026 - 9 - 25)
candidates <- setdiff(unique(scores$model), c(base_model, "dynamical"))
comparison <- bind_rows(lapply(split(scores, scores$experiment), function(d) {
  present <- intersect(candidates, unique(d$model))
  pairs <- tibble(candidate = present, reference = base_model)
  # the finer omega mesh on top of the finer xi mesh: omega's effect alone
  xi1200 <- paste0(base_model, "_mesh-xi1200")
  omega5000 <- paste0(base_model, "_mesh-omega5000")
  if (all(c(xi1200, omega5000) %in% present)) {
    pairs <- add_row(pairs, candidate = omega5000, reference = xi1200)
  }
  # xi on the omega mesh against the intermediate xi mesh
  xifull <- paste0(base_model, "_mesh-xifull")
  if (all(c(xi1200, xifull) %in% present)) {
    pairs <- add_row(pairs, candidate = xifull, reference = xi1200)
  }
  # both meshes refined, against each refined alone
  both <- paste0(base_model, "_mesh-omega5000_xi2500")
  for (single in c(omega5000, xifull)) {
    if (all(c(both, single) %in% present)) {
      pairs <- add_row(pairs, candidate = both, reference = single)
    }
  }
  bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
    compare(d, pairs$candidate[i], pairs$reference[i]) %>%
      mutate(experiment = d$experiment[1], .before = 1)
  }))
}))

if (nrow(comparison) == 0) stop("no mesh-configuration runs scored yet")

coverage <- summary_table %>%
  filter(fold == "pooled", stratum == "all",
         model %in% c(base_model, candidates)) %>%
  select(experiment, model, coverage_50, coverage_95)
# coverage is of the pooled experiment; for a candidate missing a fold it is
# over the folds it has, since the pooled summary requires every fold
comparison <- comparison %>%
  left_join(coverage, by = c("experiment", "model")) %>%
  left_join(coverage %>% filter(model == base_model) %>%
              select(experiment, coverage_50_base = coverage_50,
                     coverage_95_base = coverage_95),
            by = "experiment")

write.csv(comparison, file.path(output_dir, "mesh_comparison.csv"),
          row.names = FALSE)

cat("\neach mesh configuration of omega_xi_u against a reference, base unless stated (paired pixel bootstrap):\n")
print(as.data.frame(comparison %>%
  transmute(experiment, model = sub("^.*_mesh-", "", model),
            vs = ifelse(reference_model == base_model, "base",
                        sub("^.*_mesh-", "", reference_model)), folds, n,
            elpd = sprintf("%.4f (ref %.4f)", elpd, elpd_ref),
            d_elpd = sprintf("%+.4f [%+.4f, %+.4f]", diff_elpd,
                             diff_elpd_lower, diff_elpd_upper),
            d_crps = sprintf("%+.5f [%+.5f, %+.5f]", diff_crps,
                             diff_crps_lower, diff_crps_upper),
            explained = sprintf("%.1f (ref %.1f)", explained, explained_ref),
            d_explained = sprintf("%+.2f [%+.2f, %+.2f]", diff_explained,
                                  diff_explained_lower, diff_explained_upper),
            cov50 = round(coverage_50, 3), cov95 = round(coverage_95, 3))),
  row.names = FALSE)

# fitted ranges against the mesh edges, and cost, per configuration
fits <- read.csv(file.path(output_dir, "fit_summary.csv"),
                 stringsAsFactors = FALSE)
if (!"mesh_config" %in% names(fits)) fits$mesh_config <- NA_character_
fits <- fits %>%
  mutate(mesh = ifelse(is.na(mesh_config), "base", mesh_config)) %>%
  filter(variant == "omega_xi_u", m_ref == "mean")
fit_table <- fits %>%
  group_by(experiment, fold, mesh) %>%
  summarise(types = n(),
            omega_nodes = paste(range(mesh_nodes), collapse = "-"),
            omega_edge_km = paste(round(range(mesh_edge_km)), collapse = "-"),
            xi_nodes = paste(range(mesh_xi_nodes), collapse = "-"),
            xi_edge_km = paste(round(range(mesh_xi_edge_km)), collapse = "-"),
            omega_below_2_edges = sum(flag_range_omega == "below 2 mesh edges",
                                      na.rm = TRUE),
            eta_below_2_edges = sum(flag_range_eta == "below 2 xi-mesh edges",
                                    na.rm = TRUE),
            not_converged = sum(convergence != 0 | is.na(convergence)),
            fit_minutes = round(sum(time_fit_s) / 60, 1),
            max_fit_minutes = round(max(time_fit_s) / 60, 1),
            peak_memory_gb = round(max(peak_memory_gb), 1),
            .groups = "drop")
write.csv(fit_table, file.path(output_dir, "mesh_fit_summary.csv"),
          row.names = FALSE)
cat("\nmeshes, range flags and cost per fold and configuration:\n")
print(as.data.frame(fit_table), row.names = FALSE)

cat("\nfitted ranges (km) by type:\n")
print(as.data.frame(fits %>%
  filter(experiment == "spatial_interpolation") %>%
  select(insecticide_type, mesh, range_omega, range_eta, mesh_edge_km,
         mesh_xi_edge_km) %>%
  mutate(across(where(is.numeric), round)) %>%
  pivot_wider(names_from = mesh,
              values_from = c(range_omega, range_eta, mesh_edge_km,
                              mesh_xi_edge_km))), row.names = FALSE)
