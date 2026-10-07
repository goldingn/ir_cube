# Summaries of every fetched species run (#47; R/species_runs.R) side by side.
#
#   Rscript R/species_fit_summaries.R [<label>=<fitted_model.RData> ...]
#
# Finds the runs fetched with `irpod fetch <name>`
# (outputs/pod_jobs/<name>/temporary/fitted_model.RData), labelled as in
# species_runs (B_f, V1, ...), and any other fits given as label=file (e.g. a
# reference fit). Summarises each that has no summary newer than its fit with
# R/species_fit_summary.R, one process at a time, and writes
#   outputs/species_runs/summary.csv          parameters x fits: posterior
#                                             summaries, R-hat, ESS
#   outputs/species_runs/summary_by_mode.csv  the new parameters within each
#                                             floor mode
#   outputs/species_runs/summary_chains.csv   one row per chain
#   outputs/species_runs/summary_fits.csv     one row per fit
#   figures/species_runs/new_parameters.png   posteriors of the new
#                                             parameters by fit (usable
#                                             chains): exp(gamma_*),
#                                             exp(delta_*), the floors, the
#                                             latent smooths' sd and range,
#                                             and shear loading
# and, for the fits R/species_misfit.R has been run on,
#   outputs/species_runs/misfit_regions.csv   their region tables together
#   outputs/species_runs/loo_compare.csv      loo::loo_compare() of their
#                                             PSIS-LOO (bioassays are
#                                             clustered: the se of the
#                                             differences is optimistic)

suppressMessages({
  library(dplyr)
  library(ggplot2)
})
source("R/species_runs.R")

summary_dir <- "outputs/species_runs/summary"
figure_dir <- "figures/species_runs"
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

fits <- setNames(
  file.path("outputs/pod_jobs", species_runs$name,
            "temporary/fitted_model.RData"),
  species_runs$label)
for (argument in commandArgs(trailingOnly = TRUE)) {
  parts <- strsplit(argument, "=", fixed = TRUE)[[1]]
  stopifnot(length(parts) == 2)
  fits[[parts[1]]] <- parts[2]
}
found <- file.exists(fits)
cat(sprintf("%d of %d fits found: %s\n", sum(found), length(fits),
            toString(names(fits)[found])))
fits <- fits[found]
if (length(fits) == 0) {
  quit(save = "no")
}

for (label in names(fits)) {
  out <- file.path(summary_dir, sprintf("%s.rds", label))
  if (file.exists(out) && file.mtime(out) > file.mtime(fits[[label]])) {
    next
  }
  status <- system2("Rscript", c("R/species_fit_summary.R", fits[[label]],
                                 label))
  if (status != 0) {
    warning("the summary of ", label, " failed")
  }
}

summary_files <- file.path(summary_dir, sprintf("%s.rds", names(fits)))
fits <- fits[file.exists(summary_files)]
summaries <- lapply(summary_files[file.exists(summary_files)], readRDS)
combine <- function(part) bind_rows(lapply(summaries, `[[`, part))
parameters <- combine("parameters")
write.csv(parameters, "outputs/species_runs/summary.csv", row.names = FALSE)
write.csv(combine("parameters_by_mode"),
          "outputs/species_runs/summary_by_mode.csv", row.names = FALSE)
write.csv(combine("chains"), "outputs/species_runs/summary_chains.csv",
          row.names = FALSE)
fit_rows <- combine("fit")
write.csv(fit_rows, "outputs/species_runs/summary_fits.csv",
          row.names = FALSE)
options(width = 160)
print(as.data.frame(select(fit_rows, -file, -options)), digits = 4)

# the new parameters: medians, 50% and 95% intervals, by fit
new <- parameters %>%
  filter(group %in% c("multiplier", "floor", "smooth"),
         !grepl("^floor_intercept", parameter)) %>%
  mutate(label = factor(label, levels = rev(names(fits))))
if (nrow(new) > 0) {
  # no effect, 1, for the multipliers (an empty layer breaks the facets)
  multipliers <- distinct(filter(new, group == "multiplier"), parameter)
  p <- ggplot(new, aes(y = label)) +
    (if (nrow(multipliers) > 0) {
      geom_vline(aes(xintercept = 1), colour = "grey60", linetype = 2,
                 data = multipliers)
    }) +
    geom_linerange(aes(xmin = q2.5, xmax = q97.5), linewidth = 0.4) +
    geom_linerange(aes(xmin = q25, xmax = q75), linewidth = 1.4) +
    geom_point(aes(x = q50), size = 2) +
    facet_wrap(~ parameter, scales = "free_x") +
    labs(x = "posterior median, 50% and 95% intervals (usable chains)",
         y = NULL,
         caption = paste("exp(): multipliers of the cumulative log fitness",
                         "and the fitness cost; 1, dashed, is no effect")) +
    theme_bw(base_size = 9) +
    theme(strip.background = element_blank())
  ggsave(file.path(figure_dir, "new_parameters.png"), p, width = 9,
         height = 2 + 0.9 * ceiling(length(unique(new$parameter)) / 3),
         dpi = 150)
}
cat("written outputs/species_runs/summary*.csv and",
    file.path(figure_dir, "new_parameters.png"), "\n")

# the misfit tables and PSIS-LOO of the fits R/species_misfit.R has been run on
misfit_dir <- "outputs/species_runs/misfit"
region_files <- file.path(misfit_dir, sprintf("%s_regions.csv", names(fits)))
if (any(file.exists(region_files))) {
  write.csv(bind_rows(lapply(region_files[file.exists(region_files)],
                             read.csv)),
            "outputs/species_runs/misfit_regions.csv", row.names = FALSE)
}
loo_files <- setNames(file.path(misfit_dir, sprintf("%s_loo.rds", names(fits))),
                      names(fits))
loo_files <- loo_files[file.exists(loo_files)]
if (length(loo_files) >= 2) {
  loos <- lapply(loo_files, readRDS)
  comparison <- as.data.frame(unclass(loo::loo_compare(loos)))
  # loo_compare() names the rows model1, model2, ... or by position; the
  # label is matched by elpd_loo
  elpd <- vapply(loos, function(x) x$estimates["elpd_loo", "Estimate"],
                 numeric(1))
  comparison <- data.frame(label = names(loos)[match(comparison$elpd_loo,
                                                     elpd)],
                           comparison, row.names = NULL)
  write.csv(comparison, "outputs/species_runs/loo_compare.csv",
            row.names = FALSE)
  cat("\nPSIS-LOO, best first (bioassays are clustered: the se of the",
      "differences is optimistic):\n")
  print(comparison[, c("label", "elpd_diff", "se_diff", "elpd_loo",
                       "p_loo")], digits = 6)
}
