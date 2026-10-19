# Bioassay or map: which better estimates mortality at a pixel-year with data?
# (#30)
#
#   Rscript R/bioassay_vs_map.R
#
# Run after R/validation_metrics.R, which writes the per-record scores read
# here (outputs/cv_scores.csv). Nothing is refitted or rescored.
#
# The quantity is n*: the number of bioassays from the same pixel-year whose
# mean predicts a held-out bioassay as accurately as the map. On held-out
# pixel-years (cell x year x insecticide) with at least two assays,
#
#   D_b = mean squared difference between an assay and its partners at the
#         pixel-year, averaged per held-out assay
#   D_m = mean squared difference between the map and the held-out assay
#   W   = D_b / 2, one assay's noise variance
#   n*  = W / (D_m - W)
#
# The held-out assay's noise is in both D_b and D_m, so no rho enters (see
# excess_ratio() in R/bioassay_vs_map_functions.R). The map is the two-stage
# model's point prediction (`map`, #32): the posterior mean of
# plogis(m + omega + xi), without u and p.
#
# Writes, rows by CV experiment x insecticide class (plus all classes):
#
#   outputs/cv_bioassay_vs_map_independent.csv (Table 1)
#       n* from every held-out assay with a partner, with a bootstrap over
#       pixels: held-out pixel-years concentrate in few pixels, and those at one
#       pixel share the map's error there. Alongside the pair estimate of one
#       assay's noise (one_bioassay = D_b / 2), the data floor from rho-hat on
#       the same assays (rho_floor), which checks the floor #36 uses (#32).
#   outputs/cv_bioassay_vs_map_surveys.csv (Table 2, a robustness check)
#       n* from different-survey partners only, with a bootstrap over clusters
#       of pixel-years linked by a shared pixel or survey. A partner from the
#       held-out assay's own survey shares its survey error, which favours the
#       bioassay. Different-survey pairs are also more often both at 0% or 100%,
#       so a difference between the tables is not the survey effect alone.
#
# Rows with few pixels (Table 1) or few clusters (Table 2) are flagged in
# `note`; cluster-bootstrap intervals are unreliable with few clusters
# (Cameron & Miller 2015, J Human Resources 50:317-372).

source("R/validation_scoring.R")
source("R/bioassay_vs_map_functions.R")

suppressMessages(library(readr))

set.seed(2026 - 10 - 5)
n_bootstrap <- 2000

# flag cells below these sizes: Table 1 by pixels, Table 2 by clusters
min_pixels <- 10
min_clusters <- 20


# held-out records, with surveys and replicate partners -------------------------

# the two-stage model's per-record scores: the map and the data floor (#32)
records <- read_csv("outputs/cv_scores.csv", show_col_types = FALSE,
                    col_select = c(model, experiment, fold, insecticide_type,
                                   insecticide_class, country_name, citation,
                                   year_start, cell, observed, map,
                                   floor)) %>%
  filter(model == "two_stage") %>%
  mutate(survey = survey_id(citation, country_name, year_start),
         map_error = (observed - map) ^ 2) %>%
  group_by(experiment, cell, year_start, insecticide_type) %>%
  mutate(pixel_year = cur_group_id(), n_pixel_year = n()) %>%
  ungroup()

# every held-out assay with a partner at its pixel-year
assays <- records %>%
  filter(n_pixel_year >= 2) %>%
  group_by(pixel_year) %>%
  mutate(replicate_partners(observed, survey)) %>%
  ungroup() %>%
  # only what the tables use, which keeps the bootstraps fast
  select(experiment, insecticide_class, cell, pixel_year, survey, floor,
         map_error, pair_any, pair_other_survey, partners_other_survey)

# the rows of every table: each class, and all classes
by_class <- function(data) {
  bind_rows(data, mutate(data, insecticide_class = "all"))
}
class_order <- c(sort(unique(records$insecticide_class)), "all")


# Table 1: all partners, pixel bootstrap ----------------------------------------

table_independent <- assays %>%
  by_class() %>%
  group_by(experiment, insecticide_class) %>%
  group_modify(function(data, key) {
    ratio <- function(d) c(data = excess_ratio(d$pair_any, d$map_error))
    drawn <- cluster_bootstrap(data, ratio, n_bootstrap)
    w <- mean(data$pair_any) / 2
    bind_cols(
      tibble(
        pixel_years = n_distinct(data$pixel_year),
        pixels = n_distinct(data$cell),
        assays = nrow(data),
        one_bioassay = w,
        rho_floor = mean(data$floor, na.rm = TRUE),
        map_error = mean(data$map_error) - w),
      as_tibble(as.list(n_star_interval(ratio(data)[["data"]],
                                        drawn[, "data"])))
    )
  }) %>%
  ungroup() %>%
  mutate(note = ifelse(pixels < min_pixels,
                       sprintf("fewer than %i pixels", min_pixels), ""),
         insecticide_class = factor(insecticide_class, class_order)) %>%
  arrange(experiment, insecticide_class)

write.csv(table_independent, "outputs/cv_bioassay_vs_map_independent.csv",
          row.names = FALSE)


# Table 2: different-survey partners, cluster bootstrap --------------------------

table_surveys <- assays %>%
  filter(!is.na(pair_other_survey)) %>%
  by_class() %>%
  group_by(experiment, insecticide_class) %>%
  group_modify(function(data, key) {
    cluster <- survey_pixel_clusters(data)
    ratio <- function(d) c(data = excess_ratio(d$pair_other_survey,
                                               d$map_error))
    drawn <- cluster_bootstrap(data, ratio, n_bootstrap, cluster = cluster)
    w <- mean(data$pair_other_survey) / 2
    bind_cols(
      tibble(
        pairs = sum(data$partners_other_survey) / 2,
        surveys = n_distinct(data$survey),
        clusters = n_distinct(cluster),
        pixel_years = n_distinct(data$pixel_year),
        pixels = n_distinct(data$cell),
        assays = nrow(data),
        one_bioassay = w,
        map_error = mean(data$map_error) - w),
      as_tibble(as.list(n_star_interval(ratio(data)[["data"]],
                                        drawn[, "data"])))
    )
  }) %>%
  ungroup() %>%
  mutate(note = ifelse(clusters < min_clusters,
                       sprintf("fewer than %i clusters: interval not interpreted",
                               min_clusters), ""),
         insecticide_class = factor(insecticide_class, class_order)) %>%
  arrange(experiment, insecticide_class)

write.csv(table_surveys, "outputs/cv_bioassay_vs_map_surveys.csv",
          row.names = FALSE)


# print ---------------------------------------------------------------------------

interval <- function(estimate, lower, upper, digits = 2) {
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
          estimate, lower, upper)
}

cat("\nTable 1, all partners (pixel bootstrap):\n")
print(as.data.frame(table_independent %>%
  transmute(experiment, class = insecticide_class, pixel_years, pixels,
            one_bioassay = round(one_bioassay, 3),
            rho_floor = round(rho_floor, 3),
            map_error = round(map_error, 3),
            n_star = interval(n_star, n_star_lower, n_star_upper),
            note)), row.names = FALSE)

cat("\nTable 2, different-survey partners (cluster bootstrap):\n")
print(as.data.frame(table_surveys %>%
  transmute(experiment, class = insecticide_class, pairs, surveys,
            clusters, one_bioassay = round(one_bioassay, 3),
            map_error = round(map_error, 3),
            n_star = interval(n_star, n_star_lower, n_star_upper),
            note)), row.names = FALSE)
