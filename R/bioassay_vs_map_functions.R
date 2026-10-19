# Functions behind R/bioassay_vs_map.R (#30): survey identifiers, replicate
# pairs within held-out pixel-years, the clusters for the survey-structured
# bootstrap, and the n* summaries. Functions only: source from the repo root,
# after R/validation_scoring.R (for cluster_bootstrap()).

suppressMessages(library(dplyr))


# surveys -----------------------------------------------------------------------

# Survey identifier: citation x country x year. Restored from tag
# two-stage-options-6edcf0d (R/two_stage_correction.R), where it defined the
# survey effect removed in #29. The citation alone does not identify a study:
# aggregated sources ("Ministry of Health", "PMI 2016", "VectorBase", personal
# communications) span countries, so the country splits them (110 of 1434
# citation-years span more than one country, with 9,902 assays). Assays with no
# citation (12 in the cleaned data, all Madagascar 2023) are one survey per
# country-year.
survey_id <- function(citation, country, year) {
  paste(citation, country, year, sep = " | ")
}


# replicate pairs ---------------------------------------------------------------

# Per held-out assay in a pixel-year (cell x year x insecticide) with at least
# one other held-out assay: the mean squared difference to its partners, over
# all of them (`pair_any`) and over those from a different survey only
# (`pair_other_survey`, NA without one), and the number of each.
#
# Averaging over each assay's partners, rather than over all pairs, weights
# every held-out assay once, as the map's squared error does. Averaging over
# pairs would weight a pixel-year of n assays by n^2 in D_b but by n in D_m.
replicate_partners <- function(y, survey) {
  difference <- outer(y, y, "-") ^ 2
  other <- outer(survey, survey, "!=")
  diag(difference) <- NA
  partners_other <- rowSums(other)
  data.frame(
    partners = length(y) - 1,
    pair_any = rowMeans(difference, na.rm = TRUE),
    partners_other_survey = partners_other,
    pair_other_survey = ifelse(partners_other > 0,
                               rowSums(difference * other, na.rm = TRUE) /
                                 partners_other,
                               NA_real_)
  )
}


# clusters ----------------------------------------------------------------------

# Connected components of a bipartite graph: items (here pixel-years) linked
# through shared labels (pixels and surveys). Union-find on the labels, with
# path halving; returns a component index per item.
link_components <- function(item, label) {
  label <- as.integer(factor(label))
  parent <- seq_len(max(label))
  find <- function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i <- parent[i]
    }
    i
  }
  # each item's labels are joined to the item's first label
  for (rows in split(seq_along(item), item)) {
    root <- find(label[rows[1]])
    for (j in rows[-1]) {
      other <- find(label[j])
      if (other != root) parent[other] <- root
    }
  }
  roots <- vapply(seq_along(parent), find, integer(1))
  first_label <- label[match(unique(item), item)]
  setNames(as.integer(factor(roots[first_label])), unique(item))
}

# Clusters for the survey-structured bootstrap (Table 2): pixel-years that share
# a pixel or a survey are in one cluster. Pixel-years at one pixel share the
# map's error there, and the two assays of a pair share their surveys' errors
# with every other assay of those surveys. `assays` holds both assays of every
# pair (an assay with a different-survey partner is that partner's partner), so
# its surveys are all the pairs' surveys. Returns a cluster per row of `assays`.
survey_pixel_clusters <- function(assays) {
  links <- bind_rows(
    assays %>% distinct(pixel_year, label = paste("pixel", cell)),
    assays %>% distinct(pixel_year, label = paste("survey", survey))
  )
  components <- link_components(links$pixel_year, links$label)
  unname(components[as.character(assays$pixel_year)])
}


# n* --------------------------------------------------------------------------

# R, the map's error variance in population mortality per unit of one
# bioassay's noise variance W:
#
#   W = D_b / 2,   R = (D_m - W) / W,   n* = 1 / R
#
# with D_b the mean squared difference between two assays at a pixel-year and
# D_m the map's mean squared error against the held-out assay. The held-out
# assay's own noise W is in both, so it cancels without an estimate of rho.
# n pooled bioassays have noise variance W / n, so n* is the number whose mean
# predicts a held-out assay as accurately as the map.
#
# The interval comes from bootstrapping R, whose denominator is always
# positive, and inverting: n* in [1 / R_upper, 1 / R_lower], with an infinite
# upper limit when R_lower <= 0. Bootstrapping n* directly fails in small
# cells, where the map's net error can be <= 0 in a draw.
excess_ratio <- function(pair, map_error) {
  w <- mean(pair) / 2
  (mean(map_error) - w) / w
}

# n* from R, infinite when the map's net error is <= 0: no number of bioassays
# matches it
invert_ratio <- function(ratio) {
  ifelse(ratio > 0, 1 / ratio, Inf)
}

# n* point estimate and interval from a ratio and its bootstrap replicates
n_star_interval <- function(ratio, replicates) {
  ratio_lower <- unname(quantile(replicates, 0.025, na.rm = TRUE))
  ratio_upper <- unname(quantile(replicates, 0.975, na.rm = TRUE))
  c(n_star = invert_ratio(ratio),
    n_star_lower = invert_ratio(ratio_upper),
    n_star_upper = invert_ratio(ratio_lower))
}

