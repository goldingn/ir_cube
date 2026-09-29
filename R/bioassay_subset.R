# the subset of the collated bioassay data used in modelling, factored out of
# R/fit_model.R so that the model and the data summaries use the same filter.
# Needs R/packages.R and R/functions.R (for sample_mode()).

# the insecticide types modelled: those with at least 1000 unique places/times,
# and alpha-cypermethrin (914 unique) because of its use in LLINs
modelled_insecticides <- c("Alpha-cypermethrin",
                           "Deltamethrin",
                           "Lambda-cyhalothrin",
                           "Permethrin",
                           "Fenitrothion",
                           "Malathion",
                           "Pirimiphos-methyl",
                           "DDT",
                           "Bendiocarb")

# numbers of unique location/time records per insecticide, and the mean number
# of mosquitoes tested per location/time
count_location_years <- function(ir_africa) {
  ir_africa %>%
    group_by(insecticide_type, latitude, longitude, year_start) %>%
    summarise(
      mosquito_number = sum(mosquito_number),
      .groups = "drop"
    ) %>%
    group_by(insecticide_type) %>%
    summarise(
      n = n(),
      mosquito_number = mean(mosquito_number),
      .groups = "drop"
    ) %>%
    arrange(desc(n))
}

# subset to the modelled insecticides, the modal concentration of each, the
# modelled years, and cells inside the mask. Adds year_id (1-indexed from
# baseline_year) and the mask cell of each record.
subset_modelled_bioassays <- function(ir_africa,
                                      mask,
                                      insecticides_keep = modelled_insecticides,
                                      baseline_year = 1995,
                                      final_data_year = 2024) {
  ir_africa %>%
    filter(insecticide_type %in% insecticides_keep) %>%
    group_by(insecticide_type) %>%
    # subset to the most common concentration for each insecticide
    filter(
      concentration == sample_mode(concentration)
    ) %>%
    ungroup() %>%
    filter(
      # drop any from before the baseline
      year_start >= baseline_year,
      year_start <= final_data_year
    ) %>%
    mutate(
      # create an index to the simulation year (in 1-indexed integers)
      year_id = year_start - baseline_year + 1,
      # add on cell ids corresponding to these observations,
      cell = terra::cellFromXY(mask,
                               as.matrix(select(., longitude, latitude)))
    ) %>%
    # drop a handful of datapoints missing covariates
    filter(
      !is.na(terra::extract(mask, cell)[, 1])
    )
}
