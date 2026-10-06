# check whether the Malaria Threats Map kdr records duplicate the kdr allele
# frequency records from Vector Atlas and IR Mapper (R/prep_kdr.R)

# load packages and functions
source("R/packages.R")
source("R/functions.R")
source("R/kdr_functions.R")

# The MTM kdr data (the biochemical assay download) give only whether a kdr
# mutation was detected, not its frequency. For the African An. gambiae complex
# L1014F (995F) and L1014S (995S) records, we look for:
#   1. a site match: a kdr record of the same variant and start year, within
#      match_km and of a compatible species group;
#   2. a study match: the same study (first author's surname and year, or DOI)
#      among the kdr records, anywhere;
# and summarise the MTM records matching neither.
match_km <- 5

# MTM records ------------------------------------------------------------------

mtm_kdr <- read_xlsx("data/raw/MTM_BIOCHEMICAL_ASSAY_20261006.xlsx",
                     sheet = "Data",
                     col_types = "text") %>%
  filter(
    grepl("^kdr", TEST_TYPE)
  ) %>%
  mutate(
    across(c(LATITUDE, LONGITUDE, YEAR_START, MOSQUITO_NUMBER), as_number),
    variant = kdr_variant_name(TEST_TYPE),
    country = kdr_country_name(COUNTRY_NAME),
    region = kdr_region(country),
    # hybrids of the two species are only identified to the complex
    species_group = if_else(grepl("hybrid", `VECTOR SPECIES`),
                            "gambiae_complex",
                            kdr_species_group(`VECTOR SPECIES`)),
    species = kdr_species_label(`VECTOR SPECIES`, species_group),
    surname = normalise_surname(CITATION),
    publication_year = as_number(coalesce(
      str_match(CITATION, "\\(([0-9]{4})\\)")[, 2],
      str_extract(CITATION, "(19|20)[0-9]{2}")
    )),
    doi = coalesce(extract_doi(CITATION_URL), extract_doi(CITATION)),
    # citations in the form "Author et al. (year) title" are publications, the
    # rest name a programme or institution
    published = grepl("et al|\\([0-9]{4}\\)", CITATION)
  ) %>%
  # African records of the An. gambiae complex
  filter(
    country %in% country_region_lookup()$country_name,
    !is.na(species_group)
  )

# the unspecified kdr records are only summarised; they cannot be matched to a
# variant
mtm_unspecified <- mtm_kdr %>%
  filter(variant == "unspecified")
mtm_kdr <- mtm_kdr %>%
  filter(variant %in% c("995F", "995S"))

# kdr records and their studies -----------------------------------------------

kdr_records <- read_csv("data/clean/kdr_records.csv",
                        col_types = cols(variant = "c", .default = col_guess()))

kdr_studies <- kdr_study_lookup(kdr_records)

# 1. site matches --------------------------------------------------------------

# all pairs of an MTM record and a kdr record of the same variant, with start
# years at most max_years apart, within max_km and of compatible species
site_pairs <- function(mtm, records, max_km, max_years = 0) {
  mtm %>%
    select(ID, variant, YEAR_START, LONGITUDE, LATITUDE, species_group,
           species) %>%
    inner_join(
      records %>%
        select(record_id, variant, year_start, longitude, latitude,
               species_group, species, frequency, n_tested, source, citation),
      by = "variant",
      suffix = c("_mtm", ""),
      relationship = "many-to-many"
    ) %>%
    filter(
      abs(YEAR_START - year_start) <= max_years
    ) %>%
    mutate(
      km = haversine_km(LONGITUDE, LATITUDE, longitude, latitude)
    ) %>%
    filter(
      km <= max_km,
      compatible_species(species_group_mtm, species_group, species_mtm,
                         species)
    )
}

# the nearest matching record of each MTM record, and how many there are
nearest_site <- site_pairs(mtm_kdr, kdr_records, match_km) %>%
  group_by(ID) %>%
  mutate(site_matches = n()) %>%
  slice_min(km, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(ID, site_matches, matched_record_id = record_id,
         matched_source = source, matched_citation = citation,
         matched_frequency = frequency, matched_n_tested = n_tested,
         matched_km = km)

# a looser site match, to see why records of studies in the kdr records do not
# match: start years a year apart, and sites up to 25 km apart
loose_site_ids <- site_pairs(mtm_kdr, kdr_records, 25, max_years = 1) %>%
  pull(ID) %>%
  unique()

# 2. study matches -------------------------------------------------------------

# whether each record's study (by surname and year, or DOI) is among
# `studies`, and whether it is with the year one either side
study_match <- function(records, studies) {
  keys <- studies %>%
    filter(!is.na(surname), !is.na(publication_year)) %>%
    transmute(key = paste(surname, publication_year)) %>%
    pull(key)
  dois <- na.omit(unique(studies$doi))
  in_studies <- function(year) {
    paste(records$surname, year) %in% keys
  }
  in_studies(records$publication_year) |
    (!is.na(records$doi) & records$doi %in% dois)
}

study_years_apart <- function(records, studies) {
  study_match(records, studies) |
    study_match(mutate(records, publication_year = publication_year - 1),
                studies) |
    study_match(mutate(records, publication_year = publication_year + 1),
                studies)
}

mtm_matches <- mtm_kdr %>%
  mutate(
    study_in_records = study_match(., filter(kdr_studies, in_records)),
    study_in_sources = study_match(., kdr_studies),
    study_in_sources_year_1 = study_years_apart(., kdr_studies)
  ) %>%
  left_join(nearest_site, by = "ID") %>%
  mutate(
    site_matches = coalesce(site_matches, 0L),
    loose_site_match = ID %in% loose_site_ids,
    match = case_when(
      site_matches > 0 ~ "site",
      study_in_records ~ "study in records",
      study_in_sources ~ "study dropped from records",
      .default = "none"
    ),
    period = kdr_period(YEAR_START)
  )

mtm_matches %>%
  select(ID, country, region, LATITUDE, LONGITUDE, YEAR_START, species,
         variant, MOSQUITO_NUMBER, MECHANISM_STATUS, STAGE_ORIGIN, DATA_SOURCE,
         CITATION, CITATION_URL, published, match, study_in_records,
         study_in_sources, study_in_sources_year_1, loose_site_match,
         site_matches, matched_record_id, matched_source, matched_citation,
         matched_frequency, matched_n_tested, matched_km) %>%
  write.csv("outputs/review/mtm_kdr_matches.csv", row.names = FALSE)

# summaries -------------------------------------------------------------------

cat(sprintf(paste("African gambiae complex MTM kdr records:",
                  "%i L1014F/L1014S, %i unspecified\n"),
            nrow(mtm_kdr), nrow(mtm_unspecified)))
cat(sprintf("unspecified kdr records whose study is in the kdr records: %i\n",
            sum(study_match(mtm_unspecified,
                            filter(kdr_studies, in_records)))))

print(count(mtm_matches, match, variant) %>%
        pivot_wider(names_from = variant, values_from = n))
print(count(mtm_matches, match, loose_site_match))
print(count(mtm_matches, match, study_in_sources_year_1))

# status consistency for the site matches: whether MTM's detected or not
# detected agrees with the nearest record's frequency, and whether MTM's number
# of mosquitoes is the record's number genotyped
site_matched <- mtm_matches %>%
  filter(match == "site")
print(count(site_matched,
            MECHANISM_STATUS,
            matched_frequency = if_else(matched_frequency > 0, "> 0", "0")) %>%
        pivot_wider(names_from = matched_frequency, values_from = n,
                    values_fill = 0))
print(count(site_matched,
            same_n = case_when(
              is.na(MOSQUITO_NUMBER) | is.na(matched_n_tested) ~ "missing",
              MOSQUITO_NUMBER == matched_n_tested ~ "same",
              .default = "different"
            )))

# the records matching neither
unmatched <- mtm_matches %>%
  filter(match == "none")
print(count(unmatched, period))
print(count(unmatched, region))
print(count(unmatched, published, period) %>%
        pivot_wider(names_from = period, values_from = n))
print(count(unmatched, DATA_SOURCE, sort = TRUE), n = 20)
print(count(unmatched, CITATION, CITATION_URL, sort = TRUE) %>%
        mutate(CITATION = str_trunc(CITATION, 70)), n = 40)
cat(sprintf("unmatched records missing MOSQUITO_NUMBER: %i of %i\n",
            sum(is.na(unmatched$MOSQUITO_NUMBER)),
            nrow(unmatched)))
