# functions for preparing and checking the kdr (Vgsc) allele frequency data

na_to_zero <- function(x) {
  coalesce(x, 0)
}

as_number <- function(x) {
  suppressWarnings(as.numeric(x))
}

# append label to the semicolon-separated flags where condition is TRUE
add_flag <- function(flag, condition, label) {
  condition <- coalesce(condition, FALSE)
  case_when(
    !condition ~ flag,
    is.na(flag) ~ label,
    .default = paste(flag, label, sep = "; ")
  )
}

# name Vgsc mutations in An. gambiae numbering
kdr_variant_name <- function(mechanism) {
  case_when(
    grepl("L1014F|L995F", mechanism) ~ "995F",
    grepl("L1014S|L995S", mechanism) ~ "995S",
    grepl("L1014C|L995C", mechanism) ~ "995C",
    grepl("N1575Y|N1570Y", mechanism) ~ "1570Y",
    grepl("V402L", mechanism) ~ "402L",
    grepl("A1746S", mechanism) ~ "1746S",
    grepl("unspecified", mechanism, ignore.case = TRUE) ~ "unspecified"
  )
}

# assign species labels to groups in the An. gambiae complex, NA for species
# outside it. The S and M molecular forms are An. gambiae s.s. and An.
# coluzzii (Coetzee et al. 2013); samples of both forms go to the complex.
kdr_species_group <- function(species) {
  species <- tolower(species)
  case_when(
    grepl("coluzzii/gambiae|gambiae/an. coluzzii|s form_m form|m and s forms",
          species) ~ "gambiae_complex",
    grepl("coluzzii|m form", species) ~ "coluzzii",
    grepl("arabiensis", species) ~ "arabiensis",
    grepl("gambiae complex|gambiae sl|gambiae s.l.",
          species) ~ "gambiae_complex",
    grepl("gambiae", species) ~ "gambiae_ss",
    grepl("melas|merus|quadriannulatus|amharicus|bwambae|fontenillei",
          species) ~ "other_complex"
  )
}

# a tidy species label for each record
kdr_species_label <- function(species, species_group) {
  epithet <- species %>%
    tolower() %>%
    str_remove("^anopheles |^an\\. ") %>%
    str_extract("^[a-z]+")
  case_when(
    species_group == "gambiae_ss" ~ "Anopheles gambiae s.s.",
    species_group == "coluzzii" ~ "Anopheles coluzzii",
    species_group == "arabiensis" ~ "Anopheles arabiensis",
    species_group == "gambiae_complex" &
      grepl("coluzzii|form", tolower(species)) ~ "Anopheles gambiae/coluzzii",
    species_group == "gambiae_complex" ~ "Anopheles gambiae s.l.",
    species_group == "other_complex" ~ paste("Anopheles", epithet)
  )
}

# species groups that could be the same population: the same group (and the
# same species outside the main three), or either at complex level
compatible_species <- function(group_a, group_b, species_a, species_b) {
  group_a == "gambiae_complex" |
    group_b == "gambiae_complex" |
    (group_a == group_b & (group_a != "other_complex" | species_a == species_b))
}

# country names as in country_region_lookup()
kdr_country_name <- function(country) {
  country <- str_to_title(str_squish(country))
  case_when(
    grepl("^C[oô]te D", country) ~ "Côte d’Ivoire",
    country == "United Republic Of Tanzania" ~ "Tanzania",
    country %in% c("Democratic Republic Of The Congo", "Dr Congo") ~ "DR Congo",
    country == "Central African Republic" ~ "CAR",
    country %in% c("Sao Tome And Principe",
                   "Sao Tome & Principe") ~ "Sao Tome & Principe",
    country == "Guinea Bissau" ~ "Guinea-Bissau",
    .default = country
  )
}

# regions for the summaries: the UN geoscheme's, with the Horn of Africa (and
# Sudan) split from Eastern Africa, and the southern countries of Eastern
# Africa moved to Southern Africa
kdr_region <- function(country) {
  unsd <- country_region_lookup()
  region <- unsd$region[match(country, unsd$country_name)]
  case_when(
    country %in% c("Djibouti", "Eritrea", "Ethiopia", "Somalia",
                   "Sudan") ~ "Horn",
    country %in% c("Comoros", "Madagascar", "Malawi", "Mauritius",
                   "Mozambique", "Zambia", "Zimbabwe") ~ "Southern",
    region == "Western Africa" ~ "West",
    region == "Middle Africa" ~ "Central",
    region == "Eastern Africa" ~ "East",
    region == "Southern Africa" ~ "Southern",
    region == "Northern Africa" ~ "North"
  )
}

# a short key for a citation (first author's surname and year), to compare
# citations between sources
citation_key <- function(author, year) {
  surname <- author %>%
    stringi::stri_trans_general("Latin-ASCII") %>%
    str_extract("[A-Za-z'-]+") %>%
    tolower()
  paste(surname, year)
}

# great circle distance in km between points given in degrees
haversine_km <- function(lon_a, lat_a, lon_b, lat_b) {
  to_radians <- pi / 180
  d_lat <- (lat_b - lat_a) * to_radians
  d_lon <- (lon_b - lon_a) * to_radians
  h <- sin(d_lat / 2) ^ 2 +
    cos(lat_a * to_radians) * cos(lat_b * to_radians) * sin(d_lon / 2) ^ 2
  2 * 6371 * asin(pmin(1, sqrt(h)))
}

# Pick which of a set of records matching in pairs (rows `a` and `b` of each
# matching pair) to keep: take the most preferred record (the smallest
# `preference`, ties by row order), drop the records it matches, and repeat
# with the records left. Returns whether each record is kept, and the row that
# replaces each dropped record.
#
# Records are only dropped against a record that is kept, so a chain of
# matching records thins to its ends rather than to one record.
greedy_duplicates <- function(preference, a, b) {
  neighbours <- split(c(b, a), c(a, b))
  kept <- rep(NA, length(preference))
  replacement <- rep(NA_integer_, length(preference))
  for (row in order(preference, seq_along(preference))) {
    if (!is.na(kept[row])) {
      next
    }
    kept[row] <- TRUE
    matched <- neighbours[[as.character(row)]]
    matched <- matched[is.na(kept[matched])]
    kept[matched] <- FALSE
    replacement[matched] <- row
  }
  tibble(kept = kept, replaced_by = replacement)
}

# Resolve whether reported frequencies are proportions or percentages, for one
# citation at a time. Where the source also reports genotype counts, the
# frequency they imply (`reference`) picks the scale. Otherwise most of the
# citation's values decide it: values above 1 can only be percentages, and
# values between 0 and 1 are only percentages if the citation reports larger
# ones too (the percentage reading of a citation whose values are all below 1
# would put every allele frequency in the study under 1%). What is left -
# citations whose only nonzero value is exactly 1 - takes the scale the source
# usually uses. Values above 1 in a citation read as proportions are divided
# by 100, the only reading that is a frequency; these and the unresolved
# citations are flagged.
resolve_frequency_scale <- function(value, reference, source_scale) {
  nonzero <- value[!is.na(value) & value > 0]
  pairs <- !is.na(value) & !is.na(reference)
  if (any(pairs)) {
    percent_error <- median(abs(value[pairs] / 100 - reference[pairs]))
    proportion_error <- median(abs(value[pairs] - reference[pairs]))
    scale <- if_else(percent_error <= proportion_error, "percent", "proportion")
    resolved <- TRUE
  } else if (length(nonzero) == 0 || all(nonzero == 1)) {
    scale <- source_scale
    resolved <- length(nonzero) == 0
  } else {
    scale <- if_else(mean(nonzero > 1) > 0.5, "percent", "proportion")
    resolved <- TRUE
  }
  inconsistent <- !is.na(value) & scale == "proportion" & value > 1
  tibble(
    frequency = if_else(scale == "percent" | inconsistent,
                        value / 100,
                        value),
    scale_ambiguous = inconsistent | (!resolved & value > 0)
  )
}

# apply resolve_frequency_scale() to each citation of a table of records
scale_frequencies <- function(records, source_scale) {
  records %>%
    group_by(citation_key) %>%
    mutate(
      resolve_frequency_scale(reported_frequency,
                              genotype_frequency,
                              source_scale)
    ) %>%
    ungroup() %>%
    mutate(
      # where the source reports genotype counts but no frequency, use the
      # counts
      frequency = coalesce(frequency, genotype_frequency),
      flag = add_flag(flag, scale_ambiguous, "scale ambiguous")
    ) %>%
    select(-scale_ambiguous, -reported_frequency, -genotype_frequency)
}

# the periods the kdr figures and summaries use
kdr_period <- function(year) {
  cut(year,
      breaks = c(-Inf, 2009, 2014, 2019, Inf),
      labels = c("to 2009", "2010-2014", "2015-2019", "2020+"))
}

# the first word of a citation (the first author's surname, or the first word
# of a compound surname), without accents, hyphens or apostrophes
normalise_surname <- function(citation) {
  citation %>%
    stringi::stri_trans_general("Latin-ASCII") %>%
    str_squish() %>%
    str_extract("^[^ ,(]+") %>%
    tolower() %>%
    str_remove_all("[^a-z]")
}

# a DOI in a citation or URL
extract_doi <- function(text) {
  text %>%
    str_extract("10\\.[0-9]{4,}/[^ ?#&]+") %>%
    str_remove("[.;,)]+$") %>%
    tolower()
}

# The first author and year of each study with Vgsc data in the raw files, and
# whether any of its records are in the kdr records (prep_kdr.R drops records
# without coordinates or a year, located only to a wide area, or of survivors
# or dead mosquitoes only). The citations are those prep_kdr.R gives: a DOI
# for Vector Atlas (or the author and year where there is none), "Author et
# al. year" for the IR Mapper contact files, and a full reference for the
# analysis-ready data, so we take the authors and years from the raw files.
# `kdr_records` are the records read from data/clean/kdr_records.csv.
kdr_study_lookup <- function(kdr_records) {
  va_studies <- read_csv("data/raw/VA_FULL_DATA_20250721.csv",
                         col_types = cols(.default = "c"),
                         col_select = c(insecticide_resistance_data,
                                        citation_doi, author, publication_year,
                                        starts_with("vgsc"))) %>%
    filter(
      insecticide_resistance_data %in% c("genotypic", "both"),
      if_any(starts_with("vgsc") & !matches("method|generation|notes"),
             ~ !is.na(.x))
    ) %>%
    transmute(
      source = "vector_atlas",
      citation = coalesce(citation_doi, paste(author, publication_year)),
      surname = normalise_surname(author),
      publication_year = as_number(publication_year),
      doi = extract_doi(citation_doi)
    )

  ir_mapper_studies <- bind_rows(
    read_xlsx("data/raw/IR Mapper data August 2022.xlsx",
              col_types = "text") %>%
      filter(grepl("^Kdr", IR_Mechanism_Name, ignore.case = TRUE)) %>%
      transmute(source = "ir_mapper_2022", citation = Reference_Name),
    read_xlsx("data/raw/IR_Mapper_Anopheles_data_2021_Sept_2024.xlsx",
              col_types = "text") %>%
      filter(grepl("^Kdr", iR_Mechanism_Name, ignore.case = TRUE)) %>%
      transmute(source = "ir_mapper_2024", citation = reference_Name)
  ) %>%
    mutate(
      surname = normalise_surname(citation),
      publication_year = as_number(str_extract(citation, "(19|20)[0-9]{2}")),
      doi = NA_character_
    )

  ready_studies <- c("6_Vgsc-allele-freq_complex-subgroup.csv",
                     "7_Vgsc-allele-freq_species.csv") %>%
    map(\(file) read.csv(file.path("data/raw",
                                   "IRMapper Analysis ready data sets",
                                   file),
                         fileEncoding = "latin1",
                         check.names = FALSE)) %>%
    map(\(x) select(x, citation = `Source 1 citation`,
                    publication_year = `Publication year`)) %>%
    bind_rows() %>%
    transmute(
      source = "ir_mapper_ready",
      citation,
      surname = normalise_surname(citation),
      publication_year = as_number(publication_year),
      doi = NA_character_
    )

  bind_rows(va_studies, ir_mapper_studies, ready_studies) %>%
    distinct() %>%
    mutate(
      in_records = paste(source, citation) %in%
        paste(kdr_records$source, kdr_records$citation)
    )
}
