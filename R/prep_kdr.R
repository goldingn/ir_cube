# prep kdr (Vgsc) allele frequency data, for the latent kdr type surface

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# Combines Vgsc allele frequencies in the An. gambiae complex from Vector
# Atlas, the two IR Mapper contact files (August 2022, and 2021 to September
# 2024) and the IR Mapper analysis-ready data (Moyes et al. 2019, Sci Data
# 6:121). One row per record (a sample and a variant). Variants are in An.
# gambiae numbering: L1014F, L1014S and L1014C (Musca domestica numbering) are
# 995F, 995S and 995C, and N1575Y is 1570Y. Frequencies are proportions of
# alleles; n_tested is the number of mosquitoes genotyped.

# records outside this box (longitude, latitude) have bad coordinates
africa_bbox <- c(xmin = -26, xmax = 64, ymin = -36, ymax = 38)

# records from different sources of the same variant and start year, at sites
# within duplicate_km of each other, with compatible species groups and
# frequencies within duplicate_tolerance, are duplicates. Such matches within a
# source, and pairs up to near_miss_km apart, are listed for checking.
duplicate_km <- 5
near_miss_km <- 10
duplicate_tolerance <- 0.02

# the sources, in order of preference for the copy of a duplicate that is kept,
# and their short names for record ids
kdr_sources <- c(vector_atlas = "va",
                 ir_mapper_2024 = "irm24",
                 ir_mapper_2022 = "irm22",
                 ir_mapper_ready = "irmr")

# helper functions ------------------------------------------------------------

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
    country %in% c("Cote D Ivoire", "Cote Divoire") ~ "Côte d’Ivoire",
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

# Vector Atlas ----------------------------------------------------------------

# the Vgsc columns: the genotype counts of each variant (homozygote,
# heterozygotes) and the reported allele percentage
va_vgsc_columns <- tribble(
  ~variant, ~homozygote, ~heterozygotes, ~percent,
  "995F", "vgsc995f_vgsc995f_n",
  c("vgsc995l_vgsc995f_n", "vgsc995f_vgsc995s_n", "vgsc995f_vgsc995c_n"),
  "vgsc995f_percent",
  "995S", "vgsc995s_vgsc995s_n",
  c("vgsc995l_vgsc995s_n", "vgsc995f_vgsc995s_n"),
  "vgsc995s_percent",
  "995C", "vgsc995c_vgsc995c_n",
  c("vgsc995l_vgsc995c_n", "vgsc995f_vgsc995c_n"),
  "vgsc995c_percent",
  "1570Y", "vgsc1570y_vgsc1570y_n", "vgsc1570n_vgsc1570y_n",
  "vgsc1570y_percent"
)

# all the genotype count columns of the 995 and 1570 codons, whose sum is the
# number of mosquitoes the genotype counts describe
va_995_counts <- c("vgsc995l_vgsc995l_n", "vgsc995l_vgsc995f_n",
                   "vgsc995f_vgsc995f_n", "vgsc995l_vgsc995s_n",
                   "vgsc995s_vgsc995s_n", "vgsc995f_vgsc995s_n",
                   "vgsc995l_vgsc995c_n", "vgsc995c_vgsc995c_n",
                   "vgsc995f_vgsc995c_n")
va_1570_counts <- c("vgsc1570n_vgsc1570n_n", "vgsc1570n_vgsc1570y_n",
                    "vgsc1570y_vgsc1570y_n")

# the genotypic records of the Vector Atlas extract, one row per sample
va_raw <- read_csv(
  "data/raw/VA_FULL_DATA_20250721.csv",
  col_types = cols(.default = "c"),
  col_select = c(insecticide_resistance_data, citation_doi, author,
                 publication_year, country, site, latitude_1, longitude_1,
                 area_type, year_start, year_end, species, species_notes,
                 bioassay_subsample_used_in_genotypic_test,
                 vgsc_number_of_mosquitoes_tested,
                 all_of(c(va_995_counts, va_1570_counts)),
                 all_of(va_vgsc_columns$percent))) %>%
  filter(
    insecticide_resistance_data %in% c("genotypic", "both")
  ) %>%
  mutate(
    across(c(latitude_1, longitude_1, year_start, year_end,
             vgsc_number_of_mosquitoes_tested,
             all_of(c(va_995_counts, va_1570_counts)),
             all_of(va_vgsc_columns$percent)),
           as_number),
    sample_id = sprintf("va_%04i", row_number())
  ) %>%
  # drop samples of bioassay survivors or dead mosquitoes only: these are
  # selected, not population frequencies
  filter(
    !(bioassay_subsample_used_in_genotypic_test %in%
        c("alive only", "dead only", "Alive", "Dead")),
    !(species_notes %in% c("Alive", "Dead"))
  ) %>%
  # drop samples whose location is only known to a wide area
  filter(
    !(area_type %in% c("large polygon", "wide area", "admin unit centroid"))
  )

# the records of one variant in the Vector Atlas extract, with the allele
# frequency its genotype counts imply (over all the mosquitoes genotyped at
# that codon, so heterozygotes of two mutations count in both their variants)
va_variant_records <- function(variant, homozygote, heterozygotes, percent) {
  counts <- if (variant == "1570Y") va_1570_counts else va_995_counts
  va_raw %>%
    mutate(
      variant = variant,
      n_genotyped = rowSums(across(all_of(counts), na_to_zero)),
      alleles = 2 * na_to_zero(.data[[homozygote]]) +
        rowSums(across(all_of(heterozygotes), na_to_zero)),
      typed = if_any(all_of(c(homozygote, heterozygotes)), ~ !is.na(.x)),
      genotype_frequency = if_else(typed & n_genotyped > 0,
                                   alleles / (2 * n_genotyped),
                                   NA),
      reported_frequency = .data[[percent]]
    )
}

va_records <- va_vgsc_columns %>%
  pmap(va_variant_records) %>%
  bind_rows() %>%
  filter(
    !is.na(reported_frequency) | !is.na(genotype_frequency)
  ) %>%
  mutate(
    source = "vector_atlas",
    citation = coalesce(citation_doi, paste(author, publication_year)),
    citation_key = citation_key(author, publication_year),
    country = kdr_country_name(country),
    species_group = kdr_species_group(species),
    n_tested = case_when(
      !is.na(vgsc_number_of_mosquitoes_tested) &
        vgsc_number_of_mosquitoes_tested > 0 ~
        vgsc_number_of_mosquitoes_tested,
      n_genotyped > 0 ~ n_genotyped
    ),
    flag = NA_character_,
    flag = add_flag(flag,
                    area_type %in% c("multi-point", "multiple points"),
                    "multiple sites"),
    flag = add_flag(flag, area_type == "small polygon", "coarse location")
  ) %>%
  select(source, citation, citation_key, country, latitude = latitude_1,
         longitude = longitude_1, year_start, year_end, species, species_group,
         variant, reported_frequency, genotype_frequency, n_tested, sample_id,
         flag) %>%
  scale_frequencies(source_scale = "percent")

# IR Mapper -------------------------------------------------------------------

# Some Kdr_frequency entries of the August 2022 file are text: either two
# frequencies, one per molecular form, or one frequency followed by the
# frequency of the double mutant. Split the first into a record per form, take
# the leading number of the second, and leave the rest (empty strings, and two
# entries whose forms cannot be told apart) unparsed.
expand_form_frequencies <- function(records, frequency, forms) {
  records <- records %>%
    mutate(
      text = {{ frequency }},
      m_value = as_number(str_match(text, "M\\s*([0-9]*\\.?[0-9]+)")[, 2]),
      s_value = as_number(str_match(text, "S\\s*([0-9]*\\.?[0-9]+)")[, 2]),
      split_forms = is.na(as_number(text)) & (!is.na(m_value) | !is.na(s_value))
    )
  numeric_records <- records %>%
    filter(!split_forms) %>%
    mutate(
      reported_frequency = coalesce(
        as_number(text),
        as_number(str_match(text, "^\\s*([0-9]*\\.?[0-9]+)")[, 2])
      ),
      flag = add_flag(flag,
                      is.na(as_number(text)) & !is.na(reported_frequency),
                      "frequency text trimmed")
    )
  # one record per form, with that form's species
  form_records <- records %>%
    filter(split_forms) %>%
    pivot_longer(c(m_value, s_value),
                 names_to = "form",
                 values_to = "reported_frequency") %>%
    filter(!is.na(reported_frequency)) %>%
    mutate(
      {{ forms }} := if_else(form == "m_value", "M", "S"),
      flag = add_flag(flag, TRUE, "frequency split by molecular form")
    ) %>%
    select(-form)
  bind_rows(numeric_records, form_records) %>%
    select(-text, -split_forms, -any_of(c("m_value", "s_value")))
}

# S and M are the molecular forms of An. gambiae s.s. and An. coluzzii
# (Coetzee et al. 2013), so where IR Mapper records the forms of a sample
# identified to An. gambiae s.s. or to the complex, they identify the species.
# A sample of both forms stays at complex level.
molecular_form_group <- function(species_group, forms) {
  has_m <- grepl("\\bM\\s*[0-9&(]|\\bM\\s*$|\\bM\\)", paste0(forms, " "))
  has_s <- grepl("\\bS\\s*[0-9&(]|\\bS\\s*$|\\bS\\)", paste0(forms, " "))
  case_when(
    !species_group %in% c("gambiae_ss", "gambiae_complex") ~ species_group,
    has_m & has_s ~ "gambiae_complex",
    has_m ~ "coluzzii",
    has_s ~ "gambiae_ss",
    .default = species_group
  )
}

# The kdr rows of the August 2022 contact file. These are PCR tests, with no
# mortality or resistance status, so IR_Test_NumExposed is the number
# genotyped.
ir_mapper_2022 <- read_xlsx("data/raw/IR Mapper data August 2022.xlsx",
                            col_types = "text") %>%
  filter(
    grepl("^Kdr", IR_Mechanism_Name, ignore.case = TRUE)
  ) %>%
  mutate(
    across(c(Latitude, Longitude, IR_Test_NumExposed,
             Mosquito_collection_Year_start, Mosquito_collection_year_end),
           as_number),
    flag = NA_character_
  ) %>%
  expand_form_frequencies(Kdr_frequency, Molecular_forms) %>%
  mutate(
    source = "ir_mapper_2022",
    variant = kdr_variant_name(IR_Mechanism_Name),
    country = kdr_country_name(Country),
    species_group = kdr_species_group(Vector_species),
    species_group = molecular_form_group(species_group, Molecular_forms),
    genotype_frequency = NA_real_,
    citation_key = citation_key(Reference_Name,
                               str_extract(Reference_Name, "[0-9]{4}")),
    # each row is one variant of a sample, so group the variants typed on the
    # same sample
    sample_id = sprintf("irm22_%04i",
                        as.integer(factor(paste(citation_key, Locality,
                                                Latitude, Longitude,
                                                Mosquito_collection_Year_start,
                                                Vector_species,
                                                Molecular_forms,
                                                IR_Test_NumExposed))))
  ) %>%
  select(source, citation = Reference_Name, citation_key, country,
         latitude = Latitude, longitude = Longitude,
         year_start = Mosquito_collection_Year_start,
         year_end = Mosquito_collection_year_end, species = Vector_species,
         species_group, variant, reported_frequency, genotype_frequency,
         n_tested = IR_Test_NumExposed, sample_id, flag)

# the rows with no usable frequency, for checking
dir.create("outputs/review", showWarnings = FALSE, recursive = TRUE)
ir_mapper_2022 %>%
  filter(is.na(reported_frequency)) %>%
  write.csv("outputs/review/kdr_unparsed_frequencies.csv", row.names = FALSE)

ir_mapper_2022 <- ir_mapper_2022 %>%
  filter(!is.na(reported_frequency)) %>%
  scale_frequencies(source_scale = "proportion")

# the kdr rows of the 2021 to September 2024 contact file
ir_mapper_2024 <- read_xlsx(
  "data/raw/IR_Mapper_Anopheles_data_2021_Sept_2024.xlsx",
  col_types = "text") %>%
  filter(
    grepl("^Kdr", iR_Mechanism_Name, ignore.case = TRUE)
  ) %>%
  mutate(
    across(c(latitude, longitude, `iR_Test_NumExposed/Tested`,
             collection_Year_Start, collection_Year_End,
             resistance_Mechanism_Frequency),
           as_number),
    source = "ir_mapper_2024",
    variant = kdr_variant_name(iR_Mechanism_Name),
    country = kdr_country_name(country),
    species_group = kdr_species_group(vector_Species),
    reported_frequency = resistance_Mechanism_Frequency,
    genotype_frequency = NA_real_,
    flag = NA_character_,
    citation_key = citation_key(reference_Name,
                               str_extract(reference_Name, "[0-9]{4}")),
    sample_id = sprintf("irm24_%04i",
                        as.integer(factor(paste(citation_key, locality,
                                                latitude, longitude,
                                                collection_Year_Start,
                                                vector_Species,
                                                `iR_Test_NumExposed/Tested`))))
  ) %>%
  select(source, citation = reference_Name, citation_key, country,
         latitude, longitude, year_start = collection_Year_Start,
         year_end = collection_Year_End, species = vector_Species,
         species_group, variant, reported_frequency, genotype_frequency,
         n_tested = `iR_Test_NumExposed/Tested`, sample_id, flag) %>%
  scale_frequencies(source_scale = "proportion")

# IR Mapper analysis-ready data sets ------------------------------------------

# the genotype count columns of the 995 codon, whose sum is the number of
# mosquitoes the counts describe
ready_counts <- c("L/L (no.)", "L/F (no.)", "F/F (no.)", "L/S (no.)",
                  "S/S (no.)", "F/S (no.)")

# Read one of the analysis-ready Vgsc allele frequency files (Moyes et al.
# 2019, Sci Data 6:121), with a record per variant. These report percentages
# of the three alleles of the 995 codon, which sum to 100, and genotype counts
# for a third of the rows. Up to 20 sites are listed per row; we take the
# first, as the coordinates of the sample, and flag the rest.
read_ready_vgsc <- function(file, species_column) {
  records <- read.csv(file.path("data/raw/IRMapper Analysis ready data sets",
                                file),
                      fileEncoding = "latin1",
                      check.names = FALSE,
                      na.strings = c("", "NA", "NR", "NF"))
  names(records) <- trimws(names(records))
  records %>%
    as_tibble() %>%
    mutate(
      across(c(Latitude, Longitude, `Start year`, `End year`,
               `No. mosquitoes tested`, all_of(ready_counts),
               `L1014F (%)`, `L1014S (%)`),
             as_number),
      sites = rowSums(!is.na(across(matches("^Latitu(de|te) [0-9]+$")))),
      n_genotyped = rowSums(across(all_of(ready_counts), na_to_zero)),
      typed = if_any(all_of(ready_counts), ~ !is.na(.x)),
      f_alleles = 2 * na_to_zero(`F/F (no.)`) + na_to_zero(`L/F (no.)`) +
        na_to_zero(`F/S (no.)`),
      s_alleles = 2 * na_to_zero(`S/S (no.)`) + na_to_zero(`L/S (no.)`) +
        na_to_zero(`F/S (no.)`),
      sample_id = sprintf("irmr_%s", `kdr test ID`),
      species = .data[[species_column]],
      flag = NA_character_,
      flag = add_flag(flag, sites > 0, "multiple sites"),
      flag = add_flag(flag, `Site type` == "polygon", "coarse location")
    ) %>%
    pivot_longer(c(`L1014F (%)`, `L1014S (%)`),
                 names_to = "variant",
                 values_to = "reported_frequency") %>%
    mutate(
      variant = kdr_variant_name(variant),
      genotype_frequency = case_when(
        !typed | n_genotyped == 0 ~ NA,
        variant == "995F" ~ f_alleles / (2 * n_genotyped),
        variant == "995S" ~ s_alleles / (2 * n_genotyped)
      ),
      test_id = `kdr test ID`,
      citation = `Source 1 citation`,
      citation_key = citation_key(citation, `Publication year`),
      species_group = kdr_species_group(species),
      n_tested = coalesce(`No. mosquitoes tested`,
                          if_else(n_genotyped > 0, n_genotyped, NA))
    ) %>%
    filter(
      !is.na(reported_frequency) | !is.na(genotype_frequency)
    ) %>%
    select(test_id, citation, citation_key, latitude = Latitude,
           longitude = Longitude, year_start = `Start year`,
           year_end = `End year`, species, species_group, variant,
           reported_frequency, genotype_frequency, n_tested, sample_id, flag)
}

ready_species <- read_ready_vgsc("7_Vgsc-allele-freq_species.csv", "Species")
ready_complex <- read_ready_vgsc("6_Vgsc-allele-freq_complex-subgroup.csv",
                                 "Complex/Subgroup")

# The same kdr test appears in both files where the sample was identified to
# species: keep the species-level rows of those tests, as prep_bioassays.R does
# for the bioassays.
ir_mapper_ready <- ready_species %>%
  bind_rows(
    ready_complex %>%
      filter(!test_id %in% ready_species$test_id)
  ) %>%
  mutate(
    source = "ir_mapper_ready"
  ) %>%
  select(-test_id) %>%
  scale_frequencies(source_scale = "percent")

# combine and clean -----------------------------------------------------------

# the African country each record is in, from its coordinates (the
# analysis-ready data sets give no country), taking the nearest country for
# points just offshore
country_from_coordinates <- function(longitude, latitude) {
  borders <- readRDS("data/clean/country_borders.RDS")
  points <- st_as_sf(tibble(longitude, latitude),
                     coords = c("longitude", "latitude"),
                     crs = st_crs(borders))
  borders$country_name[as.integer(st_nearest_feature(points, borders))]
}

kdr_all <- bind_rows(va_records,
                     ir_mapper_2024,
                     ir_mapper_2022,
                     ir_mapper_ready) %>%
  # drop records without the basics
  filter(
    !is.na(frequency),
    !is.na(latitude),
    !is.na(longitude),
    !is.na(year_start),
    !is.na(variant)
  ) %>%
  # drop the bad coordinates (one Vector Atlas record has longitude -42833)
  filter(
    between(longitude, africa_bbox["xmin"], africa_bbox["xmax"]),
    between(latitude, africa_bbox["ymin"], africa_bbox["ymax"])
  ) %>%
  # drop species outside the An. gambiae complex
  filter(
    !is.na(species_group)
  ) %>%
  mutate(
    source = factor(source, levels = names(kdr_sources)),
    species = kdr_species_label(species, species_group),
    country = coalesce(country, country_from_coordinates(longitude, latitude)),
    region = kdr_region(country),
    year_end = coalesce(year_end, year_start),
    # a few records report no mosquitoes genotyped, which we treat as missing
    n_tested = if_else(n_tested == 0, NA, n_tested),
    flag = add_flag(flag, is.na(n_tested), "no sample size")
  )

# drop any frequency above 1, which is impossible on either scale
cat(sprintf("records dropped with a frequency above 1: %i\n",
            sum(kdr_all$frequency > 1)))

kdr_all <- kdr_all %>%
  filter(
    frequency <= 1
  ) %>%
  arrange(source) %>%
  mutate(
    record_id = sprintf("%s_%04i", kdr_sources[as.character(source)],
                        row_number())
  )

# deduplicate -----------------------------------------------------------------

# All pairs of records of the same variant and start year within near_miss_km
# of each other, with the distance between them and the difference in their
# frequencies.
match_kdr_pairs <- function(records, max_km) {
  keys <- records %>%
    mutate(row = row_number()) %>%
    select(row, variant, year_start, longitude, latitude, frequency)
  keys %>%
    inner_join(keys,
               by = c("variant", "year_start"),
               suffix = c("_a", "_b"),
               relationship = "many-to-many") %>%
    filter(row_a < row_b) %>%
    mutate(
      km = haversine_km(longitude_a, latitude_a, longitude_b, latitude_b),
      frequency_difference = abs(frequency_a - frequency_b)
    ) %>%
    filter(km <= max_km) %>%
    select(a = row_a, b = row_b, km, frequency_difference)
}

pairs <- match_kdr_pairs(kdr_all, near_miss_km) %>%
  mutate(
    compatible = compatible_species(kdr_all$species_group[a],
                                    kdr_all$species_group[b],
                                    kdr_all$species[a],
                                    kdr_all$species[b]),
    same_source = kdr_all$source[a] == kdr_all$source[b],
    match = compatible &
      km <= duplicate_km &
      frequency_difference <= duplicate_tolerance,
    # Only matches between sources are duplicates. Matches within a source,
    # mostly within a citation, are most likely distinct samples (neighbouring
    # villages, seasons or species, many at frequency 0 or 1), so as for the
    # bioassays (#31) they are listed for checking and kept.
    duplicate = match & !same_source,
    near_miss = compatible & !match &
      (km <= duplicate_km | frequency_difference <= duplicate_tolerance)
  )

# Of each pair of duplicates, keep the record with a sample size, then the one
# identified to species, then the one from the preferred source.
duplicates <- pairs %>%
  filter(duplicate)
kdr_all <- kdr_all %>%
  mutate(
    preference = 100 * as.integer(is.na(n_tested)) +
      10 * as.integer(species_group == "gambiae_complex") +
      as.integer(source),
    greedy_duplicates(preference, duplicates$a, duplicates$b),
    kept_record_id = record_id[replaced_by]
  )

# the matches and near misses, for checking
review_columns <- c("record_id", "source", "citation", "latitude", "longitude",
                    "year_start", "species", "variant", "frequency",
                    "n_tested", "kept")
kdr_review <- bind_cols(
  kdr_all[pairs$a, review_columns] %>%
    rename_with(~ paste0(.x, "_a")),
  kdr_all[pairs$b, review_columns] %>%
    rename_with(~ paste0(.x, "_b"))
) %>%
  mutate(
    km = pairs$km,
    frequency_difference = pairs$frequency_difference,
    match = case_when(
      pairs$duplicate ~ "duplicate",
      pairs$match ~ "within_source",
      pairs$near_miss & pairs$km > duplicate_km ~ "5-10 km apart",
      pairs$near_miss ~ "frequency differs",
      .default = "incompatible species"
    )
  ) %>%
  filter(
    match != "incompatible species"
  ) %>%
  arrange(match, desc(km))

kdr_review %>%
  write.csv("outputs/review/kdr_near_matches.csv", row.names = FALSE)

# the deduplicated records
kdr_records <- kdr_all %>%
  filter(kept) %>%
  select(record_id, source, citation, country, latitude, longitude, year_start,
         year_end, species, species_group, variant, frequency, n_tested,
         sample_id, flag)

write.csv(kdr_records, "data/clean/kdr_records.csv", row.names = FALSE)

# map -------------------------------------------------------------------------

# the periods the figure facets by
kdr_period <- function(year) {
  cut(year,
      breaks = c(-Inf, 2009, 2014, 2019, Inf),
      labels = c("to 2009", "2010-2014", "2015-2019", "2020+"))
}

borders <- readRDS("data/clean/country_borders.RDS")

map_records <- kdr_records %>%
  filter(
    variant %in% c("995F", "995S")
  ) %>%
  mutate(
    period = kdr_period(year_start)
  ) %>%
  # plot the highest frequencies last, so they are not hidden
  arrange(frequency)

kdr_map <- ggplot() +
  geom_sf(
    data = borders,
    fill = "white",
    colour = grey(0.8),
    linewidth = 0.2
  ) +
  geom_point(
    aes(
      x = longitude,
      y = latitude,
      fill = frequency
    ),
    data = map_records,
    shape = 21,
    colour = grey(0.4),
    stroke = 0.15,
    size = 1.1,
    alpha = 0.9
  ) +
  facet_grid(variant ~ period) +
  scale_fill_gradient(
    name = "allele frequency",
    low = "#deebf7",
    high = "#08306b",
    limits = c(0, 1)
  ) +
  coord_sf(xlim = c(-18, 51), ylim = c(-35, 20)) +
  labs(
    title = sprintf("kdr (Vgsc) allele frequency records, %i to %i",
                    min(map_records$year_start),
                    max(map_records$year_start))
  ) +
  theme_ir_maps() +
  theme(strip.background = element_blank())

ggsave(
  "figures/kdr_records_map.png",
  plot = kdr_map,
  bg = "white",
  width = 9,
  height = 4.5,
  dpi = 200
)

# summarise -------------------------------------------------------------------

# records by variant and source, before and after deduplication
kdr_all %>%
  count(source, variant, kept) %>%
  pivot_wider(names_from = kept, values_from = n, names_prefix = "kept_",
              values_fill = 0) %>%
  mutate(records = kept_TRUE + kept_FALSE) %>%
  select(source, variant, records, deduplicated = kept_TRUE) %>%
  arrange(variant, source) %>%
  print(n = Inf)

print(count(kdr_records, variant, period = kdr_period(year_start)) %>%
        pivot_wider(names_from = period, values_from = n, values_fill = 0))
print(count(kdr_records, region = kdr_region(country), variant) %>%
        pivot_wider(names_from = variant, values_from = n, values_fill = 0))
print(count(kdr_records, species_group, variant) %>%
        pivot_wider(names_from = variant, values_from = n, values_fill = 0))

# samples with both of the two common variants typed
cat(sprintf("samples typing both 995F and 995S: %i of %i\n",
            kdr_records %>%
              filter(variant %in% c("995F", "995S")) %>%
              distinct(sample_id, variant) %>%
              count(sample_id) %>%
              filter(n == 2) %>%
              nrow(),
            n_distinct(kdr_records$sample_id)))

cat(sprintf("records missing the number genotyped: %i of %i\n",
            sum(is.na(kdr_records$n_tested)),
            nrow(kdr_records)))

print(count(kdr_records, flag))
print(count(kdr_review, match))
