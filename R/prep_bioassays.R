# Prepare the bioassay data for modelling: read the five databases, put them on
# a common format, remove copies of the same bioassay, and save the An. gambiae
# complex records to data/clean/all_gambiae_complex_data.RDS. The supplement
# (doc/supp.qmd, "Bioassay data processing") describes the method;
# outputs/review/bioassay_attrition.csv counts the records each step removes,
# outputs/review/bioassay_paper_versions.csv lists each paper's versions and
# the one kept, outputs/review/bioassay_links.csv lists every other copy
# removed, outputs/review/bioassay_review_pairs.csv the pairs of copies
# flagged for review, and outputs/review/bioassay_country_relabels.csv the
# records whose country label disagrees with their coordinates.

source("R/packages.R")
source("R/functions.R")
source("R/bioassay_prep_functions.R")
source("R/bioassay_subset.R")

# ---- read the databases ----

# The 20261005 MTM discriminating concentration download fills missing values
# of MORTALITY_ADJUSTED with 0: it has no missing values, and all 306 records
# missing in the 20251031 download are exactly 0, 87 of them labelled
# susceptible. The papers checked report nonzero mortality for these (#27). So
# we drop them, by ID against the 20251031 download; the 232 zeros already in
# that download are kept.
mtm_filled_mortality_ids <- read_xlsx(
  "data/raw/MTM_DISCRIMINATING_CONCENTRATION_BIOASSAY_20251031.xlsx",
  sheet = "Data", col_types = "text") %>%
  filter(is.na(MORTALITY_ADJUSTED)) %>%
  pull(ID)

mtm_columns <- c(record_id = "ID", citation = "CITATION",
                 link = "CITATION_URL", site = "SITE_NAME",
                 country_name = "COUNTRY_NAME", latitude = "LATITUDE",
                 longitude = "LONGITUDE", year_start = "YEAR_START",
                 generation = "STAGE_ORIGIN",
                 insecticide_type = "INSECTICIDE_TYPE",
                 concentration = "INSECTICIDE_CONC", test_type = "TEST_TYPE",
                 mosquito_number = "MOSQUITO_NUMBER",
                 mortality_adjusted = "MORTALITY_ADJUSTED")

# An IR Mapper analysis-ready extract, as a source (below), given the column
# naming its species. The extracts are in Windows-1252, with "NR" and "NF" for
# missing values. Records pooled over several sites ("multipoint") are
# dropped.
ir_mapper_extract <- function(source, file, species_column) {
  force(file)
  list(
    source = source,
    database = "IR Mapper extracts",
    read = function() {
      read_csv(file, col_types = cols(.default = "c"),
               na = c("", "NA", "NR", "NF"),
               locale = locale(encoding = "windows-1252"))
    },
    keep = function(d) filter(d, `Site type` == "point"),
    columns = c(ir_mapper_extract_columns, species = species_column)
  )
}
ir_mapper_extract_columns <- c(record_id = "Bioassay ID",
                               citation = "Source 1 citation",
                               site = "Site name", country_name = "Country",
                               latitude = "Latitude", longitude = "Longitude",
                               year_start = "Start year",
                               year_end = "End year",
                               month = "Start month",
                               generation = "Generation",
                               insecticide_type = "Insecticide tested",
                               concentration = "Concentration (%)",
                               test_type = "Test protocol",
                               mosquito_number = "No. mosquitoes tested",
                               mortality_adjusted = "Percent mortality")

# Vector Atlas cites a study by DOI, so we build a citation from its author and
# publication year columns. These sometimes vary within a study (a fill-down
# error gives one study a year per record, 2014 to 2031), so we take the most
# common, and the earliest of equally common years.
read_vector_atlas <- function() {
  va <- read_csv("data/raw/VA_FULL_DATA_20250721.csv",
                 col_types = cols(.default = "c")) %>%
    filter(insecticide_resistance_data == "phenotypic")
  study_citations <- va %>%
    count(source_id, author, publication_year, article_title,
          journal_title) %>%
    arrange(source_id, desc(n), publication_year) %>%
    distinct(source_id, .keep_all = TRUE) %>%
    transmute(source_id,
              citation = paste0(author, " (", publication_year, ") ",
                                article_title, ". ", journal_title, "."))
  va %>%
    left_join(study_citations, by = "source_id")
}

# each source: its database, how to read it, its own filters, and its columns.
# Sources are listed in order of preference, which chooses a paper's version
# where versions tie (below) and the copy of a bioassay to keep. The IR Mapper
# extracts come first as the analysis-ready datasets of Moyes et al. (2019),
# abstracted and checked to one standard; where a test is in three or more
# databases, they are the least often the one database that disagrees with the
# rest (see the supplement).
sources <- list(
  ir_mapper_extract("mapper_dis_species",
                    "data/raw/2_standard-WHO-susc-test_species.csv",
                    "Species"),
  ir_mapper_extract("mapper_int_complex",
                    "data/raw/5_intensity-bioassays_complex-subgroup.csv",
                    "Complex/Subgroup"),
  ir_mapper_extract("mapper_dis_complex",
                    "data/raw/1_standard-WHO-susc-test_complex-subgroup.csv",
                    "Complex/Subgroup"),
  list(
    source = "va_new",
    database = "Vector Atlas",
    read = read_vector_atlas,
    # records at a point, with a positive number of mosquitoes if recorded
    keep = function(d) {
      d %>%
        filter(area_type == "point",
               is.na(mosquitoes_tested_n) |
                 as.numeric(mosquitoes_tested_n) != 0)
    },
    columns = c(citation = "citation", link = "citation_doi", site = "site",
                country_name = "country", latitude = "latitude_1",
                longitude = "longitude_1", year_start = "year_start",
                year_end = "year_end",
                month = "month_start", species = "species",
                generation = "generation",
                insecticide_type = "insecticide_tested",
                concentration = "concentration_percent",
                test_type = "test_protocol",
                mosquito_number = "mosquitoes_tested_n",
                mortality_adjusted = "percent_mortality",
                synergist = "synergist_tested",
                intensity = "intensity_multiplier")
  ),
  list(
    source = "mtm_dis",
    database = "MTM",
    read = function() {
      read_xlsx(
        "data/raw/MTM_DISCRIMINATING_CONCENTRATION_BIOASSAY_20261005.xlsx",
        sheet = "Data", col_types = "text")
    },
    # drop the filled mortalities (above) and a record in the Americas
    keep = function(d) {
      d %>%
        filter(!ID %in% mtm_filled_mortality_ids,
               as.numeric(LONGITUDE) > -50)
    },
    columns = c(mtm_columns, species = "SPECIES")
  ),
  list(
    source = "mtm_int",
    database = "MTM",
    read = function() {
      read_xlsx(
        "data/raw/MTM_INTENSITY_CONCENTRATION_BIOASSAY_20261005.xlsx",
        sheet = "Data", col_types = "text")
    },
    keep = function(d) filter(d, as.numeric(LONGITUDE) > -50),
    columns = c(mtm_columns, species = "VECTOR_SPECIES")
  ),
  list(
    source = "mapper_source_2022",
    database = "IR Mapper 2022",
    read = function() {
      read_xlsx("data/raw/IR Mapper data August 2022.xlsx", col_types = "text")
    },
    keep = identity,
    columns = c(citation = "Reference_Name", link = "URL", site = "Locality",
                country_name = "Country", latitude = "Latitude",
                longitude = "Longitude",
                year_start = "Mosquito_collection_Year_start",
                year_end = "Mosquito_collection_year_end",
                species = "Vector_species", insecticide_type = "Chemical_type",
                concentration = "Insecticide_dosage",
                test_type = "IR_Test_Method",
                mosquito_number = "IR_Test_NumExposed",
                mortality_adjusted = "IR_Test_mortality")
  ),
  list(
    source = "mapper_source_2024",
    database = "IR Mapper 2024",
    read = function() {
      read_xlsx("data/raw/IR_Mapper_Anopheles_data_2021_Sept_2024.xlsx",
                col_types = "text")
    },
    keep = identity,
    columns = c(citation = "reference_Name", link = "url", site = "locality",
                country_name = "country", latitude = "latitude",
                longitude = "longitude", year_start = "collection_Year_Start",
                year_end = "collection_Year_End",
                species = "vector_Species", insecticide_type = "chemical_Type",
                concentration = "insecticide_Dosage",
                test_type = "iR_Test_Method",
                mosquito_number = "iR_Test_NumExposed/Tested",
                mortality_adjusted = "iR_Test_Mortality")
  )
)

# A paper is labelled by its first author and year (citation_author_year()),
# but databases misspell surnames, split compound ones, or give a paper's online
# year or its preprint, and two papers can share a first author, year and
# country (e.g. Edi et al. 2014 in Parasites & Vectors and in PLoS Genetics).
# review/bioassay_paper_labels.csv gives these citations their correct labels,
# with the evidence.
paper_labels <- read_csv("review/bioassay_paper_labels.csv",
                         col_types = cols(.default = "c"))

bioassays <- lapply(sources, read_bioassay_source) %>%
  bind_rows() %>%
  mutate(
    record = row_number(),
    source_rank = match(source, map_chr(sources, "source")),
    country_name = standardise_country(country_name),
    insecticide_type = standardise_insecticide(insecticide_type),
    insecticide_class = insecticide_class_of(insecticide_type),
    species = standardise_species(species),
    # the paper a record cites, as its first author and year of publication
    # (see the study labels, below) unless a reviewer labelled it, and
    # whether the citation is a paper (is_paper_citation(); see one version
    # of each paper, below)
    paper = coalesce(
      paper_labels$paper[match(citation, paper_labels$citation)],
      citation_author_year(citation), citation),
    is_paper = is_paper_citation(citation),
    # the years a test ran: one year where no end year is given (as in MTM),
    # or where the end year is before the start (3 records)
    year_end = if_else(is.na(year_end) | year_end < year_start, year_start,
                       year_end),
    mosquito_number = round(mosquito_number),
    n_recorded = !is.na(mosquito_number)
  )

# ---- select records ----

# Each step keeps the records meeting a condition, and counts by source the
# records it removes
attrition <- count(bioassays, source, name = "read")
keep_records <- function(records, condition, step) {
  kept <- filter(records, {{ condition }})
  removed <- anti_join(records, kept, by = "record") %>%
    count(source, name = step)
  attrition <<- attrition %>%
    left_join(removed, by = "source") %>%
    mutate(across(all_of(step), ~ coalesce(.x, 0L)))
  kept
}

# Vector Atlas gives a test's synergist (e.g. PBO) and its intensity (e.g. 5 or
# 10 times the discriminating concentration) in their own columns, with the
# insecticide's base concentration, so a synergist or intensity test would
# look like a standard test. The other databases give no synergist tests, and
# give intensity tests at their own concentrations.
bioassays <- bioassays %>%
  keep_records(!is.na(country_name), "not in Africa") %>%
  keep_records(grepl("who", tolower(test_type)) &
                 !grepl("bottle|larval|synergist", tolower(test_type)),
               "not a WHO tube test") %>%
  keep_records((is.na(synergist) | tolower(synergist) == "none") &
                 (is.na(intensity) | intensity == "1"),
               "synergist or intensity test") %>%
  keep_records(!is.na(latitude) & !is.na(longitude) & !is.na(year_start) &
                 !is.na(concentration) & !is.na(mortality_adjusted) &
                 !is.na(insecticide_type),
               "missing a field") %>%
  keep_records(species %in% gambiae_complex_species,
               "not An. gambiae complex")

# map pixels, which the sums of species-level records (below) are matched on
mask <- rast("data/clean/raster_mask.tif")
bioassays <- bioassays %>%
  mutate(cell = terra::cellFromXY(mask, cbind(longitude, latitude))) %>%
  keep_records(!is.na(cell), "outside the map")

# ---- the number of mosquitoes, and proportions ----

# Where the number of mosquitoes isn't recorded we take the number closest to
# 100 that gives a whole number dead, or 100 if none between 50 and 200 does
# (infer_sample_size()). Copies are found from the recorded numbers only
# (n_recorded), below.
bioassays <- bioassays %>%
  mutate(
    mosquito_number = replace(
      mosquito_number, !n_recorded,
      infer_sample_size(mortality_adjusted[!n_recorded])),
    died = round(mosquito_number * mortality_adjusted / 100)
  ) %>%
  # a positive mortality that gives no dead mosquito is an error: a proportion
  # read as a percentage, e.g. 1.0 for 100% (#60), or a number dead, e.g. 3
  # of 3 or 4
  keep_records(!(died == 0 & mortality_adjusted > 0),
               "positive mortality below one mosquito")

# ---- copies of the same bioassay ----

# A download sometimes holds a paper more than once (imported twice), or
# repeats its rows (one row copied down for many tests): at least half of the
# paper's records in the download are repeated, identical in every column but
# their ID (row_hash). One of each is kept. Records of one paper in one
# download that are repeated otherwise are taken as the paper's separate tests:
# a paper reporting many tests at one site can give two the same result. The
# supplement gives the shares of repeated records and the checks by hand.
bioassays <- bioassays %>%
  group_by(source, citation) %>%
  add_count(row_hash, name = "row_copies") %>%
  mutate(repeated = mean(row_copies > 1) >= 0.5 & duplicated(row_hash)) %>%
  ungroup() %>%
  select(-row_copies) %>%
  keep_records(!repeated, "paper repeated in a download")

# The databases hold many of the same bioassays, in three ways, removed in
# this order:
# 1. within one version (a database's records of one paper, or of one other
#    citation): one test in two downloads of the database, or a sum of the
#    version's own species-level records;
# 2. one paper in two or more databases: we keep one database's version of
#    the paper, whole;
# 3. one test under two citations, e.g. a paper and a report that includes its
#    data, in one database or two: matched record by record.
# So a version of a paper is chosen before any of its records is matched with
# another citation's, and is never changed by another version of the same
# paper.
#
# Copies have the same country, years, insecticide, concentration and result
# (agreeing_pairs()), but not always the same coordinates: the databases
# georeference sites differently, and copies in different databases are often
# tens of km apart. So copies are looked for within max_km of each other in
# different databases, and within max_km_within in the same database.
#
# The distances are a choice. Matching on the percentage alone, where a number
# of mosquitoes is missing, is less specific, so its distance is shorter.
# Within a database, copies are mostly at the same coordinates.
max_km <- c(count = 100, percent = 50, misscaled = 100)
max_km_within <- 10

# Which copy to keep: one with a recorded number of mosquitoes, then from the
# first source in the list above. A version (`version_key`) is one database's
# records of one paper (its study label, in one country), however the
# database's downloads write the citation (e.g. MTM's "Sovi et al.. (2020)"
# and "Sovi et al. (2020)"), or of one other citation. Copies under different
# citations are matched by the citation as written (`source_key`), so a paper
# entered twice in one download under two spellings (e.g. "Djègbè" and
# "Djegbe") is still found.
bioassays <- bioassays %>%
  arrange(!n_recorded, source_rank, record) %>%
  mutate(preference = row_number(),
         record_key = paste(source, record_id),
         source_key = paste(database, coalesce(citation, "")),
         version_key = if_else(is_paper,
                               paste(database, paper, country_name),
                               source_key)) %>%
  arrange(record)

# A reviewer records whether a pair of records are copies ("copy") or
# separate tests ("different") in review/bioassay_duplicates.csv, by pair_id.
# A pair decided "different" is never taken as copies; any other decision
# (blank, "unsure") leaves the rules below to decide. The pairs flagged for
# review are written to outputs/review/bioassay_review_pairs.csv (below).
decisions <- read_csv("review/bioassay_duplicates.csv",
                      col_types = cols(.default = "c")) %>%
  filter(decision %in% c("copy", "different")) %>%
  select(pair_id, decision)

record_info <- bioassays %>%
  select(record, record_key, database, source, source_key, version_key,
         citation, site,
         month, generation, preference)
pairs <- agreeing_pairs(bioassays, max_km, max_km_within) %>%
  left_join(record_info, by = c("a" = "record")) %>%
  left_join(record_info, by = c("b" = "record"), suffix = c("_a", "_b")) %>%
  mutate(pair_id = pair_id(record_key_a, record_key_b),
         name_match = same_site_name(site_a, site_b),
         same_version = version_key_a == version_key_b,
         same_source = source_key_a == source_key_b) %>%
  left_join(select(decisions, pair_id, decision), by = "pair_id")

# Keep the preferred record of each pair of copies
drop_less_preferred <- function(pairs) {
  mutate(pairs,
         dropped = if_else(preference_a > preference_b, a, b),
         kept = if_else(dropped == a, b, a))
}

# ---- 1. copies within a version ----

# Two records of one version are one test in two downloads of the database
# (e.g. the IR Mapper extracts' species and complex files, or the MTM's
# discriminating-concentration and intensity downloads) if they are at the
# same coordinates in the same month and generation, paired one-to-one; else
# they are separate tests. Within one download, a paper's records are separate
# tests (bar a paper repeated in a download, above): we trust the authors and
# the databases' extractors. A version's complex-level records that sum its
# own species-level records are copies too (aggregated_copies()).
same_text <- function(x, y) {
  coalesce(as.character(x), "") == coalesce(as.character(y), "")
}
dropped_downloads <- pairs %>%
  filter(same_version, source_a != source_b, km < 0.01,
         same_text(month_a, month_b), same_text(generation_a, generation_b),
         coalesce(decision != "different", TRUE)) %>%
  one_to_one(between = c("source_a", "source_b")) %>%
  drop_less_preferred()
dropped_sums <- aggregated_copies(
  filter(bioassays, !record %in% dropped_downloads$dropped))
version_records <- filter(bioassays,
                          !record %in% c(dropped_downloads$dropped,
                                         dropped_sums))

# ---- 2. one version of each paper ----

# Many papers are in two or more databases, each database's version extracted
# separately: different sites and years, numbers tested, or errors. Matching
# copies record by record mixes these versions, e.g. one database's site with
# another's number tested, and an error in one version (e.g. every test
# entered as 100% mortality) is kept where its records don't match the
# others'. So for each paper (a study label from a citation that is a paper,
# is_paper_citation(), in one country) in two or more databases, we keep one
# database's whole version and drop the others.
#
# The version kept is the best on its numbers tested, places and years, then
# the one with most records the model uses (subset_modelled_bioassays()), then
# the first in the order of databases of the sources (paper_versions()).
database_order <- unique(map_chr(sources, "database"))
modelled <- version_records %>%
  mutate(year = (year_start + year_end) %/% 2) %>%
  subset_modelled_bioassays(mask) %>%
  pull(record)
versions <- paper_versions(version_records, modelled, database_order)

# Dropping a version loses nothing if the kept version holds its tests with
# the same results. So each record the model uses in a dropped version is
# checked against the kept version (unmatched_in_kept()): it is matched if the
# kept version has the same test (insecticide and concentration, a site within
# 10 km or of the same name within 100 km, and overlapping years) with a
# mortality within 10 percentage points (about two standard errors for a test
# of 100 mosquitoes at 50% mortality), or has a record agreeing with it exactly
# (the same number tested and number dead) within the distances used for copies
# between databases, as databases place one site differently (e.g. 20 km
# apart). A dropped version with 5 or more records unmatched either
# holds tests the kept version lacks or disagrees with it, and the data can't
# say which version is right, so the paper is flagged for review, as is a
# paper whose versions conflict (paper_versions()). Fewer than 5 unmatched
# records make little difference to the models.
#
# A reviewer records the version to keep in review/bioassay_paper_decisions.csv:
# a database, or "none" to drop the paper, where no version is usable (e.g.
# every version's coordinates are wrong). Every flagged paper needs a
# decision, so the prep warns of any without one, and keeps the rule's version
# meanwhile.
unmatched <- unmatched_in_kept(
  filter(version_records, is_paper, record %in% modelled),
  versions %>%
    filter(rank == 1) %>%
    select(paper, country_name, kept_database = database),
  copies = filter(pairs, !within, agreement != "misscaled"))
paper_decisions <- read_csv("review/bioassay_paper_decisions.csv",
                            col_types = cols(.default = "c")) %>%
  select(paper, country_name, decision)
versions <- versions %>%
  left_join(count(unmatched, paper, country_name, database,
                  name = "unmatched"),
            by = c("paper", "country_name", "database")) %>%
  mutate(unmatched = if_else(rank == 1, NA, coalesce(unmatched, 0L))) %>%
  left_join(paper_decisions, by = c("paper", "country_name")) %>%
  group_by(paper, country_name) %>%
  mutate(flagged = conflict | any(unmatched >= 5, na.rm = TRUE),
         kept = case_when(decision == "none" ~ FALSE,
                          decision %in% database ~ database == decision,
                          .default = rank == 1)) %>%
  ungroup()
# a decision must name a version of a paper in two or more databases, or none
bad_decisions <- paper_decisions %>%
  left_join(distinct(versions, paper, country_name, database),
            by = c("paper", "country_name")) %>%
  group_by(paper, country_name, decision) %>%
  summarise(ok = !anyNA(database) &
              (first(decision) == "none" | first(decision) %in% database),
            .groups = "drop") %>%
  filter(!ok)
if (nrow(bad_decisions) > 0) {
  stop("paper decisions that name no version: ",
       paste(bad_decisions$paper, bad_decisions$country_name,
             bad_decisions$decision, sep = ", ", collapse = "; "))
}
undecided <- nrow(distinct(filter(versions, flagged, is.na(decision)),
                           paper, country_name))
if (undecided > 0) {
  warning(undecided, " papers flagged for review have no decision; see ",
          "outputs/review/bioassay_paper_versions.csv")
}
dropped_versions <- version_records %>%
  filter(is_paper) %>%
  semi_join(filter(versions, !kept),
            by = c("paper", "country_name", "database")) %>%
  pull(record)
dropped_papers <- version_records %>%
  filter(is_paper) %>%
  semi_join(filter(versions, decision == "none"),
            by = c("paper", "country_name")) %>%
  pull(record)
dropped_versions <- setdiff(dropped_versions, dropped_papers)

# ---- 3. copies under different citations ----

# Among the records left (the kept version of each paper, and the records
# whose citation isn't a paper), copies under different citations, in one
# database or two, are paired one-to-one: a record pairs with at most one
# record of each other citation in each database (one_to_one()). Two
# citations whose copies include at least two results between 0% and 100%
# share one dataset (e.g. a paper and a report that includes its data), and
# all their copies are copies. Otherwise their copies are mostly at 0% or 100%,
# which separate tests often share: where two citations have 5 or more such
# copies, each is decided by hand; fewer make little difference to the
# models, and are taken as copies, as about three quarters of those decided by
# hand in October 2026 were.
left <- setdiff(version_records$record, c(dropped_versions, dropped_papers))
copies <- pairs %>%
  filter(!same_source, agreement != "misscaled", a %in% left, b %in% left) %>%
  one_to_one(between = c("source_key_a", "source_key_b")) %>%
  mutate(source_1 = pmin(source_key_a, source_key_b),
         source_2 = pmax(source_key_a, source_key_b)) %>%
  group_by(source_1, source_2) %>%
  mutate(shared = sum(!extreme), copies = n()) %>%
  ungroup() %>%
  mutate(rule = case_when(shared >= 2 ~ "shared results",
                          copies >= 5 ~ "review",
                          .default = "few copies"),
         copy = coalesce(decision == "copy", TRUE)) %>%
  drop_less_preferred()
undecided <- sum(copies$rule == "review" & is.na(copies$decision))
if (undecided > 0) {
  warning(undecided, " pairs flagged for review have no decision; see ",
          "outputs/review/bioassay_review_pairs.csv")
}
dropped_cross <- filter(copies, copy)$dropped

# A record on the wrong scale (its mortality, read as a number dead or a
# proportion, gives another database's number dead) is dropped, unless it
# also has a copy on its own scale in another database (then the match was a
# coincidence).
has_copy <- pairs %>%
  filter(!within, agreement != "misscaled") %>%
  select(a, b) %>%
  unlist() %>%
  unique()
misscaled <- pairs %>%
  filter(agreement == "misscaled", !a %in% has_copy, a %in% left,
         b %in% left) %>%
  one_to_one() %>%
  mutate(dropped = a, kept = b)

dropped_pooled <- pooled_copies(filter(bioassays, record %in% left))

# every copy dropped, with the record kept in its place; the versions dropped
# are listed in outputs/review/bioassay_paper_versions.csv
links <- bind_rows(
  mutate(dropped_downloads, rule = "one version: one test in two downloads"),
  tibble(dropped = dropped_sums,
         rule = "one version: sum of its species-level records"),
  copies %>%
    filter(copy) %>%
    mutate(rule = paste("different citations:", coalesce(decision, rule))),
  mutate(misscaled, rule = "different citations: wrong scale"),
  tibble(dropped = dropped_pooled,
         rule = "different citations: separate tests pooled into one record")
) %>%
  select(rule, km, dropped, kept)

before_copies <- bioassays
bioassays <- bioassays %>%
  keep_records(!record %in% dropped_downloads$dropped,
               "one test in two downloads") %>%
  keep_records(!record %in% dropped_sums, "sum of species records") %>%
  keep_records(!record %in% dropped_versions,
               "another database's version of the paper") %>%
  keep_records(!record %in% dropped_papers,
               "paper unusable in every version") %>%
  keep_records(!record %in% dropped_cross, "copy under another citation") %>%
  keep_records(!record %in% misscaled$dropped,
               "number dead given as mortality") %>%
  keep_records(!record %in% dropped_pooled, "pooled separate tests")

# the models use one year for each test: the middle of its years, rounded
# down, so a test over two years takes the first
bioassays <- mutate(bioassays, year = (year_start + year_end) %/% 2)

# ---- country labels ----

# the country labels must agree with the coordinates; see check_country_labels()
country_check <- check_country_labels(bioassays)
dir.create("outputs/review", showWarnings = FALSE, recursive = TRUE)
bioassays %>%
  select(longitude, latitude, year_start, insecticide_type, source,
         citation) %>%
  bind_cols(country_check) %>%
  filter(action != "keep") %>%
  write.csv("outputs/review/bioassay_country_relabels.csv", row.names = FALSE)
bioassays <- bioassays %>%
  mutate(
    country_name = if_else(country_check$action == "relabel",
                           country_check$country_geometry, country_name),
    longitude = if_else(country_check$action == "move",
                        country_check$new_longitude, longitude),
    latitude = if_else(country_check$action == "move",
                       country_check$new_latitude, latitude)
  ) %>%
  keep_records(country_check$action != "drop", "far from its country") %>%
  left_join(country_region_lookup(), by = "country_name")

# ---- study labels ----

# A label for each study, to model noise shared by the bioassays of a study
# (the survey is the study and year). The databases cite the same paper
# differently, e.g. "Dabire KR, Diabate A, ... 2009;103(11):1113-20.", "Dabire
# et al. 2009" and "Dabire et al. (2009) Distribution of ...", so we label a
# study by its paper label (above): the first author's surname and the
# publication year parsed from the citation, or the citation itself where none
# parses (institutional or missing citations such as "Ministry of Health").
# Labels are split by country, which stops institutional labels and personal
# communications of multi-country datasets pooling surveys across countries.
# See #53.
#
# Some studies' data are cited differently again (in a report and a paper, or
# by an institution and an author). Their labels are not merged: merging them
# chains through reports that compile many studies (e.g. the PMI and WHO
# reports) into one label for most of the data.
bioassays <- bioassays %>%
  mutate(study = paste(paper, country_name, sep = ", "),
         species_complex = "gambiae complex",
         n_imputed = !n_recorded,
         link = paper_link(link, citation))

# ---- records for review, and reports ----

# The pairs flagged for review, with what the reviewer needs: each record's
# citation, a link to the paper (its DOI, the database's URL, or a Google
# Scholar search), the fields that might distinguish two tests, and any
# decision already recorded. A reviewer copies a pair's row into
# review/bioassay_duplicates.csv and fills in `decision` ("copy" or
# "different") and `notes`. Undecided pairs are taken as copies.
review_fields <- c("source", "record_id", "citation", "link", "year_start",
                   "year_end", "site", "latitude", "longitude", "month",
                   "generation", "species", "mosquito_number", "n_recorded",
                   "died", "mortality_adjusted")
record_fields <- before_copies %>%
  mutate(link = paper_link(link, citation)) %>%
  select(record, insecticide_type, concentration,
         all_of(review_fields))
side_by_side <- function(pairs) {
  a <- record_fields %>%
    rename_with(~ paste0(.x, "_a"), all_of(review_fields))
  b <- record_fields %>%
    select(record, all_of(review_fields)) %>%
    rename_with(~ paste0(.x, "_b"), all_of(review_fields))
  pairs %>%
    select(any_of(c("pair_id", "reason", "database_a", "database_b", "km",
                    "shared", "a", "b"))) %>%
    left_join(a, by = c("a" = "record")) %>%
    left_join(b, by = c("b" = "record")) %>%
    select(any_of(c("pair_id", "reason", "database_a", "database_b", "km",
                    "shared", "insecticide_type",
                    "concentration")),
           all_of(as.vector(rbind(paste0(review_fields, "_a"),
                                  paste0(review_fields, "_b")))))
}
filter(copies, rule == "review") %>%
  mutate(reason = rule, km = round(km, 2)) %>%
  side_by_side() %>%
  left_join(decisions, by = "pair_id") %>%
  relocate(pair_id, decision) %>%
  arrange(!is.na(decision), reason, citation_a) %>%
  write_excel_csv("outputs/review/bioassay_review_pairs.csv", na = "")

# Each paper's versions, the one kept and why: the rule's rank and reason, the
# records of each dropped version the kept version doesn't match, whether the
# paper is flagged for review, and any decision. A reviewer adds a flagged
# paper's decision to review/bioassay_paper_decisions.csv.
versions %>%
  arrange(desc(flagged), paper, country_name, rank) %>%
  select(flagged, decision, paper, country_name, database, kept, rank, reason,
         conflict, numbers_vary, places, years, records, modelled_records,
         unmatched) %>%
  write_excel_csv("outputs/review/bioassay_paper_versions.csv", na = "")

links %>%
  rename(a = dropped, b = kept) %>%
  mutate(reason = rule, km = round(km, 2)) %>%
  side_by_side() %>%
  rename_with(~ sub("_a$", "_dropped", sub("_b$", "_kept", .x))) %>%
  write_excel_csv("outputs/review/bioassay_links.csv", na = "")

attrition <- attrition %>%
  left_join(count(bioassays, source, name = "kept"), by = "source") %>%
  arrange(match(source, map_chr(sources, "source")))
write.csv(attrition, "outputs/review/bioassay_attrition.csv", row.names = FALSE)
print(as.data.frame(attrition))

bioassays %>%
  select(mortality_adjusted, mosquito_number, died, n_imputed, year,
         year_start, year_end, month, latitude, longitude, site, species,
         species_complex, insecticide_type, insecticide_class, concentration,
         test_type, country_name, region, citation, link, study, source,
         database, record_id) %>%
  saveRDS("data/clean/all_gambiae_complex_data.RDS")
