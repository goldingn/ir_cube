# prep data for IR mapping

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load data from malaria threat map discrimintating concentraation bioassays and
# subset to Africa
ir_dis_mtm_africa <- read_xlsx(
  path = "data/raw/MTM_DISCRIMINATING_CONCENTRATION_BIOASSAY_20240610.xlsx",
  sheet = "Data") %>%
  # these were read in as characters and readxl is a pain to set column types
  mutate(
    across(
      c(LATITUDE,
        LONGITUDE,
        MOSQUITO_NUMBER,
        MORTALITY_ADJUSTED),
      as.numeric)
  ) %>%
  # subset to Africa
  filter(
    COUNTRY_NAME %in% mtm_africa_countries()
  ) %>%
  # throw out a dodgy coordinate
  filter(
    LONGITUDE > -50
  ) %>%
  # subset to required columns
  select(
    MOSQUITO_NUMBER,
    MORTALITY_ADJUSTED,
    YEAR_START,
    LATITUDE,
    LONGITUDE,
    SPECIES,
    INSECTICIDE_TYPE,
    INSECTICIDE_CLASS,
    INSECTICIDE_CONC,
    TEST_TYPE,
    COUNTRY_NAME
  ) %>%
  mutate(
    # impute the number of mosquitoes from mortality rates where needed
    MOSQUITO_NUMBER = case_when(
      is.na(MOSQUITO_NUMBER) ~ infer_sample_size(MORTALITY_ADJUSTED),
      .default = round(MOSQUITO_NUMBER)
    ), 
    # compute the integer number that died, for model
    DIED = round(MOSQUITO_NUMBER * MORTALITY_ADJUSTED / 100)
  ) %>%  
  # drop those where resistance isn't recorded, for some reason
  filter(
    !is.na(MORTALITY_ADJUSTED)
  ) %>%
  # remove these shouty caps
  rename_all(
    tolower
  )  %>% 
  # add database name
  mutate(source = "mtm_dis") %>%   
  # coerce concentration to numeric
  mutate(
    concentration = as.numeric(str_remove(insecticide_conc, "%")
                               )
    # note this will cause coercion NAs because some would have mu-gram units or be not recorded, we
    # could remove later if needed
  ) %>% 
  select(-insecticide_conc)


# load data from malaria threat map intensity concentration bioassays and
# subset to Africa
ir_int_mtm_africa <- read_xlsx(
  path = "data/raw/MTM_INTENSITY_CONCENTRATION_BIOASSAY_20240612.xlsx",
  sheet = "Data") %>%
  # these were read in as characters and readxl is a pain to set column types
  mutate(
    across(
      c(LATITUDE,
        LONGITUDE,
        MOSQUITO_NUMBER,
        MORTALITY_ADJUSTED),
      as.numeric)
  ) %>%
  # subset to Africa
  filter(
    COUNTRY_NAME %in% mtm_africa_countries()
  ) %>%
  # probably don't need this but won't hurt
  filter(
    LONGITUDE > -50
  ) %>%
  # subset to required columns
  select(
    MOSQUITO_NUMBER,
    MORTALITY_ADJUSTED,
    YEAR_START,
    LATITUDE,
    LONGITUDE,
    VECTOR_SPECIES,
    INSECTICIDE_TYPE,
    INSECTICIDE_CLASS,
    INSECTICIDE_CONC,
    #INSECTICIDE_INTENSITY, # don't need this because it is just a multiplier of discrimination
    #concentration
    TEST_TYPE,
    COUNTRY_NAME
  ) %>%
  mutate(
    # impute the number of mosquitoes from mortality rates where needed
    MOSQUITO_NUMBER = case_when(
      is.na(MOSQUITO_NUMBER) ~ infer_sample_size(MORTALITY_ADJUSTED),
      .default = round(MOSQUITO_NUMBER)
    ), 
    # compute the integer number that died, for model
    DIED = round(MOSQUITO_NUMBER * MORTALITY_ADJUSTED / 100)
  ) %>%  
  # drop those where resistance isn't recorded, for some reason
  filter(
    !is.na(MORTALITY_ADJUSTED)
  ) %>%
  # remove these shouty caps
  rename_all(
    tolower
  ) %>% 
  # add database name
  mutate(source = "mtm_int") %>%   
  # coerce concentration to numeric
  mutate(
    concentration = as.numeric(str_remove(insecticide_conc, "%")
    )
    # note this will cause coercion NAs because some would have mu-gram units or be not recorded, we
    # could remove later if needed
  ) %>% 
  select(-insecticide_conc) %>% 
  # standardise species column name
  rename(species = "vector_species")

# load data from IR mapper discriminating assay data at species level
# these are all Africa so no need for country subset
ir_mapper_dis_species_africa <- read_csv("data/raw/2_standard-WHO-susc-test_species.csv",
                                         na = c("", "NA","NR","NF")) %>% 
  # these were read in as characters and readxl is a pain to set column types
  mutate(
    across(
      c(Latitude,
        Longitude,
        `No. mosquitoes dead`,
        `No. mosquitoes tested`,
        `Start year`),
        as.numeric)
  ) %>%
  # subset to required columns
  select(
    `No. mosquitoes dead`,
    `No. mosquitoes tested`,
    `Start year`,
    Latitude,
    Longitude,
    Species,
    `Insecticide tested`,
    `Insecticide class`,
    `Concentration (%)`,
    `Test protocol`,
    Country,
    `Percent mortality`,
    `Site type`
  ) %>%
  # drop multipoints for now
  filter(`Site type` == "point") %>% 
  mutate(
    # impute the number of mosquitoes from mortality rates where needed
    `No. mosquitoes tested` = case_when(
      is.na(`No. mosquitoes tested`) ~ infer_sample_size(`Percent mortality`),
      .default = round(`No. mosquitoes tested`)
    ), 
    # compute the integer number that died, for model
    `No. mosquitoes dead` = round(`No. mosquitoes tested` * `Percent mortality` / 100),
    # change country name inconsistencies with mtm
    Country = if_else(Country == "C\xf4te d\x92Ivoire",
                      "Côte d’Ivoire",
                      Country),
    Country = if_else(Country == "Tanzania, United Republic of",
                      "United Republic of Tanzania",
                      Country),
    Country = if_else(Country == "The Gambia",
                      "Gambia",
                      Country)
  ) %>%  
  # drop those where resistance isn't recorded, for some reason
  filter(
    !is.na(`Percent mortality`)
  ) %>%
  # standardise names
  rename(
    died = `No. mosquitoes dead`,
    mosquito_number = `No. mosquitoes tested`,
    year_start = `Start year`,
    latitude = Latitude,
    longitude = Longitude,
    species = Species,
    insecticide_type = `Insecticide tested`,
    insecticide_class = `Insecticide class`,
    concentration = `Concentration (%)`,
    test_type = `Test protocol`,
    country_name = Country,
    mortality_adjusted = `Percent mortality`
  ) %>% 
  select(-`Site type`) %>% 
  # add database name
  mutate(source = "mapper_dis_species")

# load data from IR mapper discriminating assay data at complex level
# these are all Africa so no need for country subset
ir_mapper_dis_complex <- read_csv("data/raw/1_standard-WHO-susc-test_complex-subgroup.csv",
                                         na = c("", "NA","NR","NF")) %>% 
  #strip out invalid string in long lat
  mutate(Latitude = stringr::str_extract(
    Latitude,
    "\\d+\\.*\\d*"
  ),
  Longitude = stringr::str_extract(
    Longitude,
    "\\d+\\.*\\d*"
  )
  ) %>% 
  # these were read in as characters and readxl is a pain to set column types
  mutate(
    across(
      c(Latitude,
        Longitude,
        `No. mosquitoes dead`,
        `No. mosquitoes tested`,
        `Start year`),
      as.numeric)
  ) %>%
  # subset to required columns
  select(
    `No. mosquitoes dead`,
    `No. mosquitoes tested`,
    `Start year`,
    Latitude,
    Longitude,
    `Complex/Subgroup`,
    `Insecticide tested`,
    `Insecticide class`,
    `Concentration (%)`,
    `Test protocol`,
    Country,
    `Percent mortality`,
    `Site type`
  ) %>%
  # drop multipoints for now
  filter(`Site type` == "point") %>% 
  mutate(
    # impute the number of mosquitoes from mortality rates where needed
    `No. mosquitoes tested` = case_when(
      is.na(`No. mosquitoes tested`) ~ infer_sample_size(`Percent mortality`),
      .default = round(`No. mosquitoes tested`)
    ), 
    # compute the integer number that died, for model
    `No. mosquitoes dead` = round(`No. mosquitoes tested` * `Percent mortality` / 100),
    # change country name inconsistencies with mtm
    Country = if_else(Country == "C\xf4te d\x92Ivoire",
                      "Côte d’Ivoire",
                      Country),
    Country = if_else(Country == "Tanzania, United Republic of",
                      "United Republic of Tanzania",
                      Country),
    Country = if_else(Country == "The Gambia",
                      "Gambia",
                      Country)
  ) %>%  
  # drop those where resistance isn't recorded, for some reason
  filter(
    !is.na(`Percent mortality`)
  ) %>%
  # standardise names
  rename(
    died = `No. mosquitoes dead`,
    mosquito_number = `No. mosquitoes tested`,
    year_start = `Start year`,
    latitude = Latitude,
    longitude = Longitude,
    species = `Complex/Subgroup`,
    insecticide_type = `Insecticide tested`,
    insecticide_class = `Insecticide class`,
    concentration = `Concentration (%)`,
    test_type = `Test protocol`,
    country_name = Country,
    mortality_adjusted = `Percent mortality`
  ) %>% 
  select(-`Site type`) %>% 
  # add database name
  mutate(source = "mapper_dis_complex")

# load data from IR mapper intensity assay data at species level
# these are all Africa so no need for country subset
ir_mapper_int_complex <- read_csv("data/raw/5_intensity-bioassays_complex-subgroup.csv",
                                  na = c("", "NA","NR","NF"),
                                  locale=locale(encoding="latin1") # not sure why this is needed
                                  ) %>% 
  #strip out invalid string in long lat
  mutate(Latitude = stringr::str_extract(
    Latitude,
    "\\d+\\.*\\d*"
  ),
  Longitude = stringr::str_extract(
    Longitude,
    "\\d+\\.*\\d*"
  )
  ) %>% 
  # these were read in as characters and readxl is a pain to set column types
  mutate(
    across(
      c(Latitude,
        Longitude,
        `No. mosquitoes dead`,
        `No. mosquitoes tested`,
        `Start year`),
      as.numeric)
  ) %>%
  # subset to required columns
  select(
    `No. mosquitoes dead`,
    `No. mosquitoes tested`,
    `Start year`,
    Latitude,
    Longitude,
    `Complex/Subgroup`,
    `Insecticide tested`,
    `Insecticide class`,
    `Concentration (%)`,
    `Test protocol`,
    Country,
    `Percent mortality`,
    `Site type`
  ) %>%
  # drop multipoints for now
  filter(`Site type` == "point") %>% 
  mutate(
    # impute the number of mosquitoes from mortality rates where needed
    `No. mosquitoes tested` = case_when(
      is.na(`No. mosquitoes tested`) ~ infer_sample_size(`Percent mortality`),
      .default = round(`No. mosquitoes tested`)
    ), 
    # compute the integer number that died, for model
    `No. mosquitoes dead` = round(`No. mosquitoes tested` * `Percent mortality` / 100),
    # change country name inconsistencies with mtm
    Country = if_else(Country == "C\xf4te d\x92Ivoire",
                      "Côte d’Ivoire",
                      Country),
    Country = if_else(Country == "Côte d\u0092Ivoire",
                      "Côte d’Ivoire",
                      Country),
    Country = if_else(Country == "Tanzania, United Republic of",
                      "United Republic of Tanzania",
                      Country),
    Country = if_else(Country == "The Gambia",
                      "Gambia",
                      Country)
  ) %>%  
  # drop those where resistance isn't recorded, for some reason
  filter(
    !is.na(`Percent mortality`)
  ) %>%
  # standardise names
  rename(
    died = `No. mosquitoes dead`,
    mosquito_number = `No. mosquitoes tested`,
    year_start = `Start year`,
    latitude = Latitude,
    longitude = Longitude,
    species = `Complex/Subgroup`,
    insecticide_type = `Insecticide tested`,
    insecticide_class = `Insecticide class`,
    concentration = `Concentration (%)`,
    test_type = `Test protocol`,
    country_name = Country,
    mortality_adjusted = `Percent mortality`
  ) %>% 
  select(-`Site type`) %>% 
  # add database name
  mutate(source = "mapper_int_complex")

# load vector atlas

ir_va_africa <- read_csv("data/raw/VA Bioassay Data 240612.csv",col_select = 1:98) %>% 
  # these were read in as characters and readxl is a pain to set column types
  mutate(
    across(
      c(`latitude_ 1`,
        `longitude_1`,
        `mosquitoes tested_n`,
        `percent mortality`),
      as.numeric)
  ) %>%
  mutate(
    # change country name inconsistencies with mtm
    country = if_else(country == "cote divoire",
                      "Côte d’Ivoire",
                      country),
    country = if_else(country == "c\xf4te d'ivoire",
                      "Côte d’Ivoire",
                      country),
    country = if_else(country == "sao tome & principe",
                      "Sao Tome and Principe",
                      country),
    country = if_else(country == "tanzania",
                      "United Republic of Tanzania",
                      country),
    # fix incorrect class
    `insecticide class` = if_else(`insecticide class` == "pyrethroid" &
                                    `insecticide tested`== c("ddt"),
                                  "organochlorines",
                                  `insecticide class`),
    `insecticide class` = if_else(`insecticide class` == "pyrethroid" &
                                    `insecticide tested` == c("malathion"),
                                  "organophosphates",
                                  `insecticide class`)
  ) %>% 
  select(
    `mosquitoes tested_n`,
    `percent mortality`,
    `mosquitoes dead_n`,
    `year_start`,
    `latitude_ 1`,
    `longitude_1`,
    species,
    `insecticide tested`,
    `insecticide class`,
    `concentration_percent`,
    `test protocol`,
    country,
    `percent mortality`,
    `area type`
  ) %>%
  # drop multipoints for now
  filter(`area type` == "point") %>% 
  mutate(
    # impute the number of mosquitoes from mortality rates where needed
    `mosquitoes tested_n` = case_when(
      is.na(`mosquitoes tested_n`) ~ infer_sample_size(`percent mortality`),
      .default = round(`mosquitoes tested_n`)
    ), 
    # compute the integer number that died, for model
    `mosquitoes dead_n` = round(`mosquitoes tested_n` * `percent mortality` / 100)
    # note 350 entries actually have neither # dead or # total, so this cannot be inferred
  ) %>%  
  # drop those where resistance isn't recorded, for some reason
  filter(
    !is.na(`percent mortality`)
  ) %>%
  # standardise names
  rename(
    died = `mosquitoes dead_n`,
    mosquito_number = `mosquitoes tested_n`,
    latitude = `latitude_ 1`,
    longitude = `longitude_1`,
    species = species,
    insecticide_type = `insecticide tested`,
    insecticide_class = `insecticide class`,
    concentration = `concentration_percent`,
    test_type = `test protocol`,
    country_name = country,
    mortality_adjusted = `percent mortality`
  ) %>% 
  select(-`area type`) %>% 
  mutate(source = "va")   # add database name


# bind together data note the order matters here: species goes before complex so that when removing
# duplicates we retain the highest taxonomic resolution possible
ir_everything <- ir_dis_mtm_africa %>% 
  bind_rows(ir_int_mtm_africa,
            ir_va_africa,
            ir_mapper_dis_species_africa,
            ir_mapper_int_complex,
            ir_mapper_dis_complex)

# standardise the data
ir_everything <- ir_everything %>% 
  #standardise all non-species name characters to lower
  mutate(
    across(c("insecticide_type",
             "insecticide_class",
             "country_name"),
           tolower)
  ) %>% 
  # standardise insecticide type nomenclature
  mutate(
    insecticide_type = if_else(insecticide_type == "alpha-cypermethrin",
                               "alphacypermethrin",
                               insecticide_type),
    insecticide_type = if_else(insecticide_type == "lambda-cyhalothrin",
                               "lambdacyhalothrin",
                               insecticide_type),
    insecticide_type = if_else(insecticide_type == "pirimiphos_methyl",
                               "pirimiphos-methyl",
                               insecticide_type)
  ) %>% 
  # standardise insecticide class nomenclature
  mutate(
    insecticide_class = if_else(insecticide_class == "carbamate",
                                "carbamates",
                                insecticide_class),
    insecticide_class = if_else(insecticide_class == "organochlorine",
                                "organochlorines",
                                insecticide_class),
    insecticide_class = if_else(insecticide_class == "organophosphate",
                                "organophosphates",
                                insecticide_class),
    insecticide_class = if_else(insecticide_class == "pyrethroid",
                                "pyrethroids",
                                insecticide_class),
  ) %>% 
  # fill in missing class from va data
  mutate(
    insecticide_class = if_else(insecticide_type == "deltamethrin" & is.na(insecticide_class),
                                "pyrethroids",
                                insecticide_class)
  ) %>% 
  # final boss now - standardise species nomenclature
  # first get ride of all abbreviated genus name
  mutate(species = gsub("An. ","Anopheles ",species),
         # fix a missing genus name
         species = if_else(species == "arabiensis",
                           "Anopheles arabiensis",
                           species),
         species = if_else(species == "coluzzii (gambiae m)",
                           "Anopheles coluzzii",
                           species),
         species = if_else(species == "gambiae (s)",
                           "Anopheles gambiae s.s.",
                           species),
         # standardise complex names
         species = if_else(species %in% c("Funestus Subgroup",
                                          "funestus",
                                          "Anopheles funestus s.l."),
                           "funestus complex",
                           species),
         species = if_else(species %in% c("Gambiae Complex",
                                          "Anopheles gambiae s.l.",
                                          "Anopheles coluzzii/gambiae"),
                           "gambiae complex",
                           species)
  ) 

# filter data
ir_everything <- ir_everything %>% 
  # drop CDC bioassay
  mutate(
    test_type = tolower(test_type)
  ) %>% 
  filter(grepl("who",test_type)) %>%
  # filter out data without basic info
  # na coords (from ir_mapper)
  filter(!(is.na(latitude) | is.na(longitude))) %>% 
  # filter out those without inferrable dead or total mosquito counts
  # shouldn't be any here because of imputation, but just in case
  filter(!(is.na(mosquito_number) | is.na(died))) %>% 
  # filter out na years
  filter(!(is.na(year_start))) %>% 
  # filter out na species
  filter(!(is.na(species))) %>% 
  #filter out na concentration
  filter(!(is.na(concentration)))
  
# check concentration and insecticide mapping
type_concentrations <- ir_everything %>%
  select(insecticide_type, concentration,source) %>%
  group_by(insecticide_type) %>%
  unique() %>%
  arrange(insecticide_type)
#View(type_concentrations)

ir_everything <- ir_everything %>% 
  # assign species to complexes for the purpose of deduplication
  mutate(species_complex = if_else(species %in% c("gambiae complex",
                                                  "Anopheles gambiae s.s.",
                                                  "Anopheles arabiensis",
                                                  "Anopheles coluzzii",
                                                  "Anopheles merus",
                                                  "Anopheles gambiae",
                                                  "Anopheles quadriannulatus"
  ),
  "gambiae complex",
  NA),
  species_complex = if_else(species %in% c("funestus complex",
                                           "Anopheles funestus s.s.",
                                           "Anopheles funestus"
  ),
  "funestus complex",
  species_complex)
  )

# deduplicate based on unique criteria
ir_distinct <- ir_everything %>% 
  distinct(year_start,
           lat_round = round(latitude,digits = 2),
           lon_round = round(longitude,digits = 2),
           species_complex,
           insecticide_type,
           mosquito_number,
           died,
           country_name, 
           concentration,
           mortality_round = round(mortality_adjusted,digits = 0),
           .keep_all = TRUE) 

# map distinct records for checking if any remaining duplicates
library(mapview)
library(sf)
jitter_coord <- function(x) jitter(x,factor = 0.1)

ir_distinct_sf <- st_as_sf(ir_distinct %>% 
                             mutate(
                               across(
                                 c("longitude","latitude"),
                                 jitter_coord
                                 )
                             ),
                           coords = c("longitude","latitude"),
                           crs = st_crs(4326))


mapview(ir_distinct_sf,zcol = "insecticide_type")

# group by nearby matching points and check if nearby points might be duplicates  
ir_distinct %>%  group_by(year_start, 
                          died, 
                          insecticide_type,
                          mortality_round,
                          species_complex,
                          country_name,
                          concentration,
                          round(lat_round,1),
                          round(lon_round,1),
                          mosquito_number) %>% 
  count() %>% View
# most of the remaining close by points with identical responses have 100% mortality, they are
# likely legit separate observations

# subset to just gambiae complex
ir_distinct_gambiae <- ir_distinct %>% filter(species_complex == "gambiae complex")
# rid columns used in dedup
ir_distinct_gambiae <- ir_distinct_gambiae %>% 
  select(-c(lat_round:mortality_round))

# check contributions from each database
table(ir_distinct_gambiae$source) %>% sort(decreasing = TRUE)

# # save the diagnostic interactive map
# mapshot(mapview(ir_distinct_sf,zcol = "insecticide_type"),url = "distinct_pts.html")
# save these out as an RDS
saveRDS(ir_distinct_gambiae,
        file = "data/clean/all_gambiae_complex_data.RDS")
