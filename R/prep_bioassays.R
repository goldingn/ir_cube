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
    COUNTRY_NAME,
    CITATION
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
    COUNTRY_NAME,
    CITATION
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



# load data from public IR mapper discriminating assay data at species level
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
    `Site type`,
    `Source 1 citation`
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
    mortality_adjusted = `Percent mortality`,
    citation = `Source 1 citation`
  ) %>% 
  select(-`Site type`) %>% 
  # add database name
  mutate(source = "mapper_dis_species")

# load data from public IR mapper discriminating assay data at complex level
# these are all Africa so no need for country subset
ir_mapper_dis_complex <- read_csv("data/raw/1_standard-WHO-susc-test_complex-subgroup.csv",
                                         na = c("", "NA","NR","NF")) %>% 
  # strip out invalid string in long lat
  mutate(
    Latitude = str_remove_all(Latitude, "\xa0"),
    Longitude = str_remove_all(Longitude, "\xa0")
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
    `Site type`,
    `Source 1 citation`
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
    mortality_adjusted = `Percent mortality`,
    citation = `Source 1 citation`
  ) %>% 
  select(-`Site type`) %>% 
  # add database name
  mutate(source = "mapper_dis_complex")

# load data from public IR mapper intensity assay data at species level
# these are all Africa so no need for country subset
ir_mapper_int_complex <- read_csv("data/raw/5_intensity-bioassays_complex-subgroup.csv",
                                  na = c("", "NA","NR","NF"),
                                  locale=locale(encoding="latin1") # not sure why this is needed
                                  ) %>% 
  mutate(
    Latitude = str_remove_all(Latitude, "\xa0"),
    Longitude = str_remove_all(Longitude, "\xa0")
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
    `Site type`,
    `Source 1 citation`
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
    mortality_adjusted = `Percent mortality`,
    citation = `Source 1 citation`
  ) %>% 
  select(-`Site type`) %>% 
  # add database name
  mutate(source = "mapper_int_complex")

# load source IR mapper data, note these are in a different format and likely
# contain additional data
ir_mapper_source_2022 <- read_xlsx("data/raw/IR Mapper data August 2022.xlsx") %>% 
  # enforce numeric column types
  mutate(
    across(
      c(Latitude,
        Longitude,
        IR_Test_NumExposed,
        IR_Test_mortality,
        Mosquito_collection_Year_start,
        Insecticide_dosage # coercing concentration to numeric would cause lots of NAs for non WHO bioassay entries due to units attached, we should deal with this later
        ),
      as.numeric)
  ) %>%
  # subset to required columns
  select(
    IR_Test_NumExposed,
    Mosquito_collection_Year_start,
    Latitude,
    Longitude,
    Vector_species,
    Chemical_type,
    Chemical_class,
    Insecticide_dosage,
    IR_Test_Method,
    Country,
    IR_Test_mortality,
    Reference_Name
  ) %>%
  mutate(
    # impute the number of mosquitoes from mortality rates where needed
    IR_Test_NumExposed = case_when(
      is.na(IR_Test_NumExposed) ~ infer_sample_size(IR_Test_mortality),
      .default = round(IR_Test_NumExposed)
    ), 
    # compute the integer number that died, for model
    died = round(IR_Test_NumExposed * IR_Test_mortality / 100),
    # change country name inconsistencies with mtm
    Country = if_else(Country == "Cote d Ivoire",
                      "Côte d’Ivoire",
                      Country),
    Country = if_else(Country == "United Republic of Tanzania",
                      "Tanzania",
                      Country),
    Country = if_else(Country == "Sao tome and principe",
                      "Sao Tome & Principe",
                      Country)
  ) %>%  
  # subset to africa
  filter(Country %in% country_region_lookup()$country_name) %>% 
  # drop those where resistance isn't recorded, likely PCR assays
  filter(
    !is.na(IR_Test_mortality)
  ) %>%
  # standardise names
  rename(
    mosquito_number = IR_Test_NumExposed,
    year_start = Mosquito_collection_Year_start,
    latitude = Latitude,
    longitude = Longitude,
    species = Vector_species,
    insecticide_type = Chemical_type,
    insecticide_class = Chemical_class,
    concentration = Insecticide_dosage,
    test_type = IR_Test_Method,
    country_name = Country,
    mortality_adjusted = IR_Test_mortality,
    citation = Reference_Name
  ) %>% 
  # add database name
  mutate(source = "mapper_source_2022")

# load source IR mapper data, note these are in a different format and likely
# contain additional data, this is the bit extra since 2022
ir_mapper_source_2024 <- read_xlsx("data/raw/IR_Mapper_Anopheles_data_2021_Sept_2024.xlsx") %>% 
  # enforce numeric column types
  mutate(
    across(
      c(latitude,
        longitude,
        `iR_Test_NumExposed/Tested`,
        iR_Test_Mortality,
        collection_Year_Start,
        insecticide_Dosage # coercing concentration to numeric would cause lots of NAs for non WHO bioassay entries due to units attached, we should deal with this later
      ),
      as.numeric)
  ) %>%
  # subset to required columns
  select(
    `iR_Test_NumExposed/Tested`,
    collection_Year_Start,
    latitude,
    longitude,
    vector_Species,
    chemical_Type,
    chemical_Class,
    insecticide_Dosage,
    iR_Test_Method,
    country,
    iR_Test_Mortality,
    reference_Name
  ) %>%
  mutate(
    # impute the number of mosquitoes from mortality rates where needed
    `iR_Test_NumExposed/Tested` = case_when(
      is.na(`iR_Test_NumExposed/Tested`) ~ infer_sample_size(iR_Test_Mortality),
      .default = round(`iR_Test_NumExposed/Tested`)
    ), 
    # compute the integer number that died, for model
    died = round(`iR_Test_NumExposed/Tested` * iR_Test_Mortality / 100),
    # change country name inconsistencies with mtm
    country = if_else(country == "Cote d Ivoire",
                      "Côte d’Ivoire",
                      country),
    country = if_else(country == "United Republic of Tanzania",
                      "Tanzania",
                      country),
    country = if_else(country == "Sao tome and principe",
                      "Sao Tome & Principe",
                      country)
  ) %>%  
  # subset to africa
  filter(country %in% country_region_lookup()$country_name) %>% 
  # drop those where resistance isn't recorded, likely PCR assays
  filter(
    !is.na(iR_Test_Mortality)
  ) %>%
  # standardise names
  rename(
    mosquito_number = `iR_Test_NumExposed/Tested`,
    year_start = collection_Year_Start,
    latitude = latitude,
    longitude = longitude,
    species = vector_Species,
    insecticide_type = chemical_Type,
    insecticide_class = chemical_Class,
    concentration = insecticide_Dosage,
    test_type = iR_Test_Method,
    country_name = country,
    mortality_adjusted = iR_Test_Mortality,
    citation = reference_Name
  ) %>% 
  # add database name
  mutate(source = "mapper_source_2024")
  

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
    `area type`,
    `citation_doi`
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
    mortality_adjusted = `percent mortality`,
    citation = `citation_doi`
  ) %>% 
  select(-`area type`) %>% 
  mutate(source = "va")   # add database name



# read in updated VA data
ir_va_africa_new <- read_csv("data/raw/VA_FULL_DATA_20250716.csv") %>% 
  # these were read in as characters and readxl is a pain to set column types
  # keep only bioassays
  filter(insecticide_resistance_data == "phenotypic") %>% 
  mutate(
    across(
      c(`latitude_1`,
        `longitude_1`,
        `mosquitoes_tested_n`,
        `percent_mortality`),
      as.numeric)
  ) %>%
  mutate(
    # change country name inconsistencies with mtm
    country = if_else(country == "cote divoire",
                      "Côte d’Ivoire",
                      country),
    country = if_else(country == "central african republic",
                      "CAR",
                      country),
    country = if_else(country == "sao tome & principe",
                      "Sao Tome & Principe",
                      country),
    country = if_else(country == "dr congo",
                      "DR Congo",
                      country),
    # fix incorrect class
    `insecticide class` = if_else(`insecticide_class` == "pyrethroid" &
                                    `insecticide_tested`== c("ddt"),
                                  "organochlorines",
                                  `insecticide_class`),
    `insecticide_class` = if_else(`insecticide_class` %in% c("pyrethroid","carbamate") &
                                    `insecticide_tested` == c("malathion"),
                                  "organophosphates",
                                  `insecticide_class`),
    `insecticide_class` = if_else(is.na(`insecticide_class`) &
                                    `insecticide_tested` == c("dieldrin"),
                                  "organochlorines",
                                  `insecticide_class`),
    `insecticide_class` = if_else(`insecticide_class` == "carbamate"  &
                                    `insecticide_tested` == c("fenitrothion"),
                                  "organophosphate",
                                  `insecticide_class`),
    `insecticide_class` = if_else(`insecticide_class` == "carbamate"  &
                                    `insecticide_tested` == c("pirimiphos_methyl"),
                                  "organophosphate",
                                  `insecticide_class`),
    `insecticide_class` = if_else(is.na(`insecticide_class`)  &
                                    `insecticide_tested` == c("pirimiphos_methyl"),
                                  "organophosphate",
                                  `insecticide_class`)
  ) %>% 
  select(
    `percent_mortality`,
    `mosquitoes_tested_n`,
    mosquitoes_dead_n,
    `year_start`,
    `latitude_1`,
    `longitude_1`,
    species,
    `insecticide_tested`,
    `insecticide_class`,
    `concentration_percent`,
    `test_protocol`,
    country,
    `area_type`,
    `citation_doi`
  ) %>%
  # drop multipoints for now
  filter(`area_type` == "point") %>% 
  mutate(
    # impute the number of mosquitoes from mortality rates where needed
    `mosquitoes_tested_n` = case_when(
      is.na(`mosquitoes_tested_n`) ~ infer_sample_size(`percent_mortality`),
      .default = round(`mosquitoes_tested_n`)
    ), 
    # compute the integer number that died, for model
    `mosquitoes_dead_n` = round(`mosquitoes_tested_n` * `percent_mortality` / 100)
  ) %>%  
  # drop those where resistance isn't recorded, for some reason
  filter(
    !is.na(`percent_mortality`)
  ) %>%
  # drop those where the number of mosquitoes tested is 0 (should probably be
  # NA)
  filter(
    mosquitoes_tested_n != 0
  ) %>%
  # standardise names
  rename(
    died = `mosquitoes_dead_n`,
    mosquito_number = `mosquitoes_tested_n`,
    latitude = `latitude_1`,
    longitude = `longitude_1`,
    species = species,
    insecticide_type = `insecticide_tested`,
    concentration = `concentration_percent`,
    test_type = `test_protocol`,
    country_name = country,
    mortality_adjusted = `percent_mortality`,
    citation = `citation_doi`
  ) %>% 
  select(-`area_type`) %>% 
  mutate(source = "va_new")   # add database name



# bind together data note the order matters here: species goes before complex so that when removing
# duplicates we retain the highest taxonomic resolution possible
ir_everything <- ir_va_africa %>% 
  bind_rows(ir_va_africa_new,
            ir_dis_mtm_africa,
            ir_int_mtm_africa,
            ir_mapper_dis_species_africa,
            ir_mapper_int_complex,
            ir_mapper_dis_complex,
            ir_mapper_source_2022,
            ir_mapper_source_2024)

# standardise the data
ir_everything <- ir_everything %>% 
  # convert graphical characters in Côte d'Ivoire to the format supported here
  mutate(
    country_name = iconv(country_name, to = "", sub = "byte"),
    # align these for macOS, at least
    country_name = case_when(
      country_name == "c<f4>te d'ivoire" ~ "Côte d’Ivoire",
      country_name == "C<f4>te d<92>Ivoire" ~ "Côte d’Ivoire",
      .default = country_name
    )
  ) %>%
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
                               insecticide_type),
    insecticide_type = if_else(insecticide_type == "alphacypermethrin 5x",
                               "alphacypermethrin",
                               insecticide_type),
    insecticide_type = if_else(insecticide_type == "alphacypermethrin 10x",
                               "alphacypermethrin",
                               insecticide_type),
    insecticide_type = if_else(insecticide_type == "deltamethrin 10x",
                               "deltamethrin",
                               insecticide_type),
    insecticide_type = if_else(insecticide_type == "deltamethrin10x",
                               "deltamethrin",
                               insecticide_type),
    insecticide_type = if_else(insecticide_type == "deltamethrin5x",
                               "deltamethrin",
                               insecticide_type),
    insecticide_type = if_else(insecticide_type == "lambda-cyhalothrin5x",
                               "lambdacyhalothrin",
                               insecticide_type),
    insecticide_type = if_else(insecticide_type == "lambda-cyhalothrin10x",
                               "lambdacyhalothrin",
                               insecticide_type),
    insecticide_type = if_else(insecticide_type == "permethrin 10x",
                               "permethrin",
                               insecticide_type),
    insecticide_type = if_else(insecticide_type == "permethrin 5x",
                               "permethrin",
                               insecticide_type)
  ) %>% 
  # standardise insecticide class nomenclature
  mutate(
    insecticide_class = if_else(insecticide_class == "carbamates",
                                "carbamate",
                                insecticide_class),
    insecticide_class = if_else(insecticide_class == "organochlorines",
                                "organochlorine",
                                insecticide_class),
    insecticide_class = if_else(insecticide_class == "organophosphates",
                                "organophosphate",
                                insecticide_class),
    insecticide_class = if_else(insecticide_class == "pyrethroids",
                                "pyrethroid",
                                insecticide_class),
    insecticide_class = if_else(insecticide_class == "neonicotinoids",
                                "neonicotinoid",
                                insecticide_class)
  ) %>% 
  # fill in missing class from va data
  mutate(
    insecticide_class = if_else(insecticide_type == "deltamethrin" & is.na(insecticide_class),
                                "pyrethroid",
                                insecticide_class)
  ) %>% 
  # final boss now - standardise species nomenclature
  # first get ride of all abbreviated genus name
  mutate(species = gsub("An. ","Anopheles ",species),
         # fix a missing genus name
         species = if_else(species == "arabiensis",
                           "Anopheles arabiensis",
                           species),
         species = if_else(species %in% c("stephensi",
                                          "Anopheles stephensi sl"),
                           "Anopheles stephensi",
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
                                          "Anopheles funestus s.l.",
                                          "Anopheles funestus sl",
                                          "FUNESTUS COMPLEX"),
                           "funestus complex",
                           species),
         species = if_else(species %in% c("Gambiae Complex",
                                          "GAMBIAE COMPLEX",
                                          "gambiae",
                                          "coluzzi",
                                          "Anopheles gambiae s.l.",
                                          "Anopheles gambiae sl",
                                          "gambiae (S_M)",
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
                                                  "Anopheles quadriannulatus",
                                                  "Anopheles gambiae/Anopheles coluzzii",
                                                  "Anopheles gambiae ss"
  ),
  "gambiae complex",
  NA),
  species_complex = if_else(species %in% c("funestus complex",
                                           "Anopheles funestus s.s.",
                                           "Anopheles funestus",
                                           "Anopheles funestus ss"
  ),
  "funestus complex",
  species_complex)
  )

# fix insecticide type and country names
ir_everything <- ir_everything %>%
  mutate(
    insecticide_type = case_when(
      insecticide_type == "deltamethrin" ~ "Deltamethrin",
      insecticide_type == "permethrin" ~ "Permethrin",
      insecticide_type == "ddt" ~ "DDT",
      insecticide_type == "bendiocarb" ~ "Bendiocarb",
      insecticide_type == "lambdacyhalothrin" ~ "Lambda-cyhalothrin",
      insecticide_type == "pirimiphos-methyl" ~ "Pirimiphos-methyl",
      insecticide_type == "fenitrothion" ~ "Fenitrothion",
      insecticide_type == "alphacypermethrin" ~ "Alpha-cypermethrin",
      insecticide_type == "malathion" ~ "Malathion",
      insecticide_type == "carbosulfan" ~ "Carbosulfan",
      insecticide_type == "chlorfenapyr" ~ "Chlorfenapyr",
      insecticide_type == "propoxur" ~ "Propoxur",
      insecticide_type == "etofenprox" ~ "Etofenprox",
      insecticide_type == "dieldrin" ~ "Dieldrin",
      insecticide_type == "cyfluthrin" ~ "Cyfluthrin",
      insecticide_type == "clothianidin" ~ "Clothianidin",
      insecticide_type == "chlorpyrifos-methyl" ~ "Chlorpyrifos-methyl",
      .default = NA
    ),
    insecticide_class = case_when(
      insecticide_class == "pyrethroid" ~ "Pyrethroids",
      insecticide_class == "organophosphate" ~ "Organophosphates",
      insecticide_class == "organochlorine" ~ "Organochlorines",
      insecticide_class == "carbamate" ~ "Carbamates",
      insecticide_class == "neonicotinoid" ~ "Neonicotinoids",
      insecticide_class == "pyriproxyfen" ~ "Pyriproxyfen",
      insecticide_class == "pyrroles" ~ "Pyrroles",
      .default = NA
    ),
    # capitalise countries again, then fix up errors
    country_name = str_to_title(country_name),
    country_name = str_replace(country_name, "D’ivoire", "d’Ivoire"),
    country_name = str_replace(country_name, " And ", " and "),
    country_name = str_replace(country_name, " Of ", " of "),
    country_name = str_replace(country_name, " The ", " the "),
    # shorten some country names for plotting
    country_name = case_when(
      country_name == "United Republic of Tanzania" ~ "Tanzania",
      country_name == "Democratic Republic of the Congo" ~ "DR Congo",
      country_name == "Dr Congo" ~ "DR Congo",
      country_name == "Central African Republic" ~ "CAR",
      country_name == "Sao Tome and Principe" ~ "Sao Tome & Principe",
    .default = country_name
    )
  ) %>%
  # drop any missing the insecticide type
  filter(
    !is.na(insecticide_type)
  )

# deduplicate based on unique criteria
ir_distinct <- ir_everything %>% 
  distinct(year_start,
           lat_round = round(latitude,
                             digits = 1),
           lon_round = round(longitude,
                             digits = 1),
           species_complex,
           insecticide_type,
           mosquito_number,
           died,
           country_name, 
           concentration,
           mortality_round = round(mortality_adjusted, digits = 0),
           .keep_all = TRUE) 

# # map distinct records for checking if any remaining duplicates
# library(mapview)
# library(sf)
# jitter_coord <- function(x) jitter(x,factor = 0.1)
# 
# ir_distinct_sf <- st_as_sf(ir_distinct %>% 
#                              mutate(
#                                across(
#                                  c("longitude","latitude"),
#                                  jitter_coord
#                                  )
#                              ),
#                            coords = c("longitude","latitude"),
#                            crs = st_crs(4326))
# 
# mapview(ir_distinct_sf,
#         zcol = "insecticide_type")

# # group by nearby matching points and check if nearby points might be duplicates  
# ir_distinct %>%
#   group_by(year_start, 
#            died, 
#            insecticide_type,
#            mortality_round,
#            species_complex,
#            country_name,
#            concentration,
#            round(lat_round,1),
#            round(lon_round,1),
#            mosquito_number) %>% 
#   count() %>%
#   arrange(desc(n)) %>%
#   View()


# most of the remaining close by points with identical responses have 100% mortality, they are
# likely legit separate observations

# subset to just gambiae complex
ir_distinct_gambiae <- ir_distinct %>% filter(species_complex == "gambiae complex")
# rid columns used in dedup
ir_distinct_gambiae <- ir_distinct_gambiae %>% 
  select(-c(lat_round:mortality_round))

# add on region information
ir_distinct_gambiae <- ir_distinct_gambiae %>%
  left_join(
    country_region_lookup(),
    by = "country_name"
  )

# check contributions from each database
table(ir_distinct_gambiae$source) %>% sort(decreasing = TRUE)

# print out some diagnostic checks for key columns
summary(ir_distinct_gambiae$concentration)
table(ir_distinct_gambiae$insecticide_class %>% as.factor()) %>% sort
table(ir_distinct_gambiae$insecticide_type %>% as.factor()) %>% sort
table(ir_distinct$species %>% as.factor()) %>% sort

# tabulate canonical insecticide classification and correct errors
insecticide_class_type_canon <- table(ir_distinct_gambiae$insecticide_type,ir_distinct_gambiae$insecticide_class) %>% 
  as_tibble(.name_repair = "unique") %>% 
  rename(type = ...1, class = ...2) %>% 
  arrange(desc(n)) %>% 
  filter(!duplicated(type)) %>% 
  select(!n)
  
# override class type match based on the canonical table
ir_distinct_gambiae$insecticide_class <- insecticide_class_type_canon$class[match(ir_distinct_gambiae$insecticide_type,insecticide_class_type_canon$type)] 

# check again if every combination of type and class is unique
(table(ir_distinct_gambiae$insecticide_type,ir_distinct_gambiae$insecticide_class) %>% 
  as_tibble(.name_repair = "unique") %>% 
  rename(type = ...1, class = ...2) %>% distinct(type, class) %>% nrow()) == (table(ir_distinct_gambiae$insecticide_type,ir_distinct_gambiae$insecticide_class) %>% 
  as_tibble(.name_repair = "unique") %>% 
  rename(type = ...1, class = ...2) %>% distinct(type, class) %>% nrow())

# # save the diagnostic interactive map
# mapshot(mapview(ir_distinct_sf,zcol = "insecticide_type"),url = "distinct_pts.html")

# save these out as an RDS
saveRDS(ir_distinct_gambiae,
        file = "data/clean/all_gambiae_complex_data.RDS")

