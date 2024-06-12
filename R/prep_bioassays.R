# prep data for IR mapping

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load data from malaria threat map discrimintating concentraation bioassays and
# subset to Africa
ir_mtm_africa <- read_xlsx(
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
    ASSAY_TYPE,
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
  )

# save these out as an RDS
saveRDS(ir_mtm_africa,
        file = "data/clean/mtm_data.RDS")
