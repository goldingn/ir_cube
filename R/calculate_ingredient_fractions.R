# Compute weights for an 'effective resistance to nets' output metric

# Analyse DHS/MIS survey data to approximate the proportion of single-AI LLINs
# that have each of the major pyrethroids (deltamethrin, permethrin,
# alpha-cypermethrin) as the active ingredient.

# load packages and functions
source("R/packages.R")
source("R/functions.R")

# load the DHS/MIS data,
survey_net_data <- read_csv("data/raw/nets_by_ingredients.csv",
                            col_types = cols(
                              program = col_character(),
                              region = col_character(),
                              surveyname = col_character(),
                              country = col_character(),
                              year = col_double(),
                              brandid = col_double(),
                              brand = col_character(),
                              brand_reclass = col_character(),
                              net_type = col_character(),
                              active_ingredient = col_character(),
                              count = col_double()
                            ))

# select only those records where the active ingredient (AI) was listed. Note
# that this is a very imperfect approach, since nets are much more likely to be
# listed with insufficient information to infer the AI, e.g. listing the net
# brand name (e.g. Permanet), but not the version/generation, thereby making it
# ambiguous whether the net is single or dual action, or the AI. This approach
# assumes the information is 'missing at random', i.e. or those where it is not
# listed, it is equally like to follow these proportions
pyrethroids <- c("Deltamethrin", "Permethrin", "Alpha-cypermethrin")

# compute weights for the past 5 years
weights <- survey_net_data %>%
  filter(
    year >= 2022 - 5,
    active_ingredient %in% pyrethroids
  ) %>%
  group_by(
    active_ingredient
  ) %>%
  summarise(
    count = sum(count),
    .groups = "drop"
  ) %>%
  mutate(
    weight = count / sum(count)
  ) %>%
  select(
    -count
  ) %>%
  pivot_wider(
    names_from = active_ingredient,
    values_from = weight
  ) %>%
  as.list()

# stash these to use in prediction
saveRDS(weights,
        "temporary/ingredient_weights.RDS")

# print the percentages
perc <- round(unlist(weights) * 100)
paste(sprintf("%s: %i%s", names(perc), perc, "%"), collapse = ", ")
