# Covariate layers and design matrix shared by the cross-validation folds.
#
# Extracted from dynamic_predictive_validation.R (#10): none of this depends on
# which fold is being fitted, so it is built once and passed to fit_fold().
# Sourcing this file expects the fold definitions in validation_folds.R to have
# been sourced already, for `df`, `unique_cells`, `classes`, `types`, `regions`
# and `countries`.

source("R/dynamical_model.R")

# build covariate rasters for the proper model

nets_cube <- pre_pad_cube(nets_cube, baseline_year)
irs_cube <- pre_pad_cube(irs_cube, baseline_year)
pop_cube <- pre_pad_cube(pop_cube, baseline_year)

nets_cube <- post_pad_cube(nets_cube, final_data_year)
irs_cube <- post_pad_cube(irs_cube, final_data_year)
pop_cube <- post_pad_cube(pop_cube, final_data_year)

crops_group <- rast("data/clean/crop_group_scaled.tif")
crops_all <- rast("data/clean/crop_scaled.tif")
crops_implicated <- c(
  crops_all$cotton,
  crops_all$vegetables,
  crops_all$rice)
covs_flat <- c(crops_group, crops_implicated)


# index to the classes for each type
classes_index <- df %>%
  distinct(type_id, class_id) %>%
  arrange(type_id) %>%
  pull(class_id)

# pull out concentrations for different types
type_concentrations <- df %>%
  select(type_id, concentration) %>%
  group_by(type_id) %>%
  filter(row_number() == 1) %>%
  arrange(type_id) %>%
  pull(concentration)


# the model terms (R/dynamical_model.R): the defaults, unless the calling
# script has set model_options before sourcing this file
if (!exists("model_options")) {
  model_options <- dynamical_model_options()
}

# create design matrix at all unique cells and for all years, as the model
# options' selection design asks (R/model_covariates.R)
selection <- selection_design_matrix(unique_cells, baseline_year,
                                     final_data_year,
                                     model_options$selection_columns)
cell_years_index <- selection$cell_years_index
x_cell_years <- selection$x_cell_years
rm(selection)

# the initial-state covariates (#19) at each cell, one row per cell_id
x_cells_init <- init_covariate_matrix(unique_cells,
                                      model_options$selection_columns)

# dimensions of things in the fitting stage
n_covs <- ncol(x_cell_years)
n_obs <- nrow(df)
n_unique_cells <- length(unique_cells)
n_times <- max(df$year_start) - min(df$year_start) + 1
n_classes <- length(classes)
n_types <- length(types)
n_regions <- length(regions)
n_countries <- length(countries)
