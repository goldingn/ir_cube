# The design matrix and initial-state covariates shared by the
# cross-validation folds (#10): none of this depends on which fold is being
# fitted, so it is built once and passed to fit_fold(). Source after
# validation_folds.R, for `df`, `unique_cells`, `baseline_year` and
# `final_data_year`.

source("R/dynamical_model.R")

# the default model terms, or an R expression for them in
# IR_CUBE_MODEL_OPTIONS (docker/run_pod_job.sh), unless the calling script has
# set model_options. A fold paired with a full fit made with other options
# needs that fit's options
if (!exists("model_options")) {
  model_options <- eval(str2lang(Sys.getenv("IR_CUBE_MODEL_OPTIONS",
                                            "dynamical_model_options()")))
}

# the design matrix at all unique cells and years (R/model_covariates.R)
selection <- selection_design_matrix(unique_cells, baseline_year,
                                     final_data_year,
                                     model_options$selection_columns)
cell_years_index <- selection$cell_years_index
x_cell_years <- selection$x_cell_years
rm(selection)

# the initial-state covariates (#19) at each cell, one row per cell_id
x_cells_init <- init_covariate_matrix(unique_cells,
                                      model_options$selection_columns)

