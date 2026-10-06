# The fits of the species model (#47) and its reference, full data only. Run
# on RunPod (docker/README.md), one fit per pod.
#
#   Rscript R/species_runs.R [code ref]
#
# writes outputs/species_runs/jobs.csv, one row per fit, and
# outputs/species_runs/pods.json, the create-pod body of each
# (net_grid_pod_body(), R/net_grid.R), for the commit to run (a full commit
# id; default: HEAD, which must be on GitHub). The fits:
#   sp_ref_floor  B_f: one trajectory, the floor estimated with prior
#                 Beta(1, 4), the prior of the species floors
#   sp_v1         V1: the species model, no floors
#   sp_v1_floor   V1f: the species model, a floor for each species, each
#                 Beta(1, 4)
#   sp_v2         V2: one trajectory, with the kdr covariate (kdr_options(),
#                 R/kdr_covariate.R; the "complex" band), no floor
#   sp_v2_floor   V2f: as V2, the floor estimated with prior Beta(1, 4)
#   sp_v3         V3: the species model with the kdr covariate, each species
#                 with its own band and slopes, no floors
#   sp_v3_floor   V3f: as V3, a floor for each species, each Beta(1, 4)
#   sp_v4         V4: V2 with the kdr-dependent floor, plogis(floor_intercept
#                 + floor_kdr k(x)), the intercept's prior matched to
#                 Beta(1, 4) (kdr_options(floor = TRUE))
#   sp_v4_class   V4 class: as V4, with an intercept per insecticide class and
#                 the kdr term for the pyrethroids and DDT only
#                 (kdr_options(floor = "class"))
# Otherwise the default options (d_half 270, reversion estimated), given in
# full so that a later change of default does not change these fits. The fits
# with floors start chains 1-2 from the low-floor mode and 3-4 from the
# high-floor mode (IR_CUBE_INITS, dynamical_chain_inits(); with the species
# model, both species' floors start at the cached floor, dynamical_inits()),
# and record each chain's floors and log posterior (chain_floor_modes()); V1
# starts from the default initial values. Before creating the pods:
#   Rscript R/floor_mode_inits.R <r2_main.RData>  # temporary/inits_floor_*.RDS
#   irpod sync    # uploads them, and data/clean/arabiensis_fraction.tif
# Then create one pod per element of the json (about 3.7 h and $1.05 per fit
# on the default pod, by the gradient timing, which is the same with and
# without the species model; the full fit of 5 October took 6.2 h; state the
# price first), and when each is done, `irpod fetch <name>` and delete its
# pod (it does not delete itself).

source("R/net_grid.R")

floor_mode_inits <- paste("temporary/inits_floor_low.RDS",
                          "temporary/inits_floor_high.RDS", sep = ",")
species_floors <- "species_options(floors = TRUE, floor_prior = c(1, 4))"
species_runs <- data.frame(
  name = c("sp_ref_floor", "sp_v1", "sp_v1_floor", "sp_v2", "sp_v2_floor",
           "sp_v3", "sp_v3_floor", "sp_v4", "sp_v4_class"),
  label = c("B_f", "V1", "V1f", "V2", "V2f", "V3", "V3f", "V4", "V4_class"),
  options = sprintf("dynamical_model_options(%s)", c(
    "mortality_floor = TRUE, floor_prior = c(1, 4)",
    "species = species_options(floors = FALSE)",
    sprintf("species = %s", species_floors),
    "kdr = kdr_options()",
    "mortality_floor = TRUE, floor_prior = c(1, 4), kdr = kdr_options()",
    "species = species_options(floors = FALSE), kdr = kdr_options()",
    sprintf("species = %s, kdr = kdr_options()", species_floors),
    paste("mortality_floor = TRUE, floor_prior = c(1, 4),",
          "kdr = kdr_options(floor = TRUE)"),
    paste("mortality_floor = TRUE, floor_prior = c(1, 4),",
          "kdr = kdr_options(floor = \"class\")"))),
  inits = c(floor_mode_inits, "", floor_mode_inits, "", floor_mode_inits, "",
            floor_mode_inits, floor_mode_inits, floor_mode_inits))

# One row per pod job, as net_grid_jobs() makes them
species_run_jobs <- function(code_ref, threads = 8) {
  jobs <- species_runs
  jobs$JOB <- sprintf("%s full --threads %d", jobs$name, threads)
  jobs$CODE_REF <- code_ref
  jobs$IR_CUBE_MODEL_OPTIONS <- jobs$options
  jobs$IR_CUBE_INITS <- jobs$inits
  jobs[, c("name", "label", "JOB", "CODE_REF", "IR_CUBE_MODEL_OPTIONS",
           "IR_CUBE_INITS")]
}

if (sys.nframe() == 0) {
  arguments <- commandArgs(trailingOnly = TRUE)
  code_ref <- if (length(arguments) >= 1) arguments[1] else
    system("git rev-parse HEAD", intern = TRUE)
  stopifnot(grepl("^[0-9a-f]{40}$", code_ref))
  # each fit's options must pass the model's checks (definitions only; no
  # data or python)
  source("R/dynamical_model.R")
  for (expression in species_runs$options) {
    check_dynamical_model_options(eval(str2lang(expression)))
  }
  inputs <- c(arabiensis_fraction_file, kdr_total_file,
              unlist(strsplit(species_runs$inits[nzchar(species_runs$inits)],
                              ",")))
  missing <- inputs[!file.exists(inputs)]
  if (length(missing) > 0) {
    warning("not found here, so not uploaded by irpod sync: ",
            toString(unique(missing)))
  }
  jobs <- species_run_jobs(code_ref)
  dir.create("outputs/species_runs", showWarnings = FALSE, recursive = TRUE)
  write.csv(jobs, "outputs/species_runs/jobs.csv", row.names = FALSE)
  bodies <- lapply(seq_len(nrow(jobs)), function(i) {
    net_grid_pod_body(jobs[i, ])
  })
  jsonlite::write_json(bodies, "outputs/species_runs/pods.json",
                       auto_unbox = TRUE, pretty = TRUE)
  cat(sprintf("%d fits at %s\n", nrow(jobs), code_ref))
  print(jobs[, c("name", "label", "JOB", "IR_CUBE_MODEL_OPTIONS",
                 "IR_CUBE_INITS")], row.names = FALSE)
}
