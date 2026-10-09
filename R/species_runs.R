# The pod jobs of V5f (#47), the default model: its test fits, its
# cross-validation, the second stage on each fold, and the two-stage maps of
# its full fit, on RunPod (docker/README.md; doc/v5f_cv_runbook.md), one job
# per pod.
#
#   Rscript R/species_runs.R [code ref] [set]
#
# writes outputs/species_runs/jobs_<set>.csv, one row per job, and
# outputs/species_runs/pods_<set>.json, the create-pod body of each
# (species_run_pod_body()), for the commit to run (a full commit id;
# default: HEAD, which must be on GitHub) and a set of run_sets, below:
# "v5f" (the default), "cv5f", "cv5f_warm", "ts_cv5f", "ts_cv5f_warm" or
# "maps_v5f". The fits of #47 before V5f (the species model, the kdr
# covariate, the weighted binomial likelihood, V5, V5r, V5i and V5h; the
# sets "species", "wb", "v5r", "v5i" and "v5h") were launched by this script
# on the branch weighted-binomial, which keeps them. Before creating the
# pods, `irpod sync` uploads the initial values (the low-floor ones,
# temporary/inits_floor_low.RDS, from R/floor_mode_inits.R). Then create
# one pod per element of the json (state the price first), and when each is
# done, `irpod fetch <name>` and delete its pod.

# The set "v5f" (v5f_runs): V5f, the default model since #47, V5h with the
# floor smooth centred over the cells with bioassays of the pyrethroids or
# DDT (smooth_centre(), R/latent_smooth.R) and the floor at a flat smooth
# with a half-normal prior of scale 0.05 (smooth_options(
# floor_intercept_prior = )), its options in full (v5f_options, which must
# equal dynamical_model_options() at the commit run): the full fit and the
# temporal forecasting folds from 2018 and 2014, to test the model before
# the cross-validation, and as its forecasting folds,
#   v5f_full, v5f_fc2018, v5f_fc2014
# every chain from the low-floor initial values (each class's floor at a
# flat smooth starting at the cached floor, 0.0015; dynamical_inits()), 2,000
# warmup and 1,500 samples. On the default pod (cpu5c, 8 vCPU, $0.28/h),
# V5h's took 1.54 h (full) and 2.4 h (2018 fold, on cpu3c), with the image
# pull about 1.8-2.7 h and $0.50-0.75 a pod.
#
# The sets "cv5f" and "cv5f_warm" (cv5f_runs()): the cross-validation of
# V5f, the five folds of R/run_validation_folds.R,
#   cv5f_blocks1, cv5f_blocks2  spatial blocks 1 and 2
#   cv5f_interp                 spatial interpolation
#   v5f_fc2014, v5f_fc2018      forecasting from 2014 and from 2018: the
#                               test set's jobs, so not in "cv5f"
# from the low-floor initial values, as the test jobs; and as a fallback, if
# the test fold does not converge, all five as cv5f_warm_<fold>
# (blocks1, blocks2, interp, fc2014, fc2018), chain i from draw i of
# the full fit v5f_full (R/draw_inits.R: temporary/inits_v5f_draw1.RDS to
# draw4.RDS, one draw from each of its chains, spread along the pyrethroid
# floor at a flat smooth), so 4 chains. About 1-1.6 h a fold (V5h's 2014
# fold took 0.95 h), plus the pull: about $1.50 for the three pods of
# "cv5f", $2.50 for the five of "cv5f_warm".
#
# The sets "ts_cv5f" and "ts_cv5f_warm" (cv5f_two_stage()): the two-stage
# model on each fold of the set, with docker/two_stage_fold.sh, run by --in
# in the fold's job directory once the fold is done; and the set "maps_v5f"
# (v5f_maps): the two-stage maps of the full fit, with
# docker/two_stage_maps.sh, on three pods (maps_v5f_full_a, _b and _c), each
# with a third of the maps. Both scripts must be on the volume, in ts_v5f/,
# and take the commit the launcher is given (the second stage's, with the
# map draws of #43), not the folds'. On cpu3c pods of 16 vCPU and 32 GB
# ($0.48/h): about 20-60 minutes a fold (the October two-stage folds of PR
# #43 took 20-58), and about 2-2.5 h for each third of the maps
# (doc/v5f_cv_runbook.md).
#
# The options, the initial values and the job names are set in one place,
# below, each overridable from the environment, so that after another change
# to the model the sets are remade by this script, given the commit:
#   CV5F_OPTIONS  the options, in full (default v5f_options)
#   CV5F_NAME     the prefix of the fold jobs' names (default cv5f), which
#                 must be free on the volume
#   CV5F_FULL     the full fit: its job name in "v5f" (default v5f_full), and
#                 the fit "maps_v5f" maps and "cv5f_warm" draws from
#   CV5F_INITS    the prefix of the warm start's init files, as
#                 R/draw_inits.R writes them (default temporary/inits_v5f_draw)
v5f_options <- Sys.getenv("CV5F_OPTIONS", paste(
  "dynamical_model_options(mortality_floor = TRUE, floor_prior = c(1, 4),",
  "smooth = smooth_options(selection = TRUE, floor = \"class\",",
  "shear = FALSE, floor_intercepts = \"class\", kernel = \"se\", c = 2,",
  "range = 1.5, basis_range = 1.5,",
  "sd_prior = list(family = \"half_normal\", scale = 0.5),",
  "floor_intercept_prior = list(family = \"half_normal\", scale = 0.05)))"))
cv5f_name <- Sys.getenv("CV5F_NAME", "cv5f")
v5f_full_name <- Sys.getenv("CV5F_FULL", "v5f_full")
low_floor_inits <- "temporary/inits_floor_low.RDS"
cv5f_draw_inits <- paste(sprintf("%s%i.RDS",
                                 Sys.getenv("CV5F_INITS",
                                            "temporary/inits_v5f_draw"),
                                 1:4),
                         collapse = ",")
v5f_runs <- data.frame(
  name = c(v5f_full_name, "v5f_fc2018", "v5f_fc2014"),
  label = c("V5f_full", "V5f_fc2018", "V5f_fc2014"),
  options = v5f_options,
  inits = low_floor_inits,
  job = c("full", "fold temporal_forecasting 2018",
          "fold temporal_forecasting 2014"))
cv5f_folds <- data.frame(
  fold = c("blocks1", "blocks2", "interp", "fc2014", "fc2018"),
  job = c("fold spatial_blocks 1", "fold spatial_blocks 2",
          "fold spatial_interpolation all",
          "fold temporal_forecasting 2014", "fold temporal_forecasting 2018"))
# the forecasting folds of the cross-validation are the test set's jobs
# (v5f_fc2014 was launched with the test fits), so "cv5f" has the other
# three; the warm start, if needed, refits all five
cv5f_test_folds <- c(fc2014 = "v5f_fc2014", fc2018 = "v5f_fc2018")
cv5f_runs <- function(warm = FALSE) {
  prefix <- if (warm) paste0(cv5f_name, "_warm") else cv5f_name
  folds <- if (warm) cv5f_folds else
    cv5f_folds[!cv5f_folds$fold %in% names(cv5f_test_folds), ]
  data.frame(
    name = sprintf("%s_%s", prefix, folds$fold),
    label = sprintf("V5f%s_%s", if (warm) "_warm" else "", folds$fold),
    options = v5f_options,
    inits = if (warm) cv5f_draw_inits else low_floor_inits,
    job = folds$job,
    row.names = NULL)
}
# the fold job of each fold of the cross-validation: those of "cv5f" and the
# test set's forecasting folds, or those of "cv5f_warm"
cv5f_fold_jobs <- function(warm = FALSE) {
  if (warm) {
    return(setNames(cv5f_runs(warm = TRUE)$name, cv5f_folds$fold))
  }
  standard <- cv5f_runs()
  jobs <- setNames(standard$name,
                   cv5f_folds$fold[!cv5f_folds$fold %in% names(cv5f_test_folds)])
  c(jobs, cv5f_test_folds)[cv5f_folds$fold]
}
ts_script_dir <- "/workspace/ir_cube/ts_v5f"
cv5f_two_stage <- function(warm = FALSE) {
  jobs <- cv5f_fold_jobs(warm)
  data.frame(name = paste0("ts_", jobs),
             label = paste0("two_stage_V5f", if (warm) "_warm" else "", "_",
                            cv5f_folds$fold),
             options = "", inits = "",
             job = sprintf("--in %s -- bash %s/two_stage_fold.sh %s {ref}",
                           jobs, ts_script_dir,
                           sub("^fold ", "", cv5f_folds$job)),
             cpu = "cpu3c", vcpu = 16L, row.names = NULL)
}
# the maps on three pods, each with the map outputs of one part (TS_OUTPUTS
# of docker/two_stage_maps.sh) in a work directory of its own (TS_PART,
# maps_<full fit>_<part>/): the three LLIN types, mapped together as
# llin_effective, and the six other types in two sets of three, balanced by
# the times of their fits and maps on the Linux box in October (about 1.6 h
# each, and llin_effective's 1.7 h)
v5f_map_parts <- c(a = "llin_effective",
                   b = "Bendiocarb,Lambda-cyhalothrin,Fenitrothion",
                   c = "DDT,Pirimiphos-methyl,Malathion")
v5f_maps <- data.frame(
  name = sprintf("maps_%s_%s", v5f_full_name, names(v5f_map_parts)),
  label = sprintf("two_stage_maps_%s", names(v5f_map_parts)),
  options = "", inits = "",
  job = sprintf(paste("--in %s -- env TS_PART=%s TS_OUTPUTS=%s",
                      "bash %s/two_stage_maps.sh {ref}"),
                v5f_full_name, names(v5f_map_parts), v5f_map_parts,
                ts_script_dir),
  cpu = "cpu3c", vcpu = 16L, row.names = NULL)

# the sets of jobs, by name
run_sets <- list(v5f = v5f_runs,
                 cv5f = cv5f_runs(), cv5f_warm = cv5f_runs(warm = TRUE),
                 ts_cv5f = cv5f_two_stage(),
                 ts_cv5f_warm = cv5f_two_stage(warm = TRUE),
                 maps_v5f = v5f_maps)

# One row per pod job, with its environment (docker/README.md): the full fit,
# or the fold in the runs' column job, if they have one, with --threads; or a
# command run by --in, with the code ref in place of {ref}; and the pod's CPU
# type and vCPUs (columns cpu and vcpu of the runs, or cpu5c and 8)
species_run_jobs <- function(code_ref, threads = 8, runs = v5f_runs) {
  jobs <- runs
  if (is.null(jobs$job)) jobs$job <- "full"
  if (is.null(jobs$cpu)) jobs$cpu <- "cpu5c"
  if (is.null(jobs$vcpu)) jobs$vcpu <- 8L
  inside <- startsWith(jobs$job, "--in ")
  jobs$JOB <- ifelse(inside,
                     paste(jobs$name, gsub("{ref}", code_ref, jobs$job,
                                           fixed = TRUE)),
                     sprintf("%s %s --threads %d", jobs$name, jobs$job,
                             threads))
  jobs$CODE_REF <- code_ref
  jobs$IR_CUBE_MODEL_OPTIONS <- jobs$options
  jobs$IR_CUBE_INITS <- jobs$inits
  jobs[, c("name", "label", "JOB", "CODE_REF", "IR_CUBE_MODEL_OPTIONS",
           "IR_CUBE_INITS", "cpu", "vcpu")]
}

# The create-pod body of a job (the body of the RunPod connector's create-pod,
# or POST /v2/pods): the image, the job's pod (by default cpu5c, 8 vCPU, 16
# GB) in EU-RO-1 with the network volume at /workspace, and the job's
# environment (docker/README.md), with the code from the fork, which has the
# commits of this branch; IR_CUBE_MODEL_OPTIONS and IR_CUBE_INITS only when
# they are set (an --in job takes its fit's options, and stops if they are
# set to anything else)
species_run_pod_body <- function(job, vcpu = job$vcpu) {
  env <- list(JOB = job$JOB, CODE_REF = job$CODE_REF,
              CODE_REPO = "goldingn/ir_cube")
  if (nzchar(job$IR_CUBE_MODEL_OPTIONS)) {
    env$IR_CUBE_MODEL_OPTIONS <- job$IR_CUBE_MODEL_OPTIONS
  }
  if (nzchar(job$IR_CUBE_INITS)) {
    env$IR_CUBE_INITS <- job$IR_CUBE_INITS
  }
  list(name = gsub("_", "-", job$name),
       image = "ghcr.io/goldingn/ir_cube-runpod:latest",
       cpu = list(id = job$cpu, vcpuCount = vcpu),
       cloud = "SECURE",
       dataCenterIds = list("EU-RO-1"),
       mounts = list(network = list(list(volumeId = "0f6xjzaxch",
                                         path = "/workspace"))),
       disk = 20,
       env = env)
}

if (sys.nframe() == 0) {
  arguments <- commandArgs(trailingOnly = TRUE)
  code_ref <- if (length(arguments) >= 1) arguments[1] else
    system("git rev-parse HEAD", intern = TRUE)
  stopifnot(grepl("^[0-9a-f]{40}$", code_ref))
  set <- if (length(arguments) >= 2) arguments[2] else "v5f"
  stopifnot(set %in% names(run_sets))
  runs <- run_sets[[set]]
  suffix <- paste0("_", set)
  # each fit's options must pass the model's checks (definitions only; no
  # data or python)
  source("R/dynamical_model.R")
  for (expression in runs$options[nzchar(runs$options)]) {
    check_dynamical_model_options(eval(str2lang(expression)))
  }
  # the test fits and cross-validation of the default model have the
  # defaults in full
  if (set %in% c("v5f", "cv5f", "cv5f_warm")) {
    stopifnot(vapply(runs$options, function(expression) {
      identical(eval(str2lang(expression)), dynamical_model_options())
    }, logical(1)))
  }
  inputs <- c(character(0),
              unlist(strsplit(runs$inits[nzchar(runs$inits)], ",")))
  missing <- inputs[!file.exists(inputs)]
  if (length(missing) > 0) {
    warning("not found here, so not uploaded by irpod sync: ",
            toString(unique(missing)))
  }
  jobs <- species_run_jobs(code_ref, runs = runs)
  dir.create("outputs/species_runs", showWarnings = FALSE, recursive = TRUE)
  write.csv(jobs, sprintf("outputs/species_runs/jobs%s.csv", suffix),
            row.names = FALSE)
  bodies <- lapply(seq_len(nrow(jobs)), function(i) {
    species_run_pod_body(jobs[i, ])
  })
  jsonlite::write_json(bodies, sprintf("outputs/species_runs/pods%s.json",
                                       suffix),
                       auto_unbox = TRUE, pretty = TRUE)
  cat(sprintf("%d fits at %s\n", nrow(jobs), code_ref))
  print(jobs[, c("name", "label", "JOB", "IR_CUBE_MODEL_OPTIONS",
                 "IR_CUBE_INITS")], row.names = FALSE)
}
