# The grid of #37: the population transform's d_half x the mortality floor,
# each fitted to the full data and to every cross-validation fold, with linear
# net use. Run on RunPod (docker/README.md), one fit per pod.
#
#   Rscript R/net_grid.R [set] [code ref]
#
# writes outputs/net_grid/jobs_<set>.csv, one row per fit, and
# outputs/net_grid/pods_<set>.json, the create-pod body of each (the body of
# the RunPod connector's create-pod, or POST /v2/pods), for the commit to run
# (a full commit id; default: HEAD, which must be on GitHub). The sets
# (net_grid_sets; default "preferred"):
#   preferred  d_half 270, the floor fixed at 0: the default options (the
#              full fit run on 5 October as dh270_lin_f0_full). 6 fits
#   issue37    d_half 270 and 50, floor estimated, chains started in both
#              floor modes. 12
#   floor0     d_half 270 and 200, floor fixed at 0. 12
#   all        every variant above. 24
# Before creating pods for a set with the floor estimated:
#   Rscript R/floor_mode_inits.R    # temporary/inits_floor_{low,high}.RDS
#   irpod sync                      # uploads them with the other inputs
# Then create one pod per element of the json (about 3.7 h and $1.05 per fit
# on the default pod; the full fit of 5 October took 6.2 h; state the price
# first), and when each is done, `irpod fetch <name>` and delete its pod (it
# does not delete itself). R/net_grid_report.R summarises the fetched fits.
#
# The variants:
#   d_half       270, where the population covariate is most spread across
#                the bioassays (sd 0.342; R/pop_d_half_spread.R); 200 (sd
#                0.339); 50, the value before the refit
#   floor        "est", the estimated floor, with chains 1-2 started from the
#                low-floor mode and 3-4 from the high-floor mode
#                (IR_CUBE_INITS, dynamical_chain_inits()); or "f0", no floor,
#                from the default initial values
#   fits         the full data, and the folds of R/run_validation_folds.R:
#                spatial blocks 1 and 2, interpolation, forecasting from 2014
#                and 2018
# Each variant's options set d_half and mortality_floor explicitly, whatever
# the defaults.
# Chains: 4, the default 2,000 warmup + 3,000 samples. Each fit with a floor
# records each chain's mode and log posterior (chain_floor_modes(),
# R/chain_floor_modes.R). windowed_hmc() estimates its mass matrix from all
# chains pooled, so while chains sit in different modes it is wider than
# either mode's, and each mixes more slowly than a fit started in one mode.

net_grid_sets <- list(
  preferred = "dh270_f0",
  issue37 = c("dh270_est", "dh50_est"),
  floor0 = c("dh270_f0", "dh200_f0"))
net_grid_sets$all <- unique(unlist(net_grid_sets))

# The grid's variants: name (dh<d_half>_<est|f0>), d_half, floor, the R
# expression for dynamical_model_options() that IR_CUBE_MODEL_OPTIONS takes,
# and the initial values (IR_CUBE_INITS; empty for the default)
net_grid_variants <- function(set = "all") {
  grid <- expand.grid(d_half = c(270, 200, 50), floor = c("est", "f0"),
                      stringsAsFactors = FALSE)
  grid$variant <- sprintf("dh%g_%s", grid$d_half, grid$floor)
  grid$options <- vapply(seq_len(nrow(grid)), function(i) {
    sprintf(paste0("dynamical_model_options(selection_columns = ",
                   "selection_design(pop_d_half = %g), ",
                   "mortality_floor = %s)"),
            grid$d_half[i], grid$floor[i] == "est")
  }, character(1))
  grid$inits <- ifelse(grid$floor == "est",
                       paste("temporary/inits_floor_low.RDS",
                             "temporary/inits_floor_high.RDS", sep = ","),
                       "")
  stopifnot(set %in% names(net_grid_sets))
  grid <- grid[match(net_grid_sets[[set]], grid$variant), ]
  rownames(grid) <- NULL
  grid[, c("variant", "d_half", "floor", "options", "inits")]
}

# The fits of each variant: the full data, and the folds as run_one_fold.R
# names them
net_grid_fits <- data.frame(
  fit = c("full", "blocks1", "blocks2", "interp", "fc2014", "fc2018"),
  experiment = c(NA, "spatial_blocks", "spatial_blocks",
                 "spatial_interpolation", "temporal_forecasting",
                 "temporal_forecasting"),
  fold = c(NA, "1", "2", "all", "2014", "2018"))
net_grid_fits$job <- ifelse(
  is.na(net_grid_fits$experiment), "full",
  paste("fold", net_grid_fits$experiment, net_grid_fits$fold))

# One row per pod job: its name (grid_<variant>_<fit>), JOB, and the
# environment
net_grid_jobs <- function(code_ref, set = "all", threads = 8) {
  jobs <- merge(net_grid_variants(set), net_grid_fits, by = NULL)
  jobs$name <- sprintf("grid_%s_%s", jobs$variant, jobs$fit)
  jobs$JOB <- sprintf("%s %s --threads %d", jobs$name, jobs$job, threads)
  jobs$CODE_REF <- code_ref
  jobs$IR_CUBE_MODEL_OPTIONS <- jobs$options
  jobs$IR_CUBE_INITS <- jobs$inits
  jobs <- jobs[order(match(jobs$variant, net_grid_sets[[set]]),
                     match(jobs$fit, net_grid_fits$fit)), ]
  rownames(jobs) <- NULL
  jobs[, c("name", "variant", "fit", "experiment", "fold", "JOB", "CODE_REF",
           "IR_CUBE_MODEL_OPTIONS", "IR_CUBE_INITS")]
}

# The create-pod body of a job: the image, the default pod (cpu5c, 8 vCPU,
# 16 GB) in EU-RO-1 with the network volume at /workspace, and the job's
# environment (docker/README.md); IR_CUBE_INITS only when it is set
net_grid_pod_body <- function(job, vcpu = 8) {
  env <- list(JOB = job$JOB, CODE_REF = job$CODE_REF,
              IR_CUBE_MODEL_OPTIONS = job$IR_CUBE_MODEL_OPTIONS)
  if (nzchar(job$IR_CUBE_INITS)) {
    env$IR_CUBE_INITS <- job$IR_CUBE_INITS
  }
  list(name = gsub("_", "-", job$name),
       image = "ghcr.io/goldingn/ir_cube-runpod:latest",
       cpu = list(id = "cpu5c", vcpuCount = vcpu),
       cloud = "SECURE",
       dataCenterIds = list("EU-RO-1"),
       mounts = list(network = list(list(volumeId = "0f6xjzaxch",
                                         path = "/workspace"))),
       disk = 20,
       env = env)
}

if (sys.nframe() == 0) {
  arguments <- commandArgs(trailingOnly = TRUE)
  set <- if (length(arguments) >= 1) arguments[1] else "preferred"
  code_ref <- if (length(arguments) >= 2) arguments[2] else
    system("git rev-parse HEAD", intern = TRUE)
  stopifnot(set %in% names(net_grid_sets),
            grepl("^[0-9a-f]{40}$", code_ref))
  # each variant's selection design must build (R/model_covariates.R),
  # without data
  source("R/model_covariates.R")
  dynamical_model_options <- function(...) list(...)
  for (expression in net_grid_variants(set)$options) {
    complete_selection_design(eval(str2lang(expression))$selection_columns)
  }
  rm(dynamical_model_options)
  jobs <- net_grid_jobs(code_ref, set)
  dir.create("outputs/net_grid", showWarnings = FALSE, recursive = TRUE)
  write.csv(jobs, sprintf("outputs/net_grid/jobs_%s.csv", set),
            row.names = FALSE)
  bodies <- lapply(seq_len(nrow(jobs)), function(i) {
    net_grid_pod_body(jobs[i, ])
  })
  jsonlite::write_json(bodies, sprintf("outputs/net_grid/pods_%s.json", set),
                       auto_unbox = TRUE, pretty = TRUE)
  cat(sprintf("set %s: %d fits (%d variants x %d) at %s\n", set, nrow(jobs),
              length(net_grid_sets[[set]]), nrow(net_grid_fits), code_ref))
  print(jobs[, c("name", "JOB")], row.names = FALSE)
  cat("\n")
  print(unique(jobs[, c("variant", "IR_CUBE_MODEL_OPTIONS", "IR_CUBE_INITS")]),
        row.names = FALSE)
}
