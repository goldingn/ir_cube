# The fits of the species model (#47) and its reference, full data only, of
# V5 with a fixed range, full and forecasting folds, and the test fits and
# cross-validation of V5f, the default model, with the second stage and the
# maps. Run on RunPod (docker/README.md), one fit per pod.
#
#   Rscript R/species_runs.R [code ref] [set]
#
# writes outputs/species_runs/jobs.csv, one row per fit, and
# outputs/species_runs/pods.json, the create-pod body of each
# (species_run_pod_body()), for the commit to run (a full commit
# id; default: HEAD, which must be on GitHub), for the set "species" (the
# default, species_runs); for the other sets of run_sets, below ("wb",
# "v5r", "v5i", "v5h", "v5f", "cv5f", "cv5f_warm", "ts_cv5f", "ts_cv5f_warm"
# and "maps_v5f"), jobs_<set>.csv and pods_<set>.json. The fits:
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
#   sp_v5         V5: no kdr covariate; latent smooths (smooth_options(),
#                 R/latent_smooth.R; squared exponential, P(range < 1,500
#                 km) = 0.05, P(sd > 1) = 0.05) of the strength of selection,
#                 for every class, and of the logit floor, for the
#                 pyrethroids and DDT, with an intercept per class, its prior
#                 matched to Beta(1, 4)
#   sp_v5_shear   V5 shear: as V5, with the selection smooth u_s = v_s + b
#                 u_f, sharing the floor's u_f through one loading b ~ N(1,
#                 0.5) (smooth_options(shear = TRUE)); to run with V5
#   sp_v5_sel     V5 sel: as V5, without the smooth of the floor (a constant
#                 floor per class)
#   sp_v5_floor   V5 floor: as V5, without the smooth of selection
# sp_v5_sel and sp_v5_floor are for comparisons after V5: do not launch them
# with it.
# Otherwise the default options of the time (#37: d_half 270, reversion
# estimated, no floor, no smooths, the beta-binomial likelihood), which each
# option string sets where it does not set them itself (pre_v5h_defaults), so
# that a rerun under the later defaults (V5h, V5f; #47) is the same model;
# V5's options give the logit-normal prior of its floor intercepts
# (floor_intercept_prior = "beta_moments"), as do V5r's, V5i's and V5h's.
# A rerun of a model with the floor smooth centres it over the cells with
# bioassays of the pyrethroids or DDT, not over every modelled cell as these
# fits did (their saved bases rebuild them as they were). The fits
# with floors start chains 1-2 from the low-floor mode and 3-4 from the
# high-floor mode (IR_CUBE_INITS, dynamical_chain_inits(); with the species
# model, both species' floors start at the cached floor, and with the kdr or
# smooth floor, its intercepts; dynamical_inits()), and record each chain's
# floors and log posterior (chain_floor_modes(); with the smooths, their sd
# and range too); V1 starts from the default initial values. Before creating
# the pods:
#   Rscript R/floor_mode_inits.R <r2_main.RData>  # temporary/inits_floor_*.RDS
#   irpod sync    # uploads them, and data/clean/arabiensis_fraction.tif
# Then create one pod per element of the json (about 3.7 h and $1.05 per fit
# on the default pod, by the gradient timing, which is the same with and
# without the species model; the full fit of 5 October took 6.2 h. The fits
# with floors took 6.3-6.9 h; V5's gradient takes about 41 ms, and V5
# shear's 42, against V4_class's 38, at 4 chains and 4 threads, so about
# 6.7-7.5 h and $2 each; state the price first), and when each is done,
# `irpod fetch <name>` and delete its pod (it does not delete itself).
#
# The set "wb" (wb_runs): the candidate models with the weighted binomial
# likelihood (#47; dynamical_model_options(likelihood = "weighted_binomial"),
# R/weighted_binomial.R), rho per type fixed at the replicate estimates
# (data/clean/bioassay_rho_replicate.csv, which irpod sync uploads), with the
# default sampler of #48 (centred data-informed hierarchy, 30-60 leapfrog
# steps, 2,000 warmup and 1,500 samples):
#   wb_ref        ref_f0: the default options, no floor
#   wb_bf         B_f, as sp_ref_floor
#   wb_v3f        V3f, as sp_v3_floor
#   wb_v4         V4, as sp_v4
#   wb_v4_class   V4_class, as sp_v4_class
#   wb_v5         V5, as sp_v5
# each with the options and initial values of the fit it copies (ref_f0 from
# the default initial values), and likelihood = "weighted_binomial". Their
# gradients took 0.87-0.95 times the default beta-binomial model's (4 chains,
# 8 threads, 8 October 2026), which took 48 ms on the default pod: about
# 42-46 ms there, and for 3,500 iterations of 45 leapfrog steps, 1.8-2.0 h
# each, about 2.5 h and $0.70 a pod with setup and saving.

# The options the species and wb fits left at their defaults, at the
# defaults of the time (#37), before those of V5h (#47)
pre_v5h_defaults <- list(mortality_floor = FALSE, floor_prior = c(1, 49),
                         smooth = FALSE, likelihood = "beta_binomial")

# An option string (a call of dynamical_model_options()) with the options in
# `settings` (a named list of values) added where it does not set them, or
# with replace = TRUE, set to them whether it does or not
set_options <- function(expression, settings, replace = FALSE) {
  call <- str2lang(expression)
  stopifnot(identical(call[[1]], quote(dynamical_model_options)))
  for (name in names(settings)) {
    if (replace || is.null(call[[name]])) call[[name]] <- settings[[name]]
  }
  paste(deparse(call, width.cutoff = 500L), collapse = " ")
}

floor_mode_inits <- paste("temporary/inits_floor_low.RDS",
                          "temporary/inits_floor_high.RDS", sep = ",")
species_floors <- "species_options(floors = TRUE, floor_prior = c(1, 4))"
# the latent smooths of V5, in full, with the smooths on selection and the
# floor, and the shear, as given; the range estimated
v5_smooth <- function(selection, floor, shear = "FALSE") {
  sprintf(paste("mortality_floor = TRUE, floor_prior = c(1, 4),",
                "smooth = smooth_options(selection = %s, floor = %s,",
                "shear = %s, floor_intercepts = \"class\", kernel = \"se\",",
                "c = 2, range = NULL, basis_range = 1,",
                "range_prior = c(1.5, 0.05), sd_prior = c(1, 0.05),",
                "floor_intercept_prior = \"beta_moments\")"),
          selection, floor, shear)
}
species_runs <- data.frame(
  name = c("sp_ref_floor", "sp_v1", "sp_v1_floor", "sp_v2", "sp_v2_floor",
           "sp_v3", "sp_v3_floor", "sp_v4", "sp_v4_class", "sp_v5",
           "sp_v5_shear", "sp_v5_sel", "sp_v5_floor"),
  label = c("B_f", "V1", "V1f", "V2", "V2f", "V3", "V3f", "V4", "V4_class",
            "V5", "V5_shear", "V5_sel", "V5_floor"),
  options = vapply(sprintf("dynamical_model_options(%s)", c(
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
          "kdr = kdr_options(floor = \"class\")"),
    v5_smooth("TRUE", "\"class\""),
    v5_smooth("TRUE", "\"class\"", shear = "TRUE"),
    v5_smooth("TRUE", "FALSE"),
    v5_smooth("FALSE", "\"class\""))),
    set_options, "", settings = pre_v5h_defaults, USE.NAMES = FALSE),
  inits = c(floor_mode_inits, "", floor_mode_inits, "", floor_mode_inits, "",
            floor_mode_inits, floor_mode_inits, floor_mode_inits,
            floor_mode_inits, floor_mode_inits, floor_mode_inits,
            floor_mode_inits))

# The candidate models with the weighted binomial likelihood (the set "wb"):
# the options of each fit it copies, with likelihood = "weighted_binomial"
weighted_binomial_options <- function(expression) {
  set_options(expression, list(likelihood = "weighted_binomial"),
              replace = TRUE)
}
wb_copies <- c(wb_bf = "B_f", wb_v3f = "V3f", wb_v4 = "V4",
               wb_v4_class = "V4_class", wb_v5 = "V5")
wb_runs <- rbind(
  data.frame(name = "wb_ref", label = "ref_f0_wb",
             options = set_options(paste0("dynamical_model_options(",
                                          "likelihood = \"weighted_binomial\")"),
                                   pre_v5h_defaults),
             inits = ""),
  data.frame(name = names(wb_copies),
             label = paste0(wb_copies, "_wb"),
             options = vapply(species_runs$options[match(wb_copies,
                                                          species_runs$label)],
                              weighted_binomial_options, ""),
             inits = species_runs$inits[match(wb_copies, species_runs$label)],
             row.names = NULL))

# The set "v5r" (v5r_runs): V5 with the range of both smooths fixed at 1,500
# km (smooth_options(range = 1.5), R/latent_smooth.R; the basis set for it,
# 164 functions), each smooth's sd estimated with its PC prior, with the
# beta-binomial (bb) and the weighted binomial (wb) likelihood, the full fit
# and the temporal forecasting folds from 2014 and 2018 (fc2014, fc2018;
# R/run_one_fold.R):
#   v5r_bb_full, v5r_wb_full, v5r_bb_fc2014, v5r_wb_fc2014, v5r_bb_fc2018,
#   v5r_wb_fc2018
# The range is fixed because V5's fitted ranges (500-650 km) piled against
# the basis's limit; 1,500 km is the spatial correlation range of recent
# pyrethroid mortality (1,230 km [840, 1,800]) and of block plateaus (1,500
# km [1,100, 2,000]; R/plateau_range.R). Every chain starts from the
# low-floor initial values: in the V5 fits, the chains from the high-floor
# ones found a mode 125 log-posterior units lower. With the default sampler
# of #48. V5r's gradient took 0.98 (bb) and 0.94 (wb) times the default
# model's (4 chains, 8 threads, 9 October 2026; as wb_v5's, within the
# timing noise): about 45-47 ms on the default pod, and for 3,500 iterations
# of 45 leapfrog steps about 2.0 h (wb_ref, at 0.90, sampled in 1.9 h), so
# about 2.3-2.5 h and $0.70 a pod with setup, predictions and saving.
v5r_options <- function(likelihood) {
  sprintf(paste("dynamical_model_options(mortality_floor = TRUE,",
                "floor_prior = c(1, 4),",
                "smooth = smooth_options(selection = TRUE, floor = \"class\",",
                "shear = FALSE, floor_intercepts = \"class\", kernel = \"se\",",
                "c = 2, range = 1.5, basis_range = 1.5,",
                "floor_intercept_prior = \"beta_moments\",",
                "sd_prior = c(1, 0.05)), likelihood = \"%s\")"),
          likelihood)
}
v5r_fits <- data.frame(fit = c("full", "fc2014", "fc2018"),
                       job = c("full", "fold temporal_forecasting 2014",
                               "fold temporal_forecasting 2018"))
v5r_likelihoods <- c(bb = "beta_binomial", wb = "weighted_binomial")
v5r_runs <- merge(data.frame(short = names(v5r_likelihoods),
                             likelihood = unname(v5r_likelihoods)),
                  v5r_fits, by = NULL)
v5r_runs <- v5r_runs[order(match(v5r_runs$fit, v5r_fits$fit),
                           match(v5r_runs$short, names(v5r_likelihoods))), ]
v5r_runs <- data.frame(
  name = sprintf("v5r_%s_%s", v5r_runs$short, v5r_runs$fit),
  label = sprintf("V5r_%s_%s", v5r_runs$short, v5r_runs$fit),
  options = v5r_options(v5r_runs$likelihood),
  inits = "temporary/inits_floor_low.RDS",
  job = v5r_runs$job,
  row.names = NULL)

# The set "v5i" (v5i_runs): V5r with the beta-binomial likelihood and the
# hierarchy of regions and countries of the initial state replaced by a
# third smooth, of the logit relative initial state, with a loading per type
# (smooth_options(init = TRUE), R/latent_smooth.R; the same kernel, basis
# and fixed range, the field's sd fixed at 1 and each loading with the PC
# prior of a smooth's sd), within the same limits init_frac_min and with the
# initial-state covariates; the full fit and the temporal forecasting folds
# from 2014 and 2018:
#   v5i_bb_full, v5i_bb_fc2014, v5i_bb_fc2018
# Every chain starts from the low-floor initial values, as V5r's, the field
# flat and each loading at 0.3 (dynamical_inits()). With the default sampler
# of #48. V5i's gradient took 1.00 times V5r's (bb) and 1.02 times the
# default model's (4 chains, 8 threads, 9 October 2026): about 49 ms on the
# default pod, and for 3,500 iterations of 45 leapfrog steps about 2.1 h, so
# about 2.4-2.6 h and $0.70 a pod with setup, predictions and saving.
v5i_options <- sub("sd_prior = c(1, 0.05))",
                   "sd_prior = c(1, 0.05), init = TRUE)",
                   v5r_options("beta_binomial"), fixed = TRUE)
v5i_runs <- data.frame(
  name = sprintf("v5i_bb_%s", v5r_fits$fit),
  label = sprintf("V5i_bb_%s", v5r_fits$fit),
  options = v5i_options,
  inits = "temporary/inits_floor_low.RDS",
  job = v5r_fits$job,
  row.names = NULL)

# The set "v5h" (v5h_runs): V5r with the beta-binomial likelihood and each
# smooth's sd with the half-normal prior of scale 0.5 in place of the PC
# prior (smooth_options(sd_prior = list(family = "half_normal", scale =
# 0.5)), R/latent_smooth.R), under which V5r's sds ran to 4-13; the full fit
# and the temporal forecasting folds from 2014 and 2018:
#   v5h_bb_full, v5h_bb_fc2014, v5h_bb_fc2018
# Every chain starts from the low-floor initial values, as V5r's. With the
# default sampler of #48. The model is V5r's but for the prior of two
# scalars, so about V5r's time and cost, 2.3-2.5 h and $0.70 a pod.
v5h_options <- sub("sd_prior = c(1, 0.05)",
                   "sd_prior = list(family = \"half_normal\", scale = 0.5)",
                   v5r_options("beta_binomial"), fixed = TRUE)
v5h_runs <- data.frame(
  name = sprintf("v5h_bb_%s", v5r_fits$fit),
  label = sprintf("V5h_bb_%s", v5r_fits$fit),
  options = v5h_options,
  inits = "temporary/inits_floor_low.RDS",
  job = v5r_fits$job,
  row.names = NULL)

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
  "floor_intercept_prior = list(family = \"half_normal\", scale = 0.05)),",
  "likelihood = \"beta_binomial\")"))
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

# the sets of fits, by name
run_sets <- list(species = species_runs, wb = wb_runs, v5r = v5r_runs,
                 v5i = v5i_runs, v5h = v5h_runs, v5f = v5f_runs,
                 cv5f = cv5f_runs(), cv5f_warm = cv5f_runs(warm = TRUE),
                 ts_cv5f = cv5f_two_stage(),
                 ts_cv5f_warm = cv5f_two_stage(warm = TRUE),
                 maps_v5f = v5f_maps)

# One row per pod job, with its environment (docker/README.md): the full fit,
# or the fold in the runs' column job, if they have one, with --threads; or a
# command run by --in, with the code ref in place of {ref}; and the pod's CPU
# type and vCPUs (columns cpu and vcpu of the runs, or cpu5c and 8)
species_run_jobs <- function(code_ref, threads = 8, runs = species_runs) {
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
  set <- if (length(arguments) >= 2) arguments[2] else "species"
  stopifnot(set %in% names(run_sets))
  runs <- run_sets[[set]]
  suffix <- if (set == "species") "" else paste0("_", set)
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
  inputs <- c(arabiensis_fraction_file, kdr_total_file,
              if (set %in% c("wb", "v5r")) replicate_rho_file,
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
