#!/usr/bin/env bash
# The two-stage maps (#21, R/two_stage_maps.R) on a full fit, with the code
# of a given commit, in a work directory of its own: all of them on one pod,
#
#   run_pod_job.sh maps_<job> --in <job> -- \
#     bash /workspace/ir_cube/ts_v5f/two_stage_maps.sh <commit>
#
# or some of them, so that pods can share the maps (the three pods of
# R/species_runs.R's set "maps_v5f"):
#
#   run_pod_job.sh maps_<job>_<part> --in <job> -- env TS_PART=<part> \
#     TS_OUTPUTS=<outputs> bash /workspace/ir_cube/ts_v5f/two_stage_maps.sh <commit>
#
# Run in the full fit's job directory, as --in does (which writes only the
# step's markers and log there, named after the --in job). Works in
# /workspace/ir_cube/ts_v5f/maps_<job>/ (TS_ROOT/maps_<job>), or in
# maps_<job>_<part>/ with TS_PART, with R/ and tmb/ of the commit of
# TS_CODE_REPO (default goldingn/ir_cube), or copied from TS_CODE_DIR; the
# volume's data, the job's fit (temporary/fitted_model.RData, linked) and its
# initial values, ingredient weights and rho tables. Steps: prepare; fit, in
# two queues of one process each, the types in alternating order of size;
# map, one output at a time, each with IR_CUBE_MAP_WORKERS forked workers
# (default 4; about 5 GB each); figures.
#
# The map outputs are the six types not in LLINs and llin_effective, which
# maps Alpha-cypermethrin, Deltamethrin and Permethrin together, weighted by
# temporary/ingredient_weights.RDS. TS_OUTPUTS (comma-separated: types not in
# LLINs, by name or by index as R/two_stage_maps.R takes them, and
# llin_effective) maps only those, and fits only their types (llin_effective's
# three for it). It needs TS_PART, so that each pod has a work directory of
# its own; the work directory records its outputs (ts_outputs.txt), which a
# later run in it must repeat. Such a pod keeps only its types' dynamical
# draws (outputs/two_stage/maps/<type>/dynamical.rds, which prepare writes
# for every type), so that the pods' work directories merge without overlap,
# and does not make the figures, which need every type: they are made
# locally from the merged directories (doc/v5f_cv_runbook.md, step 3).
# TS_STEPS (a subset of "prepare fit map figures", or with TS_OUTPUTS of
# "prepare fit map") runs only those steps, to resume.
#
# The outputs stay in the work directory (outputs/two_stage/, figures/
# two_stage/): fetch them with the aws command of doc/v5f_cv_runbook.md. On
# a 16 vCPU, 32 GB pod (cpu3c), all on one pod: prepare minutes, the fits
# about 1.3 h, the maps about 4 h.
set -euo pipefail
export IR_CUBE_NO_MEMORY_WAIT=1
ref=$1
job=$(pwd); name=$(basename "$job")
part=${TS_PART:-}
outputs_wanted=${TS_OUTPUTS:-}
[ -z "$part" ] || [[ $part =~ ^[A-Za-z0-9_-]+$ ]] || {
  echo "TS_PART must be letters, digits, - or _" >&2; exit 1; }
[ -z "$outputs_wanted" ] || [ -n "$part" ] || {
  echo "TS_OUTPUTS needs TS_PART, for a work directory of its own" >&2; exit 1; }
work=${TS_ROOT:-/workspace/ir_cube/ts_v5f}/maps_$name${part:+_$part}
threads=${RUNPOD_CPU_COUNT:-$(OMP_NUM_THREADS= nproc)}
half=$(( threads / 2 > 0 ? threads / 2 : 1 ))
workers=${IR_CUBE_MAP_WORKERS:-4}
if [ -n "$outputs_wanted" ]; then
  steps=${TS_STEPS:-prepare fit map}
else
  steps=${TS_STEPS:-prepare fit map figures}
fi
say() { echo "[$(date -u '+%F %T')] $*"; }
has() { [[ " $steps " == *" $1 "* ]]; }
if [ -n "$outputs_wanted" ] && has figures; then
  echo "the figures need every type: make them from the merged work directories" >&2
  exit 1
fi

say "set up $work from $job"
mkdir -p "$work"; cd "$work"
# the outputs this work directory is for: a run for others stops here
record=ts_outputs.txt
if [ -f $record ]; then
  [ "$(cat $record)" = "${outputs_wanted:-all}" ] || {
    echo "$work is for the outputs $(cat $record), not ${outputs_wanted:-all}" >&2
    exit 1; }
else
  printf '%s\n' "${outputs_wanted:-all}" > $record
fi
if [ ! -d R ]; then
  if [ -n "${TS_CODE_DIR:-}" ]; then
    cp -r "$TS_CODE_DIR/R" "$TS_CODE_DIR/tmb" .
  else
    [[ $ref =~ ^[0-9a-f]{40}$ ]] || { echo "the commit must be a full id" >&2; exit 1; }
    curl -fsSL "https://codeload.github.com/${TS_CODE_REPO:-goldingn/ir_cube}/tar.gz/$ref" \
      | tar -xz --strip-components=1 --wildcards '*/R/' '*/tmb/'
  fi
fi
mkdir -p outputs/two_stage/maps temporary figures
ln -sfn "$(readlink -f "$job/data")" data
cp "$job"/temporary/*.RDS temporary/
ln -sfn "$job/temporary/fitted_model.RData" temporary/fitted_model.RData
cp "$job"/outputs/bioassay_rho*.csv outputs/

say "compile the correction template"
Rscript -e 'source("R/two_stage_correction.R"); invisible(load_correction_template())'

if has prepare; then
  say "prepare"
  OPENBLAS_NUM_THREADS=$threads Rscript R/two_stage_maps.R prepare
fi

if has prepare || has fit || has map; then
  # the map outputs (type indices, and llin_effective), and the types they
  # need fitted, by number of bioassays, largest first, dealt to the two fit
  # queues in the order a, b, b, a, a, b, b, a, a. With TS_OUTPUTS, the
  # other types' dynamical draws are deleted
  selection=$(TS_OUTPUTS=$outputs_wanted Rscript -e '
    maps_dir <- "outputs/two_stage/maps"
    files <- list.files(maps_dir, "dynamical.rds", recursive = TRUE,
                        full.names = TRUE)
    types <- readRDS(files[1])$parameters$types
    llin <- names(unlist(readRDS("temporary/ingredient_weights.RDS")))
    stopifnot(all(llin %in% types))
    wanted <- strsplit(Sys.getenv("TS_OUTPUTS"), "[, ]+")[[1]]
    wanted <- wanted[nzchar(wanted)]
    split <- length(wanted) > 0
    if (!split) wanted <- c(which(!types %in% llin), "llin_effective")
    outputs <- vapply(wanted, function(o) {
      if (o == "llin_effective") return(o)
      k <- if (grepl("^[0-9]+$", o)) as.integer(o) else match(o, types)
      if (is.na(k) || k < 1 || k > length(types)) {
        stop("no insecticide type ", o)
      }
      if (types[k] %in% llin) stop(types[k], " is mapped by llin_effective")
      as.character(k)
    }, "", USE.NAMES = FALSE)
    stopifnot(!anyDuplicated(outputs))
    fit <- unique(unlist(lapply(outputs, function(o) {
      if (o == "llin_effective") match(llin, types) else as.integer(o)
    })))
    if (split) {
      drop <- file.path(maps_dir, types[-fit])
      unlink(file.path(drop, "dynamical.rds"))
      for (d in drop[dir.exists(drop)]) {
        if (length(list.files(d, all.files = TRUE, no.. = TRUE)) == 0) {
          unlink(d, recursive = TRUE)
        }
      }
    }
    n <- vapply(types[fit], function(type) {
      nrow(readRDS(file.path(maps_dir, type, "dynamical.rds"))$keys)
    }, integer(1))
    order_k <- fit[order(-n)]
    a <- rep_len(c(TRUE, FALSE, FALSE, TRUE), length(order_k))
    message("map: ", paste(ifelse(outputs == "llin_effective", outputs,
                                  types[suppressWarnings(as.integer(outputs))]),
                           collapse = ", "),
            "; fit: ", paste(types[order_k], collapse = ", "))
    writeLines(c(paste(outputs, collapse = " "),
                 paste(order_k[a], collapse = ","),
                 paste(order_k[!a], collapse = ",")))
  ')
  map_outputs=$(sed -n 1p <<< "$selection")
  queue_a=$(sed -n 2p <<< "$selection")
  queue_b=$(sed -n 3p <<< "$selection")
fi

if has fit; then
  say "fit queue a: $queue_a; queue b: ${queue_b:-none}"
  run_queue() {
    local failed=0
    for k in ${1//,/ }; do
      say "start fit $k"
      if OPENBLAS_NUM_THREADS=$half Rscript R/two_stage_maps.R fit "$k" \
           > "fit_$k.log" 2>&1; then say "done fit $k"; else say "FAILED fit $k"; failed=1; fi
    done
    return $failed
  }
  set +e
  run_queue "$queue_a" > fit_queue_a.log 2>&1 & pa=$!
  run_queue "$queue_b" > fit_queue_b.log 2>&1 & pb=$!
  wait $pa; sa=$?; wait $pb; sb=$?
  set -e
  cat fit_queue_a.log fit_queue_b.log
  for f in fit_*.log; do echo "=== $f"; tail -2 "$f"; done
  [ $sa -eq 0 ] && [ $sb -eq 0 ] || { say "a fit failed"; exit 1; }
fi

if has map; then
  for o in $map_outputs; do
    say "start map $o"
    IR_CUBE_MAP_WORKERS=$workers OPENBLAS_NUM_THREADS=1 \
      Rscript R/two_stage_maps.R map "$o" > "map_$o.log" 2>&1 \
      || { say "FAILED map $o"; tail -20 "map_$o.log"; exit 1; }
    say "done map $o: $(tail -1 "map_$o.log")"
  done
fi

if has figures; then
  say "figures"
  Rscript R/two_stage_maps.R figures
fi
say "done"
