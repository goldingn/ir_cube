#!/usr/bin/env bash
# The two-stage maps (#21, R/two_stage_maps.R) on a full fit, with the code
# of a given commit, in a work directory of its own.
#
#   run_pod_job.sh maps_<job> --in <job> -- \
#     bash /workspace/ir_cube/ts_v5f/two_stage_maps.sh <commit>
#
# Run in the full fit's job directory, as --in does (which writes only the
# step's markers and log there). Works in /workspace/ir_cube/ts_v5f/maps_<job>/
# (TS_ROOT/maps_<job>), with R/ and tmb/ of the commit of TS_CODE_REPO
# (default goldingn/ir_cube), or copied from TS_CODE_DIR; the volume's data,
# the job's fit (temporary/fitted_model.RData, linked) and its initial
# values, ingredient weights and rho tables. Steps: prepare; fit 1-9, in two
# queues of one process each, the types in alternating order of size; map of
# the six types not in LLINs and of llin_effective, one at a time, each with
# IR_CUBE_MAP_WORKERS forked workers (default 4; about 5 GB each); figures.
# TS_STEPS (a subset of "prepare fit map figures") runs only those steps, to
# resume. The outputs stay in the work directory (outputs/two_stage/,
# figures/two_stage/): fetch them with the aws command of
# doc/v5f_cv_runbook.md. On a 16 vCPU, 32 GB pod (cpu3c): prepare minutes,
# the fits about 1.5 h, the maps about 4 h.
set -euo pipefail
export IR_CUBE_NO_MEMORY_WAIT=1
ref=$1
job=$(pwd); name=$(basename "$job")
work=${TS_ROOT:-/workspace/ir_cube/ts_v5f}/maps_$name
threads=${RUNPOD_CPU_COUNT:-$(OMP_NUM_THREADS= nproc)}
half=$(( threads / 2 > 0 ? threads / 2 : 1 ))
workers=${IR_CUBE_MAP_WORKERS:-4}
steps=${TS_STEPS:-prepare fit map figures}
say() { echo "[$(date -u '+%F %T')] $*"; }
has() { [[ " $steps " == *" $1 "* ]]; }

say "set up $work from $job"
mkdir -p "$work"; cd "$work"
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

if has fit; then
  # the types by number of bioassays, largest first, dealt to the two queues
  # in the order a, b, b, a, a, b, b, a, a
  read -r queue_a queue_b < <(Rscript -e '
    files <- list.files("outputs/two_stage/maps", "dynamical.rds",
                        recursive = TRUE, full.names = TRUE)
    d <- readRDS(files[1])
    types <- d$parameters$types
    n <- vapply(types, function(type) {
      nrow(readRDS(file.path("outputs/two_stage/maps", type,
                             "dynamical.rds"))$keys)
    }, integer(1))
    order_k <- order(-n)
    a <- rep_len(c(TRUE, FALSE, FALSE, TRUE), length(order_k))
    cat(paste(order_k[a], collapse = ","), paste(order_k[!a], collapse = ","),
        "\n")
  ')
  say "fit queue a: $queue_a; queue b: $queue_b"
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
  # the six types not in LLINs (temporary/ingredient_weights.RDS), then the
  # three LLIN types together
  outputs=$(Rscript -e '
    files <- list.files("outputs/two_stage/maps", "dynamical.rds",
                        recursive = TRUE, full.names = TRUE)
    types <- readRDS(files[1])$parameters$types
    llin <- names(unlist(readRDS("temporary/ingredient_weights.RDS")))
    cat(which(!types %in% llin), "llin_effective")
  ')
  for o in $outputs; do
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
