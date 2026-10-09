#!/usr/bin/env bash
# The two-stage model (#21) on one dynamical cross-validation fold, with the
# code of a given commit, in a work directory of its own: the fold's code can
# be older than the second stage's (e.g. without the map draws of #43).
#
#   run_pod_job.sh ts_<job> --in <job> -- \
#     bash /workspace/ir_cube/ts_v5f/two_stage_fold.sh <experiment> <fold> <commit>
#
# e.g. for the fold job cv5f_blocks1 (R/species_runs.R):
#   JOB="ts_cv5f_blocks1 --in cv5f_blocks1 -- bash /workspace/ir_cube/ts_v5f/two_stage_fold.sh spatial_blocks 1 <commit>"
#
# Run in the fold's job directory, as --in does. Works in
# /workspace/ir_cube/ts_v5f/<job>/ (TS_ROOT/<job>), with R/ and tmb/ of the
# commit of TS_CODE_REPO (default goldingn/ir_cube), or copied from
# TS_CODE_DIR; the volume's data, and the job's dynamical fold, initial
# values and rho tables. The steps of R/run_two_stage_folds.R: stage_one
# (which checks that the rebuilt held-out set and the recomputed draws are
# the fold's), the types in two queues of one process each, then assemble.
# The results are copied into the job directory, so that `irpod fetch <job>`
# brings them: outputs/cv_draws/two_stage__<experiment>__<fold>.rds and
# outputs/two_stage/ (fit_summary.csv, parts/, logs). The stage-one cache
# stays in the work directory. TS_TYPES (comma-separated type indices) fits
# only those types, the others falling back to the dynamical draws: for smoke
# tests. On a 16 vCPU, 32 GB pod (cpu3c), about 1 h per fold; the largest
# types peak at about 14 GB, so the queues take the types in alternating
# order of training size.
set -euo pipefail
export IR_CUBE_NO_MEMORY_WAIT=1
experiment=$1; fold=$2; ref=$3
job=$(pwd); name=$(basename "$job"); key=${experiment}__${fold}
work=${TS_ROOT:-/workspace/ir_cube/ts_v5f}/$name
threads=${RUNPOD_CPU_COUNT:-$(OMP_NUM_THREADS= nproc)}
half=$(( threads / 2 > 0 ? threads / 2 : 1 ))
say() { echo "[$(date -u '+%F %T')] $*"; }

say "set up $work from $job"
mkdir -p "$work"; cd "$work"
if [ -n "${TS_CODE_DIR:-}" ]; then
  cp -r "$TS_CODE_DIR/R" "$TS_CODE_DIR/tmb" .
else
  [[ $ref =~ ^[0-9a-f]{40}$ ]] || { echo "the commit must be a full id" >&2; exit 1; }
  curl -fsSL "https://codeload.github.com/${TS_CODE_REPO:-goldingn/ir_cube}/tar.gz/$ref" \
    | tar -xz --strip-components=1 --wildcards '*/R/' '*/tmb/'
fi
mkdir -p outputs/cv_draws outputs/two_stage/parts temporary figures
ln -sfn "$(readlink -f "$job/data")" data
cp "$job"/temporary/*.RDS temporary/
cp "$job"/outputs/bioassay_rho*.csv outputs/
ln -sfn "$job/outputs/cv_draws/dynamical__$key.rds" \
  "outputs/cv_draws/dynamical__$key.rds"

say "compile the correction template"
Rscript -e 'source("R/two_stage_correction.R"); invisible(load_correction_template())'

say "stage_one"
OPENBLAS_NUM_THREADS=$threads Rscript R/run_two_stage_folds.R "$experiment" "$fold" stage_one

# the types by training size, largest first, dealt to the two queues in the
# order a, b, b, a, a, b, b, a, a
read -r queue_a queue_b < <(Rscript -e '
  cache <- readRDS(commandArgs(TRUE)[1])
  n <- tabulate(cache$training$type_id, length(cache$types))
  only <- Sys.getenv("TS_TYPES")
  types <- order(-n)
  if (nzchar(only)) types <- types[types %in% as.integer(strsplit(only, ",")[[1]])]
  a <- rep_len(c(TRUE, FALSE, FALSE, TRUE), length(types))
  cat(paste(types[a], collapse = ","), paste(types[!a], collapse = ","))
' "outputs/two_stage/stage_one__$key.rds")
say "queue a: ${queue_a:-none}; queue b: ${queue_b:-none}"
run_queue() {
  local failed=0
  for k in ${1//,/ }; do
    say "start $k"
    if OPENBLAS_NUM_THREADS=$half Rscript R/run_two_stage_folds.R "$experiment" "$fold" "$k" \
         > "type_$k.log" 2>&1; then say "done $k"; else say "FAILED $k"; failed=1; fi
  done
  return $failed
}
set +e
run_queue "${queue_a:-}" > queue_a.log 2>&1 & pa=$!
run_queue "${queue_b:-}" > queue_b.log 2>&1 & pb=$!
wait $pa; sa=$?; wait $pb; sb=$?
set -e
cat queue_a.log queue_b.log
for f in type_*.log; do [ -e "$f" ] || continue; echo "=== $f"; tail -3 "$f"; done
[ $sa -eq 0 ] && [ $sb -eq 0 ] || { say "a type failed"; exit 1; }

say "assemble"
Rscript R/run_two_stage_folds.R "$experiment" "$fold" assemble

say "copy the results to $job"
cp "outputs/cv_draws/two_stage__$key.rds" "$job/outputs/cv_draws/"
mkdir -p "$job/outputs/two_stage"
cp -r outputs/two_stage/fit_summary.csv outputs/two_stage/parts \
  queue_a.log queue_b.log "$job/outputs/two_stage/"
cp type_*.log "$job/outputs/two_stage/" 2>/dev/null || true
say "done"
