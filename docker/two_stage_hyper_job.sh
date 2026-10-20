#!/usr/bin/env bash
# One #34 experiment (R/two_stage_hyper_folds.R) on a pod: a fold, a type and
# its methods, with the code of a given commit, in a work directory of its own.
#
#   JOB="<name> --in <fold job> -- bash /workspace/ir_cube/ts_hyper/hyper_job.sh \
#        <experiment> <fold> <type index> <methods> <commit>"
#
# hyper_job.sh is a copy of this script on the volume. <fold job> is the
# dynamical fold's job (e.g. cv5pg_fc2014), whose second stage
# (ts_v5f/<fold job>/) holds the stage-one cache. Works in
# /workspace/ir_cube/ts_hyper/<name>/. HYPER_FITS_FROM=<earlier name> starts
# from that job's saved fits (its outputs/two_stage/hyper/fits/), and
# STAGE_ONE_DIR=<dir> takes the stage-one cache from <dir> (e.g. a copy with
# the study labels of #54). The *_is
# methods evaluate their draws in one single-threaded process per vCPU; the
# other methods use all the vCPUs.
set -euo pipefail
export IR_CUBE_NO_MEMORY_WAIT=1
experiment=$1; fold=$2; k=$3; methods=$4; ref=$5
job=$(pwd); fold_job=$(basename "$job"); key=${experiment}__${fold}
name=${HYPER_NAME:?set HYPER_NAME to the job name}
work=/workspace/ir_cube/ts_hyper/$name
threads=${RUNPOD_CPU_COUNT:-$(OMP_NUM_THREADS= nproc)}
say() { echo "[$(date -u '+%F %T')] $*"; }

[[ $ref =~ ^[0-9a-f]{40}$ ]] || { echo "the commit must be a full id" >&2; exit 1; }
say "set up $work"
mkdir -p "$work"; cd "$work"
curl -fsSL "https://codeload.github.com/${TS_CODE_REPO:-goldingn/ir_cube}/tar.gz/$ref" \
  | tar -xz --strip-components=1 --wildcards '*/R/' '*/tmb/'
mkdir -p outputs/two_stage temporary
ln -sfn "$(readlink -f "$job/data")" data
cp "$job"/outputs/bioassay_rho*.csv outputs/
stage_one_dir=${STAGE_ONE_DIR:-/workspace/ir_cube/ts_v5f/$fold_job/outputs/two_stage}
ln -sfn "$stage_one_dir/stage_one__$key.rds" \
  "outputs/two_stage/stage_one__$key.rds"
[ -e "outputs/two_stage/stage_one__$key.rds" ] || { echo "no stage-one cache" >&2; exit 1; }

if [ -n "${HYPER_FITS_FROM:-}" ]; then
  say "fits from $HYPER_FITS_FROM"
  mkdir -p outputs/two_stage/hyper/fits
  cp /workspace/ir_cube/ts_hyper/$HYPER_FITS_FROM/outputs/two_stage/hyper/fits/*.rds \
    outputs/two_stage/hyper/fits/
fi
say "compile the correction template"
Rscript -e 'source("R/two_stage_correction.R"); invisible(load_correction_template())'

# the methods other than the importance sampling ones (*_is) with every
# thread; those in forked single-threaded processes
other=$(echo "$methods" | tr ',' '\n' | grep -v '_is$' | paste -sd, || true)
is_methods=$(echo "$methods" | tr ',' '\n' | grep '_is$' | paste -sd, || true)
if [ -n "$other" ]; then
  say "methods $other, $threads threads"
  OPENBLAS_NUM_THREADS=$threads Rscript R/two_stage_hyper_folds.R \
    "$experiment" "$fold" "$k" "$other"
fi
if [ -n "$is_methods" ]; then
  say "methods $is_methods, $threads processes x 1 thread"
  OPENBLAS_NUM_THREADS=1 IR_CUBE_IS_CORES=$threads Rscript \
    R/two_stage_hyper_folds.R "$experiment" "$fold" "$k" "$is_methods"
fi
say "done"
