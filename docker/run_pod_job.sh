#!/usr/bin/env bash
# Run one fit on the pod, detached, from a private copy of the checkout.
#
#   run_pod_job.sh <name> fold <experiment> <fold> [flags]
#   run_pod_job.sh <name> full [flags]
#   run_pod_job.sh <name> --in <job> -- <command ...>
#
# flags: --threads N (default nproc), --chains N, --warmup N, --samples N
# (defaults: dynamical_mcmc_settings()). Other model options are an R
# expression in IR_CUBE_MODEL_OPTIONS, read by R/fit_model.R and
# R/validation_covariates.R. A fit runs in /workspace/ir_cube/jobs/<name>/
# with its own R/, tmb/, temporary/ and outputs/. --in runs any command (e.g.
# a second stage, predict.R) in an existing job's directory, on that job's fit.
# Each writes <prefix>.log, .pid, .spec, then <prefix>.DONE, or
# <prefix>.FAILED holding the exit status; the prefix is "job" for a fit and
# <name> for --in.
set -euo pipefail
root=/workspace/ir_cube; code=${CODE:-$root/code}
name=$1; kind=$2; shift 2

# a command in an existing job's directory
detach() {  # detach <prefix> <command ...>
  local prefix=$1; shift
  setsid nohup bash -c 'p=$1; shift; "$@" > "$p.log" 2>&1; s=$?
    echo "ended: $(date -u "+%F %T UTC"), status $s" >> "$p.spec"
    if [ $s -eq 0 ]; then touch "$p.DONE"; else echo $s > "$p.FAILED"; fi' \
    job "$prefix" "$@" < /dev/null > /dev/null 2>&1 &
  echo $! > "$prefix.pid"
}
if [ "$kind" = --in ]; then
  job=$root/jobs/$1; shift
  [ "${1:-}" = -- ] && shift
  [ $# -gt 0 ] || { echo "no command" >&2; exit 1; }
  cd "$job"
  [ ! -e "$name.spec" ] || { echo "$name already run in $job" >&2; exit 1; }
  printf 'command: %s\nstarted: %s\n' "$*" "$(date -u '+%F %T UTC')" > "$name.spec"
  detach "$name" "$@"
  echo "started $name in $job"
  exit 0
fi
args=()
[ "$kind" = fold ] && { args=("$1" "$2"); shift 2; }
threads=$(OMP_NUM_THREADS= nproc); chains=""; warmup=""; samples=""
while [ $# -gt 0 ]; do
  case $1 in
    --threads) threads=$2 ;; --chains) chains=$2 ;;
    --warmup) warmup=$2 ;; --samples) samples=$2 ;;
    *) echo "unknown flag $1" >&2; exit 1 ;;
  esac
  shift 2
done
case $kind in
  # run_one_fold.R takes chains, threads, warmup, samples by position
  fold) [ -z "$samples" ] || [ -n "$warmup" ] || { echo "--samples needs --warmup" >&2; exit 1; }
        command=(Rscript R/run_one_fold.R "${args[@]}" "${chains:-default}"
                 "$threads" $warmup $samples) ;;
  full) command=(Rscript R/fit_model.R "$threads")
        s=${chains:+n_chains = $chains, }${warmup:+warmup = $warmup, }${samples:+n_samples = $samples, }
        [ -z "$s" ] || export IR_CUBE_MCMC_SETTINGS="dynamical_mcmc_settings(${s%, })" ;;
  *) echo "kind is fold or full" >&2; exit 1 ;;
esac

job=$root/jobs/$name
mkdir -p "$root/jobs"; mkdir "$job"  # fails if the name is taken
mkdir -p "$job/temporary" "$job/outputs/cv_draws" "$job/figures"
cp -r "$code/R" "$code/tmb" "$job/"
ln -s "$root/data" "$job/data"
# copies, not links: fit_model.R rewrites the inits file
cp "$root"/temporary/*.RDS "$job/temporary/"
cp "$root"/outputs/bioassay_rho*.csv "$job/outputs/"
cat > "$job/job.spec" <<EOF
code: $(git -C "$code" rev-parse HEAD)
command: ${command[*]}
IR_CUBE_MODEL_OPTIONS: ${IR_CUBE_MODEL_OPTIONS:-}
IR_CUBE_MCMC_SETTINGS: ${IR_CUBE_MCMC_SETTINGS:-}
host: $(hostname) ${RUNPOD_POD_ID:-}, $(OMP_NUM_THREADS= nproc) vCPU
started: $(date -u '+%F %T UTC')
EOF
cd "$job"
detach job "${command[@]}"
echo "started $name: $job"
