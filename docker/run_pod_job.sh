#!/usr/bin/env bash
# Run one fit on the pod, from a private copy of the code.
#
#   run_pod_job.sh <name> fold <experiment> <fold> [flags]
#   run_pod_job.sh <name> full [flags]
#   run_pod_job.sh <name> --in <job> -- <command ...>
#
# flags: --threads N (default the pod's vCPUs), --chains N, --warmup N,
# --samples N (defaults: dynamical_mcmc_settings()). Other model options are
# an R expression in IR_CUBE_MODEL_OPTIONS, read by R/fit_model.R and
# R/validation_covariates.R; IR_CUBE_INITS names the cached initial values,
# comma-separated, for a share of the chains each (dynamical_inits_files()). A fit runs in /workspace/ir_cube/jobs/<name>/
# with its own R/, tmb/, temporary/ and outputs/, the code of the commit
# CODE_REF (a full commit id) of CODE_REPO (default idem-lab/ir_cube),
# downloaded from GitHub. --in runs any command (e.g. a second stage,
# predict.R) in an existing job's directory, on that job's fit and with its
# model options. Each writes <prefix>.log, .pid, .spec, then <prefix>.DONE,
# or <prefix>.FAILED holding the exit status; the prefix is "job" for a fit
# and <name> for --in.
#
# Commands run detached, so that they outlive an ssh session, unless
# RUN_POD_JOB_FOREGROUND=1, as when the image's entrypoint runs a job given at
# pod creation; then the log also goes to the container's output.
set -euo pipefail
root=/workspace/ir_cube
name=$1; kind=$2; shift 2

# run a command, writing <prefix>.log and the markers
run() {  # run <prefix> <command ...>
  local prefix=$1; shift
  if [ "${RUN_POD_JOB_FOREGROUND:-0}" = 1 ]; then
    echo $$ > "$prefix.pid"
    set +e; "$@" 2>&1 | tee "$prefix.log"; local s=${PIPESTATUS[0]}; set -e
    echo "ended: $(date -u "+%F %T UTC"), status $s" >> "$prefix.spec"
    if [ $s -eq 0 ]; then touch "$prefix.DONE"; else echo $s > "$prefix.FAILED"; fi
    return $s
  fi
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
  # the fit's model options, so that the command rebuilds the fit's design
  fit_options=$(cat job.options)
  if [ -n "${IR_CUBE_MODEL_OPTIONS+set}" ] &&
     [ "$IR_CUBE_MODEL_OPTIONS" != "$fit_options" ]; then
    echo "IR_CUBE_MODEL_OPTIONS differs from the fit's" >&2; exit 1
  fi
  if [ -n "$fit_options" ]; then export IR_CUBE_MODEL_OPTIONS=$fit_options; fi
  printf 'command: %s\nIR_CUBE_MODEL_OPTIONS: %s\nstarted: %s\n' "$*" \
    "$fit_options" "$(date -u '+%F %T UTC')" > "$name.spec"
  echo "started $name in $job"
  run "$name" "$@"
  exit
fi
args=()
[ "$kind" = fold ] && { args=("$1" "$2"); shift 2; }
# RunPod sets RUNPOD_CPU_COUNT; nproc can report the host's cores
threads=${RUNPOD_CPU_COUNT:-$(OMP_NUM_THREADS= nproc)}; chains=""; warmup=""; samples=""
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
        [ -z "${IR_CUBE_MCMC_SETTINGS:-}" ] || {
          echo "folds take --chains, --warmup and --samples, not IR_CUBE_MCMC_SETTINGS" >&2; exit 1; }
        command=(Rscript R/run_one_fold.R "${args[@]}" "${chains:-default}"
                 "$threads" $warmup $samples) ;;
  full) command=(Rscript R/fit_model.R "$threads")
        s=${chains:+n_chains = $chains, }${warmup:+warmup = $warmup, }${samples:+n_samples = $samples, }
        [ -z "$s" ] || export IR_CUBE_MCMC_SETTINGS="dynamical_mcmc_settings(${s%, })" ;;
  *) echo "kind is fold or full" >&2; exit 1 ;;
esac

[[ ${CODE_REF:-} =~ ^[0-9a-f]{40}$ ]] || {
  echo "CODE_REF must be a full commit id" >&2; exit 1; }
repo=${CODE_REPO:-idem-lab/ir_cube}
job=$root/jobs/$name
mkdir -p "$root/jobs"; mkdir "$job"  # fails if the name is taken
mkdir -p "$job/temporary" "$job/outputs/cv_draws" "$job/figures"
curl -fsSL "https://codeload.github.com/$repo/tar.gz/$CODE_REF" \
  | tar -xz -C "$job" --strip-components=1 --wildcards '*/R/' '*/tmb/'
ln -s "$root/data" "$job/data"
# copies, not links: fit_model.R rewrites the inits file
cp "$root"/temporary/*.RDS "$job/temporary/"
cp "$root"/outputs/bioassay_rho*.csv "$job/outputs/"
printf '%s' "${IR_CUBE_MODEL_OPTIONS:-}" > "$job/job.options"
cat > "$job/job.spec" <<EOF
code: $repo@$CODE_REF
inputs (sha256): bioassays $(sha256sum < "$root/data/clean/all_gambiae_complex_data.RDS" | cut -c1-16), initial values $(sha256sum < "$root/temporary/inits_refit.RDS" | cut -c1-16)
command: ${command[*]}
IR_CUBE_MODEL_OPTIONS: ${IR_CUBE_MODEL_OPTIONS:-}
IR_CUBE_MCMC_SETTINGS: ${IR_CUBE_MCMC_SETTINGS:-}
IR_CUBE_INITS: ${IR_CUBE_INITS:-}
host: $(hostname) ${RUNPOD_POD_ID:-}, ${RUNPOD_CPU_COUNT:-$(OMP_NUM_THREADS= nproc)} vCPU
started: $(date -u '+%F %T UTC')
EOF
cd "$job"
echo "started $name: $job"
run job "${command[@]}"
