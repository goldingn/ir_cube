# Running fits on RunPod

One greta fit (the full fit or one CV fold) per RunPod CPU pod. The run recipe
(which folds, inits, post-processing) is in `doc/cv_run_plan.md`.

## Image

`docker/Dockerfile`: rocker/geospatial 4.6.1, R packages from a dated posit
snapshot, greta at `282944f`, python 3.12 with TensorFlow 2.21 and TensorFlow
Probability 0.25 (`docker/requirements.txt`), OpenBLAS, and sshd taking the
key RunPod passes in `PUBLIC_KEY` (`docker/entrypoint.sh`), with
`run_pod_job.sh`. It holds software only; the code is
downloaded per job and the data live on the network volume. Checked against the local
greta 0.6 setup: `R/check_dynamical_model.R` agreed to every digit.

The image in use is `ghcr.io/goldingn/ir_cube-runpod:latest` (public). To
rebuild, run `.github/workflows/runpod-image.yml` by hand, which pushes
`ghcr.io/<owner>/ir_cube-runpod:<tag>` and `:<short sha>`:

```bash
gh workflow run runpod-image.yml --ref <branch> -f tag=latest
```

or build it locally and push it with `docker build -f docker/Dockerfile .`
and `docker push`.

## Pods and volume

Network volume `ir-cube` (id `0f6xjzaxch`, EU-RO-1, 50 GB, about $3.50 a
month), mounted at `/workspace`. Under `/workspace/ir_cube/`: `data/`,
`temporary/` (the RDS inputs),
`outputs/` (`bioassay_rho*.csv`) and `jobs/`. Pods must be in EU-RO-1 to
mount it.

| pod | threads | one fit, 4 chains, 2000 + 3000 | cost per fit |
|---|---|---|---|
| cpu5c, 8 vCPU, 16 GB (default) | 8 | 3.7 h (measured) | about $1.05 |
| cpu5c, 16 vCPU, 32 GB | 16 | about 3 h | about $1.70 |

cpu5c capacity in EU-RO-1 can run out (32 vCPU once); cpu3c is the fallback.
More chains with fewer samples each did not pay off: warmup dominates.

## Jobs

Routine runs need no ssh. Each job is given when its pod is created, its log
is the pod's log, and its results are on the volume. Data go to and from the
volume through its S3-compatible API.

**Once per machine:** install `docker/irpod` outside the repository, so that
the command a permission rule allows can't change with an edit to the
repository, and add an S3 API key from the RunPod console (Settings, S3 API
Keys) as the profile `runpod` in `~/.aws/credentials`:

```bash
python3 -m venv ~/.local/share/irpod/venv
~/.local/share/irpod/venv/bin/pip install awscli
install -D -m 555 docker/irpod ~/.local/bin/irpod
```

**Inputs:** from the repository root, `irpod sync` uploads `data/clean/`, the
UNSD table, `temporary/*.RDS` and `outputs/bioassay_rho*.csv`. Sync again
after `fit_model.R` rewrites the cached initial values.

**A job:** create a pod with the image, the volume at `/workspace`, and these
environment variables:

| variable | value |
|---|---|
| `JOB` | the arguments of `run_pod_job.sh`: `<name> fold <experiment> <fold>`, `<name> full`, or `<name> --in <job> -- <command>`, with `--threads`, `--chains`, `--warmup`, `--samples` |
| `CODE_REF` | the commit to run (a full commit id), downloaded from GitHub |
| `CODE_REPO` | if not `idem-lab/ir_cube`, e.g. a fork for a PR branch |
| `IR_CUBE_MODEL_OPTIONS` | an R expression for `dynamical_model_options()`, if not the defaults (V5h since #47: population d½ 270, a floor per class shifted by a latent smooth, a latent smooth of selection, the beta-binomial likelihood); give it in full anyway, so that a later change of default cannot change the job, and so that a fold matches the full fit it is paired with. It may contain spaces, e.g. `dynamical_model_options(mortality_floor = FALSE, smooth = FALSE)` for the default before V5h, `dynamical_model_options(mortality_floor = FALSE, smooth = FALSE, species = species_options())` for the species model (#47), which reads `data/clean/arabiensis_fraction.tif`, or `smooth = FALSE, kdr = kdr_options()` for the kdr covariate, which reads `data/clean/kdr_total_2015.tif`. The latent smooths read no further input |
| `IR_CUBE_INITS` | cached initial values under `temporary/`, comma-separated, each for an equal share of the chains (e.g. both floor modes, `R/floor_mode_inits.R`, or one draw of a full fit per chain, `R/draw_inits.R`), if not `temporary/inits_refit.RDS` |
| `OMP_NUM_THREADS`, `IR_CUBE_NO_MEMORY_WAIT` | as the job needs |

`JOB` is split on spaces, so its arguments cannot contain spaces or quotes.
An `--in` job takes its fit's `IR_CUBE_MODEL_OPTIONS`. For example,
`JOB="main full --threads 8"`,
`JOB="blocks1 fold spatial_blocks 1 --threads 8"`, or, on a finished fit,
`JOB="ts_blocks1 --in blocks1 -- Rscript R/run_two_stage_folds.R spatial_blocks 1 stage_one"`.
The steps of a second stage and their order are in each script's header.
Loading a saved fold takes 4 to 8 GB, so on a 16 GB pod run such steps one
at a time.

The job runs in `jobs/<name>/` on the volume, with its own `R/`, `tmb/`,
`temporary/` and `outputs/`, and writes `job.log`, `job.DONE` or
`job.FAILED` (for `--in`, `<name>.log` and so on), and `job.spec`: the commit,
checksums of the bioassay data and initial values (the job reads the volume's
live `data/`, which a later sync changes), the command and the options. A name already taken
fails at once, so a restarted container does not run the job twice. Then:

```bash
irpod status <name>    # spec, markers and the end of the log
irpod fetch <name>     # results to outputs/pod_jobs/<name>/
```

A pod cannot delete itself (RunPod refuses its pod-scoped key), so it idles,
billed, until it is deleted. When `irpod status` shows `job.DONE` or
`job.FAILED`, delete the pod with the RunPod connector; its log ends "job
finished". At the end of a session, check `list-pods` for pods left running.

**Interactively**, with ssh (`~/.ssh/runpod_ed25519`, the pod's direct TCP
port 22), `run_pod_job.sh` runs a job detached, so that it outlives the
session, with the same variables.

## Pitfalls

- Over ssh, detach every job (the launcher does): an ssh or agent session
  ending kills its children.
- A new pod pulls the image first, which took 5 to 25 minutes.
- In a container `/proc/meminfo` is the host's; set `IR_CUBE_NO_MEMORY_WAIT=1`
  for scripts that call `wait_for_memory()`.
- zsh does not word-split a variable holding ssh options; wrap in `bash -c`
  or write them out.
- Pods bill until deleted, idle or not. The volume is billed monthly whether
  or not a pod is running.
