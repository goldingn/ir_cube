# Running fits on RunPod

One greta fit (the full fit or one CV fold) per RunPod CPU pod. The run recipe
(which folds, inits, post-processing) is in `doc/cv_run_plan.md`.

## Image

`docker/Dockerfile`: rocker/geospatial 4.6.1, R packages from a dated posit
snapshot, greta at `282944f`, python 3.12 with TensorFlow 2.21 and TensorFlow
Probability 0.25 (`docker/requirements.txt`), OpenBLAS, and sshd taking the
key RunPod passes in `PUBLIC_KEY` (`docker/entrypoint.sh`). It holds software
only; code and data live on the network volume. Checked against the local
greta 0.6 setup: `R/check_dynamical_model.R` agreed to every digit.

The image in use is `ghcr.io/goldingn/ir_cube-runpod:latest` (public). To
rebuild, run `.github/workflows/runpod-image.yml` by hand, which pushes
`ghcr.io/<owner>/ir_cube-runpod:<tag>` and `:<short sha>`:

```bash
gh workflow run runpod-image.yml --ref <branch> -f tag=latest
```

## Pods and volume

Network volume `ir-cube` (id `0f6xjzaxch`, EU-RO-1, 50 GB, about $3.50 a
month), mounted at `/workspace`. Under `/workspace/ir_cube/`: `code/` (a clone,
checked out at the commit to run), `data/`, `temporary/` (the RDS inputs),
`outputs/` (`bioassay_rho*.csv`), `run_pod_job.sh` and `jobs/`. Pods must be
in EU-RO-1 to mount it. Pod template: the image above, TCP port 22 exposed,
volume at `/workspace`; ssh with `~/.ssh/runpod_ed25519` to the IP and port of
the pod's direct TCP mapping.

| pod | threads | one fit, 4 chains, 2000 + 3000 | cost per fit |
|---|---|---|---|
| cpu5c, 8 vCPU, 16 GB (default) | 8 | 3.7 h (measured) | about $1.05 |
| cpu5c, 16 vCPU, 32 GB | 16 | about 3 h | about $1.70 |

cpu5c capacity in EU-RO-1 can run out (32 vCPU once); cpu3c is the fallback.
More chains with fewer samples each did not pay off: warmup dominates.

## Jobs

Sync inputs from the local repository (`-L`: `data/raw` and `temporary/` hold
symlinks), then update the volume's checkout:

```bash
rsync -avL -e "ssh -p <port> -i ~/.ssh/runpod_ed25519" data temporary \
  root@<ip>:/workspace/ir_cube/
ssh ... 'git -C /workspace/ir_cube/code fetch && git -C /workspace/ir_cube/code checkout <commit>'
```

`docker/run_pod_job.sh` runs one fit detached, in `jobs/<name>/` with its own
copy of `R/` and `tmb/`, so the checkout can change while it runs:

```bash
J=/workspace/ir_cube/run_pod_job.sh   # copy it there once
$J blocks1 fold spatial_blocks 1 --threads 8
$J main full --threads 8
IR_CUBE_MODEL_OPTIONS='dynamical_model_options(selection_columns = selection_design(net_w = 0.47))' \
  $J netw047 fold spatial_interpolation all --threads 8
```

Check `jobs/<name>/job.log`, and `job.DONE` or `job.FAILED`; stop with
`kill -- -$(cat jobs/<name>/job.pid)`.

`--in <job>` runs any command detached in an existing job's directory, on its
fit: second stages, `predict.R`, figures, scoring. The steps and their order
are in each script's header. Its files are `<name>.log`, `.pid`, `.spec`,
`.DONE` or `.FAILED`:

```bash
$J ts_blocks1 --in blocks1 -- Rscript R/run_two_stage_folds.R spatial_blocks 1 stage_one
IR_CUBE_NO_MEMORY_WAIT=1 $J ts_maps --in main -- Rscript R/two_stage_maps.R prepare
```

Loading a saved fold takes 4 to 8 GB, so on a 16 GB pod run such steps one at
a time. Fetch results with

```bash
rsync -av --exclude R/ --exclude tmb/ --exclude data \
  -e "ssh -p <port> -i ~/.ssh/runpod_ed25519" \
  root@<ip>:/workspace/ir_cube/jobs/<name> outputs/pod_jobs/
```

## Pitfalls

- Detach every job (the launcher does): an ssh or agent session ending kills
  its children.
- A new pod pulls the image first, which took 5 to 25 minutes.
- In a container `/proc/meminfo` is the host's; set `IR_CUBE_NO_MEMORY_WAIT=1`
  for scripts that call `wait_for_memory()`.
- zsh does not word-split a variable holding ssh options; wrap in `bash -c`
  or write them out.
- Agents create pods with the RunPod connector. `.claude/settings.local.json`
  needs allow rules for `mcp__claude_ai_runpod__create-pod`, `delete-pod`,
  `pod-action`, `Bash(ssh -i ~/.ssh/runpod_ed25519:*)` and `Bash(rsync:*)`.
- Terminate pods when their jobs finish. The volume is billed monthly whether
  or not a pod is running.
