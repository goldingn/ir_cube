# V5f cross-validation and two-stage maps: runbook (#47)

The test fits, the dynamical folds and the two-stage model on each fold of
V5f, the default model since #47 (V5h with the floor smooth centred over the
cells with bioassays of the pyrethroids or DDT, and a half-normal prior of
scale 0.05 on the floor at a flat smooth), and the two-stage maps of its full
fit, on RunPod (`docker/README.md`). Pod creation needs Nick's approval;
state the price first.

## Code and inputs

| what | where |
|---|---|
| fit and fold code | `goldingn/ir_cube`, branch `weighted-binomial` (`CODE_REF`: its head) |
| second-stage code | branch `v5h-cv`: `weighted-binomial` merged with `cv-update` (PR #43), for the map draws (`map_draws`, `map_before`) that the #43 scoring reads. Its dynamical code is `weighted-binomial`'s; `stage_one` checks the recomputed draws against each fold's to 1e-6 |
| options | the defaults, in full in `IR_CUBE_MODEL_OPTIONS` (`v5f_options` in `R/species_runs.R`) |
| initial values | `temporary/inits_floor_low.RDS` (each class's floor at a flat smooth starts at its floor, 0.0015). The warm start, a fallback: `temporary/inits_v5f_draw1.RDS` to `draw4.RDS`, one draw of `v5f_full` per chain (`R/draw_inits.R`). `irpod sync` uploads them |
| scripts on the volume | `ts_v5f/two_stage_fold.sh`, `ts_v5f/two_stage_maps.sh`: copies of `docker/two_stage_fold.sh` and `docker/two_stage_maps.sh` |

The pod bodies come from the launcher, which checks that the options equal
`dynamical_model_options()` at the commit it is run at:

```bash
Rscript R/species_runs.R <CODE_REF> v5f            # outputs/species_runs/pods_v5f.json: v5f_full, v5f_fc2018, v5f_fc2014
Rscript R/species_runs.R <CODE_REF> cv5f           # pods_cv5f.json: the three spatial folds
Rscript R/species_runs.R <CODE_REF> cv5f_warm      # pods_cv5f_warm.json: the fallback, after R/draw_inits.R
Rscript R/species_runs.R <V5H_CV_REF> ts_cv5f      # pods_ts_cv5f.json: the second stage on the five folds
Rscript R/species_runs.R <V5H_CV_REF> maps_v5f     # pods_maps_v5f.json: the maps of v5f_full
```

`CV5F_OPTIONS`, `CV5F_NAME`, `CV5F_FULL` and `CV5F_INITS` override the
options, the folds' name prefix, the full fit and the warm start's files
(`R/species_runs.R`), e.g. after another change to the model. Upload the
scripts (again after any change to them) with

```bash
aws=~/.local/share/irpod/venv/bin/aws
s3() { $aws --profile runpod --region eu-ro-1 --endpoint-url https://s3api-eu-ro-1.runpod.io s3 "$@"; }
s3 cp docker/two_stage_fold.sh s3://0f6xjzaxch/ir_cube/ts_v5f/two_stage_fold.sh
s3 cp docker/two_stage_maps.sh s3://0f6xjzaxch/ir_cube/ts_v5f/two_stage_maps.sh
```

## 0. Test fits, and the forecasting folds

`v5f_full` (`JOB="v5f_full full --threads 8"`), `v5f_fc2018`
(`JOB="v5f_fc2018 fold temporal_forecasting 2018 --threads 8"`) and
`v5f_fc2014` (`JOB="v5f_fc2014 fold temporal_forecasting 2014 --threads 8"`),
cpu5c, 8 vCPU, 16 GB, $0.28/h, 4 chains, 2,000 + 1,500, from the low-floor
initial values (all three launched 9 October at `933f8ae`). V5h's took 1.54 h
(full, cpu5c), 0.95 h (2014 fold, cpu5c) and 2.4 h (2018 fold, cpu3c); with
the image pull about 1.2-2.7 h, $0.35-0.75 a pod. Check the 2018 fold's
chains: per-chain means of the floors at a flat smooth (`floor_at_u0`) and
of the smooth sds, and the log posterior per chain (`chain_modes` in the
fold file), and rank Rhat. The two forecasting folds are the
cross-validation's; the full fit is the one the maps are made from (step 3),
and the warm start's draws.

## 1. Spatial folds

Three pods, one fold each, as the test jobs: cpu5c, 8 vCPU, $0.28/h.

| job | `JOB` |
|---|---|
| `cv5f_blocks1` | `cv5f_blocks1 fold spatial_blocks 1 --threads 8` |
| `cv5f_blocks2` | `cv5f_blocks2 fold spatial_blocks 2 --threads 8` |
| `cv5f_interp` | `cv5f_interp fold spatial_interpolation all --threads 8` |

each with `IR_CUBE_MODEL_OPTIONS` (the options in full) and
`IR_CUBE_INITS=temporary/inits_floor_low.RDS`. About 1-1.6 h a fold plus
the pull: about $1.50 for the three.

**If the test fold does not converge,** the warm start: when `v5f_full` is
done, `irpod fetch v5f_full`, then

```bash
Rscript R/draw_inits.R outputs/pod_jobs/v5f_full/temporary/fitted_model.RData temporary/inits_v5f_draw
irpod sync
Rscript R/species_runs.R <CODE_REF> cv5f_warm
```

which refits all five folds, `cv5f_warm_<fold>` (`blocks1`, `blocks2`,
`interp`, `fc2014`, `fc2018`), chain i at draw i of the full fit, one from
each of its chains, spread along the pyrethroid floor at a flat smooth: 4
chains, one init file each; about $2.50. Its second stage is `ts_cv5f_warm`.

When `irpod status <job>` shows `job.DONE` (or `job.FAILED`): `irpod fetch
<job>`, then delete the pod. Each fold saves
`outputs/cv_draws/dynamical__<experiment>__<fold>.rds`, with its draws, its
chains' floors and log posteriors (`chain_modes`) and Rhat
(`convergence`).

## 2. Two-stage folds

One pod per fold, once its dynamical fold is done (its pod can be deleted
first): cpu3c, 16 vCPU, 32 GB, $0.48/h, no `IR_CUBE_MODEL_OPTIONS` (`--in`
takes the fold's). The five, on the fold jobs `cv5f_blocks1`,
`cv5f_blocks2`, `cv5f_interp`, `v5f_fc2014` and `v5f_fc2018`:

| job | `JOB` |
|---|---|
| `ts_cv5f_blocks1` | `ts_cv5f_blocks1 --in cv5f_blocks1 -- bash /workspace/ir_cube/ts_v5f/two_stage_fold.sh spatial_blocks 1 <V5H_CV_REF>` |
| `ts_cv5f_blocks2` | `ts_cv5f_blocks2 --in cv5f_blocks2 -- bash /workspace/ir_cube/ts_v5f/two_stage_fold.sh spatial_blocks 2 <V5H_CV_REF>` |
| `ts_cv5f_interp` | `ts_cv5f_interp --in cv5f_interp -- bash /workspace/ir_cube/ts_v5f/two_stage_fold.sh spatial_interpolation all <V5H_CV_REF>` |
| `ts_v5f_fc2014` | `ts_v5f_fc2014 --in v5f_fc2014 -- bash /workspace/ir_cube/ts_v5f/two_stage_fold.sh temporal_forecasting 2014 <V5H_CV_REF>` |
| `ts_v5f_fc2018` | `ts_v5f_fc2018 --in v5f_fc2018 -- bash /workspace/ir_cube/ts_v5f/two_stage_fold.sh temporal_forecasting 2018 <V5H_CV_REF>` |
 The script works in `ts_v5f/<job>/` with the `v5h-cv` code:
`stage_one`, the nine types in two queues (the largest peak at about 14 GB),
then `assemble`; it copies `two_stage__<experiment>__<fold>.rds` and
`outputs/two_stage/` (`fit_summary.csv`, `parts/`, logs) into the fold's job
directory. The October two-stage folds of PR #43 took 50-58 min (blocks,
interpolation), 20 min (2014) and 35 min (2018) on this pod, plus the pull:
about $2.40 for five. When `ts_<job>.DONE` (or `.FAILED`) appears in
`jobs/<job>/` (`irpod status <job>` lists it): `irpod fetch <job>` again, then
delete the pod.

## 3. Two-stage maps of the full fit

One pod: cpu3c, 16 vCPU, 32 GB, $0.48/h, no `IR_CUBE_MODEL_OPTIONS`:

```
JOB="maps_v5f_full --in v5f_full -- bash /workspace/ir_cube/ts_v5f/two_stage_maps.sh <V5H_CV_REF>"
```

It needs only the full fit, so it can run beside the folds. The script works
in `ts_v5f/maps_v5f_full/` (the `--in` writes only its markers and log in
`jobs/v5f_full/`): `prepare` (minutes, 7 GB), `fit` of the nine types in two
queues (7-20 min and 8-14 GB each; about 1.3 h), `map` of the six types not
in LLINs and of `llin_effective`, one at a time with 4 workers (about 20 GB;
on the Linux box 22-30 min a type and 81 min for `llin_effective`; about 4
h), and `figures`. About 6 h, $3. `TS_STEPS` (e.g. `"map figures"`) resumes
from a step. The outputs add about 10 GB to the 50 GB volume (fits 4.3 GB, 720
rasters). Fetch them with

```bash
s3 sync s3://0f6xjzaxch/ir_cube/ts_v5f/maps_v5f_full/ outputs/pod_jobs/maps_v5f_full/ \
  --exclude 'R/*' --exclude 'tmb/*' --exclude 'data' --exclude 'data/*' --exclude 'temporary/*'
```

## Totals

| step | pods | wall time | cost |
|---|---|---|---|
| test fits and forecasting folds | 3 × cpu5c 8 | 1.2-2.7 h | $1.50-2 |
| spatial folds | 3 × cpu5c 8 | 1.2-1.9 h | about $1.50 |
| two-stage folds | 5 × cpu3c 16 | 0.5-1.3 h | about $2.40 |
| two-stage maps | 1 × cpu3c 16 | about 6 h | about $3 |

About $9 in all. With the spatial folds and the maps started once the test
fits are done, and each second stage once its fold is, about 9 h of wall
time from the test fits' start, most of it the maps.

## 4. Afterwards, locally

In a worktree of `v5h-cv`, from a frozen copy of `R/` (`doc/cv_run_plan.md`,
section 3), with the current data:

1. Put the five folds in `outputs/cv_draws/`: `dynamical__*.rds` and
   `two_stage__*.rds` from `outputs/pod_jobs/<job>/outputs/cv_draws/` of
   `cv5f_blocks1`, `cv5f_blocks2`, `cv5f_interp`, `v5f_fc2014` and
   `v5f_fc2018`. Do
   this before `run_validation_folds.R`, which otherwise fits the missing
   dynamical folds locally.
2. `outputs/two_stage/fit_summary.csv`: the rows of the five folds' own
   `fit_summary.csv` (each holds only its fold's), bound together.
3. For the maps: `outputs/two_stage/maps/` and `outputs/two_stage/ir_maps/`
   from the fetched work directory, and the full fit as
   `temporary/fitted_model.RData`.
4. The run order of `doc/cv_run_plan.md`, section 3 (#43):
   `fig_illustrate_bioassay_variability.R` (the rho tables),
   `run_validation_folds.R` (the nulls only, once the folds are in place),
   `validation_metrics.R`, `validation_geometry.R`, `variance_explained.R`,
   `two_stage_metrics.R`, `bioassay_vs_map.R`, `full_fit_maps_at_assays.R`,
   `validation_level.R`, `validation_change.R`, `fig_variance_explained.R`,
   `fig_predictive_validation.R`; and the published map figures
   (`fig_ir_maps.R` and those of `doc/two_stage_plan.md`, step 6), which read
   `outputs/two_stage/ir_maps/`.
