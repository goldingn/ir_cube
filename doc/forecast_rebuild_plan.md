# Plan: rebuild the temporal forecasting experiment on five-year windows

> **Prompt to start work from this file.**
>
> Read `doc/forecast_rebuild_plan.md` and carry out the tasks in §4 in order.
> Work on branch `posterior-predictive-validation` (PR #12). Do not refit
> anything not listed in §4. Run all fitting from a frozen copy of `R/` as
> described in §3, and arm monitors before launching. Ask before posting
> anything to GitHub.

---

## 1. Where things stand

Branch `posterior-predictive-validation`, PR #12. All scoring machinery is
written and working; `R/validation_metrics.R` scores whatever is in
`outputs/cv_draws/` and writes the tables the figures read.

Folds already fitted, all 4 chains / 2,000 warmup / 5,000 samples:

| experiment | fold | held out | worst Rhat |
|---|---|---|---|
| spatial_extrapolation | 6 countries | 8,694 | 1.12 – 2.33 |
| spatial_interpolation | all | 1,045 | 1.18 |
| spatial_blocks | 1, 2 | 5,382 / 3,312 | 1.104 / 1.239 |
| temporal_forecasting | all (cut 2020, 3-year windows) | 1,461 | 1.107 |

The block folds and the corrected forecasting fold are new as of 22 Sept and
already scored. `outputs/cv_draws_leaky_forecast/` holds the pre-fix forecasting
fold; `outputs/cv_draws_defunct/` holds four 2-chain folds with no `draws`
object. Both are superseded.

## 2. The decision, and why

The forecasting experiment currently uses a three-year holdout starting 2020.
That window is the worst available choice, for two reasons found by
`R/fig_change_power.R`:

- **It is the one pause in the record.** Sliding a three-year holdout against
  the three years before it, the observed change at pixels assayed in both is
  decisively negative at every origin from 2005 to 2017 and flattens only at
  2018–2020. The 2020 origin reads +0.011 [−0.017, +0.040].
- **It is the thinnest window.** 280 paired (pixel, insecticide) groups and
  1,127 assays, against 1,436 and 6,346 at the 2013 origin. Enough to bound the
  mean change, nowhere near enough to score any single site — which is why the
  direction test came out at 49%, a coin flip, and must not be reported as a
  finding.

The model is not misbehaving on the trend: it predicts −0.093 over the window,
and the slope fitted to the training years 2012–2019 implies −0.075 to −0.082.
It extrapolates the historical rate faithfully; the rate stopped.

**Five-year windows fix this.** At every origin they give 30–50% more paired
pixels and a signal roughly 5/3 larger, because the gap between window midpoints
is the window length. Critically, a five-year window spans the 2018–2020 pause
*and* the decline either side of it, so the pause no longer swallows the test:
the 2018 origin goes from a signal-to-resolvable ratio of 0.2 at three years to
3.7 at five.

Covariates end in 2022, so the latest feasible five-year origin is 2018.

**Two origins, chosen:**

| cut | training | % of data | holdout | paired pixels | observed change | ratio to resolvable |
|---|---|---|---|---|---|---|
| **2014** | 1995–2013 | 52% | 2014–2018, 9,922 assays | 1,748 | −0.084 | 11.3 |
| **2018** | 1995–2017 | 82% | 2018–2022, 4,096 assays | 881 | −0.042 | 3.7 |

Non-overlapping holdouts bar the shared endpoint, training at half and
four-fifths of the data, and a factor of two between their true rates of
decline. That contrast is the test: does the model track a slowing rate, or
carry a fixed slope forward?

**Earlier origins are out.** A 2010 origin has the strongest signal in the
record but trains on 17% of the data, so it is not the model being deployed and
its forecast skill would not transfer. Same argument kills 2011 and 2012.

The existing 2020 three-year fold is **kept and reported as a supplementary
observation** — "the most recent window shows no decline" — not as the headline
forecast test. It is already paid for.

## 3. How to run fits, and the constraints that matter

**Environment — do not upgrade greta.** greta 0.6.0 cannot sample this model
(`while_loop` shape error via `iterate_dynamic_function`;
greta-dev/greta.dynamics#45). Working pin, in `R/packages.R`:

```r
remotes::install_version("tensorflow", version = "2.16.0")
remotes::install_github("njtierney/greta@4cc989f")            # 0.5.0.9000
remotes::install_github("greta-dev/greta.dynamics@db7df31")   # 0.2.2
```

**Two ordering constraints, both encoded at the top of `R/run_one_fold.R`:**
TensorFlow will not change its thread count after initialisation, so threads
must be set through `reticulate` before python starts; and python must
initialise before `terra` or `sf` are attached, or `tensorflow_probability`
fails to import.

**Run from a frozen copy of the scripts.** R reads `--file=` incrementally, so
editing a script mid-run corrupts that run — two folds were lost this way after
62 h each. Copy `R/` to a scratch directory and symlink `data/`, `outputs/`,
`temporary/`, `figures/`, `doc/` into it, then launch from there.

**Smoke test first.** `CV_DRAWS_DIR=outputs/cv_draws_smoke Rscript
R/run_one_fold.R temporal_forecasting 2014 2 3 5 5` takes a few minutes and
exercises the whole path. `CV_DRAWS_DIR` keeps the five-sample result out of
`outputs/cv_draws/`, where it would be mistaken for a fold and skipped. Delete
the smoke directory afterwards.

**Measured costs on this machine (16 cores, 30 GB):**

| configuration | per fold |
|---|---|
| 2 concurrent, 4 threads each | **62 h** |
| 3 concurrent, 3 threads each | 95 h |

So two fits concurrently is **~62 h total**, under three days. Peak resident
size is 2.5 GB during sampling and about 8 GB during `calculate()`; two
concurrent folds are comfortable, three were tight (available memory dipped to
3 GB).

**Sampling cannot be resumed across sessions.** `calculate(values = draws)`
works on a reloaded `draws` object; `extra_samples()` does not. So every
prediction target must be requested at fitting time — including the
before-window predictions the change score needs. A longer run must be asked
for up front.

**Arm monitors before launching**, on the fold logs for completions and failure
signatures, and on the process count so a silent death is caught. Use a pgrep
pattern that cannot match your own shell — `"exec/[R] .*run_one_fold"`.

## 4. Tasks, in order

### 4.1 Parameterise the forecasting experiment by cut year and window length

`R/validation_folds.R` currently hard-codes a three-year window ending at the
covariate limit. Replace with a function of cut year and window length:

- `training` = `year_start < cut`, strictly. (This is the leak fix from the #12
  review; keep the `stopifnot`.)
- `test` = `year_start %in% cut:(cut + window - 1)`
- `before` = `year_start %in% (cut - window):(cut - 1)`

Keep the existing `temporal_forecasting` object for the 2020 three-year fold so
the already-scored result still resolves, and add the parameterised one
alongside. Window length 5 and cuts 2014 and 2018.

`R/run_one_fold.R`: dispatch `temporal_forecasting` with `fold` as the cut year,
e.g. `Rscript R/run_one_fold.R temporal_forecasting 2014 4 4`. The existing 2020
fold is `fold = "all"`; either rename it in place (file and the `fold` field
inside the `.rds`, no refit needed) or leave it and accept the inconsistent
label — renaming is cleaner.

`R/run_validation_folds.R`: add the two new folds to the fold list, set
`n_concurrent <- 2` and `threads_per_fold <- 4`.

### 4.2 Rebuild the nulls for the two new origins

Minutes, no MCMC. `run_validation_folds.R` does this before dispatching, and
skips nulls already on disk. Check that both new origins get an `intercept` and
a `nearest_neighbour` file. `n_years_prior` for the nearest-neighbour null on a
forecasting fold is 3 in the current code; consider whether it should be 5 to
match the window, and if changed, say so in the write-up.

### 4.3 Smoke test, then launch

Smoke test both origins at 5 samples per §3, confirm `p_draws_before` comes back
with one column per before-window record and aligned with `before_df`, then
delete the smoke directory and launch the real run from a frozen copy.

~62 h. Check Rhat on the first fold to land; the block folds came in at 1.104
and 1.239 and the production fit is 1.215, so anything under about 1.25 is in
family.

### 4.4 Score

`Rscript R/validation_metrics.R` — picks up whatever is in `outputs/cv_draws/`,
takes about 40 minutes for 30 folds, and needs the memory fix already in place
(`score_fold` must not retain the whole fold object; retaining `draws` for all
folds took 26 GB and put the machine into swap).

Then `Rscript R/validation_change.R`, which currently reads only the single
forecasting fold. **Generalise it to loop over origins** and emit one row per
cut, so the rolling-origin comparison is a table rather than three runs.

Then `Rscript R/validation_geometry.R` for the distance-stratified tables.

### 4.5 Figures

`R/fig_change_power.R` already produces both power figures in ggplot, on
three-year windows. Regenerate with five-year windows:

- the sliding-origin figure (`figures/CV_change_power_windows.png`) — rebuild on
  five-year windows, which is the version that justifies the design;
- the paired-years figure (`figures/CV_change_power.png`) — keep the multi-gap
  panels, since the comparison of gaps is the point, but add a 5-year facet if
  not already present and mark which origins are feasible given the covariate
  limit.

Check with Nick what he wants from "a ggplot version of the figures" — both are
already ggplot, so he may mean something else: a version with the 5-year binning
throughout, or the underlying objects saved for reuse in the paper.

`R/fig_predictive_validation.R` then needs its forecasting panels to handle
several origins rather than one.

### 4.6 Update the documentation

`doc/cv_run_plan.md` §4.1 and §6 describe the three-year forecasting design;
update to the five-year rolling origin. Keep the record of why the 2020 window
was abandoned — it is the justification for the whole change.

## 5. Results to carry forward

Scored and current as of 22 Sept, all models at the external replicate-based
overdispersion, skill anchored on the per-insecticide intercept null, pixel
cluster bootstrap:

| experiment | model | excess MSE | variance explained |
|---|---|---|---|
| spatial blocks | dynamical | 0.061 | +0.06 [−0.03, +0.13] |
| spatial blocks | nearest neighbour | 0.120 | −0.86 [−1.05, −0.68] |
| interpolation | dynamical | 0.042 | +0.19 [−0.02, +0.37] |
| interpolation | nearest neighbour | 0.027 | +0.48 [+0.32, +0.62] |
| extrapolation | dynamical | 0.087 | −0.31 [−0.42, −0.22] |
| extrapolation | nearest neighbour | 0.087 | −0.31 [−0.44, −0.19] |
| forecasting (2020, 3-year) | dynamical | 0.068 | +0.26 [+0.05, +0.45] |
| forecasting (2020, 3-year) | nearest neighbour | 0.049 | +0.47 [+0.32, +0.59] |

Three findings that do not depend on the forecasting rebuild:

- **The nearest neighbour collapses at distance.** −0.86 on the block folds,
  where held-out pixels sit a median 77 km from training data, against +0.48 on
  interpolation where the median is 36 km. Decisive in both block folds.
- **Leave-one-country-out measures the wrong thing.** Held-out bias correlates
  with the country's fitted effect at r = −0.94 and excess MSE with its
  magnitude at r = +0.91; Ethiopia and Kenya carry the whole pooled result. It
  measures the difficulty of an entirely unsampled country, which the deployed
  model never faces.
- **The dynamical model's lowest predictions are too low.** On interpolation the
  lowest predicted decile predicts 0.178 against 0.442 observed, seven times
  outside the envelope its own posterior predictive distribution would produce,
  while both nulls are inside theirs. This is issue #14, a floor on mortality in
  a fully resistant population.

## 6. Traps already hit

- `calculate(..., nsim = n)` returns an independent resample and destroys MCMC
  ordering; use `calculate(values = draws)` with no `nsim`.
- Reliability must be binned on the **prediction**, never on observed mortality:
  binning on the noisy outcome makes a perfectly calibrated model show +0.15
  bias in the lowest decile. `reliability_ppc()` supplies the correct envelope.
- The stored PIT must be one randomisation replicate, `pit[, 1]`, not the mean
  over replicates — the mean converges to the mid-P value, which is not uniform
  for discrete data.
- In `dplyr::mutate()`, a column created earlier in the call shadows a variable
  of the same name later in it. This silently reduced a distance matrix to its
  vector of minima, and silently broke two weighted means.
- `pgrep -f` matches your own shell's command line. Use a bracketed character
  class.
- A completed background job was reported as still running for six days because
  `pgrep` matched the waiting shell. Check the output, not the process.
