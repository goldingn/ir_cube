# Cross-validation: run recipe, and the rebuild after the #12 review

Written for whoever runs the next set of fits. Sections 1–3 are the recipe and
the constraints it exists to satisfy; sections 4–7 are the rebuild the review of
PR #12 asks for, and what it costs.

---

## 1. Environment: do not upgrade greta

Sampling this model fails on greta 0.6.0 — `mcmc()` raises a TensorFlow
`while_loop` shape error whenever the likelihood depends on
`iterate_dynamic_function()` output. Reported upstream as
greta-dev/greta.dynamics#45.

The working combination, verified for both master's model definition and
`fit_fold()`:

```r
remotes::install_version("tensorflow", version = "2.16.0")
remotes::install_github("njtierney/greta@4cc989f")            # 0.5.0.9000
remotes::install_github("greta-dev/greta.dynamics@db7df31")   # 0.2.2
```

Python side is a conda env built by `greta::install_greta_deps()`:
TensorFlow 2.15.1, TFP 0.23.0.

## 2. Two ordering constraints, both load-bearing

Both are encoded at the top of `R/run_one_fold.R`, and the file will fail in
confusing ways if they are disturbed:

- TensorFlow will not change its thread count after initialisation, so threads
  must be set **before** python comes up. greta exposes no interface for this;
  it must go through `reticulate::import("tensorflow")$config$threading$...`.
  The environment variables (`TF_NUM_INTRAOP_THREADS` and friends) are ignored.
- python must initialise **before** `terra` or `sf` are attached. Those load the
  system XML libraries, against which the conda environment's `pyexpat` is then
  resolved, and `tensorflow_probability` fails to import.

## 3. Running the folds

```bash
Rscript R/run_validation_folds.R        # nulls, then dispatch the model folds
Rscript R/validation_metrics.R          # score everything on disk
Rscript R/validation_geometry.R         # skill against distance and data volume
Rscript R/fig_predictive_validation.R   # figures and the table
```

`run_validation_folds.R` fits the nulls in process (minutes) and then dispatches
each model fold as a separate `Rscript R/run_one_fold.R <experiment> <fold>`,
two at a time, logging to `outputs/cv_logs/<experiment>__<fold>.log`. A fold
whose `.rds` is already in `outputs/cv_draws/` is skipped, by `run_one_fold.R`
itself, so the run resumes after an interruption.

**Run from a frozen copy of the scripts.** R reads `--file=` incrementally, so
editing a script while a fold is running corrupts that fold — two folds were
lost this way after 62 h each, both crashing at the `saveRDS` block. Copy `R/`
to a scratch directory, symlink `data/`, `outputs/` and `temporary/` into it,
and launch from there.

Smoke test the whole path before committing to a long run:

```bash
Rscript R/run_one_fold.R spatial_extrapolation Kenya 2 4 5 5
```

takes a few minutes and exercises everything. Delete the resulting `.rds`
afterwards, or the real fold will be skipped.

### Sampling settings

4 chains, 2,000 warmup, 5,000 post-warmup samples, `hmc(Lmin = 15, Lmax = 30)`,
initialised from `temporary/inits.RDS`. Roughly 62 h per fold at four threads,
two folds at a time.

Four chains rather than two because greta pools information across chains when
adapting during warmup: at two chains the Kenya fold reached Rhat 7.6, and at
four it reached 1.15. Chains cost super-linearly in this model (6.96, 15.21 and
70.61 s per iteration at 2, 4 and 8 chains), and TensorFlow threads scale poorly
beyond about four (8.39 s/iteration at 2 threads against 6.96 s using the whole
machine), which is why the budget goes into concurrent folds rather than into
threads.

The September run reached this same 5,000 samples by taking 500 and topping up
with `extra_samples()` towards a 1,000 ESS target. That target was set on the
~690 raw hierarchical parameters, whose minimum ESS was 76–96 on every fold, so
it was never reachable and the cap always bound. The loop has been removed and
the samples are asked for directly; the two are statistically equivalent, since
`extra_samples()` continues the same chains without re-adapting.

## 4. What is being changed, and why

The review of PR #12 found one blocking defect and one set of experiments that
measures something other than what it claims to. Both force refits; the rest of
the review is scoring and reporting, and has been done without refitting.

### 4.1 The temporal forecasting fold leaked post-horizon data — refit

The training set was "every record whose year is not 2020–2022". The data run to
2024 while the covariates stop in 2022, so 888 records from 2023 and 2024 stayed
in training: the model was fitted on both sides of the window it was asked to
forecast. 354 of the 1,461 held-out assays (24%) sat at pixels that also carried
post-horizon training data, and on those pixels the dynamical model's MSE was
0.060 against 0.086 elsewhere.

Fixed in `validation_folds.R` (`year_start < min(validation_years)`, with a
`stopifnot`). The old fit is in `outputs/cv_draws_leaky_forecast/` with a note;
its scores are in the git history of `outputs/cv_summary.csv` at `6b1ba1b`. The
nulls have been regenerated on the corrected split already, since they need no
MCMC; the dynamical model waits for the refit.

### 4.2 Leave-one-country-out confounds spatial skill with an unidentified initial condition — retained, reported with the caveat

A held-out country's `init_country_raw` has no data, so its initial resistant
fraction reverts to the region prior, and that error is amplified through 15–29
years of deterministic selection before the comparison year. The review argued
this from the model structure; it is now confirmed empirically. Across the six
folds, the held-out bias correlates with that country's fitted country effect in
the full-data fit at **r = −0.94**, and excess MSE against the magnitude of that
effect at **r = +0.91**. Côte d'Ivoire and Ethiopia have the two largest
negative country effects (−3.4 and −4.7 on the logit initial-fraction scale) and
are the two worst-predicted folds, both overpredicting mortality by 0.27–0.36.

These folds are **not refitted**. The Côte d'Ivoire fold's worst Rhat of 2.33
raised the possibility that the negative result was a convergence artefact; the
per-fold breakdown rules that out — five of six countries have the dynamical
model worse than the insecticide mean, the ordering tracks the country effect
rather than Rhat, and Côte d'Ivoire is in fact the one country where the
dynamical model wins. So the result is structural and the fits are informative
as they stand. They are reported as what they measure: the difficulty of
predicting an entirely unsampled country, which is not a situation the deployed
model faces, since there are bioassays in every country.

### 4.3 A sub-national spatial block design replaces the national extrapolation concept — three new fits

Holding out blocks of cells within countries keeps every country intercept
identified, so the test isolates spatial prediction rather than compounding it
with the initial condition. See §5.

### 4.4 Scoring and reporting changes, all made without refitting

- **The PIT column was the mid-P value, not a randomised PIT.** `rowMeans()` over
  100 randomisations converges to `cdf_below + 0.5 * pmf_at`, which is not
  uniform under calibration for discrete data. Now `pit[, 1]`. The headline
  uniformity statistics used the full matrix and were never affected; the
  figures were. `check_validation_functions.R` now demonstrates the distortion
  on data with a large atom at 100% mortality: coverage 0.973 at nominal 0.95,
  CvM 10.58 against 0.03.
- **Skill is anchored on the intercept null**, not the nearest neighbour, which
  pinned an informative baseline at zero by construction. `excess` (MSE above
  the noise floor, in absolute mortality² units) and `rms_p` (its square root)
  are reported alongside every ratio.
- **Every model is scored at the external replicate-based overdispersion.**
  Letting each model fit its own made coverage a comparison of dispersion rather
  than of prediction: the intercept null reached 0.96 coverage by inflating rho
  to 0.45 against an external estimate of 0.12–0.22. What each model's own
  residuals imply is kept as a diagnostic in `cv_rho_comparison.csv`.
  `score_at_external_rho <- FALSE` in `validation_metrics.R` recovers the old
  behaviour as a sensitivity. The caveat to state in the paper: the dynamical
  model's posterior on p was fitted jointly with its own rho ≈ 0.3, so the swap
  is not perfectly clean without a refit.
- **The metric surface is smaller.** Coverage, mean PIT and Cramér–von Mises are
  three functionals of one PIT distribution, which is enough; the
  Kolmogorov–Smirnov statistic and the PIT ECDF figure are gone. The CvM null
  band is gone from the reported tables, because it assumes independent PIT
  values and held-out records share a posterior, so it is too narrow — CvM is an
  ordering, not a test. The WHO-threshold block scored the sample quantity
  rather than the population quantity and averaged predictive quantiles across
  folds; gone. The pixel-year aggregation rung averaged 1.3 assays per group, so
  it was the unpooled comparison under another name; only the country-year rung
  (17–18 assays per group) is kept.
- **Per-fold and per-lead-year breakdowns are restored** (`cv_by_fold.csv`,
  `cv_by_year.csv`), which master reported and the first version of this
  pipeline dropped.
- **Skill is reported against distance and data volume**, in
  `R/validation_geometry.R`. For every held-out record: km to the nearest
  training record of the same insecticide, number of such records within 100 km,
  and years since the last observation at that pixel. This is what makes the
  arbitrary geometry of the folds matter less, and it is the answer to the
  question a user of the map actually has.

## 5. The sub-national spatial block design

**Construction** (`R/validation_blocks.R`). Within each country, its
data-bearing cells are partitioned into three contiguous blocks carrying about a
third of that country's records each. Two families of cut are tried:

- **slabs** — project the cell centroids onto a direction and cut across it at
  the record-count terciles, giving three bands each spanning the country's full
  width in the perpendicular direction;
- **sectors** — take the bearing from the country's record-weighted centroid and
  cut on that, giving three wedges meeting at the centre.

Both are contiguous and both are balanced in records by construction. Thirty-six
orientations of each are scored, and the winner is the cut that maximises the
**distance from held-out cells to the nearest cell in another block** — the
record-weighted 25th percentile of it, so that one unlucky pair either side of a
boundary cannot decide the cut.

Scoring on separation rather than on block area matters, and was learned by
getting it wrong first. Maximising the smallest block's area does not work:
rotating a cut barely changes the areas, since they are three thirds of the same
country however it is sliced, so the objective is nearly flat across
orientations and ends up choosing on density noise. What it cannot see is block
*shape*, and the cuts it picked were thin slices that put a fifth of the
held-out records within 15 km of the next block. Scoring separation directly
also picks the axis sensibly without being told to: cutting a long country
across its length gives three roughly square blocks, while cutting along its
length gives three ribbons, and the ribbons score far worse.

Sectors win in 14 countries and slabs in 20, so trying both was worth it.

No buffer is applied. The blocks are large enough that a few cells near a
boundary cannot carry the result, and a buffer would remove training records
from exactly the countries whose intercepts this design exists to keep
identified.

Countries with fewer than 12 data-bearing cells or 40 records stay wholly in the
training set for every fold: CAR, Comoros, Djibouti, Equatorial Guinea, Eritrea,
Eswatini, Gabon, Guinea-Bissau, Mauritania, Mayotte, Sao Tome & Principe and
South Sudan, 406 records in total. Thirty-four countries are split.

**What the folds look like.** 9,094 / 8,907 / 8,954 held-out assays, every
blocked record held out exactly once, no pixel in both sets of any fold. Worst
record-share imbalance 0.09 from a third (South Africa and Togo); 28 of 34
countries under 0.03. Bioassay-weighted distance from a held-out pixel to the
nearest training pixel, pooled over the folds:

| | min | 10% | 25% | median | 75% | 90% | max |
|---|---|---|---|---|---|---|---|
| block folds | 4 | 18 | 35 | **76** | 131 | 216 | 766 |
| interpolation fold, for comparison | 14 | 21 | 25 | 36 | 48 | 66 | 466 |

So this reaches roughly seven times further than the interpolation experiment
could, and spans a range wide enough for the distance-stratified reporting of
§4.4 to say something.

**The residual limitation, which cannot be designed away.** 9% of held-out
assays still sit within 15 km of a training pixel, concentrated where a country
holds most of its records in one tight cluster: South Africa (within-country
separation 6 km at the 25th percentile), Togo (10 km), Kenya (15 km — the Lake
Victoria cluster), Cameroon (19 km), Benin (25 km). No contiguous partition
carrying a third of the records each can separate a cluster that is itself more
than a third of the records. The distance-stratified reporting is what handles
it: those records are not discarded, they are read at the distance they actually
represent.

**Sampling settings.** Unchanged from §3. Pinning each country's initial
condition in every fold removes the parameter that was previously unidentified,
so convergence should be better than the national folds, not worse. Check Rhat
on the first completed fold before launching the other two. If it has not
improved, a longer warmup is available for these folds, because they are a new
experiment and do not have to match the sampling settings of the retained
national folds.

## 6. Change-based scoring for the forecasting experiment

With the leak fixed, the experiment is still substantially a spatial test: most
held-out site-years are at sites with training data one to three years earlier,
so a local method gets the level nearly free. Scoring the *change* differences
the site level out and leaves the local slope, which is what the model claims to
know.

For each group *g* — a (cell, insecticide) or (country, insecticide) pair with
data in both windows — with before window *B* (the three years preceding the
cut) and holdout window *H*:

```
delta_obs(g)  = sum_H died / sum_H tested  -  sum_B died / sum_B tested
delta_pred(g) = mean over draws of ( weighted mean of p over H
                                     - weighted mean of p over B )
```

The target is a difference of two empirical proportions, with no model
assumptions in it. The floor is the sum of the two windows' irreducible
variances, which is what `noise_floor_var_pooled()` computes (added to
`validation_functions.R`, with checks: it reduces to `noise_floor_mse()` for a
single assay, and recovers the variance of a pooled proportion to within 2% over
3,000 simulated groups). So the existing MSE-minus-floor framework applies
unchanged, anchored at "no change" — which is what the nearest neighbour null
predicts by construction, so it needs no separate treatment.

Report a sign test alongside: the proportion of groups where the model gets the
direction of change right, among groups whose |delta_obs| exceeds its own noise
standard deviation.

Run it at both scales. Cell-insecticide is the most local and the thinnest, and
the floor handles that honestly; country-insecticide has 17–18 assays per group,
so the floor falls roughly seventeen-fold and the comparison is almost purely
about the population fraction.

**What this requires at fit time.** Predictions at the before-window cell-years,
which means a second index into `dynamic_cells$all_states` and one more term in
the `calculate()` call in `fit_validation_fold.R`. Because sampling cannot be
resumed across sessions (§8), this has to be in place before the refit starts —
it cannot be added to a finished fold. The before-window records are defined in
`validation_folds.R` alongside the training and test sets.

A single temporal origin means the conclusion rests on one realisation of the
recent trend. A rolling origin would cost another full run and is not worth it;
the limitation should be stated.

## 7. Cost

| fits | wall clock, two at a time |
|---|---|
| temporal forecasting, corrected split | ~62 h |
| three sub-national block folds | ~124 h |
| **total** | **~124 h, about five days** |

Per-fold cost is set by the dynamics graph, which is solved over all unique
cells × years × types regardless of how many rows enter the likelihood, so a
smaller K does not make each fold cheaper — it runs fewer of them. Four fits two
at a time is about 124 h against 250 h for the September run.

Running all four concurrently at two threads each is worth testing on one fold
first, given how weakly threads scale; memory is the constraint, since each
process holds the full dynamics graph.

The seven retained folds in `outputs/cv_draws/` are not refitted. They can be
slimmed offline — `rho_draws` is stored expanded to `n_draws × n_test` when only
`n_draws × n_classes` is distinct, about 1.8 GB across the folds — but the
scoring path already accepts either layout, and rewriting irreplaceable 62 h
files for disk space is not obviously worth the risk.

## 8. What survives a session, and what does not

Established by experiment, not assumption:

| operation on a reloaded `draws` object | result |
|---|---|
| `calculate(target, values = draws)` | **works** — targets are recoverable from `attr(draws, "model_info")` and the graph is re-traced |
| `extra_samples(draws, n_samples = ...)` | **fails** — `"object is from previous session and is now invalid"` |

The sampler state is bound to the session that created it, and redefining the
model produces new nodes the draws cannot attach to.

Consequences for planning:

- **A longer run must be requested up front**, through `warmup` and `n_samples`.
  Sampling cannot be topped up afterwards.
- **Folds to be compared must share their sampling settings.** Refitting one
  fold with longer warmup means refitting all of them.
- **Each saved fold therefore keeps `draws` and `prediction_arrays`.** Nothing in
  the committed pipeline reads them back, but they are what lets a finished fold
  produce a new prediction target — new cell-years, new aggregations, the
  before-window predictions of §6 — without re-running 62 h of MCMC. They are
  also most of each file's size. The four defunct folds in
  `outputs/cv_draws_defunct/` have no `draws` object and so cannot be used with
  greta's prediction interface at all, which is what that costs.

## 9. Pitfalls already hit, worth not repeating

- `calculate(..., nsim = n)` returns an **independent resample** of the
  posterior. It preserves joint structure across quantities, so it is valid for
  prediction, but it destroys MCMC ordering, so effective sample size cannot be
  recovered from it. Use `calculate(values = draws)` without `nsim`.
- `coda::effectiveSize()` on 20 draws returns roughly 20. Any ESS computed from
  a short run is meaningless; several tuning conclusions were drawn from such
  numbers and had to be withdrawn.
- `future.callr` buffers worker stdout until the future resolves, so a long run
  under it is invisible. Hence one process per fold.
- Functions called from a worker must have every dependency passed explicitly.
  `codetools::findGlobals(fit_fold, merge = FALSE)$variables` catches free
  variables; `n_unique_cells` was missing this way and cost a run.
- Verify that string edits to these scripts actually applied. A silently
  non-matching replacement cost a second run.
- Arm a monitor on the logs and check that the monitor itself is alive. A run
  died and went unnoticed for 14 h because the watcher had exited days earlier.
- `pgrep -f "validation_metrics.R"` matches the shell waiting on it as well as
  the R process. A completed run was reported as still running for six days on
  the strength of that.

## 10. Superseded artefacts

- `outputs/cv_draws_defunct/` — four folds from the earlier 2-chain run, no
  `draws` object, Kenya did not converge. Delete once the rebuild is complete.
- `outputs/cv_draws_leaky_forecast/` — the forecasting fold fitted on the leaky
  split (§4.1). Keep until the corrected fold is reported, as the record of what
  the leak was worth.
- `predictive_validation.R` and `dynamic_predictive_validation.R` still hold the
  plug-in deviance path and, in the latter, three copies of the model
  definition. Left in place deliberately, pending a decision on whether to
  report the plug-in results alongside; deleting
  `dynamic_predictive_validation.R` is what makes "three copies became one"
  true.
