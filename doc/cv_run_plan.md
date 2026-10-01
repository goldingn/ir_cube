# Cross-validation: design, run recipe, and results

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
Rscript R/fig_illustrate_bioassay_variability.R  # overdispersion per type; everything downstream needs it
Rscript R/run_validation_folds.R                 # nulls, then dispatch the model folds
Rscript R/validation_metrics.R                   # score everything on disk
Rscript R/validation_change.R                    # score predicted change, per forecast origin
Rscript R/validation_geometry.R                  # fold separation, and the leak check
Rscript R/variance_explained.R                   # variance explained and the noise ceiling
Rscript R/fig_variance_explained.R               # the bar figures
Rscript R/fig_predictive_validation.R            # figures and the table
```

The overdispersion fit comes first: `rho_lookup()` stops rather than falling
back, so nothing scores until `outputs/bioassay_rho_hierarchical.csv` and
`outputs/bioassay_rho_type_draws.rds` exist.

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
Rscript R/run_one_fold.R spatial_blocks 1 2 4 5 5
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

### 4.1 The temporal forecasting design: a leak, then a window with no power — rebuilt as a rolling origin

Two separate faults, found in that order.

**The leak.** The training set was "every record whose year is not 2020–2022".
The data run to 2024 while the covariates stop in 2022, so 888 records from 2023
and 2024 stayed in training: the model was fitted on both sides of the window it
was asked to forecast. 354 of the 1,461 held-out assays (24%) sat at pixels that
also carried post-horizon training data, and on those pixels the dynamical
model's MSE was 0.060 against 0.086 elsewhere. Fixed in `validation_folds.R`
(`year_start < cut_year`, with a `stopifnot`). The old fit is in
`outputs/cv_draws_leaky_forecast/` with a note; its scores are in the git
history of `outputs/cv_summary.csv` at `6b1ba1b`.

**The window.** The corrected 2020 three-year holdout turned out to be the worst
window in the record on both counts that matter. Sliding a three-year holdout
against the three years before it, the observed change at pixels assayed in both
is decisively negative at every origin from 2005 to 2017 and flattens only at
2018–2020: the 2020 origin reads +0.011 [−0.017, +0.040]. It is also the
thinnest, at 280 paired (pixel, insecticide) groups and 1,127 assays against
1,436 and 6,346 at a 2013 origin — enough to bound the mean change, nowhere near
enough to score any single site, which is why the direction test came out at
49%, a coin flip, and must not be reported as a finding.

The model is not misbehaving on the trend. It predicts −0.093 over that window,
and the slope fitted to the training years 2012–2019 implies −0.075 to −0.082.
It extrapolates the historical rate faithfully; the rate stopped.

**Five-year windows fix it.** At every origin they give 30–50% more paired
pixels and a signal roughly 5/3 larger, because the gap between window midpoints
is the window length, and a five-year window spans the 2018–2020 pause as well
as the decline either side of it. The 2018 origin's signal-to-resolvable ratio
goes from 0.2 at three years to 3.8 at five. The power analysis measures
this; the script and its figures are in the stub branch (§6).

**Two origins, fitted:**

| cut | training | % of data | holdout | paired pixels | observed change | ratio to resolvable |
|---|---|---|---|---|---|---|
| 2014 | 1995–2013, 14,285 | 52% | 2014–2018, 9,922 assays | 1,748 | −0.084 | 12.0 |
| 2018 | 1995–2017, 22,377 | 82% | 2018–2022, 4,096 assays | 881 | −0.042 | 3.8 |

Non-overlapping holdouts bar the shared endpoint, training at half and
four-fifths of the data, and a factor of two between their true rates of
decline. That contrast is the test: does the model track a slowing rate, or
carry a fixed slope forward? Earlier origins have a stronger signal still — a
2010 origin is the strongest in the record — but train on 17% of the data, so
they are not the model being deployed and their skill would not transfer.
Covariates end in 2022, so 2018 is the latest feasible five-year origin.

The 2020 three-year fold is deleted, not reported. Its training set was drawn
with `year_start <= cut`, so the cut year appeared in both training and holdout,
and its window landed on the one pause in twenty years of decline and was also
the thinnest. Its draws are parked in `outputs/cv_draws_defunct/`.

**The two origins are not independent.** The 2018 fold's before-window, 2013 to
2017, sits inside the 2014 fold's holdout, 2014 to 2018, and 2018 itself is in
both holdouts. Pooling them is still the right way to report a single forecasting
bar - the alternative is two bars whose difference is mostly window difficulty -
but the pooled interval carries roughly one and a bit folds' worth of
information rather than two, and should not be read as though the origins were
replicates. `outputs/cv_variance_explained_by_fold.csv` has them separately.

**How it is coded.** `validation_folds.R` exposes `forecasting_fold(cut_year,
window)` and a named list `temporal_forecasting_folds`, holding 2014 and 2018.
The single-fold alias `temporal_forecasting` is gone. Dispatch is `Rscript R/run_one_fold.R
temporal_forecasting 2014 4 4`. Each origin is stored under
`experiment = "temporal_forecasting_<cut>"` so that scoring never pools two
holdout windows whose true rates of decline differ by a factor of two; the file
name keeps the plain experiment name, so the origins sit together in
`outputs/cv_draws/`.

**Two changes the five-year window forced:**

- *The nearest neighbour null's lookback.* It searches the same year and
  `n_years_prior` earlier ones, intersected with training. With a five-year
  window and `n_years_prior = 3`, held-out years at lead 4 and 5 have no valid
  training year at all — and the old code then took `sort(...)[n]` of a vector
  of `Inf`, giving a threshold of `Inf`, and selected the *entire* training set:
  a global mean wearing a nearest-neighbour label, silently. `n_years_prior` is
  now the window length, and `predict_null_fixed_nn_counts()` stops rather than
  falling back. The three-year fold escaped this by one year, so its result is
  unaffected. The nearest-neighbour null is no longer tuned at all: it is
  reported at one neighbour, as a practice baseline, and separately at its best
  k on the held-out records, as an oracle bound.
- *Prediction volume.* The 2014 fold asks for 20,522 predictions against the
  three-year fold's 5,724. `fit_fold()` now thins the stored draws to 2,000 —
  which is what the scoring and the change score thin to anyway, and ESS is
  still measured on the unthinned ordered draws — and takes the before-window
  predictions in a second `calculate()` call, so peak memory tracks the larger
  window rather than both at once. `calculate(values = draws)` with no `nsim` is
  a deterministic function of the draws, so the pairing the change score needs
  survives the split.

### 4.2 Leave-one-country-out confounds spatial skill with an unidentified initial condition — removed

A held-out country's `init_country_raw` has no data, so its initial resistant
fraction reverts to the region prior, and that error is amplified through 15–29
years of deterministic selection before the comparison year. The review argued
this from the model structure; it was confirmed empirically before the folds
were dropped. Across the six folds the held-out bias correlated with that
country's fitted country effect at **r = −0.94**, and excess MSE against the
magnitude of that effect at **r = +0.91**. Côte d'Ivoire and Ethiopia have the
two largest negative country effects and were the two worst-predicted folds.

So these folds measure the difficulty of predicting an entirely unsampled
country, which is not a situation the deployed model faces — there are bioassays
in every country. The fold definitions and every code path that scored them have
been deleted, and the draws are parked in `outputs/cv_draws_defunct/`. The
sub-national blocks of §5 test spatial prediction without the confound.

### 4.3 A sub-national spatial block design replaces the national extrapolation concept — two new fits

Holding out blocks of cells within countries keeps every country intercept
identified, so the test isolates spatial prediction rather than compounding it
with the initial condition. See §5.

### 4.4 Scoring and reporting changes, all made without refitting

**One noise floor, one definition of variance explained.** The floor is
`noise_floor_mse()`: `y(1-y) k / (1-k)`, which since `E[y(1-y)] = p(1-p)(1-k)`
is exactly unbiased for a subset's mean sampling variance and assumes nothing
about how p is distributed. An empirical Bayes floor fitted per insecticide type
was tried and dropped — it was 1.4 to 9.6% high depending on rho, and is the
more fragile of the two where p is bimodal, as it is for DDT and
Alpha-cypermethrin. Variance explained is `1 - MSE/Var(y)`, floor-free, with the
noise share shown alongside as a band rather than divided out; the
intercept-null-referenced skill that also carried that name is gone, with the
second pixel-cluster bootstrap that supported it.

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
  to 0.45 against an external estimate of 0.12–0.22. Fitting a rho per null model
  is gone entirely — the nulls earn their place through mean squared error and
  variance explained, which need point predictions only, and the dynamical
  model's calibration is judged against held-out data directly. What the
  dynamical model's own residuals imply is still reported, in
  `cv_rho_comparison.csv` and its figure: a fitted 0.28 against an external 0.16
  is the model treating its own misfit as bioassay noise. The caveat to state in
  the paper: that posterior on p was fitted jointly with that rho, so scoring at
  the external value is not perfectly clean without a refit.
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

**Which countries are split.** Only the six the national folds hold out: Côte
d'Ivoire, Ethiopia, Kenya, Nigeria, Senegal and Tanzania. `validation_folds.R`
selected those as the countries with more than ten bioassays of every one of the
nine insecticide types since 2010, and that criterion matters here for the same
reason. Keeping to them buys three things:

- the block experiment becomes directly comparable with the national one,
  differing in exactly the intended respect and no other — same countries, same
  insecticide coverage, and only whether each country's initial condition is
  identified, which is the whole mechanism under test;
- the held-out set keeps a controlled mix of insecticides, rather than becoming
  a weighted average over whatever the thinner countries happen to hold;
- each fold trains on 88% of the data rather than 67%, much closer to the
  production fit, so the result transfers to the deployed model more directly.

Splitting all 34 data-bearing countries would cost the same three fits and buy
precision that is not the binding constraint: the model differences are already
determined to a standard error of 0.003 to 0.006 in excess mean squared error.
The other 40 countries (17,588 records) stay wholly in training in every fold,
as do the other two blocks of each split country, so every country intercept is
identified throughout. That is the entire purpose of blocking rather than
holding out countries.

**Construction** (`R/validation_blocks.R`). Within each split country, its
cells are partitioned into two contiguous blocks that are **large and even in
area**, subject to a floor on the bioassays each one carries. Two families of
cut are tried — slabs (project the cells onto a direction and cut across it) and
sectors (cut on the bearing from the country's area centroid) — at 36
orientations each, and the winner maximises the smallest block's area. Slabs
win in all six countries.

Area means the number of land cells of the mask falling in the block, not the
convex hull of its data-bearing cells: the cut is applied to every land cell of
the country, so the three areas are directly comparable and sum to the country.
Because the three sum to a constant, maximising the smallest block's area is
exactly what "large and even" means — a cut that makes one block small
necessarily makes another large.

This is the third objective tried, and the reasoning behind it is worth keeping.

- *Balance the records, maximise the smallest block's area.* Does not work.
  Rotating a cut barely changes the areas, so the objective is nearly flat
  across orientations and chooses on density noise; the cuts it picked were thin
  slices, with a fifth of held-out records within 15 km of training data.
- *Balance the records, maximise separation directly.* Better, but bounded by
  the record constraint. Bioassay effort is wildly uneven in space, so cutting
  at the record terciles puts a boundary straight through the densest cluster,
  and the block holding a third of the records occupies a small area. Short
  separations then occur exactly where most of the data is.
- *Even areas with a record floor.* Lets a dense cluster sit whole inside one
  block, which keeps both blocks large and lengthens the separation for the
  sparser one. The cost is an uneven split of records between folds, which
  costs nothing overall: every record is still held out exactly once, so only the
  per-fold counts differ.

Floors, so that favouring area cannot leave a fold too thin to score: each block
must carry at least 150 bioassays or 12% of its country's records, whichever is
larger, and at least 10 data-bearing cells. The floor binds only in Kenya.

No buffer is applied. The blocks are large enough that a few cells near a
boundary cannot carry the result, and a buffer would remove training records
from exactly the countries whose intercepts this design exists to keep
identified.

**What the folds look like.** 5,382 / 3,312 held-out assays — 8,694 in total,
over the same six countries and the same "2010 and later" test window. All nine insecticide types
appear in every fold. Every blocked record from 2010 on is held out exactly
once; records at a held-out cell from before 2010 are dropped from the
experiment rather than returned to training, which is what the national folds
also do with pre-2010 records of a held-out country. No pixel is in both sets of
any fold.

Per country, with the smallest block's area, how even the three areas are, and
the worst departure of a block's record share from a third:

| country | family | smallest block | area evenness | worst record share |
|---|---|---|---|---|
| Senegal | slab | 68,900 km² | 1.00 | 0.19 |
| Côte d'Ivoire | slab | 109,300 km² | 1.00 | 0.12 |
| Kenya | slab | 130,600 km² | 0.43 | 0.41 |
| Nigeria | slab | 309,500 km² | 1.00 | 0.15 |
| Tanzania | slab | 319,100 km² | 1.00 | 0.19 |
| Ethiopia | slab | 383,700 km² | 1.00 | 0.16 |

Bioassay-weighted distance from a held-out pixel to the nearest training pixel,
pooled over the folds:

| | min | 10% | 25% | median | 75% | 90% | max |
|---|---|---|---|---|---|---|---|
| block folds | 5 | 19 | 46 | **79** | 132 | 178 | 375 |
| interpolation fold, for comparison | 14 | 21 | 25 | 36 | 48 | 66 | 466 |

7% of held-out assays sit within 15 km of training data and 14% within 25 km,
against 18% and 27% under the first objective. Separation is capped at a few
hundred km rather than the ~800 km the national folds reach, because a held-out
block is surrounded by training data in neighbouring countries. That is the
price of keeping the country intercepts identified, and it is the right price.

Per-country median separation: Nigeria 132 km, Ethiopia 98, Tanzania 93, Côte
d'Ivoire 83, Senegal 58, **Kenya 53**. Kenya is still the weakest, because it
holds the majority of its records in one cluster west of Lake Victoria; the area
objective at least keeps that cluster whole in one block, which lifted Kenya's
median separation from 26 km to 53 and halved the share within 15 km. The
distance-stratified reporting of §4.4 handles what remains: those records are
not discarded, they are read at the distance they actually represent.

One cosmetic note: a border cell is assigned to whichever country holds most of
its records, so fold 3's held-out set touches a seventh country through one such
cell.

**Sampling settings.** Unchanged from §3. Pinning each country's initial
condition in every fold removes the parameter that was previously unidentified,
so convergence should be better than the national folds, not worse. Check Rhat
on the first completed fold before launching the other two. If it has not
improved, a longer warmup is available for these folds, because they are a new
experiment and do not have to match the sampling settings of the retained
national folds.

## 6. Change-based scoring, and two analyses that are not here

`validation_change.R` scores predicted change between the before-window and the
holdout, per forecast origin, which is the quantity the forecasting experiment is
actually about. The before-window predictions it needs are saved with each fold
by `fit_validation_fold.R`; because sampling cannot be resumed across sessions
(§8), that has to be in place before a fit starts and cannot be added afterwards.

Two supporting analyses were written and are not in this branch, having served
their purpose:

- the power analysis over window length and origin, which chose the five-year
  window and showed 2014 to be the strongest feasible origin;
- skill against separation from the training data, which was a negative result:
  the joint test gave p = 0.058 for the dynamical model and 0.264 for the
  nearest survey. Distance has an intraclass correlation of 0.91 by pixel, so a
  pixel random effect absorbs the identifying contrast rather than controlling
  for it, and binned weighted least squares with a pixel-cluster bootstrap is
  what the design supports. The raw distance bins that preceded it were
  confounded in exactly that way and are also gone.

Both live in the stub branch if they are wanted again.

## 7. Cost

| fits | wall clock, two at a time |
|---|---|
| ~~temporal forecasting, corrected split~~ | done, 62 h |
| ~~two sub-national block folds~~ | done |
| forecast origins 2014 and 2018, five-year windows | ~62 h |

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

Parked in `outputs/cv_draws_defunct/`, not read by anything, kept as the record
of what was tried:

- the six leave-one-country-out folds (§4.2);
- the three-year 2020 forecasting fold, whose training set was drawn with
  `year_start <= cut`, so the cut year appeared in both training and holdout,
  and whose window landed on the one pause in twenty years of decline;
- `cv_draws_leaky_forecast/`, the same fold before the leak was fixed;
- four folds from the earlier two-chain run, with no `draws` object.

`validation_folds.R` defines `validation_experiments`, and both
`run_one_fold.R` and `validation_metrics.R` refuse anything outside it, so a
stray fold cannot be scored by accident.

Deleted on this branch: `predictive_validation.R` and
`dynamic_predictive_validation.R`, which held the plug-in deviance path and
three copies of the model definition; and `validation_metric_eval.R`, the
simulation study that supported the choice of metric.

## 11. What the validation found, and what follows

The dynamical model beats the nearest-recent-survey baseline on sub-national
spatial extrapolation (+8.8 [+3.3, +14.7] percentage points of variance
explained, paired within bootstrap replicates) but not on forecasting
(−3.6 [−8.5, +1.5]). Diagnosis, from the before-window predictions saved with
each fold:

- the model has the direction of local change right — 77–83% sign agreement at
  the 2014 origin among groups whose observed change exceeds its own noise — but
  predicts two to three times too much decline, and the signed error grows with
  forecast horizon;
- selection is linear in each covariate, and the net-use response saturates in
  the raw data (#23);
- resistance can only increase in the model: there is no fitness cost or decay
  term (#24).

Both fixes need a refit and are out of scope here.
