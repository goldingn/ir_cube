# Plan: moving cross-validation from point prediction to posterior predictive assessment

## 1. Why change

The current cross-validation scores the model on its ability to predict the value of an
individual held-out bioassay. But an individual bioassay is a noisy, overdispersed
measurement (beta-binomial, `rho` ~ 0.1) of the quantity the model actually targets:
the susceptibility fraction of the whole mosquito population at that place, time and
insecticide. Two consequences:

1. **The metrics have an unreachable floor.** Even a model that knows the true
   population fraction exactly incurs a large expected deviance, because the data are
   noisy around it. The absolute numbers in `figures/CV_deviance.png` are therefore
   uninterpretable, and genuine differences between models are compressed by a large
   constant.
2. **They answer the wrong question.** We do not care whether the model can guess that
   a particular cone test killed 63 of 100 mosquitoes. We care whether the distribution
   of bioassay results that the model *predicts* for a held-out country/year matches the
   distribution actually observed. That is a posterior predictive question.

The change is therefore: score the **full out-of-sample posterior predictive
distribution**, not a point summary of it, and report calibration and coverage
alongside accuracy.

## 2. Current state

### Fold definitions (`R/predictive_validation.R:1-395`)

Three experiments, test data restricted to `year_start >= 2010` throughout:

| Experiment | Split | Code |
|---|---|---|
| Spatial extrapolation | Leave-one-country-out over countries with >10 bioassays for each of the 9 insecticides since 2010 | `:92-183` |
| Spatial interpolation | k-means (50 centroids) on unique locations; within 5th-percentile centroid distance -> test, beyond a 1.5x buffer -> training, between -> excluded | `:186-247` |
| Temporal forecasting | Final 3 years of covariate coverage held out | `:359-390` |

### Baseline models (same file)

- **Intercept-only**: binomial GLM `cbind(died, n - died) ~ insecticide_type` (`:441-530`).
- **Nearest neighbour**: pooled `died`/`tested` over the k nearest training bioassays of
  the same insecticide within the same or previous n years, empirical-logit smoothed;
  k chosen by grid search on a 100-record internal holdout carved from training, so no
  test leakage (`:540-860`).

### Model folds (`R/dynamic_predictive_validation.R`)

Per fold: re-define the full greta dynamical model with the likelihood restricted to the
training rows, run HMC (8 chains, 2000 warmup, 1000 samples), then collapse to point
estimates at `:365-390`:

```r
predicted_test      <- calculate(population_mortality_vec_test, values = draws, nsim = 1e3)
predicted_test_mean <- apply(predicted_test[[1]], 2:3, mean) %>% as.numeric()
rho_classes_test_mean <- apply(rho_classes_test[[1]], 2:3, mean) %>% as.numeric()
```

The three experiments are three near-identical copy-pasted ~200-line blocks
(extrapolation `:204-440`, interpolation `:460-745`, forecasting `:750-985`).

### Metrics and presentation

- `betabinom_dev()`: -2 log-likelihood of a beta-binomial evaluated at the posterior
  **mean** p and rho, averaged per record.
- `bias`: `mean(predicted_mean - died / mosquito_number)`.
- Nulls scored at a **fixed `rho = 0.14`**; the dynamical model at its own fitted rho.
- `outputs/CV_result_for_plot.csv` -> `figures/CV_deviance.png`, `figures/CV_bias.png`:
  tile grids of model x (country | year), faceted by experiment, best value bolded.

### Related but disconnected code

- `R/validation_metric_eval.R`: simulation comparing deviance/RMSE/CRPS/KS/CvM under
  "truth"/"biased"/"noisy" scenarios. Already notes that CvM is Taggart's (2022) PS2 and
  has expectation 0 for a calibrated model.
- `R/fig_internal_validation.R`, `R/illustrate_validation.R`: DHARMa randomised quantile
  residuals -- but **in-sample only**, on the full fitted model.

### Defects to fix in passing

1. **Point collapse** (above) discards parameter uncertainty entirely.
2. **No noise floor**, so scores have no interpretable scale.
3. **Unfair dispersion comparison**: nulls at fixed rho vs model at fitted rho means a
   model can win on deviance by fitting rho better rather than predicting better.
4. **No out-of-sample calibration or coverage** anywhere.
5. `rmse()` at `predictive_validation.R:419` is `sqrt(mean(o - p)^2)` -- that is
   |mean error|, not RMSE. Currently unused in headline results.
6. Three copy-pasted model blocks; a change to the model must be made three times.

## 3. The central object

For each held-out bioassay *i* with `n_i` mosquitoes, the out-of-sample posterior
predictive distribution is a finite mixture over the S posterior draws from the
training-fold fit:

```
F_i(y) = (1/S) * sum_s  pbbinom(y; size = n_i, p = p_i^(s), rho = rho_{c(i)}^(s))
f_i(y) = (1/S) * sum_s  dbbinom(y; size = n_i, p = p_i^(s), rho = rho_{c(i)}^(s))
```

This is **analytic** -- `extraDistr::dbbinom`/`pbbinom` averaged over draws. No forward
simulation is needed for the pmf, cdf, quantiles, or entropy, which removes Monte Carlo
noise from the residuals and makes tie handling exact. Simulation is needed only for
CRPS and for aggregated/pooled checks.

Everything below is a functional of `F_i`. The only change needed upstream is to
**save the draws** of `population_mortality_vec_test` and `rho_classes` instead of
their means.

## 4. Metrics

### 4.1 Log score, elpd, and the KL divergence question

Per record, the log predictive density `log f_i(y_i)`; averaged over the test set this
is the **elpd**. This is the right home for the KL idea, with an important caveat:

- **Per-observation KL is not estimable.** We have one draw from the true distribution
  of bioassay *i*, so `KL(true_i || F_i)` cannot be computed pointwise.
- **Differences in elpd are exactly differences in KL.** Since
  `E[log f_true] - E[log f_model] = KL`, and `E[log f_true]` is common to all
  candidates, `elpd_A - elpd_B` estimates `KL_B - KL_A` with no bias. So model
  comparison on the KL scale is available directly, with a bootstrap or
  paired-difference standard error over test records.
- **An absolute scale via predictive entropy.** Because `y` is supported on
  `{0, ..., n_i}`, the entropy of the predictive distribution is exactly computable:
  `H_i = -sum_y f_i(y) log f_i(y)`. `-mean(H_i)` is the elpd a perfectly calibrated
  model would attain -- the noise floor. The gap `-mean(H_i) - elpd` is an estimate of
  the average KL divergence from the true data-generating distribution (exact when the
  true and predictive distributions have equal entropy; approximate otherwise). Report
  it as "divergence above the achievable minimum" with that caveat stated.
- **Skill score for readability**:
  `skill = (elpd_model - elpd_null) / (-mean(H) - elpd_null)`, so 0 = the nearest
  neighbour heuristic, 1 = the noise floor. This is the single number that makes the
  deviance results interpretable and is what should appear in the table.

### 4.2 Randomised quantile residuals / PIT

For discrete data, the randomised PIT (Dunn & Smyth 1996):

```
u_i = F_i(y_i - 1) + v_i * f_i(y_i),   v_i ~ U(0, 1)
z_i = qnorm(u_i)
```

Under a correct predictive distribution `u_i ~ U(0,1)` iid. Compute R = 100
randomisation replicates and show the band, so conclusions do not hinge on one draw of
`v`. This is what `DHARMa::createDHARMa` does in `fig_internal_validation.R:22-29`, but
computed analytically and, critically, **out of sample**.

Summaries of deviation from uniformity, each with a null band obtained by simulating
`u ~ U(0,1)` at the same sample size:

- **Cramer-von Mises** (`cvm_stat()` already in `validation_metric_eval.R:135-140`) --
  headline distributional-mismatch statistic. Equals Taggart's PS2, expectation 0 under
  calibration, and decomposes into over/under-prediction and over/under-dispersion
  components. The decomposition is the reason to prefer it over KS.
- **Kolmogorov-Smirnov** and **Anderson-Darling** as secondary (AD is more sensitive in
  the tails, which is where the overdispersion assumption is most likely to fail).

### 4.3 Coverage and calibration

- **Coverage curve**: for nominal levels 10%...95%, the empirical fraction of `y_i`
  falling inside the central posterior predictive interval, from mixture quantiles.
  Plotted against the diagonal. This is the most intuitive plot available for a
  non-statistical reader.
- **Headline number**: empirical coverage of the 95% posterior predictive interval.
- **PIT ECDF-difference plot** with simultaneous confidence bands
  (Saeilynoja, Buerkner & Vehtari 2022; `bayesplot::ppc_pit_ecdf`), overall and stratified
  by experiment, insecticide class, year, country, and predicted-mortality decile.
- **Weighted interval score / interval score**, decomposed into width +
  under-prediction penalty + over-prediction penalty. Gives sharpness and calibration in
  one number, in observation units.

### 4.4 Accuracy in interpretable units

- **CRPS** on the mortality-proportion scale, from posterior predictive draws
  (`scoringRules::crps_sample`; the existing `crps()` in
  `validation_metric_eval.R:47-92` needs adapting to draw p and rho from the posterior
  rather than fixing them). Reads as "on average the predictive distribution sits X
  mortality percentage points from the observation".
- **CRPS skill** relative to the nearest-neighbour null, same normalisation as elpd.

### 4.5 Bias, done properly

- Keep `mean(E[y_i]/n_i - y_i/n_i)` but compute it **per posterior draw**, giving a
  posterior distribution for the bias rather than a point.
- **Calibration in the large**: `mean(u_i) - 0.5`.
- **Reliability diagram**: bin test records by predicted mortality decile; plot mean
  observed proportion against mean predicted, with intervals. Because many bioassays are
  averaged within a bin, the bin mean estimates the *population-level* fraction with
  much less noise than any individual assay -- this plot is the direct answer to
  "does the model get the population quantity right", and should be a main-text panel.

### 4.6 Aggregated (population-level) validation

The most direct attack on the noise problem. Pool held-out records into groups
(country x insecticide x year, or spatial cluster x insecticide x year), and compare:

- observed pooled mortality `sum(died) / sum(n)`, against
- the posterior predictive for the pooled count (sum of independent beta-binomials
  given a draw of p and rho -- obtained by simulation, not analytically).

Report coverage and bias at increasing aggregation levels. The noise floor falls as
groups grow, so this shows how much of the residual scatter is irreducible assay noise
versus model error. Presenting the same metric at 1, ~5 and ~20 assays per group makes
the argument for the whole approach visible in one figure.

### 4.7 Out-of-sample posterior predictive checks

For each fold, compare observed test-set summary statistics against their distribution
under replicate datasets simulated from the posterior predictive:

- mean and SD of mortality,
- the full ECDF of mortality,
- fraction of bioassays below the WHO thresholds (<90% = confirmed resistance,
  90-98% = possible resistance).

The WHO-threshold version is policy-relevant and needs no statistical training to read:
"in held-out countries, the model predicted that 41% of bioassays would show confirmed
resistance; 39% did."

## 5. Fair comparison between models

All candidates must supply a **predictive distribution**, not a point prediction, or the
proper scoring rules are not comparable:

- **Intercept-only GLM**: beta-binomial with `rho` fitted by MLE on the training fold,
  and uncertainty in p from the GLM's sampling distribution (or a conjugate
  beta posterior on the pooled counts per insecticide type).
- **Nearest neighbour**: beta-binomial with p from the pooled neighbour counts and a
  beta posterior for p; `rho` fitted on the same internal holdout used to choose k.
- **Dynamical model**: the posterior predictive mixture defined in section 3.

Then report the fitted `rho` per model as a diagnostic in its own right. This removes
the current confound where the nulls are handicapped by a fixed `rho = 0.14`.

## 6. Presentation

### Main text: one figure, one table

**Figure -- "Out-of-sample predictive performance".** Columns = the three experiments
(interpolation, extrapolation, forecasting). Rows:

- (a) **Coverage curve** -- nominal vs empirical coverage, one line per model, diagonal
  reference. "When the model said there was a 95% chance, how often was it right?"
- (b) **Reliability diagram** -- predicted vs observed pooled mortality by decile.
  "When the model predicts 60% mortality, is the average outcome 60%?"
- (c) **PIT ECDF difference** with simultaneous 95% band. "Does the model reproduce the
  whole spread of results, not just the average?"

**Table -- one block per model, one row per experiment.** Columns:

| n test | 95% PPI coverage (%) | mean PIT | CRPS (mortality points) | elpd skill vs NN | CvM (null band) |
|---|---|---|---|---|---|

Give each column a plain-language subheading in the caption:

- bias / mean PIT -> "Are predictions right on average?"
- coverage -> "Are the uncertainty ranges the right width?"
- CRPS -> "How close are predictions to reality?"
- elpd skill -> "How much better than a simple spatial average? (0 = no better, 1 = as
  good as the noise allows)"
- CvM -> "Does the model reproduce the full spread of real results?"

Deliberately do **not** make elpd/KL the headline number in the abstract or main text
prose -- lead with coverage and the reliability diagram, and keep the information-theoretic
quantities in the table, normalised as skill scores.

### Supplement

- Retain the existing tile plots (`CV_deviance.png` layout is good), but with CRPS skill
  and 95% coverage replacing raw deviance, so the per-country and per-year breakdown
  survives the change.
- Aggregated-validation figure (section 4.6) at increasing group sizes.
- Out-of-sample PPC of the mortality ECDF and WHO-threshold fractions.
- PIT histograms stratified by insecticide class, region, year, predicted decile.
- The extended metric simulation study (section 7, step 7) justifying the choice of
  metrics and demonstrating the noise floor.
- A pedagogical panel: one true population fraction, the beta-binomial spread of
  bioassays it generates, and the model's posterior predictive for the same. Much of
  this exists already in `R/fig_illustrate_bioassay_variability.R`.

## 7. Implementation

### New: `R/validation_functions.R`

Shared, model-agnostic, unit-testable. All take `p_draws` (S x n) and `rho_draws`
(S x n or S x n_class + index):

- `ppd_pmf(y, n, p_draws, rho_draws)`, `ppd_cdf()`, `ppd_quantile()`, `ppd_interval()`
- `ppd_entropy(n, p_draws, rho_draws)` -- noise floor
- `ppd_pit(y, n, p_draws, rho_draws, n_rep = 100)` -- randomised PIT, returns matrix
- `log_score()`, `elpd()`, `crps_ppd()`, `interval_score()`, `wis()`
- `coverage_curve()`, `reliability_bins()`
- `cvm_stat()`, `ks_stat()`, `ad_stat()` (move `cvm_stat`/`ks_stat` here from
  `validation_metric_eval.R:112-140`), plus `pit_null_band()` by simulation
- `skill()` normalisation helper
- fixed `rmse()`

### New: `R/validation_folds.R`

Extract `predictive_validation.R:1-395` (data prep, `split_data()`, the three fold
definitions, the interpolation diagnostic plot) so that both the null-model and the
dynamical-model scripts source the fold definitions rather than the latter sourcing the
former in its entirety.

### Modify: `R/predictive_validation.R`

Reduce to the null models. Each returns a tidy tibble of test records plus a matched
`p_draws`/`rho_draws` array, saved to `outputs/cv_draws/`, so the nulls and the
dynamical model flow into an identical scoring path.

### Modify: `R/dynamic_predictive_validation.R`

1. Collapse the three copy-pasted model blocks into a single
   `fit_fold(train_df, test_df, ...)` returning the draws. This is the largest single
   maintainability win in the file and removes the risk of the blocks drifting apart.
2. Replace the `apply(..., mean)` collapse at `:379-387` with saving the raw draws:
   `outputs/cv_draws/<experiment>__<fold>.rds` holding `p_draws` (S x n_test),
   `rho_draws`, and the test-row key. Roughly 1000 x 500 x 8 bytes = 4 MB per fold,
   ~15 folds, ~60 MB total. Add `outputs/cv_draws/` to `.gitignore`.
3. Move all metric computation and plotting out of this script.

### New: `R/validation_metrics.R`

Load every `outputs/cv_draws/*.rds`, compute all metrics from
`R/validation_functions.R`, write `outputs/cv_scores.csv` (one row per model x
experiment x test record) and `outputs/cv_summary.csv` (one row per model x experiment x
stratum). All plotting reads these two files.

### New: `R/fig_predictive_validation.R`

The main-text figure and the supplementary panels of section 6. Retains the tile-plot
code from `dynamic_predictive_validation.R:1010-1150` for the supplement, re-pointed at
the new metrics.

### Extend: `R/validation_metric_eval.R`

The existing simulation has "truth", "biased" and "noisy" prediction scenarios. Add
**overdispersed** and **underdispersed** scenarios -- dispersion error is exactly what
the new approach is meant to catch and what the current deviance-only summary conflates
with location error. Also add the noise-floor demonstration (elpd of a model that knows
the true p, vs sample size), which is the evidence for the whole change of approach.

### Align: `R/fig_internal_validation.R`

Re-point the in-sample RQR calculation at `ppd_pit()` from
`R/validation_functions.R` instead of DHARMa, so in-sample and out-of-sample residuals
are computed identically and are directly comparable. The `Inf` clamping at `:36-38`
becomes unnecessary once the PIT is computed analytically rather than from 1e4
simulations.

## 8. Compute cost and ordering

The MCMC cost is unchanged -- the same fits, saving draws rather than means. But the
saved outputs from the last run hold only posterior means
(`outputs/dynamic_spatial_extrapolatoin_CV_pred.rds`, and even that is not currently in
`outputs/`), so **all folds must be re-run**: roughly 12 country folds + 1 interpolation
+ 1 forecasting fold, each an 8-chain HMC run.

Suggested order, so that the expensive step happens once and last:

1. `R/validation_functions.R` + unit checks against known cases (no MCMC).
2. Extended `R/validation_metric_eval.R` simulation -- confirms the metrics behave as
   intended on data where the truth is known, before touching the real analysis. Cheap,
   and it produces a supplementary figure.
3. `R/validation_folds.R` extraction; verify the folds are byte-identical to the current
   ones (seeded k-means at `predictive_validation.R:196`, so this is checkable).
4. Null models re-emitting predictive distributions; score them. Cheap, and exercises
   the whole scoring path end to end.
5. `fit_fold()` refactor of `dynamic_predictive_validation.R`; verify one fold
   reproduces the current posterior means before proceeding.
6. Re-run all folds, saving draws.
7. `R/validation_metrics.R` and `R/fig_predictive_validation.R`.

## 9. Decisions still open

1. **Whether to keep deviance at all.** Recommendation: drop the plug-in deviance
   entirely and replace it with elpd, which is its posterior-predictive generalisation.
   Keeping both invites the reviewer question "why do these disagree?".
2. **How to present the KL quantity.** Recommendation: as an elpd skill score bounded
   by the noise floor, with the raw elpd and its standard error in the supplement. The
   entropy-based absolute divergence is worth reporting but needs its equal-entropy
   caveat stated plainly.
3. **Aggregation grouping** for section 4.6 -- country x insecticide x year is the
   natural choice given the fold structure, but group sizes will be uneven. Consider
   fixed-size groups drawn at random within country x insecticide, repeated, so the
   noise floor is comparable across groups.
4. **Whether the nulls get uncertainty in p, or only in rho.** Giving them full
   predictive distributions is the fair comparison, but it makes them harder to describe
   in the paper. Recommendation: do it, and describe them in the supplement as
   "beta-binomial with a beta posterior on the mean".
5. Whether the Hancock et al. comparison (commented out at
   `predictive_validation.R:940-1120` and `dynamic_predictive_validation.R:421-451`) is
   being revived. It only supplies point predictions, so it can enter the reliability
   diagram, bias and CRPS-of-the-point comparisons, but not elpd, coverage or PIT.
