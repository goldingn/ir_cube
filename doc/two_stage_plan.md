# Two-stage correction: run sheet

Issue idem-lab/ir_cube#21. Stacked on the #12 branch (`posterior-predictive-validation`).

## Folds

Five outer folds, the same as #12, with no change to the fold design:

| experiment | fold | dynamical draws |
|---|---|---|
| spatial_interpolation | all | available |
| spatial_blocks | 1 | available |
| spatial_blocks | 2 | available |
| temporal_forecasting | 2014 | available |
| temporal_forecasting | 2018 | available |

Leave-one-country-out and the legacy 2020 forecasting fold are defunct and are not run.

## Steps per fold and insecticide type

1. **Stage-one predictions.** `R/dynamical_predictions.R` recomputes the dynamical-model draws at the training assays and at held-out assays, paired with the saved `p_draws`.
2. **Reference predictor.** `m_ref` is the posterior mean of the logit draws at the training assays. Option (see step 6): the approximate leave-out mean.
3. **Stage A fit** (`R/two_stage_correction.R`, `tmb/two_stage_correction.cpp`):
   - Response: empirical logit, with v inflated by the per-type rho from `outputs/bioassay_rho_hierarchical.csv`.
   - Two variants: `omega_u`, then `omega_xi_u`.
4. **Prediction.** Joint latent draws, with the cut-posterior mode shift for each dynamical draw. Draws go to `outputs/cv_draws_two_stage/<model>__<experiment>__<fold>.rds`, in the #12 format.
5. **Scoring.** Beta-binomial log score, PIT and coverage, using `R/validation_functions.R`. All models are scored at the per-type rho: the dynamical model, the nulls and both two-stage variants.
6. **Double-counting diagnostic (PSIS).** Pareto-smoothed importance sampling on the saved dynamical draws approximates leave-out predictions at the training assays, without refitting:
   - Leave out a whole pixel-year, since most assays share one.
   - Report the Pareto k values, and the shift of the leave-out logit mean from `m_ref` relative to the residual SD.
   - Also report the posterior SD of `m` relative to the residual SD. This bounds the bias from fitting the hyperparameters at `m_ref`.
   - **Potential solution:** if the shift is material and the stage-A intervals are overconfident on the held-out assays (PIT or coverage), refit stage A with the leave-out mean as `m_ref`. This is a cheap approximation to stacking; full K-fold stacking of stage one is unaffordable at 60–95 h per MCMC fit. Where k > 0.7 for a group, fall back to `m_ref` and count those groups.
   - Run on `spatial_interpolation` first.
7. **Diagnostics:**
   - Empirical-logit residuals at 0% and 100% mortality, by region and by bin of `m`. This decides whether stage B (PQL) is needed.
   - Skill gain by forecast horizon.
   - Large-scale structure in the fitted fields, and their correlation with covariates.

## Further steps

8. **Maps** (`R/two_stage_maps.R`). Fit `omega_xi_u` per type to all the data, on top of the full dynamical fit. Map, in the layout of the dynamical-model figures:
   - two-stage predicted mortality;
   - the second-stage correction (ω + ξ, logit scale);
   - the difference from the dynamical model (percentage points);
   - the correction SD.
9. **Mesh resolution.** Refine the ξ mesh (and the ω mesh) beyond the defaults, and keep the finer mesh if the cross-validation scores improve.

## Implementation notes

- **Meshes.** The ω mesh has a node at every site, merged within a cutoff that grows until there are at most 2500 nodes, with edges of at most 250 km between sites. ξ uses a separate, coarser mesh with at most 600 nodes. Putting ξ on the ω mesh gives a Cholesky factor with about 29M non-zeros, which is too slow at 2500 nodes × 30 years. Pass `mesh_xi = mesh` to use one mesh. The ξ mesh cutoff is about 190 km, so revisit its size if the fitted η range is short.
- **Compilation and BLAS.** The template must be compiled with TMBad; CppAD is about 100× slower. Runs use OpenBLAS via `LD_PRELOAD`, because the reference BLAS is about 10× slower for the Cholesky.
- **Simulation check** (`R/check_two_stage_correction.R`):
  - Hyperparameters are recovered.
  - The cut-posterior shift matches a refit to about 1e-15.
  - 95% coverage for `omega_xi_u` is 0.93 for interpolation and 0.91 for forecasting. The central intervals are slightly narrow, which is consistent with plugging in the hyperparameters.

## Stage-one results (`R/stage_one_effective_parameters.R`)

Effective number of parameters of the dynamical model on each fold's training assays. The nominal count is 698.

| fold | assays | pixel-years | pD (mean p) | pV | p_WAIC | p_loo |
|---|---|---|---|---|---|---|
| interpolation | 24318 | 19453 | 287 | 442 | 340 | 342 |
| blocks 1 | 21321 | 16861 | 285 | 480 | 338 | 340 |
| blocks 2 | 23628 | 18530 | 288 | 462 | 346 | 352 |
| forecasting 2014 | 14285 | 10863 | 222 | 397 | 267 | 270 |
| forecasting 2018 | 22377 | 17533 | 278 | 452 | 323 | 326 |

- Per type, p_WAIC and p_loo are 23–53 (7–50 on the forecasting folds).
- pD with the plug-in at the posterior mean of the *parameters* is unreliable: it is negative for Alpha-cypermethrin on interpolation. The averaged non-centred parameters, pushed through exp(beta) and the recursion, are far from the posterior mode. So use pD at the posterior mean of p, or p_loo.

Grouped PSIS, leaving out a pixel-year (and, in brackets, a whole pixel):

| fold | groups with k > 0.7 | residual SD | posterior SD / residual SD | RMS leave-out shift / residual SD | leave-out / in-sample residual SD |
|---|---|---|---|---|---|
| interpolation | 0.10% (0.32%) | 1.94 | 0.12 | 0.036 (0.052) | 1.010 |
| blocks 1 | 0.13% (0.49%) | 1.90 | 0.13 | 0.041 (0.059) | 1.012 |
| blocks 2 | 0.15% (0.55%) | 1.89 | 0.13 | 0.039 (0.059) | 1.011 |
| forecasting 2014 | 0.12% (0.33%) | 1.83 | 0.17 | 0.053 (0.062) | 1.01 |
| forecasting 2018 | 0.10% (0.33%) | 1.87 | 0.14 | 0.040 (0.056) | 1.01 |

- **Double counting is negligible.** Fitting stage A to in-sample residuals understates the residual SD by 1–4%.
- **The organophosphates are affected most, and still only slightly:** RMS shift 0.07–0.11 of the residual SD, and posterior SD 0.26–0.41 of the residual SD. Worst case: Pirimiphos-methyl on forecasting 2014 (151 training assays), RMS shift 0.20, posterior SD 0.80.
- **Decision:** keep `m_ref` as the posterior mean. `m_ref=loo` stays available in the runner as a sensitivity check.
- **Caveat on the saved draws:** 282 of the interpolation fold's 2000 paired draws are exact repeats (81 in blocks 1, 311 in blocks 2, 47 in forecasting 2014, 107 in forecasting 2018), where HMC stuck. PSIS is run on the distinct draws only; the repeats give spurious k = Inf.

## Term inclusion

Keep `xi` only if `omega_xi_u` improves the interpolation or forecasting scores over `omega_u`.

## Mesh resolution results

`omega_xi_u` only, `m_ref` = posterior mean, on interpolation and blocks 1+2 (step 9). Run with `R/run_two_stage_folds.R <experiment> <fold> mesh=<tag> variants=omega_xi_u`; the model is saved as `two_stage_omega_xi_u_mesh-<tag>`. Compared by `R/two_stage_mesh_comparison.R` (writes `outputs/two_stage/mesh_comparison.csv` and `mesh_fit_summary.csv`).

**Configurations.** Node counts are ranges over the nine types. Edge is the median data-triangle edge.

| tag | ω mesh | ω nodes | ω edge (km) | ξ mesh | ξ nodes | ξ edge (km) |
|---|---|---|---|---|---|---|
| base | cutoff 30 km, inner edge ≤ 250 km, ≤ 2500 nodes | 1231–2294 | 82–104 | ≤ 600 nodes | 475–597 | 246–313 |
| xi1200 | base | 1231–2294 | 82–104 | ≤ 1200 nodes | 1114–1199 | 101–176 |
| xifull | base | 1231–2294 | 82–104 | the ω mesh | 1231–2294 | 82–104 |
| omega5000 | cutoff 15 km, inner edge ≤ 150 km, ≤ 5000 nodes | 2220–4357 | 49–71 | ≤ 1200 nodes | 1114–1199 | 101–176 |
| omega5000_xi2500 | as omega5000 | 2220–4357 | 49–71 | ≤ 2500 nodes (= the base ω mesh) | 1231–2294 | 82–104 |

**Scores.** Pooled over each experiment's folds. Differences are from `base`, with 95% paired pixel-bootstrap intervals. The dynamical model scores −4.329 (blocks) and −3.863 (interpolation).

| experiment | tag | log score | Δ log score | Δ CRPS (×10⁻³) | variance explained | Δ explained (points) | cover 50 / 95 |
|---|---|---|---|---|---|---|---|
| blocks (n = 8694) | base | −3.779 | | | 35.3 | | 0.388 / 0.859 |
| | xi1200 | −3.742 | +0.037 [+0.029, +0.045] | −2.3 [−2.9, −1.6] | 36.9 | +1.6 [+1.2, +2.1] | 0.392 / 0.872 |
| | xifull | −3.733 | +0.046 [+0.037, +0.055] | −2.6 [−3.3, −2.0] | 36.9 | +1.7 [+1.2, +2.2] | 0.397 / 0.875 |
| | omega5000 | −3.732 | +0.048 [+0.039, +0.058] | −2.9 [−3.5, −2.2] | 37.2 | +1.9 [+1.4, +2.4] | 0.396 / 0.875 |
| | omega5000_xi2500 | **−3.724** | +0.055 [+0.046, +0.065] | −3.2 [−3.9, −2.5] | 37.4 | +2.1 [+1.6, +2.6] | 0.401 / 0.878 |
| interpolation (n = 1045) | base | −3.591 | | | 51.3 | | 0.363 / 0.859 |
| | xi1200 | −3.577 | +0.014 [+0.005, +0.023] | −1.0 [−1.8, −0.2] | 51.8 | +0.5 [−0.1, +1.0] | 0.362 / 0.861 |
| | xifull | −3.573 | +0.017 [+0.009, +0.026] | −1.3 [−2.0, −0.6] | 52.2 | +0.9 [+0.4, +1.4] | 0.365 / 0.862 |
| | omega5000 | −3.575 | +0.015 [−0.005, +0.034] | −0.5 [−2.1, +1.0] | 51.4 | +0.1 [−1.1, +1.2] | 0.364 / 0.858 |
| | omega5000_xi2500 | **−3.568** | +0.022 [+0.005, +0.042] | −1.4 [−3.0, 0.0] | 52.1 | +0.8 [−0.1, +1.8] | 0.364 / 0.859 |

Head-to-head on blocks, where the test set is 8× larger:

- The finer ξ mesh gains on top of the finer ω mesh: omega5000_xi2500 − omega5000 = +0.008 [+0.004, +0.012].
- The finer ω mesh gains on top of the finer ξ mesh: omega5000_xi2500 − xifull = +0.0095 [+0.006, +0.013].
- On interpolation, the head-to-head differences are within ±0.01, and every interval spans 0.

The blocks gain is in the pyrethroids, DDT and Bendiocarb, which have large ξ (σ_η ≈ 0.4). The organophosphates and Alpha-cypermethrin are unchanged.

**Fitted ranges** (interpolation fold, km; ω's median data-triangle edge in brackets):

| type | range ω: base (edge) | range ω: omega5000_xi2500 (edge) | range η: base | range η: xi1200 | range η: xifull |
|---|---|---|---|---|---|
| Lambda-cyhalothrin | 116 (100) | 72 (63) | 1024 | 1134 | 871 |
| Permethrin | 88 (90) | 38 (55) | 1649 | 1612 | 1624 |
| Deltamethrin | 91 (83) | 54 (51) | 1754 | 1612 | 1561 |
| Fenitrothion | 127 (102) | 85 (65) | 97 | 126 | 133 |
| Bendiocarb | 127 (90) | 124 (56) | 1347 | 1381 | 1377 |
| DDT | 91 (85) | 59 (55) | 1287 | 1354 | 1208 |
| Alpha-cypermethrin | 808 (96) | 137 (69) | 11753 | 4019 | 4043 |
| Malathion | 5438 (96) | 638 (56) | 112 | 5089 | 5042 |
| Pirimiphos-methyl | 674 (95) | 167 (63) | 94 | 6009 | 852 |

- **ω's range shrinks as the mesh is refined, and stays at about one mesh edge.** For 5–7 of 9 types, it is below two edges under every configuration. The ω signal is partly finer than any affordable mesh resolves: between-site variation that the pixel-year u does not absorb. Refining ω still helps the scores.
- **η's range for the pyrethroids, DDT and Bendiocarb is stable at 1200–1750 km**, with φ ≈ 0.03–0.2. The ξ mesh does not limit these fields, yet the finer ξ mesh still improves their scores, through a better-resolved ξ at intermediate scales.
- **The ~95 km η range of the organophosphates on the base mesh was a mesh artefact.** At 250–300 km edges it jumps to 800–6000 km once the ξ mesh is refined. Fenitrothion keeps a short η range (126–135 km, σ_η ≈ 0.09).
- All fits converged under every configuration.

**Cost.** Per fold: the sum of the nine types' fit times, and the process peak RSS, which includes reading the dynamical fold. The wall times are confounded by varying machine load: the interpolation runs of xi1200, omega5000 and xifull shared the machine with the maps job and two MCMC fits.

| tag | fits per fold (min) | slowest type (min) | peak RSS (GB) | Deltamethrin alone: nnz(L), RSS |
|---|---|---|---|---|
| base | 24–28 | 4 | 5–7 | – |
| xi1200 | 38–68 | 5–11 | 6–7 | 50M, 4.6 GB |
| xifull | 31–69 | 6–13 | 10–11 | 94M, 9.2 GB |
| omega5000 | 18–62 | 3–9 | 6–7 | – |
| omega5000_xi2500 | 27–29 | 5–6 | 10–11 | – |

**Recommendation: omega5000_xi2500.**

- It is the best configuration on both experiments.
- On blocks, it beats every alternative with intervals clear of 0: +0.055 log score and +2.1 points of variance explained over base, and 95% coverage 0.859 → 0.878.
- It costs about 30 min per fold without contention, and about 11 GB peak per process. That is affordable for one job at a time.
- If memory is tight, omega5000 (ξ ≤ 1200 nodes, about 7 GB) keeps most of the gain: +0.048.

**Adopted.** omega5000_xi2500 is the reported configuration for both variants and the runner's default; `mesh=base` reproduces the earlier runs, whose files keep their unsuffixed names. `fit_correction()`'s defaults are unchanged. `R/two_stage_maps.R` uses it too (`mesh_config`).

## Results

Reported configuration: `m_ref` = posterior mean, meshes omega5000_xi2500, both variants. Scored by `R/two_stage_metrics.R`: `outputs/two_stage/cv_headline_two_stage.csv`, `term_inclusion_two_stage.csv`, `cv_horizon_two_stage.csv`; figures in `figures/two_stage/`. All models at the per-type rho. Intervals: 95% pixel bootstrap. Forecasting pooled stacks both origins' test sets (years from 2018 appear in both). NN oracle: nearest-neighbour null at the k that minimises its own held-out error.

**Headline.** Log score: mean beta-binomial log predictive density. Variance explained ceilings (noise floor): 78–82%.

| experiment | model | log score | CRPS | variance explained (%) | cover 50 / 95 |
|---|---|---|---|---|---|
| interpolation (n = 1045) | dynamical | −3.863 | 0.1355 | 37.0 [25.2, 45.5] | 0.298 / 0.786 |
| | intercept | −4.235 | 0.1496 | 27.0 [22.2, 30.2] | 0.306 / 0.721 |
| | nearest neighbour | −4.011 | 0.1356 | 33.7 [21.0, 43.8] | 0.342 / 0.776 |
| | NN oracle | −3.614 | 0.1108 | 52.0 [42.8, 59.4] | 0.391 / 0.833 |
| | two-stage ω+u | −3.614 | 0.1243 | 42.9 [31.6, 51.1] | 0.377 / 0.882 |
| | two-stage ω+ξ+u | **−3.568** | 0.1135 | 52.1 [42.9, 59.2] | 0.364 / 0.859 |
| blocks 1+2 (n = 8694) | dynamical | −4.329 | 0.1587 | 22.7 [18.0, 27.2] | 0.292 / 0.736 |
| | intercept | −4.587 | 0.1640 | 19.3 [17.0, 21.2] | 0.304 / 0.709 |
| | nearest neighbour | −4.679 | 0.1686 | 13.8 [7.1, 19.9] | 0.293 / 0.691 |
| | NN oracle | −4.190 | 0.1415 | 34.1 [29.1, 38.6] | 0.323 / 0.752 |
| | two-stage ω+u | −3.796 | 0.1423 | 31.3 [26.9, 35.6] | 0.371 / 0.859 |
| | two-stage ω+ξ+u | **−3.724** | 0.1335 | 37.4 [32.8, 41.5] | 0.401 / 0.878 |
| forecasting 2014+2018 (n = 14018) | dynamical | −4.361 | 0.1620 | 26.2 [20.9, 31.4] | 0.311 / 0.732 |
| | intercept | −4.706 | 0.1847 | 13.4 [10.5, 15.7] | 0.268 / 0.663 |
| | nearest neighbour | −4.323 | 0.1549 | 29.8 [26.0, 33.6] | 0.310 / 0.731 |
| | NN oracle | −3.869 | 0.1281 | 48.3 [45.8, 50.8] | 0.331 / 0.788 |
| | two-stage ω+u | −4.078 | 0.1610 | 25.1 [19.4, 30.7] | 0.349 / 0.782 |
| | two-stage ω+ξ+u | **−3.744** | 0.1438 | 35.1 [30.6, 39.7] | 0.408 / 0.856 |
| forecasting 2014 (n = 9922) | dynamical | −4.323 | 0.1588 | 21.6 [15.2, 27.3] | 0.319 / 0.737 |
| | intercept | −4.623 | 0.1760 | 12.7 [9.6, 15.1] | 0.269 / 0.673 |
| | nearest neighbour | −4.314 | 0.1529 | 25.0 [20.4, 29.3] | 0.315 / 0.728 |
| | NN oracle | −3.874 | 0.1280 | 43.8 [40.8, 46.5] | 0.326 / 0.779 |
| | two-stage ω+u | −4.020 | 0.1580 | 20.2 [13.1, 26.4] | 0.356 / 0.785 |
| | two-stage ω+ξ+u | **−3.703** | 0.1432 | 29.5 [23.8, 34.5] | 0.407 / 0.855 |
| forecasting 2018 (n = 4096) | dynamical | −4.453 | 0.1698 | 31.3 [23.0, 38.4] | 0.292 / 0.722 |
| | intercept | −4.905 | 0.2058 | 10.4 [5.7, 14.5] | 0.265 / 0.639 |
| | nearest neighbour | −4.344 | 0.1598 | 35.5 [28.0, 42.1] | 0.300 / 0.738 |
| | NN oracle | −3.858 | 0.1284 | 54.5 [49.8, 58.6] | 0.345 / 0.810 |
| | two-stage ω+u | −4.217 | 0.1682 | 31.0 [22.3, 38.8] | 0.334 / 0.773 |
| | two-stage ω+ξ+u | **−3.843** | 0.1453 | 42.8 [35.1, 49.6] | 0.411 / 0.858 |

Two-stage ω+ξ+u minus dynamical: log score +0.295 [+0.188, +0.397] (interpolation), +0.605 [+0.541, +0.670] (blocks), +0.617 [+0.554, +0.678] (forecasting); variance explained +15.1 [+9.6, +21.5], +14.7 [+12.4, +17.2], +8.9 [+6.9, +10.9] points.

- ω+ξ+u has the best log score everywhere, including over the NN oracle.
- On variance explained it ties the NN oracle on interpolation and trails it on forecasting (−13 points) and blocks (+3, intervals overlap).
- 95% intervals of both two-stage variants undercover (0.86–0.88 on the spatial experiments, 0.86 forecasting), less than every other model (0.64–0.83).

**Term inclusion.** ω+ξ+u minus ω+u, same meshes:

| experiment | Δ log score | Δ CRPS (×10⁻³) | Δ variance explained (points) |
|---|---|---|---|
| interpolation | +0.046 [−0.023, +0.117] | −10.8 [−17.6, −3.9] | +9.2 [+4.3, +14.5] |
| blocks 1+2 | +0.072 [+0.049, +0.094] | −8.8 [−11.2, −6.5] | +6.0 [+4.4, +7.7] |
| forecasting 2014+2018 | +0.334 [+0.300, +0.371] | −17.1 [−19.7, −14.8] | +10.0 [+8.3, +11.8] |
| forecasting 2014 | +0.317 [+0.278, +0.355] | −14.8 [−17.5, −12.3] | +9.3 [+7.4, +11.6] |
| forecasting 2018 | +0.374 [+0.321, +0.428] | −22.9 [−26.9, −18.9] | +11.8 [+9.4, +14.4] |

**Keep ξ.** It improves forecasting on every metric, and interpolation on CRPS and variance explained (log score interval spans 0). Without ξ, the correction does not improve forecast variance explained over the dynamical model (−1.0 [−2.8, +0.7]).

**Skill by horizon** (`cv_horizon_two_stage.csv`, `figures/two_stage/skill_by_horizon.png`, `skill_gain_by_horizon.png`). Every test year of both origins is after every type's last training year (T = 2013 and 2017), so all forecasting predictions use the AR(1) forecast of ξ; horizons are 1–5 years at both origins. Both origins pooled:

| horizon (years) | n | ω+ξ+u variance explained | Δ log score vs dynamical | Δ explained vs dynamical | ω+ξ+u cover 95 | NN oracle explained |
|---|---|---|---|---|---|---|
| 1 | 4228 | 47.3 [42.3, 51.8] | +0.53 [+0.46, +0.61] | +14.7 [+11.4, +18.2] | 0.886 | 49.0 |
| 2 | 3176 | 45.3 [39.5, 50.9] | +0.50 [+0.41, +0.58] | +5.5 [+2.5, +8.5] | 0.873 | 55.5 |
| 3 | 2247 | 33.2 [21.0, 43.3] | +0.45 [+0.36, +0.55] | +4.2 [+0.8, +7.6] | 0.852 | 47.0 |
| 4 | 2222 | 17.7 [4.9, 28.1] | +0.85 [+0.70, +1.01] | +6.6 [+3.1, +9.9] | 0.825 | 45.0 |
| 5 | 2145 | 18.8 [7.3, 29.2] | +0.89 [+0.76, +1.02] | +11.4 [+8.4, +14.4] | 0.806 | 40.8 |

- Skill of all model-based forecasts falls with horizon; the two-stage gain over the dynamical model does not shrink (log score gain is largest at 4–5 years, where the dynamical model's own skill is near 0).
- 95% coverage of ω+ξ+u falls from 0.89 to 0.81 with horizon: the AR(1) forecast of ξ is overconfident at long horizons.
- The NN oracle's variance explained stays at 41–56% at every horizon.

**Stage B diagnostic** (`train_residual_extremes.csv`, reported meshes, both variants, all 5 folds). Standardised training residuals, observed vs simulated from the stage-A fit:

| mortality class | assays observed / implied by fit | mean residual observed / simulated | SD observed / simulated |
|---|---|---|---|
| 0% | 1.1–1.7× | −0.79 to −1.00 / −0.55 to −0.72 | 0.28–0.34 / 0.17–0.29 |
| interior | 0.87–0.88× | +0.01 to +0.05 / −0.10 to −0.12 | 0.75–0.80 / 0.92–0.93 |
| 100% | 1.47–1.61× | +0.34 to +0.40 / +0.24 to +0.29 | 0.26–0.32 / 0.23–0.28 |

- The same pattern on every fold and both variants: 47–61% more 100%-mortality assays than the fit implies, and residuals at both extremes further out than the selection effect explains.
- Interior residuals are narrower than v implies (SD 0.78 vs 0.93).
- So the Gaussian empirical-logit response misfits the extremes: the fitted latent is not extreme enough where mortality saturates. Stage B (PQL on the beta-binomial) is indicated; its effect on the scores is untested.
