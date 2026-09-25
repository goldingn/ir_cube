# Two-stage correction: run sheet

Issue idem-lab/ir_cube#21. Stacked on the #12 branch (`posterior-predictive-validation`).

## Folds

Five outer folds, the same as #12, with no change to the fold design:

| experiment | fold | dynamical draws |
|---|---|---|
| spatial_interpolation | all | available |
| spatial_blocks | 1 | available |
| spatial_blocks | 2 | available |
| temporal_forecasting | 2014 | still running on #12 |
| temporal_forecasting | 2018 | still running on #12 |

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

- Per type, p_WAIC and p_loo are 23–53.
- pD with the plug-in at the posterior mean of the *parameters* is unreliable: it is negative for Alpha-cypermethrin on interpolation. The averaged non-centred parameters, pushed through exp(beta) and the recursion, are far from the posterior mode. So use pD at the posterior mean of p, or p_loo.

Grouped PSIS, leaving out a pixel-year (and, in brackets, a whole pixel):

| fold | groups with k > 0.7 | residual SD | posterior SD / residual SD | RMS leave-out shift / residual SD | leave-out / in-sample residual SD |
|---|---|---|---|---|---|
| interpolation | 0.10% (0.32%) | 1.94 | 0.12 | 0.036 (0.052) | 1.010 |
| blocks 1 | 0.13% (0.49%) | 1.90 | 0.13 | 0.041 (0.059) | 1.012 |
| blocks 2 | 0.15% (0.55%) | 1.89 | 0.13 | 0.039 (0.059) | 1.011 |

- **Double counting is negligible.** Fitting stage A to in-sample residuals understates the residual SD by 1–4%.
- **The organophosphates are affected most, and still only slightly:** RMS shift 0.07–0.11 of the residual SD, and posterior SD 0.26–0.41 of the residual SD.
- **Decision:** keep `m_ref` as the posterior mean. `m_ref=loo` stays available in the runner as a sensitivity check.
- **Caveat on the saved draws:** 282 of the interpolation fold's 2000 paired draws are exact repeats (81 in blocks 1, 311 in blocks 2), where HMC stuck. PSIS is run on the distinct draws only; the repeats give spurious k = Inf.

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

The runner's default is still `base`, so earlier results reproduce. `fit_correction()`'s defaults are unchanged. The maps build their meshes explicitly, so the defaults would not reach them anyway.

To adopt omega5000_xi2500:

- Run the remaining folds with `mesh=omega5000_xi2500`.
- Rerun `omega_u` with the same ω mesh, so the term comparison stays paired.
- Change the two `build_correction_mesh()` calls in `R/two_stage_maps.R` to `cutoff = 15, max_edge_inner = 150, max_nodes = 5000` for ω, and `max_nodes = 2500` for ξ. Fitted to all the data, ω may reach the 5000-node cap, and memory may be above 11 GB.
