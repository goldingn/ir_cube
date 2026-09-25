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
