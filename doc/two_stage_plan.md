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

## Term inclusion

Keep `xi` only if `omega_xi_u` improves the interpolation or forecasting scores over `omega_u`.
