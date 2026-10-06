# Two-stage correction: run sheet

Issue idem-lab/ir_cube#21. Stacked on #12 (`posterior-predictive-validation`). Code removed during development is kept at two tags:

- `two-stage-full-e31b2b5`: the diagnostics, the rejected variants (`omega_u`, `m_ref=loo`, the other meshes) and the separate simulation checks.
- `two-stage-options-6edcf0d`: the survey effect s, joint ρ and damped ξ.

## Final model

Per insecticide type, on the logit scale, for assay i at pixel s and year t:

λ_i = m_i + ω(s) + ξ(s, t) + u(s, t) + p(s)

- m: the dynamical model's posterior mean logit at the training assays (`m_ref`).
- ω: a static Matérn field, correcting the initial conditions.
- ξ(s, t) = Σ η: the accumulation from ξ(·, t0 = 1995) = 0 of annual anomalies η, AR(1) in time (φ) with Matérn innovations. It is undamped, so a forecast holds the accumulated deviation and its mean plateaus at a rate set by φ.
- u (per pixel-year, SD τ) and p (per pixel, SD σ_p): iid.
- Meshes (`build_correction_meshes()`):
  - In the data region, ω has ≤ 5000 nodes (15 km cutoff, inner edges ≤ 150 km) and ξ has ≤ 2500 nodes.
  - Around both, a coarse outer region (edges ≤ 800 km) extends 1500 km beyond the data and the prediction mask, so every map cell is inside the mesh. Callers pass the mask's coordinates (`prediction_mask_coords()`, from `data/clean/raster_mask.tif`).
- PC priors: P(range < 50 km) = 0.05 and P(σ > 1) = 0.05 for both fields, and P(sd > 1) = 0.05 for u and p. Persistence 1/(1 − φ) is lognormal with median 5 years.
- Fitting (`fit_correction()`):
  - Stage A: the empirical logit, with its variance inflated by the external per-type ρ (`outputs/bioassay_rho_hierarchical.csv`). Hyperparameters by penalised maximum marginal likelihood in TMB.
  - Stage B: PQL on the beta-binomial counts from the stage-A fit (`R/two_stage_pql.R`), with one re-estimation of the hyperparameters.
  - The latent posterior is N(mode, H⁻¹), with the hyperparameters plugged in.
- Cut posterior: the dynamical model is never updated. Each dynamical draw shifts the latent mode by H⁻¹A′D(m_ref − m_draw), and is paired with its own latent draw.

**Target.** The prediction target is m + ω + ξ: that is what the maps, the published figures and the supplement figures show. u and p are observation-level noise. They are in the fitted model, so that this noise stays out of ω and ξ. They enter predictions only of new assays, and only as fresh draws, N(0, τ²) per pixel-year and N(0, σ_p²) per pixel, never as their fitted values. This reverses #21, which put u in the target, and an earlier version of this run sheet, which put p in it. Cross-validation scores the predictive distribution of held-out assays, so it includes fresh u and p (`predict_correction()`).

## Headline results

Scored as in #12 (`R/validation_metrics.R` → `outputs/cv_scores.csv`, `cv_summary.csv`), where the final model is `two_stage`. Variance explained, its ceiling and its paired contrasts are from `R/variance_explained.R` → `outputs/cv_variance_explained.csv`, `_by_fold.csv` and `_by_horizon.csv`. The other paired differences are from `R/two_stage_metrics.R` → `outputs/two_stage/cv_headline_two_stage.csv` and `cv_horizon_two_stage.csv`.

- Log score: mean beta-binomial log predictive density at the per-type ρ.
- CRPS: on the mortality scale.
- Variance explained: 100 (1 − MSE / Var(y)). The ceiling is the `noise_floor_mse()` bound.
- Coverage: the expected coverage of the central interval over the PIT randomisation (`expected_coverage()`), in every table. #12's tables used the mean over 100 PIT replicates (`cv_summary.csv`) or one replicate (`cv_by_fold.csv`); the switch changes #12's models' coverage by at most 0.002 in `cv_summary.csv` and 0.009 in `cv_by_fold.csv`.
- Intervals: 95% pixel bootstrap.
- NN: the nearest recent survey (k = 1). NN oracle: the nearest surveys at the k that minimises held-out error.

**Interpolation** (n = 1045; ceiling 78.0%)

| model | log score | CRPS | variance explained (%) | cover 50 / 95 |
|---|---|---|---|---|
| dynamical | −3.863 | 0.1353 | 37.0 [26.0, 46.3] | 0.298 / 0.787 |
| intercept | −4.233 | 0.1497 | 27.0 [22.7, 30.3] | 0.306 / 0.721 |
| NN | −4.010 | 0.1355 | 33.7 [20.5, 44.1] | 0.343 / 0.777 |
| NN oracle | −3.614 | 0.1107 | 52.0 [43.0, 59.7] | 0.390 / 0.833 |
| **two-stage** | **−3.372** | **0.1028** | **57.1 [49.2, 63.8]** | 0.448 / 0.921 |

**Spatial blocks 1+2** (n = 8694; ceiling 79.5%)

| model | log score | CRPS | variance explained (%) | cover 50 / 95 |
|---|---|---|---|---|
| dynamical | −4.327 | 0.1587 | 22.7 [17.9, 26.9] | 0.292 / 0.736 |
| intercept | −4.585 | 0.1640 | 19.3 [17.1, 21.3] | 0.304 / 0.709 |
| NN | −4.678 | 0.1685 | 13.8 [7.3, 19.8] | 0.293 / 0.691 |
| NN oracle | −4.189 | 0.1415 | 34.1 [29.3, 38.5] | 0.323 / 0.752 |
| **two-stage** | **−3.663** | **0.1327** | **36.9 [32.1, 41.3]** | 0.440 / 0.911 |

**Forecasting 2014+2018** (n = 14018; ceiling 82.8%)

| model | log score | CRPS | variance explained (%) | cover 50 / 95 |
|---|---|---|---|---|
| dynamical | −4.359 | 0.1619 | 26.2 [20.6, 31.4] | 0.311 / 0.732 |
| intercept | −4.703 | 0.1847 | 13.4 [10.8, 15.9] | 0.268 / 0.663 |
| NN | −4.322 | 0.1549 | 29.8 [25.7, 33.4] | 0.311 / 0.731 |
| NN oracle | −3.869 | **0.1281** | **48.3 [45.6, 50.9]** | 0.332 / 0.788 |
| **two-stage** | **−3.614** | 0.1315 | 42.7 [38.7, 46.6] | 0.457 / 0.900 |

**Two-stage minus dynamical**, paired:

| experiment | Δ log score | Δ CRPS (×10⁻³) | Δ variance explained (points) |
|---|---|---|---|
| interpolation | +0.491 [+0.366, +0.618] | −32.4 [−41.1, −23.4] | +20.1 [+13.5, +28.6] |
| blocks 1+2 | +0.664 [+0.597, +0.734] | −26.0 [−29.7, −22.3] | +14.2 [+11.6, +16.8] |
| forecasting 2014+2018 | +0.746 [+0.672, +0.825] | −30.4 [−34.6, −26.6] | +16.6 [+14.0, +19.2] |
| forecasting 2014 | +0.772 [+0.690, +0.860] | −29.7 [−34.3, −25.6] | +17.5 [+14.4, +20.8] |
| forecasting 2018 | +0.682 [+0.585, +0.784] | −32.3 [−38.4, −26.8] | +15.6 [+11.9, +19.0] |

**By forecast horizon** (both origins pooled):

| horizon (years) | n | variance explained | Δ log score vs dynamical | Δ explained vs dynamical | cover 95 | NN oracle explained |
|---|---|---|---|---|---|---|
| 1 | 4228 | 51.1 [46.0, 55.5] | +0.62 [+0.53, +0.71] | +18.5 [+14.8, +22.4] | 0.914 | 49.0 |
| 2 | 3176 | 51.9 [46.0, 57.1] | +0.61 [+0.51, +0.73] | +12.0 [+8.2, +16.0] | 0.914 | 55.5 |
| 3 | 2247 | 39.2 [27.5, 49.3] | +0.55 [+0.44, +0.66] | +10.3 [+6.2, +14.3] | 0.883 | 47.0 |
| 4 | 2222 | 30.4 [21.1, 38.6] | +1.05 [+0.86, +1.23] | +19.4 [+14.4, +24.8] | 0.882 | 45.0 |
| 5 | 2145 | 30.9 [21.1, 39.2] | +1.09 [+0.93, +1.26] | +23.5 [+19.4, +28.2] | 0.887 | 40.8 |

- The two-stage model has the best log score in every experiment.
- On variance explained it beats the NN oracle on interpolation (+5) and blocks (+3), and trails it on forecasting (−6), mostly at horizons of 3–5 years.
- Its 95% coverage is 0.89–0.92. Every other model's is 0.66–0.83.
- These are from the rerun of the five folds with the current code (meshes covering the prediction mask, `m_ref` unclamped, seeds per fold and type). Against the earlier draws, the posterior mean predictions correlate at 0.9998–0.9999 per fold (mean |difference| 0.2–0.4 points of mortality), and the headline scores move by at most 0.002 in log score, 0.2 points of variance explained and 0.002 in coverage.

## Decisions

Scores are from the cross-validation of the time. A decision's numbers can predate later choices; the tables above are the current model.

| option | result | adopted | code |
|---|---|---|---|
| ξ (`omega_xi_u` vs `omega_u`) | Forecasting +0.33 log score, +10 points explained; blocks +0.07, +6. Without ξ, forecasting gains nothing over the dynamical model in variance explained. | yes | full tag |
| Meshes omega5000_xi2500 (vs four coarser configurations) | Best everywhere. Blocks +0.055 log score and +2.1 points over the base mesh; ~11 GB per fit. | yes | full tag |
| `m_ref` = posterior mean (vs PSIS leave-out mean) | Leave-out shift 0.04–0.05 residual SD; double counting negligible. | posterior mean | full tag |
| Stage B, PQL on the beta-binomial (vs Gaussian empirical logit) | +0.06 to +0.19 log score; 95% coverage 0.86–0.88 → 0.90–0.92; removes the excess of 100% assays in the training residuals. | yes | — |
| Pixel effect p | +0.003 to +0.005 log score; about 30% fewer ω hotspots. | yes, as observation noise | — |
| Survey effect s | +0.006 to +0.013 log score, all from a wider held-out predictive distribution; none in the mean or the map. | no | options tag |
| Joint ρ | Estimated ρ is a median of 0.35 × the replicate estimate, with the noise moved into u (a PQL bias also seen in simulation); −0.04 to −0.07 log score at its own ρ. | no | options tag |
| Damped ξ (ψ) | Forecasting −0.038 log score, −0.12 and −0.16 at 4–5 years; spatial experiments unchanged. | no | options tag |
| u and p in the target | u and p are observation noise, as fresh draws for new assays only (see Final model). | target m + ω + ξ | — |
| Diagnostics: covariate and hotspot structure, training residuals at 0% and 100%, mesh comparison, stage-one effective parameters | Finished. They motivated stage B, p and the meshes. | removed | full tag |

## Open issues

- #28: the temporal panel of #12's variance-explained bar charts should score change, not level.
- #34: check the stage-B inference before the publishable run: convergence of the hyperparameters, PQL against a Laplace fit of the beta-binomial, a Gaussian latent, and plugged-in hyperparameters.
- Wide intervals away from data. ξ is close to a random walk (φ ≈ 0.01–0.07 at stage B), so its prior SD grows as σ_η√t. Far from data, and at long horizons, the two-stage SD map is wide and the posterior mean mortality is pulled towards 50%. Damping ξ did not help held-out forecasts (Decisions).
- Residual narrowness. After stage B, interior training residuals are narrower than v implies (SD 0.70 vs 0.95). A lower ρ does not explain it (joint ρ). The stage-A standardisation used by the diagnostic is a candidate.

## How to run

Use OpenBLAS (`LD_PRELOAD=.../libopenblas.so.0 OPENBLAS_NUM_THREADS=3`), because reference BLAS is ~10× slower. Compile the template with TMBad, which `load_correction_template()` does.

1. **Simulation check:** `Rscript R/check_two_stage.R`.
2. **Folds** (`R/run_two_stage_folds.R`), for each of `spatial_interpolation all`, `spatial_blocks 1|2` and `temporal_forecasting 2014|2018`:
   - `Rscript R/run_two_stage_folds.R <experiment> <fold> stage_one`: rebuilds the fold and recomputes the paired dynamical draws. Loading a fold takes 4–8 GB.
   - `... <experiment> <fold> <k>` for k = 1–9: one process per type. Peak memory is up to ~15 GB (Deltamethrin).
   - `... <experiment> <fold> assemble`: writes `outputs/cv_draws/two_stage__<experiment>__<fold>.rds` and `outputs/two_stage/fit_summary.csv`.
3. **Metrics:** `Rscript R/validation_metrics.R`, then `Rscript R/variance_explained.R`, then `Rscript R/two_stage_metrics.R`.
4. **Maps** (`R/two_stage_maps.R`):
   - `prepare`: loads `temporary/fitted_model.RData`, checks that its data and covariates are those the current scripts build, and saves per-type dynamical draws (2 min, 7 GB).
   - `fit <k>`: fits type k to all data and saves `outputs/two_stage/maps/<type>/fit.rds` (the fit without its TMB object) and `hyperparameters.csv`. 7–20 min and 8–14 GB per type (Deltamethrin 13.9).
   - `map <k>` for the six types not in LLINs, and `map llin_effective`, which maps Alpha-cypermethrin, Deltamethrin and Permethrin together with their combination, weighted draw by draw by `temporary/ingredient_weights.RDS`. Each writes the posterior mean and SD of mortality for every year 1995–2030 to `outputs/two_stage/ir_maps/<output>/ir_<year>_susceptibility(_sd).tif`, the layout of the dynamical model's `outputs/ir_maps` (`R/predict.R`), from 1000 paired draws in batches of 100 (Monte Carlo SE of a mean map: posterior SD / √1000, at most 1.3 percentage points where the SD reaches its maximum of 40). Set `IR_CUBE_MAP_WORKERS` to map groups of cells in forked workers. On the Linux box a type took 23–28 min at 3–4 workers and `llin_effective` 81 min at 3 workers; each worker holds ~5 GB, so run one map at a time.
   - `figures`: writes `figures/two_stage/`. Per type: the posterior SD of mortality (percentage points, including the dynamical model's uncertainty), the correction ω + ξ, and the difference from the dynamical model. Also the all-types correction panel for 2020 and the hyperparameters.
5. **Supplement figures:** `Rscript R/fig_two_stage_components.R`. Reads the saved Deltamethrin fit and its draws, and writes `supp_components_realisations` and `supp_components_intervals` (captions: `doc/two_stage_supplement_captions.md`).
6. **Published figures.** The figure scripts that show predictions use the two-stage model: `R/fig_ir_maps.R`, `fig_ir_map_2025.R`, `fig_ir_map_2000_2030.R` and `fig_ento_epi_impact.R` read the rasters of step 4 (`ir_map_files()` in `R/functions.R`), and `R/fig_temporal_preds_data.R`, `fig_temporal_preds_net_use.R`, `summarise_model_fit.R` and `fig_internal_validation.R` draw from the saved fits at the data cells (`R/two_stage_predictions.R`). Each of the latter takes 4–11 min; `fig_ento_epi_impact.R` takes ~1 h.
