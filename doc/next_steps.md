# Next steps, as of 25 September 2026

> **Prompt to resume from this file.** Read `doc/next_steps.md`. Work on branch
> `posterior-predictive-validation` (PR #12). §1 is running and needs no action
> until it lands. Start at whichever of §2–§6 the user asks for. Ask before
> posting anything to GitHub. When editing the bar-chart figures, render and
> **look at the PNG** with the Read tool after every change — reasoning about
> ggplot spacing from the code does not work and cost several rounds.

## 1. Running now, no action needed

Two forecast folds, launched 10:42 on 23 September from the frozen copy at
`/tmp/claude-1000/.../scratchpad/frozen`, two concurrent at four threads:

| fold | progress at 09:00, 25 Sep | expected |
|---|---|---|
| `temporal_forecasting 2018` | 3,500 / 5,000 sampling | late 25 Sep |
| `temporal_forecasting 2014` | 2,850 / 5,000 sampling | 26 Sep |

Logs in `outputs/cv_logs/temporal_forecasting__{2014,2018}.log`. Each writes
`outputs/cv_draws/dynamical__temporal_forecasting__{2014,2018}.rds`.

**Do not edit anything under `R/` that the run reads.** R reads `--file=`
incrementally, so editing a script mid-run corrupts that run; the frozen copy
exists for this reason. Two folds were lost this way at 62 h each.

## 2. Finish the forecast rebuild — `doc/forecast_rebuild_plan.md` §4.4–4.6

Once both folds are on disk:

1. `Rscript R/validation_metrics.R` — about 40 minutes, picks up whatever is in
   `outputs/cv_draws/`.
2. `Rscript R/validation_change.R` — already generalised over origins, emits one
   row per cut.
3. `Rscript R/validation_geometry.R` — distance-stratified tables.
4. `Rscript R/variance_explained.R` then `Rscript R/fig_variance_explained.R` —
   the temporal change bars currently use the superseded 2020 three-year fold,
   which is the only forecasting fit that existed. Decide whether the figure
   shows one origin or both; if both, `experiments` in `variance_explained.R`
   needs a second forecasting entry and the facet count changes.
5. `R/fig_predictive_validation.R` still needs its forecasting panels
   generalised from one origin to several. The label function handles the
   naming; the panels do not.

Check Rhat on each fold as it lands. The block folds came in at 1.104 and
1.239, the production fit at 1.215, so anything under about 1.25 is in family.

## 3. The bar-chart figures need their spacing redone

`R/fig_variance_explained.R`, producing `figures/CV_variance_explained*.png/pdf`
and `outputs/figure_variance_explained.RDS`.

The encoding is settled and correct: each bar is the whole variance in held-out
mortality, light grey; the bioassay-noise share is washed toward white from the
top with a dashed rule at the ceiling; the model's share is coloured from the
bottom, solid to the lower 95% bound and translucent to the upper with a rule
at the estimate; three rotated key labels name the regions in the right margin.
Caption material for all three figures is in comments above each `ggsave`.

**What is wrong: the spacing, in several places, and it got worse over the last
few rounds rather than better.** Redo it by rendering and looking, not by
reasoning about the code. Specific things that have gone wrong before and are
worth checking explicitly:

- Bars flush with the panel edges, because the facet strip label aligns to the
  panel edge and needs to line up with the first bar. This is why the x window
  is set with `coord_cartesian(xlim = ...)` and not with scale expansion.
- The key geoms must not train the x scale. They sit past the last bar and are
  drawn with `clip = "off"` into a right margin set in `plot.margin`. Putting a
  hard `limits` on the x scale silently drops them.
- Gaps: between the paired bars within an insecticide, and between bar groups.
  These are currently too tight in places.
- The two-line "bioassay noise" label straddles its anchor, so its offset is
  larger than the single-line labels', and it is nudged toward the top of its
  band so it does not read as one phrase with "unexplained variance".
- The per-insecticide key references the best dynamical bar in each panel, not
  the first: the first alphabetically is Alpha-cypermethrin at about -137% on
  interpolation, which pushes two of three key segments off the axis.
- Negative bars stay clipped and absent, by decision. Do not make the axis
  negative.

## 4. Swap the per-type overdispersion into the validation floor

**This is now unblocked.** The hierarchical MCMC fit converged this morning:
worst Rhat 1.040, minimum effective sample size 129, in
`outputs/bioassay_rho_hierarchical.csv`.

| insecticide | rho | 95% |
|---|---|---|
| Alpha-cypermethrin | 0.252 | 0.218 – 0.293 |
| Malathion | 0.244 | 0.196 – 0.302 |
| Pirimiphos-methyl | 0.179 | 0.116 – 0.261 |
| Permethrin | 0.176 | 0.162 – 0.191 |
| Deltamethrin | 0.154 | 0.145 – 0.163 |
| Bendiocarb | 0.125 | 0.108 – 0.145 |
| Fenitrothion | 0.121 | 0.088 – 0.163 |
| DDT | 0.118 | 0.105 – 0.131 |
| Lambda-cyhalothrin | 0.095 | 0.079 – 0.113 |

`R/variance_explained.R` already prefers this file when `worst_rhat < 1.05`, so
rerunning it will pick the per-type values up automatically — it fell back to
per-class maximum likelihood only because the fit had not converged when it last
ran. What still uses the per-class values, and should be moved over:

- `R/validation_metrics.R` — `rho_for_class()`, which sets the noise floor for
  every scored fold.
- `R/validation_change.R` — same lookup.
- `R/validation_distance.R` — same.

This matters most where the class value is furthest from the type's own:
Lambda-cyhalothrin's floor is overstated 72% by the pyrethroid value and
Fenitrothion's 78% by the organophosphate one. Fenitrothion on the spatial
blocks currently scores an excess mean squared error of **-0.0003**, i.e. below
its own floor, which is the arithmetic signature of an inflated floor.

No refitting is needed for any of this. Minutes, not hours.

## 5. Not started, and deliberately out of scope for this PR

- **Issue #20**, type-level rho inside the dynamical model. Needs a refit, so it
  belongs after the forecast folds. The implementation to lift is in
  `R/fig_illustrate_bioassay_variability.R`.
- The cheap test proposed for #20 before committing 62 h: a betabinomial GLM on
  the same covariates with class- versus type-level dispersion, to see whether
  the dispersion structure moves predictions at all.

## 6. Housekeeping

Uncommitted, and all of it wanted:

- `R/variance_explained.R`, `R/fig_variance_explained.R` — new.
- `R/fig_illustrate_bioassay_variability.R` — now fits the hierarchical rho by
  MCMC and writes `outputs/bioassay_rho_hierarchical.csv`; its own figure uses
  each example's own type-level rho.
- `R/null_models.R` — the guard against a missing neighbour count.
- `figures/bioassay_variability.png`, `figures/cluster_sampling_power.png`,
  `figures/CV_variance_explained*`, `figures/CV_distance_smooth.png`.

Untracked and disposable: `Rplots.pdf`, `distinct_pts.html`. Superseded draws:
`outputs/cv_draws_defunct/`, `outputs/cv_draws_leaky_forecast/`,
`outputs/cv_draws_superseded_nn/` (this last has a README explaining why).

Analysis scripts still living in the scratchpad rather than `R/`, worth
promoting if any of their outputs are going in the paper:
`predictable_variance.R`, `decomposition.R`, `ranking.R`,
`decomposition_by_insecticide.R`, `floor_free.R`, `shrinkage.R`,
`rho_weighting.R`, `rho_hier.R`, `rho_groups.R`.

## 7. Results settled in this conversation, for the write-up

**Variance explained out of sample, floor-free** (`outputs/cv_variance_explained.csv`),
as % of observed variance in held-out mortality, 95% pixel-cluster bootstrap:

| experiment | dynamical | nearest recent survey | best-k NN | insecticide mean | bioassay noise |
|---|---|---|---|---|---|
| spatial interpolation | 37.0 [26.2, 46.0] | 33.7 [20.7, 44.0] | 52.0 | 27.0 | 21.8 |
| spatial extrapolation | 22.7 [17.9, 26.9] | 13.8 [7.3, 19.8] | 34.1 | 19.3 | 20.6 |
| temporal change (2020) | 35.7 [21.3, 47.7] | 28.2 [14.9, 39.7] | 56.7 | 18.4 | 15.5 |

**Against the intercept null** (`outputs/cv_skill_ci.csv`) the dynamical model is
not distinguishable from it in any headline row: blocks pooled +5.5 [-2.9,
+12.8], interpolation +19.4 [-1.4, +35.6], block fold 1 -10.1 [-23.7, +1.4].
Only block fold 2 clears zero, at +24.1. It does beat the one-neighbour
baseline on blocks pooled, +14.7 [+5.6, +24.5], and loses to the best-k bound
everywhere, -19 to -29.

**Murphy decomposition** (`outputs/cv_decomposition.csv`): on the blocks the
dynamical model's discrimination equals the single nearest neighbour's — 36.9
against 37.0 — so its whole advantage over that baseline is calibration, 8.3
against 19.6. Within insecticide it ranks two same-insecticide locations
correctly 64% of the time on blocks and 69% on interpolation, against 50% for no
information, and slightly *worse* than the nearest neighbour.

**Predictions are over-dispersed, not over-shrunk**: the optimal recalibration
slope is 0.82 on blocks and 0.79 on interpolation. This refuted the proposed
mechanism by which an inflated rho would compress predictions, and it is
consistent with issue #14.

**No distance effect is resolvable** (`outputs/cv_distance.csv`, §6b of
`doc/cv_run_plan.md`): joint test p = 0.058 for the dynamical model, 0.264 for
one-neighbour, and a one-degree-of-freedom slope of +7.2 [-1.25, +15.7] % per
e-fold of distance. Closed as a negative result.

**The nearest neighbour null was never run on the block folds** before this
conversation — a missing row in `optimal_nn.csv` reduced it to Beta(0.5, 0.5) —
so the previously reported -0.86 variance explained was the score of a prior.
The null is no longer tuned at all; it is reported as a practice baseline at one
neighbour and as an oracle bound at its best k.

**The two rho estimators agree**, which validates the Laplace intervals: on the
same six most-sampled groups, maximum likelihood gives 0.1948 [0.1467, 0.2539]
and MCMC 0.1948 [0.1473, 0.2499]. On all 3,713 replicated groups the two models
separate as expected, the hierarchical marginal likelihood giving 0.1547 [0.1489,
0.1606] against 0.1769 [0.1705, 0.1839] for flat independent priors on the group
fractions — the hyperprior does real work when most groups hold two assays.
