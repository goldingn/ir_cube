# Misfit screens and the latent kdr surface: status and next steps

Written 6 October 2026 (Claude Code, for Nick), before compacting a long
session. Two threads: PR idem-lab/ir_cube#44 (the data update, nearly done)
and the #37 misfit work that has led to a plan for a latent genetic surface.

## 1. PR #44 (milestone 3, data update) — nearly done

- Branch `data-update` on goldingn (worktree `../ir_cube_data_update`), PR
  #44 against idem-lab master. Pushed head `c4a6889`; local commit `db72c06`
  (drops the cumulative net-use transform and the grid's saturating arm) is
  **not yet pushed**.
- Contents: #31 cross-database dedup; #27 October 2026 MTM, VA S/M form
  labels, the 306 zero-filled MTM records dropped; #33 MITN run 06 net use
  2000-2025 (2 Oct download); #37 pre-screen scripts and grid launcher
  (`R/net_grid.R`, now d_half x floor only, linear nets); `pop_d_half_spread.R`
  on the cleaned data (max spread at d_half = 270).
- To do before merge:
  1. Push `db72c06`.
  2. Rewrite the PR description (it was posted without Nick's review and
     still recommends the cumulative transform). Draft in chat first.
  3. Decide which #37 screen scripts from `net-transform-screen` go into #44
     (probably none; they belong to the next PR). Do NOT merge
     `net-transform-screen` again as is: its `2a07396` re-adds the cumulative
     transform that #44 removed.
  4. Draft the #37 comment (Nick approves before posting). #27/#31/#33 need
     no comments; #44 closes them.
- Follow-ups outside #44: `final_data_year` -> 2025 after the figures PR
  (adds 32 bioassays); newer IR Mapper data for the final run (#35); remove
  worktrees after merge.
- Volume: `data/clean` on the RunPod volume now holds the #44 data (the old
  copy is identical to `../ir_cube_refit/data/clean`).

## 2. What the #37 screens established

All on the RunPod fit `dh270_lin_f0_full` (d_half = 270, floor fixed at 0,
linear nets, w = 0.25, #44 data; chain 4 stuck and dropped; local copy
`outputs/pod_jobs/dh270_lin_f0_full/temporary/fitted_model.RData`; stage 2
fitted and mapped, outputs under `outputs/pod_jobs/dh270_lin_f0_full/outputs/two_stage/`).

- **Net transforms and hinges do not help.** Every approach (OLS hinges on z,
  MLE hinges, encounter transforms) prefers near-step functions of net use,
  i.e. a time trend. The z GAM's "jump at 0.2-0.3" is an artefact of z's
  standardisation (low-net pairs are early, near-100% mortality, noisy); on
  the selection scale the misfit is about linear in net use. Population:
  keep d_half = 270.
- **Net coefficient set by timing.** Within pixels and between sites at the
  same time, resistance does not rise with net use as the coefficient
  (~0.48) implies; the coefficient comes from the common timing of
  resistance and net scale-up.
- **Timing (w).** At w = 0.25 the model under-predicts resistance in
  2005-14 and over-predicts after 2018. A pseudo-likelihood screen (net
  coefficient and country levels refitted, other effects fixed;
  `R/net_w_screen.R`) improves to w ~ 3, but w > 1 is not biologically
  feasible (a cITN cannot select more than an LLIN); it signals missing
  early (pre-2010) pressure. Nick: constrain w to (0, 1]; w = 1 gives +90
  log-lik over 0.25. w does not change the late-period plateau maps.
- **Noise.** Study (citation) explains ~30% of the level residual, more than
  country or space; stage 2's pixel/pixel-year terms cover that. Species
  (arabiensis) is a partial, confounded signal, too sparse to model (left to
  stage 2). Cotton and other covariates do not explain the West/East split.
- **Plateaus.** Selection covariates can only change rates, never produce a
  plateau. The bimodal mortality floor (#14, #37) arises because plateau
  heights differ by place. Plateau heights do not follow local exposure.
- **The lead: kdr type.** The plateau map (change in level residual,
  2012-16 to 2019-25) aligns with the kdr surfaces: missing plateau
  strongest where Vgsc-995S dominates (Ethiopia ~ +2 logit, Kenya / Lake
  Victoria ~ +1), weakest where 995F is near fixation (coastal Cote
  d'Ivoire/Ghana neutral), blue where little kdr (SW Africa, Madagascar).
  Visual only so far. Figure: `figures/plateau_vs_kdr.png`.

## 3. Next: a latent genetic surface (Nick's direction)

Goal: a latent surface (or two) of resistance-mechanism type that we are
happy with, then fold it into the dynamical model, testing alternatives in
parallel on RunPod: as a modifier of selection strength, of the fitness cost
(reversion), and of the mortality floor/plateau.

1. **Quantify the lead first:** regress the plateau change (by bioassay
   location) on 995F, 995S frequency / share, restricted to areas with
   marker records.
2. **Data:** Vector Atlas genotypic records (`VA_FULL_DATA_20250721.csv`;
   995F 674 records, 995S 294; often only one variant typed per study; one
   record with longitude ~ -43000, filtered). Add MTM molecular-mechanism
   data and IR Mapper kdr data (Nick downloads; never scrape APIs). Consider
   Vgsc-1570Y and 402L (995F background); keep Ace-1 / metabolic markers
   for a possible second factor, not the first.
3. **Model for the surface:** joint model of variant frequencies with shared
   latent spatial factor(s), logit p_v(x,t) = alpha_v + beta_v (t - 2012) +
   lambda_v u(x) (v = variant; not m, which is the dynamical model). The
   single-factor version (`R/kdr_factor.R`) gave essentially 995F (lambda_S
   ~ 0 to -0.23); 995S needs its own factor or a two-factor model. Needs:
   overdispersion between studies (binomial AIC is invalid), shrinkage /
   low-rank smooth to avoid wild extrapolation (North Africa, Madagascar,
   southern Africa have no records), and a sensible prior on smoothness.
   Nick wants something simple and defensible, not a complex East African
   smooth.
4. **Model variants to fit on RunPod (later):** surface as (a) multiplier
   of net/pyrethroid selection, (b) driver of the fitness cost / reversion,
   (c) driver of the floor or plateau; each with w in (0, 1], floor
   handling as appropriate, d_half = 270. Launch with the grid launcher
   pattern (`R/net_grid.R`, `docker/README.md`); pods do NOT delete
   themselves; state prices; cpu5c 8 vCPU $0.28/h (~6 h per full fit),
   cpu3c 16 vCPU $0.48/h for stage 2.

## 4. Files (worktree `ir_cube_netscreen`, branch `net-transform-screen`)

Screens (all take the fit file and a label; outputs in
`outputs/hinge_screen/`, figures in `figures/`):

| Script | What |
|---|---|
| `R/change_residual_screen.R` | within-pixel pairs (min gap), z = qnorm(PIT), counts and yearly covariates |
| `R/change_misfit_plot.R` | misfit on the selection scale vs covariates |
| `R/hinge_select.R`, `R/hinge_select_plot.R` | OLS hinges on z (abandoned) |
| `R/hinge_mle.R`, `R/hinge_mle_plot.R` | concave hinge MLE on pairs (abandoned) |
| `R/transform_mle.R`, `R/transform_z_select.R` | encounter transforms for nets / d_half |
| `R/net_w_screen.R` | w pseudo-likelihood screen; saves logits for mock-ups |
| `R/between_within_nets.R` | level vs between-site net slopes |
| `R/species_misfit.R`, `R/cotton_misfit.R`, `R/misfit_regression.R` | candidate explanations |
| `R/misfit_maps.R`, `R/period_maps.R`, `R/plateau_map.R` | spatial misfit maps (optionally from mock-up logits) |
| `R/kdr_share.R`, `R/kdr_factor.R` | kdr surfaces (share; single latent factor) |

Superseded: `R/concave_hinge_select.R` (NNLS, buggy, abandoned).

Key outputs: `figures/plateau_map_dhalf270f0_w1.png`,
`figures/plateau_vs_kdr.png`, `figures/misfit_maps_dhalf270f0*.png`,
`figures/period_maps_dhalf270f0_2018_2025.png`,
`outputs/hinge_screen/dhalf270f0_net_w*.csv`,
`outputs/hinge_screen/dhalf270f0_misfit_regression*.txt`,
`outputs/hinge_screen/kdr_995s_share.tif`, `kdr_factor.tif`.

Sign conventions: level residual = empirical - predicted logit mortality
(positive: fit predicts too much resistance). Change misfit in
`change_misfit_plot.R` = observed - predicted fall (negative: over-predicted
selection); in `misfit_maps.R` it is flipped so red = over-prediction.
