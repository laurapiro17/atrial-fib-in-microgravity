# Spec: Upgrade to a human-atrial ionic substrate (Courtemanche–Ramirez–Nattel)

**Date:** 2026-06-14
**Status:** Approved (design phase)
**Author:** Laura Piñero Roig

## Motivation

The current model uses the dimensionless two-variable Aliev–Panfilov (AP) kinetics.
It is fast and reproduces spiral re-entry, but it cannot make quantitative,
physiologically anchored claims: parameters are abstract, conduction velocity is in
"cells per unit time", and "microgravity remodelling" is a change to an opaque
parameter. To raise the project from a hypothesis-generating toy toward a
preprint-grade, ESA-credible computational study, we replace the substrate with a
**validated human atrial ionic model** — Courtemanche, Ramirez & Nattel (1998, *Am J
Physiol* 275:H301) — while keeping the project tested, reproducible, and honest about
its limitations.

This unlocks claims in **physical units** (mV, ms, cm/s, mM), a mechanistic link via
**APD and CV restitution**, and remodelling mappings anchored to the ionic-current
changes documented in atrial electrical remodelling.

## Goals

1. A correct, vectorised CRN single-cell ionic model, validated against published
   action-potential biomarkers.
2. A model-agnostic 2D monodomain sheet supporting **anisotropic** (fibre-aligned)
   diffusion, integrated with operator splitting (Rush–Larsen for gates + explicit
   diffusion).
3. Spatially **correlated** fibrosis textures (Gaussian random field), replacing
   blocky random squares.
4. Microgravity remodelling expressed as **CRN ionic-current scalings** + geometry +
   fibrosis + anisotropy, each individually toggleable and literature-anchored.
5. Mechanistic read-outs: **APD restitution** and **CV restitution** (S1–S2), in
   addition to the existing phase-singularity metrics.
6. A statistically honest ensemble: n ≥ 25 substrate realisations per condition with
   **bootstrap confidence intervals**, plus dose-response sweeps.

## Non-goals (explicit honesty boundary)

- No 3D atrial geometry, no realistic atrial anatomy.
- Monodomain only (no bidomain / extracellular field).
- No pulmonary-vein ectopic triggers, no stretch-activated channels.
- Microgravity → ionic mappings remain hypotheses in the *direction* reported in the
  literature, not direct microgravity tissue measurements.
- We do not claim clinical prediction. The deliverable remains hypothesis-generating,
  now on a physically faithful substrate.

## Architecture

Each unit has one purpose, a clear interface, and is independently testable.

```
src/afib_microgravity/
  crn.py          # CRN ionic kinetics: state layout, 12 currents, Rush–Larsen gate
                  #   update, concentration ODEs. Fully vectorised (each state var is a
                  #   2D ndarray). Conductances in a dataclass (CRNParams) so remodelling
                  #   = scaling fields/scalars.
  cell_model.py   # Common interface: a cell model exposes (state init, reaction step
                  #   given V and dt -> new state + I_ion, V accessor). AP and CRN both
                  #   implement it so the sheet is model-agnostic.
  diffusion.py    # Anisotropic finite-volume div(D grad V) with a per-cell conductivity
                  #   tensor (Dxx, Dyy, Dxy) built from a fibre-angle field. No-flux
                  #   boundaries. Isotropic uniform field reduces to the current 5-point
                  #   Laplacian (regression-checked against existing operator).
  model.py        # AtrialSheet refactored to hold a cell_model + a diffusion operator.
                  #   step(): operator splitting — reaction (cell_model, Rush–Larsen) then
                  #   diffusion (explicit FV) under an enforced stability bound.
  fibrosis.py     # Gaussian-random-field correlated low-coupling field: smooth noise
                  #   (FFT or separable Gaussian blur of white noise) thresholded to a
                  #   target density and correlation length. Deterministic by seed.
  remodeling.py   # Spaceflight -> substrate. microgravity_params() returns CRN
                  #   conductance scalings (I_CaL down, I_to down, I_Kur down by
                  #   literature-anchored factors); make_condition() assembles geometry +
                  #   fibrosis + anisotropy + kinetics. Each mechanism toggleable.
  restitution.py  # S1–S2 protocol: APD90 restitution and CV restitution curves on a
                  #   1D strand (cheap), the mechanistic explanation for fragmentation.
  metrics.py      # Existing: phase singularities (topological charge), dominant
                  #   frequency. Added: PS density time series, sustained-AF classifier,
                  #   single-cell AP biomarkers (Vrest, Vpeak, APD90, dVdt_max).
experiments/
  validate_single_cell.py   # writes AP biomarkers + figure; asserts physiological range
  restitution_curves.py     # ground vs microgravity APD/CV restitution
  ensemble_with_ci.py       # n>=25 per condition, bootstrap CIs on PS density
  dose_response.py          # sweeps: dilation, fibrosis density, remodelling severity
tests/
  test_crn_singlecell.py    # AP biomarkers vs published CRN values (the keystone test)
  test_diffusion.py         # anisotropic operator: isotropic limit, conservation, axis CV
  test_model.py             # stability, planar propagation, operator splitting
  test_fibrosis.py          # density + correlation length of the random field
  test_metrics.py           # PS detector on a synthetic spiral; restitution monotonicity
  test_remodeling.py        # each toggle changes the intended current/field only
```

## Numerics

- **Units:** physical. V in mV, t in ms, currents in pA/pF, concentrations in mM.
  Space in mm; diffusion coefficient tuned so longitudinal CV ≈ 60 cm/s
  (≈ 0.06 cm/ms). dx ≈ 0.25 mm.
- **Time stepping:** operator splitting. Gating variables via **Rush–Larsen**
  (exponential integrator using steady-state + time constant) for stability at larger
  dt; membrane voltage and concentrations via forward Euler within the reaction step;
  diffusion via explicit finite volume under `dt <= dx^2 / (4 * Dmax)`.
- **Reaction dt** ≈ 0.02–0.05 ms (validated in the single-cell test); diffusion
  sub-stepping if the stability bound is tighter.
- Pure NumPy, vectorised over the grid. No global-environment changes — all work in an
  isolated virtual environment.

## Validation strategy (keystone)

The single largest risk is an incorrect CRN transcription. Mitigation: **physical
validation tests**, not just "it runs":

- Single-cell AP at 1 Hz pacing: resting Vm ≈ −81 mV, peak ≈ +20 mV,
  **APD90 ≈ 300 ms**, dV/dt_max in physiological range. (Published CRN benchmarks.)
- Planar CV ≈ 60 cm/s longitudinal after tuning D; anisotropy ratio ≈ 2–3× transverse.
- Restitution curves monotonic; remodelling shortens APD and steepens restitution.
- PS detector returns exactly 1 for a synthetic single spiral, ~0 for planar waves.

Numbers (APD90, CV, etc.) will be reported in `figures/results.json` and the README,
each traceable to the experiment that produced it.

## Remodelling mapping (literature-anchored)

| Spaceflight change | Mechanism | CRN representation |
|---|---|---|
| Shortened atrial refractoriness | autonomic / electrical remodelling | I_CaL ↓, I_to ↓, I_Kur ↓ (AF-type), toggleable scalings |
| Atrial dilation | headward fluid shift ↑ filling | larger sheet (more wavelengths) |
| Interstitial fibrosis | deconditioning / loading | correlated low-coupling Gaussian random field |
| Fibre structure | anatomical | anisotropic diffusion tensor |

Honesty: AF-type remodelling is used because it is the documented electrical-remodelling
pattern that shortens refractoriness — the *direction* Khine et al. (2018) report
post-spaceflight. Magnitudes are scaled hypotheses, not microgravity measurements.

## Compute budget

CRN 2D runs are ~30–40 min each (grid ~200×200, several seconds simulated). Ensemble
(n ≥ 25 × conditions) and dose-response sweeps run for hours / overnight in the isolated
venv; results are checkpointed to JSON and verified on completion.

## Phasing (for the implementation plan)

1. CRN single-cell + validation test (keystone — gate everything on this passing).
2. Anisotropic diffusion operator + isotropic-limit regression.
3. Model refactor (cell_model interface + operator splitting) + AP back-compat.
4. Correlated fibrosis field.
5. Restitution module + curves.
6. Remodelling mapping (CRN scalings) + toggles.
7. Extended metrics (PS density, sustained-AF classifier).
8. Experiments: single-cell, restitution, ensemble+CI, dose-response.
9. README rewrite with physical numbers, figures, GIFs; CI stays green.

## References (to verify before formal citation)

- Courtemanche M, Ramirez RJ, Nattel S. *Ionic mechanisms underlying human atrial
  action potential properties.* Am J Physiol. 1998;275(1):H301–21.
- Rush S, Larsen H. *A practical algorithm for solving dynamic membrane equations.*
  IEEE Trans Biomed Eng. 1978;25(4):389–92.
- Khine HW et al. *Effects of prolonged spaceflight on atrial size, atrial
  electrophysiology, and risk of atrial fibrillation.* Circ Arrhythm Electrophysiol.
  2018;11(5):e005959. (PMID 29752376 — verified.)
- Aliev RR, Panfilov AV. Chaos Solitons Fractals. 1996;7(3):293–301. (legacy model)
- Gray RA, Pertsov AM, Jalife J. Nature. 1998;392:75–78. (phase singularities)
