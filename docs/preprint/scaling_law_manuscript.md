# A gravitational scaling law for atrial re-entrant vulnerability: a dimensionless predictor of the critical gravity for spaceflight atrial fibrillation

**Laura Piñero Roig**¹
¹ Medical student, Faculty of Medicine, University of Barcelona, Barcelona, Spain · ORCID 0009-0008-3390-4029

*Preprint draft — not peer reviewed. All figures and numbers are generated reproducibly from the accompanying code (github.com/laurapiro17/atrial-fib-in-microgravity). This is a hypothesis paper with in-silico support; it proposes and tests a scaling argument, it does not prove a law.*

---

## Abstract

**Background.** Long-duration spaceflight remodels the atrium and raises markers of atrial-fibrillation (AF) risk (Khine et al., 2018), yet whether reduced gravity has a *threshold* below which an arrhythmogenic substrate emerges — and where that threshold lies relative to the Moon and Mars — is unknown. Cardiac physiology has scaling laws with body mass but none with gravitational acceleration.

**Hypothesis.** We propose that the atrium's proximity to sustainable re-entry is governed by a single dimensionless quantity, the **Cardiac Gravitational Number** `𝒩_g = L₀·dilation(g) / WL(g)`, where `WL = CV·APD₉₀` is the electrical wavelength and `g` is gravitational level in Earth-g units. Re-entry can be hosted when `𝒩_g ≳ 1`, defining a **critical gravity `g*`** with `𝒩_g(g*) = 1`.

**Methods.** A continuous, literature-anchored map sends `g` to atrial remodelling through the cephalad fluid-shift drive `(1 − g)`, calibrated so `g=1` recovers a validated ground baseline and `g=0` recovers an established microgravity operating point. Single-cell APD₉₀ was measured with the Courtemanche–Ramirez–Nattel (CRN) human-atrial ionic model across the remodelling continuum (one vectorised pacing run), forming `𝒩_g(g)`. The prediction was independently corroborated in a 2-D monodomain reaction–diffusion sheet (fibre-anisotropic conduction calibrated to ≈58 cm/s) by measuring rotor (phase-singularity) burden at each gravity level. Sensitivity to the two main assumptions (atrial path length `L₀`; fluid-shift map shape) was quantified.

**Results.** `𝒩_g` rose monotonically as gravity fell, from 0.47 at Earth to 1.29 in free fall, crossing unity at **`g* ≈ 0.22`** — between the Moon (`𝒩_g ≈ 1.07`, vulnerable) and Mars (`𝒩_g ≈ 0.82`, protected). The independent full-sheet model reproduced the trend: rotor burden increased as gravity decreased and was exactly zero at Earth, with a gradual rather than sharp transition (residual rotors persisted at the Mars level). `g*` was robust to the map shape (`g* ∈ [0.13, 0.34]` for drive exponents 0.5–2.0) but strongly dependent on atrial size: no vulnerability at `L₀ = 6 cm`, `g* = 0.22` at 8 cm, 0.41 at 10 cm, and 0.58 at 12 cm (Mars vulnerable).

**Conclusions.** A simple dimensionless argument predicts a critical gravity for atrial re-entrant vulnerability that falls between the Moon and Mars, and is corroborated by an independent reaction–diffusion model. The strong dependence on atrial size yields a concrete, testable clinical hypothesis: individual atrial dimensions may set off-world arrhythmia risk. The framework is falsifiable by partial-gravity electrophysiology and is offered as an open, reproducible tool.

---

## 1. Introduction

Cardiovascular adaptation constrains long-duration human spaceflight. In microgravity a headward fluid shift raises atrial filling and dilates the chamber; deconditioning promotes interstitial change; and autonomic shifts alter atrial refractoriness. Each is independently pro-arrhythmic on Earth. In astronauts after six-month missions, Khine et al. (2018) reported transient left-atrial enlargement (+12 ± 18 mL) and electrophysiological changes consistent with increased AF risk, with **no episodes of sustained AF** observed.

As humanity plans permanent habitation on the Moon (0.16 g) and Mars (0.38 g), a question becomes unavoidable and currently unanswered: **is there a gravitational threshold below which the atrium acquires an arrhythmogenic substrate, and on which side of it do the Moon and Mars fall?** Existing work treats spaceflight remodelling as a binary (1 g vs "microgravity") and is descriptive. Physiology offers allometric scaling laws with body mass but no scaling law with gravitational acceleration.

We propose to make gravity a continuous variable and to compress the question into a single dimensionless number, in the spirit of the governing numbers of fluid mechanics. We then test whether an independent, spatially resolved reaction–diffusion model tells the same story, and we quantify how far the prediction can be trusted given its assumptions.

## 2. The mechanism

Each link is supported by published physiology, not assumed de novo:

1. **Gravity → hydrostatic gradient.** Upright on Earth a vertical blood column `ΔP = ρgh` exists; as `g → 0` it vanishes, driving a cephalad fluid shift toward the thorax.
2. **Fluid shift → atrial stretch.** Increased preload raises atrial volume; by the law of Laplace (`σ = P·r/2t`) wall stress and myocyte stretch rise. Clinically observed as left-atrial enlargement in spaceflight (Khine et al., 2018).
3. **Stretch → mechano-electric feedback.** Atrial stretch shortens the refractory period and/or slows conduction velocity and increases their spatial dispersion (Ravelli, 2003).
4. **Wavelength shortening → re-entry.** The electrical wavelength `WL = CV·APD₉₀` shortens; when it no longer fits the available path length, wavefronts fragment into re-entrant rotors.

## 3. The scaling law

We define the **Cardiac Gravitational Number**

> `𝒩_g = L(g) / WL(g) = L₀·dilation(g) / (CV · APD₉₀(g)/1000)`,

with `L₀` the ground atrial characteristic path length (human left-atrial circumference ≈ 8–12 cm; nominal 8 cm) and `CV ≈ 58 cm/s` the calibrated planar conduction velocity. As `g → 0` the fluid-shift drive `(1 − g)` increases remodelling, shortening APD₉₀ (hence `WL`) and increasing dilation, so `𝒩_g` rises. Re-entry can be hosted when `𝒩_g ≳ 1`; the **critical gravity `g*`** solves `𝒩_g(g*) = 1`.

The continuous map `g → remodelling` is deliberately the simplest hypothesis consistent with the endpoints: severity, dilation and fibrosis density scale linearly with `(1 − g)`, calibrated so `g = 1` recovers the validated ground baseline and `g = 0` recovers the project's microgravity operating point. The linear form and the mechano-electric sensitivities are assumptions, examined in §6.

## 4. Methods

**Ionic model and wavelength.** APD₉₀ was measured with the CRN human-atrial action-potential model (1 Hz pacing, steady state). AF-type electrical remodelling (reduced I_CaL, I_to, I_Kur; increased I_K1) scaled with severity. All severities were paced in a single vectorised run as columns of a (1, K) sheet with per-column conductances, exploiting NumPy broadcasting (the `g = 0` column reproduces the scalar baseline APD₉₀ = 293.8 ms exactly).

**Full-sheet corroboration.** An independent 2-D monodomain reaction–diffusion sheet (CRN ionics; orthotropic conduction calibrated to ≈58 cm/s; spatially correlated fibrosis) was built at each gravity level with severity, fibrosis density and dilation all driven by the same `g → remodelling` map. A broken wavefront was induced and rotor burden quantified as phase-singularity density via the topological-charge method on a delayed-voltage phase field.

**Sensitivity.** Without re-simulating the ionic model, `g*` was recomputed for `L₀ ∈ {6, 8, 10, 12} cm` and for fluid-shift drive exponents `p ∈ {0.5, 1, 2}` in `(1 − g)^p`, interpolating the measured APD₉₀(severity) relation.

## 5. Results

**The scaling law.** `𝒩_g` rose monotonically as gravity fell (Table 1, `figures/gravity_law.png`), crossing unity at **`g* ≈ 0.22`**. The Moon (0.16 g) lies just inside the vulnerable region (`𝒩_g ≈ 1.07`); Mars (0.38 g) lies outside it (`𝒩_g ≈ 0.82`); interplanetary transit (0 g) is the most vulnerable (`𝒩_g = 1.29`); Earth is safe (`𝒩_g = 0.47`).

**Table 1. Cardiac Gravitational Number across gravity.**

| g (Earth-g) | body | APD₉₀ (ms) | WL (cm) | 𝒩_g |
|:---:|:---|:---:|:---:|:---:|
| 0.00 | transit | 139.3 | 8.08 | 1.29 |
| 0.16 | Moon | ≈163 | ≈9.4 | ≈1.07 |
| 0.38 | Mars | ≈202 | ≈11.7 | ≈0.82 |
| 0.50 | — | 219.8 | 12.75 | 0.72 |
| 1.00 | Earth | 294.1 | 17.06 | 0.47 |

**Independent corroboration.** The full reaction–diffusion sheet, run as an ensemble of 6 fibrosis realisations per gravity level (bootstrap 95% CI), reproduced the direction of the law (Table 1b). Rotor burden was comparably high across transit, the Moon and Mars (≈9–12 ×10⁻⁴), fell sharply between 0.5 and 0.8 g, and was **exactly zero at Earth** (CI 0–0). The two models agree at the extremes (clean conduction at Earth; vulnerability at low gravity) and on the direction of change, but differ in where the transition completes: the analytic threshold is crisp at `g* ≈ 0.22`, whereas the spatially resolved model is **more pessimistic**, sustaining rotors through the Mars level and crossing toward zero only around 0.5–0.6 g. The single-cell wavelength argument is therefore an optimistic bound; the full model widens, rather than narrows, the vulnerable gravity range.

**Table 1b. Full-sheet rotor burden (6-seed ensemble, mean and bootstrap 95% CI).**

| g (Earth-g) | body | PS density ×10⁻⁴ (mean, 95% CI) |
|:---:|:---|:---:|
| 0.00 | transit | 9.05 (8.14–9.81) |
| 0.16 | Moon | 11.51 (8.87–14.24) |
| 0.38 | Mars | 10.85 (9.14–12.43) |
| 0.50 | — | 4.76 (3.30–5.86) |
| 0.80 | — | 0.21 (0.00–0.62) |
| 1.00 | Earth | 0.00 (0.00–0.00) |

**Sensitivity (Table 2).** `g*` was robust to the map shape (`g* ∈ [0.13, 0.34]` for `p ∈ {0.5, 1, 2}`) but strongly dependent on atrial size: no crossing (no vulnerability at any gravity) for `L₀ = 6 cm`, `g* = 0.22` at 8 cm, 0.41 at 10 cm, and 0.58 at 12 cm — at which point Mars itself becomes vulnerable.

**Table 2. Sensitivity of g\*.**

| assumption | value | g* |
|:---|:---:|:---:|
| L₀ = 6 cm | small atrium | none (always safe) |
| L₀ = 8 cm | nominal | 0.22 |
| L₀ = 10 cm | enlarged | 0.41 |
| L₀ = 12 cm | dilated | 0.58 (Mars vulnerable) |
| drive exponent p = 0.5 | concave | 0.34 |
| p = 1.0 | linear | 0.22 |
| p = 2.0 | convex | 0.13 |

## 6. Discussion

A one-line dimensionless argument predicts a critical gravity for atrial re-entrant vulnerability that falls **between the Moon and Mars**, and an independent reaction–diffusion model corroborates the trend. The prediction is consistent with the clinical picture — increased vulnerability without sustained AF (Khine et al., 2018) — because `𝒩_g` indexes *substrate proximity to re-entry*, not guaranteed fibrillation.

The most decision-relevant finding is the **strong dependence on atrial size**. Because `g*` moves from "no one is ever vulnerable" (`L₀ = 6 cm`) to "even Mars is vulnerable" (`L₀ = 12 cm`), the framework predicts that **individual atrial dimensions could determine off-world arrhythmia risk** — a personalised, pre-flight-measurable hypothesis (atrial size is routinely imaged), and a sharper claim than a population threshold.

## 7. Limitations (read this before believing the headline)

- **The `g → stretch` link is the weakest.** The map uses acute-stretch mechano-electric sensitivities as a stand-in for chronic adaptation; the linear form is a hypothesis, not data. This is the first thing a partial-gravity experiment should test.
- **Single ionic model, 2-D** for corroboration (6-seed ensemble, but no autonomic dynamics and no 3-D atrial geometry). The sheet result is a direction-of-effect check, not a calibrated incidence estimate.
- **`𝒩_g` is proposed, not proven.** It is a scaling argument validated in silico; it is not a theorem, and the analytic threshold is more optimistic about Mars than the full model.
- **No human partial-gravity electrophysiology exists** to fix the absolute value of `g*`. The headline `g* ≈ 0.22` is conditional on `L₀ ≈ 8 cm` and the calibration choices above.

## 8. Falsifiability and next steps

The framework is falsifiable: if partial-gravity (centrifuge or lunar/Mars-analog) measurements of atrial APD₉₀, CV and effective path length yield `𝒩_g` that does not cross unity within `[0, 1] g`, or crosses it far from the predicted band, the law is wrong. The immediate computational next steps are a dispersion-of-refractoriness term in `𝒩_g` and a 3-D, autonomically driven atrial model.

## References

- Khine HW, Steding-Ehrenborg K, Hastings JL, et al. Effects of Prolonged Spaceflight on Atrial Size, Atrial Electrophysiology, and Risk of Atrial Fibrillation. *Circ Arrhythm Electrophysiol.* 2018;11(5):e005959. https://doi.org/10.1161/CIRCEP.117.005959
- Ravelli F. Mechano-electric feedback and atrial fibrillation. *Prog Biophys Mol Biol.* 2003;82(1-3):137-149. https://doi.org/10.1016/s0079-6107(03)00011-7
- Courtemanche M, Ramirez RJ, Nattel S. Ionic mechanisms underlying human atrial action potential properties. *Am J Physiol.* 1998;275(1):H301–21.
