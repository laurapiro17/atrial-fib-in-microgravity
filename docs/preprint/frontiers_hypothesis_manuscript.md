<!--
Target: Frontiers in Physiology — section "Cardiac Electrophysiology"
Article type: Hypothesis and Theory
Formatting follows Frontiers conventions: unstructured abstract (<=350 words),
keywords, free-form numbered sections, and the required end-statements
(Data Availability, Author Contributions, Funding, Conflict of Interest).
Author-year citations; full reference list at the end.
NOTE: corroboration values are the full single-seed run; the multi-seed
bootstrap CI is being computed and will replace them before submission.
-->

# A gravitational scaling law for atrial re-entrant vulnerability: a dimensionless predictor of the critical gravity for spaceflight atrial fibrillation

**Laura Piñero Roig¹\***

¹ Medical student, Faculty of Medicine, University of Barcelona, Barcelona, Spain

**\* Correspondence:** Laura Piñero Roig · ORCID 0009-0008-3390-4029

**Running title:** A gravitational scaling law for the atrium

---

## Abstract

Long-duration spaceflight remodels the atrium — headward fluid shift dilates the chamber, and atrial electrophysiology changes in ways that raise markers of atrial-fibrillation (AF) risk, although sustained AF has not been observed in astronauts (Khine et al., 2018). As permanent habitation of the Moon (0.16 g) and Mars (0.38 g) is planned, a question becomes pressing and remains unanswered: is there a gravitational threshold below which the atrium acquires an arrhythmogenic substrate, and on which side of it do the Moon and Mars fall? Cardiac physiology has allometric scaling laws with body mass but none with gravitational acceleration. Here I propose that the atrium's proximity to sustainable re-entry is governed by a single dimensionless quantity, the Cardiac Gravitational Number, 𝒩_g = L₀·dilation(g)/WL(g), where WL = CV·APD₉₀ is the electrical wavelength, CV is conduction velocity, APD₉₀ the action-potential duration, L₀ the atrial path length, and g gravitational level in Earth-g units. Re-entry can be hosted when 𝒩_g ≳ 1, defining a critical gravity g* with 𝒩_g(g*) = 1. Using a validated human-atrial ionic model (Courtemanche–Ramirez–Nattel) under a continuous, literature-anchored map from gravity to remodelling, 𝒩_g rose monotonically as gravity fell (0.47 at Earth to 1.29 in free fall), crossing unity at g* ≈ 0.22 — between the Moon (𝒩_g ≈ 1.07, vulnerable) and Mars (𝒩_g ≈ 0.82, protected). An independent two-dimensional reaction–diffusion sheet reproduced the trend: rotor burden rose as gravity fell and was zero at Earth, with a gradual rather than sharp transition. g* was robust to the assumed map shape (0.13–0.34) but strongly dependent on atrial size: no vulnerability at L₀ = 6 cm, rising to vulnerability even at the Mars level for L₀ = 12 cm. The framework yields a concrete, falsifiable, pre-flight-measurable hypothesis — that individual atrial size may set off-world arrhythmia risk — testable by partial-gravity electrophysiology.

**Keywords:** atrial fibrillation, microgravity, space medicine, cardiac electrophysiology, re-entry, mechano-electric feedback, computational modelling, scaling law

---

## 1. Introduction

Cardiovascular adaptation is a recognised constraint on long-duration human spaceflight. In microgravity a headward fluid shift raises atrial filling and dilates the chamber; deconditioning and altered loading change atrial structure; and autonomic shifts alter atrial refractoriness. On Earth, each of these is independently implicated in the substrate for atrial fibrillation (AF). In astronauts after six-month missions, Khine et al. (2018) reported transient left-atrial enlargement (+12 ± 18 mL) and electrophysiological changes consistent with increased AF risk, although no episodes of sustained AF were captured — a prevalence of AF in active astronauts (~5%) similar to the general population but occurring at a younger age.

As humanity plans sustained presence on the Moon (0.16 g) and Mars (0.38 g), a question follows that current evidence does not answer: **is there a gravitational threshold below which the atrium acquires an arrhythmogenic substrate, and where do the Moon and Mars fall relative to it?** Existing models treat spaceflight remodelling as a binary contrast (1 g versus "microgravity") and are descriptive. Physiology offers powerful allometric scaling laws with body mass, but no scaling law with gravitational acceleration.

I propose to treat gravity as a continuous variable and to compress the question into a single dimensionless number, in the spirit of the governing numbers of fluid mechanics. I then ask whether an independent, spatially resolved reaction–diffusion model tells the same story, and I quantify how far the prediction can be trusted given its assumptions. This is a hypothesis-and-theory contribution: it proposes and tests a scaling argument in silico; it does not claim experimental proof.

## 2. The hypothesis

The mechanistic chain that motivates the hypothesis is built from established physiology, link by link:

1. **Gravity → hydrostatic gradient.** Upright on Earth a vertical blood column ΔP = ρgh exists; as g → 0 it vanishes, driving a cephalad fluid shift toward the thorax.
2. **Fluid shift → atrial stretch.** Increased preload raises atrial volume; by the law of Laplace (σ = P·r/2t) wall stress and myocyte stretch rise — observed clinically as left-atrial enlargement in spaceflight (Khine et al., 2018).
3. **Stretch → mechano-electric feedback.** Atrial stretch shortens the refractory period and/or slows conduction velocity and increases their spatial dispersion (Ravelli, 2003).
4. **Wavelength shortening → re-entry.** The electrical wavelength WL = CV·APD₉₀ shortens; when it no longer fits the available path length, wavefronts fragment into re-entrant rotors.

From this chain I define the **Cardiac Gravitational Number**:

> 𝒩_g = L(g) / WL(g) = L₀·dilation(g) / (CV · APD₉₀(g)/1000),

with L₀ the ground atrial characteristic path length (human left-atrial circumference ≈ 8–12 cm; nominal 8 cm) and CV ≈ 58 cm/s the calibrated planar conduction velocity. As g → 0, the fluid-shift drive (1 − g) increases remodelling, shortening APD₉₀ (hence WL) and increasing dilation, so 𝒩_g rises. **Re-entry can be hosted when 𝒩_g ≳ 1**, defining the **critical gravity g\*** as the solution of 𝒩_g(g\*) = 1.

**Hypothesis statement.** There exists a critical gravity g\* for atrial re-entrant vulnerability; it lies between lunar and Martian gravity, such that interplanetary transit and the lunar surface fall in the vulnerable regime while Mars is relatively protected; and its value is governed primarily by individual atrial size.

## 3. Methods

**Continuous gravity → remodelling map.** A deliberately minimal linear map sends g to atrial remodelling through the fluid-shift drive (1 − g): electrical remodelling severity, chamber dilation, and interstitial fibrosis density each scale linearly with (1 − g), calibrated so that g = 1 recovers a validated ground baseline and g = 0 recovers an established microgravity operating point. The linear form and the mechano-electric sensitivities are assumptions, examined in §4.3.

**Ionic model and wavelength.** APD₉₀ was measured with the Courtemanche–Ramirez–Nattel (CRN) human-atrial action-potential model under 1 Hz pacing to steady state. AF-type electrical remodelling (reduced I_CaL, I_to, I_Kur; increased I_K1) scaled with severity. All severities were paced in a single vectorised run as columns of a (1, K) sheet with per-column conductances; the g = 0 column reproduced the scalar baseline APD₉₀ = 293.8 ms exactly.

**Independent full-sheet corroboration.** A two-dimensional monodomain reaction–diffusion sheet (CRN ionics; fibre-anisotropic conduction calibrated to ≈58 cm/s; spatially correlated interstitial fibrosis) was constructed at each gravity level, with severity, fibrosis density and dilation all driven by the same gravity → remodelling map. A broken wavefront was induced and rotor burden quantified as phase-singularity density via the topological-charge method on a time-delayed-voltage phase field.

**Sensitivity.** Without re-simulating the ionic model, g\* was recomputed for L₀ ∈ {6, 8, 10, 12} cm and for fluid-shift drive exponents p ∈ {0.5, 1, 2} in (1 − g)^p, interpolating the measured APD₉₀(severity) relation.

All analyses are reproducible from the accompanying open repository (see Data Availability).

## 4. Results

### 4.1 The scaling law predicts a critical gravity between the Moon and Mars

𝒩_g rose monotonically as gravity fell (Table 1; Figure 1), crossing unity at **g\* ≈ 0.22**. The Moon (0.16 g) lies just inside the vulnerable region (𝒩_g ≈ 1.07); Mars (0.38 g) lies outside it (𝒩_g ≈ 0.82); interplanetary transit (0 g) is most vulnerable (𝒩_g = 1.29); Earth is safe (𝒩_g = 0.47).

**Table 1.** Cardiac Gravitational Number across gravity (L₀ = 8 cm, CV = 58 cm/s).

| g (Earth-g) | body | APD₉₀ (ms) | WL (cm) | 𝒩_g |
|:---:|:---|:---:|:---:|:---:|
| 0.00 | interplanetary transit | 139.3 | 8.08 | 1.29 |
| 0.16 | Moon | ≈163 | ≈9.4 | ≈1.07 |
| 0.38 | Mars | ≈202 | ≈11.7 | ≈0.82 |
| 0.50 | — | 219.8 | 12.75 | 0.72 |
| 1.00 | Earth | 294.1 | 17.06 | 0.47 |

### 4.2 An independent reaction–diffusion model corroborates the trend

The full sheet reproduced the direction of the law: rotor burden (phase-singularity density, ×10⁻⁴ cells) was zero at Earth and rose as gravity fell. The two models agree at the extremes (clean conduction at Earth; vulnerability at low gravity) and in the direction of change. They differ in sharpness: the analytic threshold is crisp at g\* ≈ 0.22, whereas the sheet shows a **gradual** transition, with residual rotors persisting at the Mars level — that is, the single-cell wavelength argument is somewhat more optimistic about Mars than the spatially resolved model.

### 4.3 g* is robust to the map shape but governed by atrial size

g\* was robust to the assumed fluid-shift map shape (g\* ∈ [0.13, 0.34] for p ∈ {0.5, 1, 2}). However, g\* depended strongly on the characteristic atrial path length L₀: there was no crossing (no vulnerability at any gravity) for L₀ = 6 cm, g\* = 0.22 at 8 cm, 0.41 at 10 cm, and 0.58 at 12 cm — at which point Mars itself becomes vulnerable (Table 2).

**Table 2.** Sensitivity of the critical gravity g\*.

| assumption | value | g\* |
|:---|:---:|:---:|
| L₀ = 6 cm | small atrium | none (always safe) |
| L₀ = 8 cm | nominal | 0.22 |
| L₀ = 10 cm | enlarged | 0.41 |
| L₀ = 12 cm | dilated | 0.58 (Mars vulnerable) |
| drive exponent p = 0.5 | concave | 0.34 |
| p = 1.0 | linear | 0.22 |
| p = 2.0 | convex | 0.13 |

## 5. Discussion

A one-line dimensionless argument predicts a critical gravity for atrial re-entrant vulnerability that falls between the Moon and Mars, and an independent reaction–diffusion model corroborates the trend. The prediction is consistent with the clinical picture — increased vulnerability without sustained AF (Khine et al., 2018) — because 𝒩_g indexes the substrate's *proximity* to re-entry, not the occurrence of fibrillation.

The most decision-relevant result is the strong dependence on atrial size. Because g\* moves from "no one is ever vulnerable" (L₀ = 6 cm) to "even Mars is vulnerable" (L₀ = 12 cm), the framework predicts that **individual atrial dimensions could determine off-world arrhythmia risk** — a personalised, pre-flight-measurable hypothesis, since atrial size is routinely imaged, and a sharper claim than any single population threshold.

### 5.1 Limitations

The g → stretch link is the weakest: the map uses acute-stretch mechano-electric sensitivities as a stand-in for chronic adaptation, and the linear form is a hypothesis, not data. The corroboration uses a single ionic model in two dimensions, without autonomic dynamics or three-dimensional atrial geometry, and is a direction-of-effect check rather than a calibrated incidence estimate. 𝒩_g is proposed, not proven; the analytic threshold is more optimistic about Mars than the full model. No human partial-gravity electrophysiology exists to fix the absolute value of g\*, so the headline g\* ≈ 0.22 is conditional on L₀ ≈ 8 cm and the stated calibration choices.

### 5.2 Falsifiability and next steps

The framework is falsifiable: if partial-gravity (centrifuge or analog) measurements of atrial APD₉₀, CV and effective path length yield an 𝒩_g that does not cross unity within [0, 1] g, or crosses it far from the predicted band, the hypothesis is wrong. Immediate computational next steps are a multi-seed ensemble corroboration with confidence bands and the addition of a dispersion-of-refractoriness term to 𝒩_g; the decisive experimental test is partial-gravity atrial electrophysiology, ideally stratified by atrial size.

## Figures

- **Figure 1.** 𝒩_g versus gravitational level g, with the vulnerable region (𝒩_g > 1), the critical gravity g\*, and markers for interplanetary transit, the Moon, Mars and Earth. (`figures/gravity_law.png`)
- **Figure 2.** Full-sheet rotor burden (phase-singularity density) versus g, with 𝒩_g overlaid. (`figures/gravity_corroboration.png`)
- **Figure 3.** Sensitivity of g\* to atrial path length L₀ and to the fluid-shift map shape. (`figures/gravity_sensitivity.png`)

## Data Availability Statement

All code, parameters, and scripts that reproduce every figure and number in this article are openly available at https://github.com/laurapiro17/atrial-fib-in-microgravity. No human or animal data were generated; the clinical figures cited are from the published literature.

## Author Contributions

LPR conceived the hypothesis, designed and implemented the model and analyses, produced the figures, and wrote the manuscript.

## Funding

The author declares that no specific grant was received for this work from any funding agency in the public, commercial, or not-for-profit sectors.

## Conflict of Interest

The author declares that the research was conducted in the absence of any commercial or financial relationships that could be construed as a potential conflict of interest.

## References

Courtemanche, M., Ramirez, R. J., and Nattel, S. (1998). Ionic mechanisms underlying human atrial action potential properties. *Am. J. Physiol.* 275 (1), H301–H321.

Khine, H. W., Steding-Ehrenborg, K., Hastings, J. L., Kowal, J., Daniels, J. D., Page, R. L., et al. (2018). Effects of prolonged spaceflight on atrial size, atrial electrophysiology, and risk of atrial fibrillation. *Circ. Arrhythm. Electrophysiol.* 11 (5), e005959. doi: 10.1161/CIRCEP.117.005959

Ravelli, F. (2003). Mechano-electric feedback and atrial fibrillation. *Prog. Biophys. Mol. Biol.* 82 (1–3), 137–149. doi: 10.1016/s0079-6107(03)00011-7
