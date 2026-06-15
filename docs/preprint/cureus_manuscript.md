<!--
Target: Cureus Journal of Medical Science — Original Article
Cureus is free to publish (no APC) and indexed in PubMed/PMC.
Format: structured abstract (Background/Methods/Results/Conclusions), keywords,
Introduction, Materials And Methods, Results, Discussion, Conclusions, then the
mandatory "Additional Information" disclosures block, and numbered references in
Cureus/AMA style with DOIs.
-->


# A Gravitational Scaling Law for Atrial Re-Entrant Vulnerability: A Dimensionless Predictor of the Critical Gravity for Spaceflight Atrial Fibrillation

**Laura Piñero Roig ¹**

¹ Medical Student, Faculty of Medicine, University of Barcelona, Barcelona, ESP

**Corresponding author:** Laura Piñero Roig · ORCID 0009-0008-3390-4029

---

## Abstract

**Background:** Long-duration spaceflight remodels the atrium: a headward fluid shift dilates the chamber, and atrial electrophysiology changes in ways that raise markers of atrial fibrillation (AF) risk, although sustained AF has not been observed in astronauts. As permanent habitation of the Moon (0.16 g) and Mars (0.38 g) is planned, it is unknown whether a gravitational threshold exists below which the atrium acquires an arrhythmogenic substrate, and where the Moon and Mars lie relative to it. Cardiac physiology has scaling laws with body mass but none with gravitational acceleration.

**Methods:** Building on the classical wavelength hypothesis of re-entry, we recast the path-length-to-wavelength ratio as a gravity-dependent dimensionless quantity, the Cardiac Gravitational Number, 𝒩_g = L₀·dilation(g)/WL(g), where WL = CV·APD₉₀ is the electrical wavelength; re-entry can be hosted when 𝒩_g ≥ 1, defining a critical gravity g* with 𝒩_g(g*) = 1. The novelty is the gravity parametrisation, not the underlying criterion. A validated human-atrial ionic model (Courtemanche–Ramirez–Nattel) was driven by a continuous, literature-anchored map from gravitational level to atrial remodelling, and 𝒩_g was computed across 0–1 g. The prediction was independently corroborated in a two-dimensional reaction–diffusion sheet by measuring rotor (phase-singularity) burden, and sensitivity to atrial size and map shape was quantified.

**Results:** 𝒩_g rose monotonically as gravity fell (0.47 at Earth to 1.29 in free fall), crossing unity at g* ≈ 0.22 — between the Moon (𝒩_g ≈ 1.07, vulnerable) and Mars (𝒩_g ≈ 0.82, protected). The full sheet reproduced the trend, with rotor burden zero at Earth and rising as gravity fell (gradual, not sharp, transition). g* was robust to the assumed map shape (0.13–0.34) but strongly dependent on atrial size (no vulnerability at L₀ = 6 cm; Mars vulnerable at L₀ = 12 cm).

**Conclusions:** A simple dimensionless argument predicts a critical gravity for atrial re-entrant vulnerability between the Moon and Mars, corroborated by an independent model. Its strong dependence on atrial size yields a falsifiable, pre-flight-measurable hypothesis: individual atrial dimensions may set off-world arrhythmia risk.

**Keywords:** atrial fibrillation, microgravity, space medicine, cardiac electrophysiology, re-entry, mechano-electric feedback, in silico model, scaling law

---

## Introduction

Cardiovascular adaptation is a recognised constraint on long-duration human spaceflight, and spaceflight predisposes to cardiac arrhythmias through fluid redistribution, autonomic dysregulation, and radiation exposure [5]. In microgravity a headward fluid shift raises atrial filling and dilates the chamber; deconditioning and altered loading change atrial structure; and autonomic shifts alter atrial refractoriness. On Earth, each of these is independently implicated in the substrate for atrial fibrillation (AF). In astronauts after six-month missions, Khine et al. reported transient left-atrial enlargement (+12 ± 18 mL) and electrophysiological changes consistent with increased AF risk, although no sustained AF was captured; AF prevalence in active astronauts (~5%) resembles the general population but occurs at a younger age [1].

As humanity plans sustained presence on the Moon (0.16 g) and Mars (0.38 g), a question follows that current evidence does not answer: is there a gravitational threshold below which the atrium acquires an arrhythmogenic substrate, and where do the Moon and Mars fall relative to it? Existing models treat spaceflight remodelling as a binary contrast (1 g versus "microgravity") and are descriptive. Physiology offers powerful allometric scaling laws with body mass, but none with gravitational acceleration.

We treat gravity as a continuous variable and compress the question into a single dimensionless number, in the spirit of the governing numbers of fluid mechanics. We then test whether an independent, spatially resolved reaction–diffusion model agrees, and quantify how far the prediction can be trusted given its assumptions. This is a hypothesis-generating study: it proposes and tests a scaling argument in silico and does not claim experimental proof.

## Materials And Methods

**The mechanistic chain.** Each link rests on established physiology: (i) on Earth a vertical hydrostatic column ΔP = ρgh exists and vanishes as g → 0, driving a cephalad fluid shift; (ii) increased preload raises atrial volume and, by the law of Laplace (σ = P·r/2t), wall stress and myocyte stretch — clinically, left-atrial enlargement in spaceflight [1]; (iii) atrial stretch shortens the refractory period and/or slows conduction velocity (CV) and increases their dispersion via mechano-electric feedback [2]; (iv) the electrical wavelength WL = CV·APD₉₀ shortens until it no longer fits the available path length, and wavefronts fragment into re-entrant rotors.

**The Cardiac Gravitational Number.** We define 𝒩_g = L(g)/WL(g) = L₀·dilation(g)/(CV·APD₉₀(g)/1000), with L₀ the ground atrial path length (human left-atrial circumference ≈ 8–12 cm; nominal 8 cm) and CV ≈ 58 cm/s. Re-entry can be hosted when 𝒩_g ≥ 1; the critical gravity g* solves 𝒩_g(g*) = 1. The ratio of path length to wavelength is itself the classical wavelength hypothesis of re-entry: Rensma et al. measured a normal atrial wavelength of 14–18 cm and atrial-fibrillation induction below ≈7.8 cm [4]. The contribution here is not this criterion but its parametrisation by gravity and the resulting threshold g*.

**Continuous gravity → remodelling map.** A minimal linear map sends g to atrial remodelling through the fluid-shift drive (1 − g): electrical-remodelling severity, dilation, and fibrosis density each scale linearly with (1 − g), calibrated so g = 1 recovers a validated ground baseline and g = 0 an established microgravity operating point.

**Ionic model.** APD₉₀ was measured with the Courtemanche–Ramirez–Nattel (CRN) human-atrial action-potential model [3] at 1 Hz to steady state, with AF-type remodelling (reduced I_CaL, I_to, I_Kur; increased I_K1) scaling with severity. All severities were paced in one vectorised run; the g = 0 column reproduced the scalar baseline APD₉₀ = 293.8 ms exactly.

**Full-sheet corroboration.** A two-dimensional monodomain reaction–diffusion sheet (CRN ionics; fibre-anisotropic conduction calibrated to ≈58 cm/s; spatially correlated fibrosis) was built at each gravity level with all tissue knobs driven by the same map. A broken wavefront was induced and rotor burden quantified as phase-singularity density via the topological-charge method on a time-delayed-voltage phase field.

**Sensitivity.** Without re-simulating the ionic model, g* was recomputed for L₀ ∈ {6, 8, 10, 12} cm and fluid-shift drive exponents p ∈ {0.5, 1, 2} in (1 − g)^p. All code and parameters are openly available (see Additional Information).

## Results

𝒩_g rose monotonically as gravity fell (Table 1), crossing unity at g* ≈ 0.22. The Moon (0.16 g) lies just inside the vulnerable region (𝒩_g ≈ 1.07); Mars (0.38 g) lies outside it (𝒩_g ≈ 0.82); interplanetary transit (0 g) is most vulnerable (𝒩_g = 1.29); Earth is safe (𝒩_g = 0.47). Notably, although the model was not calibrated to wavelength data, its Earth wavelength (17.1 cm) falls within the normal range measured by Rensma et al. (14–18 cm) and its free-fall wavelength (8.1 cm) sits at their AF-induction threshold (≈7.8 cm) [4] — an unprompted, cross-species consistency check.

*Table 1: Cardiac Gravitational Number across gravitational level (L₀ = 8 cm, CV = 58 cm/s).*

| g (Earth-g) | Body | APD₉₀ (ms) | WL (cm) | 𝒩_g |
|:---:|:---|:---:|:---:|:---:|
| 0.00 | Interplanetary transit | 139.3 | 8.08 | 1.29 |
| 0.16 | Moon | ≈163 | ≈9.4 | ≈1.07 |
| 0.38 | Mars | ≈202 | ≈11.7 | ≈0.82 |
| 0.50 | — | 219.8 | 12.75 | 0.72 |
| 1.00 | Earth | 294.1 | 17.06 | 0.47 |

The full reaction–diffusion sheet (6-seed ensemble per gravity level, bootstrap 95% CI) reproduced the direction of the law (Table 1b): rotor burden was comparably high across transit, the Moon and Mars (≈9–12 ×10⁻⁴), fell sharply between 0.5 and 0.8 g, and was exactly zero at Earth (CI 0–0). The models agree at the extremes and in direction but differ in where the transition completes: the analytic threshold is crisp at g* ≈ 0.22, while the spatially resolved model is more pessimistic, sustaining rotors through the Mars level and approaching zero only near 0.5–0.6 g. The single-cell argument is thus an optimistic bound on the vulnerable gravity range.

*Table 1b: Full-sheet rotor burden (6-seed ensemble, mean and bootstrap 95% CI).*

| g (Earth-g) | Body | PS density ×10⁻⁴ (mean, 95% CI) |
|:---:|:---|:---:|
| 0.00 | Interplanetary transit | 9.05 (8.14–9.81) |
| 0.16 | Moon | 11.51 (8.87–14.24) |
| 0.38 | Mars | 10.85 (9.14–12.43) |
| 0.50 | — | 4.76 (3.30–5.86) |
| 0.80 | — | 0.21 (0.00–0.62) |
| 1.00 | Earth | 0.00 (0.00–0.00) |

g* was robust to the assumed fluid-shift map shape (g* ∈ [0.13, 0.34] for p ∈ {0.5, 1, 2}) but depended strongly on atrial path length: no crossing (no vulnerability at any gravity) for L₀ = 6 cm, g* = 0.22 at 8 cm, 0.41 at 10 cm, and 0.58 at 12 cm — at which point Mars itself becomes vulnerable (Table 2).

*Table 2: Sensitivity of the critical gravity g\*.*

| Assumption | Value | g\* |
|:---|:---:|:---:|
| L₀ = 6 cm | Small atrium | None (always safe) |
| L₀ = 8 cm | Nominal | 0.22 |
| L₀ = 10 cm | Enlarged | 0.41 |
| L₀ = 12 cm | Dilated | 0.58 (Mars vulnerable) |
| Drive exponent p = 0.5 | Concave | 0.34 |
| p = 1.0 | Linear | 0.22 |
| p = 2.0 | Convex | 0.13 |

## Discussion

A one-line dimensionless argument predicts a critical gravity for atrial re-entrant vulnerability between the Moon and Mars, and an independent reaction–diffusion model corroborates the trend. The prediction is consistent with the clinical picture — increased vulnerability without sustained AF [1] — because 𝒩_g indexes the substrate's proximity to re-entry, not the occurrence of fibrillation.

The most decision-relevant result is the strong dependence on atrial size: because g* moves from "no one is ever vulnerable" (L₀ = 6 cm) to "even Mars is vulnerable" (L₀ = 12 cm), the framework predicts that individual atrial dimensions could determine off-world arrhythmia risk — a personalised, pre-flight-measurable hypothesis, since atrial size is routinely imaged.

This study has clear limitations. The gravity → stretch link is the weakest: the map uses acute-stretch mechano-electric sensitivities as a stand-in for chronic adaptation, and the linear form is a hypothesis, not data. The corroboration uses a single ionic model in two dimensions, without autonomic dynamics or three-dimensional geometry, and is a direction-of-effect check rather than a calibrated incidence estimate. 𝒩_g is proposed, not proven; the analytic threshold is more optimistic about Mars than the full model. No human partial-gravity electrophysiology exists to fix the absolute value of g*, so the headline g* ≈ 0.22 is conditional on L₀ ≈ 8 cm and the stated calibration. The framework is nonetheless falsifiable: partial-gravity measurements of atrial APD₉₀, CV and effective path length that yield an 𝒩_g not crossing unity within [0, 1] g, or crossing it far from the predicted band, would refute the hypothesis.

## Conclusions

A dimensionless Cardiac Gravitational Number predicts a critical gravity for atrial re-entrant vulnerability that falls between the Moon and Mars and is corroborated by an independent reaction–diffusion model. The prediction depends primarily on atrial size, giving a concrete, falsifiable, pre-flight-measurable hypothesis: that individual atrial dimensions may set off-world arrhythmia risk. The decisive test is partial-gravity atrial electrophysiology, ideally stratified by atrial size.

## Additional Information

**Author Contributions:** LPR conceived the hypothesis, designed and implemented the model and analyses, produced the figures, and wrote the manuscript.

**Disclosures:** *Human subjects:* All authors have confirmed that this study did not involve human participants or tissue. *Animal subjects:* All authors have confirmed that this study did not involve animal subjects or tissue. *Conflicts of interest:* In compliance with the ICMJE uniform disclosure form, the author declares no conflicts of interest. *Financial relationships:* The author declares no financial relationships with any organizations that might have an interest in the submitted work. *Other relationships:* The author declares no other relationships or activities that could appear to have influenced the submitted work.

**Data Availability:** All code, parameters, and scripts that reproduce every figure and number are openly available at https://github.com/laurapiro17/atrial-fib-in-microgravity.

## References

1. Khine HW, Steding-Ehrenborg K, Hastings JL, et al.: Effects of prolonged spaceflight on atrial size, atrial electrophysiology, and risk of atrial fibrillation. Circ Arrhythm Electrophysiol. 2018, 11:e005959. 10.1161/CIRCEP.117.005959
2. Ravelli F: Mechano-electric feedback and atrial fibrillation. Prog Biophys Mol Biol. 2003, 82:137-149. 10.1016/s0079-6107(03)00011-7
3. Courtemanche M, Ramirez RJ, Nattel S: Ionic mechanisms underlying human atrial action potential properties: insights from a mathematical model. Am J Physiol. 1998, 275:H301-H321. 10.1152/ajpheart.1998.275.1.H301
4. Rensma PL, Allessie MA, Lammers WJ, Bonke FI, Schalij MJ: Length of excitation wave and susceptibility to reentrant atrial arrhythmias in normal conscious dogs. Circ Res. 1988, 62:395-410. 10.1161/01.res.62.2.395
5. Mircea AA, Pistritu DV, Fortner A, et al.: Space travel: the radiation and microgravity effects on the cardiovascular system. Int J Mol Sci. 2024, 25:11812. 10.3390/ijms252111812
