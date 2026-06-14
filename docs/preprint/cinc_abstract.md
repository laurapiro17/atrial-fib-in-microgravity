# In-Silico Atrial Re-entrant Vulnerability under Microgravity-Induced Remodelling

*Draft for Computing in Cardiology (CinC). Target length ~4 pages, 2-column IEEE format.
Numbers marked [TO FILL] are populated from the reproducible ensemble run.*

**Laura Piñero Roig**
Faculty of Medicine, University of Barcelona, Barcelona, Spain

## Abstract

Long-duration spaceflight remodels the atrium — chamber dilation, interstitial
fibrosis, and shortened refractoriness — and a clinical study (Khine et al., 2018)
found electrophysiological changes that may raise atrial-fibrillation (AF) risk,
without observed sustained AF. We tested, in silico, whether such remodelling
increases re-entrant vulnerability. A 2-D human-atrial monodomain model
(Courtemanche-Ramirez-Nattel ionics, orthotropic conduction calibrated to ≈58 cm/s)
was driven by an S1-S2 protocol across an ensemble of fibrosis realisations; we
quantified wavefront fragmentation as an artifact-rejected phase-singularity
(wavebreak) burden. Healthy tissue conducted cleanly; the remodelled substrate
fragmented wavefronts (burden [TO FILL] vs [TO FILL] PS·ms). A mechanism-isolation
analysis attributed the effect predominantly to [TO FILL], and the vulnerable
window widened from [TO FILL] to [TO FILL] ms. The model reproduces increased
vulnerability without sustained AF, consistent with the clinical picture, and is a
reproducible, open tool for spaceflight arrhythmia risk.

## 1. Introduction

Cardiovascular risk gates long-duration human spaceflight. Microgravity produces
documented atrial remodelling, each component independently pro-arrhythmic on
Earth. Khine et al. reported transient left-atrial enlargement (12 ± 18 mL) and
electrophysiological change after six-month missions, raising — but not
confirming — AF risk. We ask mechanistically whether the combined remodelling
lowers the re-entry threshold, and which component dominates.

## 2. Methods

**Ionic model.** Courtemanche-Ramirez-Nattel human atrial action potential
(21 states, Rush-Larsen gating), validated against published biomarkers (resting
≈ −81 mV, APD₉₀ ≈ 290–300 ms at 1 Hz, upstroke > 100 V/s).

**Tissue.** 2-D monodomain, operator-split (reaction then explicit finite-volume
diffusion), orthotropic 3:1 conduction, dx = 0.25 mm, longitudinal diffusion
calibrated to planar CV ≈ 58 cm/s, no-flux boundaries. A Numba kernel gives ≈9×
speed-up, validated identical to the reference.

**Microgravity remodelling.** AF-type ionic scaling (I_CaL −70 %, I_to −50 %,
I_Kur −50 %, I_K1 +100 %; APD₉₀ → ≈135 ms), correlated Gaussian-random-field
fibrosis, and dilation (anchored to Khine 12 ± 18 mL).

**Vulnerability read-out.** S1-S2 cross-field induction; phase singularities (PS)
counted by the topological-charge method on a delayed-V phase field, with
artifact rejection (PS ≤ 40, t ≥ 40 ms). Headline = wavebreak burden (∫PS dt).
Detector specificity validated (single spiral → 1; planar/static step → 0).
Ensembles report bootstrap 95 % CIs. Mechanism-isolation, vulnerable-window,
restitution-slope and fibrosis-threshold analyses dissect the effect.

## 3. Results

- Single-cell and CV validation: [Figure 1].
- Wavebreak burden: ground [TO FILL] PS·ms (95 % CI [TO FILL]); microgravity
  [TO FILL] PS·ms (95 % CI [TO FILL]); seeds with wavebreak [TO FILL]/N.
- Mechanism isolation: dominant driver = [TO FILL] [Figure 2].
- Vulnerable-window width: ground [TO FILL] ms vs microgravity [TO FILL] ms.
- APD restitution max slope: ground [TO FILL] vs microgravity [TO FILL].
- Critical fibrosis density for wavebreak ≈ [TO FILL].
- Activity self-terminated within a few hundred ms; no sustained AF.

## 4. Discussion and Conclusions

Microgravity-type remodelling converts a cleanly-conducting substrate into one
that fragments wavefronts into transient rotors — the first step toward re-entry —
without producing sustained AF, matching Khine et al. The mechanism-isolation and
restitution-slope analyses provide the mechanistic link. Limitations: 2-D
monodomain, directional (not measured) remodelling magnitudes, no sustained-AF
regime at laptop-feasible domains. The model and analyses are open and reproducible
(github.com/laurapiro17/atrial-fib-in-microgravity).

## References

[1] Khine HW et al. Circ Arrhythm Electrophysiol 2018;11(5):e005959.
[2] Courtemanche M, Ramirez RJ, Nattel S. Am J Physiol 1998;275:H301-21.
[3] Rush S, Larsen H. IEEE Trans Biomed Eng 1978;25:389-92.
[4] Gray RA, Pertsov AM, Jalife J. Nature 1998;392:75-78.
