# Microgravity-induced atrial remodelling increases re-entrant vulnerability: an in-silico human-atrial study

**Laura Piñero Roig**¹
¹ Faculty of Medicine, University of Barcelona, Barcelona, Spain · ORCID 0009-0008-3390-4029

*Preprint draft — not peer reviewed. Figures and ensemble statistics are generated reproducibly from the accompanying code (github.com/laurapiro17/atrial-fib-in-microgravity).*

---

## Abstract

**Background.** Long-duration spaceflight produces atrial structural and electrophysiological remodelling — chamber dilation from headward fluid shift, interstitial fibrosis, and shortened atrial refractoriness — each independently pro-arrhythmic on Earth. A clinical study of astronauts (Khine et al., 2018) reported left-atrial enlargement and electrophysiological changes that *may* raise atrial-fibrillation (AF) risk, though no sustained AF was observed. Whether these changes lower the threshold for re-entry has not been tested mechanistically.

**Methods.** We built a 2-D monodomain model of human atrial tissue using the Courtemanche–Ramirez–Nattel (CRN) ionic action-potential model, with fibre-anisotropic conduction calibrated to a physiological longitudinal conduction velocity (≈58 cm/s). Microgravity remodelling was represented as (i) AF-type ionic remodelling (reduced I_CaL, I_to, I_Kur; increased I_K1), which shortened the action-potential duration (APD₉₀) from ≈294 ms to ≈135 ms; (ii) spatially correlated interstitial fibrosis (low-coupling Gaussian-random-field patches); and (iii) chamber dilation. An S1–S2 cross-field protocol probed re-entrant vulnerability across an ensemble of fibrosis realisations; we quantified wavefront fragmentation as the burden of phase singularities (PS), with bootstrap confidence intervals.

**Results.** [TO FILL FROM ENSEMBLE: ground wavebreak burden = __ PS·ms (95% CI __–__); microgravity = __ PS·ms (95% CI __–__); seeds with wavebreak __/10 vs __/10; peak PS ground __ vs microgravity __.] Healthy ("ground") tissue conducted the induced wavefronts cleanly, with negligible phase-singularity formation, whereas the microgravity-remodelled fibrotic substrate fragmented propagating wavefronts into transient phase singularities (rotor cores) — the substrate of re-entry.

**Conclusions.** In this minimal but physiologically grounded model, microgravity-type atrial remodelling increases susceptibility to wavefront fragmentation and transient re-entry without producing sustained fibrillation — consistent with the clinical observation of increased electrophysiological vulnerability but no overt AF in astronauts. The model is a transparent, reproducible, hypothesis-generating tool for spaceflight cardiac-arrhythmia risk.

---

## 1. Introduction

Cardiovascular deconditioning is a recognised constraint on long-duration human spaceflight. In microgravity the heart undergoes documented remodelling: a headward fluid shift raises atrial filling pressure and dilates the chamber; deconditioning and altered loading promote structural remodelling and interstitial fibrosis; and autonomic shifts shorten atrial refractoriness. On Earth, each of these is independently pro-arrhythmic and is implicated in the substrate for atrial fibrillation (AF). Khine et al. (2018) reported, in astronauts after long-duration missions, transient increases in left-atrial size and changes in atrial electrophysiology consistent with an increased AF risk, although no episodes of sustained AF were captured.

The question this raises is mechanistic: does the *combination* of spaceflight-induced remodelling lower the threshold for, and sustain, fibrillatory re-entry? Direct experimental access is limited — astronaut numbers are small and invasive electrophysiology in flight is impractical. In-silico electrophysiology offers a complementary route: a tissue model lets each remodelling mechanism be toggled and the arrhythmic consequence measured.

Here we test, in a physiologically grounded 2-D human-atrial model, whether microgravity-type remodelling increases re-entrant vulnerability, and we characterise the *nature* of that vulnerability.

## 2. Methods

### 2.1 Ionic model

Membrane kinetics use the Courtemanche–Ramirez–Nattel (CRN) human atrial action-potential model — 21 state variables and 12 sarcolemmal currents plus sarcoplasmic-reticulum calcium handling. Gating variables are integrated with the Rush–Larsen exponential scheme; membrane voltage and ionic concentrations with forward Euler. The single-cell model was validated against published CRN biomarkers: resting potential ≈ −81 mV, peak ≈ +30 mV, APD₉₀ ≈ 290–300 ms at 1 Hz, maximum upstroke velocity > 100 V/s (≈220 V/s in the validation run), with physiological rate-dependent APD shortening.

### 2.2 Tissue model

Tissue is a 2-D monodomain sheet, `∂V/∂t = ∇·(D∇V) − I_ion/C_m`, integrated by operator splitting (Rush–Larsen reaction step, then an explicit finite-volume diffusion step under the diffusion-stability bound). Conduction is **orthotropic** with longitudinal:transverse ratio ≈ 3:1 (fibres aligned to the x-axis; a spatially rotating fibre field is supported by the implementation but not used here). The longitudinal diffusion coefficient was calibrated so that planar conduction velocity ≈ 58 cm/s (physiological human atrium). Space step dx = 0.25 mm; no-flux boundaries. The implementation is pure NumPy with an optional Numba-accelerated kernel (≈9× speed-up, validated identical to the reference to 1×10⁻¹³).

### 2.3 Microgravity remodelling

Three documented spaceflight changes are mapped to the substrate, each individually toggleable:

| Spaceflight change | Mechanism | Model representation |
|---|---|---|
| Shortened refractoriness | autonomic / AF-type electrical remodelling | I_CaL ↓, I_to ↓, I_Kur ↓, I_K1 ↑ → APD₉₀ ≈ 135 ms |
| Interstitial fibrosis | deconditioning / altered loading | correlated low-coupling Gaussian-random-field patches |
| Atrial dilation | headward fluid shift ↑ filling | enlarged sheet |

The ionic-remodelling pattern is the canonical AF-remodelling direction; magnitudes are scaled hypotheses chosen to reproduce the *direction* (shortened refractoriness) reported post-spaceflight, not direct microgravity tissue measurements. The dilation magnitude is anchored to Khine et al. (2018), who measured a transient left-atrial volume increase of 12 ± 18 mL (≈20%) after six months in space; the modelled enlargement is an illustrative bound, and — given that the standard deviation exceeds the mean — dilation is expected to be a minor contributor, which the mechanism-isolation analysis (Section 3.2) tests directly.

### 2.4 Re-entry induction and read-out

Re-entry was probed with an S1–S2 cross-field protocol: a planar S1 wavefront followed by an S2 stimulus in a quadrant, timed to break the recovering wavefront. Phase singularities (PS) — rotor cores — were counted over time using the topological-charge method on a phase field derived from V and a time-delayed copy of V. The sharp stimulus edges briefly produce a large spurious PS count (a wavefront artifact); these samples (PS above an artifact threshold, and the initial window) were excluded. The headline observable is the **wavebreak burden**: the artifact-rejected, time-integrated phase-singularity count (PS·ms), reflecting how strongly the substrate fragments propagating wavefronts. For each condition (ground, microgravity) we ran an ensemble of independent fibrosis realisations and report the mean with a 2000-sample bootstrap 95% confidence interval, plus the peak PS and the fraction of seeds exhibiting any wavebreak.

### 2.5 Reproducibility

All results derive from `experiments/ensemble_with_ci.py` at the accompanying repository; the model and metrics are covered by an automated test suite (single-cell biomarkers, diffusion operator, restitution monotonicity, fibrosis density, metric correctness) run in continuous integration.

## 3. Results

### 3.1 Single-cell and tissue validation

The CRN cell reproduced human atrial action-potential biomarkers (Section 2.1). Tissue conduction velocity was calibrated to ≈58 cm/s; APD restitution was monotonic and shortened with the AF-type remodelling. [Figure 1: single-cell AP and restitution curves.]

### 3.2 Re-entrant vulnerability

[TO FILL FROM `figures/results_crn.json["ensemble"]`:]
- Ground wavebreak burden: __ PS·ms (95% CI __–__); peak PS __; seeds with wavebreak __/10.
- Microgravity wavebreak burden: __ PS·ms (95% CI __–__); peak PS __; seeds with wavebreak __/10.
- Fold change / qualitative statement: __.

Healthy tissue conducted the induced wavefronts cleanly (negligible phase-singularity formation). The microgravity-remodelled fibrotic substrate fragmented wavefronts into transient phase singularities. [Figure 2: PS(t) ground vs microgravity with bootstrap band. Figure 3: late-state V snapshots — quiescent ground vs fragmented microgravity.]

### 3.3 Transience

In both conditions the induced activity self-terminated within a few hundred milliseconds; no sustained fibrillation occurred at the simulated tissue sizes. This is the expected consequence of the long human-atrial wavelength (conduction velocity × APD) relative to a laptop-feasible domain, and matches the clinical picture of increased vulnerability without observed sustained AF.

## 4. Discussion

The model supports a graded, mechanistic interpretation of spaceflight atrial risk: microgravity-type remodelling — particularly the combination of shortened refractoriness and interstitial fibrosis — converts a substrate that conducts cleanly into one that fragments wavefronts into transient rotors. This is the recognised first step toward re-entrant arrhythmia. The result aligns with Khine et al. (2018): electrophysiological vulnerability increases, but sustained AF is not an inevitable consequence.

### 4.1 Limitations

This is a hypothesis-generating model, not a quantitative clinical prediction. (i) It is a 2-D monodomain sheet — no realistic 3-D atrial geometry, no bidomain effects, no pulmonary-vein triggers. (ii) The microgravity→ionic mappings are directional hypotheses, not measured microgravity tissue relationships. (iii) Sustained fibrillation is not reproduced, partly because physiological human-atrial wavelengths approach chamber dimensions and a domain large enough to host multi-wavelet AF exceeds laptop compute; the claim is restricted to *vulnerability*, not sustained AF. (iv) Fibrosis is modelled as a low-coupling field, not discrete microstructure. References should be independently verified before formal citation.

### 4.2 Future work

Fit to additional human-atrial data; 3-D and bidomain extensions; larger GPU-scale ensembles to test for sustained multi-wavelet dynamics; coupling to measured microgravity / bedrest / dry-immersion cardiac datasets.

## Appendix A — Parameter table

**Ionic model (CRN baseline maximal conductances / fluxes).** Default Courtemanche–Ramirez–Nattel (1998) values:

| Current | Symbol | Value | Unit |
|---|---|---|---|
| Fast Na⁺ | g_Na | 7.8 | nS/pF |
| Inward rectifier K⁺ | g_K1 | 0.09 | nS/pF |
| Transient outward K⁺ | g_to | 0.1652 | nS/pF |
| Ultrarapid delayed K⁺ | g_Kur | 0.005 + 0.05/(1+e^(−(V−15)/13)) | nS/pF |
| Rapid delayed K⁺ | g_Kr | 0.0294 | nS/pF |
| Slow delayed K⁺ | g_Ks | 0.129 | nS/pF |
| L-type Ca²⁺ | g_CaL | 0.12375 | nS/pF |
| Na⁺–K⁺ pump | I_NaK,max | 0.599 | pA/pF |
| Na⁺–Ca²⁺ exchanger | I_NaCa,max | 1600 | pA/pF |
| Membrane capacitance | C_m | 100 | pF |

**Microgravity (AF-type) remodelling — scalings applied at severity = 1.0:**

| Current | Change | Rationale |
|---|---|---|
| I_CaL | −70 % | canonical AF electrical remodelling (APD/refractoriness ↓) |
| I_to | −50 % | canonical AF electrical remodelling |
| I_Kur | −50 % | canonical AF electrical remodelling |
| I_K1 | +100 % | canonical AF electrical remodelling (resting stabilisation, APD ↓) |
| I_Na, I_Kr | unchanged | conduction velocity preserved |

Net effect: APD₉₀ ≈ 294 ms (baseline) → ≈ 135 ms (remodelled). `severity = 0` recovers baseline; scalings are linear in severity.

**Tissue / numerics:**

| Parameter | Value |
|---|---|
| Space step dx | 0.25 mm |
| Longitudinal diffusion D_∥ | 0.15 mm²/ms (planar CV ≈ 58 cm/s) |
| Transverse diffusion D_⊥ | 0.05 mm²/ms (anisotropy ≈ 3:1) |
| Time step dt | 0.02 ms (operator-split; Rush–Larsen gates) |
| Dilation factor (microgravity) | 1.3× linear |
| Fibrosis area fraction | 0.30 |
| Fibrosis correlation length | 4 cells (1.0 mm) |
| Scar coupling | 0.05 × healthy |
| Boundary conditions | no-flux (Neumann) |

*All values are transcribed from the cited CRN source or are the calibrated/illustrative choices documented above; remodelling magnitudes are directional hypotheses, not microgravity tissue measurements.*

## References

1. Courtemanche M, Ramirez RJ, Nattel S. Ionic mechanisms underlying human atrial action potential properties. *Am J Physiol.* 1998;275(1):H301–21.
2. Khine HW, Steding-Ehrenborg K, Hastings JL, et al. Effects of prolonged spaceflight on atrial size, atrial electrophysiology, and risk of atrial fibrillation. *Circ Arrhythm Electrophysiol.* 2018;11(5):e005959. (PMID 29752376.)
3. Rush S, Larsen H. A practical algorithm for solving dynamic membrane equations. *IEEE Trans Biomed Eng.* 1978;25(4):389–92.
4. Gray RA, Pertsov AM, Jalife J. Spatial and temporal organization during cardiac fibrillation. *Nature.* 1998;392:75–78.
5. Garrett-Bakelman FE, et al. The NASA Twins Study. *Science.* 2019;364:eaau8650.
