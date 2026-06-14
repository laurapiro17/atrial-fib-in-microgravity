# Design spec — A gravitational scaling law of atrial arrhythmogenesis

**Date:** 2026-06-15
**Author:** Laura Piñero Roig
**Repo:** `atrial-fib-in-microgravity`
**Status:** design approved, pending implementation plan

---

## 0. One-line pitch

There may exist a **dimensionless number** that governs how close the human atrium
is to sustaining re-entry as a function of gravitational acceleration `g`, and a
**critical gravity `g*`** below which an arrhythmogenic substrate emerges. We propose
it from first principles, validate it in-silico across the `0 g → 1 g` continuum with
an existing, tested reaction–diffusion atrial model, and state a falsifiable prediction
about Moon vs Mars vs interplanetary transit.

This turns cardiac space medicine from *descriptive* ("astronauts show risk markers")
into *predictive* ("a physical law says at which gravity the substrate appears").

---

## 1. Why this is worth doing

- **Gap is real.** Allometric physiology has scaling laws with body mass (Kleiber).
  There is **no published scaling law with gravity** for cardiac electrophysiology.
- **It fits a rare skill intersection** (mathematics + physics + medicine + an existing
  validated cardiac reaction–diffusion model). Few cardiologists can derive a
  dimensionless governing number; few physicists have a calibrated atrial model.
- **It is feasible now.** The repo already contains the machinery (see §6); the novelty
  is making gravity a *continuous* variable and framing the result as a law.

**Honest scope.** This is a *hypothesis paper with in-silico support*, not a proof and
not a Nobel. The honest, defensible win is a falsifiable scaling framework + a clean
figure + a concrete prediction — a genuine, publishable contribution.

---

## 2. Central claim (falsifiable)

> Define a dimensionless **Cardiac Gravitational Number** `N_g`. As `g` decreases,
> headward fluid shift increases atrial stretch, which (via mechano-electric feedback)
> shortens the electrical wavelength `WL = CV · ERP`. Re-entry becomes sustainable when
> `N_g = L / WL(g) ≳ 1`. There exists a **critical gravity `g*`** solving `N_g(g*) = 1`.
> We predict `g*` lies such that **interplanetary transit (0 g) is below threshold**
> (vulnerable) while **Mars (0.38 g) is at or above threshold** (relatively protected);
> the Moon (0.16 g) is the decisive test case.

**How to falsify:** if the model, swept across `g` with literature-calibrated
sensitivities, shows no monotonic rise of vulnerability with decreasing `g`, or no
crossing of `N_g = 1` within `[0, 1] g`, the law is wrong. A future flight/centrifuge
experiment measuring atrial ERP and CV vs partial gravity would test `g*` directly.

---

## 3. Verified evidence base

All three pillars verified against PubMed on 2026-06-15 (not assumed):

| Pillar | Source | What it establishes | DOI |
|---|---|---|---|
| Atrial enlargement + AF risk markers in spaceflight | Khine et al. 2018, *Circ Arrhythm Electrophysiol* 11(5):e005959 | LA volume +12±18 mL after 6 months (transient); AF ~5% in astronauts at younger age; risk markers up, **no sustained AF** | https://doi.org/10.1161/CIRCEP.117.005959 |
| Mechano-electric feedback chain | Ravelli 2003, *Prog Biophys Mol Biol* 82(1-3):137-49 | Atrial stretch shortens refractory period and/or slows conduction velocity and increases spatial dispersion → re-entry | https://doi.org/10.1016/s0079-6107(03)00011-7 |
| Headward fluid shift in microgravity | spaceflight fluid-shift literature (candidate PMIDs 12434816, 9688754, 8828643, 11537424, 11537422 — to be read & one selected) | Cephalad fluid shift raises central filling | (to confirm during implementation) |

**Consistency check (important):** clinical data = increased *vulnerability markers*
without observed sustained AF. The repo already reproduces exactly this (transient
wavebreak, not sustained fibrillation). `N_g` predicts substrate vulnerability, **not**
guaranteed AF — fully consistent. No overclaim.

---

## 4. Mechanism B — the causal chain (each link literature-backed)

1. **Gravity → hydrostatic gradient.** Upright on Earth, a vertical column `ΔP = ρ g h`
   exists; as `g → 0` it vanishes → headward fluid shift (~litres toward thorax).
2. **Fluid → atrial stretch.** ↑ preload → ↑ atrial volume → by **Laplace**
   `σ = P · r / (2 t)` → ↑ wall stress → myocyte stretch.
3. **Stretch → mechano-electric feedback.** Stretch-activated channels depolarise
   resting potential → ↓ conduction velocity (CV); refractory period (ERP) shortens
   (Ravelli 2003).
4. **CV↓, ERP↓ → wavelength↓ → re-entry.** `WL = CV · ERP`. When `WL < L` (tissue
   path length), wavefronts fragment into rotors.

---

## 5. Law A — the derivation

**Cardiac Gravitational Number:**

```
N_g  =  L_atrium / WL(g)  =  L / ( CV(g) · ERP(g) )
```

- `L` = characteristic atrial path length (fixed geometry).
- `CV(g)`, `ERP(g)` = gravity-dependent via the §4 chain.

As `g → 0`: fluid shift ↑ → stretch ↑ → `CV↓, ERP↓` → `WL↓` → **`N_g ↑`**.
Re-entry sustains when `N_g ≳ 1`. Solve `N_g(g*) = 1` → **critical gravity `g*`** with
a confidence interval (from the existing ensemble bootstrap).

**The new modelling link.** Today the repo's remodelling is a binary toggle
(1 g vs "microgravity"). We replace it with a **continuous map `g → remodelling
magnitude`**:

- fluid shift `∝ (1 − g/g₀)` → ΔV atrium → (Laplace) → stretch `λ(g)`
- stretch → `CV(g)`, `ERP(g)` via mechano-electric **sensitivities taken from
  literature** (slope of CV and ERP vs % stretch). Where a numeric slope is not
  available, it is declared an **assumed parameter with a stated range**, swept in
  sensitivity analysis — never presented as measured fact.

Start simple (CV, ERP). Optional richer version adds **dispersion of refractoriness**
as a second term.

---

## 6. Computational validation plan (reuses existing machinery)

Existing assets to reuse:
- `src/afib_microgravity/` — CRN ionic model + monodomain (CV calibrated ≈ 58 cm/s) ✅
- `experiments/dose_response.py` — already sweeps a parameter ✅
- `experiments/apd_wavelength.py` — already measures wavelength ✅
- ensemble + 95% CI bootstrap; phase-singularity (PS) detector ✅

New work:
1. **`remodeling.py` extension** — add `gravity_to_remodeling(g)` returning the
   continuous parameter set (fluid shift → stretch → CV/ERP deltas).
2. **`experiments/gravity_sweep.py`** — sweep `g ∈ {0, 0.16, 0.38, 0.50, 0.80, 1.00}`,
   per `g` measure `WL(g)` and wavebreak burden (PS·ms), across the ensemble.
3. **Compute `N_g(g)`** and locate the `N_g = 1` crossing → `g*` with CI.
4. **Headline figure** — `N_g` vs `g` curve with error bars; vertical markers for
   Moon / Mars / transit; shaded "vulnerable" region where `N_g > 1`.
5. **Tests** — unit tests for `gravity_to_remodeling` monotonicity and limits
   (`g=1` recovers the validated ground baseline; `g=0` recovers existing microgravity
   case), keeping CI green.

---

## 7. Deliverables

- Extended `remodeling.py` + new `gravity_sweep.py` + tests (CI green).
- One headline figure (`N_g` vs `g`) + supporting panels.
- A short **hypothesis manuscript** framing `N_g` and `g*`, target venue e.g.
  *npj Microgravity*, *Frontiers in Physiology*, or a Letter — to decide later.
- README section + reproducibility (matches existing repo style).

---

## 8. Integrity safeguards (non-negotiable)

- **Zero fabricated citations.** Every physiological number traces to a verified
  source (checked via PubMed/arXiv) or is labelled an assumed parameter with range.
- **Derived vs calibrated vs assumed** are always distinguished in text and figure
  captions.
- **`N_g` is a proposed hypothesis, not a theorem.** Language: "we propose and test
  in-silico", never "we prove".
- **No overclaim on clinical outcome.** Substrate vulnerability ≠ guaranteed AF;
  state explicitly, consistent with Khine 2018.
- **Honest limitations section** mandatory (2-D monodomain, single ionic model,
  acute-stretch sensitivities used for chronic adaptation, no autonomic dynamics).

---

## 9. Honest risks & limitations

- The `g → stretch` map is the weakest link; acute-stretch MEF data may not equal
  chronic adaptation. Mitigated by sensitivity sweep + explicit labelling.
- `g*` may fall outside `[0,1]` or have a CI too wide to separate Moon from Mars — that
  is itself a publishable, honest negative/uncertain result.
- 2-D, single-cell-model, no fibre-level 3-D geometry — stated as limitation, not hidden.

---

## 10. Definition of done

- `gravity_sweep.py` runs end-to-end, CI green, ensemble with CI.
- `N_g` vs `g` figure produced; `g*` (or "no crossing") reported with CI.
- Manuscript draft with verified citations and an explicit limitations section.
- Every claim traceable to: first principles, a verified citation, or a labelled
  assumption.
