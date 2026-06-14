# Gravitational Scaling Law Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a continuous map from gravitational level `g` to atrial remodelling, a dimensionless Cardiac Gravitational Number `N_g`, and a sweep experiment that locates the critical gravity `g*` and produces the headline `N_g` vs `g` figure.

**Architecture:** A new pure module `gravity_law.py` holds the physics (g→remodelling map, wavelength, `N_g`, crossing interpolation) so it is unit-testable in isolation. A new experiment `gravity_sweep.py` drives the existing CRN single-cell / sheet machinery across `g ∈ {0, 0.16, 0.38, 0.5, 0.8, 1.0}`, feeds measured APD90 into `gravity_law`, and renders the figure. Endpoints are calibrated so `g=1` recovers the validated ground baseline and `g=0` recovers the project's existing microgravity operating point.

**Tech Stack:** Python, NumPy, Matplotlib, pytest, existing `afib_microgravity` package (CRN ionic model, monodomain sheet, phase-singularity metrics).

---

## File structure

- Create: `src/afib_microgravity/gravity_law.py` — pure physics of the scaling law.
- Create: `tests/test_gravity_law.py` — unit tests for the pure module.
- Create: `experiments/gravity_sweep.py` — the g-sweep driver + figure.
- Create: `tests/test_gravity_sweep_smoke.py` — smoke test that the experiment runs.
- Modify: `README.md` — add a "Gravitational scaling law" section.

Commit rule for this repo: commits are authored as Laura, **no Co-Authored-By trailer**.

---

### Task 1: The continuous g → remodelling map

**Files:**
- Create: `src/afib_microgravity/gravity_law.py`
- Test: `tests/test_gravity_law.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_gravity_law.py
import math
import pytest
from afib_microgravity.gravity_law import (
    gravity_to_remodeling, SEVERITY_MAX, DILATION_MAX, DENSITY_MAX,
)


def test_earth_recovers_baseline():
    r = gravity_to_remodeling(1.0)
    assert r.severity == pytest.approx(0.0)
    assert r.dilation == pytest.approx(1.0)
    assert r.density == pytest.approx(0.0)


def test_freefall_recovers_microgravity_operating_point():
    r = gravity_to_remodeling(0.0)
    assert r.severity == pytest.approx(SEVERITY_MAX)
    assert r.dilation == pytest.approx(DILATION_MAX)
    assert r.density == pytest.approx(DENSITY_MAX)


def test_severity_monotonic_decreasing_in_gravity():
    sev = [gravity_to_remodeling(g).severity for g in (0.0, 0.25, 0.5, 0.75, 1.0)]
    assert all(a >= b for a, b in zip(sev, sev[1:]))


def test_gravity_clamped_to_unit_interval():
    assert gravity_to_remodeling(1.5).severity == pytest.approx(0.0)
    assert gravity_to_remodeling(-0.3).severity == pytest.approx(SEVERITY_MAX)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd ~/dev/atrial-fib-in-microgravity && .venv/bin/pytest tests/test_gravity_law.py -q`
Expected: FAIL with `ModuleNotFoundError: No module named 'afib_microgravity.gravity_law'`

- [ ] **Step 3: Write minimal implementation**

```python
# src/afib_microgravity/gravity_law.py
"""Gravitational scaling law of atrial arrhythmogenesis.

A continuous map from gravitational level ``g`` (in Earth-g units: g=1.0 is
Earth, g=0 is free fall) to the microgravity remodelling knobs, plus the
dimensionless Cardiac Gravitational Number

    N_g = L(g) / WL(g) = L0 * dilation(g) / (CV * APD90(g) / 1000)

Re-entry can be hosted when N_g >= 1. The critical gravity g* solves
N_g(g*) = 1.

HONESTY: the g -> remodelling map is a deliberately simple *linear* hypothesis
in the cephalad fluid-shift drive (1 - g), calibrated so g=1 recovers the
validated ground baseline and g=0 recovers the project's existing microgravity
operating point (see remodeling.py). The linear form and the mechano-electric
sensitivities are assumptions swept in sensitivity analysis, not measured
facts. Evidence base: Khine 2018 (10.1161/CIRCEP.117.005959, atrial
enlargement + AF risk markers in spaceflight) and Ravelli 2003
(10.1016/s0079-6107(03)00011-7, stretch shortens refractoriness / slows CV).
"""
from __future__ import annotations

from dataclasses import dataclass

EARTH_G = 1.0
# Endpoints recovered at g=0, taken from remodeling.make_condition_crn("microgravity").
SEVERITY_MAX = 1.0
DILATION_MAX = 1.3
DENSITY_MAX = 0.3


@dataclass(frozen=True)
class GravityRemodeling:
    gravity: float
    severity: float
    dilation: float
    density: float


def gravity_to_remodeling(g: float) -> GravityRemodeling:
    """Map gravitational level ``g`` to remodelling knobs.

    Linear in the fluid-shift drive ``(1 - g)``; ``g`` clamped to [0, 1].
    """
    g = min(max(g, 0.0), 1.0)
    drive = 1.0 - g
    return GravityRemodeling(
        gravity=g,
        severity=SEVERITY_MAX * drive,
        dilation=1.0 + (DILATION_MAX - 1.0) * drive,
        density=DENSITY_MAX * drive,
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/pytest tests/test_gravity_law.py -q`
Expected: PASS (4 passed)

- [ ] **Step 5: Commit**

```bash
git add src/afib_microgravity/gravity_law.py tests/test_gravity_law.py
git commit -m "feat(gravity_law): continuous g->remodelling map"
```

---

### Task 2: Wavelength, N_g, and the critical-gravity crossing

**Files:**
- Modify: `src/afib_microgravity/gravity_law.py`
- Test: `tests/test_gravity_law.py`

- [ ] **Step 1: Write the failing test (append to tests/test_gravity_law.py)**

```python
from afib_microgravity.gravity_law import (
    wavelength_cm, cardiac_gravitational_number, interpolate_crossing,
)


def test_wavelength_is_cv_times_apd_in_cm():
    # 58 cm/s * 0.300 s = 17.4 cm
    assert wavelength_cm(300.0, cv_cm_s=58.0) == pytest.approx(17.4)


def test_Ng_is_one_when_path_equals_wavelength():
    # WL = 58*0.1379... ; choose apd so wl=8 cm, l0=8, dilation=1 -> N_g=1
    apd = 8.0 / 58.0 * 1000.0
    ng = cardiac_gravitational_number(apd, dilation=1.0, l0_cm=8.0, cv_cm_s=58.0)
    assert ng == pytest.approx(1.0)


def test_Ng_rises_when_wavelength_shortens():
    long_wl = cardiac_gravitational_number(300.0, dilation=1.0)
    short_wl = cardiac_gravitational_number(150.0, dilation=1.0)
    assert short_wl > long_wl


def test_interpolate_crossing_finds_linear_root():
    xs = [0.0, 0.2, 0.4, 0.6]
    ys = [2.0, 1.5, 0.5, 0.2]   # crosses 1.0 between 0.2 and 0.4
    gstar = interpolate_crossing(xs, ys, target=1.0)
    assert gstar == pytest.approx(0.3)


def test_interpolate_crossing_returns_none_without_crossing():
    assert interpolate_crossing([0.0, 1.0], [3.0, 2.0], target=1.0) is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `.venv/bin/pytest tests/test_gravity_law.py -q`
Expected: FAIL with `ImportError: cannot import name 'wavelength_cm'`

- [ ] **Step 3: Write minimal implementation (append to gravity_law.py)**

```python
def wavelength_cm(apd90_ms: float, cv_cm_s: float = 58.0) -> float:
    """Electrical wavelength WL = CV * APD90, in cm."""
    return cv_cm_s * (apd90_ms / 1000.0)


def cardiac_gravitational_number(
    apd90_ms: float,
    dilation: float,
    l0_cm: float = 8.0,
    cv_cm_s: float = 58.0,
) -> float:
    """Dimensionless N_g = L0*dilation / WL.

    l0_cm: characteristic atrial path length at ground (human left-atrial
    circumference ~8-12 cm; default 8). It dilates with the fluid shift.
    """
    wl = wavelength_cm(apd90_ms, cv_cm_s)
    return (l0_cm * dilation) / wl


def interpolate_crossing(xs, ys, target: float = 1.0):
    """First x (linear-interpolated) where y crosses ``target``.

    ``xs`` ascending. Returns None if no sign change of (y - target).
    """
    for i in range(len(xs) - 1):
        d0 = ys[i] - target
        d1 = ys[i + 1] - target
        if d0 == 0.0:
            return xs[i]
        if d0 * d1 < 0.0:
            t = (target - ys[i]) / (ys[i + 1] - ys[i])
            return xs[i] + t * (xs[i + 1] - xs[i])
    return None
```

- [ ] **Step 4: Run test to verify it passes**

Run: `.venv/bin/pytest tests/test_gravity_law.py -q`
Expected: PASS (9 passed)

- [ ] **Step 5: Commit**

```bash
git add src/afib_microgravity/gravity_law.py tests/test_gravity_law.py
git commit -m "feat(gravity_law): wavelength, N_g, critical-gravity crossing"
```

---

### Task 3: The gravity sweep experiment (APD90 → N_g curve → g*)

**Files:**
- Create: `experiments/gravity_sweep.py`

This reuses the single-cell APD90 measurement pattern from `experiments/apd_wavelength.py` (CRN cell, `ap_biomarkers`) and the `gravity_to_remodeling` severity to set `microgravity_crn_params(severity=...)`.

- [ ] **Step 1: Write the experiment**

```python
# experiments/gravity_sweep.py
"""Gravity sweep: the Cardiac Gravitational Number N_g vs gravitational level.

For each gravity level g we map g -> remodelling severity/dilation, measure the
single-cell APD90 under microgravity_crn_params(severity), form the electrical
wavelength WL = CV*APD90 and the dimensionless N_g = L0*dilation/WL, then locate
the critical gravity g* where N_g crosses 1. Markers for Moon (0.16), Mars
(0.38) and interplanetary transit (0.0) are overlaid.

Outputs:
  * figures/gravity_law.png        -- N_g vs g, with g* and body markers
  * figures/results_crn.json       -- under ["gravity_law"]

Run:  python experiments/gravity_sweep.py [--full]
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import argparse
import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from afib_microgravity.crn import CRNCell, ap_biomarkers  # noqa: E402
from afib_microgravity.remodeling import microgravity_crn_params  # noqa: E402
from afib_microgravity.gravity_law import (  # noqa: E402
    gravity_to_remodeling,
    cardiac_gravitational_number,
    interpolate_crossing,
)

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
RESULTS = os.path.join(FIG, "results_crn.json")

DT = 0.02
CV_CM_S = 58.0          # planar CV (calibrated; g_Na unchanged across conditions)
L0_CM = 8.0             # ground atrial characteristic path length
BODIES = {"transit": 0.0, "Moon": 0.16, "Mars": 0.38, "Earth": 1.0}


def apd90_ms(params, bcl=1000.0, n_beats=12):
    cell = CRNCell(shape=(1, 1), params=params, use_numba=True)
    spb = int(bcl / DT)
    ss = int(2.0 / DT)
    tl, vl = [], []
    for beat in range(n_beats):
        for s in range(spb):
            I = -2300.0 if s < ss else 0.0
            dV = cell.reaction_step(DT, I_stim=I)
            cell.set_V(cell.V + DT * dV)
            if beat == n_beats - 1:
                vl.append(float(cell.V.ravel()[0]))
                tl.append(s * DT)
    return ap_biomarkers(np.array(tl), np.array(vl))["APD90"]


def merge_results(section, payload):
    os.makedirs(FIG, exist_ok=True)
    data = {}
    if os.path.exists(RESULTS):
        with open(RESULTS) as f:
            data = json.load(f)
    data[section] = payload
    with open(RESULTS, "w") as f:
        json.dump(data, f, indent=2)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--full", action="store_true")
    args = ap.parse_args()

    gs = ([0.0, 0.16, 0.38, 0.5, 0.8, 1.0] if args.full
          else [0.0, 0.5, 1.0])

    rows = []
    for g in gs:
        rem = gravity_to_remodeling(g)
        apd = apd90_ms(microgravity_crn_params(severity=rem.severity))
        ng = cardiac_gravitational_number(apd, dilation=rem.dilation,
                                          l0_cm=L0_CM, cv_cm_s=CV_CM_S)
        rows.append({"g": g, "severity": rem.severity, "dilation": rem.dilation,
                     "apd90_ms": apd, "wavelength_cm": CV_CM_S * apd / 1000.0,
                     "N_g": ng})
        print(f"  g={g:.2f}  sev={rem.severity:.2f}  APD90={apd:.1f}ms  N_g={ng:.2f}",
              flush=True)

    gstar = interpolate_crossing([r["g"] for r in rows],
                                 [r["N_g"] for r in rows], target=1.0)

    # figure
    fig, axis = plt.subplots(figsize=(6.5, 4.2))
    axis.plot([r["g"] for r in rows], [r["N_g"] for r in rows],
              "o-", color="#1b4079", lw=2, label=r"$\mathcal{N}_g$")
    axis.axhline(1.0, color="#888", ls="--", lw=1)
    axis.fill_between([r["g"] for r in rows], 1.0,
                      [max(r["N_g"], 1.0) for r in rows],
                      color="#d65a31", alpha=0.15, label=r"$\mathcal{N}_g>1$ (vulnerable)")
    for name, gval in BODIES.items():
        axis.axvline(gval, color="#444", ls=":", lw=0.8)
        axis.text(gval, axis.get_ylim()[1], name, rotation=90,
                  va="top", ha="right", fontsize=8, color="#444")
    if gstar is not None:
        axis.plot([gstar], [1.0], "*", color="#e84545", ms=15,
                  label=fr"$g^*$={gstar:.2f}")
    axis.set_xlabel("gravitational level  g  (Earth-g units)")
    axis.set_ylabel(r"Cardiac Gravitational Number  $\mathcal{N}_g$")
    axis.set_title("Gravitational scaling of atrial re-entry vulnerability")
    axis.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "gravity_law.png"), dpi=130)
    plt.close(fig)

    merge_results("gravity_law", {
        "mode": "full" if args.full else "smoke",
        "cv_cm_s": CV_CM_S, "l0_cm": L0_CM,
        "g_star": gstar, "rows": rows,
    })
    print(f"[gravity_sweep] g*={gstar}  wrote figures/gravity_law.png", flush=True)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the experiment in smoke mode**

Run: `.venv/bin/python experiments/gravity_sweep.py`
Expected: three `g=...` lines printed, a `g*=...` line, and `figures/gravity_law.png` written. (APD90 should fall as g decreases; N_g should rise.)

- [ ] **Step 3: Commit**

```bash
git add experiments/gravity_sweep.py figures/gravity_law.png figures/results_crn.json
git commit -m "feat(experiments): gravity sweep -> N_g curve and critical gravity g*"
```

---

### Task 4: Smoke test wiring the experiment into CI

**Files:**
- Create: `tests/test_gravity_sweep_smoke.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_gravity_sweep_smoke.py
import json
import os
import subprocess
import sys

ROOT = os.path.join(os.path.dirname(__file__), "..")


def test_gravity_sweep_smoke_runs_and_writes_results():
    proc = subprocess.run(
        [sys.executable, "experiments/gravity_sweep.py"],
        cwd=ROOT, capture_output=True, text=True, timeout=600,
    )
    assert proc.returncode == 0, proc.stderr
    with open(os.path.join(ROOT, "figures", "results_crn.json")) as f:
        data = json.load(f)
    assert "gravity_law" in data
    rows = data["gravity_law"]["rows"]
    gmap = {r["g"]: r["N_g"] for r in rows}
    # vulnerability rises as gravity falls: N_g(0) > N_g(1)
    assert gmap[0.0] > gmap[1.0]
```

- [ ] **Step 2: Run test to verify it fails (red) then passes after Task 3 exists**

Run: `.venv/bin/pytest tests/test_gravity_sweep_smoke.py -q`
Expected: PASS (Task 3 already created the experiment). If `gravity_sweep.py` were absent it would fail with returncode != 0.

- [ ] **Step 3: Commit**

```bash
git add tests/test_gravity_sweep_smoke.py
git commit -m "test: smoke test for gravity sweep experiment"
```

---

### Task 5: README section + full run + verify CI

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Run the full sweep to produce the publication figure**

Run: `.venv/bin/python experiments/gravity_sweep.py --full`
Expected: six `g=...` lines (0, 0.16, 0.38, 0.5, 0.8, 1.0), a `g*` value, updated `figures/gravity_law.png`.

- [ ] **Step 2: Add a README section** (insert after the existing headline block)

```markdown
## A gravitational scaling law (N_g)

Beyond the binary ground-vs-microgravity contrast, we treat gravity as a
*continuous* variable. A simple linear map sends the gravitational level `g`
(Earth-g units) to atrial remodelling via the cephalad fluid-shift drive
`(1 - g)`, calibrated so `g=1` recovers the validated ground baseline and
`g=0` recovers the microgravity operating point. From the measured APD90 we
form the electrical wavelength `WL = CV·APD90` and a **dimensionless Cardiac
Gravitational Number**

> `N_g = L0·dilation(g) / WL(g)`

Re-entry can be hosted when `N_g ≳ 1`; the **critical gravity `g*`** solves
`N_g(g*) = 1`. `experiments/gravity_sweep.py --full` sweeps
`g ∈ {0, 0.16, 0.38, 0.5, 0.8, 1.0}` and renders `figures/gravity_law.png`,
with markers for the Moon, Mars and interplanetary transit.

**This is a falsifiable hypothesis, not a proof.** The linear `g→remodelling`
map and the mechano-electric sensitivities are assumptions; `N_g` predicts a
rise in *substrate vulnerability*, not guaranteed atrial fibrillation —
consistent with Khine et al. 2018 (risk markers up, no sustained AF observed).
```

- [ ] **Step 3: Run the whole test suite**

Run: `.venv/bin/pytest -q`
Expected: all tests pass (existing + new gravity_law + smoke).

- [ ] **Step 4: Commit**

```bash
git add README.md figures/gravity_law.png figures/results_crn.json
git commit -m "docs: README section for the gravitational scaling law"
```

---

## Self-review

- **Spec coverage:** mechanism chain (§4 spec) → encoded in `gravity_to_remodeling` + N_g (Tasks 1-2); law derivation (§5) → Tasks 1-2; validation sweep + figure + g* with bodies (§6) → Tasks 3-5; integrity safeguards (§8) → honest docstrings/README, "proposed not proven", assumptions labelled. CI green requirement → Tasks 4-5.
- **Placeholder scan:** none — every step has concrete code/commands.
- **Type consistency:** `gravity_to_remodeling` returns `GravityRemodeling(.severity/.dilation/.density)`, used identically in Task 3; `cardiac_gravitational_number(apd90_ms, dilation, l0_cm, cv_cm_s)` and `interpolate_crossing(xs, ys, target)` signatures match between Task 2 and Task 3.
- **Honesty note:** confidence-interval / ensemble band on N_g is intentionally deferred — APD90 is deterministic given params, so the headline N_g curve has no stochasticity to bootstrap. If a vulnerability band is wanted, it would come from PS-density ensembles (separate follow-up), and that limitation is stated rather than faked.
