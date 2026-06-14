# CRN Ionic Substrate Upgrade — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the dimensionless Aliev–Panfilov substrate with a validated human-atrial ionic model (Courtemanche–Ramirez–Nattel 1998), on an anisotropic 2D monodomain sheet, with correlated fibrosis, restitution read-outs, and a statistically honest ensemble — while keeping AP available behind a common interface.

**Architecture:** A `CellModel` interface lets the `AtrialSheet` integrate either AP or CRN. CRN gating uses Rush–Larsen; the sheet uses operator splitting (reaction then explicit anisotropic finite-volume diffusion). Remodelling is expressed as CRN conductance scalings + geometry + correlated-fibrosis conductivity field + fibre anisotropy. Everything is in physical units (mV, ms, mm, cm/s, mM).

**Tech Stack:** Python ≥3.9, NumPy, SciPy, Matplotlib, pytest. Pure NumPy (vectorised over the grid). Isolated virtualenv only — never touch the global environment.

---

## Ground rules for the implementer

- **Commits:** author and committer are always `Laura Piñero Roig <laurapineroroig@gmail.com>`. **Never add a `Co-Authored-By: Claude` trailer** (project owner's standing preference). Use `git -c user.name=... -c user.email=...` if needed.
- **Environment:** `python3 -m venv .venv && source .venv/bin/activate && pip install -e ".[dev]"`. Do not install anything globally.
- **TDD:** every task is red → green → commit. Run the named test after each step.
- **Units:** mV, ms, mm, mM, pA/pF. Keep them consistent; the validation tests enforce physicality.
- **Backwards-compat:** the existing AP path and its tests must keep passing throughout. CI matrix (3.9/3.11/3.12) must stay green.

## File structure (decomposition)

| File | Responsibility |
|---|---|
| `src/afib_microgravity/cell_model.py` | `CellModel` protocol/ABC + `APCell` adapter wrapping current AP kinetics |
| `src/afib_microgravity/crn.py` | CRN single-cell ionic kinetics, vectorised, Rush–Larsen gates |
| `src/afib_microgravity/diffusion.py` | Anisotropic finite-volume `div(D·∇V)` from a fibre field; isotropic limit |
| `src/afib_microgravity/model.py` | `AtrialSheet` refactor: holds a `CellModel` + diffusion op; operator-split `step()` |
| `src/afib_microgravity/fibrosis.py` | Correlated (Gaussian-random-field) low-coupling conductivity field |
| `src/afib_microgravity/remodeling.py` | Spaceflight → CRN scalings + geometry + fibrosis + anisotropy; toggles |
| `src/afib_microgravity/restitution.py` | S1–S2 APD90 and CV restitution on a 1D strand |
| `src/afib_microgravity/metrics.py` | + PS-density time series, sustained-AF classifier, AP biomarkers |
| `experiments/validate_single_cell.py` | AP biomarkers + figure; asserts physiological range |
| `experiments/restitution_curves.py` | ground vs microgravity restitution |
| `experiments/ensemble_with_ci.py` | n≥25/condition, bootstrap CIs |
| `experiments/dose_response.py` | sweeps over dilation / fibrosis density / remodelling severity |
| `tests/test_cell_model.py`, `test_crn_singlecell.py`, `test_diffusion.py`, `test_fibrosis.py`, `test_restitution.py`, `test_remodeling_crn.py`, `test_metrics.py` | unit tests |

Phases are ordered so each produces working, tested software. The keystone (Task 2, CRN validation) gates everything downstream.

---

## Task 1: `CellModel` interface + AP adapter

**Files:**
- Create: `src/afib_microgravity/cell_model.py`
- Test: `tests/test_cell_model.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_cell_model.py
import numpy as np
from afib_microgravity.cell_model import APCell

def test_apcell_roundtrip_shapes_and_rest():
    cell = APCell(shape=(8, 8))
    assert cell.V.shape == (8, 8)
    # at rest, ionic "current" (here dimensionless du-reaction) is ~0 with V,v=0
    I = cell.reaction_step(dt=0.02)
    assert I.shape == (8, 8)
    assert np.allclose(cell.V, 0.0, atol=1e-9)  # no stimulus => stays at rest
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_cell_model.py -v`
Expected: FAIL (`ModuleNotFoundError: afib_microgravity.cell_model`).

- [ ] **Step 3: Write minimal implementation**

```python
# src/afib_microgravity/cell_model.py
"""Common interface so AtrialSheet can drive either the AP or CRN kinetics.

A CellModel owns its own state arrays (the membrane variable plus any gating /
concentration variables). reaction_step advances the *local* kinetics one dt and
returns the membrane time-derivative contribution the sheet will combine with
diffusion. The sheet writes the diffused membrane field back via set_V.
"""
from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

from .model import APParams  # existing dataclass


class CellModel(ABC):
    @property
    @abstractmethod
    def V(self) -> np.ndarray: ...

    @abstractmethod
    def set_V(self, V: np.ndarray) -> None: ...

    @abstractmethod
    def reaction_step(self, dt: float) -> np.ndarray:
        """Advance local kinetics by dt; return dV/dt_reaction (same shape)."""

    @abstractmethod
    def stimulate(self, mask: np.ndarray, value: float) -> None: ...


class APCell(CellModel):
    """Aliev–Panfilov kinetics behind the CellModel interface (dimensionless)."""

    def __init__(self, shape=(200, 200), params: APParams | None = None):
        self.p = params or APParams()
        self.u = np.zeros(shape, dtype=float)
        self.v = np.zeros(shape, dtype=float)

    @property
    def V(self) -> np.ndarray:
        return self.u

    def set_V(self, V: np.ndarray) -> None:
        self.u = V

    def reaction_step(self, dt: float) -> np.ndarray:
        u, v, p = self.u, self.v, self.p
        i_ion = p.k * u * (u - p.a) * (u - 1.0) + u * v
        eps = p.eps0 + p.mu1 * v / (u + p.mu2)
        dv = eps * (-v - p.k * u * (u - p.a - 1.0))
        self.v = v + dt * dv
        return -i_ion  # reaction contribution to du/dt (diffusion added by sheet)

    def stimulate(self, mask: np.ndarray, value: float = 1.0) -> None:
        self.u[mask] = value
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_cell_model.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/afib_microgravity/cell_model.py tests/test_cell_model.py
git -c user.name="Laura Piñero Roig" -c user.email="laurapineroroig@gmail.com" \
  commit -m "feat: CellModel interface with Aliev-Panfilov adapter"
```

---

## Task 2: CRN single-cell ionic model (KEYSTONE)

**Approach (important):** Do **not** reproduce CRN rate constants from memory. Transcribe the equations from an authoritative source — the original appendix (Courtemanche, Ramirez & Nattel, *Am J Physiol* 1998;275:H301) or the Physiome/CellML CRN listing — into a vectorised NumPy module. The validation test below (published biomarkers) is the hard gate: the task is **not done** until single-cell biomarkers land in range. Build incrementally: stimulus + INa + IK1 first (gives an upstroke and a resting potential), then add repolarising currents, checking APD as you go.

**State layout (21 vars), all 2D ndarrays of grid shape:**
`V`; fast/gates: `m, h, j` (INa), `oa, oi` (Ito), `ua, ui` (IKur), `xr` (IKr), `xs` (IKs), `d, f, fca` (ICaL), `u_rel, v_rel, w_rel` (SR release); concentrations: `Nai, Ki, Cai, Ca_up, Ca_rel`.

**Gating integration — Rush–Larsen** for every gate `g` with steady state `g_inf` and time constant `tau_g`:
`g_new = g_inf - (g_inf - g) * exp(-dt / tau_g)`.
Voltage and concentrations: forward Euler within the reaction step. `I_ion` (pA/pF) is the sum of sarcolemmal currents; `dV/dt_reaction = -(I_ion + I_stim) / Cm` with `Cm = 100 pF`, capacitance per area handled consistently.

**Files:**
- Create: `src/afib_microgravity/crn.py`
- Test: `tests/test_crn_singlecell.py`

- [ ] **Step 1: Write the failing validation test (the gate)**

```python
# tests/test_crn_singlecell.py
import numpy as np
from afib_microgravity.crn import CRNCell, ap_biomarkers

def _pace_single_cell(cell, bcl=1000.0, n_beats=8, dt=0.02, stim_amp=-2800.0,
                      stim_dur=2.0):
    """Pace a 1x1 CRN cell at basic cycle length bcl (ms); return last-beat V trace."""
    steps_per_beat = int(bcl / dt)
    stim_steps = int(stim_dur / dt)
    trace, times = [], []
    for beat in range(n_beats):
        for s in range(steps_per_beat):
            I_stim = stim_amp if s < stim_steps else 0.0
            cell.reaction_step(dt, I_stim=I_stim)
            if beat == n_beats - 1:
                trace.append(float(cell.V.ravel()[0]))
                times.append(s * dt)
    return np.array(times), np.array(trace)

def test_crn_single_cell_biomarkers_in_physiological_range():
    cell = CRNCell(shape=(1, 1))
    t, V = _pace_single_cell(cell)
    bm = ap_biomarkers(t, V)
    assert -86.0 <= bm["V_rest"] <= -75.0,  bm     # resting ~ -81 mV
    assert  10.0 <= bm["V_peak"] <= 35.0,   bm     # overshoot ~ +20 mV
    assert 250.0 <= bm["APD90"]  <= 360.0,  bm     # human atrial APD90 ~ 300 ms
    assert bm["dVdt_max"] >= 100.0,          bm    # fast upstroke (V/s)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_crn_singlecell.py -v`
Expected: FAIL (`ModuleNotFoundError`).

- [ ] **Step 3: Implement CRN incrementally until the gate passes**

Create `src/afib_microgravity/crn.py` with:
- `@dataclass CRNParams` holding maximal conductances/scalings, each defaulting to the published CRN value: `g_Na, g_K1, g_to, g_Kur, g_Kr, g_Ks, g_CaL, I_NaK_max, I_NaCa_max, g_bCa, g_bNa, I_pCa_max` plus SR (`I_up_max`, `K_rel`, …) and fixed constants (`Cm=100 pF`, `R, T, F`, external `Nao=140, Ko=5.4, Cao=1.8 mM`).
- `class CRNCell(CellModel)` (import the ABC from `cell_model` — note: to avoid a cycle, `cell_model` imports only `APParams` from `model`; `crn` imports `CellModel` from `cell_model`). Implements `V`, `set_V`, `stimulate`, and `reaction_step(self, dt, I_stim=0.0)`.
- Inside `reaction_step`: compute each current from `self.V` and gates (transcribed from the CRN appendix), update gates by Rush–Larsen, update `Nai, Ki, Cai, Ca_up, Ca_rel` by forward Euler, set `dVdt = -(I_ion + I_stim) / 1.0` (currents already in pA/pF; Cm normalised), `self.V = self.V + dt*dVdt`, and `return dVdt`.
- `ap_biomarkers(t, V)` → dict with `V_rest=min(V)`, `V_peak=max(V)`, `dVdt_max=max(diff(V)/diff(t))`, and `APD90` computed at 90% repolarisation from the upstroke (V crossing `V_rest + 0.1*(V_peak - V_rest)` downward).

Build order to reach the gate: (a) stimulus + INa + IK1 → upstroke + rest; (b) add Ito, IKur, IKr, IKs, ICaL → repolarisation/APD; (c) add INaK, INaCa, background, IpCa + Ca handling → stable concentrations across beats. Re-run the test after each addition; tune nothing except fixing transcription errors — defaults are the published values.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_crn_singlecell.py -v`
Expected: PASS (biomarkers in range).

- [ ] **Step 5: Add a second guard test (rate dependence) and commit**

```python
# append to tests/test_crn_singlecell.py
def test_crn_apd_shortens_at_faster_rate():
    slow = CRNCell(shape=(1, 1)); fast = CRNCell(shape=(1, 1))
    _, Vs = _pace_single_cell(slow, bcl=1000.0)
    _, Vf = _pace_single_cell(fast, bcl=400.0)
    assert ap_biomarkers(_pace_single_cell(fast, bcl=400.0)[0], Vf)["APD90"] < \
           ap_biomarkers(_pace_single_cell(slow, bcl=1000.0)[0], Vs)["APD90"]
```

Run: `pytest tests/test_crn_singlecell.py -v` → PASS, then:

```bash
git add src/afib_microgravity/crn.py tests/test_crn_singlecell.py
git -c user.name="Laura Piñero Roig" -c user.email="laurapineroroig@gmail.com" \
  commit -m "feat: validated CRN human-atrial single-cell ionic model (Rush-Larsen)"
```

---

## Task 3: Anisotropic finite-volume diffusion

**Files:**
- Create: `src/afib_microgravity/diffusion.py`
- Test: `tests/test_diffusion.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_diffusion.py
import numpy as np
from afib_microgravity.diffusion import AnisotropicDiffusion

def test_isotropic_limit_matches_5point_laplacian():
    rng = np.random.default_rng(0)
    u = rng.standard_normal((16, 16))
    op = AnisotropicDiffusion(shape=(16,16), d_long=1.0, d_trans=1.0,
                              theta=np.zeros((16,16)), dx=1.0)
    lap = (np.roll(u,1,0)+np.roll(u,-1,0)+np.roll(u,1,1)+np.roll(u,-1,1)-4*u)
    # interior agreement (boundaries differ due to no-flux)
    assert np.allclose(op.apply(u)[1:-1,1:-1], lap[1:-1,1:-1], atol=1e-9)

def test_no_flux_conserves_total_for_uniform_field():
    u = np.ones((10, 10)) * 0.7
    op = AnisotropicDiffusion(shape=(10,10), d_long=2.0, d_trans=0.5,
                              theta=np.full((10,10), 0.3), dx=0.25)
    assert np.allclose(op.apply(u), 0.0, atol=1e-9)  # flat field => zero diffusion
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/test_diffusion.py -v` → FAIL (`ModuleNotFoundError`).

- [ ] **Step 3: Implement**

```python
# src/afib_microgravity/diffusion.py
"""Anisotropic monodomain diffusion: div(D grad V) with a per-cell conductivity
tensor built from a fibre-angle field. No-flux (Neumann) boundaries.

D = R(theta) diag(d_long, d_trans) R(theta)^T, so
  Dxx = d_long cos^2 + d_trans sin^2
  Dyy = d_long sin^2 + d_trans cos^2
  Dxy = (d_long - d_trans) sin cos
The cross (Dxy) term uses centred differences; axial terms use a harmonic-mean
finite-volume flux (robust across sharp conductivity contrasts, e.g. fibrosis).
"""
from __future__ import annotations

import numpy as np


class AnisotropicDiffusion:
    def __init__(self, shape, d_long, d_trans, theta, dx=0.25, coupling=None):
        self.ny, self.nx = shape
        self.dx = float(dx)
        c, s = np.cos(theta), np.sin(theta)
        scale = np.ones(shape) if coupling is None else np.asarray(coupling, float)
        dl = d_long * scale
        dt_ = d_trans * scale
        self.Dxx = dl * c * c + dt_ * s * s
        self.Dyy = dl * s * s + dt_ * c * c
        self.Dxy = (dl - dt_) * s * c

    @staticmethod
    def _hmean(a, b):
        s = a + b
        out = np.zeros_like(a)
        nz = s > 0
        out[nz] = 2.0 * a[nz] * b[nz] / s[nz]
        return out

    def apply(self, V):
        inv = 1.0 / (self.dx * self.dx)
        # axial x
        Dr = self._hmean(self.Dxx, np.roll(self.Dxx, -1, 1))
        fx_r = Dr * (np.roll(V, -1, 1) - V); fx_r[:, -1] = 0.0
        Dl = self._hmean(self.Dxx, np.roll(self.Dxx, 1, 1))
        fx_l = Dl * (V - np.roll(V, 1, 1)); fx_l[:, 0] = 0.0
        # axial y
        Dd = self._hmean(self.Dyy, np.roll(self.Dyy, -1, 0))
        fy_d = Dd * (np.roll(V, -1, 0) - V); fy_d[-1, :] = 0.0
        Du = self._hmean(self.Dyy, np.roll(self.Dyy, 1, 0))
        fy_u = Du * (V - np.roll(V, 1, 0)); fy_u[0, :] = 0.0
        axial = ((fx_r - fx_l) + (fy_d - fy_u)) * inv
        # cross term Dxy * d2V/dxdy (centred, zeroed at borders via no-flux proxy)
        dVdy = (np.roll(V, -1, 0) - np.roll(V, 1, 0)) * 0.5
        cross = (np.roll(self.Dxy * dVdy, -1, 1) -
                 np.roll(self.Dxy * dVdy, 1, 1)) * 0.5 * inv
        cross[0, :] = cross[-1, :] = cross[:, 0] = cross[:, -1] = 0.0
        return axial + cross

    @property
    def Dmax(self):
        return float(max(self.Dxx.max(), self.Dyy.max()))
```

- [ ] **Step 4: Run to verify it passes**

Run: `pytest tests/test_diffusion.py -v` → PASS.

- [ ] **Step 5: Commit**

```bash
git add src/afib_microgravity/diffusion.py tests/test_diffusion.py
git -c user.name="Laura Piñero Roig" -c user.email="laurapineroroig@gmail.com" \
  commit -m "feat: anisotropic finite-volume diffusion with fibre tensor"
```

---

## Task 4: Refactor `AtrialSheet` to operator splitting over a `CellModel`

**Files:**
- Modify: `src/afib_microgravity/model.py` (keep `APParams`; add `AtrialSheet2`-style path without breaking the existing class/tests)
- Test: `tests/test_model.py` (extend; existing AP tests must still pass)

Design note: keep the existing `AtrialSheet` (AP, isotropic) untouched so legacy tests pass. Add a new `MonodomainSheet` that composes any `CellModel` + `AnisotropicDiffusion`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_model.py (append)
import numpy as np
from afib_microgravity.cell_model import APCell
from afib_microgravity.diffusion import AnisotropicDiffusion
from afib_microgravity.model import MonodomainSheet

def test_monodomain_propagates_a_planar_wave():
    shape = (40, 40)
    cell = APCell(shape=shape)
    diff = AnisotropicDiffusion(shape, d_long=1.0, d_trans=1.0,
                                theta=np.zeros(shape), dx=1.0)
    sheet = MonodomainSheet(cell, diff, dt=0.01)
    cell.stimulate(np.s_[:, :3] and (np.arange(shape[1])[None, :] < 3) &
                   np.ones(shape, bool), value=1.0)
    front0 = (cell.V > 0.5).sum()
    sheet.run(400)
    assert (cell.V > 0.5).sum() > front0  # wave advanced into the sheet
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest tests/test_model.py::test_monodomain_propagates_a_planar_wave -v` → FAIL.

- [ ] **Step 3: Implement `MonodomainSheet`**

```python
# src/afib_microgravity/model.py (append; do not modify existing AtrialSheet)
class MonodomainSheet:
    """Operator-split monodomain: reaction (CellModel) then explicit diffusion."""

    def __init__(self, cell, diffusion, dt=0.01):
        self.cell = cell
        self.diff = diffusion
        self.dt = float(dt)
        self.t = 0.0
        self._check_stability()

    def _check_stability(self):
        bound = self.diff.dx ** 2 / (4.0 * max(self.diff.Dmax, 1e-12))
        if self.dt > bound:
            raise ValueError(f"dt={self.dt} exceeds diffusion bound {bound:.4g}")

    def step(self, I_stim=0.0):
        dV_react = self.cell.reaction_step(self.dt) if I_stim == 0.0 \
            else self.cell.reaction_step(self.dt)  # AP ignores I_stim
        V = self.cell.V + self.dt * dV_react
        V = V + self.dt * self.diff.apply(V)
        self.cell.set_V(V)
        self.t += self.dt

    def run(self, n_steps, I_stim=0.0):
        for _ in range(n_steps):
            self.step(I_stim)
```

(Note: for CRN, reaction already advances V internally; the sheet then applies the diffusion increment to `cell.V`. Implement `reaction_step` to return dV/dt and NOT mutate V for both models — adjust APCell/CRNCell so V mutation happens only in the sheet. Update Task 1/2 code accordingly: `reaction_step` returns dV/dt_reaction and updates only gates/concentrations.)

- [ ] **Step 4: Reconcile reaction/diffusion contract**

Make both `APCell.reaction_step` and `CRNCell.reaction_step` **return dV/dt and not write V** (gates/concentrations updated in place). The sheet performs `V += dt*(dV_react + diffusion)`. Re-run Task 1, 2 tests and this test.

Run: `pytest tests/test_cell_model.py tests/test_crn_singlecell.py tests/test_model.py -v`
Expected: PASS (adjust the single-cell pacing helper to add `dt*dVdt` to `cell.V` itself, since there is no sheet in that test).

- [ ] **Step 5: Commit**

```bash
git add src/afib_microgravity/model.py src/afib_microgravity/cell_model.py \
        src/afib_microgravity/crn.py tests/
git -c user.name="Laura Piñero Roig" -c user.email="laurapineroroig@gmail.com" \
  commit -m "feat: operator-split MonodomainSheet over CellModel interface"
```

---

## Task 5: Correlated fibrosis conductivity field

**Files:**
- Create: `src/afib_microgravity/fibrosis.py`
- Test: `tests/test_fibrosis.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_fibrosis.py
import numpy as np
from afib_microgravity.fibrosis import correlated_fibrosis

def test_density_is_approximately_hit():
    f = correlated_fibrosis((128,128), density=0.3, corr_len=4.0,
                            d_healthy=1.0, d_scar=0.05, seed=1)
    frac_scar = np.mean(f < 0.5)
    assert abs(frac_scar - 0.3) < 0.05

def test_zero_density_is_homogeneous():
    f = correlated_fibrosis((32,32), density=0.0, corr_len=4.0, seed=1)
    assert np.allclose(f, 1.0)

def test_larger_corr_len_makes_larger_patches():
    small = correlated_fibrosis((128,128), 0.3, corr_len=1.0, seed=2)
    large = correlated_fibrosis((128,128), 0.3, corr_len=8.0, seed=2)
    # autocorrelation at lag 3 is higher for the larger correlation length
    def autocorr(a, lag):
        b = (a < 0.5).astype(float); b -= b.mean()
        return np.mean(b * np.roll(b, lag, axis=0))
    assert autocorr(large, 3) > autocorr(small, 3)
```

- [ ] **Step 2: Run to verify it fails** → `pytest tests/test_fibrosis.py -v` FAIL.

- [ ] **Step 3: Implement**

```python
# src/afib_microgravity/fibrosis.py
"""Spatially correlated fibrosis as a low-coupling conductivity field.

White noise smoothed by a Gaussian kernel (correlation length corr_len, in cells)
is thresholded at the quantile that yields the requested area fraction. Smoothing
gives realistic connected patches instead of isolated squares; thresholding pins
the density exactly. Returns a coupling-scale field (d_scar inside, d_healthy out).
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter


def correlated_fibrosis(shape, density=0.0, corr_len=4.0, d_healthy=1.0,
                        d_scar=0.05, seed=0):
    field = np.full(shape, d_healthy, dtype=float)
    if density <= 0:
        return field
    rng = np.random.default_rng(seed)
    noise = rng.standard_normal(shape)
    smooth = gaussian_filter(noise, sigma=corr_len, mode="reflect")
    thresh = np.quantile(smooth, density)        # lowest `density` fraction -> scar
    field[smooth <= thresh] = d_scar
    return field
```

- [ ] **Step 4: Run to verify it passes** → PASS.

- [ ] **Step 5: Commit**

```bash
git add src/afib_microgravity/fibrosis.py tests/test_fibrosis.py
git -c user.name="Laura Piñero Roig" -c user.email="laurapineroroig@gmail.com" \
  commit -m "feat: correlated Gaussian-random-field fibrosis"
```

---

## Task 6: Restitution (S1–S2) — APD90 and CV

**Files:**
- Create: `src/afib_microgravity/restitution.py`
- Test: `tests/test_restitution.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_restitution.py
from afib_microgravity.restitution import apd_restitution

def test_apd_restitution_is_monotonic_and_shortens_with_short_di():
    di, apd = apd_restitution(diastolic_intervals=(300., 200., 120., 60.),
                              dt=0.05)
    assert all(a > 0 for a in apd)
    assert apd[0] > apd[-1]      # shorter DI -> shorter APD
    assert all(apd[i] >= apd[i+1] for i in range(len(apd)-1))
```

- [ ] **Step 2: Run to verify it fails** → FAIL.

- [ ] **Step 3: Implement** (single CRN cell; S1 train then one S2 at each DI; measure APD90 of the S2 beat using `ap_biomarkers`). CV restitution function `cv_restitution(...)` paces a 1D strand (a `(1, N)` `MonodomainSheet` with CRN cells) and times the wavefront between two probes, returning cm/s for each DI. Keep durations short (a few S1 beats) so the test runs in seconds.

- [ ] **Step 4: Run to verify it passes** → PASS.

- [ ] **Step 5: Commit**

```bash
git add src/afib_microgravity/restitution.py tests/test_restitution.py
git -c user.name="Laura Piñero Roig" -c user.email="laurapineroroig@gmail.com" \
  commit -m "feat: APD and CV restitution (S1-S2)"
```

---

## Task 7: Remodelling mapping in CRN terms

**Files:**
- Modify: `src/afib_microgravity/remodeling.py`
- Test: `tests/test_remodeling_crn.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_remodeling_crn.py
from afib_microgravity.remodeling import microgravity_crn_params, baseline_crn_params

def test_af_type_remodelling_scales_the_intended_currents_only():
    base = baseline_crn_params()
    mg = microgravity_crn_params(severity=1.0)
    assert mg.g_CaL < base.g_CaL      # ICaL down
    assert mg.g_to  < base.g_to       # Ito down
    assert mg.g_Kur < base.g_Kur      # IKur down
    assert mg.g_Na == base.g_Na       # INa untouched
    assert mg.g_Kr == base.g_Kr       # IKr untouched

def test_severity_zero_recovers_baseline():
    base = baseline_crn_params(); mg = microgravity_crn_params(severity=0.0)
    assert (mg.g_CaL, mg.g_to, mg.g_Kur) == (base.g_CaL, base.g_to, base.g_Kur)
```

- [ ] **Step 2: Run to verify it fails** → FAIL.

- [ ] **Step 3: Implement** `baseline_crn_params()` (returns default `CRNParams`) and `microgravity_crn_params(severity)` applying literature-anchored AF-type scalings, each toggleable: `g_CaL *= 1 - 0.5*severity`, `g_to *= 1 - 0.5*severity`, `g_Kur *= 1 - 0.5*severity`. Add `make_condition_crn(name, base_shape, seed)` returning `(cell, diffusion)` assembling geometry (dilation for microgravity), `correlated_fibrosis`, fibre `theta`, and the CRN params. Keep the legacy AP `make_condition` intact.

- [ ] **Step 4: Run to verify it passes** → PASS.

- [ ] **Step 5: Commit**

```bash
git add src/afib_microgravity/remodeling.py tests/test_remodeling_crn.py
git -c user.name="Laura Piñero Roig" -c user.email="laurapineroroig@gmail.com" \
  commit -m "feat: CRN-based microgravity remodelling (AF-type ionic scalings)"
```

---

## Task 8: Extended metrics — PS density series, sustained-AF classifier

**Files:**
- Modify: `src/afib_microgravity/metrics.py`
- Test: `tests/test_metrics.py` (append; keep existing PS tests)

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_metrics.py (append)
import numpy as np
from afib_microgravity.metrics import is_sustained_af, ps_density

def test_sustained_af_true_when_ps_persists():
    series = [5, 6, 4, 7, 5, 6]      # many rotors throughout
    assert is_sustained_af(series, window_frac=0.5, threshold=2) is True

def test_sustained_af_false_when_quiescent():
    series = [3, 2, 1, 0, 0, 0]      # dies out
    assert is_sustained_af(series, window_frac=0.5, threshold=2) is False

def test_ps_density_normalises_by_area():
    assert ps_density(8, area_cells=40000) == 8 / 40000 * 1e4
```

- [ ] **Step 2: Run to verify it fails** → FAIL.

- [ ] **Step 3: Implement** `ps_density(count, area_cells)` and `is_sustained_af(series, window_frac, threshold)` (mean PS over the last `window_frac` of the series ≥ threshold). Reuse existing `count_phase_singularities`/`phase_field` for V-based phase (use a CRN-appropriate phase from `(V, dVdt)` or a fixed-V0 threshold).

- [ ] **Step 4: Run to verify it passes** → PASS.

- [ ] **Step 5: Commit**

```bash
git add src/afib_microgravity/metrics.py tests/test_metrics.py
git -c user.name="Laura Piñero Roig" -c user.email="laurapineroroig@gmail.com" \
  commit -m "feat: PS density and sustained-AF classifier"
```

---

## Task 9: Experiments + README + figures

**Files:**
- Create: `experiments/validate_single_cell.py`, `experiments/restitution_curves.py`,
  `experiments/ensemble_with_ci.py`, `experiments/dose_response.py`
- Modify: `README.md`; write outputs to `figures/`

- [ ] **Step 1: `validate_single_cell.py`** — pace one CRN cell, write `figures/crn_ap.png` and biomarkers to `figures/results.json["single_cell"]`; assert physiological range (mirrors the keystone test) so running it is self-checking.

- [ ] **Step 2: `restitution_curves.py`** — APD and CV restitution for ground vs microgravity; save `figures/restitution.png`.

- [ ] **Step 3: `ensemble_with_ci.py`** — for each condition, n≥25 fibrosis seeds; for each run a `MonodomainSheet` of CRN cells with broken-wavefront induction; record PS-density (2nd-half mean) and sustained-AF flag; compute **bootstrap 95% CIs**; write per-seed + summary to `figures/results.json`. Checkpoint after each seed so an overnight run is resumable.

```python
# bootstrap helper (put in metrics.py, test it in Task 8 if preferred)
def bootstrap_ci(values, n_boot=2000, alpha=0.05, seed=0):
    import numpy as np
    rng = np.random.default_rng(seed)
    vals = np.asarray(values, float)
    boots = [rng.choice(vals, size=vals.size, replace=True).mean()
             for _ in range(n_boot)]
    lo, hi = np.quantile(boots, [alpha/2, 1-alpha/2])
    return float(vals.mean()), float(lo), float(hi)
```

- [ ] **Step 4: `dose_response.py`** — sweep dilation factor, fibrosis density, remodelling severity independently; plot PS density vs each; save `figures/dose_response.png`.

- [ ] **Step 5: Regenerate animations + GIFs** (reuse the existing ffmpeg 2-pass-palette pipeline: `fps=12, scale=300, palettegen max_colors=64`), update the README comparison table; rewrite README results with **physical numbers + CIs** and an updated honesty section (now ionic substrate, still 2D monodomain). Keep CI green.

- [ ] **Step 6: Commit and push**

```bash
git add experiments/ README.md figures/
git -c user.name="Laura Piñero Roig" -c user.email="laurapineroroig@gmail.com" \
  commit -m "feat: CRN experiments (validation, restitution, ensemble CIs, dose-response) + README"
git push origin main
```

---

## Self-review (completed by plan author)

- **Spec coverage:** CRN single-cell (Task 2) ✓; cell-model interface w/ AP kept (Task 1) ✓; anisotropic diffusion (Task 3) ✓; operator-split sheet (Task 4) ✓; correlated fibrosis (Task 5) ✓; restitution APD+CV (Task 6) ✓; CRN remodelling toggles (Task 7) ✓; PS density + sustained-AF + bootstrap CI (Tasks 8–9) ✓; experiments + dose-response + README honesty (Task 9) ✓; physical units & validation biomarkers ✓.
- **Reaction/diffusion contract:** unified in Task 4 — `reaction_step` returns dV/dt and never mutates V; the sheet owns the V update. Single-cell pacing helper updates V itself (no sheet). Consistent across APCell/CRNCell.
- **Import cycle guard:** `model` defines `APParams`; `cell_model` imports `APParams` from `model`; `crn` imports `CellModel` from `cell_model`; `model.MonodomainSheet` references cells/diffusion only via duck-typing → no cycle.
- **Placeholder scan:** CRN rate constants are intentionally delegated to authoritative-source transcription gated by the biomarker test — this is a deliberate accuracy decision, not a vague placeholder; the acceptance criteria (exact biomarker ranges) are concrete.
- **Naming consistency:** `CRNParams.g_CaL/g_to/g_Kur` used identically in Tasks 2 and 7; `ap_biomarkers`, `correlated_fibrosis`, `AnisotropicDiffusion`, `MonodomainSheet`, `is_sustained_af`, `ps_density`, `bootstrap_ci` used consistently.
