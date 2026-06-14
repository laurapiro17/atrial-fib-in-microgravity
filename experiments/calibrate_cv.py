"""Calibrate d_long for physiological planar CV, and APD90/wavelength
for ground vs microgravity substrates. Numba-accelerated single cell + strand."""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import numpy as np
from afib_microgravity.crn import CRNCell, CRNParams, ap_biomarkers
from afib_microgravity.diffusion import AnisotropicDiffusion
from afib_microgravity.model import MonodomainSheet


def apd90_single(params, bcl=1000.0, n_beats=10, dt=0.02,
                 stim_amp=-2300.0, stim_dur=2.0):
    """Pace a single CRN cell to steady state; return last-beat APD90."""
    cell = CRNCell(shape=(1, 1), params=params, use_numba=True)
    spb = int(bcl / dt); ss = int(stim_dur / dt)
    tl, vl = [], []
    for beat in range(n_beats):
        for s in range(spb):
            I = stim_amp if s < ss else 0.0
            dV = cell.reaction_step(dt, I_stim=I)
            cell.set_V(cell.V + dt * dV)
            if beat == n_beats - 1:
                vl.append(float(cell.V.ravel()[0])); tl.append(s * dt)
    return ap_biomarkers(np.array(tl), np.array(vl))


def planar_cv(params, d_long, n=160, dx=0.25, dt=0.02, d_trans=None,
              bcl=1000.0, n_s1=2):
    """Planar CV (cm/s) on a 1xN strand paced from the left end.
    Times the upstroke (V crosses -40 mV) at two interior probes on the
    final (2nd) S1 beat."""
    if d_trans is None:
        d_trans = d_long
    shape = (1, n)
    cell = CRNCell(shape=shape, params=params, use_numba=True)
    diff = AnisotropicDiffusion(shape, d_long, d_trans, theta=np.zeros(shape), dx=dx)
    sheet = MonodomainSheet(cell, diff, dt=dt)
    left = np.zeros(shape, dtype=bool); left[0, 0:3] = True
    p1, p2 = n // 4, 3 * n // 4
    dist_cm = (p2 - p1) * dx * 0.1

    def beat(record=False):
        ts, v1s, v2s = [], [], []
        t = 0.0
        ss = int(2.0 / dt)
        for s in range(int(bcl / dt)):
            if s < ss:
                cell.stimulate(left, 20.0)
            sheet.step(); t += dt
            if record:
                ts.append(t); v1s.append(float(cell.V[0, p1]))
                v2s.append(float(cell.V[0, p2]))
        return ts, v1s, v2s

    for _ in range(n_s1 - 1):
        beat(record=False)
    ts, v1s, v2s = beat(record=True)

    def upstroke(ts, vs, thr=-40.0):
        for k in range(1, len(vs)):
            if vs[k - 1] < thr <= vs[k]:
                v0, v1 = vs[k - 1], vs[k]; t0, t1 = ts[k - 1], ts[k]
                return t0 + (thr - v0) * (t1 - t0) / (v1 - v0)
        return float("nan")

    t1 = upstroke(ts, v1s); t2 = upstroke(ts, v2s)
    if not (t2 > t1):
        return float("nan")
    return dist_cm / (t2 - t1) * 1000.0


if __name__ == "__main__":
    base = CRNParams()
    print("=== CV vs d_long (baseline params, isotropic strand) ===", flush=True)
    for dl in (0.06, 0.10, 0.15, 0.20, 0.25, 0.30):
        t0 = time.perf_counter()
        cv = planar_cv(base, d_long=dl)
        print(f"  d_long={dl:.3f}: CV={cv:.1f} cm/s  ({time.perf_counter()-t0:.1f}s)", flush=True)
