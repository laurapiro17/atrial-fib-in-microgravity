"""Artifact control (audit C2): does the fibrotic substrate produce phase
singularities WITHOUT any re-entry?

We drive the microgravity fibrotic sheet with a SINGLE planar S1 wave (no S2, no
cross-field induction). A clean rotor detector must report ~0 phase singularities
once the wave has swept across and the tissue has repolarised, because there is
no rotation -- only a propagating front crossing low-coupling patches. If instead
PS persist (and co-locate with patch boundaries), the headline wavebreak count is
a detection artifact rather than a biological effect.

We report: peak PS during propagation, PS after the wave has cleared (the tail),
and the fraction of late-time PS that sit on a fibrosis boundary.

Run:  python experiments/artifact_control.py [--full]
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import argparse
import json
import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from afib_microgravity.crn import CRNCell  # noqa: E402
from afib_microgravity.diffusion import AnisotropicDiffusion  # noqa: E402
from afib_microgravity.fibrosis import correlated_fibrosis  # noqa: E402
from afib_microgravity.metrics import count_phase_singularities, phase_field  # noqa: E402
from afib_microgravity.model import MonodomainSheet  # noqa: E402
from afib_microgravity.remodeling import microgravity_crn_params  # noqa: E402

from ensemble_with_ci import (  # noqa: E402
    ARTIFACT_PS_THRESHOLD, MIN_ANALYSIS_MS, crn_phase, s1_planar,
)

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
RESULTS = os.path.join(FIG, "results_crn.json")
DT = 0.02
D_LONG, D_TRANS, DX = 0.15, 0.05, 0.25
DELAY_MS = 30.0


def boundary_mask(coupling):
    """True at cells adjacent to a fibrosis (low-coupling) boundary."""
    scar = coupling < 0.5
    nb = np.zeros_like(scar)
    nb[:-1, :] |= scar[1:, :]; nb[1:, :] |= scar[:-1, :]
    nb[:, :-1] |= scar[:, 1:]; nb[:, 1:] |= scar[:, :-1]
    return nb & ~scar | scar  # patch cells and their immediate neighbours


def run(grid, density, seed, total_ms, sample_every):
    shape = (grid, grid)
    coupling = correlated_fibrosis(shape, density=density, corr_len=4.0, seed=seed)
    cell = CRNCell(shape=shape, params=microgravity_crn_params(severity=1.0),
                   use_numba=True)
    diff = AnisotropicDiffusion(shape, D_LONG, D_TRANS, theta=np.zeros(shape),
                                dx=DX, coupling=coupling)
    sheet = MonodomainSheet(cell, diff, dt=DT)
    edge = boundary_mask(coupling)

    n_steps = round(total_ms / DT)
    delay_snaps = max(1, round(DELAY_MS / (sample_every * DT)))
    ring, ps_series, times = [], [], []
    late_on_edge, late_total = 0, 0

    s1_planar(cell)              # SINGLE planar S1 only -- NO S2, NO induction
    for i in range(n_steps):
        sheet.step()
        if i % sample_every == 0:
            ring.append(cell.V.copy())
            if len(ring) > delay_snaps + 1:
                ring.pop(0)
            ph = crn_phase(cell.V, ring[0])
            ps = count_phase_singularities(ph)
            ps_series.append(ps); times.append(i * DT)

    ps = np.asarray(ps_series, float); t = np.asarray(times, float)
    valid = (t >= MIN_ANALYSIS_MS) & (ps <= ARTIFACT_PS_THRESHOLD)
    # "tail" = last third of the (post-S1, artifact-rejected) trace: should be ~0
    vt = t[valid]; vps = ps[valid]
    tail = vps[vt >= (0.66 * total_ms)] if vt.size else np.array([])
    return {
        "peak_ps_valid": int(vps.max()) if vps.size else 0,
        "tail_mean_ps": float(tail.mean()) if tail.size else 0.0,
        "tail_max_ps": int(tail.max()) if tail.size else 0,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--grid", type=int, default=None)
    ap.add_argument("--seeds", type=int, default=None)
    args = ap.parse_args()
    if args.full:
        grid = args.grid or 160; n_seeds = args.seeds or 5; total_ms = 800.0; samp = 100
    else:
        grid = args.grid or 60; n_seeds = args.seeds or 2; total_ms = 300.0; samp = 50

    print(f"[artifact-control] microgravity fibrosis, SINGLE S1 (no induction); "
          f"grid={grid} seeds={n_seeds}")
    rows = [run(grid, 0.3, s, total_ms, samp) for s in range(n_seeds)]
    tail_means = [r["tail_mean_ps"] for r in rows]
    verdict = ("CLEAN (no spurious PS once wave clears -> headline effect is real "
               "re-entry, not a detector artifact)"
               if max(r["tail_max_ps"] for r in rows) <= 1
               else "WARNING: persistent PS without induction -> possible edge artifact")
    out = {"mode": "full" if args.full else "smoke", "grid": grid,
           "n_seeds": n_seeds, "per_seed": rows,
           "tail_mean_ps_mean": float(np.mean(tail_means)), "verdict": verdict}
    for s, r in enumerate(rows):
        print(f"  seed={s}  peak_valid={r['peak_ps_valid']}  "
              f"tail_mean={r['tail_mean_ps']:.2f}  tail_max={r['tail_max_ps']}")
    print(f"[artifact-control] VERDICT: {verdict}")

    data = json.load(open(RESULTS)) if os.path.exists(RESULTS) else {}
    data["artifact_control"] = out
    json.dump(data, open(RESULTS, "w"), indent=2)


if __name__ == "__main__":
    main()
