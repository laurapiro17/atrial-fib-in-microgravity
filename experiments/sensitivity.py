"""Sensitivity analysis: is the ground-vs-microgravity wavebreak difference robust
to the two main analysis choices -- the artifact-rejection threshold and the
conduction calibration (d_long)?

(1) Artifact-threshold robustness is cheap: we record the full PS(t) trace once
    per condition/seed and re-compute the wavebreak burden at several thresholds
    post-hoc (no re-simulation). If the ground<microgravity ordering holds across
    thresholds, the result is not an artifact of one cutoff.
(2) d_long robustness (optional, --full): re-run at a few conduction velocities.

Run:  python experiments/sensitivity.py [--full]
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
from afib_microgravity.metrics import count_phase_singularities  # noqa: E402
from afib_microgravity.model import MonodomainSheet  # noqa: E402
from afib_microgravity.remodeling import baseline_crn_params, microgravity_crn_params  # noqa: E402

from ensemble_with_ci import MIN_ANALYSIS_MS, S2_DELAY_MS, crn_phase, s1_planar, s2_quadrant  # noqa: E402

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
RESULTS = os.path.join(FIG, "results_crn.json")
DT = 0.02
D_TRANS_RATIO = 1.0 / 3.0
DX = 0.25
DELAY_MS = 30.0


def trace(cond, grid, seed, d_long, total_ms, sample_every):
    """Full PS(t) trace for one run (no thresholding here)."""
    shape = (grid, grid)
    if cond == "ground":
        params, coupling = baseline_crn_params(), None
    else:
        params = microgravity_crn_params(severity=1.0)
        coupling = correlated_fibrosis(shape, density=0.3, corr_len=4.0, seed=seed)
    cell = CRNCell(shape=shape, params=params, use_numba=True)
    diff = AnisotropicDiffusion(shape, d_long, d_long * D_TRANS_RATIO,
                                theta=np.zeros(shape), dx=DX, coupling=coupling)
    sheet = MonodomainSheet(cell, diff, dt=DT)
    n_steps = round(total_ms / DT); s2_step = round(S2_DELAY_MS / DT)
    delay_snaps = max(1, round(DELAY_MS / (sample_every * DT)))
    ring, ps, t = [], [], []
    s1_planar(cell)
    for i in range(n_steps):
        if i == s2_step:
            s2_quadrant(cell)
        sheet.step()
        if i % sample_every == 0:
            ring.append(cell.V.copy())
            if len(ring) > delay_snaps + 1:
                ring.pop(0)
            ps.append(count_phase_singularities(crn_phase(cell.V, ring[0])))
            t.append(i * DT)
    return np.asarray(ps, float), np.asarray(t, float)


def burden(ps, t, thr, total_ms):
    valid = (t >= MIN_ANALYSIS_MS) & (ps <= thr)
    if not valid.any():
        return 0.0
    w = ps[valid]
    dt_s = float(np.median(np.diff(t[valid]))) if valid.sum() > 1 else 0.0
    return float(w.sum() * dt_s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--grid", type=int, default=None)
    ap.add_argument("--seeds", type=int, default=None)
    args = ap.parse_args()
    if args.full:
        grid = args.grid or 160; n_seeds = args.seeds or 4; total_ms = 800.0; samp = 100
        d_longs = [0.10, 0.15, 0.20]
    else:
        grid = args.grid or 60; n_seeds = args.seeds or 2; total_ms = 250.0; samp = 50
        d_longs = [0.15]
    thresholds = [20, 30, 40, 60, 80]
    t0 = time.time()
    out = {"mode": "full" if args.full else "smoke", "grid": grid, "n_seeds": n_seeds,
           "thresholds": thresholds, "threshold_sweep": {}, "d_long_sweep": {}}

    # (1) record traces once at the nominal d_long=0.15, then re-threshold post-hoc
    print("[sensitivity] recording traces (d_long=0.15) ...")
    traces = {c: [trace(c, grid, s, 0.15, total_ms, samp) for s in range(n_seeds)]
              for c in ("ground", "microgravity")}
    for thr in thresholds:
        row = {}
        for c in ("ground", "microgravity"):
            bs = [burden(ps, t, thr, total_ms) for ps, t in traces[c]]
            row[c] = float(np.mean(bs))
        out["threshold_sweep"][str(thr)] = row
        print(f"  thr={thr:3d}  ground={row['ground']:.0f}  micro={row['microgravity']:.0f}  "
              f"micro>ground={row['microgravity'] > row['ground']}")

    # (2) d_long robustness
    for dl in d_longs:
        row = {}
        for c in ("ground", "microgravity"):
            bs = []
            for s in range(n_seeds):
                ps, t = (traces[c][s] if dl == 0.15 else trace(c, grid, s, dl, total_ms, samp))
                bs.append(burden(ps, t, 40, total_ms))
            row[c] = float(np.mean(bs))
        out["d_long_sweep"][str(dl)] = row
        print(f"  d_long={dl}  ground={row['ground']:.0f}  micro={row['microgravity']:.0f}")

    # figure: burden vs threshold
    fig, ax = plt.subplots(figsize=(6, 4))
    for c, col in (("ground", "#1b4079"), ("microgravity", "#c1121f")):
        ys = [out["threshold_sweep"][str(thr)][c] for thr in thresholds]
        ax.plot(thresholds, ys, "o-", color=col, label=c)
    ax.set_xlabel("artifact-rejection threshold (PS)")
    ax.set_ylabel("wavebreak burden (PS·ms)")
    ax.set_title("Robustness of the wavebreak difference to the artifact cutoff")
    ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(FIG, "sensitivity.png"), dpi=120); plt.close(fig)

    out["wall_seconds"] = round(time.time() - t0, 1)
    data = json.load(open(RESULTS)) if os.path.exists(RESULTS) else {}
    data["sensitivity"] = out
    json.dump(data, open(RESULTS, "w"), indent=2)
    print(f"[sensitivity] wall={out['wall_seconds']}s -> figures/sensitivity.png")


if __name__ == "__main__":
    main()
