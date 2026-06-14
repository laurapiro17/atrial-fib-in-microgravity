"""Critical fibrosis density for wavebreak (percolation-style).

On the microgravity (AF-remodelled) kinetics, we sweep the interstitial-fibrosis
area fraction and measure the wavebreak burden at each density. This locates the
critical density at which the substrate starts to fragment wavefronts -- a
percolation-like threshold that mirrors the clinical observation that the
*pattern and amount* of atrial fibrosis gates AF, not the remodelling alone.

Run:  python experiments/fibrosis_threshold.py [--full]
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
from afib_microgravity.metrics import bootstrap_ci, count_phase_singularities  # noqa: E402
from afib_microgravity.model import MonodomainSheet  # noqa: E402
from afib_microgravity.remodeling import microgravity_crn_params  # noqa: E402

from ensemble_with_ci import (  # noqa: E402
    S2_DELAY_MS, crn_phase, s1_planar, s2_quadrant, window_metrics,
)

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
RESULTS = os.path.join(FIG, "results_crn.json")
DT = 0.02
D_LONG, D_TRANS, DX = 0.15, 0.05, 0.25
CORR = 4.0
DELAY_MS = 30.0


def burden_at_density(grid, density, seed, total_ms, sample_every):
    shape = (grid, grid)
    params = microgravity_crn_params(severity=1.0)
    coupling = (None if density <= 0 else
                correlated_fibrosis(shape, density=density, corr_len=CORR, seed=seed))
    cell = CRNCell(shape=shape, params=params, use_numba=True)
    diff = AnisotropicDiffusion(shape, D_LONG, D_TRANS, theta=np.zeros(shape),
                                dx=DX, coupling=coupling)
    sheet = MonodomainSheet(cell, diff, dt=DT)
    n_steps = round(total_ms / DT); s2_step = round(S2_DELAY_MS / DT)
    delay_snaps = max(1, round(DELAY_MS / (sample_every * DT)))
    ring, ps_series, times = [], [], []
    s1_planar(cell)
    for i in range(n_steps):
        if i == s2_step:
            s2_quadrant(cell)
        sheet.step()
        if i % sample_every == 0:
            ring.append(cell.V.copy())
            if len(ring) > delay_snaps + 1:
                ring.pop(0)
            ps_series.append(count_phase_singularities(crn_phase(cell.V, ring[0])))
            times.append(i * DT)
    ps = np.asarray(ps_series, float); t = np.asarray(times, float)
    return window_metrics(ps, t, S2_DELAY_MS, total_ms, grid * grid)["break_burden_ps_ms"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--grid", type=int, default=None)
    ap.add_argument("--seeds", type=int, default=None)
    args = ap.parse_args()
    if args.full:
        grid = args.grid or 160; n_seeds = args.seeds or 5
        densities = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]; total_ms = 800.0; sample_every = 100
    else:
        grid = args.grid or 60; n_seeds = args.seeds or 2
        densities = [0.0, 0.3, 0.5]; total_ms = 250.0; sample_every = 50

    print(f"[fibrosis-threshold] grid={grid} seeds={n_seeds} densities={densities}")
    out = {"mode": "full" if args.full else "smoke", "grid": grid,
           "densities": densities, "points": []}
    mus, los, his = [], [], []
    t0 = time.time()
    for dens in densities:
        burdens = [burden_at_density(grid, dens, s, total_ms, sample_every)
                   for s in range(n_seeds)]
        mu, lo, hi = bootstrap_ci(burdens, n_boot=2000 if args.full else 500, seed=0)
        out["points"].append({"density": dens, "burden_mean": mu,
                              "burden_ci95": [lo, hi]})
        mus.append(mu); los.append(lo); his.append(hi)
        print(f"  density={dens:.2f}  burden={mu:.0f} PS·ms  95%CI=[{lo:.0f}, {hi:.0f}]")

    fig, ax = plt.subplots(figsize=(6, 4))
    err = [np.array(mus) - np.array(los), np.array(his) - np.array(mus)]
    ax.errorbar(densities, mus, yerr=err, fmt="o-", capsize=4, color="#c1121f")
    ax.set_xlabel("fibrosis area fraction"); ax.set_ylabel("wavebreak burden (PS·ms)")
    ax.set_title("Critical fibrosis density for wavebreak (microgravity kinetics)")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fibrosis_threshold.png"), dpi=120)
    plt.close(fig)

    out["wall_seconds"] = round(time.time() - t0, 1)
    data = json.load(open(RESULTS)) if os.path.exists(RESULTS) else {}
    data["fibrosis_threshold"] = out
    json.dump(data, open(RESULTS, "w"), indent=2)
    print(f"[fibrosis-threshold] wall={out['wall_seconds']}s -> figures/fibrosis_threshold.png")


if __name__ == "__main__":
    main()
