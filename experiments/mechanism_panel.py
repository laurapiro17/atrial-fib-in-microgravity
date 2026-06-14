"""Mechanism-isolation panel: which microgravity remodelling driver dominates
re-entrant vulnerability?

We toggle each remodelling mechanism on its own and measure the wavebreak burden
(artifact-rejected, time-integrated phase-singularity count) with the same S1-S2
induction as the main ensemble:

  * ground        -- baseline kinetics, no fibrosis, nominal size
  * electrical    -- AF-type ionic remodelling only (short APD), no fibrosis
  * fibrosis      -- correlated low-coupling patches only, baseline kinetics
  * dilation      -- enlarged sheet only, baseline kinetics
  * combined      -- all three (== the microgravity condition)

For each mechanism we run a small fibrosis-seed ensemble and report the mean
wavebreak burden with a bootstrap 95% CI. This answers the reviewer's question:
is the vulnerability driven by the electrical remodelling, the fibrosis, the
dilation, or only their combination?

Run:  python experiments/mechanism_panel.py [--full]
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
from afib_microgravity.remodeling import (  # noqa: E402
    baseline_crn_params,
    microgravity_crn_params,
)

# Reuse the exact induction + metric definitions from the headline ensemble.
from ensemble_with_ci import (  # noqa: E402
    ARTIFACT_PS_THRESHOLD,
    MIN_ANALYSIS_MS,
    S2_DELAY_MS,
    crn_phase,
    s1_planar,
    s2_quadrant,
    window_metrics,
)

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
RESULTS = os.path.join(FIG, "results_crn.json")
DT = 0.02
D_LONG, D_TRANS, DX = 0.15, 0.05, 0.25
DILATION = 1.3
FIBROSIS_DENSITY, FIBROSIS_CORR = 0.3, 4.0
DELAY_MS = 30.0

MECHANISMS = ["ground", "electrical", "fibrosis", "dilation", "combined"]


def build(mech, base_grid, seed):
    """Return (cell, diffusion) for one isolated mechanism."""
    electrical = mech in ("electrical", "combined")
    fibrotic = mech in ("fibrosis", "combined")
    dilated = mech in ("dilation", "combined")

    g = int(base_grid * DILATION) if dilated else base_grid
    shape = (g, g)
    params = microgravity_crn_params(severity=1.0) if electrical else baseline_crn_params()
    coupling = (correlated_fibrosis(shape, density=FIBROSIS_DENSITY,
                                    corr_len=FIBROSIS_CORR, seed=seed)
                if fibrotic else None)
    cell = CRNCell(shape=shape, params=params, use_numba=True)
    diff = AnisotropicDiffusion(shape, D_LONG, D_TRANS, theta=np.zeros(shape),
                                dx=DX, coupling=coupling)
    return cell, diff


def run_one(mech, base_grid, seed, total_ms, sample_every):
    cell, diff = build(mech, base_grid, seed)
    sheet = MonodomainSheet(cell, diff, dt=DT)
    n_steps = round(total_ms / DT)
    s2_step = round(S2_DELAY_MS / DT)
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
    m = window_metrics(ps, t, S2_DELAY_MS, total_ms, cell.shape[0] * cell.shape[1])
    return m["break_burden_ps_ms"], m["peak_ps"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--grid", type=int, default=None)
    ap.add_argument("--seeds", type=int, default=None)
    ap.add_argument("--total-ms", type=float, default=None, dest="total_ms")
    ap.add_argument("--sample-every", type=int, default=None, dest="sample_every")
    args = ap.parse_args()

    if args.full:
        grid = args.grid or 160; n_seeds = args.seeds or 8
        total_ms = args.total_ms or 800.0; sample_every = args.sample_every or 100
    else:
        grid = args.grid or 60; n_seeds = args.seeds or 2
        total_ms = args.total_ms or 250.0; sample_every = args.sample_every or 50

    print(f"[mechanism] grid={grid} seeds={n_seeds} total_ms={total_ms}")
    out = {"mode": "full" if args.full else "smoke", "grid": grid,
           "n_seeds": n_seeds, "total_ms": total_ms, "mechanisms": {}}
    labels, mus, los, his = [], [], [], []
    t_start = time.time()
    for mech in MECHANISMS:
        burdens = []
        for seed in range(n_seeds):
            b, pk = run_one(mech, grid, seed, total_ms, sample_every)
            burdens.append(b)
        mu, lo, hi = bootstrap_ci(burdens, n_boot=2000 if args.full else 500, seed=0)
        out["mechanisms"][mech] = {"burden_mean": mu, "burden_ci95": [lo, hi],
                                   "per_seed_burden": burdens}
        labels.append(mech); mus.append(mu); los.append(lo); his.append(hi)
        print(f"  {mech:10s} burden={mu:.0f} PS·ms  95%CI=[{lo:.0f}, {hi:.0f}]")

    # bar chart with CI error bars
    fig, ax = plt.subplots(figsize=(6.5, 4))
    x = np.arange(len(labels))
    err = [np.array(mus) - np.array(los), np.array(his) - np.array(mus)]
    ax.bar(x, mus, yerr=err, capsize=4, color="#c1121f", alpha=0.85)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15)
    ax.set_ylabel("wavebreak burden (PS·ms)")
    ax.set_title("Which remodelling driver fragments wavefronts?")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "mechanism_panel.png"), dpi=120)
    plt.close(fig)

    out["wall_seconds"] = round(time.time() - t_start, 1)
    os.makedirs(FIG, exist_ok=True)
    data = json.load(open(RESULTS)) if os.path.exists(RESULTS) else {}
    data["mechanism_panel"] = out
    json.dump(data, open(RESULTS, "w"), indent=2)
    print(f"[mechanism] wall={out['wall_seconds']}s -> figures/mechanism_panel.png")


if __name__ == "__main__":
    main()
