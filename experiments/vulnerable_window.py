"""Vulnerable-window width: the gold-standard electrophysiological measure of
re-entrant vulnerability.

We sweep the S1-S2 coupling interval and, for each value, measure whether the S2
induces wavebreak (artifact-rejected peak phase-singularity count > 0). The set of
coupling intervals that induce re-entry is the *vulnerable window*; its width
quantifies how easy the substrate is to push into re-entry. A wider window in the
microgravity substrate than on ground is direct evidence of increased vulnerability.

Run:  python experiments/vulnerable_window.py [--full]
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

from afib_microgravity.metrics import count_phase_singularities  # noqa: E402
from afib_microgravity.model import MonodomainSheet  # noqa: E402
from afib_microgravity.remodeling import make_condition_crn  # noqa: E402

from ensemble_with_ci import (  # noqa: E402
    crn_phase,
    s1_planar,
    s2_quadrant,
    window_metrics,
)

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
RESULTS = os.path.join(FIG, "results_crn.json")
DT = 0.02
DELAY_MS = 30.0


def peak_break_for_s2(condition, grid, seed, s2_delay_ms, total_ms, sample_every):
    cell, diff = make_condition_crn(condition, base_shape=(grid, grid), seed=seed)
    sheet = MonodomainSheet(cell, diff, dt=DT)
    n_steps = round(total_ms / DT)
    s2_step = round(s2_delay_ms / DT)
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
    m = window_metrics(ps, t, s2_delay_ms, total_ms, cell.shape[0] * cell.shape[1])
    return m["peak_ps"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--grid", type=int, default=None)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    if args.full:
        grid = args.grid or 160
        s2_delays = list(range(90, 241, 15))   # ms
        total_ms = 600.0; sample_every = 100
    else:
        grid = args.grid or 60
        s2_delays = [120, 160, 200]
        total_ms = 250.0; sample_every = 50

    print(f"[vuln-window] grid={grid} s2_delays={s2_delays}")
    out = {"mode": "full" if args.full else "smoke", "grid": grid,
           "s2_delays_ms": s2_delays, "conditions": {}}
    t_start = time.time()
    fig, ax = plt.subplots(figsize=(6.5, 4))
    colors = {"ground": "#1b4079", "microgravity": "#c1121f"}

    for cond in ("ground", "microgravity"):
        peaks = []
        for d in s2_delays:
            pk = peak_break_for_s2(cond, grid, args.seed, d, total_ms, sample_every)
            peaks.append(int(pk))
        vulnerable = [d for d, p in zip(s2_delays, peaks) if p > 0]
        width = (max(vulnerable) - min(vulnerable)) if vulnerable else 0
        out["conditions"][cond] = {"peak_ps_by_delay": peaks,
                                   "vulnerable_delays_ms": vulnerable,
                                   "window_width_ms": width}
        ax.plot(s2_delays, peaks, "o-", color=colors[cond], label=cond)
        print(f"  {cond:12s} peaks={peaks}  vulnerable@{vulnerable}  width={width}ms")

    ax.set_xlabel("S1-S2 coupling interval (ms)")
    ax.set_ylabel("peak wavebreak (phase singularities)")
    ax.set_title("Vulnerable window: re-entry induction vs S2 timing")
    ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "vulnerable_window.png"), dpi=120)
    plt.close(fig)

    out["wall_seconds"] = round(time.time() - t_start, 1)
    data = json.load(open(RESULTS)) if os.path.exists(RESULTS) else {}
    data["vulnerable_window"] = out
    json.dump(data, open(RESULTS, "w"), indent=2)
    gw = out["conditions"]["ground"]["window_width_ms"]
    mw = out["conditions"]["microgravity"]["window_width_ms"]
    print(f"[vuln-window] ground width={gw}ms  microgravity width={mw}ms  "
          f"wall={out['wall_seconds']}s")


if __name__ == "__main__":
    main()
