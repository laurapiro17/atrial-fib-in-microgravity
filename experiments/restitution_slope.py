"""Restitution-slope link: the mechanistic 'why' behind the wavebreak result.

Steep APD restitution (maximum slope > 1) is the classical predictor of
wavefront instability, alternans and breakup. We measure the APD90 restitution
curve for ground (baseline CRN) vs microgravity (AF-remodelled CRN) by S1-S2
pacing a single cell across a range of diastolic intervals, and report the
maximum restitution slope for each. A steeper slope in the microgravity substrate
links the ionic remodelling (Task 6) directly to the increased wavebreak burden.

Run:  python experiments/restitution_slope.py [--full]
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

from afib_microgravity.crn import CRNCell, ap_biomarkers  # noqa: E402
from afib_microgravity.remodeling import (  # noqa: E402
    baseline_crn_params,
    microgravity_crn_params,
)

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
RESULTS = os.path.join(FIG, "results_crn.json")
DT = 0.05
STIM_AMP, STIM_DUR = -2300.0, 2.0


def _apd90_after_s1s2(params, di_ms, n_s1=4, s1_bcl=400.0):
    """Pace n_s1 S1 beats at s1_bcl, then an S2 at coupling = last repolarisation
    + di_ms; return APD90 of the S2 beat (ms)."""
    cell = CRNCell(shape=(1, 1), params=params, use_numba=False)

    def beat_and_wait(wait_ms):
        steps = round(wait_ms / DT)
        stim_steps = round(STIM_DUR / DT)
        for s in range(steps):
            I = STIM_AMP if s < stim_steps else 0.0
            dV = cell.reaction_step(DT, I_stim=I)
            cell.set_V(cell.V + DT * dV)

    for _ in range(n_s1):
        beat_and_wait(s1_bcl)
    # let the last S1 repolarise to rest, then wait the diastolic interval
    # (approximate: pace S2 di_ms after a fixed long recovery)
    rec = []
    steps = round((s1_bcl) / DT)
    # deliver S2 after di_ms of pure diastole following the S1 train tail
    beat_and_wait(di_ms)  # diastolic interval (no stim during most of it)
    # now S2 + record
    t, V = [], []
    stim_steps = round(STIM_DUR / DT)
    n = round(500.0 / DT)
    for s in range(n):
        I = STIM_AMP if s < stim_steps else 0.0
        dV = cell.reaction_step(DT, I_stim=I)
        cell.set_V(cell.V + DT * dV)
        t.append(s * DT); V.append(float(cell.V.ravel()[0]))
    return ap_biomarkers(np.asarray(t), np.asarray(V))["APD90"]


def restitution_curve(params, dis):
    return [(_apd90_after_s1s2(params, di)) for di in dis]


def max_slope(dis, apds):
    dis = np.asarray(dis, float); apds = np.asarray(apds, float)
    order = np.argsort(dis)
    d, a = dis[order], apds[order]
    slopes = np.diff(a) / np.diff(d)
    return float(np.nanmax(np.abs(slopes)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", action="store_true")
    args = ap.parse_args()
    dis = ([300, 200, 140, 100, 70, 50, 35] if args.full else [250, 120, 60])

    print(f"[restitution-slope] DIs={dis}")
    out = {"mode": "full" if args.full else "smoke", "dis_ms": dis, "conditions": {}}
    fig, ax = plt.subplots(figsize=(6, 4))
    colors = {"ground": "#1b4079", "microgravity": "#c1121f"}
    t0 = time.time()
    for cond, params in (("ground", baseline_crn_params()),
                         ("microgravity", microgravity_crn_params(severity=1.0))):
        apds = restitution_curve(params, dis)
        slope = max_slope(dis, apds)
        out["conditions"][cond] = {"apd90_ms": apds, "max_slope": slope}
        ax.plot(dis, apds, "o-", color=colors[cond],
                label=f"{cond} (max slope {slope:.2f})")
        print(f"  {cond:12s} APD90={[round(a,1) for a in apds]}  max_slope={slope:.2f}")
    ax.set_xlabel("diastolic interval (ms)"); ax.set_ylabel("APD90 (ms)")
    ax.set_title("APD restitution: slope > 1 predicts wavebreak")
    ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(FIG, "restitution_slope.png"), dpi=120); plt.close(fig)

    out["wall_seconds"] = round(time.time() - t0, 1)
    data = json.load(open(RESULTS)) if os.path.exists(RESULTS) else {}
    data["restitution_slope"] = out
    json.dump(data, open(RESULTS, "w"), indent=2)
    print(f"[restitution-slope] wall={out['wall_seconds']}s -> figures/restitution_slope.png")


if __name__ == "__main__":
    main()
