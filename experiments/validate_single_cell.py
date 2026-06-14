"""Single-cell CRN validation: pace one atrial cell and check the action
potential lands in physiological ranges.

This is the smoke-test / sanity gate for the whole CRN stack: if the kinetics
do not produce a physiological action potential here, nothing downstream
(restitution, sheets, ensembles) can be trusted. We pace a single CRN cell with
a brief stimulus, record the membrane trace, compute ``ap_biomarkers`` and assert
each biomarker sits in a literature-plausible window for human atrium.

Outputs:
  * ``figures/crn_ap.png``           -- V(t) of the paced action potential
  * ``figures/results_crn.json``     -- biomarkers under ["single_cell"]

Run:  ``python experiments/validate_single_cell.py [--full]``  (seconds, headless)
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

from afib_microgravity.crn import CRNCell, CRNParams, ap_biomarkers  # noqa: E402

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
RESULTS = os.path.join(FIG, "results_crn.json")

# Physiological acceptance windows for a CRN human-atrial action potential.
RANGES = {
    "V_rest": (-86.0, -75.0),
    "V_peak": (10.0, 35.0),
    "APD90": (250.0, 360.0),
    "dVdt_max": (100.0, np.inf),
}


def _scalar(V):
    return float(np.asarray(V).flat[0])


def pace_cell(dt, stim_amp, stim_dur, settle_ms, record_ms):
    """Pace one CRN cell once and return (times, volts) of the recorded beat."""
    cell = CRNCell(shape=(1, 1), params=CRNParams())

    # let the cell relax to its resting state before the stimulus
    for _ in range(round(settle_ms / dt)):
        dVdt = cell.reaction_step(dt, I_stim=0.0)
        cell.set_V(cell.V + dt * dVdt)

    ts = [0.0]
    vs = [_scalar(cell.V)]
    t = 0.0

    # stimulus window
    for _ in range(round(stim_dur / dt)):
        dVdt = cell.reaction_step(dt, I_stim=stim_amp)
        cell.set_V(cell.V + dt * dVdt)
        t += dt
        ts.append(t)
        vs.append(_scalar(cell.V))

    # free repolarisation
    for _ in range(round((record_ms - stim_dur) / dt)):
        dVdt = cell.reaction_step(dt, I_stim=0.0)
        cell.set_V(cell.V + dt * dVdt)
        t += dt
        ts.append(t)
        vs.append(_scalar(cell.V))

    return np.asarray(ts), np.asarray(vs)


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
    ap.add_argument("--full", action="store_true",
                    help="production config (finer dt, longer record); "
                         "default is the fast smoke config")
    args = ap.parse_args()

    if args.full:
        dt, record_ms, settle_ms = 0.02, 500.0, 200.0
    else:
        dt, record_ms, settle_ms = 0.05, 450.0, 100.0
    stim_amp, stim_dur = -2300.0, 2.0  # pA, ms

    t0 = time.time()
    ts, vs = pace_cell(dt, stim_amp, stim_dur, settle_ms, record_ms)
    bm = ap_biomarkers(ts, vs)
    wall = time.time() - t0

    # --- physiological self-check ---
    failures = []
    for key, (lo, hi) in RANGES.items():
        val = bm[key]
        if not (lo <= val <= hi):
            failures.append(f"{key}={val:.2f} outside [{lo}, {hi}]")

    print(f"[single_cell] biomarkers: {bm}")
    print(f"[single_cell] wall={wall:.1f}s mode={'full' if args.full else 'smoke'}")

    # --- figure ---
    os.makedirs(FIG, exist_ok=True)
    fig, axis = plt.subplots(figsize=(6, 4))
    axis.plot(ts, vs, color="#1b4079", lw=1.5)
    axis.axhline(bm["V_rest"], color="grey", ls=":", lw=0.8)
    axis.set_xlabel("time (ms)")
    axis.set_ylabel("V (mV)")
    axis.set_title(f"CRN atrial AP  (APD90={bm['APD90']:.0f} ms, "
                   f"dV/dt_max={bm['dVdt_max']:.0f} mV/ms)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "crn_ap.png"), dpi=120)
    plt.close(fig)

    merge_results("single_cell", {
        "mode": "full" if args.full else "smoke",
        "dt": dt,
        "biomarkers": bm,
        "ranges": {k: [lo, (None if np.isinf(hi) else hi)]
                   for k, (lo, hi) in RANGES.items()},
        "passed": not failures,
        "wall_seconds": round(wall, 2),
    })
    print(f"Wrote figures/crn_ap.png + results_crn.json[single_cell]")

    if failures:
        raise AssertionError("physiological check failed: " + "; ".join(failures))
    print("[single_cell] PASS: all biomarkers physiological")


if __name__ == "__main__":
    main()
