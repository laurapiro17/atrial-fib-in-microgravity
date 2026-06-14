"""APD90 restitution curves for ground vs microgravity-remodelled atrium.

Restitution -- APD90 as a function of the preceding diastolic interval (DI) --
is the canonical single-cell read-out of fibrillatory risk: a shorter APD and a
steeper slope at short DI favour wavebreak. We compute S1-S2 APD90 restitution
for the baseline CRN cell and for the microgravity-remodelled cell (reduced
ICaL/Ito/IKur) and overlay the two curves.

``restitution.apd_restitution`` builds a fresh ``CRNCell(CRNParams())`` per DI
internally and does not accept a params argument, so to vary kinetics we run it
once with the default (baseline) params, then re-run the same protocol inline
with the microgravity params via a small local copy of the pacing loop.

Outputs:
  * ``figures/restitution.png``      -- APD90 vs DI, both conditions
  * ``figures/results_crn.json``     -- curves under ["restitution"]

Run:  ``python experiments/restitution_curves.py [--full]``
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
_REPOL_MV = -70.0


def _scalar(V):
    return float(np.asarray(V).flat[0])


def _advance(cell, dt, n_steps, i_stim=0.0):
    for _ in range(int(n_steps)):
        dVdt = cell.reaction_step(dt, I_stim=i_stim)
        cell.set_V(cell.V + dt * dVdt)


def apd_restitution_params(diastolic_intervals, params, dt=0.05, n_s1=4,
                           s1_bcl=500.0, stim_amp=-2300.0, stim_dur=2.0):
    """S1-S2 APD90 restitution for a CRN cell built with the given params.

    Mirrors ``restitution.apd_restitution`` but lets the caller supply the
    kinetic parameters (the library function hard-codes baseline CRNParams).
    """
    dis = list(diastolic_intervals)
    apds = []
    for di in dis:
        cell = CRNCell(shape=(1, 1), params=params)
        for _ in range(n_s1):
            _advance(cell, dt, round(stim_dur / dt), i_stim=stim_amp)
            _advance(cell, dt, round((s1_bcl - stim_dur) / dt), i_stim=0.0)
        for _ in range(round(600.0 / dt)):
            dVdt = cell.reaction_step(dt, I_stim=0.0)
            cell.set_V(cell.V + dt * dVdt)
            if _scalar(cell.V) < _REPOL_MV:
                break
        _advance(cell, dt, round(di / dt), i_stim=0.0)

        ts = [0.0]
        vs = [_scalar(cell.V)]
        t = 0.0
        for _ in range(round(stim_dur / dt)):
            dVdt = cell.reaction_step(dt, I_stim=stim_amp)
            cell.set_V(cell.V + dt * dVdt)
            t += dt
            ts.append(t)
            vs.append(_scalar(cell.V))
        for _ in range(round(500.0 / dt)):
            dVdt = cell.reaction_step(dt, I_stim=0.0)
            cell.set_V(cell.V + dt * dVdt)
            t += dt
            ts.append(t)
            vs.append(_scalar(cell.V))
        apds.append(ap_biomarkers(ts, vs)["APD90"])
    return dis, apds


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
                    help="6 DIs and more S1 beats; default is 3-DI smoke config")
    args = ap.parse_args()

    if args.full:
        dis = [50, 100, 150, 250, 400, 600]
        dt, n_s1 = 0.05, 5
    else:
        dis = [80, 200, 400]
        dt, n_s1 = 0.1, 2

    t0 = time.time()
    _, apd_ground = apd_restitution_params(
        dis, baseline_crn_params(), dt=dt, n_s1=n_s1)
    _, apd_ug = apd_restitution_params(
        dis, microgravity_crn_params(severity=1.0), dt=dt, n_s1=n_s1)
    wall = time.time() - t0

    print(f"[restitution] DIs={dis}")
    print(f"[restitution] ground APD90={['%.0f' % a for a in apd_ground]}")
    print(f"[restitution] microg APD90={['%.0f' % a for a in apd_ug]}")
    print(f"[restitution] wall={wall:.1f}s mode={'full' if args.full else 'smoke'}")

    os.makedirs(FIG, exist_ok=True)
    fig, axis = plt.subplots(figsize=(6, 4))
    axis.plot(dis, apd_ground, "o-", color="#1b4079", label="ground")
    axis.plot(dis, apd_ug, "s--", color="#c1121f", label="microgravity")
    axis.set_xlabel("diastolic interval (ms)")
    axis.set_ylabel("APD90 (ms)")
    axis.set_title("CRN APD90 restitution")
    axis.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "restitution.png"), dpi=120)
    plt.close(fig)

    merge_results("restitution", {
        "mode": "full" if args.full else "smoke",
        "dt": dt,
        "diastolic_intervals": list(dis),
        "apd90_ground": [float(a) for a in apd_ground],
        "apd90_microgravity": [float(a) for a in apd_ug],
        "wall_seconds": round(wall, 2),
    })
    print("Wrote figures/restitution.png + results_crn.json[restitution]")


if __name__ == "__main__":
    main()
