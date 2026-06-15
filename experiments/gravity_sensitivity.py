"""Sensitivity of the critical gravity g* to the modelling assumptions.

The headline N_g(g) curve rests on two assumptions a reviewer will attack:

  1. the characteristic atrial path length L0 (we used 8 cm);
  2. the *linear* g -> remodelling map (severity = 1 - g).

This script quantifies how much g* moves when those assumptions change, WITHOUT
re-running the ionic model. It reuses the (severity, APD90) pairs already measured
by gravity_sweep.py (stored in figures/results_crn.json) and interpolates APD90 as
a function of severity; everything else is analytic.

  * L0 sweep: g* for L0 in {6, 8, 10, 12} cm with the linear map.
  * map-shape sweep: severity = (1 - g)**p for p in {0.5, 1.0, 2.0}
    (concave / linear / convex fluid-shift drive), with L0 = 8 cm.

Outputs:
  * figures/gravity_sensitivity.png   -- two panels of g* vs assumption
  * figures/results_crn.json          -- under ["gravity_sensitivity"]

Run (after gravity_sweep.py):  python experiments/gravity_sensitivity.py
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from afib_microgravity.gravity_law import (  # noqa: E402
    DILATION_MAX,
    cardiac_gravitational_number,
    interpolate_crossing,
)

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
RESULTS = os.path.join(FIG, "results_crn.json")
CV_CM_S = 58.0


def load_apd_of_severity():
    """Return (severities, apd90s) sorted by severity, from the sweep results."""
    with open(RESULTS) as f:
        data = json.load(f)
    rows = data["gravity_law"]["rows"]
    pairs = sorted({(r["severity"], r["apd90_ms"]) for r in rows})
    sev = np.array([p[0] for p in pairs])
    apd = np.array([p[1] for p in pairs])
    return sev, apd


def ng_curve(gs, severity_of_g, apd_of_sev, l0_cm):
    """N_g over a g grid given a map g->severity and an APD90(severity) interp."""
    ng = []
    for g in gs:
        sev = severity_of_g(g)
        apd = float(np.interp(sev, apd_of_sev[0], apd_of_sev[1]))
        dilation = 1.0 + (DILATION_MAX - 1.0) * (1.0 - g)
        ng.append(cardiac_gravitational_number(apd, dilation=dilation,
                                               l0_cm=l0_cm, cv_cm_s=CV_CM_S))
    return ng


def main():
    sev_tab, apd_tab = load_apd_of_severity()
    apd_of_sev = (sev_tab, apd_tab)
    gs = list(np.linspace(0.0, 1.0, 51))

    # 1) L0 sweep (linear map)
    l0s = [6.0, 8.0, 10.0, 12.0]
    gstar_l0 = []
    for l0 in l0s:
        ng = ng_curve(gs, lambda g: 1.0 - g, apd_of_sev, l0)
        gstar_l0.append(interpolate_crossing(gs, ng, target=1.0))

    # 2) map-shape sweep (L0 = 8)
    ps = [0.5, 1.0, 2.0]
    gstar_p = []
    for p in ps:
        ng = ng_curve(gs, lambda g, p=p: (1.0 - g) ** p, apd_of_sev, 8.0)
        gstar_p.append(interpolate_crossing(gs, ng, target=1.0))

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    axes[0].plot(l0s, [g if g is not None else np.nan for g in gstar_l0],
                 "o-", color="#1b4079")
    axes[0].axhspan(0.16, 0.16, color="#999")
    for y, name in ((0.16, "Moon"), (0.38, "Mars")):
        axes[0].axhline(y, ls=":", color="#444", lw=0.8)
        axes[0].text(l0s[-1], y, name, va="bottom", ha="right", fontsize=8)
    axes[0].set_xlabel("atrial path length L0 (cm)")
    axes[0].set_ylabel("critical gravity g*")
    axes[0].set_title("g* vs L0 (linear map)")

    axes[1].plot(ps, [g if g is not None else np.nan for g in gstar_p],
                 "s-", color="#d65a31")
    for y, name in ((0.16, "Moon"), (0.38, "Mars")):
        axes[1].axhline(y, ls=":", color="#444", lw=0.8)
        axes[1].text(ps[-1], y, name, va="bottom", ha="right", fontsize=8)
    axes[1].set_xlabel("fluid-shift drive exponent p  ((1-g)^p)")
    axes[1].set_ylabel("critical gravity g*")
    axes[1].set_title("g* vs map shape (L0=8 cm)")

    fig.suptitle("Sensitivity of critical gravity g* to assumptions")
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "gravity_sensitivity.png"), dpi=130)
    plt.close(fig)

    with open(RESULTS) as f:
        data = json.load(f)
    data["gravity_sensitivity"] = {
        "l0_cm": l0s, "gstar_vs_l0": gstar_l0,
        "drive_exponent": ps, "gstar_vs_exponent": gstar_p,
        "note": "g* interpolated from APD90(severity) measured in gravity_law sweep; no re-simulation.",
    }
    with open(RESULTS, "w") as f:
        json.dump(data, f, indent=2)

    print("[gravity_sensitivity] g* vs L0:", dict(zip(l0s, gstar_l0)), flush=True)
    print("[gravity_sensitivity] g* vs p :", dict(zip(ps, gstar_p)), flush=True)
    print("Wrote figures/gravity_sensitivity.png", flush=True)


if __name__ == "__main__":
    main()
