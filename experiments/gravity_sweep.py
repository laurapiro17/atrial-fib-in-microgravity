"""Gravity sweep: the Cardiac Gravitational Number N_g vs gravitational level.

For each gravity level g we map g -> remodelling severity/dilation, measure the
single-cell APD90 under microgravity_crn_params(severity), form the electrical
wavelength WL = CV*APD90 and the dimensionless N_g = L0*dilation/WL, then locate
the critical gravity g* where N_g crosses 1. Markers for Moon (0.16), Mars
(0.38) and interplanetary transit (0.0) are overlaid.

Outputs:
  * figures/gravity_law.png        -- N_g vs g, with g* and body markers
  * figures/results_crn.json       -- under ["gravity_law"]

Run:  python experiments/gravity_sweep.py [--full]
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import argparse
import dataclasses
import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from afib_microgravity.crn import CRNCell, CRNParams, ap_biomarkers  # noqa: E402
from afib_microgravity.gravity_law import (  # noqa: E402
    gravity_to_remodeling,
    cardiac_gravitational_number,
    interpolate_crossing,
)

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
RESULTS = os.path.join(FIG, "results_crn.json")

DT = 0.02
CV_CM_S = 58.0          # planar CV (calibrated; g_Na unchanged across conditions)
L0_CM = 8.0             # ground atrial characteristic path length
BODIES = {"transit": 0.0, "Moon": 0.16, "Mars": 0.38, "Earth": 1.0}


def apd90_over_severities(severities, bcl=1000.0, n_beats=6):
    """APD90 (ms) for every severity in one vectorised pacing run.

    Each severity becomes one column of a (1, K) sheet with its own AF-remodelling
    conductances (same mapping as remodeling.microgravity_crn_params). The CRN
    NumPy reaction step broadcasts the per-column scalar conductances, so K
    severities cost essentially the same as a single cell. use_numba=False because
    the Numba kernel takes scalar conductances; the NumPy path takes arrays.
    """
    severities = list(severities)
    K = len(severities)
    sev = np.asarray(severities, dtype=float).reshape(1, K)
    base = CRNParams()
    params = dataclasses.replace(
        base,
        g_CaL=base.g_CaL * (1.0 - 0.70 * sev),
        g_to=base.g_to * (1.0 - 0.50 * sev),
        g_Kur_scale=base.g_Kur_scale * (1.0 - 0.50 * sev),
        g_K1=base.g_K1 * (1.0 + 1.0 * sev),
    )
    cell = CRNCell(shape=(1, K), params=params, use_numba=False)
    spb = int(bcl / DT)
    ss = int(2.0 / DT)
    tl, Vtr = [], []
    for beat in range(n_beats):
        for s in range(spb):
            I = -2300.0 if s < ss else 0.0
            dV = cell.reaction_step(DT, I_stim=I)
            cell.set_V(cell.V + DT * dV)
            if beat == n_beats - 1:
                Vtr.append(cell.V.reshape(-1).copy())
                tl.append(s * DT)
    t = np.asarray(tl)
    Vtr = np.asarray(Vtr)  # (T, K)
    return [ap_biomarkers(t, Vtr[:, k])["APD90"] for k in range(K)]


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
    ap.add_argument("--full", action="store_true")
    args = ap.parse_args()

    # The vectorised pacing cost is independent of the number of g-points, so we
    # always use a dense grid; --full only paces more beats for tighter steady state.
    gs = ([round(x, 3) for x in np.linspace(0.0, 1.0, 11)] if args.full
          else [0.0, 0.16, 0.38, 0.5, 0.8, 1.0])
    n_beats = 10 if args.full else 6

    severities = [gravity_to_remodeling(g).severity for g in gs]
    apds = apd90_over_severities(severities, n_beats=n_beats)

    rows = []
    for g, apd in zip(gs, apds):
        rem = gravity_to_remodeling(g)
        ng = cardiac_gravitational_number(apd, dilation=rem.dilation,
                                          l0_cm=L0_CM, cv_cm_s=CV_CM_S)
        rows.append({"g": g, "severity": rem.severity, "dilation": rem.dilation,
                     "apd90_ms": apd, "wavelength_cm": CV_CM_S * apd / 1000.0,
                     "N_g": ng})
        print(f"  g={g:.2f}  sev={rem.severity:.2f}  APD90={apd:.1f}ms  N_g={ng:.2f}",
              flush=True)

    gstar = interpolate_crossing([r["g"] for r in rows],
                                 [r["N_g"] for r in rows], target=1.0)

    # figure
    fig, axis = plt.subplots(figsize=(6.5, 4.2))
    axis.plot([r["g"] for r in rows], [r["N_g"] for r in rows],
              "o-", color="#1b4079", lw=2, label=r"$\mathcal{N}_g$")
    axis.axhline(1.0, color="#888", ls="--", lw=1)
    axis.fill_between([r["g"] for r in rows], 1.0,
                      [max(r["N_g"], 1.0) for r in rows],
                      color="#d65a31", alpha=0.15, label=r"$\mathcal{N}_g>1$ (vulnerable)")
    for name, gval in BODIES.items():
        axis.axvline(gval, color="#444", ls=":", lw=0.8)
        axis.text(gval, axis.get_ylim()[1], name, rotation=90,
                  va="top", ha="right", fontsize=8, color="#444")
    if gstar is not None:
        axis.plot([gstar], [1.0], "*", color="#e84545", ms=15,
                  label=fr"$g^*$={gstar:.2f}")
    axis.set_xlabel("gravitational level  g  (Earth-g units)")
    axis.set_ylabel(r"Cardiac Gravitational Number  $\mathcal{N}_g$")
    axis.set_title("Gravitational scaling of atrial re-entry vulnerability")
    axis.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "gravity_law.png"), dpi=130)
    plt.close(fig)

    merge_results("gravity_law", {
        "mode": "full" if args.full else "smoke",
        "cv_cm_s": CV_CM_S, "l0_cm": L0_CM,
        "g_star": gstar, "rows": rows,
    })
    print(f"[gravity_sweep] g*={gstar}  wrote figures/gravity_law.png", flush=True)


if __name__ == "__main__":
    main()
