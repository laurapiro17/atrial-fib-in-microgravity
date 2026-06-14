"""Dose-response: how re-entry burden scales with each microgravity mechanism.

The microgravity condition bundles three changes -- atrial dilation (bigger
sheet), interstitial fibrosis (more low-coupling scar), and electrical
remodelling (reduced ICaL/Ito/IKur). This experiment sweeps each knob
independently around the microgravity operating point and records the
phase-singularity density of one short CRN sheet per value, isolating the
marginal contribution of each mechanism.

Three panels:
  * dilation factor          -- grid linear scale relative to the base shape
  * fibrosis density         -- scar area fraction in the coupling field
  * remodelling severity     -- ICaL/Ito/IKur reduction (0=baseline, 1=halved)

Outputs:
  * ``figures/dose_response.png``    -- 3 panels, ps_density vs knob value
  * ``figures/results_crn.json``     -- swept values under ["dose_response"]

Run:  ``python experiments/dose_response.py [--full]``
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
from afib_microgravity.metrics import (  # noqa: E402
    count_phase_singularities,
    ps_density,
)
from afib_microgravity.model import MonodomainSheet  # noqa: E402
from afib_microgravity.remodeling import microgravity_crn_params  # noqa: E402

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
RESULTS = os.path.join(FIG, "results_crn.json")

V_DEPOL = 20.0
V0 = -40.0
DELAY_MS = 30.0

# microgravity operating point (the value held fixed while the others sweep)
BASE_DILATION = 1.3
BASE_DENSITY = 0.3
BASE_SEVERITY = 1.0
SEED = 0


def crn_phase(V_now, V_delayed):
    return np.arctan2(V_now - V0, V_delayed - V0)


def induce_broken_wavefront(cell):
    ny, nx = cell.shape
    mask = np.zeros((ny, nx), dtype=bool)
    mask[:, : nx // 2] = True
    mask[: ny // 2, : nx // 2] = False
    cell.stimulate(mask, V_DEPOL)


def run_sheet(shape, density, severity, dt, n_steps, record_every):
    """Build a CRN sheet from explicit knobs and return its 2nd-half ps_density."""
    cell = CRNCell(shape=shape, params=microgravity_crn_params(severity=severity))
    coupling = correlated_fibrosis(shape, density=density, corr_len=4.0, seed=SEED)
    diff = AnisotropicDiffusion(shape=shape, d_long=0.06, d_trans=0.02,
                                theta=np.zeros(shape), dx=0.25, coupling=coupling)
    sheet = MonodomainSheet(cell, diff, dt=dt)
    induce_broken_wavefront(cell)

    delay_steps = max(1, round(DELAY_MS / dt))
    ring = [cell.V.copy()]
    ps_series = []
    for i in range(n_steps):
        sheet.step()
        ring.append(cell.V.copy())
        if len(ring) > delay_steps + 1:
            ring.pop(0)
        if i % record_every == 0:
            ps_series.append(count_phase_singularities(crn_phase(cell.V, ring[0])))
    if not ps_series:
        return 0.0
    ps_mean = float(np.asarray(ps_series[len(ps_series) // 2:], float).mean())
    return ps_density(ps_mean, shape[0] * shape[1])


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
    ap.add_argument("--dt", type=float, default=0.02)
    args = ap.parse_args()

    if args.full:
        base_grid = 140
        n_steps, record_every = 60000, 200
        dilations = [1.0, 1.15, 1.3, 1.45, 1.6]
        densities = [0.0, 0.15, 0.3, 0.45, 0.6]
        severities = [0.0, 0.25, 0.5, 0.75, 1.0]
    else:
        base_grid = 48
        n_steps, record_every = 400, 20
        dilations = [1.0, 1.4]
        densities = [0.0, 0.4]
        severities = [0.0, 1.0]

    def grid(dil):
        n = int(base_grid * dil)
        return (n, n)

    print(f"[dose_response] mode={'full' if args.full else 'smoke'} "
          f"base_grid={base_grid} steps={n_steps}")
    t_start = time.time()

    sweeps = {}

    # 1) dilation (fibrosis + severity held at base)
    vals = []
    for dil in dilations:
        d = run_sheet(grid(dil), BASE_DENSITY, BASE_SEVERITY,
                      args.dt, n_steps, record_every)
        vals.append(d)
        print(f"  [dilation={dil}] ps_density={d:.2f}")
    sweeps["dilation"] = {"values": dilations, "ps_density_x1e4": vals}

    # 2) fibrosis density (dilation + severity held at base)
    vals = []
    for dens in densities:
        d = run_sheet(grid(BASE_DILATION), dens, BASE_SEVERITY,
                      args.dt, n_steps, record_every)
        vals.append(d)
        print(f"  [density={dens}] ps_density={d:.2f}")
    sweeps["fibrosis_density"] = {"values": densities, "ps_density_x1e4": vals}

    # 3) remodelling severity (dilation + fibrosis held at base)
    vals = []
    for sev in severities:
        d = run_sheet(grid(BASE_DILATION), BASE_DENSITY, sev,
                      args.dt, n_steps, record_every)
        vals.append(d)
        print(f"  [severity={sev}] ps_density={d:.2f}")
    sweeps["remodelling_severity"] = {"values": severities, "ps_density_x1e4": vals}

    wall = time.time() - t_start

    panels = [
        ("dilation", "dilation factor"),
        ("fibrosis_density", "fibrosis density"),
        ("remodelling_severity", "remodelling severity"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    for axis, (key, label) in zip(axes, panels):
        axis.plot(sweeps[key]["values"], sweeps[key]["ps_density_x1e4"],
                  "o-", color="#1b4079")
        axis.set_xlabel(label)
        axis.set_ylabel("PS density (/1e4 cells)")
    fig.suptitle("Dose-response of re-entry burden to microgravity mechanisms")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "dose_response.png"), dpi=120)
    plt.close(fig)

    merge_results("dose_response", {
        "mode": "full" if args.full else "smoke",
        "dt": args.dt, "base_grid": base_grid, "n_steps": n_steps,
        "base_operating_point": {
            "dilation": BASE_DILATION, "fibrosis_density": BASE_DENSITY,
            "remodelling_severity": BASE_SEVERITY},
        "sweeps": sweeps,
        "wall_seconds": round(wall, 1),
    })
    print(f"[dose_response] wall={wall:.1f}s")
    print("Wrote figures/dose_response.png + results_crn.json[dose_response]")


if __name__ == "__main__":
    main()
