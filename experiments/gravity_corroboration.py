"""Does the analytic law N_g actually predict rotors in the full sheet?

gravity_sweep.py builds N_g from a single-cell wavelength argument. This script
is the independent check: it runs the *2-D monodomain reaction-diffusion sheet*
at several gravity levels, with the tissue knobs (electrical severity, fibrosis
density, dilation) all driven by the same gravity_to_remodeling(g) map, induces a
broken wavefront, and measures the phase-singularity (rotor) density.

The law is corroborated if the measured rotor burden rises across the same
gravity region where N_g crosses 1 (i.e. high PS density for g < g*, low for
g > g*). This does NOT prove the law; it checks the single-cell argument and the
full PDE tell the same story.

Outputs:
  * figures/gravity_corroboration.png  -- PS density vs g, with N_g overlaid
  * figures/results_crn.json           -- under ["gravity_corroboration"]

Run:  python experiments/gravity_corroboration.py [--full]
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import argparse
import json
import os
import sys

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
from afib_microgravity.gravity_law import gravity_to_remodeling  # noqa: E402

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
RESULTS = os.path.join(FIG, "results_crn.json")

V_DEPOL = 20.0
V0 = -40.0
DELAY_MS = 30.0
SEED = 0


def crn_phase(V_now, V_delayed):
    return np.arctan2(V_now - V0, V_delayed - V0)


def induce_broken_wavefront(cell):
    ny, nx = cell.shape
    mask = np.zeros((ny, nx), dtype=bool)
    mask[:, : nx // 2] = True
    mask[: ny // 2, : nx // 2] = False
    cell.stimulate(mask, V_DEPOL)


def run_sheet_for_gravity(g, base_grid, dt, n_steps, record_every):
    rem = gravity_to_remodeling(g)
    n = int(base_grid * rem.dilation)
    shape = (n, n)
    cell = CRNCell(shape=shape, params=microgravity_crn_params(severity=rem.severity))
    coupling = correlated_fibrosis(shape, density=rem.density, corr_len=4.0, seed=SEED)
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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--dt", type=float, default=0.02)
    args = ap.parse_args()

    gs = [0.0, 0.16, 0.38, 0.5, 0.8, 1.0]
    if args.full:
        base_grid, n_steps, record_every = 100, 12000, 100
    else:
        base_grid, n_steps, record_every = 40, 400, 20

    psd = []
    for g in gs:
        d = run_sheet_for_gravity(g, base_grid, args.dt, n_steps, record_every)
        psd.append(d)
        print(f"  g={g:.2f}  ps_density={d:.2f}", flush=True)

    # overlay N_g from the sweep, if present
    ng = None
    if os.path.exists(RESULTS):
        with open(RESULTS) as f:
            data = json.load(f)
        if "gravity_law" in data:
            ngmap = {r["g"]: r["N_g"] for r in data["gravity_law"]["rows"]}
            ng = [ngmap.get(g) for g in gs]

    fig, ax1 = plt.subplots(figsize=(6.5, 4.2))
    ax1.plot(gs, psd, "o-", color="#d65a31", lw=2, label="rotor burden (PS density)")
    ax1.set_xlabel("gravitational level  g  (Earth-g units)")
    ax1.set_ylabel("PS density (/1e4 cells)", color="#d65a31")
    if ng is not None and all(v is not None for v in ng):
        ax2 = ax1.twinx()
        ax2.plot(gs, ng, "s--", color="#1b4079", lw=1.5, label=r"$\mathcal{N}_g$")
        ax2.axhline(1.0, color="#888", ls=":", lw=1)
        ax2.set_ylabel(r"$\mathcal{N}_g$", color="#1b4079")
    ax1.set_title("Sheet rotor burden vs the analytic law")
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, "gravity_corroboration.png"), dpi=130)
    plt.close(fig)

    data = {}
    if os.path.exists(RESULTS):
        with open(RESULTS) as f:
            data = json.load(f)
    data["gravity_corroboration"] = {
        "mode": "full" if args.full else "smoke",
        "base_grid": base_grid, "n_steps": n_steps,
        "g": gs, "ps_density_x1e4": psd,
    }
    with open(RESULTS, "w") as f:
        json.dump(data, f, indent=2)
    print(f"[gravity_corroboration] wrote figures/gravity_corroboration.png", flush=True)


if __name__ == "__main__":
    main()
