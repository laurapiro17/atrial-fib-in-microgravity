"""Fibrosis-seed ensemble of re-entry burden: ground vs microgravity, with CI.

This is the headline experiment of the CRN stack. For each condition we run an
ensemble of fibrosis realisations; in every realisation we induce re-entry with a
broken-wavefront protocol, evolve the monodomain sheet, and track the number of
phase singularities (rotors) over time. Because the microgravity substrate is
random, we report the bootstrap confidence interval of the rotor *density* across
seeds -- not a single run.

Phase singularities need two state variables that wind around the rotor. The CRN
cell exposes only V, so we build a slow "recovery proxy" from a time-delayed copy
of V: ``phase = arctan2(V_now - V0, V_delayed - V0)``. The delay (~20-40 ms)
gives V and its lagged self a quarter-cycle offset, exactly what the phase-charge
counter needs. We keep a short ring buffer of past V frames for this.

Outputs:
  * ``figures/ensemble_ps.png``      -- rotor count vs time, per condition
  * ``figures/results_crn.json``     -- per-seed stats + bootstrap CI under
                                        ["ensemble"]

Run:  ``python experiments/ensemble_with_ci.py [--full]``
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

from afib_microgravity.metrics import (  # noqa: E402
    bootstrap_ci,
    count_phase_singularities,
    is_sustained_af,
    ps_density,
)
from afib_microgravity.model import MonodomainSheet  # noqa: E402
from afib_microgravity.remodeling import make_condition_crn  # noqa: E402

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
RESULTS = os.path.join(FIG, "results_crn.json")

V_DEPOL = 20.0          # mV: value written into stimulated cells
V0 = -40.0              # mV: phase-plane origin (mid-upstroke), rotor reference
DELAY_MS = 30.0         # time-delay for the recovery proxy


def crn_phase(V_now, V_delayed):
    """Phase from V and its time-delayed copy (CRN recovery proxy)."""
    return np.arctan2(V_now - V0, V_delayed - V0)


def induce_broken_wavefront(cell):
    """Seed one broken wavefront: depolarise the left half, but blank the top-left
    quadrant so the wavefront has a free end that curls into a rotor."""
    ny, nx = cell.shape
    mask = np.zeros((ny, nx), dtype=bool)
    mask[:, : nx // 2] = True          # left half depolarised
    mask[: ny // 2, : nx // 2] = False  # ... except the top-left quadrant
    cell.stimulate(mask, V_DEPOL)


def run_seed(condition, base_shape, seed, dt, n_steps, record_every):
    """One realisation -> (ps_series, shape). PS sampled every record_every steps
    using the time-delayed-V phase proxy."""
    cell, diff = make_condition_crn(condition, base_shape=base_shape, seed=seed)
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
            V_now = cell.V
            V_delayed = ring[0]
            ps = count_phase_singularities(crn_phase(V_now, V_delayed))
            ps_series.append(ps)
    return ps_series, cell.shape


def merge_results(section, payload):
    os.makedirs(FIG, exist_ok=True)
    data = {}
    if os.path.exists(RESULTS):
        with open(RESULTS) as f:
            data = json.load(f)
    data[section] = payload
    with open(RESULTS, "w") as f:
        json.dump(data, f, indent=2)


def second_half_mean(series):
    series = np.asarray(series, dtype=float)
    if series.size == 0:
        return 0.0
    return float(series[series.size // 2:].mean())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--grid", type=int, default=None,
                    help="base grid side (ny=nx); overrides mode default")
    ap.add_argument("--seeds", type=int, default=None,
                    help="number of fibrosis seeds per condition")
    ap.add_argument("--steps", type=int, default=None,
                    help="number of time steps to evolve")
    ap.add_argument("--dt", type=float, default=0.02)
    args = ap.parse_args()

    if args.full:
        grid = args.grid or 160
        n_seeds = args.seeds or 15
        n_steps = args.steps or 120000      # ~2400 ms at dt=0.02
        record_every = 200
        n_boot = 2000
    else:
        grid = args.grid or 48
        n_seeds = args.seeds or 2
        n_steps = args.steps or 400         # ~8 ms at dt=0.02 (smoke: runs, not AF)
        record_every = 20
        n_boot = 500

    base_shape = (grid, grid)
    seeds = list(range(n_seeds))
    conditions = ["ground", "microgravity"]

    print(f"[ensemble] mode={'full' if args.full else 'smoke'} grid={base_shape} "
          f"seeds={n_seeds} steps={n_steps} dt={args.dt}")

    results = {"mode": "full" if args.full else "smoke",
               "base_shape": list(base_shape), "n_seeds": n_seeds,
               "n_steps": n_steps, "dt": args.dt, "conditions": {}}
    plot_series = {}
    t_start = time.time()

    for cond in conditions:
        per_seed = []
        densities = []
        rep_series = None
        for seed in seeds:
            t0 = time.time()
            ps_series, shape = run_seed(
                cond, base_shape, seed, args.dt, n_steps, record_every)
            area = shape[0] * shape[1]
            ps_mean = second_half_mean(ps_series)
            ps_max = int(max(ps_series)) if ps_series else 0
            dens = ps_density(ps_mean, area)
            sustained = is_sustained_af(ps_series, window_frac=0.5, threshold=2)
            per_seed.append({
                "seed": seed, "shape": list(shape),
                "ps_mean_second_half": ps_mean, "ps_max": ps_max,
                "ps_density_x1e4": dens, "sustained_af": bool(sustained),
                "wall_seconds": round(time.time() - t0, 2),
            })
            densities.append(dens)
            if rep_series is None:
                rep_series = ps_series
            print(f"  [{cond} seed={seed}] ps_mean2nd={ps_mean:.2f} "
                  f"max={ps_max} dens={dens:.2f} sustained={sustained} "
                  f"({per_seed[-1]['wall_seconds']}s)")

        mean, lo, hi = bootstrap_ci(densities, n_boot=n_boot, alpha=0.05, seed=0)
        results["conditions"][cond] = {
            "per_seed": per_seed,
            "ps_density_mean_x1e4": mean,
            "ps_density_ci95_x1e4": [lo, hi],
            "n_sustained": int(sum(d["sustained_af"] for d in per_seed)),
        }
        plot_series[cond] = rep_series
        print(f"[{cond}] ps_density={mean:.2f} 95%CI=[{lo:.2f}, {hi:.2f}] "
              f"sustained={results['conditions'][cond]['n_sustained']}/{n_seeds}")

    # fold-change headline (guard against zero ground density)
    g = results["conditions"]["ground"]["ps_density_mean_x1e4"]
    u = results["conditions"]["microgravity"]["ps_density_mean_x1e4"]
    results["density_fold_increase"] = round(u / g, 3) if g > 1e-9 else None
    results["total_wall_seconds"] = round(time.time() - t_start, 1)

    os.makedirs(FIG, exist_ok=True)
    fig, axis = plt.subplots(figsize=(6, 4))
    colors = {"ground": "#1b4079", "microgravity": "#c1121f"}
    for cond, series in plot_series.items():
        if series:
            t = np.arange(len(series)) * record_every * args.dt
            axis.plot(t, series, color=colors[cond], label=f"{cond} (seed 0)")
    axis.set_xlabel("time (ms)")
    axis.set_ylabel("phase singularities")
    axis.set_title("Re-entry burden over time (representative seed)")
    axis.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "ensemble_ps.png"), dpi=120)
    plt.close(fig)

    merge_results("ensemble", results)
    print(f"[ensemble] total wall={results['total_wall_seconds']}s")
    print("Wrote figures/ensemble_ps.png + results_crn.json[ensemble]")


if __name__ == "__main__":
    main()
