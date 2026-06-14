"""Fibrosis-seed ensemble of TRANSIENT re-entrant vulnerability: ground vs
microgravity, with bootstrap CIs.

This is the headline experiment of the CRN stack. The scientific claim is
*transient re-entrant vulnerability*, not sustained AF: an S1-S2 cross-field
stimulus seeds a broken wavefront; in the microgravity substrate the resulting
rotor lives measurably longer before self-terminating, whereas on ground it
extinguishes almost immediately. We do NOT claim sustained fibrillation.

For each condition and each fibrosis seed we:
  1. build (cell, diffusion) via ``make_condition_crn``;
  2. induce re-entry with the S1-planar + S2-quadrant cross-field protocol
     (the exact induction reused from ``experiments/sustain_probe.py``);
  3. evolve the monodomain sheet for ``total_ms``, sampling phase singularities
     (PS) via the time-delayed-V phase proxy;
  4. discard the first ~100 ms after S2 (the sharp S1/S2 wavefronts spike the
     phase proxy -- an artifact) and compute metrics on the analysis window
     [S2 + 100 ms, end].

Per-seed metrics over the analysis window: peak PS, time-to-extinction, mean PS,
peak PS density, and (for completeness) the sustained-AF flag (expected False).
Per condition we bootstrap the mean + 95% CI of mean-PS and time-to-extinction,
and report the microgravity/ground fold-increase in mean PS.

Outputs:
  * ``figures/ensemble_ps.png``            -- PS(t) mean +- band, window shaded
  * ``figures/snapshot_crn_ground.png``    -- late-window V field (quiescent)
  * ``figures/snapshot_crn_microgravity.png`` -- late-window V field (rotor)
  * ``figures/results_crn.json``           -- everything, under ["ensemble"]

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
D_LONG = 0.15           # calibrated planar CV ~58 cm/s (see make_condition_crn)

# S1-S2 cross-field protocol timing (ms). S2 is delivered s2_delay_ms after S1.
S2_DELAY_MS = 160.0
# Blank-out window after S2 before metrics start (ms): the S1/S2 wavefronts spike
# the phase proxy and must be excluded from the analysis window.
BLANK_AFTER_S2_MS = 100.0

APD_NOTE = ("microgravity = AF-remodelled CRN, APD90 ~135 ms; ground = baseline "
            "CRN, APD90 ~294 ms")


def crn_phase(V_now, V_delayed):
    """Phase from V and its time-delayed copy (CRN recovery proxy)."""
    return np.arctan2(V_now - V0, V_delayed - V0)


def s1_planar(cell, width=4):
    """S1: planar wavefront from the left edge."""
    ny, nx = cell.shape
    mask = np.zeros((ny, nx), dtype=bool)
    mask[:, :width] = True
    cell.stimulate(mask, V_DEPOL)


def s2_quadrant(cell):
    """S2 over the lower-left quadrant: it can only propagate into recovered
    tissue, breaking the wavefront -> a rotor (the sustain_probe induction)."""
    ny, nx = cell.shape
    mask = np.zeros((ny, nx), dtype=bool)
    mask[ny // 2:, : nx // 2] = True
    cell.stimulate(mask, V_DEPOL)


def run_seed(condition, base_shape, seed, dt, total_ms, sample_every):
    """One realisation. Returns ``(ps_series, times_ms, s2_ms, shape, V_final)``.

    S1 planar, wait S2_DELAY_MS, S2 quadrant, then evolve to total_ms sampling PS
    every ``sample_every`` steps with the time-delayed-V phase proxy.
    """
    cell, diff = make_condition_crn(condition, base_shape=base_shape, seed=seed)
    sheet = MonodomainSheet(cell, diff, dt=dt)

    n_steps = round(total_ms / dt)
    s2_step = round(S2_DELAY_MS / dt)
    # The phase proxy needs V and a ~DELAY_MS-old copy. We snapshot only at the
    # sample cadence and read the snapshot ``delay_snaps`` records back, so the
    # ring holds O(DELAY/sample) frames -- a small, fixed footprint.
    record_every_ms = sample_every * dt
    delay_snaps = max(1, round(DELAY_MS / record_every_ms))

    ring = []
    ps_series, times = [], []

    s1_planar(cell)
    for i in range(n_steps):
        if i == s2_step:
            s2_quadrant(cell)
        sheet.step()
        if i % sample_every == 0:
            ring.append(cell.V.copy())
            if len(ring) > delay_snaps + 1:
                ring.pop(0)
            V_delayed = ring[0]
            ps = count_phase_singularities(crn_phase(cell.V, V_delayed))
            ps_series.append(ps)
            times.append(i * dt)

    return (np.asarray(ps_series, dtype=float), np.asarray(times, dtype=float),
            S2_DELAY_MS, cell.shape, cell.V.copy())


def window_metrics(ps_series, times, s2_ms, total_ms, area):
    """Transient-vulnerability metrics over the analysis window
    [s2_ms + BLANK_AFTER_S2_MS, total_ms]."""
    win_start = s2_ms + BLANK_AFTER_S2_MS
    in_win = times >= win_start
    w_ps = ps_series[in_win]
    w_t = times[in_win]

    if w_ps.size == 0:
        return {
            "peak_ps": 0, "mean_ps_window": 0.0,
            "time_to_extinction_ms": float(total_ms - s2_ms),
            "ps_density_peak_x1e4": 0.0, "sustained_af": False,
            "window_start_ms": win_start, "window_n": 0,
        }

    peak = int(w_ps.max())
    mean_ps = float(w_ps.mean())

    # time_to_extinction: time from S2 until PS first hits 0 and stays 0 to end.
    tte = float(total_ms - s2_ms)  # never extinguishes within window
    zero = w_ps == 0
    for k in range(w_ps.size):
        if zero[k] and zero[k:].all():
            tte = float(w_t[k] - s2_ms)
            break

    return {
        "peak_ps": peak,
        "mean_ps_window": mean_ps,
        "time_to_extinction_ms": tte,
        "ps_density_peak_x1e4": ps_density(peak, area),
        "sustained_af": bool(is_sustained_af(w_ps, window_frac=0.5, threshold=2)),
        "window_start_ms": win_start,
        "window_n": int(w_ps.size),
    }


def merge_results(section, payload):
    os.makedirs(FIG, exist_ok=True)
    data = {}
    if os.path.exists(RESULTS):
        with open(RESULTS) as f:
            data = json.load(f)
    data[section] = payload
    with open(RESULTS, "w") as f:
        json.dump(data, f, indent=2)


def save_snapshot(V, condition, total_ms, fname):
    fig, ax = plt.subplots(figsize=(4.2, 4))
    im = ax.imshow(V, cmap="inferno", vmin=-85, vmax=20, origin="lower")
    ax.set_title(f"{condition}: V at t={total_ms:.0f} ms")
    fig.colorbar(im, ax=ax, label="V (mV)", shrink=0.8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, fname), dpi=120)
    plt.close(fig)


def save_ps_figure(curves, s2_ms, total_ms):
    """curves[cond] = list of ps_series arrays (one per seed), all same length."""
    os.makedirs(FIG, exist_ok=True)
    colors = {"ground": "#1b4079", "microgravity": "#c1121f"}
    fig, ax = plt.subplots(figsize=(6.5, 4))
    win_start = s2_ms + BLANK_AFTER_S2_MS
    ax.axvspan(win_start, total_ms, color="0.85", alpha=0.5,
               label=f"analysis window (>{win_start:.0f} ms)")
    ax.axvline(s2_ms, color="0.4", ls="--", lw=1, label=f"S2 ({s2_ms:.0f} ms)")
    for cond, data in curves.items():
        t = data["times"]
        stack = np.vstack(data["series"])
        mean = stack.mean(axis=0)
        ax.plot(t, mean, color=colors.get(cond, "k"), label=f"{cond} (mean)")
        if stack.shape[0] > 1:
            sd = stack.std(axis=0)
            ax.fill_between(t, mean - sd, mean + sd, color=colors.get(cond, "k"),
                            alpha=0.18)
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("phase singularities")
    ax.set_title("Transient re-entrant vulnerability: PS(t) after S1-S2")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "ensemble_ps.png"), dpi=120)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--grid", type=int, default=None,
                    help="base grid side (ny=nx); overrides mode default")
    ap.add_argument("--seeds", type=int, default=None,
                    help="number of fibrosis seeds per condition")
    ap.add_argument("--total-ms", type=float, default=None, dest="total_ms",
                    help="total ms to evolve")
    ap.add_argument("--dt", type=float, default=None)
    ap.add_argument("--sample-every", type=int, default=None, dest="sample_every",
                    help="sample PS every N steps")
    args = ap.parse_args()

    if args.full:
        grid = args.grid or 180
        n_seeds = args.seeds or 12
        total_ms = args.total_ms or 1000.0
        dt = args.dt or 0.02
        sample_every = args.sample_every or 200
        n_boot = 2000
    else:
        grid = args.grid or 64
        n_seeds = args.seeds or 2
        total_ms = args.total_ms or 300.0
        dt = args.dt or 0.02
        sample_every = args.sample_every or 50
        n_boot = 500

    base_shape = (grid, grid)
    seeds = list(range(n_seeds))
    conditions = ["ground", "microgravity"]

    print(f"[ensemble] mode={'full' if args.full else 'smoke'} grid={base_shape} "
          f"seeds={n_seeds} total_ms={total_ms} dt={dt} sample_every={sample_every}")

    results = {
        "mode": "full" if args.full else "smoke",
        "claim": "transient re-entrant vulnerability (NOT sustained AF)",
        "base_shape": list(base_shape), "n_seeds": n_seeds,
        "total_ms": total_ms, "dt": dt, "sample_every": sample_every,
        "s2_delay_ms": S2_DELAY_MS, "blank_after_s2_ms": BLANK_AFTER_S2_MS,
        "d_long": D_LONG, "apd_note": APD_NOTE,
        "conditions": {},
    }
    curves = {}
    snapshots = {"ground": "snapshot_crn_ground.png",
                 "microgravity": "snapshot_crn_microgravity.png"}
    t_start = time.time()

    for cond in conditions:
        per_seed = []
        means, ttes = [], []
        series_list, common_times = [], None
        last_V, last_shape = None, None
        for seed in seeds:
            t0 = time.time()
            ps_series, times, s2_ms, shape, V_final = run_seed(
                cond, base_shape, seed, dt, total_ms, sample_every)
            area = shape[0] * shape[1]
            m = window_metrics(ps_series, times, s2_ms, total_ms, area)
            m["seed"] = seed
            m["shape"] = list(shape)
            m["wall_seconds"] = round(time.time() - t0, 2)
            per_seed.append(m)
            means.append(m["mean_ps_window"])
            ttes.append(m["time_to_extinction_ms"])
            series_list.append(ps_series)
            common_times = times
            last_V, last_shape = V_final, shape
            print(f"  [{cond} seed={seed}] peak={m['peak_ps']} "
                  f"mean_win={m['mean_ps_window']:.2f} "
                  f"tte={m['time_to_extinction_ms']:.0f}ms "
                  f"dens_peak={m['ps_density_peak_x1e4']:.2f} "
                  f"sustained={m['sustained_af']} ({m['wall_seconds']}s)")

        mean_ps_mu, mean_ps_lo, mean_ps_hi = bootstrap_ci(
            means, n_boot=n_boot, alpha=0.05, seed=0)
        tte_mu, tte_lo, tte_hi = bootstrap_ci(
            ttes, n_boot=n_boot, alpha=0.05, seed=0)
        results["conditions"][cond] = {
            "per_seed": per_seed,
            "mean_ps_window_mean": mean_ps_mu,
            "mean_ps_window_ci95": [mean_ps_lo, mean_ps_hi],
            "time_to_extinction_ms_mean": tte_mu,
            "time_to_extinction_ms_ci95": [tte_lo, tte_hi],
            "n_sustained": int(sum(d["sustained_af"] for d in per_seed)),
        }
        curves[cond] = {"times": common_times, "series": series_list}
        save_snapshot(last_V, cond, total_ms, snapshots[cond])
        print(f"[{cond}] mean_ps_window={mean_ps_mu:.2f} "
              f"95%CI=[{mean_ps_lo:.2f}, {mean_ps_hi:.2f}]  "
              f"tte={tte_mu:.0f}ms 95%CI=[{tte_lo:.0f}, {tte_hi:.0f}]  "
              f"sustained={results['conditions'][cond]['n_sustained']}/{n_seeds}")

    # fold-increase headline on mean_ps_window (guard against zero ground value).
    g = results["conditions"]["ground"]["mean_ps_window_mean"]
    u = results["conditions"]["microgravity"]["mean_ps_window_mean"]
    if g > 1e-9:
        results["fold_increase"] = round(u / g, 3)
        results["fold_increase_note"] = "microgravity / ground mean_ps_window"
    else:
        results["fold_increase"] = None
        results["fold_increase_note"] = (
            f"ground mean_ps_window ~0 ({g:.3g}); fold undefined, "
            f"microgravity mean_ps_window={u:.3g}")
    results["total_wall_seconds"] = round(time.time() - t_start, 1)

    save_ps_figure(curves, S2_DELAY_MS, total_ms)
    merge_results("ensemble", results)
    print(f"[ensemble] total wall={results['total_wall_seconds']}s  "
          f"fold_increase={results['fold_increase']}")
    print("Wrote figures/ensemble_ps.png, snapshot_crn_{ground,microgravity}.png, "
          "results_crn.json[ensemble]")


if __name__ == "__main__":
    main()
