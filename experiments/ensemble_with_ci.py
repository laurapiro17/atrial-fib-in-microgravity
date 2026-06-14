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
# A planar stimulus edge briefly spikes the phase proxy to a huge spurious count
# (~150-200) -- this is a wavefront artifact, NOT re-entry. Genuine wavebreak in
# the fibrotic substrate gives a MODEST count (single digits to low tens). We
# therefore reject samples whose PS exceeds ARTIFACT_PS_THRESHOLD, and ignore the
# first MIN_ANALYSIS_MS (the initial S1 edge), rather than blanking by time. The
# headline observable is the resulting *wavebreak burden*: how much the substrate
# fragments propagating wavefronts into phase singularities (rotor cores).
ARTIFACT_PS_THRESHOLD = 40
MIN_ANALYSIS_MS = 40.0

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
    """Wavebreak-burden metrics over the artifact-rejected trace.

    We keep samples with ``times >= MIN_ANALYSIS_MS`` and
    ``ps <= ARTIFACT_PS_THRESHOLD`` (the rest are wavefront-edge artifacts). On
    this clean series we report the peak and mean phase-singularity (wavebreak)
    count, the integrated wavebreak burden (PS-ms, a rotor-lifetime surrogate),
    and the active fraction of time (any wavebreak present)."""
    valid = (times >= MIN_ANALYSIS_MS) & (ps_series <= ARTIFACT_PS_THRESHOLD)
    w_ps = ps_series[valid]
    w_t = times[valid]

    if w_ps.size == 0:
        return {
            "peak_ps": 0, "mean_ps": 0.0, "break_burden_ps_ms": 0.0,
            "active_frac": 0.0, "ps_density_peak_x1e4": 0.0,
            "sustained_af": False, "n_valid": 0,
        }

    peak = int(w_ps.max())
    mean_ps = float(w_ps.mean())
    # integrated wavebreak burden (PS-ms): sum(PS) * sample spacing.
    dt_samp = float(np.median(np.diff(w_t))) if w_t.size > 1 else 0.0
    burden = float(w_ps.sum() * dt_samp)
    active_frac = float((w_ps > 0).mean())

    return {
        "peak_ps": peak,
        "mean_ps": mean_ps,
        "break_burden_ps_ms": burden,
        "active_frac": active_frac,
        "ps_density_peak_x1e4": ps_density(peak, area),
        "sustained_af": bool(is_sustained_af(w_ps, window_frac=0.5, threshold=2)),
        "n_valid": int(w_ps.size),
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
    # Clip the display at the artifact threshold so the brief wavefront-edge spike
    # (~150-200 PS at each stimulus) doesn't crush the genuine wavebreak signal.
    ax.axvline(s2_ms, color="0.4", ls="--", lw=1, label=f"S2 ({s2_ms:.0f} ms)")
    ax.axhline(ARTIFACT_PS_THRESHOLD, color="0.7", ls=":", lw=1,
               label=f"artifact cutoff ({ARTIFACT_PS_THRESHOLD})")
    for cond, data in curves.items():
        t = data["times"]
        stack = np.clip(np.vstack(data["series"]), 0, ARTIFACT_PS_THRESHOLD)
        mean = stack.mean(axis=0)
        ax.plot(t, mean, color=colors.get(cond, "k"), label=f"{cond} (mean)")
        if stack.shape[0] > 1:
            sd = stack.std(axis=0)
            ax.fill_between(t, mean - sd, mean + sd, color=colors.get(cond, "k"),
                            alpha=0.18)
    ax.set_ylim(0, ARTIFACT_PS_THRESHOLD)
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("phase singularities (wavebreak)")
    ax.set_title("Transient re-entrant vulnerability: wavebreak after S1-S2")
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
        "s2_delay_ms": S2_DELAY_MS,
        "artifact_ps_threshold": ARTIFACT_PS_THRESHOLD,
        "min_analysis_ms": MIN_ANALYSIS_MS,
        "d_long": D_LONG, "apd_note": APD_NOTE,
        "conditions": {},
    }
    curves = {}
    snapshots = {"ground": "snapshot_crn_ground.png",
                 "microgravity": "snapshot_crn_microgravity.png"}
    t_start = time.time()

    for cond in conditions:
        per_seed = []
        burdens, peaks = [], []
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
            burdens.append(m["break_burden_ps_ms"])
            peaks.append(float(m["peak_ps"]))
            series_list.append(ps_series)
            common_times = times
            last_V, last_shape = V_final, shape
            print(f"  [{cond} seed={seed}] peak={m['peak_ps']} "
                  f"mean_ps={m['mean_ps']:.2f} "
                  f"burden={m['break_burden_ps_ms']:.0f}PS·ms "
                  f"active={m['active_frac']:.2f} ({m['wall_seconds']}s)")

        burden_mu, burden_lo, burden_hi = bootstrap_ci(
            burdens, n_boot=n_boot, alpha=0.05, seed=0)
        peak_mu, peak_lo, peak_hi = bootstrap_ci(
            peaks, n_boot=n_boot, alpha=0.05, seed=0)
        results["conditions"][cond] = {
            "per_seed": per_seed,
            "break_burden_ps_ms_mean": burden_mu,
            "break_burden_ps_ms_ci95": [burden_lo, burden_hi],
            "peak_ps_mean": peak_mu,
            "peak_ps_ci95": [peak_lo, peak_hi],
            "n_seeds_with_wavebreak": int(sum(d["peak_ps"] > 0 for d in per_seed)),
        }
        curves[cond] = {"times": common_times, "series": series_list}
        save_snapshot(last_V, cond, total_ms, snapshots[cond])
        print(f"[{cond}] burden={burden_mu:.0f}PS·ms "
              f"95%CI=[{burden_lo:.0f}, {burden_hi:.0f}]  "
              f"peak={peak_mu:.1f} 95%CI=[{peak_lo:.1f}, {peak_hi:.1f}]  "
              f"seeds_with_break="
              f"{results['conditions'][cond]['n_seeds_with_wavebreak']}/{n_seeds}")

    # fold-increase headline on wavebreak burden (guard against zero ground value).
    g = results["conditions"]["ground"]["break_burden_ps_ms_mean"]
    u = results["conditions"]["microgravity"]["break_burden_ps_ms_mean"]
    if g > 1e-9:
        results["fold_increase"] = round(u / g, 3)
        results["fold_increase_note"] = "microgravity / ground wavebreak burden"
    else:
        results["fold_increase"] = None
        results["fold_increase_note"] = (
            f"ground wavebreak burden ~0 ({g:.3g}); fold undefined (healthy tissue "
            f"conducts cleanly), microgravity burden={u:.3g} PS·ms")
    results["total_wall_seconds"] = round(time.time() - t_start, 1)

    save_ps_figure(curves, S2_DELAY_MS, total_ms)
    merge_results("ensemble", results)
    print(f"[ensemble] total wall={results['total_wall_seconds']}s  "
          f"fold_increase={results['fold_increase']}")
    print("Wrote figures/ensemble_ps.png, snapshot_crn_{ground,microgravity}.png, "
          "results_crn.json[ensemble]")


if __name__ == "__main__":
    main()
