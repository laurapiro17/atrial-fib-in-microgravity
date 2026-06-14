"""Sustainability probe: S1-S2 cross-field induction on a CRN monodomain sheet,
run long enough (>=1500-2000 ms) to judge whether re-entry sustains.

Reports per condition: PS time series summary (mean over the last third, whether
PS stays > 0), max PS, and a SUSTAIN verdict, plus a late-state snapshot PNG to
figures/ so a human can eyeball whether a real rotor is present.

This module also exposes ``s1s2_cross_field`` and ``run_induction`` so the
ensemble/production script can reuse the exact induction.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import argparse
import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from afib_microgravity.crn import CRNCell  # noqa: E402
from afib_microgravity.diffusion import AnisotropicDiffusion  # noqa: E402
from afib_microgravity.fibrosis import correlated_fibrosis  # noqa: E402
from afib_microgravity.metrics import count_phase_singularities  # noqa: E402
from afib_microgravity.model import MonodomainSheet  # noqa: E402
from afib_microgravity.remodeling import (  # noqa: E402
    baseline_crn_params,
    microgravity_crn_params,
)

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")

V_DEPOL = 20.0      # mV written into stimulated cells
V0 = -40.0          # mV phase-plane origin (mid-upstroke)
DELAY_MS = 30.0     # time-delay for the recovery proxy
DT = 0.02
D_LONG = 0.15       # calibrated: planar CV ~58 cm/s
D_TRANS = 0.05      # anisotropy ratio ~3:1 (atrial-like)


def crn_phase(V_now, V_delayed):
    return np.arctan2(V_now - V0, V_delayed - V0)


def s1_planar(cell, width=4):
    ny, nx = cell.shape
    mask = np.zeros((ny, nx), dtype=bool)
    mask[:, :width] = True
    cell.stimulate(mask, V_DEPOL)


def s2_quadrant(cell):
    """S2 over the lower-left quadrant: it can only propagate into recovered
    tissue, breaking the wavefront -> a rotor."""
    ny, nx = cell.shape
    mask = np.zeros((ny, nx), dtype=bool)
    mask[ny // 2:, : nx // 2] = True
    cell.stimulate(mask, V_DEPOL)


def make_cell_diff(condition, shape, seed):
    if condition == "ground":
        params = baseline_crn_params()
        coupling = None
    else:
        params = microgravity_crn_params(severity=1.0)
        coupling = correlated_fibrosis(shape, density=0.3, corr_len=4.0, seed=seed)
    cell = CRNCell(shape=shape, params=params, use_numba=True)
    diff = AnisotropicDiffusion(shape, D_LONG, D_TRANS,
                               theta=np.zeros(shape), dx=0.25, coupling=coupling)
    return cell, diff


def run_induction(condition, shape, seed, total_ms, s2_delay_ms,
                  record_every_ms=10.0, snapshot=None):
    """S1 planar, wait s2_delay_ms, S2 quadrant, then evolve to total_ms.
    Returns (ps_series, times_ms, wall_seconds)."""
    cell, diff = make_cell_diff(condition, shape, seed)
    sheet = MonodomainSheet(cell, diff, dt=DT)

    rec_steps = max(1, round(record_every_ms / DT))
    # The phase proxy needs V and a ~DELAY_MS-old copy. We snapshot only at the
    # record cadence and read the snapshot ``delay_snaps`` records back, so the
    # ring holds O(DELAY/record) frames, not O(DELAY/dt) -- a small, fixed memory
    # footprint and no per-step copies.
    delay_snaps = max(1, round(DELAY_MS / record_every_ms))
    s2_step = round(s2_delay_ms / DT)
    n_steps = round(total_ms / DT)

    ring = []          # recent V snapshots (most recent last)
    ps_series, times = [], []
    t0 = time.perf_counter()

    s1_planar(cell)
    for i in range(n_steps):
        if i == s2_step:
            s2_quadrant(cell)
        sheet.step()
        if i % rec_steps == 0:
            ring.append(cell.V.copy())
            if len(ring) > delay_snaps + 1:
                ring.pop(0)
            V_delayed = ring[0] if len(ring) > delay_snaps else ring[0]
            ps = count_phase_singularities(crn_phase(cell.V, V_delayed))
            ps_series.append(ps)
            times.append(i * DT)
    wall = time.perf_counter() - t0

    if snapshot:
        os.makedirs(FIG, exist_ok=True)
        fig, ax = plt.subplots(figsize=(4.2, 4))
        im = ax.imshow(cell.V, cmap="inferno", vmin=-85, vmax=20, origin="lower")
        ax.set_title(f"{condition}: V at t={total_ms:.0f} ms")
        fig.colorbar(im, ax=ax, label="V (mV)", shrink=0.8)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, snapshot), dpi=120)
        plt.close(fig)

    return ps_series, times, wall


def summarise(ps_series):
    s = np.asarray(ps_series, dtype=float)
    n = s.size
    if n == 0:
        return dict(max=0, last_third_mean=0.0, last_third_frac_pos=0.0, sustains=False)
    k = max(1, n // 3)
    last = s[-k:]
    return dict(
        max=int(s.max()),
        last_third_mean=float(last.mean()),
        last_third_frac_pos=float((last > 0).mean()),
        sustains=bool(last.mean() >= 1.0 and (last > 0).mean() >= 0.8),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", type=int, default=200)
    ap.add_argument("--total", type=float, default=2000.0, help="ms")
    ap.add_argument("--s2delay", type=float, default=160.0, help="ms")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--conditions", default="ground,microgravity")
    args = ap.parse_args()

    shape = (args.grid, args.grid)
    size_cm = args.grid * 0.025
    print(f"[probe] grid={shape} ({size_cm:.1f}x{size_cm:.1f} cm) total={args.total}ms "
          f"s2delay={args.s2delay}ms seed={args.seed}", flush=True)

    for cond in args.conditions.split(","):
        ps, t, wall = run_induction(
            cond, shape, args.seed, args.total, args.s2delay,
            snapshot=f"sustain_{cond}.png")
        summ = summarise(ps)
        print(f"\n=== {cond} ===", flush=True)
        print(f"  wall={wall:.1f}s  max_PS={summ['max']}  "
              f"last3rd_mean={summ['last_third_mean']:.2f}  "
              f"last3rd_frac_pos={summ['last_third_frac_pos']:.2f}  "
              f"SUSTAINS={summ['sustains']}", flush=True)
        # coarse trajectory
        if ps:
            idx = np.linspace(0, len(ps) - 1, min(10, len(ps))).astype(int)
            traj = [(round(t[i]), ps[i]) for i in idx]
            print(f"  PS(t): {traj}", flush=True)


if __name__ == "__main__":
    main()
