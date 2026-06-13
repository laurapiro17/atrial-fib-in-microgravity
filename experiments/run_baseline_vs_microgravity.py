"""Main experiment: does microgravity remodelling promote fibrillatory re-entry?

For each condition we seed a single broken wavefront (one rotor), evolve the tissue,
and track the number of phase singularities (rotors) over time.  The hypothesis:

    ground tissue stabilises a single organised rotor, whereas the microgravity-
    remodelled atrium fragments it into multiple, persistent wavelets -- fibrillation.

The quantitative headline is taken from a **fibrosis-seed ensemble** at full duration,
not a single run: the microgravity substrate is random, and its vulnerability is
strongly substrate-dependent (some realisations tip into runaway fibrillation, others
sustain only a few extra wavelets).  Reporting the ensemble mean +- std, and the range,
is the honest summary.  Rotor counts are area-normalised (the dilated atrium is larger)
so we compare rotor *density*, not raw counts.

Outputs (written to ``figures/``):
  * ``anim_<condition>.mp4``   -- the membrane-potential movie (microgravity = a
                                  high-burden representative seed, for illustration)
  * ``snapshot_<condition>.png`` -- a representative late frame
  * ``ps_timeseries.png``      -- rotor count vs time, both conditions
  * ``results.json``           -- ensemble statistics and the area-normalised headline

Run:  ``python experiments/run_baseline_vs_microgravity.py``  (a few minutes, headless)
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from afib_microgravity import (  # noqa: E402
    AtrialSheet,
    count_phase_singularities,
    phase_field,
    seed_broken_wavefront,
)
from afib_microgravity.remodeling import make_condition, microgravity_params, \
    fibrosis_field  # noqa: E402
from afib_microgravity import viz  # noqa: E402

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
DT = 0.02
N_STEPS = 9000          # ~180 model time-units
RECORD_EVERY = 60       # frame / metric cadence
ENSEMBLE_SEEDS = [7, 100, 101, 102, 103, 104]
ANIM_SEED = 7           # the representative (high-burden) seed shown in the movie


def evolve(sheet, n_steps, record_every, keep_frames):
    """Evolve a sheet, returning (ps_series, frames)."""
    ps_series, frames = [], []
    for i in range(n_steps):
        sheet.step()
        if i % record_every == 0:
            ps_series.append(count_phase_singularities(phase_field(sheet.u, sheet.v)))
            if keep_frames:
                frames.append(sheet.u.copy())
    return ps_series, frames


def _second_half_mean(series):
    return float(np.mean(series[len(series) // 2:]))


def run_ground():
    """Ground atrium is deterministic; one run suffices."""
    shape, params, d_field = make_condition("ground")
    sheet = AtrialSheet(shape=shape, params=params, D_field=d_field, dt=DT)
    seed_broken_wavefront(sheet)
    t0 = time.time()
    ps_series, frames = evolve(sheet, N_STEPS, RECORD_EVERY, keep_frames=True)
    area = shape[0] * shape[1]
    summary = {
        "condition": "ground", "shape": list(shape), "area": area,
        "ps_max": int(max(ps_series)),
        "ps_mean_second_half": _second_half_mean(ps_series),
        "ps_density_x1e4": _second_half_mean(ps_series) / area * 1e4,
        "wall_seconds": round(time.time() - t0, 1),
    }
    print(f"[ground] {summary}")
    return ps_series, frames, summary


def run_microgravity_ensemble():
    """Microgravity atrium across fibrosis realisations (random substrate)."""
    base = (220, 220)
    shape = (int(base[0] * 1.3), int(base[1] * 1.3))
    area = shape[0] * shape[1]
    params = microgravity_params(autonomic_severity=1.0)
    per_seed, anim_frames, anim_series = [], None, None
    for seed in ENSEMBLE_SEEDS:
        d_field = fibrosis_field(shape, density=0.55, seed=seed)
        sheet = AtrialSheet(shape=shape, params=params, D_field=d_field, dt=DT)
        seed_broken_wavefront(sheet)
        keep = seed == ANIM_SEED
        ps_series, frames = evolve(sheet, N_STEPS, RECORD_EVERY, keep_frames=keep)
        m = _second_half_mean(ps_series)
        per_seed.append({"seed": seed, "ps_mean_second_half": m,
                         "ps_max": int(max(ps_series))})
        print(f"  [microgravity seed={seed}] mean2nd={m:.2f} max={max(ps_series)}")
        if keep:
            anim_frames, anim_series = frames, ps_series
    means = [d["ps_mean_second_half"] for d in per_seed]
    summary = {
        "condition": "microgravity", "shape": list(shape), "area": area,
        "n_seeds": len(ENSEMBLE_SEEDS), "per_seed": per_seed,
        "ps_mean_second_half": float(np.mean(means)),
        "ps_std_second_half": float(np.std(means)),
        "ps_min_second_half": float(np.min(means)),
        "ps_max_second_half": float(np.max(means)),
        "ps_density_x1e4": float(np.mean(means)) / area * 1e4,
    }
    print(f"[microgravity] mean2nd={summary['ps_mean_second_half']:.2f}"
          f"+/-{summary['ps_std_second_half']:.2f}")
    return anim_series, anim_frames, summary


def main():
    os.makedirs(FIG, exist_ok=True)
    results = {}

    ground_series, ground_frames, results["ground"] = run_ground()
    print("Running microgravity fibrosis-seed ensemble (full duration) ...")
    ug_series, ug_frames, results["microgravity"] = run_microgravity_ensemble()

    # animations + snapshots
    viz.animate(ground_frames, os.path.join(FIG, "anim_ground.mp4"),
                fps=20, title="ground atrium")
    viz.save_snapshot(ground_frames[-1], os.path.join(FIG, "snapshot_ground.png"),
                      title="ground atrium (late)")
    viz.animate(ug_frames, os.path.join(FIG, "anim_microgravity.mp4"),
                fps=20, title=f"microgravity atrium (seed {ANIM_SEED})")
    viz.save_snapshot(ug_frames[-1], os.path.join(FIG, "snapshot_microgravity.png"),
                      title="microgravity atrium (late)")
    viz.plot_ps_timeseries(
        {"Ground atrium": ground_series,
         f"Microgravity atrium (seed {ANIM_SEED})": ug_series},
        os.path.join(FIG, "ps_timeseries.png"), dt=DT, record_every=RECORD_EVERY)

    # area-normalised headline (rotor density, not raw count)
    g_density = results["ground"]["ps_density_x1e4"]
    u_density = results["microgravity"]["ps_density_x1e4"]
    results["headline"] = {
        "metric": "phase-singularity density (rotors per 1e4 cells), 2nd-half mean",
        "ground": round(g_density, 3),
        "microgravity_mean": round(u_density, 3),
        "density_fold_increase": round(u_density / max(g_density, 1e-9), 2),
        "microgravity_seed_range_rotors": [
            round(results["microgravity"]["ps_min_second_half"], 2),
            round(results["microgravity"]["ps_max_second_half"], 2),
        ],
        "note": ("High substrate-dependent variance: a subset of fibrosis patterns "
                 "tips into sustained fibrillation while others sustain few extra "
                 "wavelets. This variance is itself the finding."),
    }

    with open(os.path.join(FIG, "results.json"), "w") as f:
        json.dump(results, f, indent=2)
    print("HEADLINE:", json.dumps(results["headline"], indent=2))
    print("Wrote figures + results.json to", os.path.abspath(FIG))


if __name__ == "__main__":
    main()
