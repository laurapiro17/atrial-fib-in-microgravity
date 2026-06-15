"""Render the CRN 'money-shot': side-by-side ground vs microgravity V(t) after the
same S1-S2 induction, as an animated GIF for the README.

Ground conducts the wavefront cleanly; the microgravity fibrotic substrate fragments
it into rotor cores. Frames are written as PNGs and assembled with ffmpeg (2-pass
palette) -> figures/anim_crn.gif.

Run:  python experiments/make_crn_animation.py [--grid 160 --total-ms 600]
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

import argparse
import os
import subprocess
import sys
import tempfile

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from afib_microgravity.model import MonodomainSheet  # noqa: E402
from afib_microgravity.remodeling import make_condition_crn  # noqa: E402
from ensemble_with_ci import S2_DELAY_MS, s1_planar, s2_quadrant  # noqa: E402

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
DT = 0.02


def run_frames(grid, total_ms, frame_every_ms, seed):
    """Return dict cond -> list of V arrays (frames)."""
    out = {}
    for cond in ("ground", "microgravity"):
        cell, diff = make_condition_crn(cond, base_shape=(grid, grid), seed=seed)
        sheet = MonodomainSheet(cell, diff, dt=DT)
        n_steps = round(total_ms / DT)
        s2_step = round(S2_DELAY_MS / DT)
        frame_every = max(1, round(frame_every_ms / DT))
        # pad ground to microgravity (dilated) shape for a uniform canvas
        frames = []
        s1_planar(cell)
        for i in range(n_steps):
            if i == s2_step:
                s2_quadrant(cell)
            sheet.step()
            if i % frame_every == 0:
                frames.append(cell.V.copy())
        out[cond] = frames
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", type=int, default=160)
    ap.add_argument("--total-ms", type=float, default=600.0, dest="total_ms")
    ap.add_argument("--frame-every-ms", type=float, default=8.0, dest="frame_every_ms")
    ap.add_argument("--seed", type=int, default=4)
    args = ap.parse_args()

    print(f"[anim] grid={args.grid} total={args.total_ms}ms")
    frames = run_frames(args.grid, args.total_ms, args.frame_every_ms, args.seed)
    n = min(len(frames["ground"]), len(frames["microgravity"]))

    tmp = tempfile.mkdtemp()
    for k in range(n):
        fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.9))
        for ax, cond in zip(axes, ("ground", "microgravity")):
            ax.imshow(frames[cond][k], cmap="inferno", vmin=-85, vmax=20, origin="lower")
            ax.set_title(f"{cond}", fontsize=11)
            ax.set_xticks([]); ax.set_yticks([])
        t_ms = k * args.frame_every_ms
        fig.suptitle(f"CRN atrium, S1-S2 induction   t = {t_ms:.0f} ms", fontsize=11)
        fig.tight_layout()
        fig.savefig(os.path.join(tmp, f"f{k:04d}.png"), dpi=90)
        plt.close(fig)

    os.makedirs(FIG, exist_ok=True)
    pal = os.path.join(tmp, "pal.png")
    gif = os.path.join(FIG, "anim_crn.gif")
    subprocess.run(["ffmpeg", "-y", "-v", "error", "-framerate", "12",
                    "-i", os.path.join(tmp, "f%04d.png"),
                    "-vf", "palettegen=max_colors=64:stats_mode=diff", pal], check=True)
    subprocess.run(["ffmpeg", "-y", "-v", "error", "-framerate", "12",
                    "-i", os.path.join(tmp, "f%04d.png"), "-i", pal,
                    "-lavfi", "paletteuse=dither=bayer:bayer_scale=4", gif], check=True)
    size_kb = os.path.getsize(gif) // 1024
    print(f"[anim] wrote {gif} ({size_kb} KB, {n} frames)")


if __name__ == "__main__":
    main()
