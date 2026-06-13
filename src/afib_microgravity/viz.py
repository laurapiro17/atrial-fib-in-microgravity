"""Rendering helpers: still snapshots and animations of the membrane potential."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless / overnight-safe

import matplotlib.animation as animation  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402


def save_snapshot(u, path, title=None):
    """Save a single ``u`` field as a PNG."""
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    ax.imshow(u, cmap="inferno", vmin=0.0, vmax=1.0, interpolation="nearest")
    ax.set_axis_off()
    if title:
        ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def animate(frames, path, fps: int = 20, title=None):
    """Write a list of ``u`` fields to an mp4 (falls back to gif if ffmpeg absent)."""
    if not frames:
        raise ValueError("no frames to animate")
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    im = ax.imshow(frames[0], cmap="inferno", vmin=0.0, vmax=1.0,
                   interpolation="nearest")
    ax.set_axis_off()
    if title:
        ax.set_title(title)
    fig.tight_layout()

    def update(k):
        im.set_data(frames[k])
        return (im,)

    anim = animation.FuncAnimation(fig, update, frames=len(frames),
                                   interval=1000 / fps, blit=True)
    try:
        anim.save(path, writer="ffmpeg", fps=fps, dpi=110)
    except Exception:
        gif_path = path.rsplit(".", 1)[0] + ".gif"
        anim.save(gif_path, writer="pillow", fps=fps)
        path = gif_path
    plt.close(fig)
    return path


def plot_ps_timeseries(series_by_label, path, dt, record_every):
    """Plot phase-singularity count vs time for each condition."""
    fig, ax = plt.subplots(figsize=(7, 4))
    for label, series in series_by_label.items():
        t = [i * record_every * dt for i in range(len(series))]
        ax.plot(t, series, label=label, lw=2)
    ax.set_xlabel("time (model units)")
    ax.set_ylabel("phase singularities (rotor count)")
    ax.set_title("Re-entrant activity: ground vs microgravity atrium")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)
