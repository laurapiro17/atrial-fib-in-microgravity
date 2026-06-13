"""Quantitative read-outs of arrhythmic activity.

The single most informative scalar is the number of **phase singularities** (PS) --
the rotation centres of spiral waves.  One PS is a single re-entrant rotor; many,
constantly created and annihilated, is the signature of fibrillation.  Counting PS
over time turns "the tissue looks chaotic" into a number we can compare between
conditions.
"""

from __future__ import annotations

import numpy as np


def phase_field(u, v, u0: float = 0.3, v0: float = 0.3):
    """Map the (u, v) state to a phase angle in (-pi, pi].

    A spiral rotor is a point around which this phase winds through a full 2*pi, so
    phase is the natural coordinate for detecting rotors.
    """
    return np.arctan2(v - v0, u - u0)


def count_phase_singularities(phase):
    """Count phase singularities via the topological-charge (loop-integral) method.

    Around each 2x2 plaquette we sum the wrapped phase increments.  A closed loop that
    encircles a rotor accumulates +-2*pi; everywhere else it sums to ~0.
    """
    def wrap(d):
        return (d + np.pi) % (2.0 * np.pi) - np.pi

    p00 = phase[:-1, :-1]
    p01 = phase[:-1, 1:]
    p11 = phase[1:, 1:]
    p10 = phase[1:, :-1]
    loop = wrap(p01 - p00) + wrap(p11 - p01) + wrap(p10 - p11) + wrap(p00 - p10)
    return int(np.sum(np.abs(loop) > np.pi))


def planar_conduction_velocity(sheet_factory, length: int = 120, dt_record: int = 1,
                               threshold: float = 0.5):
    """Measure planar conduction velocity (cells per unit time).

    Stimulates the left edge and times the wavefront crossing two probes.  Returned in
    dimensionless model units; useful as a regression check that coupling changes move
    CV in the expected direction.

    Parameters
    ----------
    sheet_factory : callable -> AtrialSheet
        Builds a fresh sheet (so the function has no side effects on caller state).
    """
    sheet = sheet_factory()
    from .protocols import planar_s1  # local import to avoid a cycle

    x1, x2 = sheet.nx // 4, 3 * sheet.nx // 4
    row = sheet.ny // 2
    t1 = t2 = None
    planar_s1(sheet)
    for _ in range(length * 50):
        sheet.step()
        if t1 is None and sheet.u[row, x1] > threshold:
            t1 = sheet.t
        if t2 is None and sheet.u[row, x2] > threshold:
            t2 = sheet.t
            break
    if t1 is None or t2 is None or t2 <= t1:
        return float("nan")
    return (x2 - x1) / (t2 - t1)


def dominant_frequency(signal, dt):
    """Dominant frequency of a single-cell ``u(t)`` trace (model units^-1)."""
    sig = np.asarray(signal, dtype=float)
    sig = sig - sig.mean()
    if sig.size < 4 or not np.any(sig):
        return float("nan")
    spectrum = np.abs(np.fft.rfft(sig))
    freqs = np.fft.rfftfreq(sig.size, d=dt)
    return float(freqs[1:][np.argmax(spectrum[1:])])
