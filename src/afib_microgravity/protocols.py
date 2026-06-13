"""Stimulation protocols used to probe the tissue.

The clinically meaningful one is the **S1-S2 cross-field protocol**: a planar wave
(S1) is launched, and a second stimulus (S2) is delivered to one half of the sheet
while the S1 wave is recovering.  The S2 wave can only propagate into already-recovered
tissue, so its wavefront is broken -- and a broken wavefront curls into a rotating
spiral, the engine of re-entrant arrhythmia.
"""

from __future__ import annotations

import numpy as np

from .model import AtrialSheet


def seed_broken_wavefront(sheet: AtrialSheet, refractory_v: float = 2.0):
    """Initialise a single broken wavefront -- the canonical robust spiral seed.

    The left half is depolarised (a wavefront sits at the midline) while the top half
    is held refractory.  The wavefront can only advance through the recovered bottom
    half, so its free upper end curls into a rotor.  This deterministic initial
    condition avoids the fragile S1-S2 timing search and gives every run exactly one
    starting rotor, so any *extra* rotors that appear are attributable to the tissue
    substrate, not the stimulus.
    """
    sheet.u[:, : sheet.nx // 2] = 1.0
    sheet.v[: sheet.ny // 2, :] = refractory_v


def planar_s1(sheet: AtrialSheet, width: int = 5):
    """Launch a planar wave from the left edge."""
    mask = np.zeros((sheet.ny, sheet.nx), dtype=bool)
    mask[:, :width] = True
    sheet.stimulate(mask)


def cross_field_s2(sheet: AtrialSheet):
    """Deliver an S2 stimulus to the lower half of the sheet."""
    mask = np.zeros((sheet.ny, sheet.nx), dtype=bool)
    mask[sheet.ny // 2:, :] = True
    sheet.stimulate(mask)


def induce_spiral(sheet: AtrialSheet, s1_s2_delay: int, record_every: int = 0,
                  n_after_s2: int = 0):
    """Run the S1-S2 protocol to seed a spiral wave.

    Parameters
    ----------
    s1_s2_delay : int
        Steps between the S1 and S2 stimuli.  This coupling interval is what decides
        whether a spiral forms; it must land in the "vulnerable window" where the S1
        wave's tail is still refractory over part of the sheet.
    record_every : int
        If > 0, store a copy of ``u`` every this-many steps and return the frames.
    n_after_s2 : int
        Steps to run after the S2 stimulus.

    Returns
    -------
    list[ndarray]
        Recorded frames (empty if ``record_every == 0``).
    """
    frames = []

    def maybe_record(i):
        if record_every and i % record_every == 0:
            frames.append(sheet.u.copy())

    planar_s1(sheet)
    for i in range(s1_s2_delay):
        sheet.step()
        maybe_record(i)
    cross_field_s2(sheet)
    for i in range(n_after_s2):
        sheet.step()
        maybe_record(i)
    return frames
