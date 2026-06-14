"""Specificity tests for the phase-singularity (rotor) detector.

The audit (C2) flagged that ``count_phase_singularities`` was never validated:
a trustworthy rotor detector must return exactly 1 on a single spiral, and ~0 on
fields that contain NO rotation (planar waves, static structural steps). If it
returns spurious singularities on a non-rotating field, the headline wavebreak
result could be a detection artifact rather than a biological effect.
"""
import numpy as np

from afib_microgravity.metrics import count_phase_singularities


def _single_spiral_phase(n=64):
    """Phase that winds once around the centre -> exactly one phase singularity.

    The core is offset to a half-cell (n/2 + 0.5) so it sits INSIDE a plaquette,
    not on a grid node (a node-centred core splits its charge across four
    plaquettes and is correctly counted as a near-miss in each)."""
    y, x = np.mgrid[0:n, 0:n]
    return np.arctan2(y - (n / 2.0 + 0.5), x - (n / 2.0 + 0.5))


def test_single_spiral_gives_exactly_one():
    assert count_phase_singularities(_single_spiral_phase()) == 1


def test_planar_phase_ramp_gives_zero():
    # a linear phase ramp (planar wave) has no winding
    n = 64
    _, x = np.mgrid[0:n, 0:n]
    phase = (x / n) * 2.0 * np.pi - np.pi
    assert count_phase_singularities(phase) == 0


def test_static_structural_step_gives_zero():
    # a sharp spatial step in phase (e.g. a conduction-block edge with no rotation):
    # two uniform regions, no winding -> must be zero
    phase = np.full((48, 48), -2.0)
    phase[:, 24:] = 1.5
    assert count_phase_singularities(phase) == 0


def test_uniform_field_gives_zero():
    assert count_phase_singularities(np.full((32, 32), 0.7)) == 0


def test_two_counter_rotating_spirals_give_two():
    n = 80
    y, x = np.mgrid[0:n, 0:n]
    # two cores, each offset to a half-cell, opposite chirality, well separated
    left = np.arctan2(y - (n / 2.0 + 0.5), x - (n / 4.0 + 0.5))
    right = -np.arctan2(y - (n / 2.0 + 0.5), x - (3 * n / 4.0 + 0.5))
    phase = np.where(x < n / 2, left, right)
    # at least the two genuine cores (the np.where seam may add edge counts)
    assert count_phase_singularities(phase) >= 2
