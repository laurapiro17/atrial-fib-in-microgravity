"""Fast, deterministic tests for the model, metrics and remodelling mappings."""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from afib_microgravity import (  # noqa: E402
    APParams,
    AtrialSheet,
    count_phase_singularities,
    fibrosis_field,
    planar_conduction_velocity,
    planar_s1,
)
from afib_microgravity.model import AtrialSheet as Sheet  # noqa: E402


def test_resting_state_is_stable():
    """With no stimulus the tissue must stay at rest (no spurious activation)."""
    s = AtrialSheet(shape=(40, 40), dt=0.02)
    s.run(500)
    assert np.allclose(s.u, 0.0, atol=1e-6)
    assert np.allclose(s.v, 0.0, atol=1e-6)


def test_state_stays_bounded():
    """u must remain physiological (roughly within [0, 1.1]) -- a stability guard."""
    s = AtrialSheet(shape=(60, 60), dt=0.02)
    planar_s1(s)
    s.run(800)
    assert s.u.max() < 1.1
    assert s.u.min() > -0.2


def test_planar_wave_propagates():
    """A left-edge stimulus must travel rightward (activate tissue beyond the stimulus).

    Conduction is slow in dimensionless units, so we track the peak activation reached
    at a mid-sheet probe over the run rather than demanding the far edge in a few steps.
    """
    s = AtrialSheet(shape=(40, 120), dt=0.02)
    planar_s1(s)
    probe = s.nx // 2
    peak = 0.0
    for _ in range(4000):
        s.step()
        peak = max(peak, s.u[20, probe])
        if peak > 0.5:
            break
    assert peak > 0.5  # the wave reached the middle of the sheet


def test_uniform_diffusion_matches_standard_laplacian():
    """On a uniform D field the FV operator equals the 5-point Laplacian."""
    s = AtrialSheet(shape=(32, 32))
    rng = np.random.default_rng(0)
    field = rng.standard_normal((32, 32))
    fv = s._div_D_grad(field)
    # standard 5-point laplacian with no-flux (edge-replicate) boundaries
    lap = (
        np.pad(field, 1, mode="edge")[2:, 1:-1]
        + np.pad(field, 1, mode="edge")[:-2, 1:-1]
        + np.pad(field, 1, mode="edge")[1:-1, 2:]
        + np.pad(field, 1, mode="edge")[1:-1, :-2]
        - 4 * field
    )
    assert np.allclose(fv, lap, atol=1e-10)


def test_conduction_velocity_drops_with_coupling():
    """Lower diffusion (fibrosis-like) must slow conduction -- a sanity regression."""
    def fast():
        return Sheet(shape=(40, 160), params=APParams(D=1.0), dt=0.02)

    def slow():
        return Sheet(shape=(40, 160), params=APParams(D=0.5),
                     dt=0.02)

    cv_fast = planar_conduction_velocity(fast)
    cv_slow = planar_conduction_velocity(slow)
    assert np.isfinite(cv_fast) and np.isfinite(cv_slow)
    assert cv_slow < cv_fast


def test_phase_singularity_counter_on_synthetic_rotor():
    """A single analytic spiral phase field must be detected as exactly one rotor."""
    n = 80
    y, x = np.mgrid[0:n, 0:n]
    # centre the rotor off-grid (half-integer) so it sits inside a plaquette, as a
    # real rotor generically does -- avoids the exact +-pi loop degeneracy.
    phase = np.arctan2(y - (n / 2 - 0.5), x - (n / 2 - 0.5))
    assert count_phase_singularities(phase) == 1


def test_phase_singularity_counter_on_flat_field():
    """A non-winding (planar gradient) phase field has no rotors."""
    n = 64
    _, x = np.mgrid[0:n, 0:n]
    phase = (x / n) * 1.0  # monotone, no winding
    assert count_phase_singularities(phase) == 0


def test_fibrosis_field_density():
    """Higher requested density yields a lower mean diffusivity."""
    homog = fibrosis_field((100, 100), density=0.0)
    fibro = fibrosis_field((100, 100), density=0.5, seed=1)
    assert homog.mean() == pytest.approx(1.0)
    assert fibro.mean() < homog.mean()


def test_broken_wavefront_seeds_one_rotor():
    """The broken-wavefront seed must produce exactly one rotor on healthy tissue."""
    from afib_microgravity import seed_broken_wavefront, phase_field
    from afib_microgravity import count_phase_singularities as count_ps

    s = AtrialSheet(shape=(160, 160), dt=0.02)
    seed_broken_wavefront(s)
    s.run(1500)  # let the free end curl into a rotor
    ps = count_ps(phase_field(s.u, s.v))
    assert ps == 1


def test_dt_stability_guard():
    """An over-large dt must be rejected up front rather than blowing up silently."""
    with pytest.raises(ValueError):
        AtrialSheet(shape=(20, 20), dt=1.0, params=APParams(D=1.0))
