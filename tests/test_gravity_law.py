import math
import pytest
from afib_microgravity.gravity_law import (
    gravity_to_remodeling, SEVERITY_MAX, DILATION_MAX, DENSITY_MAX,
)


def test_earth_recovers_baseline():
    r = gravity_to_remodeling(1.0)
    assert r.severity == pytest.approx(0.0)
    assert r.dilation == pytest.approx(1.0)
    assert r.density == pytest.approx(0.0)


def test_freefall_recovers_microgravity_operating_point():
    r = gravity_to_remodeling(0.0)
    assert r.severity == pytest.approx(SEVERITY_MAX)
    assert r.dilation == pytest.approx(DILATION_MAX)
    assert r.density == pytest.approx(DENSITY_MAX)


def test_severity_monotonic_decreasing_in_gravity():
    sev = [gravity_to_remodeling(g).severity for g in (0.0, 0.25, 0.5, 0.75, 1.0)]
    assert all(a >= b for a, b in zip(sev, sev[1:]))


def test_gravity_clamped_to_unit_interval():
    assert gravity_to_remodeling(1.5).severity == pytest.approx(0.0)
    assert gravity_to_remodeling(-0.3).severity == pytest.approx(SEVERITY_MAX)


from afib_microgravity.gravity_law import (
    wavelength_cm, cardiac_gravitational_number, interpolate_crossing,
)


def test_wavelength_is_cv_times_apd_in_cm():
    # 58 cm/s * 0.300 s = 17.4 cm
    assert wavelength_cm(300.0, cv_cm_s=58.0) == pytest.approx(17.4)


def test_Ng_is_one_when_path_equals_wavelength():
    # choose apd so wl=8 cm, l0=8, dilation=1 -> N_g=1
    apd = 8.0 / 58.0 * 1000.0
    ng = cardiac_gravitational_number(apd, dilation=1.0, l0_cm=8.0, cv_cm_s=58.0)
    assert ng == pytest.approx(1.0)


def test_Ng_rises_when_wavelength_shortens():
    long_wl = cardiac_gravitational_number(300.0, dilation=1.0)
    short_wl = cardiac_gravitational_number(150.0, dilation=1.0)
    assert short_wl > long_wl


def test_interpolate_crossing_finds_linear_root():
    xs = [0.0, 0.2, 0.4, 0.6]
    ys = [2.0, 1.5, 0.5, 0.2]   # crosses 1.0 between 0.2 and 0.4
    gstar = interpolate_crossing(xs, ys, target=1.0)
    assert gstar == pytest.approx(0.3)


def test_interpolate_crossing_returns_none_without_crossing():
    assert interpolate_crossing([0.0, 1.0], [3.0, 2.0], target=1.0) is None
